function [ strand_subscripts, strand_energies ]                                                                  ...
        = get_strand_objects( edge_subscripts, edge_energies, edge_indices_in_strands, edge_backwards_in_strands )
%% get_strand_objects SAM 4/23/19
% The purpose of this function is to loop through the strand objects, putting the edges together
% into strands as encoded by the inputs edge_indices_in_strands and edge_backwards_in_strands.  

number_of_strands = length( edge_indices_in_strands );

strand_index_range = 1 : number_of_strands ;

% count the number of edges in each strand, output into vector with length of the number of strands
numbers_of_edges_in_strands                                                                         ...
            = cellfun( @( indices_at_strand ) size( indices_at_strand, 1 ), edge_indices_in_strands );

strand_subscripts = cell( size( edge_indices_in_strands ));
strand_energies   = cell( size( edge_indices_in_strands ));

parfor strand_index = strand_index_range

    % Assembly into Strands from Edges:
    %
    % look up the edge_subscripts associated with the edges in this strand. The variable
    % edge_subscripts_at_strand_in_edges is a cell array of vector lists, one list for each edge in
    % the current strand. The edges are called in their order of appearance in the strand, but
    % within each edge the list is as likely to be backwards of the master list as it is to be in
    % the same order. So we will have to flip about half of them so that there is a preserved order
    % in listing between edges. The edges to be flipped are encoded in the edge_backwards_in_strands
    % input.
    edge_subscripts_at_strand_in_edges = edge_subscripts( edge_indices_in_strands{ strand_index });
    edge_energies_at_strand_in_edges   = edge_energies(   edge_indices_in_strands{ strand_index });    
    
    edge_index_in_strand_range = 1 : numbers_of_edges_in_strands( strand_index );
    
    % loop through the edges in the curent strand
    for edge_index_in_strand = edge_index_in_strand_range
        
        % flip this edge if needed
        if edge_backwards_in_strands{ strand_index }( edge_index_in_strand )
        
                      edge_subscripts_at_strand_in_edges{ edge_index_in_strand } ...
            = flipud( edge_subscripts_at_strand_in_edges{ edge_index_in_strand });
        
                        edge_energies_at_strand_in_edges{ edge_index_in_strand } ...
            = flipud(   edge_energies_at_strand_in_edges{ edge_index_in_strand });        
            
        end % backwards edge IF
    
        % trim the last vector of each edge except on the last edge (to avoid double listing of the
        % vertices
        if edge_index_in_strand < numbers_of_edges_in_strands( strand_index )
            
              edge_subscripts_at_strand_in_edges{ edge_index_in_strand }                 ...
            = edge_subscripts_at_strand_in_edges{ edge_index_in_strand }( 1 : end - 1, : );
        
                edge_energies_at_strand_in_edges{ edge_index_in_strand }                 ...
            =   edge_energies_at_strand_in_edges{ edge_index_in_strand }( 1 : end - 1, : );        
        
        end % not last edge IF
    end % edge in strand FOR                                   
    
    % combine all the subscripts from all the edges into a single ordered list representing the
    % current strand
    strand_subscripts{ strand_index } = cell2mat( edge_subscripts_at_strand_in_edges );
    strand_energies{   strand_index } = cell2mat(   edge_energies_at_strand_in_edges );
    
%     % erase any repeated subscripts (for later interpolation we need unique locations).  Repeated
%     % subscripts can arise if a parent/child edge doesn't end up being a bifurcation, but part of a
%     % single strand.
%     [ strand_subscripts{ strand_index }, unique_indices ]                                           ...
%                                       = unique( strand_subscripts{ strand_index }, 'rows', 'stable' );
%                                   
%     strand_energies{ strand_index } = strand_energies{ strand_index }( unique_indices );
    
    % erase any repeated subscripts (for later interpolation we need unique locations).  Repeated
    % subscripts can arise if a parent/child edge doesn't end up being a bifurcation, but part of a
    % single strand.  Since this is after some smoothing at the edge step, we must specify the
    % number of digits that two values must be identical to be considered a repeat. Choosing one
    % digit past the decimal point.
    [ ~, unique_indices ] = unique( round( 10 * strand_subscripts{ strand_index }( :, 1 : 3 )), 'rows', 'stable' );
                       
    strand_subscripts{ strand_index } = strand_subscripts{ strand_index }( unique_indices, : );   
    strand_energies{   strand_index } = strand_energies{   strand_index }( unique_indices    );
    
    
end % strand FOR
end % FUNCTION