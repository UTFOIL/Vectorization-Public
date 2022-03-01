function [ vertex_indices_in_strands, edge_indices_in_strands_out, edge_backwards_in_strands ]      ...
                                        = sort_network_V180( edges2vertices, end_vertices_of_strands, ...
                                                             edge_indices_in_strands                )
%% sort_network_V180 SAM 7/11/18
% the purpose of this function is to sort the output of the get_network_V182 function so that the
% edges are not only grouped by strand assignment (the output of get_network_V182) but are also
% ordered such that within each cell array entry of edge_indices_in_... or vertex_indices_in_... the
% objects are called in the physical order that they would appear if you started at one end of the
% strand/bifurcation and moved to the other end.
%
% Note that this function could be improved by combining strands and junctions when this will not
% lead to a bifurcation. SAM 7/13/18

%% Strand sorting
% (same as junction sorting)
number_of_strands = length( end_vertices_of_strands );

edge_indices_in_strands_out = edge_indices_in_strands ;

vertex_indices_in_strands = cell( size( edge_indices_in_strands ));
edge_backwards_in_strands = cell( size( edge_indices_in_strands ));

% strand_average_sizes = cellfun( @( x1 ) mean( cell2mat( edge_scale_subscripts( x1 ))), edge_indices_in_strands );

strand_index_range = 1 : number_of_strands ;

for strand_index = strand_index_range
    
    % find the starting vertex for the strand from the set of starting vertices from the edges that
    % are in the current strand
    vertices_of_edges_at_strand = edges2vertices( edge_indices_in_strands{ strand_index }, : );
    
    try
        
        % if this doesn't work then the strand is cyclical (has no ending vertices because it is a loop)        
        next_index_of_strand = find( vertices_of_edges_at_strand( : ) == end_vertices_of_strands{ strand_index }( 1 ), 1 );
        
    catch
        
        % start anywhere
        next_index_of_strand = 1 ;
        
    end
    
    % loop through the edges in the strand and build the sorted lists from the unsorted ones
    number_of_edges_in_strand = length( edge_indices_in_strands{ strand_index });
        
    edge_index_range = 1 : number_of_edges_in_strand ;
        
    for strand_edge_index = edge_index_range
        
        % if the index is in the second column of vertices_of_edges_at_strand, then the edge needs to be
        % flipped to fit head to tail in the strand (keep track of this with the
        % edge_backwards_in_strands variable)
        if next_index_of_strand > number_of_edges_in_strand

            next_edge_origin = 2 ; next_edge_terminus = 1 ;

            next_index_of_strand = next_index_of_strand - number_of_edges_in_strand ;

        else

            next_edge_origin = 1 ; next_edge_terminus = 2 ;

        end % IF edge backwards
        
        vertex_indices_in_strands{ strand_index }( strand_edge_index ) = vertices_of_edges_at_strand( next_index_of_strand, next_edge_origin );
                
        next_edge_index = edge_indices_in_strands{ strand_index }( next_index_of_strand );
        
        edge_indices_in_strands_out{ strand_index }( strand_edge_index ) = next_edge_index ;
        
        edge_backwards_in_strands{ strand_index }( strand_edge_index ) = next_edge_terminus == 1 ;
        
        ending_vertex_of_next_index = vertices_of_edges_at_strand( next_index_of_strand, next_edge_terminus );
        
        % zero out the vertex addresses for the current edge to avoid finding it in the next line
        vertices_of_edges_at_strand( next_index_of_strand, : ) = 0 ;
                  
        next_index_of_strand = find( vertices_of_edges_at_strand( : ) == ending_vertex_of_next_index );
        
    end % FOR edges in strands
    
    vertex_indices_in_strands{ strand_index }( strand_edge_index + 1 ) = ending_vertex_of_next_index ;
    
end % FOR strands

end % FUNCTION