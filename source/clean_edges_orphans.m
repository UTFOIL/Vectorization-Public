function [ original_edge_indices ] = clean_edges_orphans( edge_space_subscripts, size_of_image, vertex_space_subscripts )
%% clean_edge_cycles
% SAM 7/1/19

number_of_edges = size( edge_space_subscripts, 1 );

original_edge_indices = 1 : number_of_edges ;

%% Remove child (grandchild, (etc.)) edges who do not connect to any other edge
% This section function should be called before the cycle-breaking section, if edges were erased
% before that section.  In that case we would want to eliminate any orphaned connecitons from the
% adjacency matrix before removing cycles and possibly removing non-orphaned edges before an
% orphaned, destined to be removed edge.
% 
% Note: In the add vertices section we again force children to have higher (worse) energy than their
% parents.

edge_space_subscripts_matrix = double( cell2mat( edge_space_subscripts ));

vertex_space_subscripts_double = double( vertex_space_subscripts );

edge_locations   = sub2ind( size_of_image,   edge_space_subscripts_matrix( :, 1 ),   edge_space_subscripts_matrix( :, 2 ),   edge_space_subscripts_matrix( :, 3 ));
vertex_locations = sub2ind( size_of_image, vertex_space_subscripts_double( :, 1 ), vertex_space_subscripts_double( :, 2 ), vertex_space_subscripts_double( :, 3 ));

edge_lengths = cellfun( @( x ) size( x, 1 ), edge_space_subscripts );

edge_locations_cell = mat2cell( edge_locations, edge_lengths, 1 );

searching_for_orphans = true ;

while searching_for_orphans
    
    edge_locations = cell2mat( edge_locations_cell );
    
    [ ~, unique_edge_location_indices ] = unique( edge_locations );

    % throw out all of the unique ones and save the repeats
    edge_location_repeats = edge_locations ; edge_location_repeats( unique_edge_location_indices ) = [ ];
    
    % redefine the edge look up table because the number of edge list numbering has changed
    number_of_edges       = length( edge_locations_cell );    
    
    edge_index_range_cell = num2cell( 1 : number_of_edges )';
    
    edge_index_LUT      = cellfun( @( x, y ) x * ones( size( y )), edge_index_range_cell, edge_locations_cell, 'UniformOutput', false  );    
    
    exterior_edge_locations                 = cell2mat( cellfun( @( x ) x( [ 1, end  ] ), edge_locations_cell, 'UniformOutput', false ));
    exterior_edge_location_index2edge_index = cell2mat( cellfun( @( x ) x( [ 1, end  ] ), edge_index_LUT     , 'UniformOutput', false ));    
    
    % look for edge terminal locations that don't coincide with a vertex or any other edge location
    [ ~, orphan_terminal_indices ] = setdiff( exterior_edge_locations( : ), union( edge_location_repeats, vertex_locations ));
    
    edge_indices_to_remove = exterior_edge_location_index2edge_index( orphan_terminal_indices );
        
    % erase these edges
    original_edge_indices( edge_indices_to_remove    ) = [ ];    
      edge_locations_cell( edge_indices_to_remove    ) = [ ];
%            edges2vertices( edge_indices_to_remove, : ) = [ ];
%     edge_space_subscripts( edge_indices_to_remove    ) = [ ];      
%     edge_scale_subscripts( edge_indices_to_remove    ) = [ ];
%             edge_energies( edge_indices_to_remove    ) = [ ];
            
    searching_for_orphans = logical( length( edge_indices_to_remove ));
    
end % WHILE searching for orphans

end % FUNCTION clean_edges_orphans
