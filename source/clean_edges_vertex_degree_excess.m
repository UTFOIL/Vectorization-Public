function [ original_edge_indices ] = clean_edges_vertex_degree_excess( edges2vertices, number_of_edges_per_vertex )
%% clean_edge_cycles
% SAM 7/1/19

number_of_edges = size( edges2vertices, 1 );

original_edge_indices = 1 : number_of_edges ;

edges2vertices = double( edges2vertices );

number_of_vertices = max( edges2vertices( : ));

number_of_edges_to_remove = 0 ;

% The first vertex in the pair is given by the row, second by the column.
edge_lookup_table   = sparse( edges2vertices( :, 1 ),   ...
                              edges2vertices( :, 2 ),   ...
                              ( 1 : number_of_edges )', ...
                              number_of_vertices,       ...
                              number_of_vertices,       ...
                              number_of_edges           );

% make sparse adjacency matrix
adjacency_matrix    = sparse( edges2vertices( :, 1 ),                ...
                              edges2vertices( :, 2 ),                ...
                              ones( number_of_edges, 1, 'logical' ), ...
                              number_of_vertices,                    ...
                              number_of_vertices,                    ...
                              number_of_edges                        );
    
vertex_degrees = sum( adjacency_matrix, 1 ) + sum( adjacency_matrix, 2 )';

vertex_excess_degrees = vertex_degrees - number_of_edges_per_vertex ;

vertices_of_excess_degree = find( vertex_excess_degrees > 0 );

edges_to_remove = [ ];

for vertex_index = vertices_of_excess_degree
        
    % find the vertex coordinates of all edges attached to this vertex of excess degree
    vertex_pairs_to_remove_1A = find( edge_lookup_table( :, vertex_index ) );
    vertex_pairs_to_remove_2B = find( edge_lookup_table( vertex_index, : )');
    
	vertex_pairs_to_remove =   [ vertex_pairs_to_remove_1A; vertex_index * ones( size( vertex_pairs_to_remove_2B )) ]                            ...
                           + ( [ vertex_index * ones( size( vertex_pairs_to_remove_1A )); vertex_pairs_to_remove_2B ] - 1 ) * number_of_vertices ;    
    
	edges_at_vertex = full( edge_lookup_table( vertex_pairs_to_remove ));

%     [ edges_at_vertex_descending, IA ] = sort( edges_at_vertex, 'descend' );
    edges_at_vertex_descending = sort( edges_at_vertex, 'descend' );
    
%     vertex_pairs_to_remove = vertex_pairs_to_remove( IA );
    
    ordinates_of_edges_to_remove = 1 : vertex_excess_degrees( vertex_index );
        
%      adjacency_matrix( vertex_pairs_to_remove( ordinates_of_edges_to_remove )) = 0 ;
%     edge_lookup_table( vertex_pairs_to_remove( ordinates_of_edges_to_remove )) = 0 ;

%     [ I, J ] = ind2sub( number_of_vertices([ 1, 1 ]), vertex_pairs_to_remove( ordinates_of_edges_to_remove ));
% 
%     symmetric_vertex_pairs_to_remove = sub2ind( number_of_vertices([ 1, 1 ]), J, I );
% 
%     adjacency_matrix( symmetric_vertex_pairs_to_remove ) = 0 ;
        
	edges_to_remove( number_of_edges_to_remove + ordinates_of_edges_to_remove ) ...
                   = edges_at_vertex_descending( ordinates_of_edges_to_remove ) ;

    number_of_edges_to_remove = number_of_edges_to_remove + vertex_excess_degrees( vertex_index );
        
end % FOR vertices_of_excess_degree

original_edge_indices( edges_to_remove ) = [ ];    

end % FUCNCTION clean_edges_vertex_excess_degree
