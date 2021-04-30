function [ original_edge_indices ] = clean_edges_hairs( edges2vertices )
%% clean_edge_hairs
% SAM 7/18/19

number_of_edges = size( edges2vertices, 1 );

original_edge_indices = 1 : number_of_edges ;

edges2vertices = double( edges2vertices );

number_of_vertices = max( edges2vertices( : ));

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
                           
% force symmetric adjacency_matrix
adjacency_matrix = adjacency_matrix | adjacency_matrix' ;

% eliminate vertices that are endpoints and neighboring bifurcations

neighbors = sum( adjacency_matrix );

   endpoint_vertices = neighbors == 1 ;
bifurcation_vertices = neighbors >= 3 ;

adjacency_matrix_without_bifurcations = adjacency_matrix ;

adjacency_matrix_without_bifurcations( bifurcation_vertices, : ) = 0 ;
adjacency_matrix_without_bifurcations( :, bifurcation_vertices ) = 0 ;

neighbors_without_bifurcations = sum( adjacency_matrix_without_bifurcations );

no_non_bifurcation_neighbor_vertices = neighbors_without_bifurcations == 0 ;

endpoint_vertices_neighboring_bifurcations = endpoint_vertices & no_non_bifurcation_neighbor_vertices ;

% find THE associated edge to delete.
edges_to_delete = [ nonzeros( edge_lookup_table( endpoint_vertices_neighboring_bifurcations, : )); ...
                    nonzeros( edge_lookup_table( :, endpoint_vertices_neighboring_bifurcations ))  ];

original_edge_indices( edges_to_delete ) = [ ];
    
% %% eliminating cycles of three vertices by checking for cycles and removing the weakest link:
% 
% vertices_to_be_truncated = sum( adjacency_matrix ) <= 1 ;
% 
% % remove the vertices that can never be in cycles
%  adjacency_matrix( vertices_to_be_truncated, : ) = [ ];
% edge_lookup_table( vertices_to_be_truncated, : ) = [ ];
% 
%  adjacency_matrix( :, vertices_to_be_truncated ) = [ ];
% edge_lookup_table( :, vertices_to_be_truncated ) = [ ];    
% 
% checking_for_cycle = true ;
% 
% cumulative_number_of_cycle_conncomps = 0 ;
% 
% while checking_for_cycle
%     
%     number_of_truncated_vertices = size( adjacency_matrix, 1 );        
%         
%     % make two-step adjacency matrix (nodes that are connected by two edges)
%     two_step_adjacency_matrix = adjacency_matrix ^ 2 ;
%     
%     cycle_adjacency_matrix = two_step_adjacency_matrix & adjacency_matrix ;
%     
%     cycle_flag_by_vertex = sum( cycle_adjacency_matrix );    
%     
%     checking_for_cycle = any( cycle_flag_by_vertex );
%     
%     % convert the adjacency matrix to a graph to find the connected components
%     cycle_graph = graph( cycle_adjacency_matrix );    
%     
%     % The connected components of this graph (with three or more nodes, impossible to have two
%     % nodes) are the vertex groups entangled in cycles.
%     vertices_in_cycles = conncomp( cycle_graph, 'OutputForm', 'cell' );
%     
%     nodes_per_cycle = cellfun( @length, vertices_in_cycles );
%     
%     % eliminate "cycles" of one node from this removal FOR loop (these have no edges to remove)
%     vertices_in_cycles( nodes_per_cycle == 1 ) = [ ];
%     
%     number_of_cycle_conncomps = length( vertices_in_cycles );
%         
%     vertex_pairs_to_remove = zeros( number_of_cycle_conncomps, 2 );
%     
%     cycle_group_index_range = 1 : number_of_cycle_conncomps ;
%     
%     parfor cycle_group_index = cycle_group_index_range
%                 
%         possible_vertex_pairs_for_cycle                                                        ...
%                         =   vertices_in_cycles{ cycle_group_index }'                           ...
%                         + ( vertices_in_cycles{ cycle_group_index } - 1 ) * number_of_truncated_vertices ;
%         
%         edge_indices_in_cycles = nonzeros( edge_lookup_table( possible_vertex_pairs_for_cycle ));
%         
%         % Mark the worst edge from this cyclical object. (Edge indices are sorted by energy.)
%         edges_to_remove( cumulative_number_of_cycle_conncomps + cycle_group_index ) = max( edge_indices_in_cycles );
%         
% %         vertex_pairs_to_remove =   edges2vertices( edges_to_remove( number_of_edges_to_remove ), [ 1, 2 ])                            ...
% %                                + ( edges2vertices( edges_to_remove( number_of_edges_to_remove ), [ 2, 1 ]) - 1 ) * number_of_vertices ;
%         
%         [ vertex_pairs_to_remove_1, vertex_pairs_to_remove_2 ] = find( edge_lookup_table == edges_to_remove( cumulative_number_of_cycle_conncomps + cycle_group_index ), 1 );
% 
%         vertex_pairs_to_remove( cycle_group_index, : )                                 ...
%                                =  [ vertex_pairs_to_remove_1, vertex_pairs_to_remove_2 ]                                      ...
%                                + ([ vertex_pairs_to_remove_2, vertex_pairs_to_remove_1 ] - 1 ) * number_of_truncated_vertices ;
%         
%         
%     end % FOR cycle group
%     
%     % Remove the worst edge from each cycle conncomp from the adjacency matrix and edge_lookup_table
%      adjacency_matrix( vertex_pairs_to_remove ) = 0 ;
%     edge_lookup_table( vertex_pairs_to_remove ) = 0 ;
%     
%     cumulative_number_of_cycle_conncomps = cumulative_number_of_cycle_conncomps + number_of_cycle_conncomps ;
%     
%     % record the vertices that don't need to be included in this analysis because they will never
%     % again be in a cycle
%     vertices_to_be_truncated = 1 : number_of_truncated_vertices ;
%         
%     vertices_to_be_truncated( cell2mat( vertices_in_cycles )) = [ ];    
%     
%     % remove the vertices that will never be in cycles again
%      adjacency_matrix( vertices_to_be_truncated, : ) = [ ];
%     edge_lookup_table( vertices_to_be_truncated, : ) = [ ];
%     
%      adjacency_matrix( :, vertices_to_be_truncated ) = [ ];
%     edge_lookup_table( :, vertices_to_be_truncated ) = [ ];    
%     
% end % WHILE checking_for_cycle
%        
% original_edge_indices( edges_to_remove    ) = [ ];
% 
% end % FUNCTION clean_edge_cycles


