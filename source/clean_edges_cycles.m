function [ original_edge_indices ] = clean_edges_cycles( edges2vertices )
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
                           
% force symmetric adjacency_matrix
adjacency_matrix = adjacency_matrix | adjacency_matrix' ;
    
%% eliminating cycles of three vertices by checking for cycles and removing the weakest link:

vertices_to_be_truncated = sum( adjacency_matrix ) <= 1 ;

% remove the vertices that can never be in cycles
 adjacency_matrix( vertices_to_be_truncated, : ) = [ ];
edge_lookup_table( vertices_to_be_truncated, : ) = [ ];

 adjacency_matrix( :, vertices_to_be_truncated ) = [ ];
edge_lookup_table( :, vertices_to_be_truncated ) = [ ];

edges_to_remove = zeros( numel( nonzeros( edge_lookup_table )), 1 );

checking_for_cycle = true ;

cumulative_number_of_cycle_conncomps = 0 ;

while checking_for_cycle
    
    number_of_truncated_vertices = size( adjacency_matrix, 1 );        
        
    % make two-step adjacency matrix (nodes that are connected by two edges)
    two_step_adjacency_matrix = adjacency_matrix ^ 2 ;
    
    cycle_adjacency_matrix = two_step_adjacency_matrix & adjacency_matrix ;
    
    cycle_flag_by_vertex = sum( cycle_adjacency_matrix );    
    
    checking_for_cycle = any( cycle_flag_by_vertex );
    
    % convert the adjacency matrix to a graph to find the connected components
    cycle_graph = graph( cycle_adjacency_matrix );    
    
    % The connected components of this graph (with three or more nodes, impossible to have two
    % nodes) are the vertex groups entangled in cycles.
    vertices_in_cycles = conncomp( cycle_graph, 'OutputForm', 'cell' );
    
    nodes_per_cycle = cellfun( @length, vertices_in_cycles );
    
    % eliminate "cycles" of one node from this removal FOR loop (these have no edges to remove)
    vertices_in_cycles( nodes_per_cycle == 1 ) = [ ];
    
    number_of_cycle_conncomps = length( vertices_in_cycles );
        
    vertex_pairs_to_remove = zeros( number_of_cycle_conncomps, 2 );
    
    cycle_group_index_range = 1 : number_of_cycle_conncomps ;
    
    % changed from PARFOR to FOR 6/8/20
    for cycle_group_index = cycle_group_index_range
                
        possible_vertex_pairs_for_cycle                                                        ...
                        =   vertices_in_cycles{ cycle_group_index }'                           ...
                        + ( vertices_in_cycles{ cycle_group_index } - 1 ) * number_of_truncated_vertices ;
        
        edge_indices_in_cycles = nonzeros( edge_lookup_table( possible_vertex_pairs_for_cycle ));
        
        % Mark the worst edge from this cyclical object. (Edge indices are sorted by energy.)
        edges_to_remove( cumulative_number_of_cycle_conncomps + cycle_group_index ) = max( edge_indices_in_cycles );
        
%         vertex_pairs_to_remove =   edges2vertices( edges_to_remove( number_of_edges_to_remove ), [ 1, 2 ])                            ...
%                                + ( edges2vertices( edges_to_remove( number_of_edges_to_remove ), [ 2, 1 ]) - 1 ) * number_of_vertices ;
        
        [ vertex_pairs_to_remove_1, vertex_pairs_to_remove_2 ] = find( edge_lookup_table == edges_to_remove( cumulative_number_of_cycle_conncomps + cycle_group_index ), 1 );

        vertex_pairs_to_remove( cycle_group_index, : )                                 ...
                               =  [ vertex_pairs_to_remove_1, vertex_pairs_to_remove_2 ]                                      ...
                               + ([ vertex_pairs_to_remove_2, vertex_pairs_to_remove_1 ] - 1 ) * number_of_truncated_vertices ;
        
        
    end % FOR cycle group
    
    % Remove the worst edge from each cycle conncomp from the adjacency matrix and edge_lookup_table
     adjacency_matrix( vertex_pairs_to_remove ) = 0 ;
    edge_lookup_table( vertex_pairs_to_remove ) = 0 ;
    
    cumulative_number_of_cycle_conncomps = cumulative_number_of_cycle_conncomps + number_of_cycle_conncomps ;
    
    % record the vertices that don't need to be included in this analysis because they will never
    % again be in a cycle
    vertices_to_be_truncated = 1 : number_of_truncated_vertices ;
        
    vertices_to_be_truncated( cell2mat( vertices_in_cycles )) = [ ];    
    
    % remove the vertices that will never be in cycles again
     adjacency_matrix( vertices_to_be_truncated, : ) = [ ];
    edge_lookup_table( vertices_to_be_truncated, : ) = [ ];
    
     adjacency_matrix( :, vertices_to_be_truncated ) = [ ];
    edge_lookup_table( :, vertices_to_be_truncated ) = [ ];    
    
end % WHILE checking_for_cycle
       
original_edge_indices( nonzeros( edges_to_remove )) = [ ];

end % FUNCTION clean_edge_cycles
% 
% %% Remove child (grandchild, (etc.)) edges who do not connect to any other edge
% % !!! This section should be turned into a function and called before the cycle-breaking section.
% % Why? In case edges were erased before that section, we want to eliminate any orphaned connecitons
% % from the adjacency matrix before removing cycles and possibly removing non-orphaned edges before
% % an orphaned, destined to be removed edge.
% % 
% % Note: In the add vertices section we again force children to have higher (worse) energy than their
% % parents.
% 
% %    mean_edge_energies( edges_to_remove    ) = [ ];
%        edges2vertices( edges_to_remove, : ) = [ ];
% 
% edge_energies          =          edge_energies( original_edge_indices );
% edge_space_subscripts  =  edge_space_subscripts( original_edge_indices );
% edge_scale_subscripts  =  edge_scale_subscripts( original_edge_indices );
% % degrees_of_edges       =       degrees_of_edges( original_edge_indices );
% % edges            =            edges( original_edge_indices );
% 
% edges2vertices = uint32( edges2vertices );
% 
% edge_space_subscripts_matrix = double( cell2mat( edge_space_subscripts ));
% 
% vertex_space_subscripts_double = double( vertex_space_subscripts );
% 
% edge_locations   = sub2ind( size_of_image,   edge_space_subscripts_matrix( :, 1 ),   edge_space_subscripts_matrix( :, 2 ),   edge_space_subscripts_matrix( :, 3 ));
% vertex_locations = sub2ind( size_of_image, vertex_space_subscripts_double( :, 1 ), vertex_space_subscripts_double( :, 2 ), vertex_space_subscripts_double( :, 3 ));
% 
% edge_lengths = cellfun( @length, edge_scale_subscripts );
% 
% edge_locations_cell = mat2cell( edge_locations, edge_lengths, 1 );
% 
% searching_for_orphans = true ;
% 
% while searching_for_orphans
%     
%     edge_locations = cell2mat( edge_locations_cell );
%     
%     [ ~, unique_edge_location_indices ] = unique( edge_locations );
% 
%     % throw out all of the unique ones and save the repeats
%     edge_location_repeats = edge_locations ; edge_location_repeats( unique_edge_location_indices ) = [ ];
%     
%     number_of_edges       = size( edges2vertices, 1 );    
%     
%     edge_index_range_cell = num2cell( 1 : number_of_edges )';
%     
%     edge_index_LUT    = cellfun( @( x, y ) x * ones( size( y )), edge_index_range_cell, edge_scale_subscripts, 'UniformOutput', false  );    
%     
%     exterior_edge_locations                 = cell2mat( cellfun( @( x ) x( [ 1, end  ] ), edge_locations_cell, 'UniformOutput', false ));
%     exterior_edge_location_index2edge_index = cell2mat( cellfun( @( x ) x( [ 1, end  ] ), edge_index_LUT     , 'UniformOutput', false ));    
%     
%     % look for edge terminal locations that don't coincide with a vertex or any other edge location
%     [ ~, orphan_terminal_indices ] = setdiff( exterior_edge_locations( : ), union( edge_location_repeats, vertex_locations ));
%     
%     edge_indices_to_remove = exterior_edge_location_index2edge_index( orphan_terminal_indices );
%         
%     % erase these edges
%     original_edge_indices( edge_indices_to_remove    ) = [ ];    
%       edge_locations_cell( edge_indices_to_remove    ) = [ ];
%            edges2vertices( edge_indices_to_remove, : ) = [ ];
%     edge_space_subscripts( edge_indices_to_remove    ) = [ ];      
%     edge_scale_subscripts( edge_indices_to_remove    ) = [ ];
%             edge_energies( edge_indices_to_remove    ) = [ ];
%             
%     searching_for_orphans = logical( length( edge_indices_to_remove ));
%     
% end % WHILE searching for orphans
% 
% % !!!! this section has bugs: 
% %
% % 2) 
% 
% %% eliminating edges that put vertices above the maximum number_of_edges_per_vertex:
% % remove edges that put a vertex at too high a degree
% % while checking_for_excess_degree
% 
% number_of_edges = size( edges2vertices, 1 );
% 
% original_edge_indices = 1 : number_of_edges ;
% 
% edges2vertices = double( edges2vertices );
% 
% number_of_vertices = max( edges2vertices( : ));
% 
% % The first vertex in the pair is given by the row, second by the column.
% edge_lookup_table   = sparse( edges2vertices( :, 1 ),   ...
%                               edges2vertices( :, 2 ),   ...
%                               ( 1 : number_of_edges )', ...
%                               number_of_vertices,       ...
%                               number_of_vertices,       ...
%                               number_of_edges           );
% 
% % make sparse adjacency matrix
% adjacency_matrix    = sparse( edges2vertices( :, 1 ),                ...
%                               edges2vertices( :, 2 ),                ...
%                               ones( number_of_edges, 1, 'logical' ), ...
%                               number_of_vertices,                    ...
%                               number_of_vertices,                    ...
%                               number_of_edges                        );
%     
% vertex_degrees = sum( adjacency_matrix );
% 
% vertex_excess_degrees = vertex_degrees - number_of_edges_per_vertex ;
% 
% vertices_of_excess_degree = find( vertex_excess_degrees > 0 );
% 
% % vertex_sub_index = 0 ;
% 
% edge_lookup_table_temp = edge_lookup_table ;
% 
% for vertex_index = vertices_of_excess_degree
%     
% %     vertex_sub_index = vertex_sub_index + 1 ;
%     
% %     edges_at_vertex = nonzeros([ edge_lookup_table_temp( :, vertex_index ), edge_lookup_table_temp( vertex_index, : )']);
%     
%     % find the vertex coordinates of all edges attached to this vertex of excess degree
%     vertex_pairs_to_remove_1A = find( edge_lookup_table_temp( :, vertex_index ) );
%     vertex_pairs_to_remove_2B = find( edge_lookup_table_temp( vertex_index, : )');
%     
% 	vertex_pairs_to_remove =   [ vertex_pairs_to_remove_1A; vertex_index * ones( size( vertex_pairs_to_remove_2B )) ]                            ...
%                            + ( [ vertex_index * ones( size( vertex_pairs_to_remove_1A )); vertex_pairs_to_remove_2B ] - 1 ) * number_of_vertices ;    
%     
% 	edges_at_vertex = full( edge_lookup_table_temp( vertex_pairs_to_remove ));
% 
%     [ edges_at_vertex_descending, IA ] = sort( edges_at_vertex, 'descend' );
%     
%     vertex_pairs_to_remove = vertex_pairs_to_remove( IA );
%     
%     ordinates_of_edges_to_remove = 1 : vertex_excess_degrees( vertex_index );
%         
%      adjacency_matrix( vertex_pairs_to_remove( ordinates_of_edges_to_remove )) = 0 ;
%     edge_lookup_table( vertex_pairs_to_remove( ordinates_of_edges_to_remove )) = 0 ;
% 
%     [ I, J ] = ind2sub( number_of_vertices([ 1, 1 ]), vertex_pairs_to_remove( ordinates_of_edges_to_remove ));
% 
%     symmetric_vertex_pairs_to_remove = sub2ind( number_of_vertices([ 1, 1 ]), J, I );
% 
%     adjacency_matrix( symmetric_vertex_pairs_to_remove ) = 0 ;
%         
% 	edges_to_remove( number_of_edges_to_remove + ordinates_of_edges_to_remove ) ...
%                    = edges_at_vertex_descending( ordinates_of_edges_to_remove ) ;
% 
%     number_of_edges_to_remove = number_of_edges_to_remove + vertex_excess_degrees( vertex_index );
%         
% end % FOR vertices_of_excess_degree
% 
% original_edge_indices( edges_to_remove ) = [ ];    


