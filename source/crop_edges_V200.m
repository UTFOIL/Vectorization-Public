function [ excluded_edges_logical ] ...
                  = crop_edges_V200( edge_space_subscripts, edge_scale_subscripts, edge_energies,            ...
                                     lumen_radius_in_microns_range, microns_per_pixel,       ...
                                     size_of_image                                                           )
% %% choose_edges_V160
%                                    
% % this function selects the "best" edges of the possible edges found in the get_edges probabalistic
% % walk through the 4D (projected to 3D acrosss scale) energy field created in the get_energy
% % function.  For a given vertex pair A to B, the best trajectory connecting A to B is called the one
% % whose maximal energy attained along the trajectory is the least among all the trajectories from A
% % to B.  SAM 5/14/18
% %
% % V181: the previous version is combined with a function that resembles choose_vertices_V184 in that
% % it paints edges and looks at the painted image to determine conflicts.  Edges conflict with other
% % edges whenever they overlap away from the vertices that they connect.  All edges paint the same
% % Inf value. For each edge, after it paints itself it also paints the index of its start vertex
% % around that vertex and the index of its end vertex around its end.  To avoid conflict between
% % adjacent edges targeting the same node, the vertex influence should be larger than the edge
% % influence. SAM 5/30/18
% %
% % V182: edges don't paint, all that is painted is the vertices they start and stop on. The edges
% % still have an influence volume for searching for violations of vertices, but this volume could be
% % the same as that used for the vertices. SAM 5/31/18
% %
% % V190: mean edge energy is used instead of max edge energy as the summary statistic.  8/2/18 SAM
% %
% % function name changed from choose_edges_V190 to crop_edges_V200.  Keeping only the part of the
% % choose edges function that does the cropping and the part that chooses best unique trajectories
% % between vertices A and B mutually. eliminating references to sigma; instead referring to lumen
% % radius.  Splitting edge subscripts into space (uint16) and scale (uint8) cell vectors.  SAM
% % 11/14/18
% 
% %% Choosing best unique trajectories between vertices A and B mutually.
% 
% number_of_original_edges = length( edges2vertices );
% 
% original_edge_indices = 1 : number_of_original_edges ;
% 
% % % not possible to have these cases (A-A edge, A- edge, or edge with nonnegative max energy) with
% % the new deterministic edge search SAM 5/31/19
% % 
% % % remove trajectories that terminate on their original vertex
% % indices_of_non_self_referential_edges = find( edges2vertices( :, 1 ) ~= edges2vertices( :, 2 ));
% % 
% % edges2vertices          =          edges2vertices( indices_of_non_self_referential_edges, : );
% % % edge_sizes_temp       =       edge_subscripts( indices_of_non_self_referential_edges, 4 );
% % edge_energies_temp    =         edge_energies( indices_of_non_self_referential_edges    );
% % original_edge_indices = original_edge_indices( indices_of_non_self_referential_edges    );
% % 
% % % remove trajectories that don't land on any vertex
% % indices_of_terminal_edges = find( edges2vertices( :, 2 ) > 0 );
% % 
% % edges2vertices          =          edges2vertices( indices_of_terminal_edges, : );
% % % edge_sizes_temp       =       edge_sizes_temp( indices_of_terminal_edges    );
% % edge_energies_temp    =    edge_energies_temp( indices_of_terminal_edges    );
% % original_edge_indices = original_edge_indices( indices_of_terminal_edges    );
% % 
% % % remove edges that ever attained a negative energy on the trajectory
% % max_edge_energies = cellfun( @max,  edge_energies_temp );
% % 
% % indices_of_negative_energy_edges = find( max_edge_energies < 0 );
% 
% [ mean_edge_energies ] = get_edge_metric( edge_energies );
% 
% % 10/4/18 Possible Improvement: weight the edge means by a factor that is low when a certain
% % variance in the size coordinate of each trajectory is large.  This variance is first adjusted so
% % that any variance explained by a linear fit is removed:
% %
% % characteristic_size_stdev = scales_per_octave ; % needs to be input to this function
% %
% % % subtract off the constant term
% % edge_sizes_temp = cellfun( @(v) v - mean(v), edge_sizes_temp, 'UniformOutput', false )
% %
% % subtract off the (orthogonal) linear fits through the origin 
% % edge_sizes_temp = cellfun( @(v) v -       ((1:length(v))-(length(v)-1)/2)'    ...
% %                                     *     ((1:length(v))-(length(v)-1)/2)     ...
% %                                     * v/( ((1:length(v))-(length(v)-1)/2)     ...
% %                                          *((1:length(v))-(length(v)-1)/2)' ), ...
% %                             edge_sizes_temp, 'UniformOutput', false           )
% % 
% % % calculate the variances of the sizes, take negative exponent, and multiply this factor by the
% % % mean energy
% % mean_edge_energies = cellfun( @mean,  edge_energies_temp ) .* exp( - cellfun( @stdev, edge_sizes_temp ) / characteristic_size_stdev, edge_sizes_temp ));
% 
% % mean_edge_energies    =     mean_edge_energies( indices_of_negative_energy_edges    );
% % original_edge_indices = original_edge_indices(  indices_of_negative_energy_edges    );
% % edges2vertices        =   edges2vertices(       indices_of_negative_energy_edges, : );
% 
% % pre-sorting trajectories by length from shortest to longest, so that any A->B and B->A with the
% % same energy will go to the shorter one.  This is because the longer one may be double counting
% % some stretch.  Also, if you can make the trajectory in fewer steps with the same energy, why not?
% degrees_of_edges = cellfun( @length, edge_scale_subscripts );
% 
% [ ~, indices_sorted_by_length ] = sort( degrees_of_edges );
% 
% mean_edge_energies    =    mean_edge_energies( indices_sorted_by_length    );
% original_edge_indices = original_edge_indices( indices_sorted_by_length    );
% edges2vertices        =        edges2vertices( indices_sorted_by_length, : );
% 
% % sorting the trajectories by activation energy in ascending order, so lowest (best) energies at the
% % top of the list.
% [ mean_edge_energies, indices_sorted_by_mean_energy ] = sort( mean_edge_energies );
% 
% original_edge_indices = original_edge_indices(  indices_sorted_by_mean_energy    );
% edges2vertices        =        edges2vertices(  indices_sorted_by_mean_energy, : );
% 
% % % ! no longer possible to have multiple trajectories from A to B with deterministic edge search. !
% % 
% % % choosing only the best trajectory from vertex A to B for all A, B in the set of all vertices.
% % % UNIQUE chooses the first of the multiple instances, so top of the list/lowest energy is preferred.
% % [ edges2vertices, indices_of_unique_edges ] = unique( edges2vertices, 'rows', 'stable' );
% % 
% % original_edge_indices = original_edge_indices(  indices_of_unique_edges );
% % mean_edge_energies    =    mean_edge_energies( indices_of_unique_edges );
% 
% % choosing the worse of the two directions of edges in the cases where we have pairs of mutual
% % trajectories (from A to B and B to A).  Choosing the worse of the two, because the edge search is
% % deterministic.  There are no bad edges.  The one that is "worse" is going to be a child as opposed
% % to a parent edge in the cases where that is possible.  Favoring child edges because those will not
% % have overlapping edge locations.
% [ ~, mutual_edge_indices, reverse_mutual_edge_indices ]            ...
%     = intersect([ edges2vertices( :, 1 ), edges2vertices( :, 2 )], ...
%                 [ edges2vertices( :, 2 ), edges2vertices( :, 1 )], ...
%                 'rows', 'stable'                                   );
% 
% is_mutual_edge_pair_to_be_erased                        ...
%     = mutual_edge_indices < reverse_mutual_edge_indices ;  
% 
% original_edge_indices( mutual_edge_indices( is_mutual_edge_pair_to_be_erased )    ) = [ ];
%    mean_edge_energies( mutual_edge_indices( is_mutual_edge_pair_to_be_erased )    ) = [ ];
%        edges2vertices( mutual_edge_indices( is_mutual_edge_pair_to_be_erased ), : ) = [ ];      
%          
% %% eliminating cycles of three vertices by checking for cycles and removing the weakest link:
% 
% number_of_vertices = max( edges2vertices( : ));
% 
% number_of_edges = size( edges2vertices, 1 );
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
% % force symmetric adjacency_matrix
% adjacency_matrix = adjacency_matrix | adjacency_matrix' ;
% 
% checking_for_cycle = true ;
% 
% number_of_edges_to_remove = 0 ;
% 
% while checking_for_cycle
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
%     vertices_in_cycles( nodes_per_cycle == 1 ) = [ ];
%     
%     for cycle_group_index = 1 : length( vertices_in_cycles )
%         
%         number_of_edges_to_remove = number_of_edges_to_remove + 1 ;
%         
%         possible_vertex_pairs_for_cycle                                                        ...
%                         =   vertices_in_cycles{ cycle_group_index }'                           ...
%                         + ( vertices_in_cycles{ cycle_group_index } - 1 ) * number_of_vertices ;
%         
%         edge_indices_in_cycles = nonzeros( edge_lookup_table( possible_vertex_pairs_for_cycle ));
%         
%         % Mark the worst edge from this cyclical object. (Edge indices are sorted by energy.)
%         edges_to_remove( number_of_edges_to_remove ) = max( edge_indices_in_cycles );
%         
%         vertex_pairs_to_remove =   edges2vertices( edges_to_remove( number_of_edges_to_remove ), [ 1, 2 ])                            ...
%                                + ( edges2vertices( edges_to_remove( number_of_edges_to_remove ), [ 2, 1 ]) - 1 ) * number_of_vertices ;
%         
%         % Remove this edge from the adjacency matrix and the edge_lookup_table
%          adjacency_matrix( vertex_pairs_to_remove ) = 0 ;
%         edge_lookup_table( vertex_pairs_to_remove ) = 0 ;
%         
%     end % FOR cycle group
% end % WHILE checking_for_cycle
%        
% original_edge_indices( edges_to_remove    ) = [ ];
%    mean_edge_energies( edges_to_remove    ) = [ ];
%        edges2vertices( edges_to_remove, : ) = [ ];
%        
% 
% edge_energies          =          edge_energies( original_edge_indices );
% edge_space_subscripts  =  edge_space_subscripts( original_edge_indices );
% edge_scale_subscripts  =  edge_scale_subscripts( original_edge_indices );
% % degrees_of_edges       =       degrees_of_edges( original_edge_indices );
% % edges            =            edges( original_edge_indices );
% 
% %% add vertices where children edges meet their parents.  
% % Remove edges whose parents no longer exist. Work from best energy to worst energy (start at the
% % top of the list), so that the oldest parents are removed first, and all of their children,
% % grandchildren, etc.
% 
% 
% 
% 
% % end % function
% % 
% 
% 
% % !!!! move cropping to its own function, move above sections to clean_edges( ) to be executed
% % BEFORE smoothing
                                 
%% cropping 

mean_edge_energies = get_edge_metric( edge_energies );

degrees_of_edges = cellfun( @length, edge_scale_subscripts );

% pull the edge subscripts out of the cell array and concatenate into matrix.  Round them so that
% they can serve as indices into the image and the size lookup table.
edge_space_subscripts_matrix = uint16( cell2mat( edge_space_subscripts ));
edge_scale_subscripts_matrix = uint8(  cell2mat( edge_scale_subscripts ));

radii_in_pixels =  uint16( lumen_radius_in_microns_range( edge_scale_subscripts_matrix( :, 1 )) ./ microns_per_pixel );

subscript_maxs = edge_space_subscripts_matrix( :, 1 : 3 ) + radii_in_pixels ;
subscript_mins = edge_space_subscripts_matrix( :, 1 : 3 ) - radii_in_pixels ;

y_is_over_cell  = subscript_maxs( :, 1 ) > size_of_image( 1 );
x_is_over_cell  = subscript_maxs( :, 2 ) > size_of_image( 2 );
z_is_over_cell  = subscript_maxs( :, 3 ) > size_of_image( 3 );

y_is_under_cell = subscript_mins( :, 1 ) <         1         ;
x_is_under_cell = subscript_mins( :, 2 ) <         1         ;
z_is_under_cell = subscript_mins( :, 3 ) <         1         ;

z_has_max_scale = edge_scale_subscripts_matrix( : ) == length( lumen_radius_in_microns_range );
z_has_min_scale = edge_scale_subscripts_matrix( : ) ==                      1                 ;

% bug results from passing uint16 degrees of edges to mat2cell because it accumulates in
% the same data type and maxes out before reaching the sum of the degrees of edges

degrees_of_edges_uint_64 = uint64( degrees_of_edges );

y_is_over_cell       = mat2cell( y_is_over_cell,  degrees_of_edges_uint_64, 1 );
x_is_over_cell       = mat2cell( x_is_over_cell,  degrees_of_edges_uint_64, 1 );
z_is_over_cell       = mat2cell( z_is_over_cell,  degrees_of_edges_uint_64, 1 );

y_is_under_cell      = mat2cell( y_is_under_cell, degrees_of_edges_uint_64, 1 );
x_is_under_cell      = mat2cell( x_is_under_cell, degrees_of_edges_uint_64, 1 );
z_is_under_cell      = mat2cell( z_is_under_cell, degrees_of_edges_uint_64, 1 );

has_max_scale_cell = mat2cell( z_has_max_scale, degrees_of_edges_uint_64, 1 );
has_min_scale_cell = mat2cell( z_has_min_scale, degrees_of_edges_uint_64, 1 );

y_is_over    = cellfun( @max, y_is_over_cell );
x_is_over    = cellfun( @max, x_is_over_cell );
z_is_over    = cellfun( @max, z_is_over_cell ); 

y_is_under   = cellfun( @max, y_is_under_cell );
x_is_under   = cellfun( @max, x_is_under_cell );
z_is_under   = cellfun( @max, z_is_under_cell );

scale_is_max = cellfun( @max, has_max_scale_cell );
scale_is_min = cellfun( @max, has_min_scale_cell );

excluded_edges_logical =  y_is_over   |  x_is_over   | z_is_over  ...
                       |  y_is_under  |  x_is_under  | z_is_under ;...
        ...               | scale_is_max | scale_is_min ; % ... removed SAM 2/5/22

%  edge_space_subscripts( excluded_edges_logical    ) = [ ];
%  edge_scale_subscripts( excluded_edges_logical    ) = [ ];
%          edge_energies( excluded_edges_logical    ) = [ ];
% %       degrees_of_edges( excluded_edges_logical    ) = [ ];
%         edges2vertices( excluded_edges_logical, : ) = [ ];
%     mean_edge_energies( excluded_edges_logical    ) = [ ];
           
end % function