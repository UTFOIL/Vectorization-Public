function [ point_coordinates, strand_points ] = casx2vmv( point_coordinates, arc_connectivity, arc_diameters )
% SAM Sep. 10th, 2020 
% !!!!!!!!!!!!!!!! unfinished function! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

% ! point_coordinates as an output has an extra column for the radius

 
% point_uniques = unique( arc_connectivity );
% 
% % sparse logical array
% adjacency_matrix = sparse( arc_connectivity( :, 1 ), arc_connectivity( :, 2 ), ones( numel( point_uniques ), 1, 'logical' ));
% 
% unique2num_neighbors = sum( adjacency_matrix, 1 )  ...
%                      + sum( adjacency_matrix, 2 )' ;
% 
% vertex_uniques = find( unique2num_neighbors == 1 | unique2num_neighbors >= 3 );
% 
% adjacency_matrix_sans_vertices = adjacency_matrix ; 
% 
% adjacency_matrix_sans_vertices( vertex_uniques, : ) = 0 ;
% adjacency_matrix_sans_vertices( :, vertex_uniques ) = 0 ;
% 
% graph_of_strand_interiors = graph( adjacency_matrix_sans_vertices );
% 
% point_uniques_in_strands = conncomp( graph_of_strand_interiors, 'Type', 'weak', 'OutputForm', 'cells' );

[ bifurcation_point_uniques, point_uniques_in_strands,                                               ...
           arc_indices_in_strands, end_point_uniques_in_strands ]                                       ...
                                                           = get_network_V190( arc_connectivity );

[ point_uniques_in_strands, arc_indices_in_strands, arc_backwards_in_strands ]      ...
                                        = sort_network_V180( arc_connectivity, end_point_uniques_in_strands, ...
                                                             arc_indices_in_strands                     );

end

