function [ vertex_space_subscripts, vertex_scale_subscripts, bifurcation_vertices, strands2vertices, strand_space_subscripts, strand_scale_subscripts, microns_per_voxel, lumen_radius_in_microns_range, size_of_image ] = input_from_LPPD( point_coordinates, arc_connectivity, arc_diameters )
%% input_from_LPPD, SAM 6/25/19
% This function imports strand and vertex objects by converting them from the format that the LPPD
% lab in Chicago, IL specified (found in the .casX file format documented in the Andreas
% Linninger/Grant Hartung directory).

%% check for redundancy in arc_connectivity matrix

% choose the first of multiple arces from A to B
[ arc_connectivity_temp, indices_of_unique_arcs ] = unique( arc_connectivity, 'rows', 'stable' );

if length( arc_connectivity_temp ) < length( arc_connectivity )

    warning( 'Redundant arc(s) detected in the arc_connectivity matrix. Erasing the later instance(s) and its(their) associated arc_diameter(s)' )

end

arc_connectivity = arc_connectivity( indices_of_unique_arcs, : );
arc_diameters    = arc_diameters(    indices_of_unique_arcs    );

% choosing the first of two directions of edges in the cases where we have pairs of mutual
% trajectories (from A to B and B to A).
[ ~, mutual_arc_indices, reverse_mutual_arc_indices ]                  ...
    = intersect([ arc_connectivity( :, 1 ), arc_connectivity( :, 2 )], ...
                [ arc_connectivity( :, 2 ), arc_connectivity( :, 1 )], ...
                'rows', 'stable'                                       );

is_mutual_arc_pair_to_be_erased = mutual_arc_indices > reverse_mutual_arc_indices ;

if ~ isempty( is_mutual_arc_pair_to_be_erased ) 
    
    warning( 'Mutual arc pair(s) detected in the arc_connectivity matrix. Erasing the later instance(s) and its(their) associated arc_diameter(s)' )

end

arc_connectivity( mutual_arc_indices( is_mutual_arc_pair_to_be_erased ), : ) = [ ];
arc_diameters(    mutual_arc_indices( is_mutual_arc_pair_to_be_erased )    ) = [ ];

%% convert arc diameters to point diameters
% casX defines diameters along each arc, whereas the strand format has radii at each position

point_indices = unique( arc_connectivity( : ))';

% check that all point indices are positive integers (MATLAB convention)
if any( point_indices <= 0 ),                      error( 'Some point index(-ices) is(are) non-positive in arc_connectivity input' ), end

if any( round( point_indices ) ~= point_indices ), error( 'Some point index(-ices) is(are) non-integer in arc_connectivity input'  ), end

number_of_points = length( point_indices );

points2radii = sparse( point_indices, ones( number_of_points, 1 ), ones( number_of_points, 1 ));

number_of_arcs = size( arc_connectivity, 1 );

for point_index = point_indices
    
    % find all arcs connected to this point and ...
    point2arc_or_arcs = mod( find( arc_connectivity( : ) == point_index ) - 1, number_of_arcs ) + 1 ;
    
    % ... average their diameters and convert to radius ( / 2 )
    points2radii( point_index ) = geomean( arc_diameters( point2arc_or_arcs )) / 2 ;
    
end % FOR point_indices

%% Decide output space and size resolutions

% decide on a spatial resolution to represent these vectors in
microns_per_voxel = [ 1, 1, 1 ];

radii_min = full( min( points2radii ));

radii_max = full( max( points2radii ));

scales_per_octave = 7 / 2 ;

% decide on a scale resolution for displaying these vectors
number_of_scales = ceil( scales_per_octave * log( radii_max / radii_min ) / log( 2 ));

% build the look up table that holds the display scales
lumen_radius_in_microns_range = logspace( log( radii_min ) / log( 10 ), log( radii_max ) / log( 10 ), number_of_scales );

% interpolate the LUT switching r_i and i, to go from radius to index
points2scale_subscripts_nonzeros = interp1( log( lumen_radius_in_microns_range ), 1 : number_of_scales, log( nonzeros( points2radii )));

points2scale_subscripts = sparse( point_indices, ones( number_of_points, 1 ), points2scale_subscripts_nonzeros );

min_object_coordinates = min( point_coordinates - points2radii );
max_object_coordinates = max( point_coordinates + points2radii );

% adjust coordinates to start at 1 (for subscripts into the image
point_coordinates  = point_coordinates - min_object_coordinates + 1 ;

% create the smallest image that captures all of the objects
size_of_image = max_object_coordinates - min_object_coordinates + 1 ;

% % excerpt from output_to_LPPD
% arc_radii = exp( interp1( log( lumen_radius_in_microns_range ), arc_radius_indices ));

%% find the vertices from the arc_connectivity matrix
% vertices are points that are either bifurcations or endpoints of vessels.

max_point_index = max( arc_connectivity( : ));

adjacency_matrix         = sparse( arc_connectivity( :, 1 ),   ...
                                   arc_connectivity( :, 2 ),   ...
                       ones( size( arc_connectivity( :, 1 ))), ...
                         max_point_index, max_point_index      );
                       
points2number_of_neighbors = sum( adjacency_matrix, 1 ) + sum( adjacency_matrix, 2 )';

% ??? if a point_coordinate isn't referenced by any arc_connectivity pair, does it even it exist ???

% points2number_of_neighbors == 0 shouldn't exist
% points2number_of_neighbors == 1 is an endpoint,
% points2number_of_neighbors == 2 is a strand interior position,
% points2number_of_neighbors  > 2 is a bifurcation

is_point_a_vertex = points2number_of_neighbors > 2 | points2number_of_neighbors == 1 ;

vertices2points = find( is_point_a_vertex );

vertex_space_subscripts =       point_coordinates(       vertices2points, : ) ;
vertex_scale_subscripts = full( points2scale_subscripts( vertices2points    ));

is_vertex_a_bifurcation = points2number_of_neighbors( vertices2points ) > 2 ;

bifurcation_vertices = find( is_vertex_a_bifurcation );

%% convert points to positions in strands 
% A strand connects two vertices with interior strand positions.  This is done by starting at each
% vertex and tracing out from there.  The points are referred to as follows:
% 
% point_A/point_C_0 -> point_B/point_C_1 -> point_C_2 -> ... -> point_C_N
%
% where N is one less than the number of positions in this example strand.

number_of_vertices = length( vertices2points );

points2vertices = sparse( vertices2points, ones( number_of_vertices, 1 ), 1 : number_of_vertices ); 

% points2vertices( vertices2points ) = 1 : number_of_vertices ;

% points will appear in more than one strand if that point is a bifurcation
number_of_repeated_points = sum( points2number_of_neighbors( vertices2points ) - 1 );

number_of_strands = sum( points2number_of_neighbors( vertices2points )) / 2 ;

number_of_positions_per_strand = zeros( number_of_strands, 1 );
strands2vertices               = zeros( number_of_strands, 2 );

total_number_of_positions = number_of_points + number_of_repeated_points ;

% this is a list of all of the points arranged in strands so that a cell array can be formed
positions2points_mat = zeros( total_number_of_positions, 1 );

strand_index = 0 ;

cumulative_number_of_positions = 0 ;

% vertex_indices = 1 : number_of_vertices ;

number_of_vertices = 0 ;

for point_A_index = vertices2points

    number_of_vertices = number_of_vertices + 1 ;
    
    % find all point B's that are connected to this point_A. Each one of these point_B's is the
    % second point of a strand that starts at point_A.
    point_B_indices = find( adjacency_matrix( :, point_A_index )' + adjacency_matrix( point_A_index, : ));
    
    for point_B_index = point_B_indices
    
        % IF point_A connected back to point_A then some point_B will be erased. Skip those
        if adjacency_matrix( point_B_index, point_A_index ) || adjacency_matrix( point_A_index, point_B_index )
        
            strand_index = strand_index + 1 ;

            % record point_A before tracing out the interior points of the strand
    %         strands2vertices( strand_index, 1 ) = number_of_vertices ;
            strands2vertices( strand_index, 1 ) = full( points2vertices( point_A_index ));

            number_of_positions_at_strand = 1 ;

            positions2points_mat( cumulative_number_of_positions + number_of_positions_at_strand ) = point_A_index ;

            % initialize the trace at point_B
            tracing = true ;

            point_C_index = point_B_index ;

            while tracing

                number_of_positions_at_strand = number_of_positions_at_strand + 1 ;

                positions2points_mat( cumulative_number_of_positions + number_of_positions_at_strand ) = point_C_index ;

                tracing = ~ is_point_a_vertex( point_C_index );

                % If this point is a vertex, don't zero it out yet.  Let it have a chance to be a
                % point_A of its own in the outer FOR loop.  If we get to that vertex and it has no
                % neighbors anymore, then there it will have no point_B's thus no strands of its own.

                % IF point_C is not a vertex
                if tracing

                    if point_C_index == point_B_index

                        % we don't want to find the point_A index and go back where we started
                        point_C_index_temp = find( adjacency_matrix( :, point_C_index )' + adjacency_matrix( point_C_index, : ));

                        point_C_index_temp( point_C_index_temp == point_A_index ) = [ ];

                    else

                        point_C_index_temp = find( adjacency_matrix( :, point_C_index )' + adjacency_matrix( point_C_index, : ), 1 );

                    end

                    adjacency_matrix( :, point_C_index ) = 0 ;
                    adjacency_matrix( point_C_index, : ) = 0 ;

                    point_C_index = point_C_index_temp ;

                end % IF tracing            

            end % WHILE tracing

            strands2vertices( strand_index, 2 ) = full( points2vertices( point_C_index ));

            number_of_positions_per_strand( strand_index ) = number_of_positions_at_strand ;

            cumulative_number_of_positions = cumulative_number_of_positions + number_of_positions_at_strand ;

        end % IF point_B exists
        
    end % FOR point_B_index

    % zero out the original vertex because we never want to see a strand connected here again
    adjacency_matrix( :, point_A_index ) = 0 ;
    adjacency_matrix( point_A_index, : ) = 0 ;
    
end % FOR vertices2points

% Do unit conversions

strand_space_subscripts = mat2cell(       point_coordinates(       positions2points_mat, : ), number_of_positions_per_strand, 3 );
strand_scale_subscripts = mat2cell( full( points2scale_subscripts( positions2points_mat    )), number_of_positions_per_strand, 1 );

end % FUNTION