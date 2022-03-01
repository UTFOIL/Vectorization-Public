function [ point_coordinates, arc_connectivity, arc_diameters ] = strand2casx( strand_subscripts, microns_per_voxel, lumen_radius_in_microns_range )
%% output_to_LPPD, SAM 6/10/19
% This function takes the strand and vertex objects and converts them into the format that the LPPD
% lab in Chicago, IL is expecting.  Format specified in the .casX file format documented in the
% Andreas Linninger/Grant Hartung directory.

% [ ~, bifurcation_vertices_LUT ] = sort( bifurcation_vertices );

% vertex_space_subscripts = double( vertex_space_subscripts );

number_of_strands = length( strand_subscripts );

% [ unique_terminal_vertices, ~, strands2unique_vertices ] = unique( strands2vertices( : ));
% %   [C,IA,IC] = UNIQUE(A) also returns index vectors IA and IC such that
% %   C = A(IA) and A = C(IC) (or A(:) = C(IC), if A is a matrix or array).

strand_endpoint_subscripts = cell2mat([ cellfun( @( x ) x(  1 , 1 : 3 ), strand_subscripts, 'UniformOutput', false ); ...
                                        cellfun( @( x ) x( end, 1 : 3 ), strand_subscripts, 'UniformOutput', false )  ]);

[ unique_terminal_points, ~, strands2unique_vertices ] = unique( strand_endpoint_subscripts, 'rows' );
%   [C,IA,IC] = UNIQUE(A,'rows') also returns index vectors IA and IC such
%   that C = A(IA,:) and A = C(IC,:). 

number_of_unique_terminal_vertices = size( unique_terminal_points, 1 );

number_of_points_in_strands = cellfun( @( x ) size( x, 1 ), strand_subscripts );

%            number_of_arcs_in_strands = number_of_points_in_strands - 1 ;
number_of_interior_points_in_strands = number_of_points_in_strands - 2 ;

total_number_of_points_in_strands = sum( number_of_points_in_strands );

total_number_of_arcs                   = total_number_of_points_in_strands -     number_of_strands ;
total_number_of_strand_interior_points = total_number_of_points_in_strands - 2 * number_of_strands ;

number_of_unique_points = number_of_unique_terminal_vertices + total_number_of_strand_interior_points ;

point_coordinates   = zeros( number_of_unique_points, 3 );
arc_connectivity    = zeros( total_number_of_arcs   , 2 );
arc_radius_indices  = zeros( total_number_of_arcs   , 1 );

% place all of the terminal vertices in the point_coordinates list first, so that we don't add the
% bifurcations in multiple times when looping through the strands
point_coordinates( 1 : number_of_unique_terminal_vertices, : ) = unique_terminal_points ;

% loop through every arc in each strands
strand_index_range = 1 : number_of_strands ;

number_of_arcs         =                                  0 ;
number_of_added_points = number_of_unique_terminal_vertices ;

for strand_index = strand_index_range
    
    % add the first half of the arc from the terminal vertex to the first interior vertex
    number_of_arcs = number_of_arcs + 1 ;

    arc_connectivity( number_of_arcs, 1 ) = strands2unique_vertices( strand_index );
    
    arc_radius_indices( number_of_arcs ) = strand_subscripts{ strand_index }( 1, 4 ) / 2 ; 
    
    for strand_point_index = 1 : number_of_interior_points_in_strands( strand_index )
        
        % add the two halves from the arcs before and after this interior point
        number_of_arcs         = number_of_arcs         + 1 ;
        number_of_added_points = number_of_added_points + 1 ;
        
        point_coordinates( number_of_added_points, : ) = strand_subscripts{ strand_index }( strand_point_index + 1, 1 : 3 );
        
        arc_connectivity( number_of_arcs - 1, 2 ) = number_of_added_points ;
        arc_connectivity( number_of_arcs    , 1 ) = number_of_added_points ;
        
        arc_radius_indices( number_of_arcs - 1 ) = arc_radius_indices( number_of_arcs - 1 ) + strand_subscripts{ strand_index }( strand_point_index + 1, 4 ) / 2 ;
        arc_radius_indices( number_of_arcs     ) =                                            strand_subscripts{ strand_index }( strand_point_index + 2, 4 ) / 2 ;
    
    end % FOR interior_points
    
    % add the last half of the arc from the last interior vertex to the other terminal vertex.
    arc_connectivity( number_of_arcs, 2 ) = strands2unique_vertices( strand_index + number_of_strands );
    
    arc_radius_indices( number_of_arcs ) = arc_radius_indices( number_of_arcs ) + strand_subscripts{ strand_index }( end, 4 ) / 2 ;
    
end % FOR strand

% switch x and y coordinates of the points.  LPPD wants x, y, z coordinates.
point_coordinates = point_coordinates( :, [ 2, 1, 3 ]);

% rotate so z is down
point_coordinates( :, [ 2, 3 ]) = - point_coordinates( :, [ 2, 3 ]);

% convert to microns:
point_coordinates = point_coordinates .* microns_per_voxel ;

arc_radii = exp( interp1( log( lumen_radius_in_microns_range ), arc_radius_indices ));

% LPPD wants diameter not radius
arc_diameters = 2 * arc_radii ;

end % FUNCTION


