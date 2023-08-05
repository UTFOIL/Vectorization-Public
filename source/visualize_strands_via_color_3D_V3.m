function visualize_strands_via_color_3D_V3( strand_subscripts, vessel_directions,                 ...
                                            microns_per_voxel, lumen_radius_in_microns_range,     ...
                                            resolution_factor, vector_directory, ROI_name, network_handle, color_code, varargin )
%% SAM 12/12/17 
% adapted from the function with the same name in the folder AA
%
% V2, in which subpixel scales and locations are passed and the ellipsoids are
% sized and placed appropriately in the image.
%
% V02 in which the inputs are changed to be ammenable to the vectorize_V02 script 3/5/18 SAM
%
% V10 in which the negative laplacian is used (instead of just the value 1) as the fill value for
% each vertex sphere SAM 4/4/18
%
% function name changed from visualize_vertices to visualize_edges
% 
% V11 which is exactly the same function as visualize_vertices_V10
%
% V12 in which the total mean_edge_energies are passed instead of the edge energies at every
% position along the trajectory.  Trajectories are sorted by their total_mean_edge_energies from
% worst to best, and the better trajectories overwrite the worse trajectories in the resulting
% image, essentially plotting the maximum of the total_mean_edge_energies at every point in the
% space that belongs to the volume of a trajectory.  Also the positions is now input as a cell array
% with entries corresponding to trajectories.  SAM 5/5/18
%
% V160 in which the max_edge_energies are passed instead of the "total mean_edge_energies" (inspired
% by thinking about the trajectories as reaction pathways and the max being like the activation
% energy).  Also combing the centerlines into this same function. % SAM 5/14/18
%
% V161 in which the position elements are constructed once into a look up table.  SAM 5/14/18
%
% function name changed from visualize_edges to visualize_depth_via_color in which the trajectories
% are plotted in the xy plane only and the z coordinate is encoded as the color.  The brightness
% encodes the contrast.  The limits of the z coordinate and the contrast values are inputs to this
% function. The x and y limits are also inputs. The original image is also passed so that the model
% can overlay it.  The histogram limits of the original intensity image then also need to be passed.
% The sizing is also demonstrated by overlaying the spheres image, except doing an erosion filter on
% it and differencing so that we only show a mask of that image corresponding to the vessel
% boundaries.  SAM 5/18/18
%
% V2 where the sizing is left out SAM 5/21/18
%
% visualize_depth_via_color_V2( edge_subscripts, max_edge_energies, pixels_per_sigma_range,     ...
% sigma_per_size, directories, original_data_handle, [ 1, 512 ], [ 1, 512 ], [ 30, 60 ], [ -5000, -1000 ], [ - 3000, 3000 ])
%
% SAM 5/21/18
%
% V3 where the objects included in the crop is expanded to include not just those whose centers make
% it into the crop but whose volumes make it into the crop.  Also the color scheme is changed from a
% divergent one to a sequential one (which makes it impossible to also show the contrast with the
% color intensity, so now there is just a lower contrast limit) SAM 5/28/18
%
% vectors from: '180521-220607' or '180528-123246' (these are the same but the earlier one was from
% before the strands phase was created, although that won't affect this figure)
%
% visualize_depth_via_color_V3( edge_subscripts, max_edge_energies, pixels_per_sigma_range,     ...
% sigma_per_size, directories, original_data_handle, [ 1, 512 ], [ 1, 512 ], [ 50, 80 ], -3500, [ - 4500, 3000 ])
%
% 5/28/18
%
% function name changed from visualize_depth_via_color to visualize_strands_via_color SAM 6/12/18
%
% the new function assigns a random color to each strand to visualize the continuity and
% connectivity of the network. SAM 6/12/18
%
% function name changed from visualize_strands_via_color to visualize_strands_via_color_3D 
%
% In this function, we make two 3D images from the vectorization at the ROI specified at the
% resolution specified:  One image is the regular volume-filling image where the volumes of the
% vessels are filled in with the contrast value for each edge.  The other image has the strand
% assignment index as the volume-filling value.  We will look for isosurfaces at the contrast
% threshold specified in the first image and use the second image to determine the color at each
% surface. SAM 6/14/18
%
% V2 in which we interpolate the size index and spatial coordinates along edge index in order to
% create more objects and avoid unfilled spaces between vectors after resolution enhancement. SAM
% 7/18/18
% 
% 11/14/22
%
% V3 in which the vectors are directly converted to an stl style format (no voxelized rendering) % SAM 11/14/22

dimensionality = 3 ;

if ~ isempty( varargin ), strand_r = varargin{ 1 }; if ~ isempty( strand_r ), dimensionality = size( strand_r{ 1 }, 2 ); ...
                                                  else,           strand_r = vessel_directions;                          end; end

%% flip z dimension and swap x and y
% z_limits = - z_limits([ 2, 1 ]);
% 
% v_limits = y_limits ;
% y_limits = x_limits ;
% x_limits = v_limits ;

coefs = [ 1, 1, -1 ] ;
dims = 1 : dimensionality ;
dim_order = [ 2, 1, 3 ];

% lumen_radius_in_pixels_range = lumen_radius_in_pixels_range( :, [ 2, 1, 3    ]);
microns_per_voxel =                           microns_per_voxel( dim_order );

vessel_directions = cellfun( @( x )   coefs          .* x( :,   dim_order        ), vessel_directions, 'UniformOutput', false );
strand_subscripts = cellfun( @( x ) [ coefs,   1   ] .* x( :, [ dim_order, 4 ]   ), strand_subscripts, 'UniformOutput', false );
strand_r          = cellfun( @( x )   coefs( dims )  .* x( :,   dim_order( dims )),          strand_r, 'UniformOutput', false );

is_box_in_center = false ;

% % erasing strands that are above the upper energy limit
% logical_strands_below_upper_energy_limit = max_strand_energies < contrast_limit;
% 
% max_strand_energies       =       max_strand_energies( logical_strands_below_upper_energy_limit );
% strand_subscripts         =         strand_subscripts( logical_strands_below_upper_energy_limit );
% 
% % erasing strands that don't lie in the crop at least partially
% max_subscripts = cell2mat( cellfun( @( x ) max( x, [ ], 1 ),  strand_subscripts, 'UniformOutput', false ));
% min_subscripts = cell2mat( cellfun( @( x ) min( x, [ ], 1 ),  strand_subscripts, 'UniformOutput', false ));
% 
% radii_in_pixels = lumen_radius_in_pixels_range( round( max_subscripts( :, 4 )), : );
% 
% logical_strands_in_crop = max_subscripts( :, 1 ) + radii_in_pixels( :, 1 ) >= y_limits( 1 ) ...
%                         & max_subscripts( :, 2 ) + radii_in_pixels( :, 2 ) >= x_limits( 1 ) ...
%                         & max_subscripts( :, 3 ) + radii_in_pixels( :, 3 ) >= z_limits( 1 ) ...
%                         & min_subscripts( :, 1 ) - radii_in_pixels( :, 1 ) <= y_limits( 2 ) ...
%                         & min_subscripts( :, 2 ) - radii_in_pixels( :, 2 ) <= x_limits( 2 ) ...
%                         & min_subscripts( :, 3 ) - radii_in_pixels( :, 3 ) <= z_limits( 2 ) ;
%                       
% max_strand_energies       =       max_strand_energies( logical_strands_in_crop );
% strand_subscripts         =         strand_subscripts( logical_strands_in_crop );
% 
% % sorting trajectories by max energy in descending order (most negative at end)
% [ max_strand_energies, indices_sorted_by_mean ] = sort( max_strand_energies, 'descend' );
% 
% strand_subscripts         =         strand_subscripts( indices_sorted_by_mean );
% vessel_directions         =         vessel_directions( indices_sorted_by_mean );
% 
% % histogram( max_strand_energies );
% 
% % set the max_strand_energies less than the contrast limits to zero and those above to 255
% % max_strand_energies_binary = max_strand_energies ;
% % 
% % max_strand_energies_binary( max_strand_energies >  contrast_limit ) = 0   ;
% % max_strand_energies_binary( max_strand_energies <= contrast_limit ) = 255 ;

% max_strand_energies( max_strand_energies >=  contrast_limit ) = 0   ;
% max_strand_energies( max_strand_energies <   contrast_limit ) = 1   ;

% interpolate vectors to be isotropic consistent with the dimension with the best resolution and
% then scale by the resolution_factor
% resolution_factors = resolution_factor     * max( lumen_radius_in_pixels_range( 1, : )) ./ lumen_radius_in_pixels_range( 1, : );
% microns_per_voxel  = microns_per_pixel_min * max( lumen_radius_in_pixels_range( 1, : )) ./ lumen_radius_in_pixels_range( 1, : ); % SAM 4/23/21

microns_per_pixel_min             = min( microns_per_voxel );

resolution_factors = resolution_factor * microns_per_voxel ...
                                     ./ microns_per_pixel_min ;

% double the number of entries in the size look up table ( resolution_factor( 4 ) == 2 )
resolution_factors = [ resolution_factors, 2 ];

lumen_radius_in_pixels_range = lumen_radius_in_microns_range ./ microns_per_voxel ;

normed_space = 'L^2' ;

[ ~, ~, ~,        vessel_directions ] = resample_vectors( lumen_radius_in_pixels_range, resolution_factors, strand_subscripts, 0, vessel_directions, normed_space );
[ ~, ~, strand_subscripts, strand_r ] = resample_vectors( lumen_radius_in_pixels_range, resolution_factors, strand_subscripts, 0,          strand_r, normed_space );


microns_per_voxel = microns_per_voxel ./ resolution_factors( 1 : 3 );

number_of_scales = size( lumen_radius_in_microns_range, 1 );

increased_number_of_scales = ceil( resolution_factors( 4 ) * ( number_of_scales - 1 )) / resolution_factors( 4 );

scale_sampling_range = 0 : 1 / resolution_factors( 4 ) : increased_number_of_scales ;

lumen_radius_in_microns_range = exp( interp1(( 0 : number_of_scales - 1 ),           ...
                                            log( lumen_radius_in_microns_range ),    ...
                                            scale_sampling_range                 ))' ;

% % is_decomposing_strands = true
% % 
% % if is_decomposing_strands
% 
%     strand_subscripts_mat = cell2mat( strand_subscripts );
%     strand_subscripts = mat2cell(     strand_subscripts_mat, ...
%                           ones( size( strand_subscripts_mat, 1 ), 1 ), 4 );
% 
%     vessel_directions_mat = cell2mat( vessel_directions );
%     vessel_directions = mat2cell(     vessel_directions_mat, ...
%                           ones( size( vessel_directions_mat, 1 ), 1 ), 3 );
% 
% %     max_strand_energies = strand_subscripts_mat( :, 1 ) * 0 + 1 ; % set all to +1
% 
% % end

% strand_positions = cellfun( @( x ) [                 microns_per_voxel .*  x( :, [ 1, 2, 3 ]),    ...
%                                      lumen_radius_in_microns_range( round( x( :,      4     )))], ...
%                             strand_subscripts,                             'UniformOutput', false );


base  = lumen_radius_in_microns_range( 1 );

ratio = lumen_radius_in_microns_range( 2 ) ...
      / lumen_radius_in_microns_range( 1 );

strand_positions = cellfun( @( x ) [ microns_per_voxel .*  x( :, [ 1, 2, 3 ]),    ...
                                         base * ratio .^ ( x( :,      4     ) - 1 )], ...
                            strand_subscripts,                              'UniformOutput', false    );

% strand_subscripts = cellfun( @round, strand_subscripts, 'UniformOutput', false );

% strand_unit_directions_abs = abs( cell2mat( cellfun( @( x )    sum( x, 1 ) ...
%                                                             / size( x, 1 ), vessel_directions, 'UniformOutput', false )));

% loop through the strands to build the contrast and strand indexed strands images
number_of_strands     = numel( strand_positions );

strand_index_range = 1 : number_of_strands ;

vertices = cell( number_of_strands, 1 );
 normals = cell( number_of_strands, 1 );
  axials = cell( number_of_strands, 1 );
  indics = cell( number_of_strands, 1 );
   radii = cell( number_of_strands, 1 );
   faces = cell( number_of_strands, 1 );
r_vector = cell( number_of_strands, 1 );

for strand_index = strand_index_range
    
    positions_at_strand = strand_positions{ strand_index };

    radius_at_strand = mean( positions_at_strand( :, 4 ));
%     radius_at_strand = max( positions_at_strand( :, 4 ));

    number_of_vertices_per_circle = max( ceil( 2 * pi * radius_at_strand / ( 2 * 3 ^ 0.5 / 3 * microns_per_voxel( 1 ))), 3 ); % target: equilateral triangle patches (will deviate along bends (more samples on insides of turns)
                
%     angle_of_circle_sampling = 2 * pi / number_of_vertices_per_circle ; % radians
% 
%     angles_of_circle_sampling = ( 0 : angle_of_circle_sampling : 2 * pi - angle_of_circle_sampling / 2 )' ;

    angles_of_circle_sampling = 2 * pi * ( 0 : number_of_vertices_per_circle - 1 )' / number_of_vertices_per_circle ;

    number_of_strand_positions = size( positions_at_strand, 1 );

    number_of_vertices_at_strand =       number_of_strand_positions       * number_of_vertices_per_circle ;
    number_of____faces_at_strand = 2 * ( number_of_strand_positions - 1 ) * number_of_vertices_per_circle ;

    vertices{ strand_index } = zeros( number_of_vertices_at_strand, 3 ); %   y, x, z  microns
     normals{ strand_index } = zeros( number_of_vertices_at_strand, 3 ); %  Ny,Nx,Nz   unit-N
      axials{ strand_index } = zeros( number_of_vertices_at_strand, 3 ); %  Ty,Tx,Tz   unit-T 
      indics{ strand_index } = zeros( number_of_vertices_at_strand, 1 ); % custom_colormap( indics, : ) is [ G, R, B ] colors % 
       radii{ strand_index } = zeros( number_of_vertices_at_strand, 1 ); %    radius  microns
       faces{ strand_index } = zeros( number_of____faces_at_strand, 3 ); % v1, v2, v3    idxs
    r_vector{ strand_index } = zeros( number_of_vertices_at_strand, dimensionality ); %   y, x, z  microns

    vertex_idx_range         = ( 0 :     number_of_vertices_per_circle - 1 )' ;
      face_idx_range         = ( 0 : 2 * number_of_vertices_per_circle - 1 )' ;

    vertex_idx_range_opp     =      vertex_idx_range +    number_of_vertices_per_circle  ;
    vertex_idx_range_aft     = mod( vertex_idx_range + 1, number_of_vertices_per_circle );
    vertex_idx_range_opp_bef = mod( vertex_idx_range - 1, number_of_vertices_per_circle )...
                             +                            number_of_vertices_per_circle  ;

    master_direction = rand( 1, 3 ) - 0.5 ; % choose random direction that is orthogonal to the vessel direction (?!?!?! this is not uniform in angle !?!?!) 

%     circle_direction_2 =      circle_direction_2 ...
%                        / dot( circle_direction_2, ...
%                               circle_direction_2  ) .^ 0.5 ;

    for strand_position_index = 1 : number_of_strand_positions
    
        % project the old master direction onto the plane perpendicular to the new vessel direction at this strand_position
%         strand_unit_directions_abs( strand_index )
        direction   = vessel_directions{ strand_index }( strand_position_index, 1 : 3 );
        direction_r =          strand_r{ strand_index }( strand_position_index,   :   );

        circle_direction_1 =   cross( master_direction  ,  direction); % a direction on the plane perpindicular to the vessel and the master_direction
        circle_direction_2 = - cross( circle_direction_1,  direction); % (unnormalized) projection of the master direction onto the plane perpindicular to the vessel 

        circle_direction_2 =      circle_direction_2 ...
                           / dot( circle_direction_2, ...
                                  circle_direction_2  ) .^ 0.5 ; % normalizing to unit length

        master_direction  = circle_direction_2 ;

%         master_direction = ... dot( master_direction, circle_direction_1 ) * circle_direction_1 ... % this dot product is zero
%                          + dot( master_direction, circle_direction_2 ) * circle_direction_2 ; % project the master direction onto the plane perpindicular to vessel
% 
%         master_direction =      master_direction  ...
%                          / dot( master_direction, ...
%                                 master_direction  ) .^ 0.5 ;

%         % rotate the master direction around the vessel direction by an angle equal to - 1/2 the
%         % equilateral triangle length (to make equilateral triangles)
%         master_direction  =                             master_direction   *       cos( angles_of_circle_sampling( 2 ) / 2 ) ...
%                           +           cross( direction, master_direction)  *       sin( angles_of_circle_sampling( 2 ) / 2 ) ...
%                           + direction * dot( direction, master_direction)  * ( 1 - cos( angles_of_circle_sampling( 2 ) / 2 )); % Rodrigues' rotation Formula

%         % rotate the master direction around the vessel direction by an angle equal to the
%         % equilateral triangle length opposite for odd positions (to make zigzag pattern (instead of spiral))
%         if rem( strand_position_index, 2 ) == 1
% 
%             master_direction  =                             master_direction   *       cos( angles_of_circle_sampling( 2 )) ...
%                               -           cross( direction, master_direction)  *       sin( angles_of_circle_sampling( 2 )) ...
%                               + direction * dot( direction, master_direction)  * ( 1 - cos( angles_of_circle_sampling( 2 ))); % Rodrigues' rotation Formula
% 
%         end

        % copy and rotate the master direction around the vessel direction to sample the circle directions
        circle_directions =                             master_direction  .*       cos( angles_of_circle_sampling          ) ...
                          +           cross( direction, master_direction) .*       sin( angles_of_circle_sampling          ) ...
                          + direction * dot( direction, master_direction) .* ( 1 - cos( angles_of_circle_sampling          )); % Rodrigues' rotation Formula

        vertex_ids =     number_of_vertices_per_circle * ( strand_position_index - 1) + 1 + vertex_idx_range ;
          face_ids = 2 * number_of_vertices_per_circle * ( strand_position_index - 1) + 1 +   face_idx_range ;

          face_ids_1 = face_ids(           1 : end / 2 );
          face_ids_2 = face_ids( end / 2 + 1 : end     );

        vertex_id_sum  = number_of_vertices_per_circle * ( strand_position_index - 1) + 1 ;

        vertices{ strand_index }( vertex_ids, : ) =                     positions_at_strand( strand_position_index, 1 : 3 ) ...
                                                  + circle_directions * positions_at_strand( strand_position_index,   4   ) ;
         normals{ strand_index }( vertex_ids, : ) = circle_directions                                                       ;
          axials{ strand_index }( vertex_ids, : ) =        direction  ...
                                     .* ones( size( circle_directions ))                                                    ;
        r_vector{ strand_index }( vertex_ids, : ) =        direction_r  ...
                                     .* ones( size( circle_directions, 1 ), dimensionality )                                                    ;
          indics{ strand_index }( vertex_ids, 1 ) =        vertex_ids                                                       ; 
           radii{ strand_index }( vertex_ids, 1 ) =                     positions_at_strand( strand_position_index,   4   ) ;
         
         

        if strand_position_index < number_of_strand_positions

           faces{ strand_index }( face_ids_1, : )        ...
                      =   vertex_id_sum            ...
                      + [ vertex_idx_range,        ...
                          vertex_idx_range_aft,    ...
                          vertex_idx_range_opp     ];

           faces{ strand_index }( face_ids_2, : )        ...
                      =   vertex_id_sum            ...
                      + [ vertex_idx_range_opp,    ...
                          vertex_idx_range_opp_bef, ...
                          vertex_idx_range          ];
           
        end
%         faces_at_strand( ) = ;


%         % find the linear index of the center pixel of the sphere that defines the strand at this
%         % position
%         strand_position_linear_index = sub2ind( size_of_outer_crop,                           ...
%                                               positions_at_strand( strand_position_index, 1 ), ...
%                                               positions_at_strand( strand_position_index, 2 ), ...
%                                               positions_at_strand( strand_position_index, 3 )  );        
%         
%         % label the spheres and centerlines in an overwriting fashion so that only the lowest energy
%         % strands will shine through in a multiply labeled region.  (Remember that we sorted by max
%         % energy attained before this for loop so the later strands to be written are lower in
%         % energy).  
%         
%                strands_image( max( min(   strand_position_linear_index                                                             ...
%                                       + structuring_element_linear_indexing{ positions_at_strand( strand_position_index, 4 )},    ...
%                                       number_of_image_voxels ),                                                                ...
%                                  1                                                                                          )) ...
%                                                                                        =       max_strand_energies( strand_index );        
%         
%         strand_index_image( max( min(   strand_position_linear_index                                                             ...
%                                       +    coloring_element_linear_indexing{ positions_at_strand( strand_position_index, 4 )},    ...
%                                       number_of_image_voxels ),                                                                ...
%                                  1                                                                                          )) ...
%                                                                                        =  max_strand_energies( strand_index )                                                                      ...
%                                                                                        * strand_index ;
%         
%         % note consider replacing the strand index image with a contrast image so that the colormap
%         % depends on the contrast
        
    end % for position
end % for trajectory

temp_cell = cellfun( @(a) a(:,1), vertices, 'UniformOutput',false );

vertex_s = num2cell([ 0; cumsum( cellfun( @(     x        )               ...
                                       size(     x   , 1  ), vertices( 1 ...
                                                               : end - 1 ...
                                                                     )))]);
vertices =     cell2mat(                                     vertices   ) ;
r_vector =     cell2mat(                                     r_vector   ) ;
 normals =     cell2mat(                                      normals   ) ;
  axials =     cell2mat(                                       axials   ) ;
   radii =     cell2mat(                                       radii   ) ;
   faces =     cell2mat(         cellfun( @(   v,  s      )             ...
                                               v + s       ,    faces,  ...
                                                             vertex_s,  ...
                                               'UniformOutput', false  ) ) ;
switch color_code 
    case 'strands', ones_cell =     cellfun(@(a) ones(size(a)), temp_cell, 'UniformOutput',false ); 
                   index_cell = mat2cell(        ( 1:length(    temp_cell ))', ...
                                                 ones(size(     temp_cell ))                     );
  indics = cell2mat(cellfun(@(a,b) a*b,ones_cell,index_cell, 'UniformOutput',false ));
    otherwise
  indics =                       ( 1 : size( vertices, 1  ))'; 
end
%   colors =     cell2mat(         cellfun( @(   v,  s      )             ...
%                                                v + s       ,   colors,  ...
%                                                              vertex_s,  ...
%                                                'UniformOutput', false  ) ) ;   

if is_box_in_center

    % y_limits_box = round([ 512 / 3, 2 * 512 / 3 ]);
    % x_limits_box = round([ 512 / 3, 2 * 512 / 3 ]);
    % z_limits_box = round([ 512 / 3, 2 * 512 / 3 ]);
    % 
    % y_limits_box = y_limits_box + outer_crop_margin ;
    % x_limits_box = x_limits_box + outer_crop_margin ;
    % z_limits_box = z_limits_box + outer_crop_margin ;
    % 
    % center_box_energy = min(max_strand_energies);
    % 
    % indexing issues !!!!!!!
    % 
    % strand_index_image( y_limits_box(1) : y_limits_box(2), x_limits_box(1) : x_limits_box(2), z_limits_box(1)                  ) = center_box_energy ;
    % strand_index_image( y_limits_box(1) : y_limits_box(2), x_limits_box(1) : x_limits_box(2), z_limits_box(2)                  ) = center_box_energy ;
    % strand_index_image( y_limits_box(1) : y_limits_box(2), x_limits_box(1)                  , z_limits_box(1) : z_limits_box(2)) = center_box_energy ;
    % strand_index_image( y_limits_box(1) : y_limits_box(2), x_limits_box(2)                  , z_limits_box(1) : z_limits_box(2)) = center_box_energy ;
    % strand_index_image( y_limits_box(1)                  , x_limits_box(1) : x_limits_box(2), z_limits_box(1) : z_limits_box(2)) = center_box_energy ;
    % strand_index_image( y_limits_box(1)                  , x_limits_box(1) : x_limits_box(2), z_limits_box(1) : z_limits_box(2)) = center_box_energy ;

    % Define the vertexes of the unit cubic

    ver = 0.5     ...
        - [ 1 1 0  ;
            0 1 0  ;
            0 1 1  ;
            1 1 1  ;
            0 0 1  ;
            1 0 1  ;
            1 0 0  ;
            0 0 0 ];

    num_vers = size( ver, 1 );

    % %  Define the faces of the unit cubic
    % fac = [ 1 2 3 4  ;
    %         4 3 5 6  ;
    %         6 7 8 5  ;
    %         1 2 8 7  ;
    %         6 7 1 4  ;
    %         2 3 5 8 ];

    fac = [ 1 2 3  ;
            1 4 3  ;
            4 3 5  ;
            4 6 5  ;
            6 7 8  ;
            6 5 8  ;
            1 2 8  ;
            1 7 8  ;
            6 7 1  ; 
            6 4 1  ;
            2 3 5  ;
            2 8 5 ];

    num_fac = size( fac, 1 );

    origin = ( size_of_inner_crop - 1 ) / 2 + 1 ;
    width  = ( size_of_inner_crop - 1 ) / 3     ;

    % ver = [ ver(:,1)*width(1)+origin(1), ...
    %         ver(:,2)*width(2)+origin(2), ...
    %         ver(:,3)*width(3)+origin(3)  ];

    ver = ver .* width + origin ;

    box_index_as_strand = 1 ;

    color = box_index_as_strand * ones( num_fac, 1 );

    % p_box = patch('Faces',fac,'Vertices',ver,'FaceColor',0.1*[1,1,1],'FaceAlpha',0.25);

end

% extract the colored isosurface from the crop and render it in 3D with lighting

% % flip the z dimension to be consistent with the normal orientation of the brain
% strands_image      = flip(      strands_image, 3 );
% strand_index_image = flip( strand_index_image, 3 );
% 
% % flip the y dimension to be consistent with the other z projected special outputs of vectorize()
% strands_image      = flip(      strands_image, 1 );
% strand_index_image = flip( strand_index_image, 1 );

% strands_image      = flip(      strands_image, 2 );
% strand_index_image = flip( strand_index_image, 2 );

number_of_strands     = numel( strand_positions );

if is_box_in_center

    % !!!!!! this code incompatible with non default color_code options, modularize !!!!!  % SAM 11/12/22
    custom_colormap = 0.1 + 0.8 * rand( number_of_strands + 1, 3 ); % default matlab colormaps have 64 values

    % % set the first entry to gray and reserve it for the box
    custom_colormap( 1, : ) = 0.9 ;

    strand_index_image = strand_index_image + 1 ;

else
    
    custom_alphamap = ones( size( axials( :, 1 )));

    switch color_code
    
        case 'strands'
            % !!!!!!!! 0.1 -> 0.2 % SAM 11/12/22
            % !!!!!!!!!!!!!!!!! num strands -> num total positions
            rng('default') % do this for consistency in random strand color between the 2D visualization FXN as well

            custom_colormap = 0.1 + 0.8 * rand( number_of_strands, 3 ); % default matlab colormaps have 64 values

        case {   'directions', ...
               'z-directions'  }
 
% %             custom_colormap = cell2mat( vessel_directions );
%             custom_colormap = abs( axials( :, dim_order ));
            custom_colormap = abs( axials );

            switch color_code
%                 case    'directions' % default case, do notjhing
                case   'z-directions'

%                     custom_colormap( :, [ 1, 2 ]) = sum( custom_colormap( :, [ 1, 2 ]) .^ 2, 2 ) .^ 0.5 ;
                    custom_colormap( :, 1 ) = sum( custom_colormap( :, [ 1, 2 ]) .^ 2, 2 ) .^ 0.5 ;
                    custom_colormap( :, 2 ) ...
                  = custom_colormap( :, 1 ); 
            end

        case { 'r-directions', ...
              'xy-directions' }
                
%                     custom_colormap( :, [ 1, 2 ]) = sum( custom_colormap( :, [ 1, 2 ]) .^ 2, 2 ) .^ 0.5 ;
%                     strand_r

            dim_order = dim_order( dims ); % ??? is dim_order effecting the following operation at all >???? SAM 3/30/23

%                 % sqt(y2 + x2), 1, proj( r, sqt(y2 + x2))  -> R, G, B
            strand_r_component = abs( dot(   axials( :, dim_order ), r_vector( :, dim_order ), 2 )) ...
                              ./      dot( r_vector( :, dim_order ), r_vector( :, dim_order ), 2 ) .^ 0.5 ;

            custom_colormap  = [ 0,  0 , 1 ] .*       strand_r_component ...
                             + [ 0, 1/2, 0 ] ...
                             + [ 1,  0 , 0 ] .* ( 1 - strand_r_component );
            
        case 'depth' % !!!!!!!!!!!! wishlist

%             custom_colormap = depth ;

        case 'radii' % !!!!!!!!!!!! wishlist
            
            max_radius  = 15 ; % mircons
            min_radius  = 3  ; % microns

%             custom_colormap = max( min(   ( log(   radii   ) .* ones( size( axials )) - log( min_radius )) ...
%                                         / ( log( max_radius )                         - log( min_radius )), 1 ), 0 );
            radius_coordinate = max( min(   ( log(   radii   )   - log( min_radius )) ...
                                          / ( log( max_radius )  - log( min_radius )), 1 ), 0 );

            custom_colormap = [ 0.25, 0, .75 ] .* ( 1 - radius_coordinate ) ...
                            + [ 0.75, 1,   0 ] .*       radius_coordinate   ;
% %             custom_alphamap =   1        - 0.5  *       radius_coordinate   ;
% %             custom_alphamap =                           radius_coordinate   ;
%             custom_alphamap =                     1 -   radius_coordinate   ;
% %             custom_alphamap =           1 * ones( size( radius_coordinate ));
% %             custom_alphamap =           0.5 * ones( size( radius_coordinate ));

%             alphamap( custom_alphamap );

        otherwise

            custom_colormap = cell2mat( varargin{ 1 }) ;

            color_code_name = color_code ;

    end
end

switch color_code

    case       'strands', color_code_name =       'Strands' ;
    case         'radii', color_code_name =        'Radius' ;
    case         'depth', color_code_name =         'Depth' ;
    case    'directions', color_code_name =    'Directions' ;
    case  'z-directions', color_code_name =  'Z-Directions' ;
    case  'r-directions', color_code_name =  'R-Directions' ;
    case 'xy-directions', color_code_name = 'XY-Directions' ;

end

is_overlaying = contains( ROI_name, 'cropped_out' ) ;
% is_overlaying = false ;

if ~ is_overlaying

    fh = figure ;

else

%     fh = gcf ;

    % find last fingure (by figure number) that matches the name of the currrent figure for overlay plotting
    fh = findobj( 'Type', 'Figure' );

    fh_numbers = arrayfun( @( x ) x.Number, fh );

    [ ~, sort_indices ] = sort( fh_numbers );

    last_match_index = find( arrayfun( @( x )    contains( x.Name, [ ' 3D ', color_code_name ])               ...
                                              && contains( x.Name,           network_handle   ), fh( sort_indices )), 1, 'last' );

    fh = fh( sort_indices( last_match_index ));

end

visualization_name = [ network_handle, '_', ROI_name, ' 3D ', color_code_name ];

set( fh, 'Name', visualization_name )

% colormap( custom_colormap )

% strands_image = strands_image( 1 + outer_crop_margin : size_of_outer_crop( 1 ) - outer_crop_margin, ...
%                                1 + outer_crop_margin : size_of_outer_crop( 2 ) - outer_crop_margin, ...
%                                1 + outer_crop_margin : size_of_outer_crop( 3 ) - outer_crop_margin  );
                       
% strand_index_image ...
%      = strand_index_image( 1 + outer_crop_margin : size_of_outer_crop( 1 ) - outer_crop_margin, ...
%                            1 + outer_crop_margin : size_of_outer_crop( 2 ) - outer_crop_margin, ...
%                            1 + outer_crop_margin : size_of_outer_crop( 3 ) - outer_crop_margin  );

% strand_index_image_nnz_idxs                = find( strand_index_image );

% color_image = zeros([ 3,                     size( strand_index_image          )]);
% color_image(          :,                           strand_index_image_nnz_idxs ) ...
%     = custom_colormap(         strand_index_image( strand_index_image_nnz_idxs ), : )';

% % begin citation: https://www.mathworks.com/help/matlab/ref/isosurface.html SAM 6/14/18
% [ F, V ] = isosurface( strands_image, 0.5 );
% 
% % convert color coding from vertices (for interp) to faces (for 'flat' FaceColor)
% 
% % do nearest neighbor interpolation (instead of whatever isosurface does without option)
% C = interp3( strand_index_image, V( :, 1 ), ...
%                                  V( :, 2 ), ...
%                                  V( :, 3 ), 'nearest' );
% % !!!!!!!!!!!!!!!! put this interp3 in a for loop across x, y, z and replace the index image with
% % the color directly (and then go from color to index into the color map if needed by nearest neighbor
% % interpolation of the 3D map index function over R, G, B color space), need color map to uniformly sample the
% % surface of a 1/8 color sphere with R,G, and B at the x,y, and z poles, respectively
% 
% R = interp3( squeeze( color_image( 1, :, :, : )),    V( :, 1 ), ...
%                                                      V( :, 2 ), ...
%                                                      V( :, 3 ) );
% G = interp3( squeeze( color_image( 2, :, :, : )),    V( :, 1 ), ...
%                                                      V( :, 2 ), ...
%                                                      V( :, 3 )  );
% B = interp3( squeeze( color_image( 3, :, :, : )),    V( :, 1 ), ...
%                                                      V( :, 2 ), ...
%                                                      V( :, 3 )  );
% 
% 
% % C = C( F( :, 1 ));    
% C = median( C( F ), 2 );    

if is_box_in_center
    
    % adjust vertex indices in the faces array to accomodate the added vertices at the start
    F = F + num_vers ;
    
end




if is_box_in_center

    % p = patch('Vertices',[ ver; V ],'Faces',[ fac; F ],'FaceVertexCData', [ color; C ],...
    %           'FaceColor','flat','FaceAlpha', 'flat', ...
    %           'AlphaData', [ 0.1 * ones( size( color )); ones( size( C ))], ...
    %           'AlphaDataMapping', 'none', 'EdgeColor','none');


    %     p = patch('Vertices',[ ver; V ],'Faces',[ fac; F ],'FaceVertexCData', [ color; C ],...
    %               'FaceColor','flat', 'EdgeColor','none');
    colormap( custom_colormap )
    % make objects transparent
    p = patch('Vertices',[ ver; V ],'Faces',[ fac; F ],'FaceVertexCData', [ indics; C ],...
              'FaceColor', 'flat',        'EdgeColor', 'none', ...
              'FaceAlpha', 'flat', 'AlphaDataMapping', 'none', ...
              'FaceVertexAlphaData', [ 0.25 * ones( size( indics )); 1 * ones( size( C ))]);
    
              % % unifrom transparency
    %     p = patch('Vertices',[ ver; V ],'Faces',[ fac; F ],'FaceVertexCData', [ color; C ],...
    %               'FaceColor','flat', 'EdgeColor','none', ...
    %               'FaceAlpha', 0.5 );          
else
    
    if is_overlaying

        FaceAlphaValue = 0.33 ;

    else

%         FaceAlphaValue = 'flat' ;
        FaceAlphaValue = 1 ;

    end
%     p = patch('Vertices',V ,'Faces',F ,'FaceVertexCData', C ,...
%               'FaceColor','flat', 'EdgeColor','none');

    p = patch('Vertices'      , vertices( :, [ 2, 1, 3 ]), ...
               'VertexNormals',  normals( :, [ 2, 1, 3 ]), ...
          'FaceVertexCData'  ,   custom_colormap( indics, : ), ...
...          'FaceVertexAlphaData', indics, ...
           'FaceVertexAlphaData', custom_alphamap,             ...
              'AlphaDataMapping', 'none', ...
                   'Faces'      ,   faces, ...
                 'FaceColor'  , 'interp',    ...
                 'FaceAlpha'  , FaceAlphaValue , ...
  ...   'FaceColor'  , [ 0.5, 0.5, 0.5 ], ...
                        'EdgeColor'  ,'none');
%                    'EdgeColor'  ,'black');

end

%% write ply file with surface mesh, color, alpha, material
file_path = [ vector_directory, visualization_name, '.ply'];
disp(['Writing 3D graphics file (Polygon File Format) to ', file_path ])
% Convert to PLY format
ply_header = sprintf( 'ply\nformat ascii 1.0\ncomment Created by the SLAVV V3 software, SAM 1/31/23 (source: github.com/UTFOIL/Vectorization-Public)\nelement vertex %d\nproperty float x\nproperty float y\nproperty float z\nproperty float nx\nproperty float ny\nproperty float nz\nproperty uchar red\nproperty uchar green\nproperty uchar blue\nelement face %d\nproperty list uchar uint vertex_indices\nend_header\n',size(vertices,1),size(faces,1));
ply_data = sprintf('%f %f %f %f %f %f %d %d %d\n',[vertices( :, [ 2, 1, 3 ]), ...
                                                    normals( :, [ 2, 1, 3 ]), round( 255 * custom_colormap( indics, : ))]');
ply_data = [ply_data sprintf('3 %d %d %d\n',faces'-1)];
% Write the PLY file
fid = fopen( file_path, 'w' );
fprintf(fid,'%s',ply_header);
fprintf(fid,'%s',ply_data);
fclose(fid);

%% render in matlab figure

% % isonormals( strands_image, p )
% p.VertexNormals = normals ;

% % convert the vertex coordinates to microns
% p.Vertices      = ( p.Vertices      - 1 ) / resolution_factor * microns_per_pixel_min ;
% % p.VertexNormals =   p.VertexNormals       * resolution_factor / microns_per_pixel_xy ; % transpose of the inverse transformation is done to normals (simplifies to a purely scaling transformation)
% % p.VertexNormals = p.VertexNormals ./ sum( p.VertexNormals .^ 2, 2 ) .^ 0.5 ;

% p.FaceColor = 'red';
% p.EdgeColor = 'none';
daspect([1 1 1])
camproj perspective
campos
% view(0,0); 
view(80,45);
axis tight
camlight right
% camlight left
camlight headlight
material shiny
%   MATERIAL SHINY makes the objects shiny.
%   MATERIAL DULL makes the objects dull.
%   MATERIAL METAL makes the objects metallic.

if is_overlaying, material dull, end
lighting gouraud
% end citation: https://www.mathworks.com/help/matlab/ref/isosurface.html SAM 6/14/18


%% write a file containing camera position and lighting as well

end % function