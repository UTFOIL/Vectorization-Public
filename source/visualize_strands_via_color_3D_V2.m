function visualize_strands_via_color_3D_V2( strand_subscripts, max_strand_energies,                 ...
                                            microns_per_pixel_xy, lumen_radius_in_pixels_range,     ...
                                            resolution_factor, y_limits, x_limits, z_limits, contrast_limit )
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

is_box_in_center = false ;

% erasing strands that are above the upper energy limit
logical_strands_below_upper_energy_limit = max_strand_energies < contrast_limit;

max_strand_energies       =       max_strand_energies( logical_strands_below_upper_energy_limit );
strand_subscripts         =         strand_subscripts( logical_strands_below_upper_energy_limit );

% erasing strands that don't lie in the crop at least partially
max_subscripts = cell2mat( cellfun( @max,  strand_subscripts, 'UniformOutput', false ));
min_subscripts = cell2mat( cellfun( @min,  strand_subscripts, 'UniformOutput', false ));

radii_in_pixels = lumen_radius_in_pixels_range( round( max_subscripts( :, 4 )), : );

logical_strands_in_crop = max_subscripts( :, 1 ) + radii_in_pixels( :, 1 ) >= y_limits( 1 ) ...
                        & max_subscripts( :, 2 ) + radii_in_pixels( :, 2 ) >= x_limits( 1 ) ...
                        & max_subscripts( :, 3 ) + radii_in_pixels( :, 3 ) >= z_limits( 1 ) ...
                        & min_subscripts( :, 1 ) - radii_in_pixels( :, 1 ) <= y_limits( 2 ) ...
                        & min_subscripts( :, 2 ) - radii_in_pixels( :, 2 ) <= x_limits( 2 ) ...
                        & min_subscripts( :, 3 ) - radii_in_pixels( :, 3 ) <= z_limits( 2 ) ;
                      
max_strand_energies       =       max_strand_energies( logical_strands_in_crop );
strand_subscripts         =         strand_subscripts( logical_strands_in_crop );

% sorting trajectories by max energy in descending order (most negative at end)
[ max_strand_energies, indices_sorted_by_mean ] = sort( max_strand_energies, 'descend' );

strand_subscripts         =         strand_subscripts( indices_sorted_by_mean );

% histogram( max_strand_energies );

% set the max_strand_energies less than the contrast limits to zero and those above to 255
% max_strand_energies_binary = max_strand_energies ;
% 
% max_strand_energies_binary( max_strand_energies >  contrast_limit ) = 0   ;
% max_strand_energies_binary( max_strand_energies <= contrast_limit ) = 255 ;

max_strand_energies( max_strand_energies >   contrast_limit ) = 0   ;
max_strand_energies( max_strand_energies <=  contrast_limit ) = 1   ;

% interpolate vectors in z to be consistent with xy
resolution_factors = resolution_factor * lumen_radius_in_pixels_range( 1, 1 ) ./ lumen_radius_in_pixels_range( 1, : );

y_limits = round(( y_limits - 1 ) * resolution_factors( 1 )) + 1 ;
x_limits = round(( x_limits - 1 ) * resolution_factors( 2 )) + 1 ;
z_limits = round(( z_limits - 1 ) * resolution_factors( 3 )) + 1 ;

% double the number of entries in the size look up table ( resolution_factor( 4 ) == 2 )
[ ~, lumen_radius_in_pixels_range, strand_subscripts ] ...
                                    = resample_vectors( lumen_radius_in_pixels_range, [ resolution_factors, 2 ], strand_subscripts, 0 );
                                
strand_subscripts = cellfun( @round, strand_subscripts, 'UniformOutput', false );

% precacalculate the spherical elements
number_of_scales = size( lumen_radius_in_pixels_range, 1 );

structuring_element_linear_indexing = cell( number_of_scales, 1 );
coloring_element_linear_indexing    = cell( number_of_scales, 1 );

size_of_inner_crop = [ y_limits( 2 ) - y_limits( 1 ) + 1, ...
                       x_limits( 2 ) - x_limits( 1 ) + 1, ...
                       z_limits( 2 ) - z_limits( 1 ) + 1  ];
                   
outer_crop_margin = round( lumen_radius_in_pixels_range( end, 1 ));                   
                   
size_of_outer_crop = size_of_inner_crop + 2 * outer_crop_margin ;

for scale_index = 1 : number_of_scales

    % find all pixel locations within the ellipsoid radii from the vertex position
    structuring_element_linear_indexing{ scale_index }                                                        ...
        = construct_structuring_element_V190( lumen_radius_in_pixels_range( scale_index, : ),     size_of_outer_crop );
    
    coloring_element_linear_indexing{ scale_index }                                                           ...
        = construct_structuring_element_V190( lumen_radius_in_pixels_range( scale_index, : ) + 1, size_of_outer_crop );
            
end % scale FOR

% loop through the strands to build the contrast and strand indexed strands images
number_of_strands     = numel( strand_subscripts );

strand_index_range = 1 : number_of_strands ;

number_of_image_voxels = prod( size_of_outer_crop );

strands_image           = zeros( size_of_outer_crop );
strand_index_image      = zeros( size_of_outer_crop );

for strand_index = strand_index_range
    
    subscripts_at_strand = strand_subscripts{ strand_index };
    
    subscripts_at_strand( :, 1 ) = subscripts_at_strand( :, 1 ) - y_limits( 1 ) + 1 + outer_crop_margin ;
    subscripts_at_strand( :, 2 ) = subscripts_at_strand( :, 2 ) - x_limits( 1 ) + 1 + outer_crop_margin ;
    subscripts_at_strand( :, 3 ) = subscripts_at_strand( :, 3 ) - z_limits( 1 ) + 1 + outer_crop_margin ;
    
    % trim the vectors whose center positions leave the crop
    within_bounds_strands_logical = subscripts_at_strand( :, 1 ) >= 1                       + outer_crop_margin ...
                                  & subscripts_at_strand( :, 1 ) <= size_of_outer_crop( 1 ) - outer_crop_margin ...
                                  & subscripts_at_strand( :, 2 ) >= 1                       + outer_crop_margin ...
                                  & subscripts_at_strand( :, 2 ) <= size_of_outer_crop( 2 ) - outer_crop_margin ...
                                  & subscripts_at_strand( :, 3 ) >= 1                       + outer_crop_margin ...
                                  & subscripts_at_strand( :, 3 ) <= size_of_outer_crop( 3 ) - outer_crop_margin ;

%     % trim the vectors that lie completely outside the crop
% 
%     radii_in_pixels_at_strand = radii_in_pixels_range( subscripts_at_strand( :, 4 ), : );
% 
%     within_bounds_strands_logical = subscripts_at_strand( :, 1 ) >= 1                 - radii_in_pixels_at_strand( :, 1 ) ...
%                                 & subscripts_at_strand( :, 1 ) <= size_of_crop( 1 ) + radii_in_pixels_at_strand( :, 1 ) ...
%                                 & subscripts_at_strand( :, 2 ) >= 1                 - radii_in_pixels_at_strand( :, 2 ) ...
%                                 & subscripts_at_strand( :, 2 ) <= size_of_crop( 2 ) + radii_in_pixels_at_strand( :, 2 ) ...
%                                 & subscripts_at_strand( :, 3 ) >= 1                 - radii_in_pixels_at_strand( :, 3 ) ...
%                                 & subscripts_at_strand( :, 3 ) <= size_of_crop( 3 ) + radii_in_pixels_at_strand( :, 3 ) ;
    
    subscripts_at_strand = subscripts_at_strand( within_bounds_strands_logical, : );
                                             
    number_of_strand_positions = size( subscripts_at_strand, 1 );
    
    for strand_position_index = 1 : number_of_strand_positions
    
        % find the linear index of the center pixel of the sphere that defines the strand at this
        % position
        strand_position_linear_index = sub2ind( size_of_outer_crop,                           ...
                                              subscripts_at_strand( strand_position_index, 1 ), ...
                                              subscripts_at_strand( strand_position_index, 2 ), ...
                                              subscripts_at_strand( strand_position_index, 3 )  );        
        
        % label the spheres and centerlines in an overwriting fashion so that only the lowest energy
        % strands will shine through in a multiply labeled region.  (Remember that we sorted by max
        % energy attained before this for loop so the later strands to be written are lower in
        % energy).  
        
               strands_image( max( min(   strand_position_linear_index                                                             ...
                                      + structuring_element_linear_indexing{ subscripts_at_strand( strand_position_index, 4 )},    ...
                                      number_of_image_voxels ),                                                                ...
                                 1                                                                                          )) ...
                                                                                       =       max_strand_energies( strand_index );        
        
        strand_index_image( max( min(   strand_position_linear_index                                                             ...
                                      +    coloring_element_linear_indexing{ subscripts_at_strand( strand_position_index, 4 )},    ...
                                      number_of_image_voxels ),                                                                ...
                                 1                                                                                          )) ...
                                                                                       =  max_strand_energies( strand_index )                                                                      ...
                                                                                       * strand_index ;
        
        % note consider replacing the strand index image with a contrast image so that the colormap
        % depends on the contrast
        
    end % for position
end % for trajectory

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

% flip the z dimension to be consistent with the normal orientation of the brain
strands_image      = flip(      strands_image, 3 );
strand_index_image = flip( strand_index_image, 3 );

% flip the y dimension to be consistent with the other z projected special outputs of vectorize()
strands_image      = flip(      strands_image, 1 );
strand_index_image = flip( strand_index_image, 1 );

% strands_image      = flip(      strands_image, 2 );
% strand_index_image = flip( strand_index_image, 2 );

number_of_strands     = numel( strand_subscripts );

if is_box_in_center

    random_colormap = 0.1 + 0.8 * rand( number_of_strands + 1, 3 ); % default matlab colormaps have 64 values

    % % set the first entry to gray and reserve it for the box
    random_colormap( 1, : ) = 0.9 ;

    strand_index_image = strand_index_image + 1 ;

else
    
    random_colormap = 0.1 + 0.8 * rand( number_of_strands, 3 ); % default matlab colormaps have 64 values

end

figure

colormap( random_colormap )

strands_image = strands_image( 1 + outer_crop_margin : size_of_outer_crop( 1 ) - outer_crop_margin, ...
                               1 + outer_crop_margin : size_of_outer_crop( 2 ) - outer_crop_margin, ...
                               1 + outer_crop_margin : size_of_outer_crop( 3 ) - outer_crop_margin  );
                       
strand_index_image ...
     = strand_index_image( 1 + outer_crop_margin : size_of_outer_crop( 1 ) - outer_crop_margin, ...
                           1 + outer_crop_margin : size_of_outer_crop( 2 ) - outer_crop_margin, ...
                           1 + outer_crop_margin : size_of_outer_crop( 3 ) - outer_crop_margin  );

% begin citation: https://www.mathworks.com/help/matlab/ref/isosurface.html SAM 6/14/18
[ F, V ] = isosurface( strands_image, 0.5 );

% convert color coding from vertices (for interp) to faces (for 'flat' FaceColor)

% do nearest neighbor interpolation (instead of whatever isosurface does without option)
C = interp3( strand_index_image, V( :, 1 ), V( :, 2 ), V( :, 3 ), 'nearest' );

% C = C( F( :, 1 ));    
C = median( C( F ), 2 );    

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

    % make objects transparent
    p = patch('Vertices',[ ver; V ],'Faces',[ fac; F ],'FaceVertexCData', [ color; C ],...
              'FaceColor','flat', 'EdgeColor','none', ...
              'FaceAlpha', 'flat', 'AlphaDataMapping', 'none', ...
              'FaceVertexAlphaData', [ 0.25 * ones( size( color )); 1 * ones( size( C ))]);

          % % unifrom transparency
%     p = patch('Vertices',[ ver; V ],'Faces',[ fac; F ],'FaceVertexCData', [ color; C ],...
%               'FaceColor','flat', 'EdgeColor','none', ...
%               'FaceAlpha', 0.5 );          
else
    
%     p = patch('Vertices',V ,'Faces',F ,'FaceVertexCData', C ,...
%               'FaceColor','flat', 'EdgeColor','none');

    p = patch('Vertices',V ,'Faces',F ,'FaceVertexCData', C ,...
              'FaceColor', 'flat', 'FaceAlpha',0.75, 'EdgeColor','none');

end

isonormals( strands_image,p)

% convert the vertex coordinates to microns
p.Vertices      = ( p.Vertices      - 1 ) / resolution_factor * microns_per_pixel_xy ;
% p.VertexNormals =   p.VertexNormals       * resolution_factor / microns_per_pixel_xy ; % transpose of the inverse transformation is done to normals (simplifies to a purely scaling transformation)
% p.VertexNormals = p.VertexNormals ./ sum( p.VertexNormals .^ 2, 2 ) .^ 0.5 ;

% p.FaceColor = 'red';
% p.EdgeColor = 'none';
daspect([1 1 1])
view(0,0); 
axis tight
camlight right
material dull
lighting gouraud
% end citation: https://www.mathworks.com/help/matlab/ref/isosurface.html SAM 6/14/18

end % function