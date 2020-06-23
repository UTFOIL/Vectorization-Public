function visualize_strands_via_color_V200( strand_subscripts, vessel_directions, max_strand_energies, ...
                                           pixels_per_sigma_range, sigma_per_size, data_directory,          ...
                                           original_handle, y_limits, x_limits, z_limits,                ...
                                           intensity_limits, is_inverted_original, are_vectors_opaque, color_code, microns_per_voxel )
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
% V2, fixed the cropping bug that changed the size of the edges_subscripts array during cropping
% (but didn't change the pixels_per_size array, leading to a size mismatch at time of comparison).
% SAM 8/14/17
%
% V200: replaced directories input with data_directory SAM 12/7/18
 
% round the strand subscripts to integers 
strand_subscripts = cellfun( @( v ) round( v ), strand_subscripts, 'UniformOutput', false );

original_file = [ data_directory, original_handle ];  

% original_file_info = h5info( original_file );
% 
% size_of_image = original_file_info.Datasets.Dataspace.Size ;

% % !!!!!!!!! kludgey
z_intra  = 2 ;
z_offset = 0 ;

z_limits = z_limits + z_offset ;

% z_intra  = 0 ;
% z_offset = 0 ;

size_of_crop = [ y_limits( 2 ) - y_limits( 1 ) + 1, ...
                 x_limits( 2 ) - x_limits( 1 ) + 1, ...
                 z_limits( 2 ) - z_limits( 1 ) + 1  ];
                          
% original_image = h52mat( original_file,       ...
%                          [ y_limits( 1 ),     ...
%                            x_limits( 1 ),     ...
%                            z_limits( 1 ) ],   ...
%                          [ size_of_crop( 1 ), ...
%                            size_of_crop( 2 ), ...
%                            size_of_crop( 3 )  ]);

original_image = h52mat( original_file,       ...
                         [ y_limits( 1 ),     ...
                           x_limits( 1 ),     ...
                           z_limits( 1 ) + z_intra + z_offset ],   ...
                         [ size_of_crop( 1 ), ...
                           size_of_crop( 2 ), ...
                           size_of_crop( 3 ) - 2 * z_intra ]); 

original_MIP = double( max( original_image, [ ], 3 ));

% histogram( original_image_max_intensity_projected( : ))

% scale intensity values to 0 : 255
if is_inverted_original

    original_MIP                                                     ...
        = 255                                                                 ...
          - ( original_MIP - intensity_limits( 1 )) ...
          / (         intensity_limits( 2 )          - intensity_limits( 1 )) ...
          * 255 ;                                                  
      
else
    
    original_MIP                                                     ...
        =   ( original_MIP - intensity_limits( 1 )) ...
          / (         intensity_limits( 2 )          - intensity_limits( 1 )) ...
          * 255 ;                                                  
    
end

% round the MIP intensities
original_MIP( original_MIP < 0   ) = 0   ;
original_MIP( original_MIP > 255 ) = 255 ;
                    
centerlines_image = zeros([ 3, size_of_crop([ 1, 2 ])]);
%     spheres_image = zeros([ 3, size_of_crop([ 1, 2 ])]);
indicator_image   = zeros([ 1, size_of_crop([ 1, 2 ])]) + eps ;

number_of_image_pixels = prod( size_of_crop([ 1, 2 ]));

% % erasing strands that are above the upper energy limit
% logical_strands_below_upper_energy_limit = max_strand_energies <= contrast_limit;
% 
% max_strand_energies       =       max_strand_energies( logical_strands_below_upper_energy_limit );
% strand_subscripts         =         strand_subscripts( logical_strands_below_upper_energy_limit );
% vessel_directions         =         vessel_directions( logical_strands_below_upper_energy_limit );

% erasing strands that don't lie in the crop at least partially
max_subscripts = cell2mat( cellfun( @max,  strand_subscripts, 'UniformOutput', false ));
min_subscripts = cell2mat( cellfun( @min,  strand_subscripts, 'UniformOutput', false ));

radii_in_pixels_range = pixels_per_sigma_range * sigma_per_size ;

max_radii_in_pixels = radii_in_pixels_range( round( max_subscripts( :, 4 )), : );

logical_strands_in_crop = max_subscripts( :, 1 ) + max_radii_in_pixels( :, 1 ) >= y_limits( 1 ) ...
                        & max_subscripts( :, 2 ) + max_radii_in_pixels( :, 2 ) >= x_limits( 1 ) ...
                        & max_subscripts( :, 3 ) + max_radii_in_pixels( :, 3 ) >= z_limits( 1 ) ...
                        & min_subscripts( :, 1 ) - max_radii_in_pixels( :, 1 ) <= y_limits( 2 ) ...
                        & min_subscripts( :, 2 ) - max_radii_in_pixels( :, 2 ) <= x_limits( 2 ) ...
                        & min_subscripts( :, 3 ) - max_radii_in_pixels( :, 3 ) <= z_limits( 2 ) ;
                      
max_strand_energies       =       max_strand_energies( logical_strands_in_crop );
strand_subscripts         =         strand_subscripts( logical_strands_in_crop );
vessel_directions         =         vessel_directions( logical_strands_in_crop );

% % sorting trajectories by max energy in descending order (most negative at end)
% [ max_strand_energies, indices_sorted_by_mean ] = sort( max_strand_energies, 'descend' );

% sorting trajectories by average radius in ascending order (biggest at end)
average_radii = cellfun( @( x ) mean( x( :, 4 )), strand_subscripts );

[ ~, indices_sorted_by_size ] = sort( average_radii, 'ascend' );

max_strand_energies       =       max_strand_energies( indices_sorted_by_size );
strand_subscripts         =         strand_subscripts( indices_sorted_by_size );
vessel_directions         =         vessel_directions( indices_sorted_by_size );

% histogram( max_strand_energies );

number_of_strands     = numel( strand_subscripts );

% colors_at_strand( :, 1 ) = 1 ;

% % set the max_strand_energies less than the contrast limits to zero and those above to 255
% max_strand_energies( max_strand_energies >  contrast_limit ) = 0   ;
% max_strand_energies( max_strand_energies <= contrast_limit ) = 255 ;
  
number_of_scales = size( pixels_per_sigma_range, 1 );

structuring_element_linear_indexing = cell( number_of_scales, 1 );

% z-project all structuring elements:
radii_in_pixels_range( :, 3 ) = eps ;

for scale_index = 1 : number_of_scales

    % find all pixel locations within the ellipsoid radii from the vertex position
    structuring_element_linear_indexing{ scale_index }                                              ...
        = construct_structuring_element_V190( radii_in_pixels_range( scale_index, : ), size_of_crop );
    
%     % erase redundant indices due to the smushing of the 3rd dimension
%     structuring_element_linear_indexing{ scale_index }               ...
%         = unique( structuring_element_linear_indexing{ scale_index });
        
end % scale FOR

structuring_element_numels = cellfun( @numel, structuring_element_linear_indexing );

strand_index_range = 1 : number_of_strands ;

for strand_index = strand_index_range
    
    subscripts_at_strand = strand_subscripts{ strand_index };
    
    subscripts_at_strand( :, 1 ) = subscripts_at_strand( :, 1 ) - y_limits( 1 ) + 1 ;
    subscripts_at_strand( :, 2 ) = subscripts_at_strand( :, 2 ) - x_limits( 1 ) + 1 ;
    subscripts_at_strand( :, 3 ) = subscripts_at_strand( :, 3 ) - z_limits( 1 ) + 1 ;    
    
%     radii_in_pixels_at_strand = radii_in_pixels_range( subscripts_at_strand( :, 4 ), : );
                           
    % trim the vectors whose center positions leave the crop
    within_bounds_strands_logical = subscripts_at_strand( :, 1 ) >= 1                 ...
                                  & subscripts_at_strand( :, 1 ) <= size_of_crop( 1 ) ...
                                  & subscripts_at_strand( :, 2 ) >= 1                 ...
                                  & subscripts_at_strand( :, 2 ) <= size_of_crop( 2 ) ...
                                  & subscripts_at_strand( :, 3 ) >= 1                 ...
                                  & subscripts_at_strand( :, 3 ) <= size_of_crop( 3 ) ;

    subscripts_at_strand =              subscripts_at_strand( within_bounds_strands_logical, : );
    directions_at_strand = vessel_directions{ strand_index }( within_bounds_strands_logical, : );
    
    switch color_code

        case 'strands'

            colors_at_strand =  ones( size( subscripts_at_strand( :, 1 )' )) ...
                             .* ( 0.1 + 0.8 * rand( 3, 1 ));            

        case 'depth'

        %     colors_at_strand ...
        %         = [   ( subscripts_at_strand( :, 3 )' -       1          ) / ( size_of_crop( 3 ) - 1 ); ...
        %             - ( subscripts_at_strand( :, 3 )' - size_of_crop( 3 )) / ( size_of_crop( 3 ) - 1 ); ...
        %                                  zeros( size( subscripts_at_strand( :, 3 )' ))                  ];

%             colors_at_strand                                                          ...
%                 = [                     ones( size( subscripts_at_strand( :, 1 )' )); ...  
%                    ( subscripts_at_strand( :, 3 )' - 1 ) / ( size_of_crop( 3 ) - 1 ); ...
%                    ( subscripts_at_strand( :, 3 )' - 1 ) / ( size_of_crop( 3 ) - 1 )  ];

            colors_at_strand                                                                       ...
                = [             ( subscripts_at_strand( :, 3 )' - 1 ) / ( size_of_crop( 3 ) - 1 ); ...  
                     zeros( size( subscripts_at_strand( :, 1 )' ));                                ...
                            1 - ( subscripts_at_strand( :, 3 )' - 1 ) / ( size_of_crop( 3 ) - 1 )  ];

            % round each color channel to be in [0, 1]
            colors_at_strand = max( colors_at_strand, zeros( size( colors_at_strand )));
            colors_at_strand = min( colors_at_strand,  ones( size( colors_at_strand )));    

        case 'directions'
            
            % y, x, z -> G, R, B
            colors_at_strand = abs( directions_at_strand( :, [ 2, 1, 3 ])' );
                        
        case 'z-directions'
            
            % y, x, z -> R, R, B
            colors_at_strand = abs( directions_at_strand( :, [ 2, 1, 3 ])' );
            
            colors_at_strand( 1, : ) = sum( colors_at_strand([ 1, 2 ], : ) .^ 2, 1 ) .^ 0.5 ;

            colors_at_strand( 2, : ) = zeros( size( subscripts_at_strand( :, 1 )' ));

    end

    number_of_strand_positions = size( subscripts_at_strand, 1 );
    
    % aspect_ratio must be a row vector here:
    
    for strand_position_index = 1 : number_of_strand_positions
    
        % find the linear index of the center pixel of the sphere that defines the strand at this
        % position
        strand_position_linear_index = sub2ind( size_of_crop([ 1, 2 ]),                       ...
                                              subscripts_at_strand( strand_position_index, 1 ), ...
                                              subscripts_at_strand( strand_position_index, 2 )  );        
        
        % label the spheres and centerlines in an overwriting fashion so that only the largest
        % strands will shine through in a multiply labeled region.  (starnds are sorted by size
        % before this FOR).  Truncating the spheres at the boundaries of the 1D indexing only (not
        % any 3D considerations, so expect wraparound effects at most image boundaries).
        
%         centerlines_image( :, strand_position_linear_index )                 ...
%             = max_strand_energies( strand_index )                            ...        
%             * colors_at_strand( :, strand_position_index );
        
        centerlines_image( :, max( min( strand_position_linear_index                                                                ...
                                    + structuring_element_linear_indexing{ subscripts_at_strand( strand_position_index, 4 )},       ...
                                 number_of_image_pixels                                                                       ),    ...
                            1                                                                                                    )) ...
      = centerlines_image( :, max( min( strand_position_linear_index                                                                ...
                                    + structuring_element_linear_indexing{ subscripts_at_strand( strand_position_index, 4 )},       ...
                                 number_of_image_pixels                                                                       ),    ...
                            1                                                                                                    )) ...
      + colors_at_strand( :, strand_position_index ) .* ones( 3, structuring_element_numels( subscripts_at_strand( strand_position_index, 4 )));
        
          indicator_image(    max( min( strand_position_linear_index                                                                ...
                                    + structuring_element_linear_indexing{ subscripts_at_strand( strand_position_index, 4 )},       ...
                                 number_of_image_pixels                                                                       ),    ...
                            1                                                                                                    )) ...
      =   indicator_image(    max( min( strand_position_linear_index                                                                ...
                                    + structuring_element_linear_indexing{ subscripts_at_strand( strand_position_index, 4 )},       ...
                                 number_of_image_pixels                                                                       ),    ...
                            1                                                                                                    )) ...
      +                                                 ones(    structuring_element_numels( subscripts_at_strand( strand_position_index, 4 )), 1 );
        
        
    end % for position
end % for strand

if strcmp( color_code, 'directions' )
    
    figure

    % from https://www.mathworks.com/help/matlab/ref/isocolors.html
    radius_of_sphere = 20 ;

    data_range_1D = linspace( 0, radius_of_sphere, radius_of_sphere );
    [x,y,z] = meshgrid(data_range_1D,data_range_1D,data_range_1D);
    data = sqrt(x.^2 + y.^2 + z.^2);
    p = patch(isosurface(x,y,z,data,radius_of_sphere));
%             isonormals(x,y,z,data,p)
    color_range_1D = 0.5 + linspace( 0, 0.5, radius_of_sphere );
    [r,g,b] = meshgrid(color_range_1D,color_range_1D,color_range_1D);
    isocolors(x,y,z,r,g,b,p)
    p.FaceColor = 'interp';
    p.EdgeColor = 'none';
    view(135,45) 
    daspect([1 1 1])
    
%     xticks([radius_of_sphere])
%     yticks([radius_of_sphere])
%     zticks([radius_of_sphere])
    
%     xticklabels(['x'])
%     yticklabels(['y'])
%     zticklabels(['z'])
    
    set(gca,'Visible','off')
    
    radius_of_text_sphere = radius_of_sphere * ( 1 + 0.1 ) ;
    
    text(radius_of_text_sphere, 0                    , 0                    , 'x', 'HorizontalAlignment', 'center', 'FontSize',18 )
    text(0                    , radius_of_text_sphere, 0                    , 'y', 'HorizontalAlignment', 'center', 'FontSize',18 )
    text(0                    , 0                    , radius_of_text_sphere, 'z', 'HorizontalAlignment', 'center', 'FontSize',18 )
        
%     xlabel('x')
%     ylabel('y')
%     zlabel('z')
    
%             camlight 
%             lighting gouraud            

end

centerlines_image = centerlines_image ./ indicator_image ;

% get the initial locations of vessels from the logical of the 2D spheres image in the first or
% second color coordinate.
% centerlines_mask = squeeze( or( centerlines_image( 1, :, : ), centerlines_image( 2, :, : )));
% 
% original_MIP( centerlines_mask ) = 0 ;

centerlines_image = permute( centerlines_image, [ 2, 3, 1 ]);

are_vectors_bright = false ;

indicator_image = indicator_image - eps ;

if are_vectors_opaque
    
    % make vectors opaque
    original_MIP( logical( indicator_image )) = 0 ; 

    composite_image = uint8(( centerlines_image * 255 + original_MIP ));

else

	weight = 0.5 ;    
    
    if are_vectors_bright
        
        composite_image = uint8(( weight + ( 1 - weight ) * centerlines_image ) .* original_MIP );

    else
       
        weighted_centerlines_image = ( 1 - weight ) * centerlines_image .* original_MIP ;        
        
        original_MIP( logical( indicator_image ))          ...
      = original_MIP( logical( indicator_image )) * weight ;
  
        
        composite_image = uint8( original_MIP + weighted_centerlines_image );

    end
        
end

figure

image( composite_image )

tick_spacing = 50 ; % microns

set( gca, 'XTick',      0 : tick_spacing / microns_per_voxel( 2 ) :   size_of_crop( 2 )                         );
set( gca, 'XtickLabel', 0 : tick_spacing                          :   size_of_crop( 2 ) * microns_per_voxel( 2 ));

set( gca, 'YTick',      0 : tick_spacing / microns_per_voxel( 1 ) :   size_of_crop( 1 )                         );
set( gca, 'YtickLabel', 0 : tick_spacing                          :   size_of_crop( 1 ) * microns_per_voxel( 1 ));

set( gca, 'TickLength', [ 0, 0 ])

xlabel('microns')

daspect([ 1 1 1 ])

switch color_code
    
    case 'depth'
        % make color bar
        map = 0.5                                                                                  ...
            + [ 0.5 *       ( linspace(1,size_of_crop( 3 ),64) - 1 )' / ( size_of_crop( 3 ) - 1 ), ...  
                 zeros( size( linspace(1,size_of_crop( 3 ),64)' )),                                ...
                0.5 * ( 1 - ( linspace(1,size_of_crop( 3 ),64) - 1 )' / ( size_of_crop( 3 ) - 1 )) ];
                    
        colormap(map)
        
        ticks_unitless = linspace( 0, 1, 4 );
        
        ticks_mincrons = - microns_per_voxel( 3 ) * ( z_limits( 1 ) - 1 + ( z_limits( 2 ) - z_limits( 1 )) * ticks_unitless )';
        
        ticks_microns_string = num2str( round( ticks_mincrons ));
        
        h = colorbar( 'Direction', 'reverse', 'Ticks',ticks_unitless,                                                 ...
                  'TickLabels',mat2cell( ticks_microns_string,ones(length(ticks_mincrons),1),size(ticks_microns_string,2 )));
             
%         ylabel( h, 'Depth' )
                            
    case  'z-directions'

        map( 1 : 64, 1 : 3 ) = linspace( 0, 1, 64 )';

        map( 1 : 64, 1     ) = sum( map( :, [ 1, 2 ]) .^ 2, 2 ) .^ 0.5 ;

        map( 1 : 64, 2     ) = 0 ;
        
        colormap(map)
        
        colorbar        
        
    case 'directions'
        % make color 1/8 sphere 
%         figure
%         
%         % make color bar
%         [x,y,z] = sphere(20);               %# Makes a 21-by-21 point sphere
%         x = x(11:end,:);                %# Keep top 11 x points
%         y = y(11:end,:);                %# Keep top 11 y points
%         z = z(11:end,:);                %# Keep top 11 z points
%         r = 1;                          %# A radius value
%         hs = surf(r.*x,r.*y,r.*z);      %# Plot the surface
%         direction = [0 1 0];            % Specify Direction
%         rotate(hs, direction, 90)       % Rotate The Object (Hemisphere) in ‘Direction’ By 90°
%         axis equal;                     %# Make the scaling on the x, y, and z axes equal
%         Ax = get(gca);                  % Axes Handle
%         XD = Ax.Children.XData;         % Get ‘XData’
%         Ax.Children.XData = XD + 0.5;   % Add 0.5 To ‘XData’ To Shift It To All > 0

%         map = 0.5 ...
%             + [ 0.5 *       ( linspace(1,size_of_crop( 3 ),64) - 1 )' / ( size_of_crop( 3 ) - 1 ), ...  
%                  zeros( size( linspace(1,size_of_crop( 3 ),64)' )),                                ...
%                 0.5 *   1 - ( linspace(1,size_of_crop( 3 ),64) - 1 )' / ( size_of_crop( 3 ) - 1 )  ];
%         
%         colormap(map)
%         
%         colorbar
%         
%         map = [ 0.5 *       ( linspace(1,size_of_crop( 3 ),64) - 1 )' / ( size_of_crop( 3 ) - 1 ), ...  
%                  zeros( size( linspace(1,size_of_crop( 3 ),64)' )),                                ...
%                 0.5 * ( 1 - ( linspace(1,size_of_crop( 3 ),64) - 1 )' / ( size_of_crop( 3 ) - 1 )) ];
%         
%         colormap(map)
%         
%         colorbar        
        
end

set(gca,'FontSize',16)

end % function
