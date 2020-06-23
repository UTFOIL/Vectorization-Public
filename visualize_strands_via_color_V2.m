function visualize_strands_via_color( edge_subscripts, max_edge_energies, edge_strand_assignmnets,  ...
                                      pixels_per_sigma_range, sigma_per_size, directories,          ...
                                      original_handle, y_limits, x_limits, z_limits,                ...
                                      contrast_limit, intensity_limits                              )
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
 
% round the edge subscripts to integers 
edge_subscripts = cellfun( @( v ) round( v ), edge_subscripts, 'UniformOutput', false );

original_file = [ directories{ 4 }, original_handle ];  

% original_file_info = h5info( original_file );
% 
% size_of_image = original_file_info.Datasets.Dataspace.Size ;

size_of_crop = [ y_limits( 2 ) - y_limits( 1 ) + 1, ...
                 x_limits( 2 ) - x_limits( 1 ) + 1, ...
                 z_limits( 2 ) - z_limits( 1 ) + 1  ];
             
original_image = h52mat( original_file,       ...
                         [ y_limits( 1 ),     ...
                           x_limits( 1 ),     ...
                           z_limits( 1 ) ],   ...
                         [ size_of_crop( 1 ), ...
                           size_of_crop( 2 ), ...
                           size_of_crop( 3 )  ]); 

original_image_max_intensity_projected = double( max( original_image, [ ], 3 ));

% histogram( original_image_max_intensity_projected( : ))

% scale intensity values to 0 : 255
inverted_original_MIP                                                     ...
    = 255                                                                 ...
      - ( original_image_max_intensity_projected - intensity_limits( 1 )) ...
      / (         intensity_limits( 2 )          - intensity_limits( 1 )) ...
      * 255 ;                                                  

% round the MIP intensities
inverted_original_MIP( inverted_original_MIP < 0   ) = 0   ;
inverted_original_MIP( inverted_original_MIP > 255 ) = 255 ;
                    
centerlines_image = zeros([ 3, size_of_crop([ 1, 2 ])]);
%     spheres_image = zeros([ 3, size_of_crop([ 1, 2 ])]);

% erasing edges that are above the upper energy limit
logical_edges_below_upper_energy_limit = max_edge_energies < contrast_limit;

max_edge_energies       =       max_edge_energies( logical_edges_below_upper_energy_limit );
edge_subscripts         =         edge_subscripts( logical_edges_below_upper_energy_limit );
edge_strand_assignmnets = edge_strand_assignmnets( logical_edges_below_upper_energy_limit );

% erasing edges that don't lie in the crop at least partially
max_subscripts = cell2mat( cellfun( @max,  edge_subscripts, 'UniformOutput', false ));
min_subscripts = cell2mat( cellfun( @min,  edge_subscripts, 'UniformOutput', false ));

radii_in_pixels_range = pixels_per_sigma_range * sigma_per_size ;

radii_in_pixels = radii_in_pixels_range( round( max_subscripts( :, 4 )), : );

logical_edges_in_crop = max_subscripts( :, 1 ) + radii_in_pixels( :, 1 ) >= y_limits( 1 ) ...
                      & max_subscripts( :, 2 ) + radii_in_pixels( :, 2 ) >= x_limits( 1 ) ...
                      & max_subscripts( :, 3 ) + radii_in_pixels( :, 3 ) >= z_limits( 1 ) ...
                      & min_subscripts( :, 1 ) - radii_in_pixels( :, 1 ) <= y_limits( 2 ) ...
                      & min_subscripts( :, 2 ) - radii_in_pixels( :, 2 ) <= x_limits( 2 ) ...
                      & min_subscripts( :, 3 ) - radii_in_pixels( :, 3 ) <= z_limits( 2 ) ;
                      
max_edge_energies       =       max_edge_energies( logical_edges_in_crop );
edge_subscripts         =         edge_subscripts( logical_edges_in_crop );
edge_strand_assignmnets = edge_strand_assignmnets( logical_edges_in_crop );

% sorting trajectories by max energy in descending order (most negative at end)
[ max_edge_energies, indices_sorted_by_mean ] = sort( max_edge_energies, 'descend' );

edge_subscripts         =         edge_subscripts( indices_sorted_by_mean );
edge_strand_assignmnets = edge_strand_assignmnets( indices_sorted_by_mean );

% histogram( max_edge_energies );

[ unique_strand_indices, ~, indices_of_unique_strands ] = unique( edge_strand_assignmnets );

number_of_strands     = length( unique_strand_indices );

% offset the indexing so that it starts at one instead of zero for the junctions
edge_strand_assignmnets = edge_strand_assignmnets + 1 ;

colors_at_strand = 0.1 + 0.8 * rand( 3, number_of_strands );

colors_at_strand( :, 1 ) = 1 ;

% set the max_edge_energies less than the contrast limits to zero and those above to 255
max_edge_energies( max_edge_energies >  contrast_limit ) = 0   ;
max_edge_energies( max_edge_energies <= contrast_limit ) = 255 ;
  
number_of_scales = size( pixels_per_sigma_range, 1 );

structuring_element_linear_indexing = cell( number_of_scales, 1 );

for scale_index = 1 : number_of_scales

    % find all pixel locations within the ellipsoid radii from the vertex position
    [ ~, structuring_element_linear_indexing{ scale_index }]                                        ...
        = construct_structuring_element_V180( radii_in_pixels_range( scale_index, : ), size_of_crop );
    
    % erase redundant indices due to the smushing of the 3rd dimension
    structuring_element_linear_indexing{ scale_index }               ...
        = unique( structuring_element_linear_indexing{ scale_index });
        
end % scale FOR

number_of_edges = length( max_edge_energies );

edge_index_range = 1 : number_of_edges ;

for edge_index = edge_index_range
    
    subscripts_at_edge = edge_subscripts{ edge_index };
    
    subscripts_at_edge( :, 1 ) = subscripts_at_edge( :, 1 ) - y_limits( 1 ) + 1 ;
    subscripts_at_edge( :, 2 ) = subscripts_at_edge( :, 2 ) - x_limits( 1 ) + 1 ;
    subscripts_at_edge( :, 3 ) = subscripts_at_edge( :, 3 ) - z_limits( 1 ) + 1 ;    
    
%     radii_in_pixels_at_edge = radii_in_pixels_range( subscripts_at_edge( :, 4 ), : );
                           
    % trim the vectors whose center positions leave the crop
    within_bounds_edges_logical = subscripts_at_edge( :, 1 ) >= 1                 ...
                                & subscripts_at_edge( :, 1 ) <= size_of_crop( 1 ) ...
                                & subscripts_at_edge( :, 2 ) >= 1                 ...
                                & subscripts_at_edge( :, 2 ) <= size_of_crop( 2 ) ...
                                & subscripts_at_edge( :, 3 ) >= 1                 ...
                                & subscripts_at_edge( :, 3 ) <= size_of_crop( 3 ) ;

    subscripts_at_edge = subscripts_at_edge( within_bounds_edges_logical, : );
    
%     colors_at_edge ...
%         = [   ( subscripts_at_edge( :, 3 )' -       1          ) / ( size_of_crop( 3 ) - 1 ); ...
%             - ( subscripts_at_edge( :, 3 )' - size_of_crop( 3 )) / ( size_of_crop( 3 ) - 1 ); ...
%                                  zeros( size( subscripts_at_edge( :, 3 )' ))                  ];
%                              
%     colors_at_edge                                                          ...
%         = [                     ones( size( subscripts_at_edge( :, 3 )' )); ...  
%            ( subscripts_at_edge( :, 3 )' - 1 ) / ( size_of_crop( 3 ) - 1 ); ...
%            ( subscripts_at_edge( :, 3 )' - 1 ) / ( size_of_crop( 3 ) - 1 )  ];
%                              
%     % round each color channel to be in [0, 1]
%     colors_at_edge = max( colors_at_edge, zeros( size( colors_at_edge )));
%     colors_at_edge = min( colors_at_edge,  ones( size( colors_at_edge )));    

%     colors_at_edge = ones( size( subscripts_at_edge( :, 3 )' )) .* colors_at_strand( :, edge_strand_assignmnets( edge_index ));

    number_of_edge_positions = size( subscripts_at_edge, 1 );
    
    % aspect_ratio must be a row vector here:
    
    for edge_position_index = 1 : number_of_edge_positions
    
        % find the linear index of the center pixel of the sphere that defines the edge at this
        % position
        edge_position_linear_index = sub2ind( size_of_crop([ 1, 2 ]),                       ...
                                              subscripts_at_edge( edge_position_index, 1 ), ...
                                              subscripts_at_edge( edge_position_index, 2 )  );        
        
        % label the spheres and centerlines in an overwriting fashion so that only the lowest energy
        % edges will shine through in a multiply labeled region.  (Remember that we sorted by max
        % energy attained before this for loop so the later edges to be written are lower in
        % energy).  Truncating the spheres at the boundaries of the 1D indexing only (not any 3D
        % considerations, so expect wraparound effects at most image boundaries).
%         centerlines_image( :, edge_position_linear_index ) ...
%             = max_edge_energies( edge_index )              ...
%             * colors_at_edge( :, edge_position_index );
        
        centerlines_image( :, edge_position_linear_index )               ...
            = max_edge_energies( edge_index )                            ...        
            * colors_at_strand( :, indices_of_unique_strands( edge_index ));
        
    end % for position
end % for trajectory

% get the initial locations of vessels from the logical of the 2D spheres image in the first or
% second color coordinate.
% centerlines_mask = squeeze( or( centerlines_image( 1, :, : ), centerlines_image( 2, :, : )));
% 
% inverted_original_MIP( centerlines_mask ) = 0 ;

centerlines_image = permute( centerlines_image, [ 2, 3, 1 ]);

composite_image = uint8( centerlines_image + inverted_original_MIP );

figure

image( composite_image )

end % function
