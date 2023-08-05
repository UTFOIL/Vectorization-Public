function get_energy_V202( matching_kernel_string, lumen_radius_in_microns_range,                    ...
                          vessel_wall_thickness_in_microns, microns_per_voxel,    ...
                          pixels_per_sigma_PSF, max_voxels_per_node, data_directory, original_handle,  ...
                          energy_handle, gaussian_to_ideal_ratio, spherical_to_annular_ratio )
%% Gaussian blur at many scales, chunk-wise and in parallel
%
% V02 in which the inputs are changed to be ammenable to the vectorize_V02 script 3/5/18 SAM
%
% V10: octaves structure is being erased.  no downsampling here.  Perhaps including the PSF and its
% deconvolution from the blurred images, see the note next to the gaussian filter function. SAM
% 3/28/18
%
% Function name changed from blur_V10 to get_energy_V140 SAM 5/5/18.  To be used with
% vectorize_V14. Computing the derivatives in the fourier space instead of the finite differencing
% in the spatial (voxel) domain.
%
% V141 in which the gradient images are also computed, from which the direction of the gradient
% vector iamge is created, all of the principal curvatures are then projected onto the plane that is
% orthogonal to the gradient before summing them up.  Goes with energy_filter_V142 for starters.
% SAM 5/7/18. Also made the chunk overlap a vector quantity and updated the get_starts_and_counts
% function to V140.
%
% V150 in which the energy is min projected across the scale dimension before saving into h5 file
% SAM 5/7/18
%
% V151 in which the symmetry factor ratio is an input to the energy filter function herein. SAM
% 5/9/18
% 
% V160 in which the symmetry factor ratio is adjusted based on the pixels per sigma, so it will
% change for different scales and will be anisotropic vector quantity too.  SAM 5/17/18
%
% V161 to be used with energy_filter V161 SAM 5/18/18
%
% V162 to be used with energy_filter V162 SAM 5/21/18  Added the feature that bottom and top scales
% are set to infinite energy SAM 5/28/18  Switched that feature back 5/30/18 0150 SAM Switched again
% 8/7/18
%
% V191 to be used with energy_filter V191.  The shape of the matching filter kernel is now an input
% to the energy filter.  SAM 8/10/18
%
% Also adjusted for the volume interpretation of the octave SAM 8/14/18 (see the construction of the
% lumen_radius_in_microns_range variable).
%
% V192 two copies of the energy function are made.  One for the vertices and one for the edges.  The
% vertices version is different in that the spots of largest and smallest scales are set to
% infinity.  This version is stored as the third component in the fourth dimension of that variable.
% SAM 8/26/18
%
% V193:  Back to one copy of the enrgy function image, we will kill the vertices at the size
% extremes after the size exclusion/curation steps.  Like before V192 but the order is different:
% Index is the first image and energy is the second. SAM 11/2/18
%
% V201: Imposes a maximum resolution for the larger scales, downsamples the data to that resolution
% and then interpolates back after the main computations.  The point is to save memory. SAM 12/8/18
%
% V202: inputs to match V200. SAM 12/11/18

% starting_resolutions = microns_per_pixel / lumen_radius_in_microns_range( 1 ); % unitless um/um
% 
% best_starting_resolution = min( starting_resolutions ) ; % unitless um/um
% 
% best_resolution_allowed  = best_starting_resolution / 10 ;

worst_resolution_to_downsample = 1 / 2.5 ;
% worst_resolution_to_downsample = 1 / 3.5 ;
% worst_resolution_to_downsample = 1 / 4.5 ;
% best_resolution_allowed = 1 / 9 ;
% best_resolution_allowed = 1 / 6 ;

% construct the directories for storing the blurred data (and the chunked intermediary data).
number_of_scales = length( lumen_radius_in_microns_range );

scales_per_octave = log(                                    2                                   ) ...
                  / log( lumen_radius_in_microns_range( 2 ) / lumen_radius_in_microns_range( 1 )) ...
                  / 3 ; % divide by three for the volume interperetation of octave

scale_subscripts_range = 1 : number_of_scales ;

octave_at_scales = ceil( scale_subscripts_range / scales_per_octave / 3 );  % times by three for the volume interperetation of octave
      
octave_range = unique( octave_at_scales );

original_file = [ data_directory, original_handle ];  
  energy_file = [ data_directory,   energy_handle ];
     
original_file_info = h5info( original_file );

size_of_image = original_file_info.Datasets.Dataspace.Size ;

size_of_energy_image = [ size_of_image, 2 ];

h5create( energy_file, '/d', size_of_energy_image )

pixels_per_radius_range         = lumen_radius_in_microns_range    ./ microns_per_voxel ;
vessel_wall_thickness_in_pixels = vessel_wall_thickness_in_microns ./ microns_per_voxel ;

% current_octave = ceil( s_subscript / scales_per_octave / 3 );
% 
% if current_octave > previous_octave % !!!!! move this IF outside the chunk PARFOR as an octave FOR !!! rechunk at every octave !! will be much more efficient and accurate for the filters at larger scales with more downsampling !! SAM 6/4/22
% 
%     previous_octave = current_octave ;
% 
% end

%% decrease octave list to just the ones with unique downsampling schemes
resolution_factors              = zeros( length( octave_range ), 3 );
largest_scale_at_current_octave = zeros(   size( octave_range )    );

for current_octave = octave_range

     largest_scale_at_current_octave( current_octave ) = min( number_of_scales, floor(  current_octave       * scales_per_octave * 3 )     );
    smallest_scale_at_current_octave                   = min( number_of_scales, floor(( current_octave - 1 ) * scales_per_octave * 3 ) + 1 );
    
    resolutions_at_octave = min( microns_per_voxel / lumen_radius_in_microns_range( smallest_scale_at_current_octave ), worst_resolution_to_downsample * ones( size( microns_per_voxel ))); % unitless um/um

    resolution_factors( current_octave, : ) = round( worst_resolution_to_downsample ./ resolutions_at_octave );    
    
end % octave FOR

[ ~, IA, IC ] = unique( resolution_factors, 'rows', 'last' );

octave_range = octave_range( IA );

octave_at_scales = octave_range( IC( octave_at_scales ));

%% octave (downsampling) FOR
% at the start of each octave (doubling of the size of object) recompute the dft of the data
% with a different overhang (offset) based on the largest size object at the current octave
for current_octave = octave_range

    chunk_directory = [ data_directory, energy_handle, '_chunks_octave_', num2str( current_octave ), '_of_', num2str( octave_range( end )), filesep ];

    mkdir( chunk_directory );    
    
%             largest_pixels_per_radius_at_octave = pixels_per_radius_range( largest_scale_at_current_octave, : ) ./ resolution_factors ;
    largest_pixels_per_radius_at_octave = pixels_per_radius_range( largest_scale_at_current_octave( current_octave ), : ); % SAM 9/24/19

    pixels_per_sigma_PSF_at_octave = pixels_per_sigma_PSF ./ resolution_factors( current_octave, : ) ;    
       microns_per_pixel_at_octave =  microns_per_voxel   .* resolution_factors( current_octave, : ) ;

	approx_size_of_image_at_octave = round( size_of_image ./ resolution_factors( current_octave, : ) );
      voxel_aspect_ratio_at_octave =                    1 ./ microns_per_pixel_at_octave ;
      
%       minimum_Gaussian_kernel = 2 * geomean( microns_per_pixel_at_octave ); % average radius of 3 voxels in each dimension
%        
% % 	number_of_voxels_in_largest_vessel_at_scale = 4 / 3 * pi * prod( largest_pixels_per_radius_at_octave );
% 	number_of_voxels_in_largest_scale_at_octave = prod( 2 * largest_pixels_per_radius_at_octave );
      
% 	max_voxels_per_node = ( 10 ) ^ 3 * number_of_voxels_in_largest_scale_at_octave ; % 10^3 = ten vessel diameters in each dimension
% 	max_voxels_per_node = 5e6 ; 
      
    [ chunk_lattice_dimensions, number_of_chunks ] ...
             = get_chunking_lattice_V190( voxel_aspect_ratio_at_octave, max_voxels_per_node, approx_size_of_image_at_octave );     

    chunk_index_range = 1 : number_of_chunks ;     

    [ ~, ~, ~, ~, ~, ~,                                                                         ...
      y_writing_starts, x_writing_starts, z_writing_starts,                                     ...
      y_writing_counts, x_writing_counts, z_writing_counts   ]                                  ...
            = get_starts_and_counts_V200( chunk_lattice_dimensions, [ 0, 0, 0 ], size_of_image, [ 1, 1, 1 ]);

    switch matching_kernel_string

        case { '3D gaussian', '3D gaussian conv spherical pulse', '3D gaussian conv annular pulse' }

 ...           chunk_overlap_vector = ceil( 10 * (   pixels_per_sigma_PSF .^ 2                   ... weird overlap artifact observed in the energy image. this solved this issue, but takes forever. 12/7/21 SAM    ...
                     chunk_overlap_vector = ceil( 6 * (   pixels_per_sigma_PSF .^ 2                   ... 12/9/21 !!!!!!!!! consider restructuring the for loop here for efficiency. THe larger scales should have different chunking lattices. permute the PARFOR and the FOR loops. Writing a file for every octave separately
...           chunk_overlap_vector = ceil( 6 * (   pixels_per_sigma_PSF .^ 2                         ...
                                               + largest_pixels_per_radius_at_octave .^ 2 ) .^ 0.5 );

        case 'spherical pulse'

            chunk_overlap_vector = ceil(   6 * pixels_per_sigma_PSF                ...
                                         + 2 * largest_pixels_per_radius_at_octave );

        case { 'annular pulse', 'annular pulse V2' }

            chunk_overlap_vector = ceil(   6 * pixels_per_sigma_PSF                ...
                                         + 6 * largest_pixels_per_radius_at_octave ...
                                         +    vessel_wall_thickness_in_pixels      ); 

        case { '3D annular gaussian', '3D annular gaussian V2' }

            chunk_overlap_vector = ceil(   6 * pixels_per_sigma_PSF                     ...
                                         + 4 * (    largest_pixels_per_radius_at_octave ...
                                                 +  vessel_wall_thickness_in_pixels     )); 

        case 'radial gaussian'

            chunk_overlap_vector = ceil(   6   * pixels_per_sigma_PSF                ...
                                         + 2   *  largest_pixels_per_radius_at_octave ...
                                         + 1.5 * vessel_wall_thickness_in_pixels     );

    end

    [ y_reading_starts, x_reading_starts, z_reading_starts,         ...
      y_reading_counts, x_reading_counts, z_reading_counts,         ...
      ~, ~, ~, ~, ~, ~,                                             ...
      y_offset,         x_offset,         z_offset,          ]      ...
            = get_starts_and_counts_V200( chunk_lattice_dimensions, ...
                                          chunk_overlap_vector,     ...
                                          size_of_image,            ...
                                          resolution_factors( current_octave, : )        );
    
    %% Main chunk PARFOR: 
    % Reads a chunk from the interpolated file that is larger than the volume of interest. Does
    % gaussian filtering at many scales.  Writes a smaller chunk to the blurred.
    parfor chunk_index = chunk_index_range

        [ y_chunk_index, x_chunk_index, z_chunk_index ] = ind2sub( chunk_lattice_dimensions, chunk_index );
        
        reading_counts  = [ y_reading_counts( y_chunk_index ), ...
                            x_reading_counts( x_chunk_index ), ...
                            z_reading_counts( z_chunk_index )  ];

        % Read original image chunk and interpolate it to downsample at higher octaves
        original_chunk   = h52mat( original_file,                         ...
                                   [ y_reading_starts( y_chunk_index ),   ...
                                     x_reading_starts( x_chunk_index ),   ...
                                     z_reading_starts( z_chunk_index ) ], ...
                                   1 + floor(( reading_counts - 1 ) ./ resolution_factors( current_octave, : ) ),  ...
                                   resolution_factors( current_octave, : )                     ); 

        chunk_dft = fourier_transform_V2( original_chunk );

    %             y_local_range = y_offset( y_chunk_index ) + ( 1 : y_writing_counts( y_chunk_index ));
    %             x_local_range = x_offset( x_chunk_index ) + ( 1 : x_writing_counts( x_chunk_index ));
    %             z_local_range = z_offset( z_chunk_index ) + ( 1 : z_writing_counts( z_chunk_index ));

        % extract a local range that is guaranteed to include the range to be written
        y_local_range = 1 + floor( y_offset( y_chunk_index )                                           / resolution_factors( current_octave, 1 )) ...
                      : 1 + ceil(( y_offset( y_chunk_index ) + y_writing_counts( y_chunk_index ) - 1 ) / resolution_factors( current_octave, 1 )) ;
        x_local_range = 1 + floor( x_offset( x_chunk_index )                                           / resolution_factors( current_octave, 2 )) ...
                      : 1 + ceil(( x_offset( x_chunk_index ) + x_writing_counts( x_chunk_index ) - 1 ) / resolution_factors( current_octave, 2 )) ;
        z_local_range = 1 + floor( z_offset( z_chunk_index )                                           / resolution_factors( current_octave, 3 )) ...
                      : 1 + ceil(( z_offset( z_chunk_index ) + z_writing_counts( z_chunk_index ) - 1 ) / resolution_factors( current_octave, 3 )) ;

    %             y_local_range = 1 + floor(  y_offset( y_chunk_index )                                      / resolution_factors( 1 )) ...
    %                           :     floor(( y_offset( y_chunk_index ) + y_writing_counts( y_chunk_index )) / resolution_factors( 1 )) ;
    %             x_local_range = 1 + floor(  x_offset( x_chunk_index )                                      / resolution_factors( 2 )) ...
    %                           :     floor(( x_offset( x_chunk_index ) + x_writing_counts( x_chunk_index )) / resolution_factors( 2 )) ;
    %             z_local_range = 1 + floor(  z_offset( z_chunk_index )                                      / resolution_factors( 3 )) ...
    %                           :     floor(( z_offset( z_chunk_index ) + z_writing_counts( z_chunk_index )) / resolution_factors( 3 )) ;

        local_ranges = { y_local_range, ...
                         x_local_range, ...
                         z_local_range  };            

        [ mesh_Y, mesh_X, mesh_Z ] =   ndgrid( linspace( 1 + mod( y_offset( y_chunk_index ), resolution_factors( current_octave, 1 )) / resolution_factors( current_octave, 1 ), ...
                                                         1 + mod( y_offset( y_chunk_index ), resolution_factors( current_octave, 1 )) / resolution_factors( current_octave, 1 )  ...
                                                                                          + ( y_writing_counts( y_chunk_index ) - 1 ) / resolution_factors( current_octave, 1 ), ...
                                                                                              y_writing_counts( y_chunk_index )),                                ... 
                                               linspace( 1 + mod( x_offset( x_chunk_index ), resolution_factors( current_octave, 2 )) / resolution_factors( current_octave, 2 ), ...
                                                         1 + mod( x_offset( x_chunk_index ), resolution_factors( current_octave, 2 )) / resolution_factors( current_octave, 2 )  ...
                                                                                          + ( x_writing_counts( x_chunk_index ) - 1 ) / resolution_factors( current_octave, 2 ), ...
                                                                                              x_writing_counts( x_chunk_index )),                                ... 
                                               linspace( 1 + mod( z_offset( z_chunk_index ), resolution_factors( current_octave, 3 )) / resolution_factors( current_octave, 3 ), ...
                                                         1 + mod( z_offset( z_chunk_index ), resolution_factors( current_octave, 3 )) / resolution_factors( current_octave, 3 )  ...
                                                                                          + ( z_writing_counts( z_chunk_index ) - 1 ) / resolution_factors( current_octave, 3 ), ...
                                                                                              z_writing_counts( z_chunk_index ))                                 );

        writing_counts  = [ y_writing_counts( y_chunk_index ), ...
                            x_writing_counts( x_chunk_index ), ...
                            z_writing_counts( z_chunk_index ), ...
                            2                                  ];

        scale_subscripts_range_at_octave = find( octave_at_scales == current_octave );
        
        number_of_scales_at_octave = length( scale_subscripts_range_at_octave );

        energy_chunk_4D = zeros([ writing_counts( 1 : 3 ), number_of_scales_at_octave ]);

        s_subscript_at_octave = 0 ;
        
        %% scale FOR
        for s_subscript = scale_subscripts_range_at_octave

            s_subscript_at_octave = s_subscript_at_octave + 1 ;
            
%             gaussian_to_ideal_ratio = min( 1, minimum_Gaussian_kernel / lumen_radius_in_microns_range( s_subscript ));
            
            energy_chunk_4D( :, :, :, s_subscript_at_octave )                                                                                ...
                = interp3( energy_filter_V200( chunk_dft, matching_kernel_string, lumen_radius_in_microns_range( s_subscript ),    ...
                                               vessel_wall_thickness_in_microns, microns_per_pixel_at_octave,                       ...
                                               pixels_per_sigma_PSF_at_octave, local_ranges, gaussian_to_ideal_ratio, spherical_to_annular_ratio, scales_per_octave ), ...
                                                                                                            mesh_X, mesh_Y, mesh_Z ); %#ok<PFBNS>

    %         scaling_factor = ( lumen_radius_in_microns_range( s_subscript ) ^ 3 / prod( microns_per_pixel_at_scale )) ^ 0.5 ; % scaling for standard error of the mean considerations (replacing standard deviation with the energy value, i.e. the LoG filtered image. The number of samples is proportional to the  radius cubed, while the image is being downsampled. Assume that each location in the LoG image is an upper estimate of the true value) SAM 11/5/21
    %         scaling_factor = ( lumen_radius_in_microns_range( s_subscript ) ^ 3 / prod( microns_per_pixel_at_scale )) ^ ( 1 / 3 ) ; % Assume a 2D shape % !!!!!!!!!!!!! not yet attempted
    %         scaling_factor = ( lumen_radius_in_microns_range( s_subscript ) ^ 3 / prod( microns_per_pixel_at_scale )) ^ ( 1 / 6 ) ; % Assume a 1D shape %  pretty good output. Over estimates for smaller vessels, robust to noise, slightly over for larger vessels
    %         scaling_factor = log( lumen_radius_in_microns_range( s_subscript ) ^ 3 / prod( microns_per_pixel_at_scale ));    
    %         scaling_factor = log( 4 / 3 * pi * lumen_radius_in_microns_range( s_subscript ) ^ 3 / prod( microns_per_pixel_at_scale ));      
    %         scaling_factor = log( 1 + 4 / 3 * pi * lumen_radius_in_microns_range( s_subscript ) ^ 3 / prod( microns_per_pixel_at_scale ));  % best so far 11/6/21
    %         scaling_factor = log( 1 + 4 / 3 * pi * ( lumen_radius_in_microns_range( s_subscript ) ^ 3 / prod( microns_per_pixel_at_scale )) ^ ( 2 / 3 ));  
    %         scaling_factor = log( 1 + 4 / 3 * pi * ( lumen_radius_in_microns_range( s_subscript ) ^ 3 / prod( microns_per_pixel_at_scale )) ^ ( 2.5 / 3 ));  

    % % no more scalnig SAM 12/8/21
    %         scaling_factor = log( 4 / 3 * pi * ( 2 * lumen_radius_in_microns_range( s_subscript )) ^ 3 / prod( microns_per_pixel_at_scale ));  % 2 times the radius for sample size considerations
    %         energy_chunk_4D( :, :, :, s_subscript ) = energy_chunk_4D( :, :, :, s_subscript ) * scaling_factor ;

        end % end scale FOR

        energy_chunk_name = [ int2str( chunk_index ), ' of ', int2str( number_of_chunks )];

        energy_chunk_file = [ chunk_directory, energy_chunk_name ];

        h5create( energy_chunk_file, '/d', writing_counts )

    % %     projection_method = 'min' ; % SAM 9/3/19
    %     projection_method = 'sum' ;
    %     
    %     switch projection_method
    %         
    %         case 'min'
    %             
    %             [ energy_chunk, energy_chunk_scale_indices ] = min( energy_chunk_4D, [ ], 4 );
    %             
    %         case 'sum'
    %             
    %             energy_chunk_4D( energy_chunk_4D( : ) > 0 ) = 0 ;
    %             
    %             energy_chunk = sum( energy_chunk_4D, 4 );            
    % 
    % %             energy_chunk_4D( energy_chunk_4D > 0 ) = 0 ;
    %             
    %             energy_chunk_scale_indices                                                       ...
    %                                  = sum(    energy_chunk_4D                                   ...
    %                                         .* permute( 1 : number_of_scales, [ 1, 3, 4, 2 ]),   ...
    %                                         4                                                 )  ...
    %                                    ./ sum( energy_chunk_4D, 4 )                              ;  
    %             
    %             energy_chunk_scale_indices( isnan( energy_chunk_scale_indices )) = 1 ;
    %             
    %             energy_chunk = energy_chunk ./ number_of_scales ;
    %             
    %         case 'conv'
    %             
    %             % blur the energy chunk along the size dimension before taking minimum
    %             octaves_per_sigma = 1 ;
    %             
    %             energy_chunk_4D( energy_chunk_4D( : ) > 0 ) = 0 ;
    %             
    %             scales_per_sigma = octaves_per_sigma * scales_per_octave ;
    %                         
    %             kernel_scale_domain = ( 1 : number_of_scales )' ;
    %             
    %             kernel_scale_domain_offset = 1 : number_of_scales ;
    %             
    %             kernel_scale_domains = kernel_scale_domain        ...
    %                                  - kernel_scale_domain_offset ;
    %             
    %             kernel_sigma_domains = kernel_scale_domains / scales_per_sigma ;
    %             
    %             indicator = ones( 1, 1, 1, number_of_scales );            
    %             
    %             gaussian_kernels = permute( exp( - kernel_sigma_domains .^ 2 / 2 ), ...
    %                                         [ 3, 4, 5, 1, 2 ]                       );
    %                             
    %             scale_mesh_1D = permute( kernel_scale_domain, [ 2, 3, 4, 1 ]);
    %              
    %        size_conv_energy_conv_kernel = squeeze( sum( scale_mesh_1D .* energy_chunk_4D .* gaussian_kernels, 4 ));
    %                  energy_conv_kernel = squeeze( sum(                  energy_chunk_4D .* gaussian_kernels, 4 ));
    %               indicator_conv_kernel = squeeze( sum(        indicator                 .* gaussian_kernels, 4 ));
    %               
    %               indicator_conv_kernel = permute( indicator_conv_kernel, [ 2, 3, 4, 1 ]);
    %               
    %             expected_scales_4D = size_conv_energy_conv_kernel ./ energy_conv_kernel ;
    % 
    %             expected_scales_4D( isnan( expected_scales_4D )) = 0 ;                                   
    %                                                
    %             energy_chunk_4D_4_smoothed = energy_conv_kernel ./ indicator_conv_kernel ;
    %             
    %             [ energy_chunk, energy_chunk_scale_indices ] = min( energy_chunk_4D_4_smoothed, [ ], 4 );
    % 
    % %             energy_chunk_scale_indices                                                              ...
    % %                 = reshape( round( expected_scales_4D( ( 1 : number_of_voxels_in_chunk )'            ...
    % %                                                      +      number_of_voxels_in_chunk               ...
    % %                                                        * ( energy_chunk_scale_indices( : ) - 1 ))), ...
    % %                            writing_counts( 1 : 3 )                                                  );
    %                         
    % %             % uncomment to save inter-scale information SAM 4/3/19             
    %             energy_chunk_scale_indices                                                              ...
    %                        = reshape( expected_scales_4D( ( 1 : number_of_voxels_in_chunk )'            ...
    %                                                      +      number_of_voxels_in_chunk               ...
    %                                                         * ( energy_chunk_scale_indices( : ) - 1 )), ...
    %                                   writing_counts( 1 : 3 )                                           );
    %                         
    %     end

        [ energy_chunk_min, energy_chunk_scale_min ] = min( energy_chunk_4D, [ ], 4 );

    %     energy_chunk_4D( energy_chunk_4D( : ) > 0 ) = 0 ;
    % 
    %     energy_chunk_sum = sum( energy_chunk_4D, 4 );            
    % %     energy_chunk_mean = mean( energy_chunk_4D, 4 );            
    % 
    %     energy_chunk_scale_sum                                                       ...
    %                          = sum(    energy_chunk_4D                                   ...
    %                                 .* permute( 1 : number_of_scales, [ 1, 3, 4, 2 ]),   ...
    %                                 4                                                 )  ...
    %                            ./ energy_chunk_sum                                       ;  
    % 
    %     energy_chunk_scale_sum( isnan( energy_chunk_scale_sum )) = 1 ;
    % 
    %     energy_chunk_sum = energy_chunk_sum ./ scales_per_octave ; % SAM 1/11/22, assuming it matches good up to +/- volume octave / 2 (15 % of radius octave)
    % %     energy_chunk_sum = energy_chunk_sum ./ number_of_scales ; % SAM before 1/11/22
    %     
    %     % min projection method goes with annular pulse signal and sum method goes with spherical signal
    %     energy_chunk               = spherical_to_annular_ratio * energy_chunk_mean              + ( 1 - spherical_to_annular_ratio ) * energy_chunk_min               ;
    %     energy_chunk_scale_indices = spherical_to_annular_ratio * energy_chunk_scale_indices_sum + ( 1 - spherical_to_annular_ratio ) * energy_chunk_scale_indices_min ;
    % 
        energy_chunk               = energy_chunk_min               ; 
        energy_chunk_scale_indices = energy_chunk_scale_min ...; % SAM before 1/11/22
                                   + sum( octave_at_scales < current_octave ); % accounting for previous octaves after switching order of nested FOR loops % SAM 6/6/22

    %     % min projection method goes with annular pulse and 3D-gaussian signals, and sum method goes with spherical pulse signal
    %     energy_chunk               = gaussian_to_ideal_ratio ^ 2 * energy_chunk_min       + ( 1 - gaussian_to_ideal_ratio ^ 2 ) * ( spherical_to_annular_ratio * energy_chunk_sum       + ( 1 - spherical_to_annular_ratio ) * energy_chunk_min       );
    %     energy_chunk_scale_indices = gaussian_to_ideal_ratio ^ 2 * energy_chunk_scale_min + ( 1 - gaussian_to_ideal_ratio ^ 2 ) * ( spherical_to_annular_ratio * energy_chunk_scale_sum + ( 1 - spherical_to_annular_ratio ) * energy_chunk_scale_min );

    %     smoothing_sizes_temp = false ;
    %     
    %     if smoothing_sizes_temp
    %         
    %         % make translatable strels to calculate gaussians within        
    %         linear_strel   = cell( number_of_scales, 1 );
    %         gaussian_strel = cell( number_of_scales, 1 );
    %         
    %         for scale_index = 1 : number_of_scales
    %                         
    %             strel_radii = round( 1 * pixels_per_radius_range( scale_index, : ));            
    %             
    %             linear_strel{ scale_index } =                                    ( - strel_radii( 1 ) : strel_radii( 1 ))' + ( - strel_radii( 2 ) : strel_radii( 2 )) * writing_counts( 1 );
    % 
    %             linear_strel{ scale_index } = linear_strel{ scale_index }( : ) + ( - strel_radii( 3 ) : strel_radii( 3 )) * prod( writing_counts( 1 : 2 ));
    % 
    %             linear_strel{ scale_index } = linear_strel{ scale_index }( : );
    %             
    %             [ mesh_Y, mesh_X, mesh_Z ] = ndgrid( - strel_radii( 1 ) : strel_radii( 1 ), ...
    %                                                  - strel_radii( 2 ) : strel_radii( 2 ), ...
    %                                                  - strel_radii( 3 ) : strel_radii( 3 )  );
    % 
    %             mesh_YXZ = cat( 4, mesh_Y, mesh_X, mesh_Z );
    %                                     
    %             gaussian_strel{ scale_index } = exp( - sum( mesh_YXZ .^ 2, 4 ) / lumen_radius_in_microns_range( scale_index ) ^ 2 ) / lumen_radius_in_microns_range( scale_index ) ;
    %             
    %             voxels_in_strel = sum( mesh_YXZ .^ 2, 4 ) / lumen_radius_in_microns_range( scale_index ) ^ 2 <= 1 ;
    %             
    %               linear_strel{ scale_index } =   linear_strel{ scale_index }( voxels_in_strel );
    %             gaussian_strel{ scale_index } = gaussian_strel{ scale_index }( voxels_in_strel );
    %             
    %         end
    %         
    %         voxel_index_range = 1 : number_of_voxels_in_chunk ;
    %         
    %         % only do the FOR loop for negative energy voxels that make the cutoff.  The other voxels
    %         % will have weights and sizes only generated from themselves (as if they had gaussian that
    %         % only rendered on their single voxel).
    %         voxel_index_range( energy_chunk > quantile( energy_chunk( energy_chunk < 0 ), 0.01 )) = [ ];
    %         
    % %                  weighting_image = ones( writing_counts( 1 : 3 )) ;
    %                  weighting_image =                          1 ./ lumen_radius_in_microns_range( round( energy_chunk_scale_indices )) ;
    %         size_and_weighting_image = energy_chunk_scale_indices ./ lumen_radius_in_microns_range( round( energy_chunk_scale_indices )) ;
    %         
    %                  weighting_image( voxel_index_range ) = 0 ;
    %         size_and_weighting_image( voxel_index_range ) = 0 ;
    %         
    % %         voxel_subscripts = ones( 1, 3 );
    %         
    %         for  voxel_index = voxel_index_range
    %                     
    %             scale_index_2 = round( energy_chunk_scale_indices( voxel_index ));
    %             
    %             % accept wraparound artifact for expedited execution
    %             valid_voxels = voxel_index + linear_strel{ scale_index_2 } >= 1 & voxel_index + linear_strel{ scale_index_2 } <= number_of_voxels_in_chunk ;
    %             
    %             current_strel = voxel_index + linear_strel{ scale_index_2 }( valid_voxels );
    %             
    %                      weighting_image( current_strel ) =          weighting_image( current_strel ) +                                             gaussian_strel{ scale_index_2 }( valid_voxels );
    %             size_and_weighting_image( current_strel ) = size_and_weighting_image( current_strel ) + energy_chunk_scale_indices( voxel_index ) * gaussian_strel{ scale_index_2 }( valid_voxels );
    %             
    % %             voxel_subscripts( 1 ) = voxel_subscripts( 1 ) + 1 ;
    % %             
    % %             if voxel_subscripts( 1 ) > size_of_image( 1 )
    % %                 
    % %                 voxel_subscripts( 1 ) = 1 ;
    % %                 
    % %                 voxel_subscripts( 2 ) = voxel_subscripts( 2 ) + 1 ;
    % %                 
    % %                 if voxel_subscripts( 2 ) > size_of_image( 2 )
    % %                     
    % %                     voxel_subscripts( 2 ) = 1 ;
    % %                     
    % %                     voxel_subscripts( 3 ) = voxel_subscripts( 3 ) + 1 ;
    % %                     
    % %                 end
    % %             end
    %         end % FOR voxel
    %         
    %         energy_chunk_scale_indices = size_and_weighting_image ./ weighting_image ;
    %         
    %     end

    %     smoothing_sizes = false ;    
    %     
    %     if smoothing_sizes
    % 
    %         % !!!! this has an edge effect at the chunk border
    %         linear_strel = ( -1 : 1 )' + ( -1 : 1 ) * writing_counts( 1 );
    %         
    %         linear_strel = linear_strel( : ) + ( -1 : 1 ) * prod( writing_counts( 1 : 2 ));
    %         
    %         linear_strel = linear_strel( : );
    %         
    %         linear_inner_image = ( 2 : writing_counts( 1 ) - 1 )' + ( 1 : writing_counts( 2 ) - 2 ) * writing_counts( 1 );
    %         
    %         linear_inner_image = linear_inner_image( : ) + ( 1 : writing_counts( 3 ) - 2 ) * prod( writing_counts( 1 : 2 ));
    %         
    %         linear_inner_image = linear_inner_image( : )' ;
    %         
    %         [ ~, min_energy_linear_strel_index ] = min( energy_chunk( linear_inner_image + linear_strel ));
    %         
    %         energy_chunk_scale_indices( linear_inner_image ) = energy_chunk_scale_indices( linear_inner_image + linear_strel( min_energy_linear_strel_index )');
    %         
    %     end % IF smoothing sizes

    %     is_adjusting_size = true ;
    %     
    %     if is_adjusting_size
    %         
    %         microns_per_sigma_PSF = pixels_per_sigma_PSF .* microns_per_voxel ;
    %        
    %         % subtracting the size (averaged over x y z dimensions) from the apparent radius of vessel
    %         % at each scale. May not be smaller than the PSF in any dimension
    %         lumen_radius_in_microns_range_adjusted =  geomean( max( 0, lumen_radius_in_microns_range .^ 2 - microns_per_sigma_PSF .^ 2 ) .^ 0.5, 2 );
    %         
    %         % looking up scale indices from the unadjusted size LUT using the adjusted sizes at each
    %         % image location
    %         energy_chunk_scale_indices = interp1( log( lumen_radius_in_microns_range ), ...
    %                                         1 : numel( lumen_radius_in_microns_range ), ...      
    %                                               log( lumen_radius_in_microns_range_adjusted( energy_chunk_scale_indices )), ...
    %                                               'linear', 1 );
    %         
    %     end


    %     force_real_size_energies = true ;
    % %     force_real_size_energies = false ; % SAM 9/3/19
    %     
    %     if force_real_size_energies
    %         
    %         energy_chunk( : ) = energy_chunk_4D(( 1 : prod( writing_counts( 1 : 3 ))) + round( energy_chunk_scale_indices( : )' - 1 ) * prod( writing_counts( 1 : 3 )));
    %         
    %     end % IF force_real_size_energies

        energy_chunk( energy_chunk >= 0 ) = 0 ;

    %     energy_chunk( energy_chunk_scale_indices ==         1        ) = 0 ;
    %     energy_chunk( energy_chunk_scale_indices == number_of_scales ) = 0 ; # SAM 1*13*22, this resulted in a worse output on noise sensitibity study

        index_and_energy = cat( 4, energy_chunk_scale_indices, energy_chunk );

        % Write h5 files to individual chunk files (to be combined later)   
        mat2h5( energy_chunk_file, index_and_energy );

    end % chunk PARFOR

    %% add to the tiled master file
    % combine the chunk files into a master file outside of a parfor to avoid simulataneous writing
    % issues
    for chunk_index = chunk_index_range

        [ y_chunk_index, ...
          x_chunk_index, ...
          z_chunk_index ] = ind2sub( chunk_lattice_dimensions, chunk_index );

        energy_chunk_name = [ int2str( chunk_index ), ' of ', int2str( number_of_chunks )];

        writing_starts = [ y_writing_starts( y_chunk_index ), ...
                           x_writing_starts( x_chunk_index ), ...
                           z_writing_starts( z_chunk_index ), ...
                           1                                  ];

        writing_counts = [ y_writing_counts( y_chunk_index ), ...
                           x_writing_counts( x_chunk_index ), ...
                           z_writing_counts( z_chunk_index ), ...
                           2                                  ];

        energy_chunk_file = [ chunk_directory, energy_chunk_name ];
        
        % only overwrite energy and scale images where the energy is better
        energy_chunk_new = h52mat( energy_chunk_file                                             );
        energy_chunk_old = h52mat( energy_file,                   writing_starts, writing_counts );
        
                                      is_better_voxel_new    ...
                            = energy_chunk_new( :, :, :, 2 ) ...
                            < energy_chunk_old( :, :, :, 2 ) ;
                                                voxels_to_overwrite                                      ...
                            =   find( is_better_voxel_new )      +     [ 0, prod( writing_counts( 1 : 3 ))];
                  
                              energy_chunk_old( voxels_to_overwrite ) ...
                            = energy_chunk_new( voxels_to_overwrite ) ;

                           mat2h5( energy_file, energy_chunk_old, writing_starts, writing_counts );

    end % chunk PARFOR

    % Remove Chunk File directories that have been tiled and will not be referenced again
    try rmdir( chunk_directory, 's' ), catch, end
    
end % octave FOR

end % FUNCTION