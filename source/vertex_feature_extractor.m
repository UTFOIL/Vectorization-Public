function [ laplacian, vesselness, energy, gradient_1, gradient_2, gradient_3, log_radius, is_cropped, axial_comp, blobness ] = vertex_feature_extractor( vertex_space_subscripts, vertex_scale_subscripts, path_to_original_data, lumen_radius_in_microns_range, microns_per_voxel, gaussian_to_ideal_ratio, spherical_to_annular_ratio, pixels_per_sigma_PSF )
%% vertex_info_extractor_V2
% This function extracts features for the vertex object using the original image for the purposes of
% training, testing, and using a machine learner to replace the manual vertex curation step.
%
% features:
%
%   curvature: laplacian (CNR), shape factors
%
%   gradient: (3-space basis chosen as 1st, 2nd, and 3rd principal curvature directions
%
% %   intensity: CNR, globally-normalized intensity
%
% normalization: all features will be unitless, and most will be locally normalized
% 
% downsampling: for larger objects over 10,000 voxels
% SAM 5/30/22
%
% NOTE: this function uses x,y,z<->1,2,3 instead of y,x,z<->1,2,3 found elsewhere in SLAVV software

vertex_space_subscripts = double( vertex_space_subscripts );

%% parameters
best_resolution_allowed = 1 / 4.5 ;

number_of_scales = length( lumen_radius_in_microns_range );

scales_per_octave = log(                                    2                                   ) ...
                  / log( lumen_radius_in_microns_range( 2 ) / lumen_radius_in_microns_range( 1 )) ...
                  / 3 ; % divide by three for the volume interperetation of octave

scale_subscripts_range = 1 : number_of_scales ;
     
original_file_info = h5info( path_to_original_data );

size_of_image = original_file_info.Datasets.Dataspace.Size ;

lumen_radius_in_voxels_range         = lumen_radius_in_microns_range    ./ microns_per_voxel ;

zero_vector = zeros( size( vertex_scale_subscripts ));

vertex_linear_indices = zero_vector ;

 laplacian = zero_vector ;
vesselness = zero_vector ;
  blobness = zero_vector ;
    energy = zero_vector ;
gradient_1 = zero_vector ;
gradient_2 = zero_vector ;
gradient_3 = zero_vector ;
log_radius = zero_vector ;
is_cropped = zero_vector ;
axial_comp = zero_vector ;

%% scale FOR (downsampling)
previous_octave = 0 ;

for s_subscript = scale_subscripts_range

    current_octave = ceil( s_subscript / scales_per_octave / 3 );
    
    if current_octave > previous_octave

        previous_octave = current_octave ;
        previous_scale  =   s_subscript  ;

        largest_scale_at_current_octave ...
                                    = min( number_of_scales, floor( current_octave * scales_per_octave * 3 ));

        resolutions_at_scale = min( microns_per_voxel / lumen_radius_in_microns_range( s_subscript ), best_resolution_allowed * ones( size( microns_per_voxel ))); % unitless um/um

        resolution_factors = round( best_resolution_allowed ./ resolutions_at_scale );

        largest_pixels_per_radius_at_octave = lumen_radius_in_voxels_range( largest_scale_at_current_octave, : ) ./ resolution_factors ;
             pixels_per_sigma_PSF_at_octave =                        pixels_per_sigma_PSF                        ./ resolution_factors ;        

        % Read original image and interpolate it to downsample at higher octaves
        image = double( h52mat( path_to_original_data,                 ...
                                  [ 1, 1, 1 ],                         ...
                                    1 +    floor(( size_of_image - 1 ) ...
                                        ./ resolution_factors ),       ...
                                           resolution_factors          ));

        % vertex strel generation
        size_of_image_at_octave      = size( image )     ;

        cumprod_size_of_image_at_octave = cumprod( size_of_image_at_octave );
        
        strel_apothems = round(  sqrt(( 1 - gaussian_to_ideal_ratio ^ 2 )) * largest_pixels_per_radius_at_octave               ...
                               + 3 * (      gaussian_to_ideal_ratio ^ 2    * largest_pixels_per_radius_at_octave .^ 2          ...
                                           +                                      pixels_per_sigma_PSF_at_octave .^ 2 ) .^ 0.5 ); % factor of 3 because the strel may contain a Gaussian with sigma equal to the vertex radius. z-score of 3 is good enough

        window_linear_template = calculate_linear_window( strel_apothems, cumprod_size_of_image_at_octave );
        
        % matched kernel generation

        microns_per_voxel_at_cctave = microns_per_voxel .* resolution_factors ;

        vertices_at_octave = find(   vertex_scale_subscripts >= previous_scale                   ...
                                   & vertex_scale_subscripts <=  largest_scale_at_current_octave );
        
                      vertex_space_subscripts( vertices_at_octave, : )       ...
        = round((     vertex_space_subscripts( vertices_at_octave, : ) - 1 ) ...
                  ./                 resolution_factors                + 1   );

        vertex_linear_indices(       vertices_at_octave ) ...
        =         sub2ind(      size_of_image_at_octave,  ...
           vertex_space_subscripts(  vertices_at_octave, 1 ), ...
           vertex_space_subscripts(  vertices_at_octave, 2 ), ...
           vertex_space_subscripts(  vertices_at_octave, 3 )  );
       
        dx = 1 ;
        dy = cumprod_size_of_image_at_octave( 1 );
        dz = cumprod_size_of_image_at_octave( 2 );

    end
    
    size_of_template = size( window_linear_template );
    
    [ matching_kernel_template, derivative_weights_from_blurring ] = calculate_matching_kernel( lumen_radius_in_microns_range( s_subscript ), size_of_template, microns_per_voxel_at_cctave, pixels_per_sigma_PSF_at_octave, spherical_to_annular_ratio, gaussian_to_ideal_ratio );
    
    vertices_at_scale = find( vertex_scale_subscripts == s_subscript )';
      
    for vertex_index = vertices_at_scale

        matching_kernel = matching_kernel_template ;        
        
        vertex_linear_index = vertex_linear_indices( vertex_index );
        
        vertex_window = window_linear_template ...
                      + vertex_linear_index    ;
                  
        % Cropping image by two voxels in each dimentions in anticipation of taking derivatives
            % Spatial derivaties add one on each side in each dimension
        lower_inset = max( 2 + strel_apothems - vertex_space_subscripts( vertex_index, : )                          , 0 ); 
        upper_inset = max( 1 + strel_apothems + vertex_space_subscripts( vertex_index, : ) - size_of_image_at_octave, 0 );

        matching_kernel = matching_kernel( 1 + lower_inset( 1 ) : end - upper_inset( 1 ), ...
                                           1 + lower_inset( 2 ) : end - upper_inset( 2 ), ...
                                           1 + lower_inset( 3 ) : end - upper_inset( 3 )  );
        
          vertex_window =   vertex_window( 1 + lower_inset( 1 ) : end - upper_inset( 1 ), ...
                                           1 + lower_inset( 2 ) : end - upper_inset( 2 ), ...
                                           1 + lower_inset( 3 ) : end - upper_inset( 3 )  );

        % normalization for differences in intesntiy variation (noise + signal) in the matched original kernel
        % normalize matching kernel in L1
        matching_kernel = matching_kernel                / sum( matching_kernel,                'all' )        ;
        % baseline correction (weighted mean of zero)
        window = image( vertex_window ); weighted_mean   = sum( matching_kernel .* window     , 'all' )        ;
        window =     window    ... % total variation (standard deviation)
                 - weighted_mean ;    standard_deviation = sum( matching_kernel .* window .^ 2, 'all' ) .^ 0.5 ;
        matching_kernel = matching_kernel ...
                                    / standard_deviation ;

        % extract image windows for computing the derivatives
        window_dx   = image( vertex_window + dx      )                                                                       - image( vertex_window - dx      );
        window_dy   = image( vertex_window + dy      )                                                                       - image( vertex_window - dy      );
        window_dz   = image( vertex_window + dz      )                                                                       - image( vertex_window - dz      );
        
        window_dxdx = image( vertex_window + dx      )               - 2 * image( vertex_window )                            + image( vertex_window - dx      );
        window_dydy = image( vertex_window + dy      )               - 2 * image( vertex_window )                            + image( vertex_window - dy      );
        window_dzdz = image( vertex_window + dz      )               - 2 * image( vertex_window )                            + image( vertex_window - dz      );
        
        window_dxdy = image( vertex_window + dx + dy ) + image( vertex_window - dx - dy ) - image( vertex_window + dx - dy ) - image( vertex_window - dx + dy );
        window_dydz = image( vertex_window + dy + dz ) + image( vertex_window - dy - dz ) - image( vertex_window + dy - dz ) - image( vertex_window - dy + dz );
        window_dzdx = image( vertex_window + dz + dx ) + image( vertex_window - dz - dx ) - image( vertex_window + dz - dx ) - image( vertex_window - dz + dx );
        
        % correct for difference in finite difference approximation length
        window_dx   = window_dx   / 2 ;
        window_dy   = window_dy   / 2 ;
        window_dz   = window_dz   / 2 ;
        
        window_dxdy = window_dxdy / 4 ;
        window_dydz = window_dydz / 4 ;
        window_dzdx = window_dzdx / 4 ;
        
        % image window (finite differences)
%         matched_original      = matching_kernel .* window      ;
        matched_original_dx   = sum( matching_kernel .* window_dx  , 'all' );
        matched_original_dy   = sum( matching_kernel .* window_dy  , 'all' );
        matched_original_dz   = sum( matching_kernel .* window_dz  , 'all' );
        matched_original_dxdx = sum( matching_kernel .* window_dxdx, 'all' );
        matched_original_dydy = sum( matching_kernel .* window_dydy, 'all' );
        matched_original_dzdz = sum( matching_kernel .* window_dzdz, 'all' );
        matched_original_dxdy = sum( matching_kernel .* window_dxdy, 'all' );
        matched_original_dydz = sum( matching_kernel .* window_dydz, 'all' );
        matched_original_dzdx = sum( matching_kernel .* window_dzdx, 'all' );
        
        % normalization for differences in real space [microns] between voxels        
            % curvature
        matched_original_dxdx = prod( derivative_weights_from_blurring([ 1, 1 ])) * matched_original_dxdx ;
        matched_original_dydy = prod( derivative_weights_from_blurring([ 2, 2 ])) * matched_original_dydy ;
        matched_original_dzdz = prod( derivative_weights_from_blurring([ 3, 3 ])) * matched_original_dzdz ;
        matched_original_dxdy = prod( derivative_weights_from_blurring([ 1, 2 ])) * matched_original_dxdy ;
        matched_original_dydz = prod( derivative_weights_from_blurring([ 2, 3 ])) * matched_original_dydz ;
        matched_original_dzdx = prod( derivative_weights_from_blurring([ 3, 1 ])) * matched_original_dzdx ;
        
            % gradient
        matched_original_dx   =       derivative_weights_from_blurring(    1    ) * matched_original_dx   ;
        matched_original_dy   =       derivative_weights_from_blurring(    2    ) * matched_original_dy   ;
        matched_original_dz   =       derivative_weights_from_blurring(    3    ) * matched_original_dz   ;
        
        % feature definitions:
            % principal curvatures:
        [ principal_directions,   ...
          principal_curvatures  ] ... % directions in COLUMNS, values in DIAGONAL
                       = eig([ matched_original_dxdx, matched_original_dxdy, matched_original_dzdx; ...
                               matched_original_dxdy, matched_original_dydy, matched_original_dydz; ...
                               matched_original_dzdx, matched_original_dydz, matched_original_dzdz  ]);
            
            % gradient (projections onto principal directions)
        principal_gradients( 1 ) = abs([ matched_original_dx, matched_original_dy, matched_original_dz ] * principal_directions( :, 1 ));
        principal_gradients( 2 ) = abs([ matched_original_dx, matched_original_dy, matched_original_dz ] * principal_directions( :, 2 ));
        principal_gradients( 3 ) = abs([ matched_original_dx, matched_original_dy, matched_original_dz ] * principal_directions( :, 3 ));

%             % gradient (projections onto principal directions and further normalization by curvatures
%         gradient_1 = abs([ matched_original_dx, matched_original_dy, matched_original_dz ] * principal_directions( :, 1 )) / principal_curvatures( 1 );
%         gradient_2 = abs([ matched_original_dx, matched_original_dy, matched_original_dz ] * principal_directions( :, 2 )) / principal_curvatures( 5 );
%         gradient_3 = abs([ matched_original_dx, matched_original_dy, matched_original_dz ] * principal_directions( :, 3 ))                            ; % ! last principal curvature will be near zero !
        
%         vesselness( vertex_index ) = max(   prod(                principal_curvatures([ 1, 5    ])     )        ...
%                                           /  sum(                principal_curvatures([ 1, 5    ]) / 2 ) ^ 2, 0 ); % inside [0,1]

%         vesselness( vertex_index ) =   max( prod(        principal_curvatures([ 1, 5    ])), 0 ) ^ 0.5   ...
%                                    / ( min(  sum(        principal_curvatures([ 1, 5, 9 ])), 0 ) / 2   ) ; % inside [0,1] (mostly, could be larger if 3rd component is positive) % SAM 10.17.22
        vesselness( vertex_index ) =   max( prod(        principal_curvatures([ 1, 5 ])), 0 ) ^ 0.5   ...
                                   / ( min(  sum(        principal_curvatures([ 1, 5 ])), 0 ) / 2 + eps ) ; % inside [0,1] % SAM 2.7.23

        blobness(  vertex_index ) =    max( prod(        principal_curvatures([ 1, 9 ])), 0 ) ^ 0.5   ...
                                   / ( min(  sum(        principal_curvatures([ 1, 9 ])), 0 ) / 2 + eps ) ; % inside [0,1] % SAM 2.7.23

        laplacian( vertex_index )  =         sum(        principal_curvatures([ 1, 5, 9 ])     )         ;
           energy( vertex_index )  = sum(                principal_curvatures([ 1, 5, 9 ])       ...
                                          .* exp( - (    principal_gradients                     ...
                                                      ./ principal_curvatures([ 1, 5, 9 ])) .^ 2 ));
        % consider summing only the first two principal components only for above two features

        gradient_1( vertex_index ) =                             principal_gradients( 1 ) ;
        gradient_2( vertex_index ) =                             principal_gradients( 2 ) ;
        gradient_3( vertex_index ) =                             principal_gradients( 3 ) ;

        axial_comp( vertex_index ) =                        abs( principal_directions( end )); % z component of third P.C.

        log_radius( vertex_index ) = log( lumen_radius_in_microns_range( vertex_scale_subscripts( vertex_index )) / 3 ); % normalize radius to 3 um, approx. size of capillary
        
        is_cropped( vertex_index ) = 1 - numel( matching_kernel ) / numel( matching_kernel_template );
        
    end % END vertex FOR    
end % END scale FOR


% intensity

end % END main FUNCTION

%% Auxillary functions
function [ matching_kernel_template, derivative_weights_from_blurring ] = calculate_matching_kernel( radius_of_lumen_in_microns, size_of_kernel, microns_per_pixel, pixels_per_sigma_PSF, spherical_to_annular_ratio, gaussian_to_ideal_ratio )

    [ y_pixel_freq_mesh, x_pixel_freq_mesh, z_pixel_freq_mesh ]                                                          ...
        = ndgrid([ 0 : ( size_of_kernel( 1 ) - 1 ) / 2, - ( size_of_kernel( 1 ) - 1 ) / 2 : - 1 ] / size_of_kernel( 1 ), ...
                 [ 0 : ( size_of_kernel( 2 ) - 1 ) / 2, - ( size_of_kernel( 2 ) - 1 ) / 2 : - 1 ] / size_of_kernel( 2 ), ...
                 [ 0 : ( size_of_kernel( 3 ) - 1 ) / 2, - ( size_of_kernel( 3 ) - 1 ) / 2 : - 1 ] / size_of_kernel( 3 )  ); 

    y_micron_freq_mesh = y_pixel_freq_mesh / microns_per_pixel( 1 );
    x_micron_freq_mesh = x_pixel_freq_mesh / microns_per_pixel( 2 );
    z_micron_freq_mesh = z_pixel_freq_mesh / microns_per_pixel( 3 );

    microns_per_sigma_PSF = pixels_per_sigma_PSF .* microns_per_pixel ;

         Gaussian_lengths         =       gaussian_to_ideal_ratio       * ( radius_of_lumen_in_microns ^ 2 + microns_per_sigma_PSF .^ 2 ) .^ 0.5 ; % SAM 12/22/21
    annular_pulse_lengths_squared = ( 1 - gaussian_to_ideal_ratio ^ 2 ) * ( radius_of_lumen_in_microns ^ 2 + microns_per_sigma_PSF .^ 2 )        ; % SAM 12/22/21

    sphere_pulse_lengths_squared = annular_pulse_lengths_squared ; % SAM 12/16/21

    y_radial_freq_mesh_squared = y_micron_freq_mesh .^ 2 * sphere_pulse_lengths_squared( 1 );
    x_radial_freq_mesh_squared = x_micron_freq_mesh .^ 2 * sphere_pulse_lengths_squared( 2 );
    z_radial_freq_mesh_squared = z_micron_freq_mesh .^ 2 * sphere_pulse_lengths_squared( 3 );

    radial_freq_mesh_sphere_pulse = (   y_radial_freq_mesh_squared          ...
                                      + x_radial_freq_mesh_squared          ...
                                      + z_radial_freq_mesh_squared ) .^ 0.5 ;

    radial_angular_freq_mesh = 2 * pi * radial_freq_mesh_sphere_pulse ;

    spherical_pulse_kernel_dft =  ( pi / 2 ./ radial_angular_freq_mesh ) .^ 0.5 ...
                               .* (   besselj( 2.5, radial_angular_freq_mesh )  ...
                                    + besselj( 0.5, radial_angular_freq_mesh )) ;

    spherical_pulse_kernel_dft( radial_angular_freq_mesh == 0 ) = 1 ;

    % ... spherical pulse addition end of section.

    y_radial_freq_mesh = y_micron_freq_mesh * Gaussian_lengths( 1 );
    x_radial_freq_mesh = x_micron_freq_mesh * Gaussian_lengths( 2 );
    z_radial_freq_mesh = z_micron_freq_mesh * Gaussian_lengths( 3 );

    %         gaussian_PSF_kernel_dft = ones( size( y_radial_freq_mesh ));

    radial_freq_mesh_gaussian     = (   y_radial_freq_mesh .^ 2          ...
                                      + x_radial_freq_mesh .^ 2          ...
                                      + z_radial_freq_mesh .^ 2 ) .^ 0.5 ;

    sigma_freq_mesh          =          radial_freq_mesh_gaussian     ;

    y_radial_freq_mesh_squared = y_micron_freq_mesh .^ 2 * annular_pulse_lengths_squared( 1 );
    x_radial_freq_mesh_squared = x_micron_freq_mesh .^ 2 * annular_pulse_lengths_squared( 2 );
    z_radial_freq_mesh_squared = z_micron_freq_mesh .^ 2 * annular_pulse_lengths_squared( 3 );

    radial_freq_mesh_sphere_pulse = (   y_radial_freq_mesh_squared          ...
                                      + x_radial_freq_mesh_squared          ...
                                      + z_radial_freq_mesh_squared ) .^ 0.5 ;        

    radial_angular_freq_mesh = 2 * pi * radial_freq_mesh_sphere_pulse ;        gaussian_kernel_dft = exp( - pi ^ 2 * 2 * sigma_freq_mesh .^ 2 );
    annular_pulse_kernel_dft = cos( radial_angular_freq_mesh );
    matching_kernel_dft = gaussian_kernel_dft .* (( 1 - spherical_to_annular_ratio ) * annular_pulse_kernel_dft + spherical_to_annular_ratio * spherical_pulse_kernel_dft );

%     matching_kernel_dft = matching_kernel_dft .* exp( -1i * y_pixel_freq_mesh * ( size_of_kernel( 1 ) - 1 ) / 2    ...
%                                                       -1i * x_pixel_freq_mesh * ( size_of_kernel( 2 ) - 1 ) / 2    ...
%                                                       -1i * z_pixel_freq_mesh * ( size_of_kernel( 3 ) - 1 ) / 2 ) .^ ( 2 * pi);

%     matching_kernel_template = ifftn( matching_kernel_dft, 'symmetric' );
    matching_kernel_template = fftshift( ifftn( matching_kernel_dft, 'symmetric' ));
    
    derivative_weights_from_blurring = Gaussian_lengths ./ microns_per_pixel ; % SAM 12/14/21
    
end

function window_linear_template = calculate_linear_window( strel_apothems, cumprod_size_of_image_at_octave )

    local_subscripts_range_1 = - strel_apothems( 1 ) : strel_apothems( 1 );
    local_subscripts_range_2 = - strel_apothems( 2 ) : strel_apothems( 2 );
    local_subscripts_range_3 = - strel_apothems( 3 ) : strel_apothems( 3 );

    window_linear_template = local_subscripts_range_1                                                                     ;
    window_linear_template = local_subscripts_range_2 * cumprod_size_of_image_at_octave( 1 ) + window_linear_template( : );
    window_linear_template = local_subscripts_range_3 * cumprod_size_of_image_at_octave( 2 ) + window_linear_template( : );
    
    window_linear_template = reshape( window_linear_template( : ), length( local_subscripts_range_1 ), ...
                                                                   length( local_subscripts_range_2 ), ...
                                                                   length( local_subscripts_range_3 )  );

end