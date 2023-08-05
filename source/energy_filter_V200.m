function [ energy_chunk ]                                                                             ...
         = energy_filter_V200( chunk_dft, matching_kernel_string, radius_of_lumen_in_microns,         ...
                               vessel_wall_thickness_in_microns, microns_per_pixel,                   ...
                               pixels_per_sigma_PSF, local_ranges, gaussian_to_ideal_ratio, spherical_to_annular_ratio, scales_per_octave )
%% energy_filter
% SAM 9/12/17
% V2 in which the convolutions are done in the frequency domain, SAM 11/12/17
%
% V3, in which the function is expecting the image matrix already Fourier
% transformed.  SAM 12/6/17
%
% V10 in which the PSF is deconvolved and the filters are analytically rendered in the Fourier
% domain before operating with the data (this is to improve the numerical stability of the filter
% construction and the deconvolution steps, which would be unstable at small blurring scales. SAM
% 3/31/18
%
% Function name changed from gaussian_filter_V10 to energy_filter with sum of the principal
% curvatures as the energy function % 5/5/18 SAM
%
% V140, energy = sum of least two principal curvatures
%
% V141, energy = least principal curvature
%
% V142 in which the gradient images are also computed, from which the direction of the gradient
% vector iamge is created, all of the principal curvatures are then projected onto the plane that is
% orthogonal to the gradient before summing them up. SAM 5/7/18
%
% V150 in which the local ranges are input to this function so that we may save time in trimming the
% edges of the chunk that will not be read before calculating energies (principal curvatures in a
% FOR loop) on those redundant image regions. SAM 5/7/18
%
% V151 in which the gradient vector is projected onto the three principal curvatures and the
% component in each is divided by the value of the principal curvature value in that direction.
% We then take the exponent of the negative absolute value of this ratio, and use this as the
% scaling factor for each principal curvature before summing them up.  SAM 5/8/18
%
% V152 in which the ratios of the gradient components to the principal curvatures is scaled by a
% factor.  This scaling factor is an input, the symmetry_ratio_factor. SAM 5/9/18
%
% V153 in which the symmetry ratio is squared instead of absolute valued. SAM 5/9/18
%
% V160 in whcih the symmetry ratio factor is 3D vector quantitied.  SAM 5/17/18
%
% V161, the symmetry ratio is back to absolute valued (before V153) SAM 5/18/18
%
% V162, in which the pixels per sigma is a variable in the smmetry ratio, (since V160 it has been
% pixels per sigma squared in there) SAM 5/21/18
%
% V190: The gaussian kernel is replace by a 3D sphere with radius of pixels_per_radius that is
% blurred by the PSF. SAM 8/2/18  This function has turned into a sandbox.  See V191 for this
% modification. SAM 8/9/18
%
% V191: the intended outcome of V190. See V190.  Also adding a string input to switch between
% different types of matching techniques (spherical, annular, gaussian).  Not all inputs will be
% used when you call this function.  The wall thickness only applies to the annular pulse kernel and
% the radius per sigma only applies to gaurrian kernel. SAM 9/10/18
%
% V200 in which the output is not cropped. SAM 12/10/18

is_only_computing_Laplacian = false ; % SAM 12/10/21
% is_only_computing_Laplacian = true ; % the rest of teh time SAM
% is_only_computing_Laplacian = false ; % 191202 SAM

if ~ is_only_computing_Laplacian
    
%     is_doing_eigenvalue_decomp = false ;
    is_doing_eigenvalue_decomp = true ; % SAM 12/10/21

%     symmetry_ratio_factor = 1 ; 
        symmetry_ratio_factor = 2 ^ 0.5 ; % SAM 12/20/21
%     symmetry_ratio_factor = 2 ; % SAM before 12/20/21
%     symmetry_ratio_factor = 4 ; 

end

% scales_per_octave_radius = scales_per_octave * 3 ;

% mesh generation for building the energy filter in the Fourier domain
[ size_of_chunk_dft( 1 ), size_of_chunk_dft( 2 ), size_of_chunk_dft( 3 )] = size( chunk_dft );

% numel_chunk = numel( chunk_dft );

[ y_pixel_freq_mesh, x_pixel_freq_mesh, z_pixel_freq_mesh ]                                                       ...
    = ndgrid([ 0 : size_of_chunk_dft( 1 ) / 2 - 1, - size_of_chunk_dft( 1 ) / 2 : - 1 ] / size_of_chunk_dft( 1 ), ...
             [ 0 : size_of_chunk_dft( 2 ) / 2 - 1, - size_of_chunk_dft( 2 ) / 2 : - 1 ] / size_of_chunk_dft( 2 ), ...
             [ 0 : size_of_chunk_dft( 3 ) / 2 - 1, - size_of_chunk_dft( 3 ) / 2 : - 1 ] / size_of_chunk_dft( 3 )  );

% MATLAB documentation on the FFT suggests a mesh like the following, but the preceding is the one
% that works in practice. SAM 8/10/18
%
% [ y_pixel_freq_mesh,                                                                        ...
%   x_pixel_freq_mesh,                                                                        ... 
%   z_pixel_freq_mesh ] = ndgrid(( 0 : size_of_chunk_dft( 1 ) - 1 ) / size_of_chunk_dft( 1 ), ...
%                                ( 0 : size_of_chunk_dft( 2 ) - 1 ) / size_of_chunk_dft( 2 ), ...
%                                ( 0 : size_of_chunk_dft( 3 ) - 1 ) / size_of_chunk_dft( 3 )  );
         
y_micron_freq_mesh = y_pixel_freq_mesh / microns_per_pixel( 1 );
x_micron_freq_mesh = x_pixel_freq_mesh / microns_per_pixel( 2 );
z_micron_freq_mesh = z_pixel_freq_mesh / microns_per_pixel( 3 );

% y_radial_freq_mesh = y_micron_freq_mesh * radius_of_lumen_in_microns ;
% x_radial_freq_mesh = x_micron_freq_mesh * radius_of_lumen_in_microns ;
% z_radial_freq_mesh = z_micron_freq_mesh * radius_of_lumen_in_microns ;
         
% y_micron_freq_mesh_squared = y_micron_freq_mesh .^ 2 ;
% x_micron_freq_mesh_squared = x_micron_freq_mesh .^ 2 ;
% z_micron_freq_mesh_squared = z_micron_freq_mesh .^ 2 ;
%          
% micron_freq_mesh = (   y_micron_freq_mesh_squared          ...
%                      + x_micron_freq_mesh_squared          ...
%                      + z_micron_freq_mesh_squared ) .^ 0.5 ;
%                      
% radial_freq_mesh = micron_freq_mesh * radius_of_lumen_in_microns ;
%                  
% sigma_PSF_freq_mesh_squared = ( y_pixel_freq_mesh * pixels_per_sigma_PSF( 1 )) .^ 2 ...
%                             + ( x_pixel_freq_mesh * pixels_per_sigma_PSF( 2 )) .^ 2 ...
%                             + ( z_pixel_freq_mesh * pixels_per_sigma_PSF( 3 )) .^ 2 ;
                        
% gaussian_PSF_kernel_dft = exp( - pi ^ 2 * 2 * sigma_PSF_freq_mesh_squared );


microns_per_sigma_PSF = pixels_per_sigma_PSF .* microns_per_pixel ;
                                                                  
% radius_of_lumen_in_voxels = radius_of_lumen_in_microns ./ microns_per_pixel ;

%% matching kernel
%     
%     case '3D gaussian'
%                 
% %         sigma_freq_mesh = radial_freq_mesh ;
%         
% %         % add PSF
% %         y_radial_freq_mesh = y_micron_freq_mesh * ( radius_of_lumen_in_microns ^ 2 +     microns_per_sigma_PSF( 1 ) ^ 2 ) ^ 0.5 ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * ( radius_of_lumen_in_microns ^ 2 +     microns_per_sigma_PSF( 2 ) ^ 2 ) ^ 0.5 ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * ( radius_of_lumen_in_microns ^ 2 +     microns_per_sigma_PSF( 3 ) ^ 2 ) ^ 0.5 ;
% %         
% %         radius_of_lumen_in_voxels = ( radius_of_lumen_in_microns ^ 2 + 2 * microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ;
%         
%         % stably deconvolve PSF
%         y_radial_freq_mesh_squared = y_micron_freq_mesh .^ 2 * max( radius_of_lumen_in_microns ^ 2 -     microns_per_sigma_PSF( 1 ) ^ 2, microns_per_sigma_PSF( 1 ) ^ 2 );
%         x_radial_freq_mesh_squared = x_micron_freq_mesh .^ 2 * max( radius_of_lumen_in_microns ^ 2 -     microns_per_sigma_PSF( 2 ) ^ 2, microns_per_sigma_PSF( 2 ) ^ 2 );
%         z_radial_freq_mesh_squared = z_micron_freq_mesh .^ 2 * max( radius_of_lumen_in_microns ^ 2 -     microns_per_sigma_PSF( 3 ) ^ 2, microns_per_sigma_PSF( 3 ) ^ 2 );
% 
% %         % stably deconvolve PSF
% %         y_radial_freq_mesh_squared = y_micron_freq_mesh .^ 2 * max( radius_of_lumen_in_microns ^ 2 -     microns_per_sigma_PSF( 1 ) ^ 2, max( microns_per_pixel( 1 ) / 2, microns_per_sigma_PSF( 1 )) ^ 2 );
% %         x_radial_freq_mesh_squared = x_micron_freq_mesh .^ 2 * max( radius_of_lumen_in_microns ^ 2 -     microns_per_sigma_PSF( 2 ) ^ 2, max( microns_per_pixel( 2 ) / 2, microns_per_sigma_PSF( 2 )) ^ 2 );
% %         z_radial_freq_mesh_squared = z_micron_freq_mesh .^ 2 * max( radius_of_lumen_in_microns ^ 2 -     microns_per_sigma_PSF( 3 ) ^ 2, max( microns_per_pixel( 3 ) / 2, microns_per_sigma_PSF( 3 )) ^ 2 );
% 
% %         % current best
% %         radius_of_lumen_in_voxels = max( radius_of_lumen_in_microns ^ 2, microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ;
% 
% % %         radius_of_lumen_in_voxels = max( radius_of_lumen_in_microns ^ 2, 2 * microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ;
% 
%         radius_of_lumen_in_voxels = max( radius_of_lumen_in_microns ^ 2 - microns_per_sigma_PSF .^ 2, microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ;
% 
% %         radius_of_lumen_in_voxels( 1 ) = real( radius_of_lumen_in_voxels_temp( 1 )) + imag( sum( radius_of_lumen_in_voxels_temp([ 2, 3 ]))) / 2 ;
% %         radius_of_lumen_in_voxels( 2 ) = real( radius_of_lumen_in_voxels_temp( 2 )) + imag( sum( radius_of_lumen_in_voxels_temp([ 3, 1 ]))) / 2 ;
% %         radius_of_lumen_in_voxels( 3 ) = real( radius_of_lumen_in_voxels_temp( 3 )) + imag( sum( radius_of_lumen_in_voxels_temp([ 1, 2 ]))) / 2 ;
% 
% %         % deconvolve PSF
% %         y_radial_freq_mesh_squared = y_micron_freq_mesh .^ 2 * ( radius_of_lumen_in_microns ^ 2 -     microns_per_sigma_PSF( 1 ) ^ 2 );
% %         x_radial_freq_mesh_squared = x_micron_freq_mesh .^ 2 * ( radius_of_lumen_in_microns ^ 2 -     microns_per_sigma_PSF( 2 ) ^ 2 );
% %         z_radial_freq_mesh_squared = z_micron_freq_mesh .^ 2 * ( radius_of_lumen_in_microns ^ 2 -     microns_per_sigma_PSF( 3 ) ^ 2 );
% 
% %         % ignore PSF
% %         y_radial_freq_mesh = y_micron_freq_mesh * radius_of_lumen_in_microns ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * radius_of_lumen_in_microns ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * radius_of_lumen_in_microns ;
% 
% %         radius_of_lumen_in_voxels = ( radius_of_lumen_in_microns ^ 2 + microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ;
% 
% %         % never search below PSF
% %         y_radial_freq_mesh_squared = y_micron_freq_mesh .^ 2 * max( radius_of_lumen_in_microns ^ 2, microns_per_sigma_PSF( 1 ) ^ 2 );
% %         x_radial_freq_mesh_squared = x_micron_freq_mesh .^ 2 * max( radius_of_lumen_in_microns ^ 2, microns_per_sigma_PSF( 2 ) ^ 2 );
% %         z_radial_freq_mesh_squared = z_micron_freq_mesh .^ 2 * max( radius_of_lumen_in_microns ^ 2, microns_per_sigma_PSF( 3 ) ^ 2 );
% % 
% %         radius_of_lumen_in_voxels = max( radius_of_lumen_in_microns ^ 2, microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ;
%         
% %         gaussian_PSF_kernel_dft = ones( size( y_radial_freq_mesh ));
%         
% %         squared_radial_freq_mesh_gaussian     = (   y_radial_freq_mesh .^ 2  ...
% %                                                   + x_radial_freq_mesh .^ 2  ...
% %                                                   + z_radial_freq_mesh .^ 2 );
% 
%         squared_radial_freq_mesh_gaussian     = (   y_radial_freq_mesh_squared  ...
%                                                   + x_radial_freq_mesh_squared  ...
%                                                   + z_radial_freq_mesh_squared );
%                                       
%         matching_kernel_dft = exp( - pi ^ 2 * 2 * squared_radial_freq_mesh_gaussian );
%         
% %         y_radial_freq_mesh = y_micron_freq_mesh * ( radius_of_lumen_in_microns ^ 2 + 2 * microns_per_sigma_PSF( 1 ) ^ 2 ) ^ 0.5 ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * ( radius_of_lumen_in_microns ^ 2 + 2 * microns_per_sigma_PSF( 2 ) ^ 2 ) ^ 0.5 ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * ( radius_of_lumen_in_microns ^ 2 + 2 * microns_per_sigma_PSF( 3 ) ^ 2 ) ^ 0.5 ;
% 
% %         % ignore PSF
% %         y_radial_freq_mesh = y_micron_freq_mesh * radius_of_lumen_in_microns ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * radius_of_lumen_in_microns ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * radius_of_lumen_in_microns ;
% 
% %         y_radial_freq_mesh = y_micron_freq_mesh * ( radius_of_lumen_in_microns ^ 2 +     microns_per_sigma_PSF( 1 ) ^ 2 ) ^ 0.5 ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * ( radius_of_lumen_in_microns ^ 2 +     microns_per_sigma_PSF( 2 ) ^ 2 ) ^ 0.5 ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * ( radius_of_lumen_in_microns ^ 2 +     microns_per_sigma_PSF( 3 ) ^ 2 ) ^ 0.5 ;
% 
% %         % add PSF
% %         radius_of_lumen_in_voxels = ( radius_of_lumen_in_microns ^ 2 + microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ;
%     
%     case 'spherical pulse'
%         
%         radial_angular_freq_mesh = 2 * pi * radial_freq_mesh  ;
% 
%         % this one will have a total sum of 1
%         matching_kernel_dft = 3 ./ ( radial_angular_freq_mesh .^ 2 )                           ...
%                         .* (   sin( radial_angular_freq_mesh ) ./ ( radial_angular_freq_mesh ) ...
%                              - cos( radial_angular_freq_mesh )                                 );
% 
%         matching_kernel_dft( 1 ) = 1 ;
% 
% %         % uncomment to have value inside sphere of approximately 1 (instead of sum to 1)
% %         spherical_volume_in_cubic_pxls = 4 / 3 * pi * radius_of_lumen_in_microns ^ 3 / prod( microns_per_pixel );   
% % 
% %         matching_kernel_dft = matching_kernel_dft * spherical_volume_in_cubic_pxls ;
%         
%     case '3D gaussian conv spherical pulse'
%         
% %         radius_of_lumen_in_voxels = radius_of_lumen_in_microns ./ microns_per_pixel ;
% 
% %         A = 0.95 ; % percent of radius made up of Gaussian to that made of spherical pulse        
% % 
% %         A = 0.8 ; % percent of radius made up of Gaussian to that made of spherical pulse        
% %         A = 0.5 ; % percent of radius made up of Gaussian to that made of spherical pulse
% %         A = 0.2 ; % percent of radius made up of Gaussian to that made of spherical pulse
% 
% 
% %         maximum_Gaussian_length_in_voxels = 2 * [ 1, 1, 1 ];
% %         maximum_Gaussian_length_in_voxels = Inf * [ 1, 1, 1 ];
% %         maximum_Gaussian_length_in_voxels = [ 0, 0, 0 ];
% %         maximum_Gaussian_length_in_voxels = 3 * [ 1, 1, 1 ];
% %         maximum_Gaussian_length_in_voxels = [ 1, 1, 1 ];
% %         maximum_Gaussian_length_in_voxels = 5 * [ 1, 1, 1 ];
% %         maximum_Gaussian_length_in_voxels = 4 * [ 1, 1, 1 ];
% 
% %         maximum_Gaussian_length_in_voxels = [ 3, 3, 2 ];
%         
% %         maximum_Gaussian_length_in_voxels = 2 * max( microns_per_pixel ) ./ microns_per_pixel ;
% %         maximum_Gaussian_length_in_voxels = max( microns_per_pixel ) ./ microns_per_pixel ;
%         
% %         maximum_Gaussian_length = max( 1.5 * microns_per_pixel );
% 
% %         maximum_Gaussian_length = 0 ;
%         
% %         maximum_Gaussian_length = microns_per_pixel ;
% 
% %         % current best
% %         maximum_Gaussian_length = max( 1.5 * microns_per_pixel + microns_per_sigma_PSF );
%         
% %                 maximum_Gaussian_length_in_voxels = 2 ;
%         
%         % !!! add a minimum Gaussian_length_in_voxels calculated as (requires extra input to fxn) !!
% %         minimum_Gaussian_length       = (   radius_of_lumen_in_microns( current  )  ...
% %                                           - radius_of_lumen_in_microns( previous )) ;
% % 
% %         maximum_Gaussian_length_in_voxels = max( maximum_Gaussian_length_in_voxels,           ...
% %                                                  minimum_Gaussian_length ./ microns_per_pixel );
% 
% %         Gaussian_lengths = max(( microns_per_pixel .^ 2 - microns_per_sigma_PSF .^ 2 ) .^ 0.5, [ 0, 0, 0 ]);
%         
% %         Gaussian_lengths =   max( radius_of_lumen_in_voxels,         ...
% %                                                  minimum_Gaussian_length_in_microns  );    
% 
% %         Gaussian_lengths =   min( radius_of_lumen_in_microns,     ...
% %                                                  Gaussian_lengths );
%         
% %         Gaussian_lengths = min( radius_of_lumen_in_microns,                             ...
% %                                                maximum_Gaussian_length_in_voxels .* microns_per_pixel );
% 
% %         Gaussian_length = min( radius_of_lumen_in_microns, maximum_Gaussian_length );
% 
%         radius_of_lumen_at_next_scale = radius_of_lumen_in_microns * 2 ^ ( 1 / scales_per_octave_radius );
%         
% %         Gaussian_lengths = max( microns_per_sigma_PSF, ( min( gaussian_to_ideal_ratio ^ 2 * ( radius_of_lumen_at_next_scale ^ 2 - radius_of_lumen_in_microns ^ 2 ), radius_of_lumen_in_microns ^ 2 ) - microns_per_sigma_PSF .^ 2 ) .^ 0.5 );
%         Gaussian_lengths = max( microns_per_sigma_PSF, (( gaussian_to_ideal_ratio * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 )) ^ 2 + microns_per_sigma_PSF .^ 2 ) .^ 0.5 );
% 
% %         sphere_pulse_lengths = radius_of_lumen_in_microns - Gaussian_lengths ;
% 
%         % assuming that the squared length of the combined kernel is the sum of the squared lengths
%         % of the two kernels being convolved.
% %         sphere_pulse_lengths_squared = max(( radius_of_lumen_in_microns ^ 2 - Gaussian_lengths .^ 2 - microns_per_sigma_PSF .^ 2 ), 0 );        
%         sphere_pulse_lengths_squared = max(( radius_of_lumen_in_microns ^ 2 - Gaussian_lengths .^ 2 + microns_per_sigma_PSF .^ 2 ), 0 );        
% %         sphere_pulse_length = max(( radius_of_lumen_in_microns ^ 2 - Gaussian_length ^ 2 ) ^ 0.5, 0 );
%                                
% %         y_radial_freq_mesh = y_micron_freq_mesh * sphere_pulse_length ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * sphere_pulse_length ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * sphere_pulse_length ;
% 
%         y_radial_freq_mesh_squared = y_micron_freq_mesh .^ 2 * sphere_pulse_lengths_squared( 1 );
%         x_radial_freq_mesh_squared = x_micron_freq_mesh .^ 2 * sphere_pulse_lengths_squared( 2 );
%         z_radial_freq_mesh_squared = z_micron_freq_mesh .^ 2 * sphere_pulse_lengths_squared( 3 );
% 
% %         radial_freq_mesh_sphere_pulse = (   y_radial_freq_mesh .^ 2          ...
% %                                           + x_radial_freq_mesh .^ 2          ...
% %                                           + z_radial_freq_mesh .^ 2 ) .^ 0.5 ;
% 
%         radial_freq_mesh_sphere_pulse = (   y_radial_freq_mesh_squared          ...
%                                           + x_radial_freq_mesh_squared          ...
%                                           + z_radial_freq_mesh_squared ) .^ 0.5 ;
%                                       
% %         % stably deconvolve PSF
% %         deconvolved_Gaussian_lengths = max( Gaussian_length ^ 2 - microns_per_sigma_PSF .^ 2, 0 ) .^ 0.5 ;
% %         
% %         y_radial_freq_mesh = y_micron_freq_mesh * deconvolved_Gaussian_lengths( 1 );
% %         x_radial_freq_mesh = x_micron_freq_mesh * deconvolved_Gaussian_lengths( 2 );
% %         z_radial_freq_mesh = z_micron_freq_mesh * deconvolved_Gaussian_lengths( 3 );
% 
% %         % add PSF
% %         y_radial_freq_mesh = y_micron_freq_mesh * ( Gaussian_length ^ 2 +     microns_per_sigma_PSF( 1 ) ^ 2 ) ^ 0.5 ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * ( Gaussian_length ^ 2 +     microns_per_sigma_PSF( 2 ) ^ 2 ) ^ 0.5 ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * ( Gaussian_length ^ 2 +     microns_per_sigma_PSF( 3 ) ^ 2 ) ^ 0.5 ;
%         
% %         % ignore PSF
% %         y_radial_freq_mesh = y_micron_freq_mesh * Gaussian_length ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * Gaussian_length ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * Gaussian_length ;
% 
%         y_radial_freq_mesh = y_micron_freq_mesh * Gaussian_lengths( 1 );
%         x_radial_freq_mesh = x_micron_freq_mesh * Gaussian_lengths( 2 );
%         z_radial_freq_mesh = z_micron_freq_mesh * Gaussian_lengths( 3 );
% 
% %         gaussian_PSF_kernel_dft = ones( size( y_radial_freq_mesh ));
%         
%         radial_freq_mesh_gaussian     = (   y_radial_freq_mesh .^ 2          ...
%                                           + x_radial_freq_mesh .^ 2          ...
%                                           + z_radial_freq_mesh .^ 2 ) .^ 0.5 ;
%                                                                             
% %         y_radial_freq_mesh = y_micron_freq_mesh * ( Gaussian_lengths( 1 ) ^ 2 + 2 * microns_per_sigma_PSF( 1 ) ^ 2 ) ^ 0.5 ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * ( Gaussian_lengths( 2 ) ^ 2 + 2 * microns_per_sigma_PSF( 2 ) ^ 2 ) ^ 0.5 ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * ( Gaussian_lengths( 3 ) ^ 2 + 2 * microns_per_sigma_PSF( 3 ) ^ 2 ) ^ 0.5 ;
% 
% %         y_radial_freq_mesh = y_micron_freq_mesh * ( sphere_pulse_lengths( 1 ) ^ 2 + Gaussian_lengths( 1 ) ^ 2 + 2 * microns_per_sigma_PSF( 1 ) ^ 2 ) ^ 0.5 ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * ( sphere_pulse_lengths( 2 ) ^ 2 + Gaussian_lengths( 2 ) ^ 2 + 2 * microns_per_sigma_PSF( 2 ) ^ 2 ) ^ 0.5 ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * ( sphere_pulse_lengths( 3 ) ^ 2 + Gaussian_lengths( 3 ) ^ 2 + 2 * microns_per_sigma_PSF( 3 ) ^ 2 ) ^ 0.5 ;
% 
% %         % add PSF
% %         y_radial_freq_mesh = y_micron_freq_mesh * ( sphere_pulse_length ^ 2 + Gaussian_length ^ 2 +     microns_per_sigma_PSF( 1 ) ^ 2 ) ^ 0.5 ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * ( sphere_pulse_length ^ 2 + Gaussian_length ^ 2 +     microns_per_sigma_PSF( 2 ) ^ 2 ) ^ 0.5 ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * ( sphere_pulse_length ^ 2 + Gaussian_length ^ 2 +     microns_per_sigma_PSF( 3 ) ^ 2 ) ^ 0.5 ;
% 
% %         % add PSF, derivative kernels at half sigma from origin
% %         y_radial_freq_mesh = y_micron_freq_mesh * ( sphere_pulse_length ^ 2 + Gaussian_length ^ 2 +     microns_per_sigma_PSF( 1 ) ^ 2 ) ^ 0.5 / 2 ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * ( sphere_pulse_length ^ 2 + Gaussian_length ^ 2 +     microns_per_sigma_PSF( 2 ) ^ 2 ) ^ 0.5 / 2 ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * ( sphere_pulse_length ^ 2 + Gaussian_length ^ 2 +     microns_per_sigma_PSF( 3 ) ^ 2 ) ^ 0.5 / 2 ;
% 
% %         % stably deconvolve PSF
% %         y_radial_freq_mesh = y_micron_freq_mesh * ( sphere_pulse_length ^ 2 + deconvolved_Gaussian_lengths ^ 2 ) ^ 0.5 ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * ( sphere_pulse_length ^ 2 + deconvolved_Gaussian_lengths ^ 2 ) ^ 0.5 ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * ( sphere_pulse_length ^ 2 + deconvolved_Gaussian_lengths ^ 2 ) ^ 0.5 ;
% 
% %         % ignore PSF
% %         y_radial_freq_mesh = y_micron_freq_mesh * ( sphere_pulse_length ^ 2 + Gaussian_length ^ 2 ) ^ 0.5 ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * ( sphere_pulse_length ^ 2 + Gaussian_length ^ 2 ) ^ 0.5 ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * ( sphere_pulse_length ^ 2 + Gaussian_length ^ 2 ) ^ 0.5 ;
% 
% %         % variance of a rectangular distribution = length^2/12 = radius^2/3 
% %         y_radial_freq_mesh = y_micron_freq_mesh * ( sphere_pulse_lengths( 1 ) ^ 2 / 3 + Gaussian_lengths( 1 ) ^ 2 + 2 * microns_per_sigma_PSF( 1 ) ^ 2 ) ^ 0.5 ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * ( sphere_pulse_lengths( 2 ) ^ 2 / 3 + Gaussian_lengths( 2 ) ^ 2 + 2 * microns_per_sigma_PSF( 2 ) ^ 2 ) ^ 0.5 ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * ( sphere_pulse_lengths( 3 ) ^ 2 / 3 + Gaussian_lengths( 3 ) ^ 2 + 2 * microns_per_sigma_PSF( 3 ) ^ 2 ) ^ 0.5 ;
% 
% %         y_radial_freq_mesh = y_micron_freq_mesh * ( sphere_pulse_length ^ 2 / 3 + Gaussian_length ^ 2 +     microns_per_sigma_PSF( 1 ) ^ 2 ) ^ 0.5 ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * ( sphere_pulse_length ^ 2 / 3 + Gaussian_length ^ 2 +     microns_per_sigma_PSF( 2 ) ^ 2 ) ^ 0.5 ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * ( sphere_pulse_length ^ 2 / 3 + Gaussian_length ^ 2 +     microns_per_sigma_PSF( 3 ) ^ 2 ) ^ 0.5 ;
% 
%         % only count blurring toward derivative weights
% %         y_radial_freq_mesh = y_micron_freq_mesh * Gaussian_length ;
% %         x_radial_freq_mesh = x_micron_freq_mesh * Gaussian_length ;
% %         z_radial_freq_mesh = z_micron_freq_mesh * Gaussian_length ;
%         
% %         radial_freq_mesh = radius_of_lumen_in_microns * micron_freq_mesh ;
% %         
% %         radial_freq_mesh_gaussian     =       A   * radial_freq_mesh ;
% %         
% %         radial_freq_mesh_sphere_pulse = ( 1 - A ) * radial_freq_mesh ;
% 
% %         radial_freq_mesh_gaussian     = Gaussian_length * micron_freq_mesh ;
% %         
% %         radial_freq_mesh_sphere_pulse = sphere_pulse_length * micron_freq_mesh ;
% 
%         sigma_freq_mesh          =          radial_freq_mesh_gaussian     ;
%         
%         radial_angular_freq_mesh = 2 * pi * radial_freq_mesh_sphere_pulse ;
%         
%         % these will have total sums of 1        
%         gaussian_kernel_dft = exp( - pi ^ 2 * 2 * sigma_freq_mesh .^ 2 );
% 
% %         spherical_pulse_kernel_dft                                                             ...
% %                     = 3 ./ ( radial_angular_freq_mesh .^ 2 )                                   ...
% %                         .* (   sin( radial_angular_freq_mesh ) ./ ( radial_angular_freq_mesh ) ...
% %                              - cos( radial_angular_freq_mesh )                                 );
% 
%         spherical_pulse_kernel_dft =  ( pi / 2 ./ radial_angular_freq_mesh ) .^ 0.5 ...
%                                    .* (   besselj( 2.5, radial_angular_freq_mesh )  ...
%                                         + besselj( 0.5, radial_angular_freq_mesh )) ;
%                                    
%         spherical_pulse_kernel_dft( radial_angular_freq_mesh == 0 ) = 1 ;
% 
% %         spherical_pulse_kernel_dft = ones( size( y_radial_freq_mesh ));
%         
%         matching_kernel_dft = gaussian_kernel_dft .* spherical_pulse_kernel_dft ;
%                 
%         % only count blurring toward derivative weights
% %         radius_of_lumen_in_voxels = ( Gaussian_length ^ 2 + microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ;
% %         radius_of_lumen_in_voxels = Gaussian_length ./ microns_per_pixel ;
% %         radius_of_lumen_in_voxels = deconvolved_Gaussian_lengths ./ microns_per_pixel ;
%         
% %         radius_of_lumen_in_voxels = max( radius_of_lumen_in_microns ^ 2, microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ;
% %         radius_of_lumen_in_voxels = Gaussian_lengths ./ microns_per_pixel ;
%         radius_of_lumen_in_voxels = max( Gaussian_lengths .^ 2 - microns_per_sigma_PSF .^ 2, microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ;
% 
% %% filters for endogenous signal images (vessel walls (endothelial cells) are lit instead of lumens)
% 
%     case '3D gaussian conv annular pulse'
        
%         radius_of_lumen_at_next_scale = radius_of_lumen_in_microns * 2 ^ ( 1 / scales_per_octave_radius );
        
%         Gaussian_lengths = max( microns_per_sigma_PSF, ( min( symmetry_ratio_factor ^ 2 * (   ( radius_of_lumen_at_next_scale + vessel_wall_thickness_in_microns / 2 ) ^ 2             ...
%                                                                                             - ( radius_of_lumen_in_microns    + vessel_wall_thickness_in_microns / 2 ) ^ 2 ),          ...
%                                                               ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ^ 2 ) - microns_per_sigma_PSF .^ 2        ) .^ 0.5 );

%         Gaussian_lengths = max( microns_per_sigma_PSF, ( min( symmetry_ratio_factor ^ 2 * ( radius_of_lumen_at_next_scale - radius_of_lumen_in_microns ) ^ 2 ,          ...
%                                                               ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ^ 2 ) - microns_per_sigma_PSF .^ 2        ) .^ 0.5 );

        

%         gaussian_to_ideal_ratio = max(          gaussian_to_ideal_ratio, ... % SAM 6/10/22
%                             microns_per_pixel > gaussian_to_ideal_ratio .* ( radius_of_lumen_in_microns ^ 2 + microns_per_sigma_PSF .^ 2 ) .^ 0.5 ); % at least one voxel of blurring in each dimension before any non-gaussian kernel is used

% %         Gaussian_lengths = max( microns_per_sigma_PSF, (( symmetry_ratio_factor * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 )) ^ 2 - microns_per_sigma_PSF .^ 2 ) .^ 0.5 );
% %         Gaussian_lengths = max( microns_per_sigma_PSF, (( gaussian_to_ideal_ratio * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 )) ^ 2 + microns_per_sigma_PSF .^ 2 ) .^ 0.5 );
% %         Gaussian_lengths = ( gaussian_to_ideal_ratio ^ 2 * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ^ 2 + microns_per_sigma_PSF .^ 2 ) .^ 0.5 ; % SAM 12/8/21, worse than 11/8/21 but unbiased radius estimation
%         Gaussian_lengths = gaussian_to_ideal_ratio * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) + [ 0, 0, 0 ]; % SAM 11/8/21 best to date, needs radius adjustmnet
%         Gaussian_lengths = gaussian_to_ideal_ratio * (( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ^ 2 + microns_per_sigma_PSF .^ 2 ) .^ 0.5 ; % SAM 12/22/21

%         Gaussian_lengths = gaussian_to_ideal_ratio .* ( radius_of_lumen_in_microns ^ 2 + microns_per_sigma_PSF .^ 2 ) .^ 0.5 ; % SAM 6/10/22
        Gaussian_lengths = gaussian_to_ideal_ratio * radius_of_lumen_in_microns + [ 0, 0, 0 ]; % SAM 7/11/22

%         target_microns_per_sigma_PSF_squared = mean( microns_per_sigma_PSF .^ 2 ); % correct half by (stable) deconvolution, blur the other half to achieve spherical blurring (and de-noising)
% %         target_microns_per_sigma_PSF_squared = max( microns_per_sigma_PSF .^ 2 );
%         
%         corrective_microns_per_sigma_PSF_squared = target_microns_per_sigma_PSF_squared ... 
%                                                  -                microns_per_sigma_PSF .^ 2 ;
        
%         Gaussian_lengths = max(( gaussian_to_ideal_ratio ^ 2 * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ^ 2 + corrective_microns_per_sigma_PSF_squared ), 0 ) .^ 0.5 ; % SAM 12/9/21

% %         Gaussian_lengths = max( [ 0, 0, 0 ], (( gaussian_to_ideal_ratio * (
% %         radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 )) ^ 2 - microns_per_sigma_PSF .^ 2 ) .^ 0.5 ); % miserable fail
% %         Gaussian_lengths = max([ 1, 1, 1 ], (( gaussian_to_ideal_ratio * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 )) ^ 2 - microns_per_sigma_PSF .^ 2 )) .^ 0.5 ; % 11/ 9 /21
% %         Gaussian_lengths = max( microns_per_sigma_PSF .^ 2, (( gaussian_to_ideal_ratio * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 )) ^ 2 - microns_per_sigma_PSF .^ 2 )) .^ 0.5 ;
        
        % assuming that the squared length of the combined kernel is the sum of the squared lengths
        % of the two kernels being convolved.
%         annular_pulse_lengths_squared = max((( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ^ 2 - Gaussian_lengths .^ 2 - microns_per_sigma_PSF .^ 2 ), 0 );        
%         annular_pulse_lengths_squared = max((( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ^ 2 - Gaussian_lengths .^ 2 + microns_per_sigma_PSF .^ 2 ), 0 );        
%         annular_pulse_lengths_squared = max(([ 1, 1, 1 ] * ( 1 - gaussian_to_ideal_ratio ^ 2 ) * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ^ 2 ), 0 );     % PREVIOUS VERSION for all of time   
%         annular_pulse_lengths_squared = max(([ 1, 1, 1 ] * ( 1 - gaussian_to_ideal_ratio ^ 2 ) * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ^ 2 - corrective_microns_per_sigma_PSF_squared ), 0 );        
%         Gaussian_lengths              = max(                     ( gaussian_to_ideal_ratio ^ 2   * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ^ 2 - microns_per_sigma_PSF .^ 2 ), 0 ) .^ 0.5 ; % SAM 12/9/21
%         annular_pulse_lengths_squared =       ([ 1, 1, 1 ] * ( 1 - gaussian_to_ideal_ratio ^ 2 ) * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ^ 2 + microns_per_sigma_PSF .^ 2 )             ;        
%         annular_pulse_lengths_squared =                      ( 1 - gaussian_to_ideal_ratio ^ 2 ) * (( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ^ 2 - microns_per_sigma_PSF .^ 2 )  ; % SAM 12/22/21
%         annular_pulse_lengths_squared =                      ( 1 - gaussian_to_ideal_ratio .^ 2 ) .* ( radius_of_lumen_in_microns ^ 2 + microns_per_sigma_PSF .^ 2 ); % SAM 6/10/22
        annular_pulse_lengths_squared =                      ( 1 - gaussian_to_ideal_ratio ^ 2 ) * radius_of_lumen_in_microns ^ 2 + microns_per_sigma_PSF .^ 2 ; % SAM 7/11/22

        % spherical pulse addition section:
        
%         % spherical component weight ( 1 is same as annular component, 2 is twice as, ... )
% %         A = 6 ; % 7/30/19 empirical datasets with spherical signal (sometimes spherical+ annular)

%         A = 0 ; % noise study for annular signal
        
% %         A = 1 ; % 7/31/19 empirical dataset with annular signal (sometimes annular + spherical) !!!! didn't work that well, perhaps too much annulus, perhaps more smoothing should have been done on input (and then added to the "PSF")
% %         A = 2 ;
% %         A = 0.5 ; % !! best so far
% %         A = 0.25 ;
% %         A = 0.75 ;

%         % fraction of spherical component of total ( 0.5 is same as annular component, ... )
%         A = 6 / 7 ; % 7/30/19 empirical datasets with spherical signal (sometimes spherical+ annular)
        
        % michaels empirical dataset with annular signal
%         A = 1 ;
%         A = 0 ;
%         A = 1/2 ;
%         A = 1/2 ; % best so far, to be used with 50% gauss
%         A = 2/3 ;
%         A = 1/3 ;
%         A = 0.1 ;

%         A = 1 ; % noise study for spherical signal

%         sphere_pulse_lengths_squared = max((  radius_of_lumen_in_microns                                          ^ 2 - Gaussian_lengths .^ 2 - microns_per_sigma_PSF .^ 2 ), 0 );        
%         sphere_pulse_lengths_squared = max((  radius_of_lumen_in_microns                                          ^ 2 - Gaussian_lengths .^ 2 + microns_per_sigma_PSF .^ 2 ), 0 );       
%         sphere_pulse_lengths_squared = ( 1 - gaussian_to_ideal_ratio ^ 2 ) * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ^ 2 ...
%                                      * [ 1, 1, 1 ];
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
        
        
        radial_angular_freq_mesh = 2 * pi * radial_freq_mesh_sphere_pulse ;
        
%         y_radial_freq_mesh_squared_B = y_micron_freq_mesh .^ 2 * sphere_pulse_lengths_squared_B( 1 );
%         x_radial_freq_mesh_squared_B = x_micron_freq_mesh .^ 2 * sphere_pulse_lengths_squared_B( 2 );
%         z_radial_freq_mesh_squared_B = z_micron_freq_mesh .^ 2 * sphere_pulse_lengths_squared_B( 3 );
% 
%         radial_freq_mesh_sphere_pulse_B = (   y_radial_freq_mesh_squared_B          ...
%                                             + x_radial_freq_mesh_squared_B          ...
%                                             + z_radial_freq_mesh_squared_B ) .^ 0.5 ;        
%         
%         
%         radial_angular_freq_mesh_B = 2 * pi * radial_freq_mesh_sphere_pulse_B ;        
%         
        % these will have total sums of 1        
        gaussian_kernel_dft = exp( - pi ^ 2 * 2 * sigma_freq_mesh .^ 2 );
% 
%         spherical_pulse_kernel_dft =  (   pi  /  2  ./  radial_angular_freq_mesh ) .^ 0.5 ...
%                                    .* (   besselj( 2.5, radial_angular_freq_mesh )        ...
%                                         + besselj( 0.5, radial_angular_freq_mesh ))       ;
%                                    
%         spherical_pulse_kernel_dft( radial_angular_freq_mesh == 0 ) = 1 ;
%         
%         spherical_pulse_kernel_dft_B =  (   pi  /  2  ./  radial_angular_freq_mesh_B ) .^ 0.5 ...
%                                      .* (   besselj( 2.5, radial_angular_freq_mesh_B )  ...
%                                           + besselj( 0.5, radial_angular_freq_mesh_B )) ;
%                                    
%         spherical_pulse_kernel_dft_B( radial_angular_freq_mesh == 0 ) = 1 ;
%         
%         spherical_volume_small_in_cubic_pxls = 4 / 3 * pi * sphere_pulse_lengths_squared_B ^ ( 3 / 2 ) / prod( microns_per_pixel );                              
%         spherical_volume_large_in_cubic_pxls = 4 / 3 * pi * sphere_pulse_lengths_squared   ^ ( 3 / 2 ) / prod( microns_per_pixel );
%         
%         annular_volume_in_cubic_pxls = 4 / 3 * pi                           ...
%                                         * (   radius_large_in_microns ^ 3   ...
%                                             - radius_small_in_microns ^ 3 ) ...
%                                         / prod( microns_per_pixel );
%         
%         annular_pulse_kernel_dft = spherical_pulse_kernel_dft - spherical_pulse_kernel_dft_B ;
                              
        % consider integrating (and averaging) cos() over different values for the radius to span
        % the vessel wall (instead of idealizing it as an infinitely thin source at a single radius )
        annular_pulse_kernel_dft = cos( radial_angular_freq_mesh );
        
%         matching_kernel_dft = gaussian_kernel_dft .* annular_pulse_kernel_dft ;
%         matching_kernel_dft = gaussian_kernel_dft .* ( annular_pulse_kernel_dft + A * spherical_pulse_kernel_dft ) / ( 1 + A );
        matching_kernel_dft = gaussian_kernel_dft .* (( 1 - spherical_to_annular_ratio ) * annular_pulse_kernel_dft + spherical_to_annular_ratio * spherical_pulse_kernel_dft );
        
%         % BEGIN PSF correction section
%         % correct for PSF without blurring by adding a spherical pulse with size equal to one sigma
%         % of PSF
% %         sphere_pulse_lengths_squared = microns_per_sigma_PSF .^ 2 ;
%         sphere_pulse_lengths_squared = [ 1, 1, 1 ] * max( microns_per_sigma_PSF ) ^ 2 ;
%                                
%         y_radial_freq_mesh_squared = y_micron_freq_mesh .^ 2 * sphere_pulse_lengths_squared( 1 );
%         x_radial_freq_mesh_squared = x_micron_freq_mesh .^ 2 * sphere_pulse_lengths_squared( 2 );
%         z_radial_freq_mesh_squared = z_micron_freq_mesh .^ 2 * sphere_pulse_lengths_squared( 3 );
% 
%         radial_freq_mesh_sphere_pulse = (   y_radial_freq_mesh_squared          ...
%                                           + x_radial_freq_mesh_squared          ...
%                                           + z_radial_freq_mesh_squared ) .^ 0.5 ;
%                      
%         radial_angular_freq_mesh = 2 * pi * radial_freq_mesh_sphere_pulse ;
%         
%         PSF_correction_dft = ( pi / 2 ./ radial_angular_freq_mesh ) .^ 0.5 ...
%                           .* (   besselj( 2.5, radial_angular_freq_mesh )  ...
%                                + besselj( 0.5, radial_angular_freq_mesh )) ;
%                                    
%         PSF_correction_dft( radial_angular_freq_mesh == 0 ) = 1 ;
%        % END PSF correctino section
        
%         matching_kernel_dft = matching_kernel_dft .* PSF_correction_dft ;
        
%         % only count blurring toward derivative weights
%         radius_of_lumen_in_voxels = Gaussian_lengths ./ microns_per_pixel ; % SAM 11/8/21 best to date, needs radius adjustmnet
        derivative_weights_from_blurring = Gaussian_lengths ./ microns_per_pixel ; % SAM 12/14/21

%             Gaussian_lengths_in_pixels =             Gaussian_lengths        ./ microns_per_pixel ;
%         
%         sphere_pulse_lengths_in_pixels = sphere_pulse_lengths_squared .^ 0.5 ./ microns_per_pixel ;
% 
%         derivative_weights_from_blurring = (                Gaussian_lengths_in_pixels                                              ...
%                                              .* (    [     Gaussian_lengths_in_pixels( 2 ) *     Gaussian_lengths_in_pixels( 3 ), ...
%                                                            Gaussian_lengths_in_pixels( 3 ) *     Gaussian_lengths_in_pixels( 1 ), ...
%                                                            Gaussian_lengths_in_pixels( 1 ) *     Gaussian_lengths_in_pixels( 2 )  ]           ... 
%                                                   +  [ sphere_pulse_lengths_in_pixels( 2 ) * sphere_pulse_lengths_in_pixels( 3 ), ...
%                                                        sphere_pulse_lengths_in_pixels( 3 ) * sphere_pulse_lengths_in_pixels( 1 ), ...
%                                                        sphere_pulse_lengths_in_pixels( 1 ) * sphere_pulse_lengths_in_pixels( 2 )  ]   ) .^ 0.5 ) .^ 0.5 ; % SAM 12/20/21

%         %  only count non-corrective blurring toward the derivative weighting (make it think the
%         %  vessel is bigger than it is.) will correct vessel sizes after edge extraction. 12/9/21
%         derivative_weights_from_blurring = gaussian_to_ideal_ratio * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 );
%         derivative_weights_from_blurring = ...
%         derivative_weights_from_blurring ./ microns_per_pixel ;

%         % only count (non-PSF matching) blurring toward derivative weights
%         radius_of_lumen_in_voxels = (( Gaussian_lengths .^ 2 - microns_per_sigma_PSF .^ 2 ) .^ 0.5 ) ./ microns_per_pixel ; % SAM 12/8/21
%         radius_of_lumen_in_voxels = max( Gaussian_lengths .^ 2 - microns_per_sigma_PSF .^ 2, microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ;
%         radius_of_lumen_in_voxels = radius_of_lumen_in_microns ./ microns_per_pixel ; % never attempted SAM 11/6/21
%         radius_of_lumen_in_voxels = (( gaussian_to_ideal_ratio * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 )) ^ 2 + 2 * microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ; % SAM 11/6/21
%         radius_of_lumen_in_voxels = (( gaussian_to_ideal_ratio * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 )) ^ 2 + microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ; % SAM 11/6/21
%         radius_of_lumen_in_voxels = gaussian_to_ideal_ratio * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ./ microns_per_pixel ; % SAM 11/6/21, same as: % SAM 11/8/21 best to date, needs radius adjustmnet
        
%     case 'annular pulse'
% 
%         
%         % do a difference of spherical pulses each with constant value approximately one inside
%         % their spheres. Normalize the difference to have this kernel sum to 1.
%         radius_small_in_microns = radius_of_lumen_in_microns ;
%         radius_large_in_microns = radius_of_lumen_in_microns + vessel_wall_thickness_in_microns ;
% 
%         radial_small_angular_freq_mesh = 2 * pi * micron_freq_mesh * radius_small_in_microns ;
%         radial_large_angular_freq_mesh = 2 * pi * micron_freq_mesh * radius_large_in_microns ;
% 
% 
%         spherical_volume_small_in_cubic_pxls = 4 / 3 * pi * radius_small_in_microns ^ 3 / prod( microns_per_pixel );                              
%         spherical_volume_large_in_cubic_pxls = 4 / 3 * pi * radius_large_in_microns ^ 3 / prod( microns_per_pixel );
% 
%         sphere_pulse_small_dft = 3 ./ ( radial_small_angular_freq_mesh .^ 2 )                                     ...
%                                .* (   sin( radial_small_angular_freq_mesh ) ./ ( radial_small_angular_freq_mesh ) ...
%                                - cos( radial_small_angular_freq_mesh )                                            );
% 
%         sphere_pulse_small_dft( 1 ) = 1 ;    
% 
%         sphere_pulse_small_dft = sphere_pulse_small_dft * spherical_volume_small_in_cubic_pxls ;
% 
%         sphere_pulse_large_dft = 3 ./ ( radial_large_angular_freq_mesh .^ 2 )                                     ...
%                                .* (   sin( radial_large_angular_freq_mesh ) ./ ( radial_large_angular_freq_mesh ) ...
%                                - cos( radial_large_angular_freq_mesh )                                            );
% 
%         sphere_pulse_large_dft( 1 ) = 1 ;    
% 
%         sphere_pulse_large_dft = sphere_pulse_large_dft * spherical_volume_large_in_cubic_pxls ;
% 
%         % normalize to have sum of 1        
%         matching_kernel_dft = sphere_pulse_large_dft - sphere_pulse_small_dft ;
% 
%         % comment out these two line to have instead a constant value of 1 (not a total sum of 1).        
%         annular_volume_in_cubic_pxls = 4 / 3 * pi                           ...
%                                         * (   radius_large_in_microns ^ 3   ...
%                                             - radius_small_in_microns ^ 3 ) ...
%                                         / prod( microns_per_pixel );        
%         
%         matching_kernel_dft = matching_kernel_dft / annular_volume_in_cubic_pxls ;
%         
%     case 'annular pulse V2'
%         
%         % do a difference of spherical pulses, one with a constant value approximately one inside,
%         % the other with a constant value A inside (where 0 < A < 1 ). Normalize their difference to
%         % have the resulting kernel sum to 1.
%         
% %         A = 0.5 ; % 8/16/18
% %         A = 0.25 ; % 8/16/18
% %         A = 0.1 ; % 8/16/18
%         spherical_to_annular_ratio = 0.25 ; % 8/16/18
%         
%         radius_small_in_microns = radius_of_lumen_in_microns ;
%         radius_large_in_microns = radius_of_lumen_in_microns + vessel_wall_thickness_in_microns ;
% 
%         radial_small_angular_freq_mesh = 2 * pi * micron_freq_mesh * radius_small_in_microns ;
%         radial_large_angular_freq_mesh = 2 * pi * micron_freq_mesh * radius_large_in_microns ;
% 
% 
%         spherical_volume_small_in_cubic_pxls = spherical_to_annular_ratio * 4 / 3 * pi * radius_small_in_microns ^ 3 / prod( microns_per_pixel );                              
%         spherical_volume_large_in_cubic_pxls =     4 / 3 * pi * radius_large_in_microns ^ 3 / prod( microns_per_pixel );
% 
%         % this kernel has average value of 1
%         sphere_pulse_small_dft = 3 ./ ( radial_small_angular_freq_mesh .^ 2 )                                     ...
%                                .* (   sin( radial_small_angular_freq_mesh ) ./ ( radial_small_angular_freq_mesh ) ...
%                                - cos( radial_small_angular_freq_mesh )                                            );
% 
%         sphere_pulse_small_dft( 1 ) = 1 ;    
% 
%         % this kernel has approxiamtely constant value of A inside the sphere
%         sphere_pulse_small_dft = sphere_pulse_small_dft * spherical_volume_small_in_cubic_pxls ;
% 
%         % this kernel has average value of 1
%         sphere_pulse_large_dft = 3 ./ ( radial_large_angular_freq_mesh .^ 2 )                                     ...
%                                .* (   sin( radial_large_angular_freq_mesh ) ./ ( radial_large_angular_freq_mesh ) ...
%                                - cos( radial_large_angular_freq_mesh )                                            );
% 
%         sphere_pulse_large_dft( 1 ) = 1 ;    
% 
%         % this kernel has approxiamtely constant value of 1 inside the sphere        
%         sphere_pulse_large_dft = sphere_pulse_large_dft * spherical_volume_large_in_cubic_pxls ;
% 
%         % this kernel has max value of 1
%         matching_kernel_dft = sphere_pulse_large_dft - sphere_pulse_small_dft ;
% 
%         % normalizing to have sum of 1  
%         %
%         % comment out these two line to have instead a max value of 1 (not a total sum of 1).        
%         annular_volume_in_cubic_pxls =    4 / 3 * pi                            ...
%                                         * (       radius_large_in_microns ^ 3   ...
%                                             - spherical_to_annular_ratio * radius_small_in_microns ^ 3 ) ...
%                                         /   prod( microns_per_pixel );        
%         
%         matching_kernel_dft = matching_kernel_dft / annular_volume_in_cubic_pxls ;
%         
%     case '3D annular gaussian'
%         
% %         A = 0.25 ; % 8/16/18        
% %         A = 1 ; % 8/17/18        
% %         A = 0.75 ; % 8/17/18
% %         A = 0.1 ; % 8/17/18 112200
% %         A = 0.5 ; % 8/17/18 113500       
%         
% %         radius_in_microns = 5 ; % approximate max radius of a capillary 8/17/18 
% %         radius_in_microns = 10 ; % 8/17/18 
%         radius_in_microns = 7.5 ; % 8/17/18 135000
% %         radius_in_microns = 5 ;
% 
%         % assuming that the loss of illumination at the on axis points in the vessel is an
%         % exponentially decreasing function of the cross-sectional area of the vessel
%         spherical_to_annular_ratio = 1 - exp( - ( radius_of_lumen_in_microns / radius_in_microns ) ^ 2 );
% 
% %         % assuming that the loss of illumination at the on axis points in the vessel is an
% %         % exponentially decreasing function of the radius of the vessel
% %         A = 1 - exp( - radius_of_lumen_in_microns / radius_in_microns ); % SAM 8/17/18 151200
% 
%         
%         sigma_small_freq_mesh = radial_freq_mesh ;
%         
%         sigma_large_freq_mesh = micron_freq_mesh                                                  ...
%                               * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns ) ;
%         
% %       % these kernels have average values of 1 
%         gaussian_kernel_small_dft = exp( - pi ^ 2 * 2 * sigma_small_freq_mesh .^ 2 );
%         gaussian_kernel_large_dft = exp( - pi ^ 2 * 2 * sigma_large_freq_mesh .^ 2 );
%         
%         max_of_small_gaussian_kernel = ( 2 * pi )  .^ - 0.5 /     radius_of_lumen_in_microns ;
%         max_of_large_gaussian_kernel = ( 2 * pi )  .^ - 0.5 / (   radius_of_lumen_in_microns       ...
%                                                                 + vessel_wall_thickness_in_microns );
%                                                     
%         % these kernels have max values of ( A and 1 ) * B, for some constant B that depends on the
%         % pixel spacing
%         gaussian_kernel_small_dft = spherical_to_annular_ratio * gaussian_kernel_small_dft / max_of_small_gaussian_kernel ;
%         gaussian_kernel_large_dft = 	gaussian_kernel_large_dft / max_of_large_gaussian_kernel ;
%         
%         % this kernel has average value of 
%         %
%         % ( 1 / max_of_large_gaussian_kernel - A / max_of_small_gaussian_kernel )
%         matching_kernel_dft = gaussian_kernel_large_dft - gaussian_kernel_small_dft ;
%         
%         average_value_of_differenced_gaussians = ( 1 / max_of_large_gaussian_kernel - spherical_to_annular_ratio / max_of_small_gaussian_kernel );
%         
%         % this kernel has average value of 1
%         matching_kernel_dft = matching_kernel_dft / average_value_of_differenced_gaussians ;
%         
%     case '3D annular gaussian V2'
%         
%         % 0 < A < 1
% %         A = 0.1 ; % 8/24/18
% %         A = 0.2 ; % 8/25/18
% %         A = 0.5 ; % 8/25/18
%         spherical_to_annular_ratio = 0.3 ; % 8/25/18
% %         A = 0.15 ; % 8/25/18
% 
% 
% %         radius_in_microns = 7.5 ; % 8/24/18
% % 
% %         % assuming that the portion of average intensity on the wall that contributes to the average
% %         % intensity in the lumen is an exponentially decreasing function of the radius of the vessel
% %         A = 1 - exp( - ( radius_of_lumen_in_microns / radius_in_microns ));
% 
% 
%         sigma_large_freq_mesh = micron_freq_mesh                                                  ...
%                               * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns ) ;
%         
%         sigma_small_freq_mesh = radial_freq_mesh ;                          
%                           
%         % these kernels have average values of 1 and A
%         gaussian_kernel_large_dft = exp( - pi ^ 2 * 2 * sigma_large_freq_mesh .^ 2 );
%         gaussian_kernel_small_dft = exp( - pi ^ 2 * 2 * sigma_small_freq_mesh .^ 2 );        
%                                                                     
%         matching_kernel_dft = ( gaussian_kernel_large_dft - spherical_to_annular_ratio * gaussian_kernel_small_dft )         ...
%                                                                                           / ( 1 - spherical_to_annular_ratio );
%                 
%                 
%     case 'radial gaussian'
%         
%         spherical_to_annular_ratio = 2 ;
%         
%         
%         % vessel wall is approximated as a radial delta function convolved with a radial Gaussian
%         sigma_vessel_wall_freq_mesh = micron_freq_mesh * vessel_wall_thickness_in_microns / 2 ;
%                 
%         radial_lumen_angular_freq_mesh = 2 * pi * micron_freq_mesh * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ;
%                 
%         % these kernels have average values of 1 
%         gaussian_kernel_vessel_wall_dft = exp( - pi ^ 2 * 2 * sigma_vessel_wall_freq_mesh .^ 2 );
%         
%         radial_delta_fxn_lumen_dft      =  sin( radial_lumen_angular_freq_mesh ) ...
%                                         ./      radial_lumen_angular_freq_mesh ;
%                                     
%         spherical_pulse_dft =    3 ./ ( radial_lumen_angular_freq_mesh .^ 2 )                           ...
%                             .* (   sin( radial_lumen_angular_freq_mesh ) ./ ( radial_lumen_angular_freq_mesh ) ...
%                             -      cos( radial_lumen_angular_freq_mesh )                                 );            
%                                     
%                spherical_pulse_dft( 1 ) = 1 ;
%         radial_delta_fxn_lumen_dft( 1 ) = 1 ;
%         
%         % their spatial convolution is a radial gaussian function
%         matching_kernel_dft =   gaussian_kernel_vessel_wall_dft ...
%                             .* ( radial_delta_fxn_lumen_dft + spherical_to_annular_ratio * spherical_pulse_dft ) ...
%                             / ( 1 + spherical_to_annular_ratio );    
%                                                         
% end % matching kernel selection SWITCH

%         % uncomment to inspect the kernel in the spatial domain
% 
%         matching_kernel_image = ifftn( matching_kernel_dft, 'symmetric' );
% 
%         matching_kernel_image( 1 : 5, 1 : 5, 1 : 5 )
% 
%         sum( matching_kernel_image( : ))
                                                                
% matching_kernel_dft = matching_kernel_dft .* gaussian_PSF_kernel_dft ;

% sum( matching_kernel_dft( : ) .^ 2 ) ^ 0.5

%         % uncomment to inspect the blurred kernel in the spatial domain
% 
%         blurred_matching_kernel_image = ifftn( matching_kernel_dft, 'symmetric' );
% 
%         blurred_matching_kernel_image( 1 : 5, 1 : 5, 1 : 5 )
%         
%         sum( blurred_matching_kernel_image( : ))
              
%% energy calcs

if ~ is_only_computing_Laplacian

    % derivative filter constructions:

    gradient_kernels_dft   = zeros([ 3, size_of_chunk_dft ]);

    if is_doing_eigenvalue_decomp
        
        curvatures_kernels_dft = zeros([ 6, size_of_chunk_dft ]);        

    else
        
        curvatures_kernels_dft = zeros([ 3, size_of_chunk_dft ]);     
        
    end

    % % derivatives with respect to radius are taken by multiplication by 2 pi i rho where rho is the
    % % radius frequency
    %                         
    %   gradient_kernels_dft( 1, :, :, : ) =  2 * pi * 1i * pixels_per_radius( 1 ) * y_freq_mesh ;    
    %   gradient_kernels_dft( 2, :, :, : ) =  2 * pi * 1i * pixels_per_radius( 2 ) * x_freq_mesh ;    
    %   gradient_kernels_dft( 3, :, :, : ) =  2 * pi * 1i * pixels_per_radius( 3 ) * z_freq_mesh ;    
    % 
    % curvatures_kernels_dft( 1, :, :, : ) = - 4 * pi ^ 2 *       pixels_per_radius_squared( 1 ) * y_freq_mesh .^ 2           ;
    % curvatures_kernels_dft( 2, :, :, : ) = - 4 * pi ^ 2 *       pixels_per_radius_squared( 2 ) * x_freq_mesh .^ 2           ;
    % curvatures_kernels_dft( 3, :, :, : ) = - 4 * pi ^ 2 *       pixels_per_radius_squared( 3 ) * z_freq_mesh .^ 2           ;
    % curvatures_kernels_dft( 4, :, :, : ) = - 4 * pi ^ 2 * prod( pixels_per_radius([ 1, 2 ])  ) * y_freq_mesh .* x_freq_mesh ;
    % curvatures_kernels_dft( 5, :, :, : ) = - 4 * pi ^ 2 * prod( pixels_per_radius([ 2, 3 ])  ) * x_freq_mesh .* z_freq_mesh ;
    % curvatures_kernels_dft( 6, :, :, : ) = - 4 * pi ^ 2 * prod( pixels_per_radius([ 3, 1 ])  ) * z_freq_mesh .* y_freq_mesh ;

    % derivatives with respect to microns in a certain dimension are taken by multiplication by 2 pi i
    % rho where rho is the micron frequency in a that dimension.  This is different than the previous
    % aproach by giving more weight to the smaller sizes (no length of radius weighting factor given to
    % larger sizes) SAM circa 8/10/18

    % % derivatives with respect to radius (of the gaussian kernel).  SAM 8/13/18

    % curvatures_kernels_dft( 1, :, :, : ) = - 4 * pi ^ 2 * y_radial_freq_mesh .^ 2 ;
    % curvatures_kernels_dft( 2, :, :, : ) = - 4 * pi ^ 2 * x_radial_freq_mesh .^ 2 ;
    % curvatures_kernels_dft( 3, :, :, : ) = - 4 * pi ^ 2 * z_radial_freq_mesh .^ 2 ;
    % curvatures_kernels_dft( 4, :, :, : ) = - 4 * pi ^ 2 * y_radial_freq_mesh .* x_radial_freq_mesh ;
    % curvatures_kernels_dft( 5, :, :, : ) = - 4 * pi ^ 2 * x_radial_freq_mesh .* z_radial_freq_mesh ;
    % curvatures_kernels_dft( 6, :, :, : ) = - 4 * pi ^ 2 * z_radial_freq_mesh .* y_radial_freq_mesh ;

    % curvatures_kernels_dft( 1, :, :, : ) = 2 * ( cos( 2 * pi *      y_radial_freq_mesh                                ) - 1 );
    % curvatures_kernels_dft( 2, :, :, : ) = 2 * ( cos( 2 * pi *      x_radial_freq_mesh                                ) - 1 );
    % curvatures_kernels_dft( 3, :, :, : ) = 2 * ( cos( 2 * pi *      z_radial_freq_mesh                                ) - 1 );
    % curvatures_kernels_dft( 4, :, :, : ) = 2 * ( cos( 2 * pi * abs( y_radial_freq_mesh .* x_radial_freq_mesh ) .^ 0.5 ) - 1 ) .* sign( y_radial_freq_mesh .* x_radial_freq_mesh );
    % curvatures_kernels_dft( 5, :, :, : ) = 2 * ( cos( 2 * pi * abs( x_radial_freq_mesh .* z_radial_freq_mesh ) .^ 0.5 ) - 1 ) .* sign( x_radial_freq_mesh .* z_radial_freq_mesh );
    % curvatures_kernels_dft( 6, :, :, : ) = 2 * ( cos( 2 * pi * abs( z_radial_freq_mesh .* y_radial_freq_mesh ) .^ 0.5 ) - 1 ) .* sign( z_radial_freq_mesh .* y_radial_freq_mesh );

    % % neighboring voxel derivatives:
    curvatures_kernels_dft( 1, :, :, : ) = prod( derivative_weights_from_blurring([ 1, 1 ])) * ( cos( 2 * pi *      y_pixel_freq_mesh                               ) - 1 );
    curvatures_kernels_dft( 2, :, :, : ) = prod( derivative_weights_from_blurring([ 2, 2 ])) * ( cos( 2 * pi *      x_pixel_freq_mesh                               ) - 1 );
    curvatures_kernels_dft( 3, :, :, : ) = prod( derivative_weights_from_blurring([ 3, 3 ])) * ( cos( 2 * pi *      z_pixel_freq_mesh                               ) - 1 );
    
    if is_doing_eigenvalue_decomp
    
    % !!!!! are these mixed derivatives off by a factor 4 ??????? SAM 6/3/22 !!!!!!!!!!!!!!!!!!!
    curvatures_kernels_dft( 4, :, :, : ) = prod( derivative_weights_from_blurring([ 1, 2 ])) * ( cos( 2 * pi * abs( y_pixel_freq_mesh .* x_pixel_freq_mesh ) .^ 0.5 ) - 1 ) .* sign( y_pixel_freq_mesh .* x_pixel_freq_mesh ) / 4 ;
    curvatures_kernels_dft( 5, :, :, : ) = prod( derivative_weights_from_blurring([ 2, 3 ])) * ( cos( 2 * pi * abs( x_pixel_freq_mesh .* z_pixel_freq_mesh ) .^ 0.5 ) - 1 ) .* sign( x_pixel_freq_mesh .* z_pixel_freq_mesh ) / 4 ;
    curvatures_kernels_dft( 6, :, :, : ) = prod( derivative_weights_from_blurring([ 3, 1 ])) * ( cos( 2 * pi * abs( z_pixel_freq_mesh .* y_pixel_freq_mesh ) .^ 0.5 ) - 1 ) .* sign( z_pixel_freq_mesh .* y_pixel_freq_mesh ) / 4 ;

    end
    
    % curvatures_kernels_dft( 1, :, :, : ) = - y_radial_freq_mesh .^ 2 ;
    % curvatures_kernels_dft( 2, :, :, : ) = - x_radial_freq_mesh .^ 2 ;
    % curvatures_kernels_dft( 3, :, :, : ) = - z_radial_freq_mesh .^ 2 ;
    % curvatures_kernels_dft( 4, :, :, : ) = - y_radial_freq_mesh .* x_radial_freq_mesh ;
    % curvatures_kernels_dft( 5, :, :, : ) = - x_radial_freq_mesh .* z_radial_freq_mesh ;
    % curvatures_kernels_dft( 6, :, :, : ) = - z_radial_freq_mesh .* y_radial_freq_mesh ;

end

if is_only_computing_Laplacian

    % laplacian kernel:
    Laplacian_kernel_dft = 2 * derivative_weights_from_blurring( 1 ) ^ 2 * ( cos( 2 * pi * y_pixel_freq_mesh ) - 1 ) ...
                         + 2 * derivative_weights_from_blurring( 2 ) ^ 2 * ( cos( 2 * pi * x_pixel_freq_mesh ) - 1 ) ...
                         + 2 * derivative_weights_from_blurring( 3 ) ^ 2 * ( cos( 2 * pi * z_pixel_freq_mesh ) - 1 ) ;
                     
    % r_pixel_freq_mesh = (   y_pixel_freq_mesh .^ 2          ...
    %                       + x_pixel_freq_mesh .^ 2          ...
    %                       + z_pixel_freq_mesh .^ 2 ) .^ 0.5 ;
    % 
    % radius_of_lumen_in_voxels_squared = (   y_pixel_freq_mesh .^ 2 * radius_of_lumen_in_voxels( 1 ) ^ 2   ...
    %                                       + x_pixel_freq_mesh .^ 2 * radius_of_lumen_in_voxels( 2 ) ^ 2   ...
    %                                       + z_pixel_freq_mesh .^ 2 * radius_of_lumen_in_voxels( 3 ) ^ 2 ) ...
    %                                  ./ r_pixel_freq_mesh                                                 ;
    %                                   
    % radius_of_lumen_in_voxels_squared( 1 ) = 0 ;
    %                               
    % % Laplacian_kernel_dft = 2 * ( cos( 3 ^ 0.5 * 2 * pi * r_pixel_freq_mesh ) - 1 ) .* radius_of_lumen_in_voxels_squared ;
    % % Laplacian_kernel_dft = 2 * ( cos( 2 * pi * r_pixel_freq_mesh ) - 1 ) .* radius_of_lumen_in_voxels_squared ;
    % 
    % Laplacian_kernel_dft = 2 * radius_of_lumen_in_voxels_squared .* ( cos( 3 ^ 0.5 * 2 * pi * r_pixel_freq_mesh ) - 1 ) ...
    %                      + 2 * radius_of_lumen_in_voxels( 1 ) ^ 2 * ( cos(           2 * pi * y_pixel_freq_mesh ) - 1 ) ...
    %                      + 2 * radius_of_lumen_in_voxels( 2 ) ^ 2 * ( cos(           2 * pi * x_pixel_freq_mesh ) - 1 ) ...
    %                      + 2 * radius_of_lumen_in_voxels( 3 ) ^ 2 * ( cos(           2 * pi * z_pixel_freq_mesh ) - 1 ) ;

end

if ~ is_only_computing_Laplacian

    % gradient_kernels_dft( 1, :, :, : ) =  2 * pi * 1i * y_radial_freq_mesh ;
    % gradient_kernels_dft( 2, :, :, : ) =  2 * pi * 1i * x_radial_freq_mesh ;
    % gradient_kernels_dft( 3, :, :, : ) =  2 * pi * 1i * z_radial_freq_mesh ;

    % gradient_kernels_dft( 1, :, :, : ) = - 8 * pi ^ 3 * 1i * y_radial_freq_mesh .^ 3 ;
    % gradient_kernels_dft( 2, :, :, : ) = - 8 * pi ^ 3 * 1i * x_radial_freq_mesh .^ 3 ;
    % gradient_kernels_dft( 3, :, :, : ) = - 8 * pi ^ 3 * 1i * z_radial_freq_mesh .^ 3 ;

    %   gradient_kernels_dft( 1, :, :, : ) =  2 * pi * 1i * y_radial_freq_mesh - 8 * pi ^ 3 * 1i * y_radial_freq_mesh .^ 3 ;
    %   gradient_kernels_dft( 2, :, :, : ) =  2 * pi * 1i * x_radial_freq_mesh - 8 * pi ^ 3 * 1i * x_radial_freq_mesh .^ 3 ;
    %   gradient_kernels_dft( 3, :, :, : ) =  2 * pi * 1i * z_radial_freq_mesh - 8 * pi ^ 3 * 1i * z_radial_freq_mesh .^ 3 ;

    % gradient_kernels_dft( 1, :, :, : ) =  1i * ( y_radial_freq_mesh - y_radial_freq_mesh .^ 3 );
    % gradient_kernels_dft( 2, :, :, : ) =  1i * ( x_radial_freq_mesh - x_radial_freq_mesh .^ 3 );
    % gradient_kernels_dft( 3, :, :, : ) =  1i * ( z_radial_freq_mesh - z_radial_freq_mesh .^ 3 );

    % gradient_kernels_dft( 1, :, :, : ) =  2 * pi * 1i * y_radial_freq_mesh - 8 * pi ^ 3 * 1i * y_radial_freq_mesh .^ 3 / 6 ;
    % gradient_kernels_dft( 2, :, :, : ) =  2 * pi * 1i * x_radial_freq_mesh - 8 * pi ^ 3 * 1i * x_radial_freq_mesh .^ 3 / 6 ;
    % gradient_kernels_dft( 3, :, :, : ) =  2 * pi * 1i * z_radial_freq_mesh - 8 * pi ^ 3 * 1i * z_radial_freq_mesh .^ 3 / 6 ;

    % gradient_kernels_dft( 1, :, :, : ) =  2 * pi * 1i * y_radial_freq_mesh + 8 * pi ^ 3 * 1i * y_radial_freq_mesh .^ 3 / 6 ;
    % gradient_kernels_dft( 2, :, :, : ) =  2 * pi * 1i * x_radial_freq_mesh + 8 * pi ^ 3 * 1i * x_radial_freq_mesh .^ 3 / 6 ;
    % gradient_kernels_dft( 3, :, :, : ) =  2 * pi * 1i * z_radial_freq_mesh + 8 * pi ^ 3 * 1i * z_radial_freq_mesh .^ 3 / 6 ;

    % gradient_kernels_dft( 1, :, :, : ) =  1i * sin( 2 * pi * y_radial_freq_mesh );
    % gradient_kernels_dft( 2, :, :, : ) =  1i * sin( 2 * pi * x_radial_freq_mesh );
    % gradient_kernels_dft( 3, :, :, : ) =  1i * sin( 2 * pi * z_radial_freq_mesh );

    % % neighboring voxel derivatives:
    % !!!!! are these derivatives off by a factor 2 ??????? SAM 6/3/22 !!!!!!!!!!!!!!!!!!!
    gradient_kernels_dft( 1, :, :, : ) =  1i * derivative_weights_from_blurring( 1 ) * sin( 2 * pi * y_pixel_freq_mesh ) / 2 ;
    gradient_kernels_dft( 2, :, :, : ) =  1i * derivative_weights_from_blurring( 2 ) * sin( 2 * pi * x_pixel_freq_mesh ) / 2 ;
    gradient_kernels_dft( 3, :, :, : ) =  1i * derivative_weights_from_blurring( 3 ) * sin( 2 * pi * z_pixel_freq_mesh ) / 2 ;

    % gradient_kernels_dft( 1, :, :, : ) = 1i * sin( 2 * pi * y_pixel_freq_mesh );
    % gradient_kernels_dft( 2, :, :, : ) = 1i * sin( 2 * pi * x_pixel_freq_mesh );
    % gradient_kernels_dft( 3, :, :, : ) = 1i * sin( 2 * pi * z_pixel_freq_mesh );

    % gradient_kernels_dft( 1, :, :, : ) = 1i * sin( 2 * pi * y_pixel_freq_mesh / 2 );
    % gradient_kernels_dft( 2, :, :, : ) = 1i * sin( 2 * pi * x_pixel_freq_mesh / 2 );
    % gradient_kernels_dft( 3, :, :, : ) = 1i * sin( 2 * pi * z_pixel_freq_mesh / 2 );

    % gradient_kernels_dft( 1, :, :, : ) = 1 - exp( - 1i * pi * y_pixel_freq_mesh );
    % gradient_kernels_dft( 2, :, :, : ) = 1 - exp( - 1i * pi * x_pixel_freq_mesh );
    % gradient_kernels_dft( 3, :, :, : ) = 1 - exp( - 1i * pi * z_pixel_freq_mesh );

    % gradient_kernels_dft( 1, :, :, : ) = 1 - exp( - 1i * 2 * pi * y_pixel_freq_mesh );
    % gradient_kernels_dft( 2, :, :, : ) = 1 - exp( - 1i * 2 * pi * x_pixel_freq_mesh );
    % gradient_kernels_dft( 3, :, :, : ) = 1 - exp( - 1i * 2 * pi * z_pixel_freq_mesh );

%     compute spatial convolultions of derivative and gaussian kernels in frequency space with .*
    curvatures_gaussian_kernel_dft = curvatures_kernels_dft .* reshape( matching_kernel_dft, [ 1, size_of_chunk_dft ]);
      gradient_gaussian_kernel_dft =   gradient_kernels_dft .* reshape( matching_kernel_dft, [ 1, size_of_chunk_dft ]);

end

if is_only_computing_Laplacian
    
    Laplacian_matching_kernel_dft = Laplacian_kernel_dft .* matching_kernel_dft ;

end

% if ~ is_only_computing_Laplacian
% 
%     % % normalize each derivative kernel by its variance about zero
%     % curvatures_gaussian_kernel_dft = numel_chunk * radius_of_lumen_in_microns ^ 2 * curvatures_gaussian_kernel_dft ./ sum( curvatures_gaussian_kernel_dft( 1 : 6, : ) .^ 2, 2 ) .^ 0.5 ;
%     %   gradient_gaussian_kernel_dft = numel_chunk * radius_of_lumen_in_microns     *   gradient_gaussian_kernel_dft ./ sum(   gradient_gaussian_kernel_dft(   :  , : ) .^ 2, 2 ) .^ 0.5 ;
% 
%     % curvatures_gaussian_kernel_dft( 1, :, :, : ) = radius_of_lumen_in_voxels( 1 ) ^ 2                              * curvatures_gaussian_kernel_dft( 1, :, :, : ) / sum( curvatures_gaussian_kernel_dft( 1, : ) .^ 2 ) .^ 0.5 ;
%     % curvatures_gaussian_kernel_dft( 2, :, :, : ) = radius_of_lumen_in_voxels( 2 ) ^ 2                              * curvatures_gaussian_kernel_dft( 2, :, :, : ) / sum( curvatures_gaussian_kernel_dft( 2, : ) .^ 2 ) .^ 0.5 ;
%     % curvatures_gaussian_kernel_dft( 3, :, :, : ) = radius_of_lumen_in_voxels( 3 ) ^ 2                              * curvatures_gaussian_kernel_dft( 3, :, :, : ) / sum( curvatures_gaussian_kernel_dft( 3, : ) .^ 2 ) .^ 0.5 ;
%     % curvatures_gaussian_kernel_dft( 4, :, :, : ) = radius_of_lumen_in_voxels( 1 ) * radius_of_lumen_in_voxels( 2 ) * curvatures_gaussian_kernel_dft( 4, :, :, : ) / sum( curvatures_gaussian_kernel_dft( 4, : ) .^ 2 ) .^ 0.5 ;
%     % curvatures_gaussian_kernel_dft( 5, :, :, : ) = radius_of_lumen_in_voxels( 2 ) * radius_of_lumen_in_voxels( 3 ) * curvatures_gaussian_kernel_dft( 5, :, :, : ) / sum( curvatures_gaussian_kernel_dft( 5, : ) .^ 2 ) .^ 0.5 ;
%     % curvatures_gaussian_kernel_dft( 6, :, :, : ) = radius_of_lumen_in_voxels( 3 ) * radius_of_lumen_in_voxels( 1 ) * curvatures_gaussian_kernel_dft( 6, :, :, : ) / sum( curvatures_gaussian_kernel_dft( 6, : ) .^ 2 ) .^ 0.5 ;
%     % 
%     %   gradient_gaussian_kernel_dft( 1, :, :, : ) = radius_of_lumen_in_voxels( 1 )                                  *   gradient_gaussian_kernel_dft( 1, :, :, : ) / sum(   gradient_gaussian_kernel_dft( 1, : ) .^ 2 ) .^ 0.5 ;
%     %   gradient_gaussian_kernel_dft( 2, :, :, : ) = radius_of_lumen_in_voxels( 2 )                                  *   gradient_gaussian_kernel_dft( 2, :, :, : ) / sum(   gradient_gaussian_kernel_dft( 2, : ) .^ 2 ) .^ 0.5 ;
%     %   gradient_gaussian_kernel_dft( 3, :, :, : ) = radius_of_lumen_in_voxels( 3 )                                  *   gradient_gaussian_kernel_dft( 3, :, :, : ) / sum(   gradient_gaussian_kernel_dft( 3, : ) .^ 2 ) .^ 0.5 ;
% 
% end

if is_only_computing_Laplacian

    Laplacian_chunk_dft = Laplacian_matching_kernel_dft .* chunk_dft ;

    Laplacian_chunk = ifftn( Laplacian_chunk_dft, 'symmetric' );
    
end

if ~ is_only_computing_Laplacian

    curvatures_chunk_dft = curvatures_gaussian_kernel_dft .* reshape( chunk_dft, [ 1, size_of_chunk_dft ]);
      gradient_chunk_dft =   gradient_gaussian_kernel_dft .* reshape( chunk_dft, [ 1, size_of_chunk_dft ]);
    
    curvatures_chunk = zeros( size( curvatures_chunk_dft ));
      gradient_chunk = zeros( size(   gradient_chunk_dft ));
    
    if is_doing_eigenvalue_decomp, curvature_index_range = 1 : 6 ; 
    else,                          curvature_index_range = 1 : 3 ; end
      
    % inverse fourier transforms to compute curvautures and gradients
    for curvature_index = curvature_index_range
    
        curvatures_chunk( curvature_index, :, :, : ) ...
            = ifftn( curvatures_chunk_dft( curvature_index, :, :, : ), 'symmetric' );
    
    end

    for gradient_index = 1 : 3
        
        gradient_chunk( gradient_index, :, :, : )    ...
            = ifftn(   gradient_chunk_dft(  gradient_index, :, :, : ), 'symmetric' );
    
    end

    curvatures_chunk = curvatures_chunk(        :         , ...
                                         local_ranges{ 1 }, ...
                                         local_ranges{ 2 }, ...
                                         local_ranges{ 3 }  );

      gradient_chunk =   gradient_chunk(        :         , ...
                                         local_ranges{ 1 }, ...
                                         local_ranges{ 2 }, ...
                                         local_ranges{ 3 }  );

%     gradient_direction = gradient_chunk ./ sum( gradient_chunk .^ 2, 1 ) .^ 0.5 ;

    if is_doing_eigenvalue_decomp

        size_of_chunk = cellfun( 'length', local_ranges );                
        
        number_of_voxels_in_chunk = numel( curvatures_chunk( 1, :, :, : ));

        voxel_index_full_range = 1 : number_of_voxels_in_chunk ;

        % Laplacian is the trace of the Hessian
        Laplacian_chunk = sum( curvatures_chunk(( 1 : 3 )' + 6 * ( voxel_index_full_range - 1 )), 1 );
        
        energy_chunk = reshape( Laplacian_chunk, size_of_chunk );        

    end
    
    if ~ is_doing_eigenvalue_decomp
        
        Laplacian_chunk = squeeze( sum( curvatures_chunk, 1 ));
        
    end
end

if is_only_computing_Laplacian
    
    energy_chunk = Laplacian_chunk( local_ranges{ 1 }, ...
                                    local_ranges{ 2 }, ...
                                    local_ranges{ 3 }  );
        
end

if ~ is_only_computing_Laplacian

    valid_voxels = Laplacian_chunk < 0 ;
    
    % Before 6/5/19
    % 
    % SAM 4/1/19
    % pixels_per_radius = radius_of_lumen_in_microns ./ microns_per_pixel ;
    % 
    % symmetry_ratio_factor_vector = symmetry_ratio_factor * pixels_per_radius' ;

    % % Before: 4/1/19 SAM 
    % symmetry_ratio_factor_vector = symmetry_ratio_factor ./ microns_per_pixel' * microns_per_pixel( 1 );                                  

    if is_doing_eigenvalue_decomp
            
%         symmetry_ratio_factor_vector = symmetry_ratio_factor * [ 1; 1; 1 ];           
        
        curvature_indices_in_hessian = [ 1, 4, 6; ...
                                         4, 2, 5; ...
                                         6, 5, 3  ];

        % only compute principal curvatures where the laplacian is negative, set the rest to infinity:                                                              
        voxel_index_range = voxel_index_full_range( valid_voxels );

        curvatures_chunk_columnated                                                                    ...
                              = curvatures_chunk(   curvature_indices_in_hessian                       ...
                                                  + 6 * ( permute( voxel_index_range, [ 1, 3, 2 ]) - 1 ));

        [ principal_curvature_vectors_cell, principal_curvature_values_cell ] = cellfun( @eig, num2cell( curvatures_chunk_columnated, [ 1, 2 ]), 'UniformOutput', false );

        principal_curvature_vectors = zeros([ 3, 3, size_of_chunk ]);
        principal_curvature_values  = zeros([ 3, 3, size_of_chunk ]);

        principal_curvature_vectors( :, :, voxel_index_range ) = cell2mat( principal_curvature_vectors_cell );
        principal_curvature_values(  :, :, voxel_index_range ) = cell2mat( principal_curvature_values_cell  );

        principal_curvature_values = reshape( principal_curvature_values, [ 9, size_of_chunk ]);

%         % % sum up each principal curvature weighted by the exponent of the negative projection of the
%         % % gradient onto that principal curvature direction( - abs( cos( theta )) * magnitude of gradient
%         % % / curvature value )
%         % energy_chunk                                                                                                        ...
%         %     = squeeze( sum(    principal_curvature_values([ 1, 5, 9 ], :, :, : )                                            ...
%         %                     .* exp( - abs(    squeeze( sum( reshape( gradient_chunk, [ 3, 1, size_of_chunk ])               ...
%         %                                                     .* principal_curvature_vectors                                  ...
%         %                                                     .* symmetry_ratio_factor_vector                  , 1 ))         ...
%         %                                       ./ principal_curvature_values([ 1, 5, 9 ], :, :, : )                  )), 1 ));
% 
%         % sum up each principal curvature weighted by the exponent of the negative projection of the
%         % gradient onto that principal curvature direction( - abs( cos( theta )) * magnitude of gradient
%         % / curvature value ).  
%         % 
%         % Then multiply by a shape factor ( sum( curvatures ) ^ 2 / sum( curvatures ^ 2 )) SAM 4/2/19
%         % energy_chunk                                                                                                        ...
%         %     = squeeze( sum(    principal_curvature_values([ 1, 5, 9 ], :, :, : )                                            ...
%         %                     .* exp( - abs(    squeeze( sum( reshape( gradient_chunk, [ 3, 1, size_of_chunk ])               ...
%         %                                                     .* principal_curvature_vectors                                  ...
%         %                                                     .* symmetry_ratio_factor_vector                  , 1 ))         ...
%         %                                       ./ principal_curvature_values([ 1, 5, 9 ], :, :, : )                  )), 1 ))...
%         %     .* reshape( laplaican_chunk, size_of_chunk ) .^ 2                                                               ...
%         %     ./ squeeze( sum( principal_curvature_values([ 1, 5, 9 ], :, :, : ) .^ 2, 1 ));
% 
%         % principal_energy_values                                                                     ...
%         %     =    principal_curvature_values([ 1, 5, 9 ], :, :, : )                                  ...
%         %     .* exp( - abs(    squeeze( sum( reshape( gradient_chunk, [ 3, 1, size_of_chunk ])       ...
%         %                                     .* principal_curvature_vectors                          ...
%         %                                     .* symmetry_ratio_factor_vector                  , 1 )) ...
%         %                       ./ principal_curvature_values([ 1, 5, 9 ], :, :, : )                  ));
% 
% %         % square instead of absolute value in the exponent argument
% %         principal_energy_values                                                                         ...
% %             =    principal_curvature_values([ 1, 5, 9 ], :, :, : )                                      ...
% %             .* exp( - (    squeeze( sum( reshape( gradient_chunk, [ 3, 1, size_of_chunk ])              ...
% %                                          .* principal_curvature_vectors                                 ...
% %                                          .* symmetry_ratio_factor_vector                  , 1 ))        ...
% %                            ./ principal_curvature_values([ 1, 5, 9 ], :, :, : )                  ) .^ 2 );

%         % caught math-concept mistake % SAM 12/20/21, results in smoother output to do L2 norm (but preserving signs) of
%         % components instead of L1 as above during this sum( ___, 1 ) call
%         principal_energy_values                                                                         ...
%    ...         =       principal_curvature_values([ 1, 5, 9 ], :, :, : ) .^ 2                             ...
%       ...     .* sign( principal_curvature_values([ 1, 5, 9 ], :, :, : ))                                 ...
% ...            =        abs( principal_curvature_values([ 1, 5, 9 ], :, :, : ))                                   ...
% ...           .*             principal_curvature_values([ 1, 5, 9 ], :, :, : )                                    ...
%             =                principal_curvature_values([ 1, 5, 9 ], :, :, : )                               ...
%            .* exp( - abs(    squeeze( sum(  reshape( abs( gradient_chunk )                               ...
%                                                       .*  gradient_chunk   , [ 3, 1, size_of_chunk ])   ...
%                                          .* abs( principal_curvature_vectors )                          ... %%% !!!! why the abs() call. why choose the dot product of the squares over the square of the dot product SAM 6/3/22 ??????
%                                          .*      principal_curvature_vectors  ,                         ...
%                            ...               .* symmetry_ratio_factor_vector                            , ...
%                                          1                                                           )) ...
%                           ./ principal_curvature_values([ 1, 5, 9 ], :, :, : ) .^ 2 / 2 ... % factor of 1/2 added 1/14/22 ... % and removed 1/14/22 % * 2                      ... % factor of 2 added 1/13/22
%                      )); ...    

       principal_energy_values                                                                               ...
           = principal_curvature_values([ 1, 5, 9 ], :, :, : )                                               ...
          .* exp( - ( squeeze(     sum(    reshape( gradient_chunk,     [ 3, 1, size_of_chunk ])             ... !!! this is how it should look 6/3/22 SAM<  
                             	       .* principal_curvature_vectors,    1                     ))           ...
                               ./ principal_curvature_values([ 1, 5, 9 ], :, :, : )               ) .^ 2 / 2 ); % factor = 1 % factor of 2 inside square: sigma is a P.Graident/P.Curvature of 0.5
                     
                     
%         % % only the least two principal curvatures are considered
%         % principal_energy_values                                                                     ...
%         %     =     principal_curvature_values([ 1, 5 ], :, :, : )                                    ...
%         %     .* exp( - abs(    squeeze( sum( reshape( gradient_chunk, [ 3, 1, size_of_chunk ])       ...
%         %                                     .* principal_curvature_vectors( :, 1 : 2, :, :, : )     ...
%         %                                     .* symmetry_ratio_factor_vector                  , 1 )) ...
%         %                       ./ principal_curvature_values([ 1, 5 ], :, :, : )                     ));
% 
%         % % square instead of absolute value in the exponent argument
%         % % only the least two principal curvatures are considered
%         % principal_energy_values                                                                         ...
%         %     =    principal_curvature_values([ 1, 5 ], :, :, : )                                      ...
%         %     .* exp( - (    squeeze( sum( reshape( gradient_chunk, [ 3, 1, size_of_chunk ])              ...
%         %                                  .* principal_curvature_vectors( :, 1 : 2, :, :, : )            ...
%         %                                  .* symmetry_ratio_factor_vector                  , 1 ))        ...
%         %                    ./ principal_curvature_values([ 1, 5 ], :, :, : )                     ) .^ 2 );
% 
% 
%         % % only the least principal curvature is considered
%         % principal_energy_values                                                                     ...
%         %     =     principal_curvature_values( 1, :, :, : )                                    ...
%         %     .* exp( - abs(    squeeze( sum( reshape( gradient_chunk, [ 3, 1, size_of_chunk ])       ...
%         %                                     .* principal_curvature_vectors( :, 1, :, :, : )     ...
%         %                                     .* symmetry_ratio_factor_vector                  , 1 )) ...
%         %                       ./ principal_curvature_values( 1, :, :, : )                     ));
% 
% 
%         % % shape factor is also considered
%         %     shape_factor_exponent = 2 ;
%         % 
%         % laplacian_factor_exponent = 1 ;
%         %                   
%         % energy_chunk = squeeze(  - abs( sum( principal_energy_values, 1 )) .^ ( shape_factor_exponent + laplacian_factor_exponent ) ...
%         %                         ./ sum( abs( principal_energy_values ) .^ shape_factor_exponent, 1 )        );
% 
%         % % remove points that dont have a negative energy before forcing the third component to be negative
%         % principal_energy_values( 1, squeeze( sum( principal_energy_values( :, : ), 1 )) >= 0 ) = Inf ;
%         % 
%         % % forcing the third component to be negative (a tube whose axial trace is a local min of intensity
%         % % (a tube between bright spots) will appear the same as one whose axial trace is a max (a bright
%         % % blob).
%         % principal_energy_values( 3, : ) = - abs( principal_energy_values( 3, : ));

        % do not hurt the objective function with unfavorable sign in the third principal component,
        % because this will be unfavorable wherever dimmer vessel approaches brighter.
        principal_energy_values( 3, principal_energy_values( 3, : ) > 0 ) = 0 ; % SAM 12/15/21   % SAM 12/23/21
      
        energy_chunk = squeeze( sum( principal_energy_values, 1 )); % previous best, noisy around bifurcations, oversizing consistently SAM 12/15/21 % SAM 12/23/21

%         % only keep least two principal curvatures
%         energy_chunk = squeeze( sum( principal_energy_values( 1 : 2, :, :, : ), 1 )); % SAM 12/14/21 % SAM 12/18/21 % SAM 12/20/21
        
%         % only keep least principal curvature
%         energy_chunk = squeeze( principal_energy_values( 1, :, :, : )); % SAM 12/18/21  % SAM 12/19/21 Joe's bday !

%         % only keep middle principal curvature
%         energy_chunk = squeeze( principal_energy_values( 2, :, :, : )); % SAM 12/19/21

%         % weighting first principal curvature most, then second, then third % SAM 12/15/21 % SAM 12/23/21
%         energy_chunk = squeeze( sum([ 1.5; 1; 1/2 ] .* principal_energy_values, 1 )); 

        
    %     % least two principal curvatures are subject to 1st derivative symmetry adjustment
    %     energy_chunk = squeeze( sum( principal_energy_values( 1 : 2, :, :, : ), 1 ) + principal_curvature_values( 9, :, :, : ));

        % only the least principal curvature is subject to 1st derivative symmetry adjustment
    %     energy_chunk = squeeze( principal_energy_values( 1, :, :, : ) + sum( principal_curvature_values([ 5, 9 ], :, :, : ), 1 ));

        % % include shape factor: 1 for two eignevalues have all the weight, 0 for one or three holding weight
        % shape_factor_chunk = 1 - abs( squeeze(    sum( principal_energy_values, 1 ) .^ 2          ...
        %                                        ./ sum( principal_energy_values .^ 2, 1   )) - 2 ) ;

    % %     include shape factor: 1 for two eignevalues have equal weight, 0 for one holding all the weight
    %     shape_factor_chunk = abs( squeeze(    sum( principal_energy_values( 1 : 2, :, :, : ), 1 ) .^ 2          ...
    %                                        ./ sum( principal_energy_values( 1 : 2, :, :, : ) .^ 2, 1   )) - 1 ) ;
    %     shape_factor_chunk = 1 - abs( squeeze(    sum( principal_energy_values( 1 : 2, :, :, : ), 1 ) .^ 2               ...
    %                                            ./ sum( principal_energy_values( 1 : 2, :, :, : ) .^ 2, 1   )) - 2 );

    %     shape_factor_chunk = abs( squeeze(    min( principal_curvature_values( 5, :, :, : ), 0 ) ...
    %                                        ./      principal_curvature_values( 1, :, :, : )      ));

    %     shape_factor_chunk = abs( squeeze(    min( sum( principal_curvature_values([ 5, 9 ], :, :, : ), 1 ), 0 ) ...
    %                                        ./           principal_curvature_values( 1, :, :, : )                 ));

    %     shape_factor_chunk = abs( squeeze(    sum( min( principal_curvature_values([ 5, 9 ], :, :, : ), 0 ), 1 ) ...
    %                                        ./           principal_curvature_values( 1, :, :, : )                 ));

    % 	shape_factor_chunk( shape_factor_chunk == Inf ) = 0 ;
    %     
    %     energy_chunk = energy_chunk .* shape_factor_chunk ;
    
    end
    
    if ~ is_doing_eigenvalue_decomp
        
        energy_chunk                                                         ...
            = squeeze( sum(    curvatures_chunk                              ...
                            .* exp( - (    symmetry_ratio_factor             ...
                                         * gradient_chunk                    ...
                                        ./ curvatures_chunk      ) .^ 2 ), 1 ));
                               
    end

    % % only least principal curvature is considered
    % energy_chunk                                                                                  ...
    %     = squeeze( principal_curvature_values( 1, :, :, : ))                                      ...
    %       .* exp( - abs(    squeeze( sum( reshape( gradient_chunk, [ 3, 1, size_of_chunk ])       ...
    %                                       .* principal_curvature_vectors( :, 1, :, :, : )         ...
    %                                       .* symmetry_ratio_factor_vector                  , 1 )) ...
    %                                  ./ squeeze( principal_curvature_values( 1, :, :, : ))        ));

    % % only the least two principal curvatures are considered
    % energy_chunk                                                                                                        ...
    %     = squeeze( sum(    principal_curvature_values([ 1, 5 ], :, :, : )                                               ...
    %                     .* exp( - abs(    squeeze( sum( reshape( gradient_chunk, [ 3, 1, size_of_chunk ])               ...
    %                                                     .* principal_curvature_vectors( :, 1 : 2, :, :, : )             ...
    %                                                     .* symmetry_ratio_factor_vector                  , 1 ))         ...
    %                                       ./ principal_curvature_values([ 1, 5 ], :, :, : )                     )), 1 ));
    % 

    % the voxels that are not valid will have an Inf for their energy value
    energy_chunk( ~ valid_voxels ) = Inf ;

    energy_chunk( isnan( energy_chunk )) = Inf ;
    
end % ELSE energy function

energy_chunk( energy_chunk >= 0 ) = Inf ;

end % function
