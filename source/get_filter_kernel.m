%function [ Laplacian_matching_kernel_image ]   ...
function [ matching_kernel_image ]...  
         = get_filter_kernel( ...%chunk_dft, matching_kernel_string, radius_of_lumen_in_microns,         ...
                               apothem_per_radius, matching_kernel_string, radius_of_lumen_in_microns, ...
                               vessel_wall_thickness_in_microns, microns_per_pixel,                    ...
                               ...%pixels_per_sigma_PSF, local_ranges, gaussian_to_ideal_ratio, spherical_to_annular_ratio)
                               pixels_per_sigma_PSF, gaussian_to_ideal_ratio, spherical_to_annular_ratio)
                           
% scales_per_octave_radius = scales_per_octave * 3 ;

radius_of_lumen_in_voxels = radius_of_lumen_in_microns ./ microns_per_pixel ;

% mesh generation for building the energy filter in the Fourier domain
% [ size_of_chunk_dft( 1 ), size_of_chunk_dft( 2 ), size_of_chunk_dft( 3 )] = size( chunk_dft );
size_of_chunk_dft = 2 * ceil(radius_of_lumen_in_voxels * apothem_per_radius) + 1;


% size_of_chunk_dft = [20 20 5];


% numel_chunk = numel( chunk_dft );

% [ y_pixel_freq_mesh,                                                                        ...
%   x_pixel_freq_mesh,                                                                        ... 
%   z_pixel_freq_mesh ] = ndgrid(( 0 : size_of_chunk_dft( 1 ) - 1 ) / size_of_chunk_dft( 1 ), ...
%                                ( 0 : size_of_chunk_dft( 2 ) - 1 ) / size_of_chunk_dft( 2 ), ...
%                                ( 0 : size_of_chunk_dft( 3 ) - 1 ) / size_of_chunk_dft( 3 )  );
%    

[ y_pixel_freq_mesh, x_pixel_freq_mesh, z_pixel_freq_mesh ]                                                   ...
    = ndgrid([ 0 : size_of_chunk_dft( 1 ) / 2, - (size_of_chunk_dft( 1 ) - 1) / 2 : - 1 ] / size_of_chunk_dft( 1 ), ...
             [ 0 : size_of_chunk_dft( 2 ) / 2, - (size_of_chunk_dft( 2 ) - 1) / 2 : - 1 ] / size_of_chunk_dft( 2 ), ...
             [ 0 : size_of_chunk_dft( 3 ) / 2, - (size_of_chunk_dft( 3 ) - 1) / 2 : - 1 ] / size_of_chunk_dft( 3 )  );
       
y_micron_freq_mesh = y_pixel_freq_mesh / microns_per_pixel( 1 );
x_micron_freq_mesh = x_pixel_freq_mesh / microns_per_pixel( 2 );
z_micron_freq_mesh = z_pixel_freq_mesh / microns_per_pixel( 3 );

y_radial_freq_mesh = y_micron_freq_mesh * radius_of_lumen_in_microns ;
x_radial_freq_mesh = x_micron_freq_mesh * radius_of_lumen_in_microns ;
z_radial_freq_mesh = z_micron_freq_mesh * radius_of_lumen_in_microns ;
         
y_micron_freq_mesh_squared = y_micron_freq_mesh .^ 2 ;
x_micron_freq_mesh_squared = x_micron_freq_mesh .^ 2 ;
z_micron_freq_mesh_squared = z_micron_freq_mesh .^ 2 ;
         
micron_freq_mesh = (   y_micron_freq_mesh_squared          ...
                     + x_micron_freq_mesh_squared          ...
                     + z_micron_freq_mesh_squared ) .^ 0.5 ;
                     
radial_freq_mesh = micron_freq_mesh * radius_of_lumen_in_microns ;
                 
sigma_PSF_freq_mesh_squared = ( y_pixel_freq_mesh * pixels_per_sigma_PSF( 1 )) .^ 2 ...
                            + ( x_pixel_freq_mesh * pixels_per_sigma_PSF( 2 )) .^ 2 ...
                            + ( z_pixel_freq_mesh * pixels_per_sigma_PSF( 3 )) .^ 2 ;
                     
microns_per_sigma_PSF = pixels_per_sigma_PSF .* microns_per_pixel ;
                                                                 
switch matching_kernel_string
    
    case '3D gaussian'
                  
        % stably deconvolve PSF
        y_radial_freq_mesh_squared = y_micron_freq_mesh .^ 2 * max( radius_of_lumen_in_microns ^ 2 -     microns_per_sigma_PSF( 1 ) ^ 2, microns_per_sigma_PSF( 1 ) ^ 2 );
        x_radial_freq_mesh_squared = x_micron_freq_mesh .^ 2 * max( radius_of_lumen_in_microns ^ 2 -     microns_per_sigma_PSF( 2 ) ^ 2, microns_per_sigma_PSF( 2 ) ^ 2 );
        z_radial_freq_mesh_squared = z_micron_freq_mesh .^ 2 * max( radius_of_lumen_in_microns ^ 2 -     microns_per_sigma_PSF( 3 ) ^ 2, microns_per_sigma_PSF( 3 ) ^ 2 );

        radius_of_lumen_in_voxels = max( radius_of_lumen_in_microns ^ 2 - microns_per_sigma_PSF .^ 2, microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ;

        squared_radial_freq_mesh_gaussian     = (   y_radial_freq_mesh_squared  ...
                                                  + x_radial_freq_mesh_squared  ...
                                                  + z_radial_freq_mesh_squared );
                                      
        matching_kernel_dft = exp( - pi ^ 2 * 2 * squared_radial_freq_mesh_gaussian );
    
    case 'spherical pulse'
        
        radial_angular_freq_mesh = 2 * pi * radial_freq_mesh  ;

        % this one will have a total sum of 1
        matching_kernel_dft = 3 ./ ( radial_angular_freq_mesh .^ 2 )                           ...
                        .* (   sin( radial_angular_freq_mesh ) ./ ( radial_angular_freq_mesh ) ...
                             - cos( radial_angular_freq_mesh )                                 );

        matching_kernel_dft( 1 ) = 1 ;

%         % uncomment to have value inside sphere of approximately 1 (instead of sum to 1)
%         spherical_volume_in_cubic_pxls = 4 / 3 * pi * radius_of_lumen_in_microns ^ 3 / prod( microns_per_pixel );   
% 
%         matching_kernel_dft = matching_kernel_dft * spherical_volume_in_cubic_pxls ;
        
    case '3D gaussian conv spherical pulse'
        Gaussian_lengths = max( microns_per_sigma_PSF, (( gaussian_to_ideal_ratio * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 )) ^ 2 + microns_per_sigma_PSF .^ 2 ) .^ 0.5 );

        % assuming that the squared length of the combined kernel is the sum of the squared lengths
        % of the two kernels being convolved.
        sphere_pulse_lengths_squared = max(( radius_of_lumen_in_microns ^ 2 - Gaussian_lengths .^ 2 + microns_per_sigma_PSF .^ 2 ), 0 );        

        y_radial_freq_mesh_squared = y_micron_freq_mesh .^ 2 * sphere_pulse_lengths_squared( 1 );
        x_radial_freq_mesh_squared = x_micron_freq_mesh .^ 2 * sphere_pulse_lengths_squared( 2 );
        z_radial_freq_mesh_squared = z_micron_freq_mesh .^ 2 * sphere_pulse_lengths_squared( 3 );

        radial_freq_mesh_sphere_pulse = (   y_radial_freq_mesh_squared          ...
                                          + x_radial_freq_mesh_squared          ...
                                          + z_radial_freq_mesh_squared ) .^ 0.5 ;
                                      

        y_radial_freq_mesh = y_micron_freq_mesh * Gaussian_lengths( 1 );
        x_radial_freq_mesh = x_micron_freq_mesh * Gaussian_lengths( 2 );
        z_radial_freq_mesh = z_micron_freq_mesh * Gaussian_lengths( 3 );

%         gaussian_PSF_kernel_dft = ones( size( y_radial_freq_mesh ));
        
        radial_freq_mesh_gaussian     = (   y_radial_freq_mesh .^ 2          ...
                                          + x_radial_freq_mesh .^ 2          ...
                                          + z_radial_freq_mesh .^ 2 ) .^ 0.5 ;
                                                                            
        sigma_freq_mesh          =          radial_freq_mesh_gaussian     ;
        
        radial_angular_freq_mesh = 2 * pi * radial_freq_mesh_sphere_pulse ;
        
        % these will have total sums of 1        
        gaussian_kernel_dft = exp( - pi ^ 2 * 2 * sigma_freq_mesh .^ 2 );

        spherical_pulse_kernel_dft =  ( pi / 2 ./ radial_angular_freq_mesh ) .^ 0.5 ...
                                   .* (   besselj( 2.5, radial_angular_freq_mesh )  ...
                                        + besselj( 0.5, radial_angular_freq_mesh )) ;
                                   
        spherical_pulse_kernel_dft( radial_angular_freq_mesh == 0 ) = 1 ;
       
        matching_kernel_dft = gaussian_kernel_dft .* spherical_pulse_kernel_dft ;
                
        radius_of_lumen_in_voxels = max( Gaussian_lengths .^ 2 - microns_per_sigma_PSF .^ 2, microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ;

%% filters for endogenous signal images (vessel walls (endothelial cells) are lit instead of lumens)

    case '3D gaussian conv annular pulse'
        
        Gaussian_lengths = max( microns_per_sigma_PSF, (( gaussian_to_ideal_ratio * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 )) ^ 2 + microns_per_sigma_PSF .^ 2 ) .^ 0.5 );
                                                          
        % assuming that the squared length of the combined kernel is the sum of the squared lengths
        % of the two kernels being convolved.
        annular_pulse_lengths_squared = max((( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ^ 2 - Gaussian_lengths .^ 2 + microns_per_sigma_PSF .^ 2 ), 0 );        

        sphere_pulse_lengths_squared = max((  radius_of_lumen_in_microns                                          ^ 2 - Gaussian_lengths .^ 2 + microns_per_sigma_PSF .^ 2 ), 0 );        
                               
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

        % these will have total sums of 1        
        gaussian_kernel_dft = exp( - pi ^ 2 * 2 * sigma_freq_mesh .^ 2 );
                              
        % consider integrating (and averaging) cos() over different values for the radius to span
        % the vessel wall (instead of idealizing it as an infinitely thin source at a single radius )
        annular_pulse_kernel_dft = cos( radial_angular_freq_mesh );
        
        matching_kernel_dft = gaussian_kernel_dft .* (( 1 - spherical_to_annular_ratio ) * annular_pulse_kernel_dft + spherical_to_annular_ratio * spherical_pulse_kernel_dft );
        
%         % only count blurring toward derivative weights
%         radius_of_lumen_in_voxels = Gaussian_lengths ./ microns_per_pixel ;

        % only count (non-PSF matching) blurring toward derivative weights
        radius_of_lumen_in_voxels = max( Gaussian_lengths .^ 2 - microns_per_sigma_PSF .^ 2, microns_per_sigma_PSF .^ 2 ) .^ 0.5 ./ microns_per_pixel ;
        
    case 'annular pulse'

        
        % do a difference of spherical pulses each with constant value approximately one inside
        % their spheres. Normalize the difference to have this kernel sum to 1.
        radius_small_in_microns = radius_of_lumen_in_microns ;
        radius_large_in_microns = radius_of_lumen_in_microns + vessel_wall_thickness_in_microns ;

        radial_small_angular_freq_mesh = 2 * pi * micron_freq_mesh * radius_small_in_microns ;
        radial_large_angular_freq_mesh = 2 * pi * micron_freq_mesh * radius_large_in_microns ;


        spherical_volume_small_in_cubic_pxls = 4 / 3 * pi * radius_small_in_microns ^ 3 / prod( microns_per_pixel );                              
        spherical_volume_large_in_cubic_pxls = 4 / 3 * pi * radius_large_in_microns ^ 3 / prod( microns_per_pixel );

        sphere_pulse_small_dft = 3 ./ ( radial_small_angular_freq_mesh .^ 2 )                                     ...
                               .* (   sin( radial_small_angular_freq_mesh ) ./ ( radial_small_angular_freq_mesh ) ...
                               - cos( radial_small_angular_freq_mesh )                                            );

        sphere_pulse_small_dft( 1 ) = 1 ;    

        sphere_pulse_small_dft = sphere_pulse_small_dft * spherical_volume_small_in_cubic_pxls ;

        sphere_pulse_large_dft = 3 ./ ( radial_large_angular_freq_mesh .^ 2 )                                     ...
                               .* (   sin( radial_large_angular_freq_mesh ) ./ ( radial_large_angular_freq_mesh ) ...
                               - cos( radial_large_angular_freq_mesh )                                            );

        sphere_pulse_large_dft( 1 ) = 1 ;    

        sphere_pulse_large_dft = sphere_pulse_large_dft * spherical_volume_large_in_cubic_pxls ;

        % normalize to have sum of 1        
        matching_kernel_dft = sphere_pulse_large_dft - sphere_pulse_small_dft ;

        % comment out these two line to have instead a constant value of 1 (not a total sum of 1).        
        annular_volume_in_cubic_pxls = 4 / 3 * pi                           ...
                                        * (   radius_large_in_microns ^ 3   ...
                                            - radius_small_in_microns ^ 3 ) ...
                                        / prod( microns_per_pixel );        
        
        matching_kernel_dft = matching_kernel_dft / annular_volume_in_cubic_pxls ;
        
    case 'annular pulse V2'
        
        spherical_to_annular_ratio = 0.25 ; % 8/16/18
        
        radius_small_in_microns = radius_of_lumen_in_microns ;
        radius_large_in_microns = radius_of_lumen_in_microns + vessel_wall_thickness_in_microns ;

        radial_small_angular_freq_mesh = 2 * pi * micron_freq_mesh * radius_small_in_microns ;
        radial_large_angular_freq_mesh = 2 * pi * micron_freq_mesh * radius_large_in_microns ;


        spherical_volume_small_in_cubic_pxls = spherical_to_annular_ratio * 4 / 3 * pi * radius_small_in_microns ^ 3 / prod( microns_per_pixel );                              
        spherical_volume_large_in_cubic_pxls =     4 / 3 * pi * radius_large_in_microns ^ 3 / prod( microns_per_pixel );

        % this kernel has average value of 1
        sphere_pulse_small_dft = 3 ./ ( radial_small_angular_freq_mesh .^ 2 )                                     ...
                               .* (   sin( radial_small_angular_freq_mesh ) ./ ( radial_small_angular_freq_mesh ) ...
                               - cos( radial_small_angular_freq_mesh )                                            );

        sphere_pulse_small_dft( 1 ) = 1 ;    

        % this kernel has approxiamtely constant value of A inside the sphere
        sphere_pulse_small_dft = sphere_pulse_small_dft * spherical_volume_small_in_cubic_pxls ;

        % this kernel has average value of 1
        sphere_pulse_large_dft = 3 ./ ( radial_large_angular_freq_mesh .^ 2 )                                     ...
                               .* (   sin( radial_large_angular_freq_mesh ) ./ ( radial_large_angular_freq_mesh ) ...
                               - cos( radial_large_angular_freq_mesh )                                            );

        sphere_pulse_large_dft( 1 ) = 1 ;    

        % this kernel has approxiamtely constant value of 1 inside the sphere        
        sphere_pulse_large_dft = sphere_pulse_large_dft * spherical_volume_large_in_cubic_pxls ;

        % this kernel has max value of 1
        matching_kernel_dft = sphere_pulse_large_dft - sphere_pulse_small_dft ;

        % normalizing to have sum of 1  
        %
        % comment out these two line to have instead a max value of 1 (not a total sum of 1).        
        annular_volume_in_cubic_pxls =    4 / 3 * pi                            ...
                                        * (       radius_large_in_microns ^ 3   ...
                                            - spherical_to_annular_ratio * radius_small_in_microns ^ 3 ) ...
                                        /   prod( microns_per_pixel );        
        
        matching_kernel_dft = matching_kernel_dft / annular_volume_in_cubic_pxls ;
        
    case '3D annular gaussian'
        radius_in_microns = 7.5 ; % 8/17/18 135000

        % assuming that the loss of illumination at the on axis points in the vessel is an
        % exponentially decreasing function of the cross-sectional area of the vessel
        spherical_to_annular_ratio = 1 - exp( - ( radius_of_lumen_in_microns / radius_in_microns ) ^ 2 );
     
        sigma_small_freq_mesh = radial_freq_mesh ;
        
        sigma_large_freq_mesh = micron_freq_mesh                                                  ...
                              * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns ) ;
        
%       % these kernels have average values of 1 
        gaussian_kernel_small_dft = exp( - pi ^ 2 * 2 * sigma_small_freq_mesh .^ 2 );
        gaussian_kernel_large_dft = exp( - pi ^ 2 * 2 * sigma_large_freq_mesh .^ 2 );
        
        max_of_small_gaussian_kernel = ( 2 * pi )  .^ - 0.5 /     radius_of_lumen_in_microns ;
        max_of_large_gaussian_kernel = ( 2 * pi )  .^ - 0.5 / (   radius_of_lumen_in_microns       ...
                                                                + vessel_wall_thickness_in_microns );
                                                    
        % these kernels have max values of ( A and 1 ) * B, for some constant B that depends on the
        % pixel spacing
        gaussian_kernel_small_dft = spherical_to_annular_ratio * gaussian_kernel_small_dft / max_of_small_gaussian_kernel ;
        gaussian_kernel_large_dft = 	gaussian_kernel_large_dft / max_of_large_gaussian_kernel ;
        
        % this kernel has average value of 
        %
        % ( 1 / max_of_large_gaussian_kernel - A / max_of_small_gaussian_kernel )
        matching_kernel_dft = gaussian_kernel_large_dft - gaussian_kernel_small_dft ;
        
        average_value_of_differenced_gaussians = ( 1 / max_of_large_gaussian_kernel - spherical_to_annular_ratio / max_of_small_gaussian_kernel );
        
        % this kernel has average value of 1
        matching_kernel_dft = matching_kernel_dft / average_value_of_differenced_gaussians ;
        
    case '3D annular gaussian V2'
        
        spherical_to_annular_ratio = 0.3 ; % 8/25/18

        sigma_large_freq_mesh = micron_freq_mesh                                                  ...
                              * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns ) ;
        
        sigma_small_freq_mesh = radial_freq_mesh ;                          
                          
        % these kernels have average values of 1 and A
        gaussian_kernel_large_dft = exp( - pi ^ 2 * 2 * sigma_large_freq_mesh .^ 2 );
        gaussian_kernel_small_dft = exp( - pi ^ 2 * 2 * sigma_small_freq_mesh .^ 2 );        
                                                                    
        matching_kernel_dft = ( gaussian_kernel_large_dft - spherical_to_annular_ratio * gaussian_kernel_small_dft )         ...
                                                                                          / ( 1 - spherical_to_annular_ratio );              
    case 'radial gaussian'
        
        spherical_to_annular_ratio = 2 ;
        
        
        % vessel wall is approximated as a radial delta function convolved with a radial Gaussian
        sigma_vessel_wall_freq_mesh = micron_freq_mesh * vessel_wall_thickness_in_microns / 2 ;
                
        radial_lumen_angular_freq_mesh = 2 * pi * micron_freq_mesh * ( radius_of_lumen_in_microns + vessel_wall_thickness_in_microns / 2 ) ;
                
        % these kernels have average values of 1 
        gaussian_kernel_vessel_wall_dft = exp( - pi ^ 2 * 2 * sigma_vessel_wall_freq_mesh .^ 2 );
        
        radial_delta_fxn_lumen_dft      =  sin( radial_lumen_angular_freq_mesh ) ...
                                        ./      radial_lumen_angular_freq_mesh ;
                                    
        spherical_pulse_dft =    3 ./ ( radial_lumen_angular_freq_mesh .^ 2 )                           ...
                            .* (   sin( radial_lumen_angular_freq_mesh ) ./ ( radial_lumen_angular_freq_mesh ) ...
                            -      cos( radial_lumen_angular_freq_mesh )                                 );            
                                    
               spherical_pulse_dft( 1 ) = 1 ;
        radial_delta_fxn_lumen_dft( 1 ) = 1 ;
        
        % their spatial convolution is a radial gaussian function
        matching_kernel_dft =   gaussian_kernel_vessel_wall_dft ...
                            .* ( radial_delta_fxn_lumen_dft + spherical_to_annular_ratio * spherical_pulse_dft ) ...
                            / ( 1 + spherical_to_annular_ratio );    
                                                        
end % matching kernel selection SWITCH

        matching_kernel_dft = matching_kernel_dft .* exp(-1i*y_pixel_freq_mesh*(size_of_chunk_dft(1)-1)/2 ...
                                                       + -1i*x_pixel_freq_mesh*(size_of_chunk_dft(2)-1)/2 ...
                                                       + -1i*z_pixel_freq_mesh*(size_of_chunk_dft(3)-1)/2).^(2*pi);
        
        matching_kernel_image = ifftn( matching_kernel_dft, 'symmetric' );

%         matching_kernel_image( 1 : 5, 1 : 5, 1 : 5 )

%         sum( matching_kernel_image( : ))
                                                                
% blurred_matching_kernel_dft = matching_kernel_dft .* gaussian_PSF_kernel_dft ;

% sum( blurred_matching_kernel_dft( : ) .^ 2 ) ^ 0.5

        % uncomment to inspect the blurred kernel in the spatial domain

%         blurred_matching_kernel_image = ifftn( blurred_matching_kernel_dft, 'symmetric' );

%         blurred_matching_kernel_image( 1 : 5, 1 : 5, 1 : 5 )
        
%         sum( blurred_matching_kernel_image( : ))
             
% laplacian kernel:
% Laplacian_kernel_dft = 2 * radius_of_lumen_in_voxels( 1 ) ^ 2 * ( cos( 2 * pi * y_pixel_freq_mesh ) - 1 ) ...
%                      + 2 * radius_of_lumen_in_voxels( 2 ) ^ 2 * ( cos( 2 * pi * x_pixel_freq_mesh ) - 1 ) ...
%                      + 2 * radius_of_lumen_in_voxels( 3 ) ^ 2 * ( cos( 2 * pi * z_pixel_freq_mesh ) - 1 ) ;
% 
% Laplacian_matching_kernel_dft = Laplacian_kernel_dft .* matching_kernel_dft ;
% 
% Laplacian_matching_kernel_image = ifftn( Laplacian_matching_kernel_dft, 'symmetric' );

end % function
