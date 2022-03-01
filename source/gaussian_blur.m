function [ blurred_image ] = gaussian_blur( image, pixels_per_sigma_gaussian )
%% gaussian_blur
% SAM 4/8/19

size_of_image = size( image );

% pad_lengths = 2 * ceil( pixels_per_sigma_gaussian * 3 / 2 ) + mod( size_of_image, 2 );
% 
% padded_image = padarray( image, pad_lengths, mean( image( : )), 'post' ); % changed padding value from default (0) to the mean image value. SAM 210909

% changed padding to do mean value by slice SAM 210914
pixels_per_sigma_gaussian( 3 ) = pixels_per_sigma_gaussian( 3 ) * 2 ; % do temporary change to z-Gaussian. put padding above and below in z

pad_lengths = 2 * ceil( pixels_per_sigma_gaussian * 3 / 2 ) + mod( size_of_image, 2 );

pixels_per_sigma_gaussian( 3 ) = pixels_per_sigma_gaussian( 3 ) / 2 ; % undo temporary change

padded_image = zeros( size_of_image + pad_lengths );

z_pad_length_pre  =  ceil( pad_lengths( 3 ) / 2 );
z_pad_length_post = floor( pad_lengths( 3 ) / 2 );

padded_image( :, :,   1 : z_pad_length_pre            ) = median( image( :, :,  1  ), 'all' );
padded_image( :, :, end - z_pad_length_post + 1 : end ) = median( image( :, :, end ), 'all' );

z_range = 1 : size_of_image( 3 );

for z_slice = z_range
    
                     padded_image( :, :, z_slice + z_pad_length_pre )              ...
        = padarray(         image( :, :, z_slice ),  pad_lengths([ 1, 2 ]),        ...
                    median( image( :, :, z_slice ), 'all' ),                'post' );
    
end

image_dft = fftn( padded_image );

size_of_image_dft = size_of_image + pad_lengths ;

% mesh generation for building the gaussian blurring filter in the Fourier domain
[ y_pixel_freq_mesh, x_pixel_freq_mesh, z_pixel_freq_mesh ]                                                       ...
    = ndgrid([ 0 : size_of_image_dft( 1 ) / 2 - 1, - size_of_image_dft( 1 ) / 2 : - 1 ] / size_of_image_dft( 1 ), ...
             [ 0 : size_of_image_dft( 2 ) / 2 - 1, - size_of_image_dft( 2 ) / 2 : - 1 ] / size_of_image_dft( 2 ), ...
             [ 0 : size_of_image_dft( 3 ) / 2 - 1, - size_of_image_dft( 3 ) / 2 : - 1 ] / size_of_image_dft( 3 )  );

sigma_freq_mesh_squared = ( y_pixel_freq_mesh * pixels_per_sigma_gaussian( 1 )) .^ 2 ...
                        + ( x_pixel_freq_mesh * pixels_per_sigma_gaussian( 2 )) .^ 2 ...
                        + ( z_pixel_freq_mesh * pixels_per_sigma_gaussian( 3 )) .^ 2 ;
                           
gaussian_kernel_dft = exp( - pi ^ 2 * 2 * ( sigma_freq_mesh_squared ));

%         % uncomment to inspect the gaussian blurring kernel in the spatial domain
% 
%         gaussian_kernel_image = ifftn( gaussian_kernel_dft, 'symmetric' );
% 
%         sum( gaussian_kernel_image( : ))

blurred_image_dft = image_dft .* gaussian_kernel_dft ;

% inverse fourier transforms to compute blurring operation
padded_blurred_image = ifftn( blurred_image_dft, 'symmetric' );

blurred_image = padded_blurred_image( 1 : size_of_image( 1 ), ...
                                      1 : size_of_image( 2 ), ...
                 z_pad_length_pre + ( 1 : size_of_image( 3 )) );

end
