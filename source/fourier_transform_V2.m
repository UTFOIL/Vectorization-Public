function [ image_dft ] = fourier_transform_V2( image )
% SAM 12/6/17, to be used with gaussian_filter_V3
%
% V2, in which the image dimensions are forced to be even (in contrast to
% being forced to a power of 2 in V1).  SAM 12/15/17

size_of_image = size( image );

next_even_image_size = 2 * ceil(( size_of_image + 1 ) / 2 );

padded_image = padarray( image, next_even_image_size - size_of_image, 'symmetric', 'post' );

image_dft = fftn( padded_image );

end