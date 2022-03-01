function image_matrix = tif2mat( tif_path )
% SAM 9/6/17

info = imfinfo( tif_path );

number_of_slices = numel( info );

slice_height = info( 1 ).Height;
slice_width  = info( 1 ).Width;

bits_per_sample = info( 1 ).BitsPerSample ;

image_matrix = zeros( slice_height, slice_width, number_of_slices );

for slice_index = 1 : number_of_slices
    
    image_matrix( :, :, slice_index ) = imread( tif_path, slice_index, 'Info', info );

end

image_matrix = image_matrix - min( image_matrix( : ));

image_matrix = eval([ 'uint', num2str( bits_per_sample ), '( image_matrix )' ]);

end % FUNCTION 