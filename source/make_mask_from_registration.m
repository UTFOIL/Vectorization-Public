function mask_image = make_mask_from_registration( path_to_registration_txt_file, size_of_image )
% uses the text file output from the registration of tiled stacks to create a mask showing where the
% tiles are not (the parts of the image that are blank). SAM 4/18/21

%             % tiling.registered.txt: 
% # Define the number of dimensions we are working on
% dim = 3
% 
% # Define the image coordinates
% 01.tif; ; (0.0, 0.0, 0.0)
% 02.tif; ; (311.4368796107142, -13.309113355605662, 0.1813471004971551)
% 03.tif; ; (8.421602933444063, 330.40924376814775, 2.7599095302603893)
% 04.tif; ; (317.8987120241055, 304.55879694341024, 2.9836940890739685)

[ starts, dims, number_of_images ] = registration_txt2mat( path_to_registration_txt_file );

% starts( 1, : ) = [0.0, 0.0, 0.0];
% starts( 2, : ) = [311.4368796107142, -13.309113355605662, 0.1813471004971551];
% starts( 3, : ) = [8.421602933444063, 330.40924376814775, 2.7599095302603893];
% starts( 4, : ) = [317.8987120241055, 304.55879694341024, 2.9836940890739685];

starts = 1 - min( starts ) + starts ; % min( starts ) -> [ 1, 1, 1 ] + min( starts )

% !!!!! note that imageJ output x,y,z triplets, but all of the vectorization code uses y,x,z order
starts = starts( :, [ 2, 1, 3 ]);

% encode these numbers in the registration text file and read them in

% dims( 1, : ) = [ 512, 512, 221 ];
% dims( 2, : ) = [ 512, 512, 221 ];
% dims( 3, : ) = [ 512, 512, 221 ];
% dims( 4, : ) = [ 512, 512, 221 ];

mask_image = ones( size_of_image );

for im_idx = 1 : number_of_images

    y_range = ceil( starts( im_idx, 1 )) : floor( starts( im_idx, 1 )) + dims( im_idx, 1 ) - 1 ;
    x_range = ceil( starts( im_idx, 2 )) : floor( starts( im_idx, 2 )) + dims( im_idx, 2 ) - 1 ;
    z_range = ceil( starts( im_idx, 3 )) : floor( starts( im_idx, 3 )) + dims( im_idx, 3 ) - 1 ;

    mask_image( y_range, ...
                x_range, ...
                z_range ) = 0 ;

end


