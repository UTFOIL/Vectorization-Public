function [ linear_indexing, subscripts ] = construct_structuring_element( radii, image_dimensions )
% SAM 9/18/17
%
% V160, the cropping of the template at the find line was modified to allow for images whose size in
% the thid dimension is 1. SAM 5/21/18
%
% V180, the buggy error when the requested element has a dimension that is larger than the image
% dimensions is eliminated by adopting a linear_indexing calculation that does not require the
% object to be contained inside the image.  SAM 5/30/18
%
% V190 only outputs the linear indexing scheme SAM 8/23/18

% template_center_pxl     =     round( radii - 0.5 ) + 1 ;
% template_pxl_dimensions = 2 * round( radii - 0.5 ) + 1 ;

% forced to be centered on a single voxel
template_center_pxl     =     round( radii ) + 1 ;
template_pxl_dimensions = 2 * round( radii ) + 1 ;

[ y_pxl_mesh,                                               ...
  x_pxl_mesh,                                               ...
  z_pxl_mesh  ] = ndgrid( 1 : template_pxl_dimensions( 1 ), ...
                          1 : template_pxl_dimensions( 2 ), ...
                          1 : template_pxl_dimensions( 3 )  );

radial_distances_squared = ( y_pxl_mesh - template_center_pxl( 1 )) .^ 2 / radii( 1 ) .^ 2 ...
                         + ( x_pxl_mesh - template_center_pxl( 2 )) .^ 2 / radii( 2 ) .^ 2 ...
                         + ( z_pxl_mesh - template_center_pxl( 3 )) .^ 2 / radii( 3 ) .^ 2 ;
                     
[ absolute_subscript_y,                                                                     ...
  absolute_subscript_x,                                                                     ...
  absolute_subscript_z ] = ind2sub( size( y_pxl_mesh ), find( radial_distances_squared <= 1 ));

subscripts = [ absolute_subscript_y - template_center_pxl( 1 ), ...
               absolute_subscript_x - template_center_pxl( 2 ), ...
               absolute_subscript_z - template_center_pxl( 3 )  ];

linear_indexing = subscripts( :, 1 )                                                     ...
                + subscripts( :, 2 ) * image_dimensions( 1 )                       ...
                + subscripts( :, 3 ) * image_dimensions( 1 ) * image_dimensions( 2 );
            
end % FUNCTION