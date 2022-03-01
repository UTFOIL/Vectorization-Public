function visualize_vertices_V200( vertex_space_subscripts, vertex_scale_subscripts, vertex_energies, ...
                                  pixels_per_radius_range, size_of_image, visual_file )
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
% V190 in which the stucturing elements are computed beforehand into a LUT.  The input positions
% variable is replaced with vertex_space_subscripts and vertex_scale_subscripts variables. SAM
% 7/19/18
%
% V200 eliminating sigma_per_size SAM + WAS 12/5/18
    [ vertices_image ] = paint_vertex_image( vertex_space_subscripts, vertex_scale_subscripts, vertex_energies, ...
                                             pixels_per_radius_range, size_of_image                             );
                                         
    mat2tif( uint16( vertices_image ), visual_file )

end