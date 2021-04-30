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

vertices_image = zeros( size_of_image );

number_of_image_voxels = prod( size_of_image );

vertex_scale_subscripts = uint8( vertex_scale_subscripts );

number_of_vertices = length( vertex_scale_subscripts );

% pre-calculating all of the ellipsoidal structuring elements to be used to paint in the image
number_of_scales = size( pixels_per_radius_range, 1 );

structuring_element_linear_indexing = cell( number_of_scales, 1 );

for scale_index = 1 : number_of_scales

    % find all pixel locations within the ellipsoid radii from the vertex position
    structuring_element_linear_indexing{ scale_index }                                                 ...
        = construct_structuring_element_V190( pixels_per_radius_range( scale_index, : ), size_of_image );
        
end % scale FOR

% loop through the vertices
for vertex_index = 1 : number_of_vertices
    
        % find the linear index of the center pixel of the current vertex
        vertex_position_linear_index = sub2ind( size_of_image,                              ...
                                                vertex_space_subscripts( vertex_index, 1 ), ...
                                                vertex_space_subscripts( vertex_index, 2 ), ...
                                                vertex_space_subscripts( vertex_index, 3 )  );        
        
        vertices_image( max( min( vertex_position_linear_index                                                        ...
                                     + structuring_element_linear_indexing{ vertex_scale_subscripts( vertex_index )}, ...
                                  number_of_image_voxels                                            ),                ...
                             1                                                                         ))             ...
                                                                                    = - vertex_energies( vertex_index );        
end

mat2tif( uint16( vertices_image ), visual_file )

end