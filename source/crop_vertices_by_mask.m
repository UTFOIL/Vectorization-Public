function [chosen_vertex_indices] = crop_vertices_by_mask( vertex_space_subscripts, vertex_scale_subscripts, vertex_energies, ...
                                                      pixels_per_radius_range, size_of_image, mask_image)
%crop_vertices_by_mask SAM 4/18/21
%   renders vertices and checks for intersection with the mask. Any vertex intersecting the boundary
%   is marked as zero in the output vector (same length as input vertex list). The mask should have
%   ones where the vertices are not allowed to be.


    is_vertex_in_crop = ones( size( vertex_scale_subscripts ));

    number_of_image_voxels = prod( size_of_image );

    vertex_scale_subscripts = uint8( vertex_scale_subscripts );

    number_of_vertices = length( vertex_scale_subscripts );

    % pre-calculating all of the ellipsoidal structuring elements to be used to paint in the image
    number_of_scales = size( pixels_per_radius_range, 1 );

    structuring_element_linear_indexing = cell( number_of_scales, 1 );

    for scale_index = 1 : number_of_scales

        % find all pixel locations within the ellipsoid radii from the vertex position
        structuring_element_linear_indexing{ scale_index }                                                 ...
            = construct_structuring_element( pixels_per_radius_range( scale_index, : ), size_of_image );

    end % scale FOR

    % loop through the vertices
    for vertex_index = 1 : number_of_vertices

            % find the linear index of the center pixel of the current vertex
            vertex_position_linear_index = sub2ind( size_of_image,                              ...
                                                    vertex_space_subscripts( vertex_index, 1 ), ...
                                                    vertex_space_subscripts( vertex_index, 2 ), ...
                                                    vertex_space_subscripts( vertex_index, 3 )  );        
    
            is_vertex_in_crop( vertex_index )                                                            ...
                = ~ any( mask_image( max( min(   vertex_position_linear_index                     ...
                                                      + structuring_element_linear_indexing{ ...
                                                            vertex_scale_subscripts( vertex_index )},    ...
                                                 number_of_image_voxels                               ), ...
                                            1                                                           )));
    
    end
    
    chosen_vertex_indices = find( is_vertex_in_crop );
    
end

