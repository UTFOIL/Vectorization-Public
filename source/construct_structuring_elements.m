function strel_linear_indexing_templates = construct_structuring_elements( lumen_radius_in_microns_range, microns_per_pixel, size_of_image )

number_of_scales = length( lumen_radius_in_microns_range );

scale_subscript_range = 1 : number_of_scales ;

strel_linear_indexing_templates = cell( number_of_scales, 1 );

radii_in_pixels_range = lumen_radius_in_microns_range ./ microns_per_pixel ;

for scale_subscript = scale_subscript_range

    % find all pixel locations within the ellipsoid radii from the vertex position    
    strel_linear_indexing_templates{ scale_subscript }                                                          ...
        = int64( construct_structuring_element( radii_in_pixels_range( scale_subscript, : ), size_of_image ));

end % constructing relative elements FOR scale

end % FUNCTION structuring_element_linear_indexing_templates
