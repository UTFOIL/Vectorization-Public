function [ strel_linear_LUT_range, numel_of_strel_range, cum_prod_image_dims, local_subscripts_range, strel_distance_LUT_range, strel_unit_vectors_LUT_range ] = calculate_linear_strel_range( size_of_image, microns_per_voxel, lumen_radius_in_microns_range )

number_of_scales = length( lumen_radius_in_microns_range );

scale_subscript_range = 1 : number_of_scales ;

      strel_linear_LUT_range = cell( number_of_scales, 1 );
      local_subscripts_range = cell( number_of_scales, 1 );
    strel_distance_LUT_range = cell( number_of_scales, 1 );
strel_unit_vectors_LUT_range = cell( number_of_scales, 1 );

numel_of_strel_range = zeros( number_of_scales, 1 );

radii_in_pixels_range = max( lumen_radius_in_microns_range ./ microns_per_voxel, 1 );

cum_prod_image_dims = cumprod( size_of_image );

for scale_subscript = scale_subscript_range

    radii = radii_in_pixels_range( scale_subscript, : );
    
    template_center_pxl     =     round( radii ) + 1 ;
    template_pxl_dimensions = 2 * round( radii ) + 1 ;

    [ y_pxl_mesh,                                               ...
      x_pxl_mesh,                                               ...
      z_pxl_mesh  ] = ndgrid( 1 : template_pxl_dimensions( 1 ), ...
                              1 : template_pxl_dimensions( 2 ), ...
                              1 : template_pxl_dimensions( 3 )  );

    voxel_Linf_distances  = max( max( abs( y_pxl_mesh - template_center_pxl( 1 )),  ...
                                      abs( x_pxl_mesh - template_center_pxl( 2 ))), ...
                                      abs( z_pxl_mesh - template_center_pxl( 3 ))  );

                          
    radial_L2_distances_squared = ( y_pxl_mesh - template_center_pxl( 1 )) .^ 2 / radii( 1 ) .^ 2 ...
                                + ( x_pxl_mesh - template_center_pxl( 2 )) .^ 2 / radii( 2 ) .^ 2 ...
                                + ( z_pxl_mesh - template_center_pxl( 3 )) .^ 2 / radii( 3 ) .^ 2 ;
                         
	is_inside_sphere         =  radial_L2_distances_squared <= 1 ;
    is_inside_27_element_box = voxel_Linf_distances         <= 1 ;

    [ absolute_subscript_y,                                              ...
      absolute_subscript_x,                                              ...
      absolute_subscript_z ] = ind2sub( size( y_pxl_mesh ),              ...
                                        find(   is_inside_sphere         ...
                                              | is_inside_27_element_box ));

    local_subscripts_range{ scale_subscript } = [ absolute_subscript_y - template_center_pxl( 1 ), ...
                                                  absolute_subscript_x - template_center_pxl( 2 ), ...
                                                  absolute_subscript_z - template_center_pxl( 3 )  ];

      strel_linear_LUT_range{ scale_subscript }                                    ...
    = local_subscripts_range{ scale_subscript }( :, 1 )                            ...
    + local_subscripts_range{ scale_subscript }( :, 2 ) * cum_prod_image_dims( 1 ) ...
    + local_subscripts_range{ scale_subscript }( :, 3 ) * cum_prod_image_dims( 2 ) ;

    strel_linear_LUT_range{ scale_subscript } = int64( strel_linear_LUT_range{ scale_subscript });

%     [ strel_y,                                     ...
%       strel_x,                                     ...
%       strel_z  ] = ndgrid( local_subscripts_range{ scale_subscript }, ...
%                            local_subscripts_range{ scale_subscript }, ...
%                            local_subscripts_range{ scale_subscript }  );
% 
%     strel_yxz = [ strel_y( : ), strel_x( : ), strel_z( : )];
% 
%     strel_distance_yxz = double( strel_yxz ).* microns_per_voxel ;
% 
%     strel_distance_LUT_range = sum( strel_distance_yxz .^ 2, 2 ) .^ 0.5 ;

    strel_relative_distances = local_subscripts_range{ scale_subscript } .* microns_per_voxel ;

        strel_distance_LUT_range{ scale_subscript } = sum( strel_relative_distances .^ 2, 2 ) .^ 0.5 ;

    strel_unit_vectors_LUT_range{ scale_subscript }                            ...
                                   = strel_relative_distances                  ...
                                  ./ strel_distance_LUT_range{ scale_subscript };
                              
	strel_unit_vectors_LUT_range{ scale_subscript }( isnan( strel_unit_vectors_LUT_range{ scale_subscript })) = 0 ; % replace undefined 0/0 with 0
    
    numel_of_strel_range( scale_subscript ) = numel( strel_distance_LUT_range{ scale_subscript });
    
end % constructing relative elements FOR scale

% strel_length = numel( strel_distance_LUT_range ) .^ ( 1 / 3 );
% 
% strel_distance_LUT_range = reshape( strel_distance_LUT_range, strel_length, strel_length, strel_length );
% 
% microns_per_voxel = [ strel_distance_LUT_range(( end + 1 ) / 2 + 1, ( end + 1 ) / 2    , ( end + 1 ) / 2     ),...
%                       strel_distance_LUT_range(( end + 1 ) / 2    , ( end + 1 ) / 2 + 1, ( end + 1 ) / 2     ),...
%                       strel_distance_LUT_range(( end + 1 ) / 2    , ( end + 1 ) / 2    , ( end + 1 ) / 2 + 1 ) ];
%                   
% strel_displacement_LUT = strel_distance_LUT_range ;     strel_displacement_LUT( 1 : ( end + 1 ) / 2 )...
%                                               = - strel_displacement_LUT( 1 : ( end + 1 ) / 2 );
%                   
% strel_unit_vectors_LUT_range ...
%     = cat( 4, strel_displacement_LUT( :, ( end + 1 ) / 2,    ( end + 1 ) / 2    ) .* ones( size( strel_distance_LUT_range )),...
%               strel_displacement_LUT(    ( end + 1 ) / 2, :, ( end + 1 ) / 2    ) .* ones( size( strel_distance_LUT_range )),...
%               strel_displacement_LUT(    ( end + 1 ) / 2,    ( end + 1 ) / 2, : ) .* ones( size( strel_distance_LUT_range )) );
%           
% strel_unit_vectors_LUT_range =          strel_unit_vectors_LUT_range    ./   strel_distance_LUT_range       ;
% 
% strel_unit_vectors_LUT_range = reshape( strel_unit_vectors_LUT_range, numel( strel_distance_LUT_range ), 3 );
% 
% strel_unit_vectors_LUT_range(( end + 1 ) / 2, : ) = 0 ; % replace undefined 0/0 with 0
% 
% strel_distance_LUT_range = strel_distance_LUT_range( : );


end % FUNCTION
