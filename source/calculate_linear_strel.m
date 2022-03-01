function [ strel_linear_LUT, numel_of_strel, cum_prod_image_dims, local_subscripts_range ] = calculate_linear_strel( size_of_image, strel_apothem )

cum_prod_image_dims = cumprod( size_of_image );

% numel_image         = prod( image_dims( 1 : 3 ));
% numel_image_yx      = prod( image_dims( 1 : 2 ));

% strel_apothem = uint16( strel_apothem );

strel_width   = 2 * strel_apothem + 1 ;

% estimate of number of voxels per edge location:
% numel_of_strel = strel_width ^ 3 ;
numel_of_strel = round( strel_width ^ 3 / 2.5 );

% % strel_dims  = uint8( strel_width * [ 1, 1, 1 ]);
% strel_dims  = [ strel_width, strel_width, strel_width ];
% middle_of_strel = ( numel_of_strel + 1 ) / 2 ;

% % underestimate of number of voxels per edge location:
% numel_of_strel_cross_section = strel_width ^ 2 ;

local_subscripts_range = - strel_apothem : strel_apothem ;

strel_linear_LUT = local_subscripts_range                                                    ;
strel_linear_LUT = local_subscripts_range * cum_prod_image_dims( 1 ) + strel_linear_LUT( : ) ;
strel_linear_LUT = local_subscripts_range * cum_prod_image_dims( 2 ) + strel_linear_LUT( : ) ;
strel_linear_LUT =                                              int64( strel_linear_LUT( : ));

end % FUNCTION
