function voxel_position = index2position( current_index, cum_prod_image_dims )
    % convert the linear indexing to spatial subscripts                        
        
    current_index = reshape( current_index, 1, numel( current_index ));
    
    pos_xy  = rem( current_index - 1, cum_prod_image_dims( 2 )) + 1 ;

    pos_z   = ( current_index - pos_xy ) / cum_prod_image_dims( 2 ) + 1 ;

    pos_y   = rem( pos_xy - 1, cum_prod_image_dims( 1 )) + 1 ;

    pos_x   = ( pos_xy - pos_y ) / cum_prod_image_dims( 1 ) + 1 ;            

    voxel_position = [ pos_y; ...
                       pos_x; ...
                       pos_z  ];

end
