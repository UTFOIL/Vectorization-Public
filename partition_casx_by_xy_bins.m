function [ point_coordinates_cell, arc_connectivity_cell, arc_diameters_cell ] = partition_casx_by_xy_bins(num_bins_X,num_bins_Y, point_coordinates, arc_connectivity, arc_diameters )

size_of_image = max( point_coordinates ) - min( point_coordinates );

border_locationsX = min( point_coordinates( :, 1 )) + ( 0 : num_bins_X ) * size_of_image( 1 ) / num_bins_X ;
border_locationsY = min( point_coordinates( :, 2 )) + ( 0 : num_bins_Y ) * size_of_image( 2 ) / num_bins_Y ;

point_coordinates_cell = cell( num_bins_X, num_bins_Y );
arc_connectivity_cell  = cell( num_bins_X, num_bins_Y ); 
arc_diameters_cell     = cell( num_bins_X, num_bins_Y );

for bin_idx_X = 1 : num_bins_X
    
    for bin_idx_Y = 1 : num_bins_Y
        
        % isolate points inside the partition
        is_point_in_partition = point_coordinates( :, 1 ) < border_locationsX( 1 + bin_idx_X ) ...
                              & point_coordinates( :, 1 ) > border_locationsX(     bin_idx_X ) ...
                              & point_coordinates( :, 2 ) < border_locationsY( 1 + bin_idx_Y ) ...
                              & point_coordinates( :, 2 ) > border_locationsY(     bin_idx_Y ) ;
                
        point_coordinates_cell{ bin_idx_X, bin_idx_Y } = point_coordinates( is_point_in_partition, : );
                   
        % isolate arcs inside the partition        
        is_arc_in_partition = is_point_in_partition( arc_connectivity( :, 1 )) ...
                            & is_point_in_partition( arc_connectivity( :, 2 )) ;
        
        arc_connectivity_cell{ bin_idx_X, bin_idx_Y } = arc_connectivity( is_arc_in_partition, : );
           arc_diameters_cell{ bin_idx_X, bin_idx_Y } = arc_diameters(    is_arc_in_partition    );
        
        % reassign the point indices inside the arc connectivity matrix, because the list of points
        % and thus their indices have changed
        indices_in_partition = find( is_point_in_partition );
        
        num_points_in_partition = numel( indices_in_partition );
        
        point_indices_tiled2partitioned = sparse( indices_in_partition, ones( num_points_in_partition, 1 ), 1 : num_points_in_partition );

        arc_connectivity_cell{ bin_idx_X, bin_idx_Y } = full( point_indices_tiled2partitioned( arc_connectivity_cell{ bin_idx_X, bin_idx_Y }));  
        
        file_path = [ 'tile_', num2str( bin_idx_X ), '_' num2str( bin_idx_Y )];
        
        % write casx file for this bin in the partition
        casx_mat2file([ file_path, '.casx' ], point_coordinates_cell{ bin_idx_X, bin_idx_Y }, ...
                                               arc_connectivity_cell{ bin_idx_X, bin_idx_Y }, ...
                                                  arc_diameters_cell{ bin_idx_X, bin_idx_Y }  );
                                       
        % write vmv file for this bin in the partition
        [ ~, ~, ~, ~, strand_space_subscripts, strand_scale_subscripts,   ...
                    microns_per_voxel, lumen_radius_in_microns_range, ~ ] ...
                               = casx2strand( point_coordinates_cell{ bin_idx_X, bin_idx_Y }, ...
                                               arc_connectivity_cell{ bin_idx_X, bin_idx_Y }, ...
                                                  arc_diameters_cell{ bin_idx_X, bin_idx_Y }  );
        
        strand_subscripts = cellfun( @( x, y ) [ x, y ], strand_space_subscripts, strand_scale_subscripts, 'UniformOutput', false );
        
        [ point_coordinates2, strand_points ] = strand2vmv( strand_subscripts, microns_per_voxel, lumen_radius_in_microns_range );
        
        vmv_mat2file([ file_path, '.vmv' ], point_coordinates2, strand_points )
        
    end
end

