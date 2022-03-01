function [ vertex_center_image, vertex_locations, reading_box_apothem, strel_linear_LUT,    ...
                   max_edge_length_in_microns_range, numel_of_strel,             ...
            cum_prod_image_dims, size_of_image, path_to_energy_data, strel_distance_LUT ]            ...
  = generate_reference_image(  lumen_radius_in_microns_range, microns_per_voxel, ...
                               vertex_space_subscripts, ...
                               strel_apothem, max_edge_length_per_origin_radius, ...
                                                  data_directory, energy_handle  )


path_to_energy_data = [ data_directory, energy_handle ];

energy_file_info    = h5info( path_to_energy_data );

image_pair_dims     = energy_file_info.Datasets.Dataspace.Size ;

size_of_image       = int64( image_pair_dims( 1 : 3 ));

[ strel_linear_LUT, numel_of_strel, cum_prod_image_dims, local_subscripts_range ] = calculate_linear_strel( size_of_image, strel_apothem );

[ strel_y,                                     ...
  strel_x,                                     ...
  strel_z  ] = ndgrid( local_subscripts_range, ...
                       local_subscripts_range, ...
                       local_subscripts_range  );

strel_yxz = [ strel_y( : ), strel_x( : ), strel_z( : )];

strel_distance_yxz = double( strel_yxz ).* microns_per_voxel ;

strel_distance_LUT = sum( strel_distance_yxz .^ 2, 2 ) .^ 0.5 ;

% average_radius_in_pixels_range = geomean( lumen_radius_in_pixels_range, 2 );
% 
% max_edge_length_in_scales = max_edge_length_per_radius .* average_radius_in_pixels_range ;

max_edge_length_in_microns_range = max_edge_length_per_origin_radius .* lumen_radius_in_microns_range ;

number_of_vertices  = size( vertex_space_subscripts, 1 ); 

%     case { 'add_vertex_to_edge', 'extend_dead_end_edge' }
%         
%         vertex_unique_range = varargin{ 1 };

% end

% normal small cache is 4 mB, double is 8 B, so max box volume is 80^3 voxels
voxel_characteristic_apothem = 10 ; 

cube_in_microns_for_reading_box_apothem = true ;

if cube_in_microns_for_reading_box_apothem

    reading_box_apothem = max( int64( voxel_characteristic_apothem ./ microns_per_voxel * geomean( microns_per_voxel )'), strel_apothem ); 

else % cube in voxels

    reading_box_apothem = voxel_characteristic_apothem * int64([ 1; 1; 1 ]);

end

% strel_apothem = uint16( strel_apothem ); SAM + WAS 12/5/18

vertex_space_subscripts = int64( vertex_space_subscripts );

vertex_locations =   vertex_space_subscripts( :, 1 )                                  ...
                 + ( vertex_space_subscripts( :, 2 ) - 1 ) * cum_prod_image_dims( 1 ) ...
                 + ( vertex_space_subscripts( :, 3 ) - 1 ) * cum_prod_image_dims( 2 );

% clear( 'vertex_space_subscripts' )

vertex_unique_range = uint32( 1 : number_of_vertices );

vertex_center_image = sparse( double( vertex_locations ), 1, double( vertex_unique_range ), double( cum_prod_image_dims( 3 )), 1, number_of_vertices );


% % spherical structuring element (strel) templates
% number_of_scales = size( lumen_radius_in_microns_range, 1 );
% 
% lumen_radius_in_voxels_range = lumen_radius_in_microns_range ./ microns_per_voxel ;
% 
% dilated_lumen_radius_in_voxels_range = lumen_radius_in_voxels_range * length_dilation_ratio ;
% 
% strel_templates = cell( number_of_scales, 1 );
% 
% for scale_index = 1 : number_of_scales
% 
%     % find all pixel locations within the ellipsoid radii from the vertex position
%     
%     strel_templates{ scale_index }                                                                  ...
%         = construct_structuring_element( dilated_lumen_radius_in_voxels_range( scale_index, : ), ...
%                                               image_dims                                              );
%         
% end % scale FOR


% % find all pixel locations within the ellipsoid radii from the vertex position
% 
% % only the largest subscript template is saved
% [ strel_template, subscript_template ]                                                                 ...
%     = construct_structuring_element( dilated_lumen_radius_in_voxels_range( number_of_scales, : ), ...
%                                           image_dims                                                   );


% numel_of_origin_vertex_strel = round( 4 / 3 * pi * lumen_radius_in_microns_range( end ) .^ 3 / prod( microns_per_voxel ));

% vertex_index_list = [ 1353 ];

end % FUNCTION generate_reference_image