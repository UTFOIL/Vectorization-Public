%% test network data conversion between strand form and casX.
% SAM 6/28/19

path_to_batch = 'E:\2P imaging\20170802_TxRed_Chronic\Processed_Images_Stack04\PMT01_Red_Images\batch_190531-164319\' ;

path_to_energy_settings = [ path_to_batch, 'settings\energy_190603-172534'                          ];
path_to_edges           = [ path_to_batch, 'vectors\edges_190612-144547_PMT01_Red_Raw_Data_16bit'   ];
path_to_network         = [ path_to_batch, 'vectors\network_190620-224054_PMT01_Red_Raw_Data_16bit' ];

load( path_to_energy_settings );
load( path_to_edges           );
load( path_to_network         );

network_statistics = calculate_network_statistics( strand_subscripts, bifurcation_vertices, lumen_radius_in_microns_range, microns_per_voxel, size_of_image );

network_histogram_plotter( network_statistics )
network_histogram_plotter( network_statistics, 'strand_z_direction' )

[ point_coordinates, arc_connectivity, arc_diameters ] = strand2casx( vertex_space_subscripts, strands2vertices, strand_subscripts, microns_per_voxel, lumen_radius_in_microns_range );

[ vertex_space_subscripts_2, vertex_scale_subscripts_2, bifurcation_vertices_2, strands2vertices_2, strand_space_subscripts_2, strand_scale_subscripts_2, microns_per_voxel_2, lumen_radius_in_microns_range_2, size_of_image_2 ] = casx2strand( point_coordinates, arc_connectivity, arc_diameters );

strand_subscripts_2 = cellfun( @( x, y ) [ x, y ], strand_space_subscripts_2, strand_scale_subscripts_2, 'UniformOutput', false );

% % placeholder vessel_directions input to not throw error
% % vessel_directions_placeholder = strand_space_subscripts ; % vessel_directions_placeholder = cellfun( @( x ) x( 1 : end - 1, : ), vessel_directions_placeholder, 'UniformOutput', false );
% [ vessel_directions_2 ] = get_vessel_directions_V3( strand_space_subscripts, microns_per_voxel );

network_statistics_2 = calculate_network_statistics( strand_subscripts_2, bifurcation_vertices_2, lumen_radius_in_microns_range_2, microns_per_voxel_2, size_of_image_2 );

network_histogram_plotter( network_statistics_2 )