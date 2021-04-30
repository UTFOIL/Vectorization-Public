%% test random anatomy generation
% SAM 7/1/19

path_to_batch = 'E:\2P imaging\20170802_TxRed_Chronic\Processed_Images_Stack04\PMT01_Red_Images\batch_190531-164319\' ;

ROI_name = 'PMT01_Red_Raw_Data_16bit' ;

network_handle = 'network_190620-224054' ;

path_to_energy_settings = [ path_to_batch, 'settings\energy_190603-172534'           ];
path_to_edges           = [ path_to_batch, 'vectors\edges_190612-144547_',  ROI_name ];
path_to_network         = [ path_to_batch, 'vectors\', network_handle, '_', ROI_name ];
visual_vector_directory = [ path_to_batch, 'visual_vectors\'                         ];

strands_visual_spheres_file        = [ visual_vector_directory,  'simulated_', network_handle, '_strands',        '_spheres', ROI_name, '.tif' ];
strands_visual_centerline_file     = [ visual_vector_directory,  'simulated_', network_handle, '_strands',    '_centerlines', ROI_name, '.tif' ];

load( path_to_energy_settings );
load( path_to_edges           );
load( path_to_network         );

network_histogram_plotter( network_statistics )

strand_space_subscripts = cellfun( @( x ) x( :, 1 : 3 ), strand_subscripts, 'UniformOutput', false );
strand_scale_subscripts = cellfun( @( x ) x( :,   4   ), strand_subscripts, 'UniformOutput', false );

size_of_image = [ 512, 512, 121 ];

radii_threshold = 5.5 ; % microns

[ strand_space_subscripts_2, strand_scale_subscripts_2 ] = randomize_anatomy( strand_space_subscripts, strand_scale_subscripts, microns_per_voxel, lumen_radius_in_microns_range, size_of_image, network_statistics, radii_threshold );

strand_subscripts_2 = cellfun( @( x, y ) [ x, y ], strand_space_subscripts_2, strand_scale_subscripts_2, 'UniformOutput', false );

% % placeholder vessel_directions input to not throw error
% % vessel_directions_placeholder = strand_space_subscripts ; % vessel_directions_placeholder = cellfun( @( x ) x( 1 : end - 1, : ), vessel_directions_placeholder, 'UniformOutput', false );
% [ vessel_directions ] = get_vessel_directions_V3( strand_space_subscripts, microns_per_voxel );

network_statistics_2 = calculate_network_statistics( strand_subscripts_2, bifurcation_vertices, lumen_radius_in_microns_range, microns_per_voxel, size_of_image );

network_histogram_plotter( network_statistics_2 )

visualize_edges_V180( strand_subscripts_2, mean_strand_energies,    ...
                      lumen_radius_in_pixels_range,            ...
                      size_of_image, strands_visual_spheres_file, ...
                      strands_visual_centerline_file              )
                  
%% render flow field

vertex_subscripts = [ double( vertex_space_subscripts ), double( vertex_scale_subscripts )];
%             edge_subscripts   = cellfun( @( x, y ) [ double( x ), double( y )], edge_space_subscripts, edge_scale_subscripts, 'UniformOutput', false );
microns_per_pixel_xy = microns_per_voxel( 1 );
z_per_xy_length_of_pxl_ratio = microns_per_voxel( 3 ) / microns_per_voxel( 1 );
%             sigma_per_size = 1 ;
%             pixels_per_sigma_range = lumen_radius_in_pixels_range ;
%             smoothing_kernel_sigma_to_lumen_radius_ratio = sigma_strand_smoothing ;
file_name = 'randomized empirical' ;

path_to_flow_field_export = 'flow_field_export_randomized_empirical' ;

%             % a cutoff at every detected scale
%             tissue_type_cutoffs_in_microns = lumen_radius_in_pixels_range( :, 1 )' * microns_per_pixel_xy ;
tissue_type_cutoffs_in_microns = [ 5.5 ];
save( path_to_flow_field_export,                       ...
...            'vertex_indices_in_strands',     ...
...              'edge_indices_in_strands',     ...
...            'edge_backwards_in_strands',     ...
'bifurcation_vertices',          ...
...            'edge_subscripts',               ...
...            'edge_energies',                 ...
...            'mean_edge_energies',            ...
'strand_space_subscripts',  ...
'strand_scale_subscripts',  ...
'strand_energies',          ...
'vessel_directions',  ...
'lumen_radius_in_pixels_range',  ...
'lumen_radius_in_microns_range', ...
...   'smoothing_kernel_sigma_to_lumen_radius_ratio',            ...
'size_of_image',                 ...
'microns_per_pixel_xy',          ...
'microns_per_voxel',             ...
'z_per_xy_length_of_pxl_ratio',  ...
'vertex_subscripts',             ...
'file_name',                     ...
'tissue_type_cutoffs_in_microns' )

flow_field_subroutine( path_to_flow_field_export )
                  