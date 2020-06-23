%% import data from LPPD.
% SAM 7/23/19

%% load casX file

casX_file_path = 'C:\Users\sam7343\Box Sync\AA\Andreas Linninger\Grant Hartung\SD1.203.tortuous.casx' ; % 7/23/19

[ point_coordinates, arc_connectivity, arc_diameters ] = casx_file2mat( casX_file_path );

% mirror in z
point_coordinates( :, 3 ) = - point_coordinates( :, 3 );

[ vertex_space_subscripts, vertex_scale_subscripts, bifurcation_vertices, strands2vertices, strand_space_subscripts, strand_scale_subscripts, microns_per_voxel, lumen_radius_in_microns_range, size_of_image ] = casx2strand( point_coordinates, arc_connectivity, arc_diameters );

%% view network histograms

strand_subscripts = cellfun( @( x, y ) [ x, y ], strand_space_subscripts, strand_scale_subscripts, 'UniformOutput', false );

% % placeholder vessel_directions input to not throw error
% % vessel_directions_placeholder = strand_space_subscripts ; % vessel_directions_placeholder = cellfun( @( x ) x( 1 : end - 1, : ), vessel_directions_placeholder, 'UniformOutput', false );
% [ vessel_directions ] = get_vessel_directions_V3( strand_space_subscripts, microns_per_voxel );

network_statistics = calculate_network_statistics( strand_subscripts, bifurcation_vertices, lumen_radius_in_microns_range, microns_per_voxel, size_of_image );

network_histogram_plotter( network_statistics )
network_histogram_plotter( network_statistics, 'strand_z_direction' )

%% replace large surface vessels with originals

% this section "replace large surface vessels with originals" only works with one tissue type cutoff
tissue_type_cutoffs_in_microns = [ 5.5 ];

microns_per_pixel_xy = microns_per_voxel( 1 );

% lumen_radius_in_pixels_range = lumen_radius_in_microns_range ./ microns_per_voxel ;
%     
% tissue_type_cutoffs        = find(   tissue_type_cutoffs_in_microns        ...
%                                    / microns_per_pixel_xy                  ...
%                                    > lumen_radius_in_pixels_range( :, 1 ), ...
%                                    1, 'last'                               );

tissue_type_cutoffs        = find(   tissue_type_cutoffs_in_microns ...
                                   > lumen_radius_in_microns_range, ...
                                   1, 'last'                               );

is_vessel_large = cellfun( @( x ) tissue_type_cutoffs <= median( x ), strand_scale_subscripts );

strand_space_subscripts( is_vessel_large ) = [ ];
strand_scale_subscripts( is_vessel_large ) = [ ];

save( 'temp_workspace' )

clear

load('E:\2P imaging\20170802_TxRed_Chronic\Processed_Images_Stack04\PMT01_Red_Images\synthetic data for MC sim\flow_field_export_empirical.mat')

% tissue_type_cutoffs_in_microns = [ 5.5 ];
% 
% microns_per_pixel_xy = microns_per_voxel( 1 );
% 
% lumen_radius_in_pixels_range = lumen_radius_in_microns_range ./ microns_per_voxel ;
    
tissue_type_cutoffs        = find(   tissue_type_cutoffs_in_microns ...
                                   > lumen_radius_in_microns_range, ...
                                   1, 'last'                               );

is_vessel_large = cellfun( @( x ) tissue_type_cutoffs <= median( x ), strand_scale_subscripts );

strand_space_subscripts_large_original = strand_space_subscripts( is_vessel_large );
strand_scale_subscripts_large_original = strand_scale_subscripts( is_vessel_large );

strand_space_coordinates_large_original = cellfun( @( x ) ( x - 1 ) .* microns_per_voxel,              strand_space_subscripts_large_original, 'UniformOutput', false );
strand_radii_large_original             = cellfun( @( x ) interp1( lumen_radius_in_microns_range, x ), strand_scale_subscripts_large_original, 'UniformOutput', false );

load( 'temp_workspace' )

% extend the size LUT if needed
common_ratio = lumen_radius_in_microns_range( 2 ) / lumen_radius_in_microns_range( 1 );

while lumen_radius_in_microns_range( end ) <= max( cellfun( @max, strand_radii_large_original ))

    lumen_radius_in_microns_range( end + 1 ) = lumen_radius_in_microns_range( end ) * common_ratio ; 
    
end

lumen_radius_in_pixels_range = lumen_radius_in_microns_range ./ microns_per_voxel ;

scale_subscript_range = ( 1 : length( lumen_radius_in_microns_range ))';

strand_space_subscripts_large_original = cellfun( @( x ) x ./ microns_per_voxel + 1,                                         strand_space_coordinates_large_original, 'UniformOutput', false );
strand_scale_subscripts_large_original = cellfun( @( x ) interp1( lumen_radius_in_microns_range, scale_subscript_range, x ), strand_radii_large_original, 'UniformOutput', false );

strand_space_subscripts = [ strand_space_subscripts; strand_space_subscripts_large_original ];
strand_scale_subscripts = [ strand_scale_subscripts; strand_scale_subscripts_large_original ];

% extend the image if needed (microns and voxels are interchangeable due to our choice of microns_per_voxel = [ 1, 1, 1 ]
size_of_image = ceil( max( cell2mat( cellfun( @( x, r ) max( x + lumen_radius_in_microns_range( ceil( r ))), strand_space_subscripts, strand_scale_subscripts, 'UniformOutput', false ))));

%% view network histograms

strand_subscripts = cellfun( @( x, y ) [ x, y ], strand_space_subscripts, strand_scale_subscripts, 'UniformOutput', false );

% % placeholder vessel_directions input to not throw error
% % vessel_directions_placeholder = strand_space_subscripts ; % vessel_directions_placeholder = cellfun( @( x ) x( 1 : end - 1, : ), vessel_directions_placeholder, 'UniformOutput', false );
% [ vessel_directions ] = get_vessel_directions_V3( strand_space_subscripts, microns_per_voxel );

network_statistics = calculate_network_statistics( strand_subscripts, bifurcation_vertices, lumen_radius_in_microns_range, microns_per_voxel, size_of_image );

network_histogram_plotter( network_statistics )
network_histogram_plotter( network_statistics, 'strand_z_direction' )

%% render flow field

path_to_flow_field_export = 'flow_field_export_synthetic';

strand_energies = cellfun( @(x) -Inf * x, strand_scale_subscripts, 'UniformOutput', false );

vertex_subscripts = [ double( vertex_space_subscripts ), double( vertex_scale_subscripts )];

% microns_per_pixel_xy = microns_per_voxel( 1 );
z_per_xy_length_of_pxl_ratio = microns_per_voxel( 3 ) / microns_per_voxel( 1 );

file_name = 'synthetic' ;

% tissue_type_cutoffs_in_microns = [ 5.5 ];

save( path_to_flow_field_export,                       ...
        'bifurcation_vertices',          ...
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

flow_field_subroutine( path_to_flow_field_export );