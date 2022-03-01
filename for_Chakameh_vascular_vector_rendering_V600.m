%% for_Chakameh_vascular_vector_rendering_V600
%
% V600: this script calls the flow_field_subroutine function with two different inputs to report
% back two different images ( one for mc3dp with just two thresholds and another for DLS-MC with
% many more thresholds) SAM 3/13/19

path_to_flow_field_export_name = 'network_190313-010401_Fused_T123_T456_Raw_pad30' ;

path_to_flow_field_export = [ 'flow_field_export_', path_to_flow_field_export_name, '.mat' ];

load( path_to_flow_field_export )

time_stamp = char( datetime('now', 'TimeZone', 'local', 'Format', 'yyMMdd-HHmmss' )); 

path_to_flow_field_rendering_input = [ 'flow_field_input_', path_to_flow_field_export_name, time_stamp ];

%% amount of smoothing to be performed optional input
% smoothing_kernel_sigma_to_lumen_radius_ratio = 1 ;

%% mc3dp: two thresholds

file_name = 'mc3dp_' ;

mc3dp_path_to_flow_field_rendering_input = [ file_name, path_to_flow_field_rendering_input ];

% two cutoffs: capillaries, mid_size, large surface vessels
tissue_type_cutoffs_in_microns = [ 5.5, 15 ];            

save_inputs( mc3dp_path_to_flow_field_rendering_input )

flow_field_subroutine( mc3dp_path_to_flow_field_rendering_input )

%% dls_mc: thresholds at every recorded scale

file_name = 'dls_mc_' ;

% a cutoff at every detected scale
tissue_type_cutoffs_in_microns = lumen_radius_in_pixels_range( :, 1 )' * microns_per_pixel_xy ;

dls_mc_path_to_flow_field_rendering_input = [ file_name, path_to_flow_field_rendering_input ];

save_inputs( dls_mc_path_to_flow_field_rendering_input )

flow_field_subroutine( dls_mc_path_to_flow_field_rendering_input )

%% save_inputs
function save_inputs( path )

    save( path,                                                ...
                              'vertex_indices_in_strands',     ...
                                'edge_indices_in_strands',     ...
                              'edge_backwards_in_strands',     ...
                              'bifurcation_vertices',          ...
                              'edge_subscripts',               ...
                              'mean_edge_energies',            ...
                              'lumen_radius_in_pixels_range',  ...
    'smoothing_kernel_sigma_to_lumen_radius_ratio',            ...                                  
                              'size_of_image',                 ...
                              'microns_per_pixel_xy',          ...
                              'z_per_xy_length_of_pxl_ratio',  ...
                              'vertex_subscripts',             ...
                              'file_name',                     ...
                              'tissue_type_cutoffs_in_microns' )

end % save_inputs