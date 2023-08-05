%% Shaun vectorization
% 4/8/2021 SAM

%% Input stacks:
% OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Control_C57\CC22703\A\210311_Session_1\' ; input_names{ 1 } = 'A_210311_tiled_medFilt_bgAdj.tif' ; 
% OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Control_C57\CC22703\A\210326_Session_2\' ; input_names{ 1 } = 'A_210326_tiled_medFilt_bgAdj.tif' ; 
% OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Control_C57\CC22703\C\210311_Session_1\' ; input_names{ 1 } = 'C_210311_tiled_medFilt_bgAdj.tif' ; 
% OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Control_C57\CC22703\C\210325_Session_2\' ; input_names{ 1 } = 'C_210325_tiled_medFilt_bgAdj.tif' ; 

% original_image = tif2mat( 'E:\2P imaging\2021_Chronic_Imaging\Control_C57\CC22703\C\210325_Session_2\C_210325_tiled_medFilt_bgAdj.tif' );

% mat2tif( 2^15*log( double(original_image )+ 1 )/log( double(max(original_image(:)) )+ 1 ), 'E:\2P imaging\2021_Chronic_Imaging\Control_C57\CC22703\C\210325_Session_2\C_210325_tiled_medFilt_bgAdj_log.tif' );
% OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Control_C57\CC22703\C\210325_Session_2\' ; input_names{ 1 } = 'C_210325_tiled_medFilt_bgAdj_log.tif' ; 

% mat2tif( 2^15*log( double(original_image )+ median(double(original_image(:))) )/log( double(max(original_image(:)) )+median(double(original_image(:))) ), 'E:\2P imaging\2021_Chronic_Imaging\Control_C57\CC22703\C\210325_Session_2\C_210325_tiled_medFilt_bgAdj_log_stable.tif' );
% OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Control_C57\CC22703\C\210325_Session_2\' ; input_names{ 1 } = 'C_210325_tiled_medFilt_bgAdj_log_stable.tif' ; % rough edge curation performed SAM 5/17/21

%% input stacks 210604
% original_name = 'Eddy_Fused_Raw_w1' ;
% 
% OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Eddy week 1\' ;
% 
% input_names{ 1 } = [ original_name, '_log_stable.tif' ]; 
% 
% original_image = tif2mat([ OutputDirectory, original_name, '.tif' ]);
% 
% original_image = original_image - min( original_image( : ));
% 
% transformed_image_path = [ OutputDirectory, input_names{ 1 }];
% 
% if ~ exist( transformed_image_path,'file')
% 
%     mat2tif( 2^15*log( double(original_image )+ median(double(original_image(:))))/log( double(max(original_image(:)) )+median(double(original_image(:))) ), transformed_image_path );
% 
% end
% % moved this log-transform pre-processing to inside the vectorization function to occure whenever a mask is applied - SAM 6/16/21

input_names = { };

% input_names{ 1 } = 'Eddy_Fused_Raw_w3.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Eddy week 3\' ; % SAM 6/16/21 %% subvolume for closeup SAM 10/31/22
% input_names{ 1 } = 'Eddy_Fused_Raw_w5.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Eddy week 5\' ; % SAM 6/18/21 %% subvolume for closeup SAM 10/31/22
% input_names{ 1 } = 'Doug_Fused_Raw_w1.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Doug week 1\' ; % SAM 6/26/21  %% SAM 3/17/22 % removed the log transform
% % input_names{ 1 } = 'Dan_w1_01.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Dan week 1\' ; % SAM 8/6/21 SAM
% % input_names{ 2 } = 'Dan_w1_02.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Dan week 1\' ; % SAM 8/10/21 SAM
% % for idx = 3 : 6, input_names{ idx - 2 } = [ 'Dan_w1_0', num2str(idx), '.tif' ]; end, OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Dan week 1\' ; % SAM 8/10/21 SAM
% input_names{ 1 } = 'Dan_Fused_Raw_w1.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Dan week 1\fused\' ; % SAM 9/2/21 % intermouse comparison 8/18/2022
% input_names{ 1 } = 'Elmer_Fused_Raw_w1.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Elmer week 1\' ; % SAM 5/10/22 % removed the log transform % intermouse comparison 8/18/2022
% % input_names{ 1 } = 'Doug_Fused_Raw_w3.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Doug week 3\' ; % SAM 5/10/22 % removed the log transform 
% input_names{ 1 } = 'Doug_Fused_Raw_w5.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Doug week 5\' ; % SAM 5/18/22 % removed the log transform 
% input_names{ 1 } = 'Elmer_Fused_Raw_w3.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Elmer week 3\' ; % SAM 5/23/22 % removed the log transform 
% input_names{ 1 } = 'Eddy_Fused_Raw_w1.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Eddy week 1\' ; % SAM 5.31.22 for plot generation % intermouse comparison 8/18/2022 %% subvolume for closeup SAM 10/31/22
% input_names{ 1 } = 'Brett_Fused_Raw_w5.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 5\' ; % SAM 6/14/22 % put the log transform back in
input_names{ 1 } = 'Brett_Fused_Raw_w1.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 1\truncated\' ; % SAM 6/30/22 % used subset of dataset that corresponds to wk 5
% input_names{ 1 } = 'Brett_Fused_Raw_w4.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 4\' ; 
% input_names{ 1 } = 'Brett_Fused_Raw_w3.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 3\truncated\' ; 
% input_names{ 1 } = 'Brett_Fused_Raw_w2.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 2\New folder\Truncated\' ; 
% input_names{ 1 } = 'Ant_Fused_Raw_w5.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Ant week 5\New folder\' ; 
% input_names{ 1 } = 'Ant_Fused_Raw_w4.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Ant week 4\New folder\' ; 

% input_names{ 1 } = 'Dwight_Fused_Raw_w5.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Dwight week 5\' ;
% input_names{ 1 } = 'Dwight_Fused_Raw_w1.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Dwight week 1\' ;

% input_names{ 1 } = 'Ant_Fused_Raw_w1.tif' ; OutputDirectory = 'E:\2P imaging\2021_Chronic_Imaging\Ant week 1\' ; 

microns_per_voxel = [ ];
% microns_per_voxel( 1, 1 : 3 ) = [1.34, 1.36, 3]; % [ y, x, z ] == microns_per_voxel - SAM 3/19/21
microns_per_voxel( 1, 1 : 3 ) = [1.37, 1.37, 3]; % [ y, x, z ] == microns_per_voxel - SAM 6/16/21

number_of_stacks = length( input_names );

stackID_range = 1 : number_of_stacks ;        

%% Workflow:

% start_workflow = 'energy' ;
% start_workflow = 'vertices' ;
% start_workflow = 'edges' ;
% start_workflow = 'network' ;
start_workflow = 'none' ;

% Settings:

% PSF_fudge_factor = 2 ;
PSF_fudge_factor = 1 ;

% % volume octaves (doubling of spherical volume)
% scales_per_octave = 3 ;
% scales_per_octave = 1 / 6 ;
% scales_per_octave = 1 / 3 ; % 1 per radius octave
% scales_per_octave = 2 / 3 ; % 2 per radius octave
scales_per_octave = 1 ; % 3 per radius octave

for stackID = stackID_range
    
    switch start_workflow
    
        case 'energy'
    
            name_value_pair_inputs = {  'OutputDirectory',                   OutputDirectory, ...
                                        'PreviousBatch',                    'none',               ...
                                        'PreviousWorkflow',                 'none',            ...
                                        'StartWorkflow',                    'energy',                           ...
                                       'FinalWorkflow',                    'network',                              ...
                                   ...     'FinalWorkflow',                    'energy',                              ...                    
                                        'Visual',                           'productive',                       ...
                                        'SpecialOutput', { 'histograms', 'vmv', 'depth', 'strands', 'directions' },        ...                                        
                                        'NewBatch',                         'yes',                    ...
                                        'Presumptive',                       true,                              ...
                                        'matching_kernel_string',            '3D gaussian conv annular pulse',            ... % !!!! depricated parameter !!!! built-in parmater to energy_filter, A = 6 !!!!! make this a recorded input
                 ...                       'symmetry_ratio_factor',                    1.5,         ... !!! depricated parameter.  
...                                        'gaussian_to_ideal_ratio',               0.6,         ...
...                                        'gaussian_to_ideal_ratio',               2^-(scales_per_octave*3/2),         ... % half the radius is due to blurring when the radius is sampled at every doubling
                                        'gaussian_to_ideal_ratio',               ( 1 - 2 ^ - ( 2 / scales_per_octave / 3 )) .^ 0.5 ,         ... % pythagorean half (root(2)/2) of the radius is due to blurring when the radius is sampled at every doubling
                                     'spherical_to_annular_ratio',               0.8, ...
    ...                                    'gaussian_to_ideal_ratio',               0.5,         ...
    ...                                 'spherical_to_annular_ratio',               0.5, ...
...                                        'microns_per_voxel',                        [0.78, 0.78, 5],          ... 
                                        'microns_per_voxel',                        microns_per_voxel( stackID, : ),          ... 
                                       'radius_of_smallest_vessel_in_microns',      1.5,                         ...
                                        'radius_of_largest_vessel_in_microns',      60,                         ...
                                        'approximating_PSF',                        true,                       ...
                                        'excitation_wavelength',                    1.3 * PSF_fudge_factor,     ...
                                        'scales_per_octave',                        scales_per_octave,                        ... 
                                        'max_voxels_per_node_energy',               1e6,                    ...
                                        'vessel_wall_thickness_in_microns',         0                       , ...
                                        ...
                                        'VertexCuration',                   'machine-manual' };

            time_stamp = vectorize_V200([ OutputDirectory, input_names{ stackID }], name_value_pair_inputs{ 1, : }); 

        
        otherwise

            name_value_pair_inputs = {  'OutputDirectory',                   OutputDirectory,  ...
                                        'PreviousBatch',                    'prompt',   ...
                      ...                  'PreviousBatch',                    'recent',   ...
                     ...                   'PreviousBatch',                    '210810-022653', ...
                    ...                    'PreviousBatch',                    '220126-010215', ...
                                        'PreviousWorkflow',                 'recent',          ...
                                        'StartWorkflow',                    start_workflow,    ...
                                   ...     'FinalWorkflow',                    'network',         ...
    ...                                    'Visual',                           'productive',      ...
                     ...                   'Visual',                           'vertices',      ...
                                        'NewBatch',                         'no',              ...
                                        'Presumptive',                       true,             ...
                            ...            'VertexCuration',                   'manual',          ...
                            ...            'VertexCuration',                   'machine-manual',          ...
                            ...            'VertexCuration',                   'none',          ...
                           ...             'VertexCuration',                   'mutual edges',          ...
                           ...          'VertexCuration',                   'auto',          ...
                                        'space_strel_apothem_edges',        int64(1),           ...
                                        'number_of_edges_per_vertex',       2,               ... 
                            ...            'number_of_edges_per_vertex',       4,               ...
                                        'sigma_edge_smoothing',             0.5,               ...
                              ...          'EdgeCuration',                     'manual', ...
                           ...             'EdgeCuration',                     'none', ...
    ...                                    'SpecialOutput', { 'depth', 'strands', 'directions', 'upsampled', '3D-strands' }};
    ...                                    'SpecialOutput', { '3D-strands' }};
...                                        'SpecialOutput', { 'depth-stats' }};
              ...                           'is_combining_strands', false,'SpecialOutput', { 'upsampled' }};    
   ...                                      'SpecialOutput', { 'depth', 'directions' }};
        ...                                 'SpecialOutput', { 'depth-stats' }};
  ...                                     'is_combining_strands', false,'SpecialOutput', { 'histograms', 'vmv', 'depth', 'strands', 'directions' }};    
    ...                                   'is_combining_strands', false,'SpecialOutput', { 'histograms', 'vmv', 'depth', 'directions', '3D-directions' }};   
                                       'is_combining_strands', false,'SpecialOutput', { 'histograms', 'vmv', '3D-directions' }};   
                 ...                     'is_combining_strands', false,'SpecialOutput', { 'histograms', 'vmv', '3D-directions', 'directions' }};   
...                                      'is_combining_strands', true,'SpecialOutput', { 'histograms', 'vmv', '3D-radii' }};               
        ...                              'is_combining_strands', false,'SpecialOutput', { 'histograms', 'vmv', '3D-radii' }};               
            ...                                          'is_combining_strands', true, 'SpecialOutput', { 'histograms', 'flow-field' }};
   ...                                     'SpecialOutput', {'histograms'}};
      ...                                  'SpecialOutput', {'vmv'}};
           ...                               'SpecialOutput', { 'vmv', 'directions' }};
            time_stamp = vectorize_V200( name_value_pair_inputs{ 1, : }); 
    end            
end
