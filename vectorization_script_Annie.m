%% Annie vectorization
% 2/27/20 SAM

% Input stacks:

% OutputDirectory = 'E:\Annie\ScanImage\' ;
OutputDirectory = 'E:\Annie\ScanImage\200302 vascular compare distortion rg vs gg\' ; 

microns_per_voxel = [ ];

% input_names{ 1 } = '200207_RG_14p0MHz_1e4_00001_2p0x_16frAvg_00001.tif' ;
input_names{ 1 } = 'resStack.tif'   ; microns_per_voxel( 1, 1 : 3 ) = [ 1.3, 0.9, 5 ];
% input_names{ 2 } = 'galvoStack.tif' ; microns_per_voxel( 2, 1 : 3 ) = [ 1.1, 1.1, 5 ];

number_of_stacks = length( input_names );

stackID_range = 1 : number_of_stacks ;        

% Workflow:

% start_workflow = 'energy' ;
% start_workflow = 'vertices' ;
% start_workflow = 'edges' ;
start_workflow = 'network' ;
% start_workflow = 'none' ;

% Settings:

PSF_fudge_factor = 2 ;

for stackID = stackID_range
    
    switch start_workflow
    
        case 'energy'
    
            name_value_pair_inputs = {  'OutputDirectory',                   OutputDirectory, ...
                                        'PreviousBatch',                    'none',               ...
                                        'PreviousWorkflow',                 'none',            ...
                                        'StartWorkflow',                    'energy',                           ...
                                        'FinalWorkflow',                    'network',                              ...
                            ...            'FinalWorkflow',                    'energy',                              ...                    
                                        'Visual',                           'productive',                       ...
                                        'NewBatch',                         'yes',                    ...
                                        'Presumptive',                       true,                              ...
                                        'matching_kernel_string',            '3D gaussian conv annular pulse',            ... % !!!! depricated parameter !!!! built-in parmater to energy_filter, A = 6 !!!!! make this a recorded input
                 ...                       'symmetry_ratio_factor',                    1.5,         ... !!! depricated parameter.  
                                        'gaussian_to_ideal_ratio',                  0.8,         ...
                                        'spherical_to_annular_ratio',               0.8, ...
...                                        'microns_per_voxel',                        [0.78, 0.78, 5],          ... 
                                        'microns_per_voxel',                        microns_per_voxel( stackID, : ),          ... 
                                        'radius_of_smallest_vessel_in_microns',     1.5,                          ...
                                        'radius_of_largest_vessel_in_microns',      60,                         ...
                                        'approximating_PSF',                        true,                       ...
                                        'excitation_wavelength',                    1.3 * PSF_fudge_factor,     ...
                                        'scales_per_octave',                        3,                        ... 
                                        'max_voxels_per_node_energy',               1e6,                    ...
                                        'vessel_wall_thickness_in_microns',         0                           };
        

            time_stamp = vectorize_V200([ OutputDirectory, input_names{ stackID }], name_value_pair_inputs{ 1, : }); 

        
        otherwise

            name_value_pair_inputs = {  'OutputDirectory',                   OutputDirectory,  ...
                                        'PreviousBatch',                    'recent',   ...
                                        'PreviousWorkflow',                 'recent',          ...
                                        'StartWorkflow',                    start_workflow,    ...
                                   ...     'FinalWorkflow',                    'network',         ...
    ...                                    'Visual',                           'productive',      ...
                     ...                   'Visual',                           'vertices',      ...
                                        'NewBatch',                         'no',              ...
                                        'Presumptive',                       true,             ...
                           ...             'VertexCuration',                   'manual',          ...
                            ...            'VertexCuration',                   'none',          ...
                            ...            'VertexCuration',                   'mutual edges',          ...
                           ...          'VertexCuration',                   'auto',          ...
                                        'space_strel_apothem_edges',        int64(1),           ...
                                        'length_dilation_ratio_vertices',   2,               ... new default
                                        'number_of_edges_per_vertex',       3,               ... new default 11/6/19
                            ...            'number_of_edges_per_vertex',       4,               ...
                                        'sigma_edge_smoothing',             0.5,               ...
                          ...              'EdgeCuration',                     'manual', ...
                           ...             'EdgeCuration',                     'none', ...
    ...                                    'SpecialOutput', { 'depth', 'strands', 'directions', 'upsampled', '3D-strands' }};
    ...                                    'SpecialOutput', { '3D-strands' }};
                                        'SpecialOutput', { 'depth-stats' }};
    ...                                    'SpecialOutput', { 'upsampled' }};    
    ...                                     'SpecialOutput', { 'depth', 'directions' }};
...                                        'SpecialOutput', 'histograms' };    

            time_stamp = vectorize_V200( name_value_pair_inputs{ 1, : }); 
    end            
end
