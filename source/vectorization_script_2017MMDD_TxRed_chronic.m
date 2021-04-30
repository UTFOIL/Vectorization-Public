%% 2017MMDD_TX_RD chronic vectorization
% SAM 8/2/2019
% summarizes vectorizations done earlier this week, and sets a standard for future vectorizations.
% This script was used for the following datasets:

% % 20170802_TxRed_Chronic
% imaging_session_and_ROI_name = '20170802_TxRed_Chronic\Processed_Images_Stack03\' ;
% imaging_session_and_ROI_name = '20170802_TxRed_Chronic\Processed_Images_Stack05\' ;
% imaging_session_and_ROI_name = '20170802_TxRed_Chronic\Processed_Images_Stack06\' ;

% 20170809_TxRed_Chronic
% imaging_session_and_ROI_name = '20170809_TxRed_Chronic\Processed_Images_Stack01\' ;
% imaging_session_and_ROI_name = '20170809_TxRed_Chronic\Processed_Images_Stack02\' ;
%
% imaging_session_and_ROI_nam = '20170809_TxRed_Chronic\Processed_Images_Stack0' ;

% imaging_session_and_ROI_nam = '20170901_TxRed_Chronic\Processed_Images_Stack0' ;
% is_tiled = false

% 10/30/19 SAM
imaging_session_and_ROI_nam = '20170809_TxRed_Chronic\tiling 180628\' ;
input_name = '20170809_txRed_chronic_tiled.tif' ;
is_tiled = true ;

if is_tiled

    OutputDirectory = [ 'E:\2P imaging\', imaging_session_and_ROI_nam ];

    stackID_range = 1 ;        

else
    
    % stackID_range = 3 : 6
    % stackID_range = 1 : 6
    
end

% start_workflow = 'energy' ;
% start_workflow = 'vertices' ;
% start_workflow = 'edges' ;
% start_workflow = 'network' ;
start_workflow = 'none' ;

PSF_fudge_factor = 2 ;

for stackID = stackID_range
    
    if ~ is_tiled
        
        imaging_session_and_ROI_name = [ imaging_session_and_ROI_nam, num2str( stackID ), '\' ];

        OutputDirectory = [ 'E:\2P imaging\', imaging_session_and_ROI_name, 'PMT01_Red_Images\' ];

    end
    
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
                                        'gaussian_to_ideal_ratio',                  0.5,         ...
                                        'spherical_to_annular_ratio',               6/7, ...
                                        'microns_per_voxel',                        [ 1.07, 1.07, 5 ],          ...
                                        'radius_of_smallest_vessel_in_microns',     1.5,                          ...
                                        'radius_of_largest_vessel_in_microns',      60,                         ...
                                        'approximating_PSF',                        true,                       ...
                                        'excitation_wavelength',                    1.3 * PSF_fudge_factor,     ...
                                        'scales_per_octave',                        3,                        ... 
                                        'max_voxels_per_node_energy',               1e6,                    ...
                                        'vessel_wall_thickness_in_microns',         0                           };
        
        if is_tiled

            time_stamp = vectorize_V200([ OutputDirectory, input_name ], name_value_pair_inputs{ 1, : }); 

        else

            time_stamp = vectorize_V200([ OutputDirectory, 'PMT01_Red_Raw_Data_16bit.tif' ], name_value_pair_inputs{ 1, : }); 

        end
        
        otherwise

            name_value_pair_inputs = {  'OutputDirectory',                   OutputDirectory,  ...
                                        'PreviousBatch',                    'recent',   ...
                                        'PreviousWorkflow',                 'recent',          ...
                                        'StartWorkflow',                    start_workflow,    ...
                                   ...     'FinalWorkflow',                    'network',         ...
    ...                                    'Visual',                           'productive',      ...
                                        'Visual',                           'network',      ...
                                        'NewBatch',                         'no',              ...
                                        'Presumptive',                       true,             ...
                                 ...       'VertexCuration',                   'manual',          ...
                                        'VertexCuration',                   'none',          ...
                                        'length_dilation_ratio_vertices',   2,               ... new default
                                        'number_of_edges_per_vertex',       2,               ... new default 11/6/19
                            ...            'number_of_edges_per_vertex',       4,               ...
                                        'sigma_edge_smoothing',             0.5,               ...
                            ...            'EdgeCuration',                     'manual'           };
                                        'EdgeCuration',                     'none', ...
    ...                                    'SpecialOutput', { 'depth', 'strands', 'directions', 'upsampled', '3D-strands' }};
    ...                                    'SpecialOutput', { '3D-strands' }};
                                        'SpecialOutput', { 'depth-stats' }};
    ...                                    'SpecialOutput', { 'upsampled' }};    
    ...                                     'SpecialOutput', { 'depth', 'directions' }};
...                                        'SpecialOutput', 'histograms' };    

            time_stamp = vectorize_V200( name_value_pair_inputs{ 1, : }); 
    end            
end
