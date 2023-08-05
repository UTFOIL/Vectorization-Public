%% register chronic vectors script

vector_set_A_workflow_mat_file = 'E:\2P imaging\2021_Chronic_Imaging\Eddy week 1\batch_210604-203819\settings\workflow_220818-110430' ;
vector_set_B_workflow_mat_file = 'E:\2P imaging\2021_Chronic_Imaging\Eddy week 3\batch_210616-194614\settings\workflow_220818-103236' ;

register_strands_wrapper( vector_set_A_workflow_mat_file, ...
                          vector_set_B_workflow_mat_file )

%% functions
function register_strands_wrapper( vector_set_A_workflow_mat_file, ...
                                   vector_set_B_workflow_mat_file )
    % SAM 9/21/22
    
    %% load strand vector sets A and B
    
    [ strand_subscripts_A, microns_per_voxel_A, lumen_radius_in_microns_range_A, ...
       production_times_A,   vector_set_batch_folder_A,   vector_set_ROI_name_A  ] = load_strands( vector_set_A_workflow_mat_file );
    
    [ strand_subscripts_B, microns_per_voxel_B, lumen_radius_in_microns_range_B, ...
       production_times_B,   vector_set_batch_folder_B,   vector_set_ROI_name_B  ] = load_strands( vector_set_B_workflow_mat_file );
    
    %             strand_subscripts_A =             strand_subscripts ;
    %             microns_per_voxel_A =             microns_per_voxel ;
    % lumen_radius_in_microns_range_A = lumen_radius_in_microns_range ;
    
    %% register the strand vectors
    [ strand_subscripts_B_registered, strand_registration_scores_B, ...
      B_strands_to_A_strands, transformation_matrix                 ] = register_strands( strand_subscripts_A, microns_per_voxel_A, lumen_radius_in_microns_range_A, ...
                                                                                          strand_subscripts_B, microns_per_voxel_B, lumen_radius_in_microns_range_B  );
                                                                     
    %% save registered vectors in folder A
    save_strands( production_times_B,   vector_set_batch_folder_B,   vector_set_ROI_name_B,                                      ...
                  strand_subscripts_B_registered, B_strands_to_A_strands, transformation_matrix, vector_set_A_workflow_mat_file, ...
                                                                                                 vector_set_B_workflow_mat_file  )
    
    %% identify largest mismatches
    
    %% manually inspect the largest mismatches and decide ground truth
    [ strand_subscripts_A, strand_subscripts_B ] = curate_strand_registration( strand_subscripts_A, original_image_A, ...
                                                                               strand_subscripts_B, original_image_B, strand_registration_scores_B );
                                                  
    %% save curated vector sets in respective folders
            

end

function [ strand_subscripts, microns_per_voxel, lumen_radius_in_microns_range, ...
           production_times,  vector_set_batch_folder,     vector_set_ROI_name, registration_path, curation_path  ] = load_strands( vector_set_workflow_mat_file )

    filesep_indcs = regexp( vector_set_workflow_mat_file, filesep );

    vector_set_batch_folder = vector_set_workflow_mat_file( 1 : filesep_indcs( end - 2 )                        );
%         vector_set_ROI_name     = vector_set_workflow_mat_file(     filesep_indcs( end     ) + 8 + 1 + 13 : end - 4 ); % '...\workflow_YYMMDD-hhmmss' : end  - '.mat'

    load( vector_set_workflow_mat_file, 'production_times', 'ROI_names' )

    vector_set_network_________file = [ vector_set_batch_folder, filesep,  'vectors', filesep, 'network_' , production_times( 5 ), ROI_names{ 1 }, '.mat' ]; % !!!! only supports one ROI in the ROI_names variable currently ? how to interpret multiple ROI's?
    vector_set_energy_settings_file = [ vector_set_batch_folder, filesep, 'settings', filesep,  'energy_' , production_times( 2 ), ROI_names{ 1 }, '.mat' ]; % !!!! only supports one ROI in the ROI_names variable currently ? how to interpret multiple ROI's?

    load( vector_set_network_________file, 'strand_subscripts'                                  )
    load( vector_set_energy_settings_file, 'microns_per_voxel', 'lumen_radius_in_microns_range' )
    
    % ??? make registratino dir ????
    registration_path = [ vector_set_batch_folder, '' ];

end

function save_strands( production_times,   vector_set_batch_folder,   vector_set_ROI_name, varargin )

%         save( )
    
end