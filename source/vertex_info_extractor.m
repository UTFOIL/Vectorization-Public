function vertex_info_extractor(vertex_set_list, extraction_type_list, root_directories_and_batch_names, target_directory, varargin)
% extracting_for_training determines if script will attempt to extract
% features of curated vertex data
overallStartTime = tic;

target_directory = char(target_directory);

extracting_for_training = logical(numel(intersect(vertex_set_list,{'training'}))); 
vertex_set_list(strcmp(vertex_set_list,'training')) = [];

extracting_derivatives = logical(numel(intersect(extraction_type_list,{'derivatives'})));
extracting_intensity_statistics = logical(numel(intersect(extraction_type_list,{'intensity'})));
extracting_spatial_and_energy_info = logical(numel(intersect(extraction_type_list,{'original'})));
if logical(numel(intersect(extraction_type_list,{'all'})))
    extracting_derivatives = true;
    extracting_intensity_statistics = true; % !!! no features defined later %SAM 220530
    extracting_spatial_and_energy_info = true; % !!! no features defined later
end

%% This scripts gets feature info of noncurated (and curated for training) vertices for machine learning curation.
% WAS 6/30/2019. 
%
% Gets and saves: all six curvatures, all three gradients

%% Get directory information. Ripped and adapted from vectorize_V200.m
for root_directory_and_batch_name_index = 1:numel(root_directories_and_batch_names)
    root_directory_and_batch_name = char(root_directories_and_batch_names{root_directory_and_batch_name_index});         

batch_directory = [ root_directory_and_batch_name, filesep ];
    
settings_directory = [ batch_directory, 'settings', filesep ];
    
batch_settings     = [ settings_directory, 'batch' ];   

% load( previous_batch_settings )
load( batch_settings )
  
    [                       ~, ...
               data_directory, ...
        visual_data_directory, ...
             vector_directory, ...
      visual_vector_directory, ...
           curation_directory, ...
           settings_directory  ]    = get_directories( batch_directory );

workflow_listing = dir([ root_directory_and_batch_name, 'settings', filesep, 'workflow_*.mat' ]); 
  
%Most recent timestamp is default
time_stamp = workflow_listing( end ).name( 10 : 22 );

%If specific time stamp
if ~isempty(varargin)
    time_stamp = varargin{2};
end
       
% time_stamp = root_directory_and_batch_name( end-13 : end );
       
workflow_handle = [ 'workflow_', time_stamp ];

path_to_workflow_settings = [ settings_directory, workflow_handle ];

load( path_to_workflow_settings, 'production_times' )  

original_data_handle = 'original' ;
 
           energy_handle = [  'energy_', production_times{ 2 }];
         vertices_handle = ['vertices_', production_times{ 3 }];
        
         path_to_energy_settings   = [ settings_directory,            energy_handle ];
         path_to_vertices_settings = [ settings_directory,          vertices_handle ];

if ~ isempty( varargin ), ROI_index_range = varargin{ 1 }; end;
         
for ROI_index = ROI_index_range

%     path_to_energy_data             = [     data_directory,               energy_handle, ROI_names{ ROI_index }]; % mat file path        
    path_to_original_data           = [     data_directory,        original_data_handle, ROI_names{ ROI_index }]; % h5  file path
    path_to_vertices                = [   vector_directory,             vertices_handle, ROI_names{ ROI_index }]; %  vectors path
    path_to_curated_vertices        = [   vector_directory, 'curated_', vertices_handle, ROI_names{ ROI_index }]; %  vectors path    

    path_to_saved_curation          = [ curation_directory,             vertices_handle, ROI_names{ ROI_index }]; % logicals path

    %new ones 
    path_to_vertex_features          = [ target_directory, '\uncurated\', vertices_handle, '_featureSet_',          ROI_names{ ROI_index }];
    path_to_curated_vertex_features  = [ target_directory, '\curated\',   vertices_handle, '_featureSet_CURATED_',  ROI_names{ ROI_index }];
    path_to_training_vertex_features = [ target_directory, '\training\',  vertices_handle, '_featureSet_TRAINING_', ROI_names{ ROI_index }];
     
    load(path_to_energy_settings)
     
    if ~ isempty( varargin )
        path_to_vertex_features = target_directory;
    end

%% get vertex derivative info. workflow code: 'derivatives'
if extracting_derivatives

    load(path_to_vertices_settings)
    
    disp([newline 'Extracting Derivative Information'])
    start_time = tic;
    
    apothem_per_radius = 3;
    padding_size = ceil((apothem_per_radius - 1) * lumen_radius_in_pixels_range(end,:) + 2);  
    
%     orig = double(h52mat(path_to_original_data));
%      size_of_image = size(orig);
%     orig = orig - min(orig(:)) + eps;
%     orig = orig ./ max(orig(:));
%     orig = padarray(orig,padding_size,'replicate','both');
    
    for vertex_set_string = vertex_set_list
        vertex_set_string = char(vertex_set_string);
        
        disp(['   - Running ', vertex_set_string, ' set'])
        
        switch vertex_set_string
            case 'uncurated'
                current_feature_path = path_to_vertex_features;
                current_vertex_vectors_path = path_to_vertices;
            case 'curated'
                current_feature_path = path_to_curated_vertex_features;    
                current_vertex_vectors_path = path_to_curated_vertices;
        end
        
        vertexFeaturePool = load_vertex_feature_pool(current_feature_path);  
        load(current_vertex_vectors_path)
        
        original_space_subscripts = double(vertex_space_subscripts);
%         vertex_space_subscripts = double(vertex_space_subscripts) + padding_size;
        
        save(current_feature_path, 'vertexFeaturePool')
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        best_resolution_allowed = 1 / 4.5 ;

        % construct the directories for storing the blurred data (and the chunked intermediary data).
        number_of_scales = length( lumen_radius_in_microns_range );

        scales_per_octave = log(                                    2                                   ) ...
                          / log( lumen_radius_in_microns_range( 2 ) / lumen_radius_in_microns_range( 1 )) ...
                          / 3 ; % divide by three for the volume interperetation of octave

        original_file_info = h5info( path_to_original_data );

        size_of_image = original_file_info.Datasets.Dataspace.Size ;

        pixels_per_radius_range         = lumen_radius_in_microns_range    ./ microns_per_voxel ;
        vessel_wall_thickness_in_pixels = vessel_wall_thickness_in_microns ./ microns_per_voxel ;
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        center_of_original_image = ((size_of_image - 1) / 2) + 1;
        
%        vertex_space_subscripts = double(vertex_space_subscripts) - padding_size;
        
        vertexFeaturePool.vertex_space_subscripts = max((double(vertex_space_subscripts) - center_of_original_image) / (size_of_image - 1),2);
        vertexFeaturePool.vertex_radii            = exp(interp1(log(lumen_radius_in_microns_range), vertex_scale_subscripts)); %in microns        
     
        [vertexFeaturePool.vertex_derivative_subscripts, std_across_kernels, ~, mean_across_kernels] ...
                        = getVertexDerivatives(original_space_subscripts, vertex_scale_subscripts, best_resolution_allowed, ...
                          number_of_scales, lumen_radius_in_microns_range, scales_per_octave,      ...
                          size_of_image, pixels_per_radius_range, vessel_wall_thickness_in_pixels, ...
                          vertexFeaturePool.vertex_radii, ...
                          vertex_set_string, microns_per_voxel, pixels_per_sigma_PSF, path_to_original_data, ...
                          padding_size, apothem_per_radius, matching_kernel_string, vessel_wall_thickness_in_microns,...
                          start_time);
                                                                          
        vertexFeaturePool.vertex_energies         = vertex_energies ./ std_across_kernels;  
        
        %average vertex intensity z-score
        vertexFeaturePool.vertex_intensity = (mean_across_kernels - mean(mean_across_kernels)) ./ std_across_kernels;
        
        save(current_feature_path, 'vertexFeaturePool')
        
    end
end %derivative info

if extracting_spatial_and_energy_info
%     
%     disp([newline 'Extracting Vertex Spatial and Energy Info']);
%     tic
%     
%     for vertex_set_string = vertex_set_list
%         vertex_set_string = char(vertex_set_string);
%         
%         disp(['   - Running ', vertex_set_string, ' set']);
%         
%         switch vertex_set_string
%             case 'uncurated'
%                 current_feature_path = path_to_vertex_features;
%                 current_vertex_vectors_path = path_to_vertices;
%             case 'curated'
%                 current_feature_path = path_to_curated_vertex_features;    
%                 current_vertex_vectors_path = path_to_curated_vertices;
%         end
%         
%         load_vertex_feature_pool(current_feature_path)
%         load(current_vertex_vectors_path)
% 
%         vertexFeaturePool.vertex_space_subscripts = double(vertex_space_subscripts) ./ max(double(vertex_space_subscripts));
%         vertexFeaturePool.vertex_scale_subscripts = vertex_scale_subscripts;       
%         vertexFeaturePool.vertex_energies = vertex_energies;       
%        
%         save(current_feature_path,'vertexFeaturePool');
%        
%     end
%     
%     fprintf('     '); toc
    
end %spatial info

%% synthesize intenity info
if extracting_intensity_statistics
% 
%     disp([newline 'Extracting Inensity Statistics Across Vertex Volumes'])
%     
%     for vertex_set_string = vertex_set_list
%         vertex_set_string = char(vertex_set_string);
%         
%         disp(['   - Running ', vertex_set_string, ' set'])
%         tic
%         
%         switch vertex_set_string
%             case 'uncurated'
%                 current_feature_path = path_to_vertex_features;
%                 current_vertex_vectors_path = path_to_vertices;
%             case 'curated'
%                 current_feature_path = path_to_curated_vertex_features;    
%                 current_vertex_vectors_path = path_to_curated_vertices;
%         end
%         
%         load_vertex_feature_pool(current_feature_path)   
%         load(current_vertex_vectors_path)
% 
%         size_of_image = h5info(path_to_original_data);
%         size_of_image = size_of_image.Datasets.Dataspace.Size;
% 
%         number_of_scales = length( lumen_radius_in_microns_range );
% 
%         scale_subscript_range = 1 : number_of_scales ;
% 
%         structuring_etlement_linear_indexing_templates = cell( number_of_scales, 1 );
% 
%         radii_in_pixels_range = lumen_radius_in_microns_range ./ microns_per_voxel ;
% 
%         for scale_subscript = scale_subscript_range
% 
%             % find all pixel locations within the ellipsoid radii from the vertex position    
%             structuring_element_linear_indexing_templates{ scale_subscript }                                            ...
%                 = int64( construct_structuring_element( radii_in_pixels_range( scale_subscript, : ), size_of_image ));
% 
%         end % constructing relative elements FOR scale
% 
%         % init structuring elements
%         space_subscripts_int64 = int64( vertex_space_subscripts );
% 
%         vertex_position_linear_indices =   space_subscripts_int64( :, 1 )                                               ...
%                                        + ( space_subscripts_int64( :, 2 ) - 1 ) * size_of_image( 1 )                    ...
%                                        + ( space_subscripts_int64( :, 3 ) - 1 ) * size_of_image( 1 ) * size_of_image( 2 );
% 
%         vertex_position_linear_indices_cell = num2cell( vertex_position_linear_indices );
% 
%         structuring_element_linear_indexing = structuring_element_linear_indexing_templates( round( vertex_scale_subscripts ))';
% 
%         vertex_structure_positions_linear_indexing = cellfun( @plus, structuring_element_linear_indexing,   ...
%                                                                        vertex_position_linear_indices_cell, ...
%                                                                                      'UniformOutput', false );   
%         % end initialize_structuring_elements
% 
%         orig = double(h52mat(path_to_original_data));
%         orig = orig - min(orig(:)) + eps;
%         orig = orig ./ max(orig(:));
%         orig = reshape(orig,[],1);
% 
%         average_intensity_across_original =  mean(orig,'all');
%          std_of_intensity_across_original = std(orig,0,'all');
         
%         vertexFeaturePool.std_of_intensity_across_vertices =  cellfun(@(x)  std(orig(x)), vertex_structure_positions_linear_indexing, 'UniformOutput', true);
%         vertexFeaturePool.average_intensity_across_vertices = cellfun(@(x) mean(orig(x))./ std(orig(x)), vertex_structure_positions_linear_indexing, 'UniformOutput', true);

        
%          vertexFeaturePool.normalized_mean_intensity_across_vertices   = average_intensity_across_vertices ./ average_intensity_across_original;
%          vertexFeaturePool.normalized_std_of_intensity_across_vertices =  std_of_intensity_across_vertices ./  std_of_intensity_across_original;

%         save(current_feature_path, 'vertexFeaturePool');
%         
%     end
%     
%     fprintf('     '); toc

end %intensity info

%% Only for extracting training data      
if extracting_for_training
    disp([newline 'Extracting Vertex Training Data'])
    
    tic
    
    load(path_to_saved_curation)
    
    vertex_truth_table = true_vertices & displayed_vertices & ~ deleted_vertices;
    vertex_truth_table = vertex_truth_table(max_number_of_added_vertices + 2:end);
    
    displayed_vertices_for_extraction = displayed_vertices(max_number_of_added_vertices + 2:end);
    
    load(path_to_vertex_features)
    
    vertexFeaturePool.truth_table = vertex_truth_table;
    vertexFeaturePool = structfun(@(x) trimFeature(x,displayed_vertices_for_extraction), vertexFeaturePool, 'UniformOutput', false);
    
    save(path_to_training_vertex_features, 'vertexFeaturePool');
    
    fprintf('     '); toc
    
end

%% keep these next to each other. end of the script
end % for roi_index range
end % for directory_and_batch_index

disp('Done with vertex feature extraction... ');
toc(overallStartTime);

end % titular function 

%% Helper Functions

% trims vertex lists for training
function output = trimFeature(x,y) 
    output = x(y,:);
end


function [      batch_directory, ...
                 data_directory, ...
          visual_data_directory, ...
               vector_directory, ...
        visual_vector_directory, ...
             curation_directory, ...
             settings_directory  ]    = get_directories( batch_directory )
        
         data_directory = [ batch_directory, 'data',           filesep ];
  visual_data_directory = [ batch_directory, 'visual_data',    filesep ];
       vector_directory = [ batch_directory, 'vectors',        filesep ];
visual_vector_directory = [ batch_directory, 'visual_vectors', filesep ];
     curation_directory = [ batch_directory, 'curations',      filesep ];
     settings_directory = [ batch_directory, 'settings',       filesep ];

end % FUNCTION get_directories

function vertexFeaturePool = load_vertex_feature_pool(path_to_features)
    if isfile(path_to_features)
        load(path_to_features)
    else
        vertexFeaturePool = struct();
    end
end %load_vertex_feature_pool
