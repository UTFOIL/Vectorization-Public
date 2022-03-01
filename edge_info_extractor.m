function edge_info_extractor(edge_set_list, edge_type_list, root_directories_and_batch_names, target_directory, varargin)
% extracting_for_training determines if script will attempt to extract
% features of curated edge data

extracting_for_training = logical(numel(intersect(edge_set_list,{'training'}))); 
edge_set_list(strcmp(edge_set_list,'training')) = [];

extracting_intensity_statistics = logical(numel(intersect(edge_type_list,{'intensity'})));
extracting_spatial_and_energy_info = logical(numel(intersect(edge_type_list,{'original'})));
if logical(numel(intersect(edge_type_list,{'all'})))
    extracting_intensity_statistics = true;
    extracting_spatial_and_energy_info = true;
end

%% This scripts gets feature info of noncurated (and curated for training) edges for machine learning curation.
% WAS 6/30/2019. 

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

% time_stamp = root_directory_and_batch_name( end-12 : end );
       
workflow_handle = [ 'workflow_', time_stamp ];

path_to_workflow_settings = [ settings_directory, workflow_handle ];

load( path_to_workflow_settings, 'production_times' )  

original_data_handle = 'original' ;
 
           energy_handle = [  'energy_', production_times{ 2 }];
         vertices_handle = ['vertices_', production_times{ 3 }];
            edges_handle = ['edges_',    production_times{ 4 }];
        
         path_to_energy_settings   = [ settings_directory,            energy_handle ];
         path_to_vertices_settings = [ settings_directory,          vertices_handle ];
         path_to_edges_settings    = [ settings_directory,             edges_handle ];

if ~ isempty( varargin ), ROI_index_range = varargin{ 1 }; end
         
for ROI_index = ROI_index_range

    path_to_energy_data             = [     data_directory,               energy_handle, ROI_names{ ROI_index }]; % mat file path        
    path_to_original_data           = [     data_directory,        original_data_handle, ROI_names{ ROI_index }]; % h5  file path
    path_to_vertices                = [   vector_directory,             vertices_handle, ROI_names{ ROI_index }]; %  vectors path
    path_to_curated_vertices        = [   vector_directory, 'curated_', vertices_handle, ROI_names{ ROI_index }]; %  vectors path    
    path_to_edges                   = [   vector_directory,                edges_handle, ROI_names{ ROI_index }];  
    path_to_curated_edges           = [   vector_directory, 'curated_',    edges_handle, ROI_names{ ROI_index }]; %  vectors path    

    path_to_saved_curation          = [ curation_directory,                edges_handle, ROI_names{ ROI_index }]; % logicals path

    %new ones 
    path_to_edge_features          = [ target_directory, '\edges\uncurated\', edges_handle, '_featureSet_',          ROI_names{ ROI_index }];
    path_to_curated_edge_features  = [ target_directory, '\edges\curated\',   edges_handle, '_featureSet_CURATED_',  ROI_names{ ROI_index }];
    path_to_training_edge_features = [ target_directory, '\edges\training\',  edges_handle, '_featureSet_TRAINING_', ROI_names{ ROI_index }];
    
    load(path_to_energy_settings)
    load( path_to_edges_settings )
    
    if ~ isempty( varargin )
        path_to_edge_features = target_directory;
    end

%% spatial info
if extracting_spatial_and_energy_info
%     
%     disp([newline 'Extracting Edge Spatial and Energy Info']);
%     tic
%     
%     for edge_set_string = edge_set_list
%         edge_set_string = char(edge_set_string);
%         
%         disp(['   - Running ', edge_set_string, ' set']);
%         
%         switch edge_set_string
%             case 'uncurated'
%                 current_feature_path = path_to_edge_features;
%                 current_edge_vectors_path = path_to_edges;
%             case 'curated'
%                 current_feature_path = path_to_curated_edge_features;    
%                 current_edge_vectors_path = path_to_curated_edges;
%         end
%         
%         load_edge_feature_pool(current_feature_path)
%         load(current_edge_vectors_path)
% 
%         size_of_original_image = size(orig);
%         center_of_original_image = ((size_of_original_image - 1) / 2) + 1;
%         mean_edge_space_subscripts = cellfun(@(x) mean(x,1),edge_space_subscripts,'UniformOutput', false);
%         mean_edge_space_subscripts = cell2mat(mean_edge_space_subscripts);
%         mean_edge_space_subscripts = max((mean_edge_space_subscripts - center_of_original_image) / (size_of_original_image - 1),2);
%         
%         mean_edge_scale_subscripts = cellfun(@(x) mean(x,1),edge_scale_subscripts);
%         mean_edge_scale_subscripts = exp(interp1(log(lumen_radius_in_microns_range), mean_edge_scale_subscripts));
% %         std_edge_scale_subscripts = cellfun(@(x) std(x,1),edge_scale_subscripts);
%         
%         edgeFeaturePool.largest_scale_subscript = cellfun(@(x) max(x(:)),edge_scale_subscripts);
%         edgeFeaturePool.edge_scale_residuals = cellfun(@(x) [x(1), x(end)],edge_scale_subscripts,'UniformOutput', false);
%         edgeFeaturePool.edge_scale_residuals = cell2mat(edgeFeaturePool.edge_scale_residuals);
%         edgeFeaturePool.edge_scale_residuals = edgeFeaturePool.largest_scale_subscript - edgeFeaturePool.edge_scale_residuals;
%         edgeFeaturePool.edge_scale_residuals = sort(edgeFeaturePool.edge_scale_residuals,2,'descend');
%         
%         edgeFeaturePool.mean_edge_space_subscripts = mean_edge_space_subscripts;
%         edgeFeaturePool.mean_edge_scale_subscripts = mean_edge_scale_subscripts;
% %         edgeFeaturePool.std_edge_scale_subscripts = std_edge_scale_subscripts;
%         edgeFeaturePool.mean_edge_energies = mean_edge_energies;
%        
%         save(current_feature_path,'edgeFeaturePool');
%        
%     end
%     
%     fprintf('     '); toc
%     
end %spatial info

%% synthesize intenity info
if extracting_intensity_statistics

    disp([newline 'Extracting Edge Inensity Statistics Across Edge Volumes'])
    
    for edge_set_string = edge_set_list
        edge_set_string = char(edge_set_string);
        
        disp(['   - Running ', edge_set_string, ' set'])
        tic
        
        switch edge_set_string
            case 'uncurated'
                current_feature_path = path_to_edge_features;
                current_edge_vectors_path = path_to_edges;
            case 'curated'
                current_feature_path = path_to_curated_edge_features;    
                current_edge_vectors_path = path_to_curated_edges;
        end
        
        load_edge_feature_pool(current_feature_path)   
        load(current_edge_vectors_path)

        size_of_image = h5info(path_to_original_data);
        size_of_image = size_of_image.Datasets.Dataspace.Size;
        

        number_of_scales = length( lumen_radius_in_microns_range );

%         scale_subscript_range = 1 : number_of_scales ;
% 
%         vertex_element_linear_indexing_templates = cell( number_of_scales, 1 );
%           edge_element_linear_indexing_templates = cell( number_of_scales, 1 );
% 
%         vertex_radii_in_pixels_range =  1 ...
%                                      .* lumen_radius_in_microns_range  ...
%                                      ./ microns_per_voxel              ;
% 
%           edge_radii_in_pixels_range =  1   ...
%                                      .* lumen_radius_in_microns_range ...
%                                      ./ microns_per_voxel             ;
%                   
%         for scale_subscript = scale_subscript_range
%       
%             % find all pixel locations within the ellipsoid radii from the element position
%             vertex_element_linear_indexing_templates{ scale_subscript }                                                        ...
%                 = int64( construct_structuring_element( vertex_radii_in_pixels_range( scale_subscript, : ), size_of_image ));
% 
%               edge_element_linear_indexing_templates{ scale_subscript }                                                        ...
%                 = int64( construct_structuring_element(   edge_radii_in_pixels_range( scale_subscript, : ), size_of_image ));
%      
%         end % constructing edge and vertex element templates FOR scale
% 
%         
%         degrees_of_edges_uint_32 = uint32( cellfun( @length, edge_scale_subscripts ));
% 
%         edge_space_subscripts_mat =         cell2mat( edge_space_subscripts );
%         edge_scale_subscripts_mat = uint8(  cell2mat( edge_scale_subscripts ));        
%         
%         space_subscripts_int64 = int64( edge_space_subscripts_mat );
% 
%         position_linear_indices =   space_subscripts_int64( :, 1 )                                               ...
%                                 + ( space_subscripts_int64( :, 2 ) - 1 ) * size_of_image( 1 )                    ...
%                                 + ( space_subscripts_int64( :, 3 ) - 1 ) * size_of_image( 1 ) * size_of_image( 2 );
% 
%         sphere_position_linear_indices_cell = num2cell( position_linear_indices );
% 
%         structuring_element_linear_indexing = edge_element_linear_indexing_templates( edge_scale_subscripts_mat );
% 
%         sphere_structure_positions_linear_indexing = cellfun( @plus, structuring_element_linear_indexing, ...
%                                                                      sphere_position_linear_indices_cell, ...
%                                                                                   'UniformOutput', false  ); 
% 
% %         vertex_structure_positions_linear_indexing = structuring_element_linear_indexing_templates( edge_scale_subscripts_mat );
% 
%         sphere_structure_positions_linear_indexing_edge_cell = mat2cell( sphere_structure_positions_linear_indexing,  degrees_of_edges_uint_32, 1 );
% 
% %         end_vertices_structure_positions_linear_indexing_1 = cellfun( @( x ) x{  1  }, sphere_structure_positions_linear_indexing_edge_cell, 'UniformOutput', false );
% %         end_vertices_structure_positions_linear_indexing_2 = cellfun( @( x ) x{ end }, sphere_structure_positions_linear_indexing_edge_cell, 'UniformOutput', false );
% 
%         space_subscripts_int64 = int64( vertex_space_subscripts );
% 
%         position_linear_indices =   space_subscripts_int64( :, 1 )                                               ...
%                                 + ( space_subscripts_int64( :, 2 ) - 1 ) * size_of_image( 1 )                    ...
%                                 + ( space_subscripts_int64( :, 3 ) - 1 ) * size_of_image( 1 ) * size_of_image( 2 );
% 
%         vertex_position_linear_indices_cell = num2cell( position_linear_indices );
% 
%         structuring_element_linear_indexing = vertex_element_linear_indexing_templates( round( vertex_scale_subscripts ));
% 
%         vertex_structure_positions_linear_indexing = cellfun( @plus, structuring_element_linear_indexing,   ...
%                                                                        vertex_position_linear_indices_cell, ...
%                                                                                      'UniformOutput', false ); 
% 
% %         vertex_structure_positions_linear_indexing = [ end_vertices_structure_positions_linear_indexing_1, ...
% %                                                        end_vertices_structure_positions_linear_indexing_2  ];
% 
%         edge_structure_positions_linear_indexing = cellfun( @cell2mat, sphere_structure_positions_linear_indexing_edge_cell, 'UniformOutput', false );
% 
%         edge_structure_positions_linear_indexing = cellfun( @( x ) unique( x, 'rows' ), edge_structure_positions_linear_indexing, 'UniformOutput', false );
% 
%         edge_structure_positions_subscript_xy = cellfun( @( x )      1 + mod( x - 1,  prod( size_of_image([ 1, 2 ]))), edge_structure_positions_linear_indexing, 'UniformOutput', false );
%         
%         edge_structure_positions_subscript_z  = cellfun( @( x ) floor( double( x ) ./ prod( size_of_image([ 1, 2 ]))), edge_structure_positions_linear_indexing, 'UniformOutput', false );
%         
%         edge_structure_positions_subscript_y  = cellfun( @( x )      1 + mod( x - 1,        size_of_image( 1 )),       edge_structure_positions_subscript_xy,    'UniformOutput', false );
%         
%         edge_structure_positions_subscript_x  = cellfun( @( x ) floor( double( x ) ./       size_of_image( 1 )),       edge_structure_positions_subscript_xy,    'UniformOutput', false );        
%                 
%         edge_structure_positions_subscripts = cellfun( @( y, x, z ) [ y, x, z ], edge_structure_positions_subscript_y,                        ...
%                                                                                  edge_structure_positions_subscript_x,                        ...
%                                                                                  edge_structure_positions_subscript_z, 'UniformOutput', false );
%                                             
%         edge_centers = cellfun( @( x ) round( mean( x, 1 )), edge_structure_positions_subscripts, 'UniformOutput', false );
% 
%         edge_centers = cell2mat( edge_centers );  
%                 % end initialize_structuring_elements
% 
% %         orig = double(h52mat(path_to_original_data));
% %         orig = orig - min(orig(:)) + eps('double');
% %         orig = orig ./ max(orig(:));
% %         size_of_original_image = size(orig);
% %         orig = reshape(orig,[],1);
        orig = double(h52mat(path_to_original_data));
        orig = orig - min(orig(:)) + eps('double');
        size_of_original_image = size(orig);
        clear orig
%         orig = orig ./ max(orig(:));
%         orig = reshape(orig,[],1);
% 
% %        average_intensity_across_original =  mean(orig,'all');
% %         std_of_intensity_across_original = std(orig,0,'all');
% 
% %        average_intensity_across_edges = cellfun(@(x) mean(orig(x)), edge_structure_positions_linear_indexing, 'UniformOutput', true);
%          std_of_intensity_across_edges = cellfun(@(x)  std(orig(x)), edge_structure_positions_linear_indexing, 'UniformOutput', true);
% 
% %        edgeFeaturePool.normalized_mean_intensity_across_edges   = average_intensity_across_edges ./ average_intensity_across_original;
% %        edgeFeaturePool.normalized_std_of_intensity_across_edges =  std_of_intensity_across_edges ./  std_of_intensity_across_original;

   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        center_of_original_image = ((size_of_original_image - 1) / 2) + 1;
        mean_edge_space_subscripts = cellfun(@(x) mean(x,1),edge_space_subscripts,'UniformOutput', false);
        mean_edge_space_subscripts = cell2mat(mean_edge_space_subscripts);
        mean_edge_space_subscripts = max((mean_edge_space_subscripts - center_of_original_image) / (size_of_original_image - 1),2);
        
        mean_edge_scale_subscripts = cellfun(@(x) mean(x,1),edge_scale_subscripts);
        mean_edge_scale_subscripts = exp(interp1(log(lumen_radius_in_microns_range), mean_edge_scale_subscripts));
        
        max_edge_scale_subscripts = cellfun(@(x) max(x(:)),edge_scale_subscripts);
        max_edge_radii = exp(interp1(log(lumen_radius_in_microns_range), max_edge_scale_subscripts));
        min_edge_scale_subscripts = cellfun(@(x) min(x(:)),edge_scale_subscripts);
        min_edge_radii = exp(interp1(log(lumen_radius_in_microns_range), min_edge_scale_subscripts));
        median_edge_scale_subscripts = cellfun(@(x) median(x(:)),edge_scale_subscripts);
        median_edge_radii = exp(interp1(log(lumen_radius_in_microns_range), median_edge_scale_subscripts));
        edgeFeaturePool.edge_scale_difference = (max_edge_radii - min_edge_radii) ./ median_edge_radii;
        
%         std_edge_scale_subscripts = cellfun(@(x) std(x,1),edge_scale_subscripts);
        
%         edgeFeaturePool.largest_scale_subscript = cellfun(@(x) max(x(:)),edge_scale_subscripts);
%         edgeFeaturePool.edge_scale_residuals = cellfun(@(x) [x(1), x(end)],edge_scale_subscripts,'UniformOutput', false);
%         edgeFeaturePool.edge_scale_residuals = cell2mat(edgeFeaturePool.edge_scale_residuals);
%         edgeFeaturePool.edge_scale_residuals = edgeFeaturePool.largest_scale_subscript - edgeFeaturePool.edge_scale_residuals;
%         edgeFeaturePool.edge_scale_residuals = sort(edgeFeaturePool.edge_scale_residuals,2,'descend');
        
%         edgeFeaturePool.mean_edge_space_subscripts = mean_edge_space_subscripts;
%         edgeFeaturePool.mean_edge_scale_subscripts = exp(interp1(log(lumen_radius_in_microns_range), mean_edge_scale_subscript));
%         edgeFeaturePool.std_edge_scale_subscripts = std_edge_scale_subscripts;

        edgeFeaturePool.mean_edge_energies = mean_edge_energies;
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
        if iscell(current_feature_path)
            current_feature_path = char(cell2mat(current_feature_path));
        end
    
        save(current_feature_path,'edgeFeaturePool');
            
%% Derivatives and Original 

        %Get vessel directions
        sigma_edge_smoothing_post = sqrt(max( 0, 0.5 - sigma_edge_smoothing^2));
        if sigma_edge_smoothing_post
            [ edge_space_subscripts_smoothed, ~ , ~ ]                               ...
                                         = smooth_edges_V2( edge_space_subscripts, edge_scale_subscripts, ...
                                                            edge_energies,                                ...
                                                            sigma_edge_smoothing_post, ...
                                                            lumen_radius_in_microns_range, microns_per_voxel );
        else
            edge_space_subscripts_smoothed = edge_space_subscripts;
        end
        edge_space_positions_in_microns = cellfun(@(x) x .* microns_per_voxel, edge_space_subscripts_smoothed, 'UniformOutput', false);
        edge_directions = cellfun(@get_edge_directions , edge_space_positions_in_microns, 'UniformOutput', false);
        edge_directions = cell2mat(edge_directions);
        
        edge_vertex_counts = cellfun(@(x) numel(x), edge_scale_subscripts);
        
        edge_vertex_scale_subscripts = cell2mat(edge_scale_subscripts);
        edge_vertex_space_subscripts = cell2mat(edge_space_subscripts);

%        edge_vertex_space_subscripts = max((edge_vertex_space_subscripts - center_of_original_image) / (size_of_image - 1),2);
%         edge_vertex_radii            = exp(interp1(log(lumen_radius_in_microns_range), edge_vertex_scale_subscripts)); %in microns        

        disp([newline 'Extracting Derivative Information'])
        start_time = tic;
    
        apothem_per_radius = 3;
        padding_size = ceil((apothem_per_radius - 1) * lumen_radius_in_pixels_range(end,:) + 2);  
        vertex_set_string = 'edge-vertex';
        
        original_space_subscripts = double(edge_vertex_space_subscripts);
    %     vertex_space_subscripts = double(vertex_space_subscripts) + padding_size;

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

%         center_of_original_image = ((size_of_image - 1) / 2) + 1;

%        edge_vertex_space_subscripts = double(edge_vertex_space_subscripts) - padding_size;

%         edge_vertex_space_subscripts = max((edge_vertex_space_subscripts - center_of_original_image) / (size_of_image - 1),2);
        edge_vertex_radii            = exp(interp1(log(lumen_radius_in_microns_range), edge_vertex_scale_subscripts)); %in microns        

        [edge_vertex_derivatives, std_across_kernels, kernel_maxs_by_scale, mean_across_kernels] = getVertexDerivatives(original_space_subscripts, edge_vertex_scale_subscripts, best_resolution_allowed, ...
                                                                              number_of_scales, lumen_radius_in_microns_range, scales_per_octave,      ...
                                                                              size_of_image, pixels_per_radius_range, vessel_wall_thickness_in_pixels, ...
                                                                              edge_vertex_radii, ...
                                                                              vertex_set_string, microns_per_voxel, pixels_per_sigma_PSF, path_to_original_data, ...
                                                                              padding_size, apothem_per_radius, matching_kernel_string, vessel_wall_thickness_in_microns,...
                                                                              start_time);

        variances_across_kernels = std_across_kernels .^ 2;
        weighted_std_across_kernels = variances_across_kernels ./ kernel_maxs_by_scale(round(edge_vertex_scale_subscripts))';
        weighted_std_across_edges = mat2cell(weighted_std_across_kernels,edge_vertex_counts);
        kernel_maxs_across_edges  = mat2cell(kernel_maxs_by_scale(round(edge_vertex_scale_subscripts))',edge_vertex_counts);
        std_across_edges = cellfun(@(x,y) sum(x ./ y), weighted_std_across_edges, kernel_maxs_across_edges);
        
        % Projection of derivative vectors onto cross section plane of lumen
        cross_products_of_vessel_directions_and_gradients = cross(edge_directions,edge_vertex_derivatives(:,1:3));
        transverse_component_of_gradient = sqrt(sum(cross_products_of_vessel_directions_and_gradients .^ 2, 2));
        transverse_gradients = mat2cell(transverse_component_of_gradient, edge_vertex_counts);
        edgeFeaturePool.mean_transverse_gradients = cellfun(@(x) mean(x,1), transverse_gradients) ./ std_across_edges;
        
%        edge_vertex_derivatives = mat2cell(edge_vertex_derivatives, edge_vertex_counts, 9 );
        
        % Curvature projections
        
        % Rotating basis to align z to vessel axis.
        % Get axis of rotation as a cross product of z unit vector and
        % vessel direction 
        axes_of_rotation = [-edge_directions(:,2), edge_directions(:,1), zeros(size(edge_directions,1),1)];
        % Make unit length
        axes_of_rotation = axes_of_rotation ./ (sqrt(sum(axes_of_rotation(:,1:2).^2,2)) + eps());
        
        % Get angles of rotation using dot product and cosine
        angles_of_rotation = acos(edge_directions(:,3));
        
        axes_of_rotation = mat2cell(axes_of_rotation,ones(size(axes_of_rotation,1),1));
        rotation_matrices = cellfun(@(x,y) cos(y)*eye(3)                                         ... 
                                           + sin(y)*([0 -x(3) x(2); x(3) 0 -x(1); -x(2) x(1) 0]) ...
                                           + (1 - cos(y))*(x' * x),                              ...
                                   axes_of_rotation, num2cell(angles_of_rotation),               ...
                                   'UniformOutput',false                                         );
                               
        transverse_curvatures = cellfun(@(p,g) p' * g([1 4 5; 4 2 6; 5 6 3]) * p, rotation_matrices, mat2cell(edge_vertex_derivatives(:,4:9),ones(size(edge_vertex_derivatives,1),1)), 'UniformOutput', false);

        transverse_curvatures = cell2mat( cellfun(@(x) x([1 2 3 4 5 6 7 8 9]), transverse_curvatures, 'UniformOutput', false));

%         rotated_curvatures = cat(3,rotated_curvatures{:});
        
        transverse_curvatures = mat2cell(transverse_curvatures,edge_vertex_counts);
        
        edgeFeaturePool.mean_transverse_curvatures = cellfun(@(x) mean(x,1), transverse_curvatures, 'UniformOutput', false);
        edgeFeaturePool.mean_transverse_curvatures = cell2mat(edgeFeaturePool.mean_transverse_curvatures)  ./ std_across_edges;
        
        %normalized intensity range across edges
        vertex_means_across_edges = mat2cell(mean_across_kernels,edge_vertex_counts);
        edgeFeaturePool.intensity_range = cellfun(@(x) (max(x) - min(x) ./ median(x)), vertex_means_across_edges, 'UniformOutput', false);
        edgeFeaturePool.intensity_range = cell2mat(edgeFeaturePool.intensity_range);
        
        if strcmp( edge_set_string, 'uncurated' )
            edgeFeaturePool = rmfield(edgeFeaturePool,'mean_edge_energies');
        end
        
        save(current_feature_path,'edgeFeaturePool');
        
        end

        fprintf('     '); toc

end %intensity info

%% Only for extracting training data      
if extracting_for_training
    disp([newline 'Extracting Edge Training Data'])
    
    tic
    
    load(path_to_saved_curation)
    
    edge_truth_table = true_edges & displayed_edges & ~ deleted_edges;
    edge_truth_table = edge_truth_table(max_number_of_added_edges + 2:end);
     
    displayed_edges_for_extraction = displayed_edges(max_number_of_added_edges + 2:end);
    
    if iscell(path_to_edge_features)
        path_to_edge_features = char(cell2mat(path_to_edge_features));
    end
    load(path_to_edge_features)
    
    edgeFeaturePool.truth_table = edge_truth_table;
    edgeFeaturePool = structfun(@(x) x(displayed_edges_for_extraction,:), edgeFeaturePool, 'UniformOutput', false);
    
    infinite_energy_indices = edgeFeaturePool.mean_edge_energies ~= -Inf;
    edgeFeaturePool = structfun(@(x) x(infinite_energy_indices,:), edgeFeaturePool, 'UniformOutput', false);
    
    edgeFeaturePool = rmfield(edgeFeaturePool,'mean_edge_energies');
    
    if iscell(path_to_training_edge_features)
        path_to_training_edge_features = char(cell2mat(path_to_training_edge_features));
    end
    save(path_to_training_edge_features, 'edgeFeaturePool');
    
    fprintf('     '); toc
    
end

%% keep these next to each other. end of the script
end % for roi_index range
end % for directory_and_batch_index
end % titular function 

%% Helper Functions

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

function load_edge_feature_pool(path_to_features)
    if isfile(path_to_features)
        load(path_to_features)
    else
        edgeFeaturePool = struct();
    end
end %load_edge_feature_pool

function edge_directions = get_edge_directions( edge_space_subs )

    [size_of_input(1), size_of_input(2)] = size( edge_space_subs );

    if size_of_input(1) >= 3
        edge_directions = edge_space_subs( 3 : end , : ) - edge_space_subs( 1 : end - 2 , : ) ;
        edge_directions = [ edge_directions( 1 , : ) ; edge_directions ; edge_directions( end , : ) ]; 
    else
        edge_directions = edge_space_subs( 2 , : ) - edge_space_subs( 1 , : ) ; 
        edge_directions = [edge_directions; edge_directions]; 
    end
    
    %Ensure unit length of direction vectors
    edge_directions = edge_directions ./ sqrt(sum(edge_directions .^ 2, 2));

end