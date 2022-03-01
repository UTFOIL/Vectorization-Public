function [feature_array, truth_table] = getTrainingArray(path_to_training_data, object_name)
   
    listing = dir([path_to_training_data, filesep, '*.mat']);

    feature_array = [];
    truth_table = [];

    for i = 1:length(listing)
        load([path_to_training_data, filesep, listing(i).name],[ object_name, 'FeaturePool'])
        
        feature_list = fieldnames(vertexFeaturePool)';

        feature_array_mini = [];
        truth_table_mini = [];
        
        for current_feature_string = feature_list
           
            if ~strcmp(current_feature_string,'truth_table')

                toCat = eval([object_name, 'FeaturePool.', char(current_feature_string)]);
                if iscell(toCat)
                    toCat = cell2mat(toCat);
                end

                if numel(feature_array_mini) > 0
                    feature_array_mini = cat(2,feature_array_mini,toCat);
                else
                    feature_array_mini = toCat;
                end
                
            else
                
                toCat = eval([object_name, 'FeaturePool.', char(current_feature_string)]);
                if iscell(toCat)
                    toCat = cell2mat(toCat);
                end

                if numel(truth_table_mini) > 0
                    truth_table_mini = cat(2,truth_table_mini,toCat);
                else
                    truth_table_mini = toCat;
                end
                
            end

        end 
        
        feature_array_mini = double(feature_array_mini);
        
        feature_array = [feature_array; feature_array_mini];
        truth_table = [truth_table; truth_table_mini];
        
    end
    
end
