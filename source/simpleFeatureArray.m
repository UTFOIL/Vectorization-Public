function [feature_array] = simpleFeatureArray(featurePool)
    
    feature_list = fieldnames(featurePool)';
         
    feature_array = [];
    
    for current_feature_string = feature_list

        toCat = eval(['featurePool.', char(current_feature_string)]);
        if iscell(toCat)
            toCat = cell2mat(toCat);
        end
            
        if numel(feature_array) > 0
            feature_array = cat(2,feature_array,toCat);
        else
            feature_array = toCat;
        end
                
    end  
    
    feature_array = double(feature_array)';
    
end