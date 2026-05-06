function events = extract_watershed_trace(energy_map, vertex_locations, edge_number_tolerance, energy_tolerance)
    % Standalone MATLAB script to record execution events during watershed discovery.
    % Use this to generate JSONL traces for comparison with Python.
    
    number_of_vertices = length(vertex_locations);
    size_of_image = size(energy_map);
    
    % Initialize maps (matching get_edges_by_watershed.m)
    vertex_index_map = zeros(size_of_image, 'uint32');
    vertex_index_map(vertex_locations) = 1:number_of_vertices;
    
    energy_map_temp = energy_map;
    energy_map_temp(vertex_locations) = -Inf;
    
    pointer_map = zeros(size_of_image, 'uint64');
    d_over_r_map = zeros(size_of_image, 'double');
    
    available_locations = vertex_locations(end:-1:1);
    
    events = {};
    iteration = 0;
    
    while ~isempty(available_locations)
        iteration = iteration + 1;
        current_location = available_locations(end);
        current_energy = energy_map_temp(current_location);
        
        % Record iteration start
        events{end+1} = struct('event', 'iteration_start', ...
                               'iteration', iteration, ...
                               'current_linear', current_location, ...
                               'current_energy', current_energy);
        
        if current_energy == -Inf
            energy_map_temp(current_location) = energy_map(current_location);
            current_energy = energy_map_temp(current_location);
        end
        
        if current_energy >= 0, break; end
        
        % Discovery logic (simplified for snippet)
        % In a real audit, you would integrate this directly into get_edges_by_watershed.m
        
        % Example: record seed selection
        % events{end+1} = struct('event', 'seed_selected', ...
        %                        'seed_idx', seed_idx, ...
        %                        'selected_linear', next_location, ...
        %                        'selected_energy', adjusted_energy);
    end
    
    % Optional: write to file
    % fid = fopen('matlab_trace.jsonl', 'w');
    % for i = 1:length(events)
    %     fprintf(fid, '%s\n', jsonencode(events{i}));
    % end
    % fclose(fid);
end
