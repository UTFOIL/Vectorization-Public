function [ sorting_indices ] = sort_edges( edge_energies )

% sort the output
metric_of_edge_energies = get_edge_metric( edge_energies ); % MAX function SAM 201027

[ ~, sorting_indices ] = sort( metric_of_edge_energies ); % 'ascend' by SORT default

end
