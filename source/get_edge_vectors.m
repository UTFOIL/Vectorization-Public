function [ edges2vertices, edge_lengths, edge_space_subscripts, edge_scale_subscripts, edge_energies ] ...
       = get_edge_vectors( number_of_edges_per_vertex, edge_indices_temp, edges2vertices, path_to_energy_data, cum_prod_image_dims )

% load the energy data and associated scale index data (two 3D images concatenated into 4D array )
energy_and_index_data = h52mat( path_to_energy_data );

energy_data = energy_and_index_data( :, :, :, 2 );

 index_data = energy_and_index_data( :, :, :, 1 );
 
clear( 'energy_and_index_data' )

% extracting outputs from the edges2vertices and edge_indices cell arrays (whose format was
% temporarilly different to allow for the parrallel FOR to run):
edges2vertices = cell2mat( edges2vertices );

% looking for nonzero in the first position (means an edge was found for that place holder)
found_edge_index_range = find( edges2vertices( :, 1 ))';

% removing placeholders that don't hold edges
edges2vertices = edges2vertices( found_edge_index_range, : );

number_of_edges = numel( found_edge_index_range );

edge_indices = cell( number_of_edges, 1 );

vertex_indices = floor(( found_edge_index_range - 1 ) / number_of_edges_per_vertex ) + 1 ;

edge_indices_at_vertices = found_edge_index_range - ( vertex_indices - 1 ) * number_of_edges_per_vertex ;

edge_index_range = 1 : number_of_edges ;

for edge_index = edge_index_range
    
    edge_indices{ edge_index } = int64([]);
    
    edge_indices{ edge_index } = nonzeros( edge_indices_temp{ vertex_indices( edge_index )}( edge_indices_at_vertices( edge_index ), : ));
    
end % FOR found edge

clear( 'edge_indices_temp' )

edge_lengths = cellfun( @( x ) uint16( numel( x )), edge_indices ); % number of positions visited (units of voxels)

% convert the linear indexing to spatial subscripts
% !!!! replace with single call to index2position
% edge_space_subscripts_xy  = cellfun( @( t )  rem( t - 1 ,    cum_prod_image_dims( 2 )) + 1, edge_indices,                                      'UniformOutput', false );
% 
% edge_space_subscripts_z   = cellfun( @( t, xy ) ( t - xy ) / cum_prod_image_dims( 2 )  + 1, edge_indices, edge_space_subscripts_xy,            'UniformOutput', false );
% 
% edge_space_subscripts_y   = cellfun( @( xy ) rem( xy - 1,    cum_prod_image_dims( 1 )) + 1, edge_space_subscripts_xy,                          'UniformOutput', false );
% 
% edge_space_subscripts_x   = cellfun( @( xy, y ) ( xy - y ) / cum_prod_image_dims( 1 )  + 1, edge_space_subscripts_xy, edge_space_subscripts_y, 'UniformOutput', false );
% 
% clear( 'edge_space_subscripts_xy' )
% 
% edge_space_subscripts     = cellfun( @( y, x, z ) uint16([ y, x, z ]), edge_space_subscripts_y, edge_space_subscripts_x, edge_space_subscripts_z, 'UniformOutput', false );

edge_space_subscripts     = cellfun( @( t ) uint16( index2position( t, cum_prod_image_dims )'), edge_indices, 'UniformOutput', false );

% clear( 'edge_space_subscripts_y', 'edge_space_subscripts_x', 'edge_space_subscripts_z' )

edge_scale_subscripts     = cellfun( @( x )  index_data( x ), edge_indices, 'UniformOutput', false );
               
edge_energies             = cellfun( @( x ) energy_data( x ), edge_indices, 'UniformOutput', false );

end % FUNCTION