function [ edge_lengths, edge_space_subscripts, edge_scale_subscripts, edge_energies ] ...
       = get_edge_vectors_V300( edge_locations, energy_map, index_map, cum_prod_image_dims )

edge_lengths              = cellfun( @( x ) numel( x ),                                         edge_locations ); % number of positions visited (units of voxels)

edge_space_subscripts     = cellfun( @( t ) uint16( index2position( t, cum_prod_image_dims )'), edge_locations, 'UniformOutput', false );

edge_scale_subscripts     = cellfun( @( x )  index_map( x ),                                    edge_locations, 'UniformOutput', false );
edge_energies             = cellfun( @( x ) energy_map( x ),                                    edge_locations, 'UniformOutput', false );

end % FUNCTION