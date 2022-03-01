function [ edge_space_subscripts, edge_scale_subscripts, edge_energies ] = smooth_edges( edge_space_subscripts, edge_scale_subscripts, edge_energies, edge_lengths, smoothing_kernel_sigma_to_lumen_radius_ratio, lumen_radius_in_pixels_range )
%% smooth_edges
% smooths edges using a gaussian kernel and energy weighting.  Outputs double precision inter pixel
% and inter scale edge coordinates at equallyl spaced intervals in real space.  SAM 4/17/19

% get the average radius of each edge in units of a characteristic voxel length

% take geo mean of the voxel side lengths to get a characteristic voxel length.
lumen_radius_in_voxel_lengths_range = prod( lumen_radius_in_pixels_range, 2 ) .^ ( 1 / 3 );    

edge_scale_subscript_averages = cellfun( @mean, edge_scale_subscripts );

edge_average_lumen_radii_in_voxel_lengths = exp( interp1( log( lumen_radius_in_voxel_lengths_range ), edge_scale_subscript_averages ));

indices_per_sigma                                                                             ...
        = smoothing_kernel_sigma_to_lumen_radius_ratio * edge_average_lumen_radii_in_voxel_lengths ;
    
radius_of_smoothing_kernels = round( 3 * indices_per_sigma );

padded_edge_lengths = 2 * radius_of_smoothing_kernels + double( edge_lengths );

edge_index_range = 1 : size( edge_scale_subscripts, 1 );

edge_subscripts = cellfun( @( u, v ) [ double( u ), double( v )], edge_space_subscripts, edge_scale_subscripts, 'UniformOutput', false );

% smooth along the edges
for edge_index = edge_index_range

    kernel_index_domain        = ( 1 : padded_edge_lengths( edge_index ))' ;
    
    kernel_index_domain_offset = ( 1 : padded_edge_lengths( edge_index ))  ;
    
    kernel_index_domains = kernel_index_domain        ...
                         - kernel_index_domain_offset ;
                     
    kernel_sigma_domains = kernel_index_domains / indices_per_sigma( edge_index );

    % Pad the edge by repeating the first and last entries before convolution. Pad with the radius
    % of the smoothing kernel, so that at the start of the convolution, the result is at least 99.7
    % % guaranteed to coincide with the origin vertex.
    edge_subscripts_padded = zeros( padded_edge_lengths( edge_index ), 4 );

    edge_subscripts_padded( 1 :   radius_of_smoothing_kernels( edge_index ),     : ) =  edge_subscripts{ edge_index }(  1, :  )             ...
                                                                                     .* ones( radius_of_smoothing_kernels( edge_index ), 4 );

    edge_subscripts_padded(       radius_of_smoothing_kernels( edge_index ) + 1                                                             ...
                          : end - radius_of_smoothing_kernels( edge_index )    , : ) =  edge_subscripts{ edge_index } ;

    edge_subscripts_padded( end - radius_of_smoothing_kernels( edge_index ) + 1                                                             ...
                              : end                                            , : ) =  edge_subscripts{ edge_index }( end, : )             ...
                                                                                     .* ones( radius_of_smoothing_kernels( edge_index ), 4 );
                                                                                         
    edge_energies_padded = zeros( padded_edge_lengths( edge_index ), 1 );
    
    edge_energies_padded( 1 : radius_of_smoothing_kernels( edge_index ))       = edge_energies{ edge_index }( 1 );
    
    edge_energies_padded(     radius_of_smoothing_kernels( edge_index ) + 1                                ...
                      : end - radius_of_smoothing_kernels( edge_index )     )  = edge_energies{ edge_index };
                      
    edge_energies_padded( end - radius_of_smoothing_kernels( edge_index ) + 1                                     ...
                          : end                                              ) = edge_energies{ edge_index }( end );
                                                                                                     
    gaussian_kernels = exp( - kernel_sigma_domains .^ 2 / 2 );

               energy_conv_kernel = sum( edge_energies_padded .* gaussian_kernels, 1 );
               
%             indicator_conv_kernel = sum(                         gaussian_kernels, 1 );
            
% energy weigthed average of the energy (more represented of the energy of the subscripts that are
% contributing to each new smoothed (energy-weighted) subscript location
   energy_conv_energy_conv_kernel = sum( edge_energies_padded .^ 2 .* gaussian_kernels, 1 );
            
    gaussian_kernels = permute( gaussian_kernels, [ 1, 3, 2 ]);        
    
subscript_conv_energy_conv_kernel = squeeze( sum( edge_subscripts_padded .* edge_energies_padded .* gaussian_kernels, 1 ));    

    edge_subscripts_at_edge_smoothed = subscript_conv_energy_conv_kernel ./ energy_conv_kernel ;
    
    edge_subscripts_at_edge_smoothed = edge_subscripts_at_edge_smoothed' ;

%     edge_energies_at_edge_smoothed = energy_conv_kernel ./ indicator_conv_kernel ;

    edge_energies_at_edge_smoothed = energy_conv_energy_conv_kernel ./ energy_conv_kernel ;
    
    edge_energies_at_edge_smoothed = edge_energies_at_edge_smoothed' ;    
    
    % guarantee that the first location coincides with the first vertex
    edge_subscripts_at_edge_smoothed(  1 , : ) = edge_subscripts{ edge_index }(  1 , : );
    edge_subscripts_at_edge_smoothed( end, : ) = edge_subscripts{ edge_index }( end, : );
    
    edge_energies_at_edge_smoothed(  1 , : ) = edge_energies{ edge_index }(  1  );
    edge_energies_at_edge_smoothed( end, : ) = edge_energies{ edge_index }( end );    
        
    % Interpolate back to the original number of vectors as the original edge

%     % convert to a common real distance unit in each spatial dimension
%     edge_subscripts_at_edge_smoothed( :, 1 : 3 ) = edge_subscripts_at_edge_smoothed( :, 1 : 3 ) ...
%                                                   ./ lumen_radius_in_pixels_range( 1, : );

    % precalculate cumulative distances covered by the vectors in each edge for later interpolation
%     edge_cumulative_length                                                                                      ...
%         = [ 0; cumsum( sum((   edge_subscripts_at_edge_smoothed( 1 + 1 : end    , 1 : 3 )                       ...
%                              - edge_subscripts_at_edge_smoothed( 1     : end - 1, 1 : 3 )) .^ 2, 2 ) .^ 0.5 )]; ...  

    edge_cumulative_lengths                                                                                 ...
        = [ 0; cumsum( max( abs(   edge_subscripts_at_edge_smoothed( 1 + 1 : end    , 1 : 3 )               ...
                                 - edge_subscripts_at_edge_smoothed( 1     : end - 1, 1 : 3 )), [ ], 2 ))]; ...  

                             
    edge_sample_lengths = linspace( 0, edge_cumulative_lengths( end ), edge_lengths( edge_index ))' ;

    edge_subscripts_at_edge_smoothed                 ...
        = interp1( edge_cumulative_lengths,          ...
                   edge_subscripts_at_edge_smoothed, ...
                   edge_sample_lengths               );
                               
    edge_energies_at_edge_smoothed                   ...
        = interp1( edge_cumulative_lengths,          ...
                   edge_energies_at_edge_smoothed,   ...
                   edge_sample_lengths               );

%     % convert back to units of voxels
%     edge_subscripts_at_edge_smoothed( :, 1 : 3 ) = edge_subscripts_at_edge_smoothed( :, 1 : 3 ) ...
%                                                   .* lumen_radius_in_pixels_range( 1, : );
                                              
    edge_space_subscripts{ edge_index } = edge_subscripts_at_edge_smoothed( :, 1 : 3 );
    edge_scale_subscripts{ edge_index } = edge_subscripts_at_edge_smoothed( :,   4   );
    edge_energies{         edge_index } = edge_energies_at_edge_smoothed              ;
    
end % FOR edge index

end % FUNCTION