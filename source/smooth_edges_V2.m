function [ edge_space_subscripts, edge_scale_subscripts, edge_energies ] = smooth_edges_V2( edge_space_subscripts, edge_scale_subscripts, edge_energies, smoothing_kernel_sigma_to_lumen_radius_ratio, lumen_radius_in_microns_range, microns_per_voxel )
%% smooth_edges
% smooths edges using a gaussian kernel and energy weighting.  Outputs double precision inter pixel
% and inter scale edge coordinates about one voxel apart (anisotropically in real space).  SAM
% 4/17/19

% get the (energy-weighted) average radius of each edge in units of microns

edge_index_range = 1 : size( edge_scale_subscripts, 1 );

% changed from PARFOR to FOR 6/8/20
for edge_index = edge_index_range
    
    % set the energy of the manually added portion to be equal to the best energy elsewhere
    is_inf_position = edge_energies{ edge_index } == - Inf ;
        
    if any( ~ is_inf_position )
        
        edge_energies{ edge_index }( is_inf_position ) = min( edge_energies{ edge_index }( ~ is_inf_position ));        

    end
end

edge_scale_subscript_averages = cellfun( @( x, y ) sum( x .* y ) / sum( y ), edge_scale_subscripts, edge_energies );

edge_average_lumen_radii = exp( interp1( log( lumen_radius_in_microns_range ), edge_scale_subscript_averages ));

microns_per_sigma = smoothing_kernel_sigma_to_lumen_radius_ratio * edge_average_lumen_radii ;
    
edge_subscripts = cellfun( @( u, v ) [ double( u ), double( v )], edge_space_subscripts, edge_scale_subscripts, 'UniformOutput', false );

% lumen_radius_in_pixels_range = lumen_radius_in_microns_range ./ microns_per_voxel ;

% (PAR)FOR smooth along the edges
for edge_index = edge_index_range

    % IF this edge was not entirely manually added
    if any( edge_energies{ edge_index } > - Inf )
        
        edge_subscripts_at_edge = edge_subscripts{ edge_index };
        edge_energies_at_edge   = edge_energies{   edge_index };

        % build gaussian kernels using real spatial distances      

        % convert to microns in each spatial dimension
        edge_microns_at_edge = edge_subscripts_at_edge( :, 1 : 3 ) .* microns_per_voxel ;

        % precalculate cumulative distances covered by the vectors in each edge for later interpolation
        edge_cumulative_length                                                                          ...
            = [ 0; cumsum( sum((   edge_microns_at_edge( 1 + 1 : end    , 1 : 3 )                       ...
                                 - edge_microns_at_edge( 1     : end - 1, 1 : 3 )) .^ 2, 2 ) .^ 0.5 )]; ...                        

        kernel_micron_domains = edge_cumulative_length  ...
                              - edge_cumulative_length' ;

        kernel_sigma_domains = kernel_micron_domains / microns_per_sigma( edge_index );                      

        gaussian_kernels = exp( - kernel_sigma_domains .^ 2 / 2 );

        energy_conv_kernel = sum( edge_energies_at_edge .* gaussian_kernels, 1 );

    %  indicator_conv_kernel = sum(                         gaussian_kernels, 1 );

    % calculate energy weigthed average of the energy (more representative of the energy of the
    % subscripts that contribute to each new smoothed (energy-weighted) subscript location
        energy_conv_energy_conv_kernel = sum( edge_energies_at_edge .^ 2 .* gaussian_kernels, 1 );

        gaussian_kernels = permute( gaussian_kernels, [ 1, 3, 2 ]);        

        subscript_conv_energy_conv_kernel = squeeze( sum( edge_subscripts_at_edge .* edge_energies_at_edge .* gaussian_kernels, 1 ));    

        edge_subscripts_at_edge_smoothed = subscript_conv_energy_conv_kernel ./ energy_conv_kernel ;

        edge_subscripts_at_edge_smoothed = edge_subscripts_at_edge_smoothed' ;

    %     edge_energies_at_edge_smoothed = energy_conv_kernel ./ indicator_conv_kernel ;

        edge_energies_at_edge_smoothed = energy_conv_energy_conv_kernel ./ energy_conv_kernel ;

        edge_energies_at_edge_smoothed = edge_energies_at_edge_smoothed' ;    

        % Force the first location to its original (before smoothing)
        edge_subscripts_at_edge_smoothed(  1 , : ) = edge_subscripts{ edge_index }(  1 , : );
        edge_subscripts_at_edge_smoothed( end, : ) = edge_subscripts{ edge_index }( end, : );

        edge_energies_at_edge_smoothed(  1 , : ) = edge_energies{ edge_index }(  1  );
        edge_energies_at_edge_smoothed( end, : ) = edge_energies{ edge_index }( end );    

        % Interpolate without changing the number of total points.  This is to have regular sampling at
        % about 1 voxel length between locations (anisotroptic in real space).  This is the best
        % sampling strategy for rendering purposes in the current resolution.  To change the resoution,
        % use the function RESAMPLE_VECTORS.
        edge_cumulative_lengths                                                                                 ...
            = [ 0; cumsum( max( abs(   edge_subscripts_at_edge_smoothed( 1 + 1 : end    , 1 : 3 )               ...
                                     - edge_subscripts_at_edge_smoothed( 1     : end - 1, 1 : 3 )), [ ], 2 ))]; ...  


        edge_sample_lengths = linspace( 0, edge_cumulative_lengths( end ), numel( edge_cumulative_lengths ))' ;

        edge_subscripts_at_edge_smoothed                 ...
            = interp1( edge_cumulative_lengths,          ...
                       edge_subscripts_at_edge_smoothed, ...
                       edge_sample_lengths               );

        edge_energies_at_edge_smoothed                   ...
            = interp1( edge_cumulative_lengths,          ...
                       edge_energies_at_edge_smoothed,   ...
                       edge_sample_lengths               );

        edge_space_subscripts{ edge_index } = edge_subscripts_at_edge_smoothed( :, 1 : 3 );
        edge_scale_subscripts{ edge_index } = edge_subscripts_at_edge_smoothed( :,   4   );
        edge_energies{         edge_index } = edge_energies_at_edge_smoothed              ;
        
    
    end % IF not manually added edge
end % FOR edge index

end % FUNCTION