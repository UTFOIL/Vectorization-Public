function [ edge_space_subscripts, edge_scale_subscripts, edge_energies ] = smooth_edges_V2( edge_space_subscripts, edge_scale_subscripts, edge_energies, smoothing_kernel_sigma_to_lumen_radius_ratio, lumen_radius_in_microns_range, microns_per_voxel )
%% smooth_edges
% smooths edges using a gaussian kernel and energy weighting.  Outputs double precision inter pixel
% and inter scale edge coordinates about one voxel apart (anisotropically in real space).  SAM
% 4/17/19

is_smoothing_to_underlying_image = ischar( edge_energies ); % smoothing the edge positions to the underlying energy image % SAM 2/6/22

edge_index_range = 1 : size( edge_scale_subscripts, 1 );

% get the (energy-weighted) average radius of each edge in units of microns
if is_smoothing_to_underlying_image
    
%     path_to_energy_data = [ data_directory, energy_handle ];
    path_to_energy_data = edge_energies ;
    
    edge_energies = cell( size( edge_space_subscripts ));
    
    energy_and_size = h52mat( path_to_energy_data );
        
    energy_image = energy_and_size( :, :, :, 2 );
      size_image = energy_and_size( :, :, :, 1 );
    
    size_of_image = size( energy_image );      
      
    edge_space_subscripts_mat = round( double( cell2mat( edge_space_subscripts )));
    edge_scale_subscripts_mat =                cell2mat( edge_scale_subscripts );
    
    edge_numels        = cellfun( @numel, edge_scale_subscripts );
    
    edge_energies_mat      = zeros( size( edge_scale_subscripts_mat ));
      
    pos_index_range        =   1 : numel( edge_scale_subscripts_mat ) ;
    
    edge_locations                                    ...
        = sub2ind( size_of_image,                     ...
                   edge_space_subscripts_mat( :, 1 ), ...
                   edge_space_subscripts_mat( :, 2 ), ...
                   edge_space_subscripts_mat( :, 3 ));
               
    [ strel_linear_LUT_range, ~, ...
        cum_prod_image_dims, ...
        local_subscripts_range, ...
        strel_r_over_R_LUT_range, ...
        strel_unit_vectors_LUT_range ] ...
       = calculate_linear_strel_range( size_of_image, microns_per_voxel, lumen_radius_in_microns_range );

    % (parfor)
    for pos_index = pos_index_range
    
        edge_scale = round( size_image( edge_locations( pos_index )));        
        
        current_strel =                 edge_locations( pos_index ) ...
                      + strel_linear_LUT_range{ edge_scale };

        is_current_strel_in_map = current_strel >= 1                      ...
                                & current_strel <= cum_prod_image_dims( 3 );

        current_strel( ~ is_current_strel_in_map ) = edge_locations( pos_index ); % to be removed later

       [ current_subscript_y,                             ...
         current_subscript_x,                             ...
         current_subscript_z ] = ind2sub(  size_of_image, ...
                                       edge_locations( pos_index ));      current_edge_subscripts ...
                                                                      = [ current_subscript_y,    ...
                                                                          current_subscript_x,    ...
                                                                          current_subscript_z     ];
       
        absolute_subscripts_range = local_subscripts_range{ edge_scale } + current_edge_subscripts ;

        is_current_strel_in_map = logical( prod(   absolute_subscripts_range >= [ 1, 1, 1 ]             ...
                                                 & absolute_subscripts_range <= size_of_image, 2 ));

        current_strel( ~ is_current_strel_in_map ) = [ ];
        
        current_strel_unit_vectors       = strel_unit_vectors_LUT_range{ edge_scale };
        current_strel_r_over_R_LUT_range =     strel_r_over_R_LUT_range{ edge_scale };

        current_strel_unit_vectors(        ~ is_current_strel_in_map, : ) = [ ];
        current_strel_r_over_R_LUT_range(  ~ is_current_strel_in_map    ) = [ ];
        
        weighting______strel =          exp( - current_strel_r_over_R_LUT_range .^ 2 / 2 ) ... % Gaussian weighted at              current location (prior distribution)
                            .* energy_image(   current_strel                             ) ;   % energy   weighted in nbrhd around current location (post. distribution)
        
        weighting_pos__strel = lumen_radius_in_microns_range( edge_scale ) ...    %       radius estimate          at current position 
                            ./                 microns_per_voxel           ...       in units of voxel length in x, y, and z directions %
                            .*                 current_strel_unit_vectors       ... % unit      3-space  direction from current position
                            .*                 current_strel_r_over_R_LUT_range ... % unitless total (L2) distance from current position
                            .* weighting______strel ;
                       
        
        weighting_size_strel = size_image( current_strel ) ... % scale idx (log-dist. in radii) of current location
                            .*      weighting______strel   ;
                        
                edge_energies_mat( pos_index    ) = sum( energy_image( current_strel )  ) / sum( weighting______strel );
        edge_scale_subscripts_mat( pos_index    ) = sum(        weighting_size_strel    ) / sum( weighting______strel );
        edge_space_subscripts_mat( pos_index, : ) = sum(        weighting_pos__strel, 1 ) / sum( weighting______strel ) ...
                                                  + edge_space_subscripts_mat( pos_index, : );
                                              
    end % FOR pos index
    
%     smoothed_edge_energies         = mat2cell(         edge_energies_mat, edge_numels    );
%     smoothed_edge_scale_subscripts = mat2cell( edge_scale_subscripts_mat, edge_numels    );
%     smoothed_edge_space_subscripts = mat2cell( edge_space_subscripts_mat, edge_numels, 3 );
%     
%     edge_energies         = cellfun( @( x, y ) [ x( 1    ), y( 2 : end - 1    ), x( end    )], edge_energies,         smoothed_edge_energies,         'UniformOutput', false );
%     edge_scale_subscripts = cellfun( @( x, y ) [ x( 1    ), y( 2 : end - 1    ), x( end    )], edge_scale_subscripts, smoothed_edge_scale_subscripts, 'UniformOutput', false );
%     edge_space_subscripts = cellfun( @( x, y ) [ x( 1, : ), y( 2 : end - 1, : ), x( end, : )], edge_space_subscripts, smoothed_edge_space_subscripts, 'UniformOutput', false );
    
    edge_energies         = mat2cell(         edge_energies_mat, edge_numels    );
    edge_scale_subscripts = mat2cell( edge_scale_subscripts_mat, edge_numels    );
    edge_space_subscripts = mat2cell( edge_space_subscripts_mat, edge_numels, 3 );
    
else                       % smoothing the edge positions to neighboring positions 

    % (parfor)
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

    current_edge_subscripts = cellfun( @( u, v ) [ double( u ), double( v )], edge_space_subscripts, edge_scale_subscripts, 'UniformOutput', false );

    % lumen_radius_in_pixels_range = lumen_radius_in_microns_range ./ microns_per_voxel ;

    % (PAR)FOR smooth along the edges
    for edge_index = edge_index_range

        % IF this edge was not entirely manually added
        if any( edge_energies{ edge_index } > - Inf )

            edge_subscripts_at_edge = current_edge_subscripts{ edge_index };
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
            edge_subscripts_at_edge_smoothed(  1 , : ) = current_edge_subscripts{ edge_index }(  1 , : );
            edge_subscripts_at_edge_smoothed( end, : ) = current_edge_subscripts{ edge_index }( end, : );

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
end

end % FUNCTION