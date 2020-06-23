function [ z_bin_ave_z, z_bin_length_densities, z_bin_SA_density, z_bin_vol_densities, z_bin_ave_radius, z_bin_ave_inclination ] = calculate_depth_statistics( strand_subscripts, lumen_radius_in_microns_range, microns_per_voxel, size_of_image, number_of_bins )
%% calculate_network_statistics SAM 5/25/19
% This function calculates depth statistics for the input vectorized network to yield depth-resolved
% statistics. All strands are transformed to an inter-position representation with associated
% lengths, positions, etc. These inter-positions are combined into a single list and sorted by
% depth.  Inter-positions are then binned according to linearly spaced points in depth.  Bins are
% normalized by the volume of image they represent. All statistics are in dimensions of microns or
% some power of microns. Averaged quanties like radius and inclination are calculated as a surface
% area weighted average.

%% inner-position strand quantities

lumen_radius_in_pixels_range = lumen_radius_in_microns_range ./ microns_per_voxel ;

% interpolate vectors in z to be consistent with xy
resolution_factors = lumen_radius_in_pixels_range( 1, 1 ) ./ lumen_radius_in_pixels_range( 1, : );

[ size_of_image, ~, strand_subscripts ] ...
                                            = resample_vectors( lumen_radius_in_pixels_range, [ resolution_factors, 1 ], strand_subscripts, size_of_image );
                                        
microns_per_voxel = microns_per_voxel([ 1, 1, 1 ]);

% length by strand
strand_delta_lengths        = cellfun( @( x ) sum(((   x( 2 : end    , 1 : 3 )                                           ...
                                                     - x( 1 : end - 1, 1 : 3 )) .* microns_per_voxel ) .^ 2, 2 ) .^ 0.5, ...
                                       strand_subscripts, 'UniformOutput', false                                         );
                            
% strand_shortest_length  = cellfun( @( x ) sum(((   x(         end, 1 : 3 )                                           ...
%                                                  - x(         1  , 1 : 3 )) .* microns_per_voxel ) .^ 2    ) .^ 0.5, ...
%                                    strand_subscripts                                                                 );

strand_radii                = cellfun( @( x ) exp( interp1( log( lumen_radius_in_microns_range ), x( :, 4 ))), ...
                                       strand_subscripts, 'UniformOutput', false                               );

% strand inner-position radii
strand_inner_pos_radii      = cellfun( @( x ) (   x( 2 : end     )          ...
                                                + x( 1 : end - 1 )) / 2,    ...
                                       strand_radii, 'UniformOutput', false );
                               
% strand inner-position radii squared
strand_inner_pos_radii_squ  = cellfun( @( x ) (   x( 2 : end     ) .^ 2        ...
                                                + x( 1 : end - 1 ) .^ 2 ) / 2, ...
                                       strand_radii, 'UniformOutput', false    );                               
                               
% strand inner-position z-position
strand_inner_pos_z          = cellfun( @( x ) (   x( 2 : end    , 3 )                                 ...
                                                + x( 1 : end - 1, 3 )) / 2 .* microns_per_voxel( 3 ), ...
                                       strand_subscripts, 'UniformOutput', false                      );

% strand inner-position inclination
strand_inner_pos_inclinations                                                                          ...
                            = cellfun( @( x, delta_l )                                                 ...
                                    asind( abs(   x( 2 : end    , 3 )                                  ...
                                                - x( 1 : end - 1, 3 ))     .* microns_per_voxel( 3 )   ...
                                                                           ./ delta_l               ), ...
                                       strand_subscripts, strand_delta_lengths, 'UniformOutput', false );
                               
xy_area = prod( size_of_image([ 1, 2 ]) .* microns_per_voxel([ 1, 2 ]));
                               
% !!!!!! input vertex locations and count how many bifurcations in each z-bin

%% decomposing strands into a single list of inner-positions, sorted by z position

% decomposing:
delta_lengths           = cell2mat( strand_delta_lengths           );
inner_pos_radii         = cell2mat( strand_inner_pos_radii         );
inner_pos_radii_squ     = cell2mat( strand_inner_pos_radii_squ     );
inner_pos_z             = cell2mat( strand_inner_pos_z             );
inner_pos_inclinations  = cell2mat( strand_inner_pos_inclinations  );

% sorting:
[ inner_pos_z, sorted_indices ] = sort( inner_pos_z );

delta_lengths           =           delta_lengths( sorted_indices );
inner_pos_radii         =         inner_pos_radii( sorted_indices );
inner_pos_radii_squ     =     inner_pos_radii_squ( sorted_indices );
inner_pos_inclinations  =  inner_pos_inclinations( sorted_indices );

%% binning in z
number_of_inner_positions = length( delta_lengths );

% bin_quantile_limits = linspace( 0, 1, number_of_bins + 1 );
% 
% bin_index_limits = round( bin_quantile_limits * number_of_inner_positions + 1 ); % left limit is included, right limit is not

bin_limits = linspace( inner_pos_z( 1 ), inner_pos_z( end ), number_of_bins + 1 );

bin_index_limits = zeros( 1, number_of_bins + 1 );

for ik = 1 : number_of_bins

    bin_index_limits( ik ) = find( inner_pos_z >= bin_limits( ik ), 1, 'first' );

end

bin_index_limits( number_of_bins + 1 ) = numel( inner_pos_z ) + 1 ; % left limit is included, right limit is not

bin_numels = bin_index_limits( 2 : end     ) ...
           - bin_index_limits( 1 : end - 1 );

% bin_index_range = 1 : number_of_bins ;

bin_z_limits = inner_pos_z( min( bin_index_limits, number_of_inner_positions ));

bin_delta_z = bin_z_limits( 2 : end     ) ...
            - bin_z_limits( 1 : end - 1 );
        
      delta_lengths_bin_cell = mat2cell(       delta_lengths, bin_numels );
    inner_pos_radii_bin_cell = mat2cell(     inner_pos_radii, bin_numels );
inner_pos_radii_squ_bin_cell = mat2cell( inner_pos_radii_squ, bin_numels );
 inner_inclinations_bin_cell = mat2cell(  inner_pos_inclinations, bin_numels );
        inner_pos_z_bin_cell = mat2cell(         inner_pos_z, bin_numels );

z_bin_length_densities  = cellfun( @( delta_l )           sum( delta_l              ),                       delta_lengths_bin_cell )                                        ./ bin_delta_z / xy_area ;
z_bin_SA_density        = cellfun( @( delta_l, r )        sum( delta_l .* r         ),                       delta_lengths_bin_cell, inner_pos_radii_bin_cell )     * 2 * pi ./ bin_delta_z / xy_area ;
z_bin_vol_densities     = cellfun( @( delta_l, r_squ )    sum( delta_l .* r_squ     ),                       delta_lengths_bin_cell, inner_pos_radii_squ_bin_cell ) *     pi ./ bin_delta_z / xy_area ;
z_bin_ave_z             = cellfun( @( delta_l, z )        sum( delta_l .* z         ) / sum( delta_l )     , delta_lengths_bin_cell, inner_pos_z_bin_cell                                  );
z_bin_ave_radius        = cellfun( @( delta_l, r )        sum( delta_l .* r .^ 2    ) / sum( r .* delta_l ), delta_lengths_bin_cell, inner_pos_radii_bin_cell                              );
z_bin_ave_inclination   = cellfun( @( delta_l, r, incl )  sum( delta_l .* r .* incl ) / sum( r .* delta_l ), delta_lengths_bin_cell, inner_pos_radii_bin_cell, inner_inclinations_bin_cell );

%% plotting

y = zeros( length( z_bin_ave_z ));

x = bin_z_limits( 1 : end - 1 );

dx = bin_delta_z ;

figure

subplot( 5, 1, 1 ), hold on, ylabel({ 'length density', '[um/um^3]'   }), dy = z_bin_length_densities ; 
box_plotter( x, y, dx, dy )
xlim([ bin_z_limits( 1 ), bin_z_limits( end )])

subplot( 5, 1, 2 ), hold on, ylabel({   'area density', '[um^2/um^3]' }), dy = z_bin_SA_density ;
box_plotter( x, y, dx, dy )
xlim([ bin_z_limits( 1 ), bin_z_limits( end )])

subplot( 5, 1, 3 ), hold on, ylabel({ 'volume density', '[um^3/um^3]' }), dy = z_bin_vol_densities ;
box_plotter( x, y, dx, dy )
xlim([ bin_z_limits( 1 ), bin_z_limits( end )])

subplot( 5, 1, 4 ), hold on, ylabel({         'radius', '[um]'        }), dy = z_bin_ave_radius ;
box_plotter( x, y, dx, dy )
xlim([ bin_z_limits( 1 ), bin_z_limits( end )])

subplot( 5, 1, 5 ), hold on, ylabel({    'inclination', '[degrees]'   }), dy = z_bin_ave_inclination ;
box_plotter( x, y, dx, dy )
xlim([ bin_z_limits( 1 ), bin_z_limits( end )])


xlabel( 'depth [um]' )

    function box_plotter( x, y, dx, dy )
        
        for ii=1:length(x), rectangle('position',[x(ii) y(ii) dx(ii) dy(ii)]), end
        
        set( gca, 'FontSize', 10 )
        
    end % FUNCTION box_plotter

end % FUNCTION