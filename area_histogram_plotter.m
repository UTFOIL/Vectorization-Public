function length_histogram_plotter( strand_subscripts, lumen_radius_in_microns_range, microns_per_voxel, size_of_image, number_of_bins )
%% length_histogram_plotter SAM 12/13/19
% This function calculates lenght-weighted histograms for radius, depth, and inlincation statistics.
% All strands are transformed to an inter-position representation with associated lengths,
% positions, etc. These inter-positions are combined into a single list and sorted by depth.
% Inter-positions are then binned according to statistic value. The total length of vessels in each
% bin are summed. Radius and depth are in microns, inclination in units of degrees.

%% inner-position strand quantities

lumen_radius_in_pixels_range = lumen_radius_in_microns_range ./ microns_per_voxel ;

z_pos_limits = [ 0, size_of_image( 3 ) * microns_per_voxel( 3 )];

% interpolate vectors in z to be consistent with xy
resolution_factors = lumen_radius_in_pixels_range( 1, 1 ) ./ lumen_radius_in_pixels_range( 1, : );

[ size_of_image, ~, strand_subscripts ] = resample_vectors( lumen_radius_in_pixels_range, [ resolution_factors, 1 ], strand_subscripts, size_of_image );
                                        
microns_per_voxel = microns_per_voxel([ 1, 1, 1 ]);

% length by strand
strand_delta_lengths        = cellfun( @( x ) sum(((   x( 2 : end    , 1 : 3 )                                           ...
                                                     - x( 1 : end - 1, 1 : 3 )) .* microns_per_voxel ) .^ 2, 2 ) .^ 0.5, ...
                                       strand_subscripts, 'UniformOutput', false                                         );
                            
strand_radii                = cellfun( @( x ) exp( interp1( log( lumen_radius_in_microns_range ), x( :, 4 ))), ...
                                       strand_subscripts, 'UniformOutput', false                               );

% strand inner-position radii
strand_inner_pos_radii      = cellfun( @( x ) (   x( 2 : end     )          ...
                                                + x( 1 : end - 1 )) / 2,    ...
                                       strand_radii, 'UniformOutput', false );
                                                              
% strand inner-position z-position
strand_inner_pos_z          = cellfun( @( x ) (   x( 2 : end    , 3 )                                 ...
                                                + x( 1 : end - 1, 3 )) / 2 .* microns_per_voxel( 3 ), ...
                                       strand_subscripts, 'UniformOutput', false                      );
                                   
% strand inner-position area
strand_inner_pos_area      = cellfun( @( x, y ) 2 * pi * x .* y,                                            ...
                                       strand_inner_pos_radii, strand_delta_lengths, 'UniformOutput', false );
                               
% strand inner-position inclination
strand_inner_pos_inclinations                                                                          ...
                            = cellfun( @( x, delta_l )                                                 ...
                                    asind( abs(   x( 2 : end    , 3 )                                  ...
                                                - x( 1 : end - 1, 3 ))     .* microns_per_voxel( 3 )   ...
                                                                           ./ delta_l               ), ...
                                       strand_subscripts, strand_delta_lengths, 'UniformOutput', false );
                                                              
%% decomposing strands into a single list of inner-positions, sorted by z position

% decomposing:
% delta_lengths           = cell2mat( strand_delta_lengths           );
inner_pos_area          = cell2mat( strand_inner_pos_area          );
inner_pos_radii         = cell2mat( strand_inner_pos_radii         );
inner_pos_z             = cell2mat( strand_inner_pos_z             );
inner_pos_inclinations  = cell2mat( strand_inner_pos_inclinations  );

radius_limits = [ min( inner_pos_radii ), max( inner_pos_radii )];

inclin_limits = [ 0, 90 ];

%% sorting, binning
[ x1, dx1, dy1 ] = sort_and_bin_by_statistic( inner_pos_z,            inner_pos_area, number_of_bins,  z_pos_limits, false );
[ x2, dx2, dy2 ] = sort_and_bin_by_statistic( inner_pos_radii,        inner_pos_area, number_of_bins, radius_limits, true  );
[ x3, dx3, dy3 ] = sort_and_bin_by_statistic( inner_pos_inclinations, inner_pos_area, number_of_bins, inclin_limits, false );

 x = [  x1;  x2;  x3 ];
dx = [ dx1; dx2; dx3 ];
dy = [ dy1; dy2; dy3 ];

function [ x, dx, dy ] = sort_and_bin_by_statistic( inner_statistics, inner_pos_area, number_of_bins, statistic_limits, is_log_dist_stat )

    [ inner_statistics, sorted_indices ] = sort( inner_statistics );

    inner_pos_area = inner_pos_area( sorted_indices );
    

    if is_log_dist_stat
    
        bin_limits = exp( linspace( log( statistic_limits( 1 )), log( statistic_limits( 2 )), number_of_bins + 1 ));
                
    else
        
    %     bin_limits = linspace( min( inner_statistics ), max( inner_statistics ), number_of_bins + 1 );
        bin_limits = linspace( statistic_limits( 1 ), statistic_limits( 2 ), number_of_bins + 1 );
    
    end
    
    bin_index_limits = zeros( 1, number_of_bins + 1 );
    
    for ik = 1 : number_of_bins
        
        bin_index_limits( ik ) = find( inner_statistics >= bin_limits( ik ), 1, 'first' );
        
    end
    
    bin_index_limits( number_of_bins + 1 ) = numel( inner_statistics ) + 1 ; % left limit is included, right limit is not
    
    bin_numels = bin_index_limits( 2 : end     ) ...
               - bin_index_limits( 1 : end - 1 );

    bin_deltas = bin_limits( 2 : end     ) ...
               - bin_limits( 1 : end - 1 );

    inner_pos_area_bin_cell = mat2cell( inner_pos_area, bin_numels )';      

    % summing lengths
    bin_inner_pos_area_sums = cellfun( @sum, inner_pos_area_bin_cell );   
    
    x = bin_limits( 1 : end - 1 );

    dx = bin_deltas ;
    
    dy = bin_inner_pos_area_sums ;    
    
%     dy = dy / 1000 ; % convert to milimeters

end
%% plotting

figure

y = zeros( 1, length( x ));

for stat_i = 1 : 3 

    subplot( 1, 3, stat_i ), hold on, 

    if stat_i == 1, ylabel( 'area [um^2]' ); end

    box_plotter( x( stat_i, : ), y, dx( stat_i, : ), dy( stat_i, : ))
    
    xlim([ x( stat_i, 1 ), x( stat_i, end ) + dx( stat_i, end )])
    
    ylim([ 0, max( dy( : ))])
    
    switch stat_i
        case 1, xlabel( 'depth [um]' )
        case 2, xlabel( 'radius [um]' )
        case 3, xlabel( 'inclination [degrees]' )
    end
end

% subplot( 5, 2, 1 ), hold on, ylabel({   'area density', '[um^2/um^3]' }), dy = z_bin_SA_density ;
% box_plotter( x, y, dx, dy )
% xlim([ bin_limits( 1 ), bin_limits( end )])
% 
% subplot( 5, 3, 1 ), hold on, ylabel({ 'volume density', '[um^3/um^3]' }), dy = z_bin_vol_densities ;
% box_plotter( x, y, dx, dy )
% xlim([ bin_limits( 1 ), bin_limits( end )])
% 
% subplot( 5, 1, 4 ), hold on, ylabel({         'radius', '[um]'        }), dy = z_bin_ave_radius ;
% box_plotter( x, y, dx, dy )
% xlim([ bin_limits( 1 ), bin_limits( end )])
% 
% subplot( 5, 1, 5 ), hold on, ylabel({    'inclination', '[degrees]'   }), dy = z_bin_ave_inclination ;
% box_plotter( x, y, dx, dy )
% xlim([ bin_limits( 1 ), bin_limits( end )])



    function box_plotter( x, y, dx, dy )
        
        for ii=1:length(x), rectangle('position',[x(ii) y(ii) dx(ii) dy(ii)]), end
        
        set( gca, 'FontSize', 10 )
        
    end % FUNCTION box_plotter

end % FUNCTION