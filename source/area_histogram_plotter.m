function area_histogram_plotter( strand_subscripts, lumen_radius_in_microns_range, microns_per_voxel, size_of_image, number_of_bins, ROI_name, network_handle )
%% length_histogram_plotter SAM 12/13/19
% This function calculates lenght-weighted histograms for radius, depth, and inlincation statistics.
% All strands are transformed to an inter-position representation with associated lengths,
% positions, etc. These inter-positions are combined into a single list and sorted by depth.
% Inter-positions are then binned according to statistic value. The total length of vessels in each
% bin are summed. Radius and depth are in microns, inclination in units of degrees.

is_smoothing = true ;

% % depth_bins = {[0, Inf]};
% depth_bins = {[0, Inf], [0, 100], [100, 200], [200, Inf]};
% 
% number_depth_bins = length( depth_bins );
% 
% colors = [ 0,   0,    0;
%            0.5, 0.75, 1;
%            1, 0.5, 0.75;
%            0.75, 1, 0.5 ];

%% inner-position strand quantities

lumen_radius_in_pixels_range = lumen_radius_in_microns_range ./ microns_per_voxel ;

% z_pos_limits = [ 0, size_of_image( 3 ) * microns_per_voxel( 3 )];
strand_subscripts_mat = cell2mat( strand_subscripts );

% z_pos_limits = [ min( strand_subscripts_mat( :, 3 )), max( strand_subscripts_mat( :, 3 ))] ...
%              * microns_per_voxel( 3 );

z_pos_limits = [ 0, max( strand_subscripts_mat( :, 3 ))] ...
             * microns_per_voxel( 3 );
         
clear strand_subscripts_mat

% interpolate vectors to be isotropic consistent with the dimension with the best resolution
resolution_factors = max( lumen_radius_in_pixels_range( 1, : )) ./ lumen_radius_in_pixels_range( 1, : );

[ size_of_image, ~, strand_subscripts ] = resample_vectors( lumen_radius_in_pixels_range, [ resolution_factors, 1 ], strand_subscripts, size_of_image );
                                        
microns_per_voxel = [ 1, 1, 1 ] * min( microns_per_voxel([ 1, 2 ]));

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
                                                + x( 1 : end - 1, 3 )) / 2  * microns_per_voxel( 3 ), ...
                                       strand_subscripts, 'UniformOutput', false                      );
                                   
% strand inner-position area
strand_inner_pos_area      = cellfun( @( x, y ) 2 * pi * x .* y,                                            ...
                                       strand_inner_pos_radii, strand_delta_lengths, 'UniformOutput', false );
                               
% strand inner-position inclination
strand_inner_pos_inclinations                                                                          ...
                            = cellfun( @( x, delta_l )                                                 ...
...                                    asind( abs(   x( 2 : end    , 3 )                                  ...
                                           abs(   x( 2 : end    , 3 )                                  ...
                                                - x( 1 : end - 1, 3 ))     .* microns_per_voxel( 3 )   ...
                                                                           ./ delta_l               ,  ... ),
                                       strand_subscripts, strand_delta_lengths, 'UniformOutput', false );
                                                              
%% decomposing strands into a single list of inner-positions, sorted by z position

% decomposing:
% delta_lengths           = cell2mat( strand_delta_lengths           );
inner_pos_area          = cell2mat( strand_inner_pos_area          );
inner_pos_radius         = cell2mat( strand_inner_pos_radii         );
inner_pos_z_pos             = cell2mat( strand_inner_pos_z             );
inner_pos_inclin  = cell2mat( strand_inner_pos_inclinations  );

%% sorting, binning
if is_smoothing, inner_pos_radius = log( inner_pos_radius ) / log( 10 ); end

radius_limits = [ min( inner_pos_radius ) - eps, ...
                  max( inner_pos_radius )        ]; % - eps becuase error was thrown when the smallest radius object didn't make it into the 1st bin because of rounding error during exp/log transform SAM 12/1/21
inclin_limits = [ 0, 1 ];

if is_smoothing
    
    z_pos_limits = [ min( inner_pos_z_pos ), ...
                     max( inner_pos_z_pos )  ];
    
     z_pos_resolution = (  z_pos_limits( 2 ) -  z_pos_limits( 1 )) / number_of_bins ;
    radius_resolution = ( radius_limits( 2 ) - radius_limits( 1 )) / number_of_bins ;
    inclin_resolution = ( inclin_limits( 2 ) - inclin_limits( 1 )) / number_of_bins ;    
    
    [ x1, pdf1, cdf1 ] = smooth_hist( inner_pos_z_pos ,  z_pos_resolution, inner_pos_area );
    [ x2, pdf2, cdf2 ] = smooth_hist( inner_pos_radius, radius_resolution, inner_pos_area );
    [ x3, pdf3, cdf3 ] = smooth_hist( inner_pos_inclin, inclin_resolution, inner_pos_area );
     
    pdf = [  pdf1,  pdf2,  pdf3 ]';
    cdf = [  cdf1,  cdf2,  cdf3 ]';
      x = [    x1,    x2,    x3 ]';
      y = cdf ;
     dy = pdf ;
     dx = zeros( size( x ));     
    
    else % binning
    
    [ x1, y1, dx1, dy1 ] = sort_and_bin_by_statistic( inner_pos_z_pos,  inner_pos_area, number_of_bins,  z_pos_limits, false );
    [ x2, y2, dx2, dy2 ] = sort_and_bin_by_statistic( inner_pos_radius, inner_pos_area, number_of_bins, radius_limits, true  );
    [ x3, y3, dx3, dy3 ] = sort_and_bin_by_statistic( inner_pos_inclin, inner_pos_area, number_of_bins, inclin_limits, false );

     x = [  x1;  x2;  x3 ];
     y = [  y1;  y2;  y3 ];
    dx = [ dx1; dx2; dx3 ];
    dy = [ dy1; dy2; dy3 ];

end


function [ x, y, dx, dy ] = sort_and_bin_by_statistic( inner_statistics, inner_pos_area, number_of_bins, statistic_limits, is_log_dist_stat )

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
    
    bin_index_limits( 1 ) = 1 ; % force left limit to be the first index
    
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
    
    y = cumsum( bin_inner_pos_area_sums );
    
%     dy = dy / 1000 ; % convert to milimeters

end
%% plotting

h = figure;

set( h, 'Name', [ network_handle, '_', ROI_name, ' Area Histograms' ])

y_zeros = zeros( 1, length( x ));

for stat_i = 1 : 3 

    for histogram_style_index = 1 : 2    

%         subplot( 1, 3, stat_i ), hold on, 
        subplot( 2, 3, ( histogram_style_index - 1 ) * 3 + stat_i )

        if is_smoothing
            if stat_i == 2 % log x axis
                fxn_plotter( 10 .^ x( stat_i, : ), pdf( stat_i, : ), cdf( stat_i, : ), histogram_style_index )
            else % linear x axis
                fxn_plotter(       x( stat_i, : ), pdf( stat_i, : ), cdf( stat_i, : ), histogram_style_index )
            end
        end
            
        if histogram_style_index == 2 % count

            if ~ is_smoothing
                if stat_i == 1, ylabel( 'area [um^2]' ); end
            else
                if stat_i == 1, ylabel( 'area [um^2] / abscissa (dB) units' ); end
            end

            if ~ is_smoothing, box_plotter( x( stat_i, : ), y_zeros, dx( stat_i, : ),  dy( stat_i, : )), end

            if ~ is_smoothing, ylim([ 0, max( dy( stat_i, : ))]), end

        else % cumulative

            if stat_i == 1, ylabel( 'area [um^2]' ); end                        

            if ~ is_smoothing, box_plotter( x( stat_i, : ), y_zeros, dx( stat_i, : ),  y( stat_i, : )), end

            ylim([ 0, max(  y( : ))])            

        end 
            
        switch stat_i
            case 1, xlabel( 'depth [um]' )
            case 2

                if is_smoothing, x( stat_i, : ) = 10 .^ x( stat_i, : ); end
                
                x_limits = [ x( stat_i, 1 ), x( stat_i, end ) + dx( stat_i, end )];

                xlim( x_limits )

                max_x_log_10 =      floor( log( x_limits( 2 )) / log( 10 ))                    ;    
                min_x_log_10 = min(  ceil( log( x_limits( 1 )) / log( 10 )), max_x_log_10 - 1 );

                x_limits_log_10 = [ min_x_log_10, max_x_log_10 ];

                x_ticks = 10 .^ x_limits_log_10 ;    
                
                xlabel( 'radius [um]' ) ...

                set( gca,'xscale','log',                                                ...
                     'XTickLabel', {[ '10^{', num2str( x_limits_log_10( 1 )), '}' ],    ...
                                    [ '10^{', num2str( x_limits_log_10( 2 )), '}' ]  }, ...                    
                     'XTick',      x_ticks,                                             ...
                     'Xlim',       [ min( x_ticks( 1 ), x_limits( 1 )), x_limits( 2 )]   )
                 
                if is_smoothing, x( stat_i, : ) = log( x( stat_i, : )) / log( 10 ); end

            case 3, xlabel( 'inclination [z component]' ), set( gca, 'Xlim', [ 0, 1 ])
        end
        
        set( gca, 'TickDir', 'both' )
        
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

function fxn_plotter( x, pdf, cdf, histogram_style_index )
    if histogram_style_index == 1, plot( x, cdf )
    else,                          plot( x, pdf )
...                      hold on, H = area( x, pdf ); set(H(1),'FaceColor',[0.5 0.75 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5 );
                      hold on, H = area( x, pdf ); set(H(1),'FaceColor',[0.5 0.75 1], 'EdgeColor', 'none', 'FaceAlpha', 0.5 );                      
    end
end

end % FUNCTION