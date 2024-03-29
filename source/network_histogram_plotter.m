function network_histogram_plotter( network_statistics, statistics_to_plot, number_of_bins, ROI_name, network_handle, plot_color ) %#ok<INUSL>

is_box_histo = false ; % IF is_box_histo == false, is_approx_PDF = true ; end

histogram_styles = { 'cumcount', 'count' };

if ~ isempty( statistics_to_plot )

    fields = statistics_to_plot ;

else % ELSE use default
    
    if isfield( network_statistics, 'strand_ave_d' ) % stats in reference to an origin (or axis)

        fields = { 'strand_lengths', ...
     ...              'strand_areas', ...
     ...              'strand_volumes', ...
                   'strand_tortuosity', ...
                   'strand_ave_d', ... % distance
                   'strand_ave_radii', ...
                   'strand_ave_r_component' };

    else % ELSE use depth and z direction stats

        fields = { 'strand_lengths', ...
...                   'strand_areas', ...
   ...                'strand_volumes', ...
                   'strand_tortuosity', ...
                   'strand_depths', ...
                   'strand_ave_radii', ...
                   'strand_z_direction'  };

    end
end % IF varargin is nonempty

number_of_statistics = length( fields );

statistic_indices = 1 : number_of_statistics ;

% h = figure;

set( gcf, 'Name', [ network_handle, '_', ROI_name, ' Strand Histograms' ])

for statistic_index = statistic_indices

%     network_statistic = eval([ 'network_statistics.', fields{ statistic_index }, ';' ]);
    network_statistic = network_statistics.( fields{ statistic_index });
    
    if strcmp( fields{ statistic_index }, 'strand_tortuosity' ), network_statistic = 1 ./ network_statistic ; end
            
    for histogram_style_index = 1 : 2

        subplot( 2, number_of_statistics, ( histogram_style_index - 1 ) * number_of_statistics + statistic_index )
    
        min_stat = min( network_statistic );
        max_stat = max( network_statistic );        
        
        switch fields{ statistic_index }

            case { 'strand_depths', ...
                   'strand_z_direction', ...
                   'strand_tortuosity', ...
                   'strand_ave_d', ...
                   'strand_ave_r_component' } % statistics that have linear scale on x axis
%             case 'strand_z_direction' % statistics that have linear scale on x axis
                
                if is_box_histo

                    edges = linspace( 0, 1, round( number_of_bins ));

                    histogram( network_statistic, edges, 'Normalization', histogram_styles{ histogram_style_index })
                    
                else % is_approx_PDF

                    [ x, pdf, cdf ] = smooth_hist( network_statistic, ( max_stat - min_stat ) / round( number_of_bins ));

                    if histogram_style_index == 1, plot( x, cdf,                 'Color', plot_color ), hold on
                    else,                          plot( x, pdf,                 'Color', plot_color ), hold on
                                               H = area( x, pdf ); set(H(1), 'FaceColor', plot_color, 'EdgeColor', 'none', 'FaceAlpha', 0.1 );
                    end

                end
                
                if histogram_style_index == 1, set( gca, 'XTickLabel', [ ]), end
                                
%                 set( gca, 'xlim',       [ min_stat, max_stat ])                
                
                switch fields{ statistic_index }

                    case { 'strand_depths', ...
                           'strand_ave_d'   }, set( gca, 'xlim', [ min_stat, max_stat ])
                    otherwise,                 set( gca, 'xlim', [        0, 1        ])

                end

            otherwise % log scale on x axis

%                 min_stat = min( network_statistic );
%                 max_stat = max( network_statistic );

                        max_stat_log_10 = floor( log( max_stat ) / log( 10 ));
                        min_stat_log_10 =  ceil( log( min_stat ) / log( 10 ));                
                   
                if max_stat_log_10 - min_stat_log_10 == -1

                        max_stat_log_10 =  ceil( log( max_stat ) / log( 10 ));
                        min_stat_log_10 = floor( log( min_stat ) / log( 10 ));

                end

                if max_stat_log_10 - min_stat_log_10 == 0

                    if max_stat - max_stat_log_10 ...
                     < min_stat_log_10 - min_stat

                        max_stat_log_10 =  ceil( log( max_stat ) / log( 10 ));
                        min_stat_log_10 =  ceil( log( min_stat ) / log( 10 ));

                    else

                        max_stat_log_10 = floor( log( max_stat ) / log( 10 ));
                        min_stat_log_10 = floor( log( min_stat ) / log( 10 ));

                    end

                end

                X_ticks_log_10 = min_stat_log_10 : max_stat_log_10 ;
                
                X_ticks   = 10 .^ X_ticks_log_10 ;
                
                if is_box_histo

                    [ ~, edges ] = histcounts( log( network_statistic )/log( 10 ), number_of_bins );

                    histogram( network_statistic, 10 .^ edges, 'Normalization', histogram_styles{ histogram_style_index })
                
                else % is_approx_PDF
                
                    [ x, pdf, cdf ] = smooth_hist( log( network_statistic ) / log( 10 ), ( log( max_stat ) - log( min_stat )) / log( 10 ) / number_of_bins );
                    
                    if histogram_style_index == 1, plot( 10 .^ x, cdf,                 'Color', plot_color ), hold on
                    else,                          plot( 10 .^ x, pdf,                 'Color', plot_color ), hold on
                                               H = area( 10 .^ x, pdf ); set(H(1), 'FaceColor', plot_color, 'EdgeColor', 'none', 'FaceAlpha', 0.1 );
                    end
                        
                    edges = x ;

                end
                
                XTickLabels = { };
                
                if histogram_style_index == 2
                
                    for i = 1 : numel( X_ticks )

                        XTickLabels{ end + 1 } = [ '10^{', num2str( X_ticks_log_10( i )), '}' ];

                    end
                end
                                
                set( gca, 'xscale',     'log',                                                 ...
                          'XTickLabel', XTickLabels, ...                    
                          'xlim',       [ min( 10 ^ edges(  1  ), X_ticks(  1  )), ...
                                          max( 10 ^ edges( end ), X_ticks( end ))], ...
                          'XTick',      X_ticks, ...
                          'XMinorTick','off')

        end % SWITCH field by which ones should have a log-distributed x axis
        
        set(gca,'fontname','times')

        switch histogram_styles{ histogram_style_index }
            
            case 'cumcount'
        
                ylim([ 0, numel( network_statistic )])
                
                if statistic_index == 1,                  ylabel( 'Number of Strands' ),             end

                if statistic_index ~= 1, set( gca, 'YTickLabel', [ ]), end                
                
            case 'count'
                
                if statistic_index == 1, if is_box_histo, ylabel( 'Number of Strands' ),                  ...
                                        else,             ylabel({'Number of Strands','/ Abscissa (dB) Unit'}),  end, end

                switch fields{ statistic_index }
                    
                    case 'strand_depths',          xlabel({    'Depth','[\mum]'       })
                    case 'strand_lengths',         xlabel({   'Length','[\mum]'       })
                    case 'strand_ave_radii',       xlabel({   'Radius','[\mum]'       })
                    case 'strand_areas',           xlabel({     'Area','[\mum^{2}]'   })
                    case 'strand_volumes',         xlabel({   'Volume','[\mum^{3}]'   })
                    case 'strand_tortuosity',      xlabel({        ' ','1 / Tortuosity' })
                    case 'strand_z_direction',     xlabel({'Alignment with Z-Axis','[Component]'})                        
                    case 'strand_ave_r_component'
                        if network_statistics.origin_dimensionality == 3
                                                   xlabel({'Orientation to Infarct Center','[Component]'      })
                        else,                      xlabel({'Orientation to Center Vessel' ,'[Component]'}) % origin_dimensionality == 2
                        end
                    case 'strand_ave_d'
                        if network_statistics.origin_dimensionality == 3
                                                   xlabel({'Distance from Infarct'      ,'[\mum]'})
                        else,                      xlabel({'Distance from Center Vessel','[\mum]'}) % origin_dimensionality == 2
                        end
                end % SWITCH field
        end % SWITCH histo style
        
        set( gca, 'TickDir', 'both' )
        
    end % FOR histo style 
end % FOR statistic

%     % Old function:
%     subplot( 2, 4, 1 )
% %         histogram( network_statistics.strand_lengths, 'Normalization', 'cdf' )
%     [ ~, edges ] = histcounts( log10( network_statistics.strand_lengths ));
%     histogram( network_statistics.strand_lengths, 10 .^ edges, 'Normalization', 'cdf' )
%     set( gca, 'xscale', 'log' )
%     ylim([ 0, 1 ])        
% 
%     subplot( 2, 4, 5 )
% %         histogram( network_statistics.strand_lengths, 'Normalization', 'pdf' )   
%     [ ~, edges ] = histcounts( log10( network_statistics.strand_lengths ));
%     histogram( network_statistics.strand_lengths, 10 .^ edges, 'Normalization', 'pdf' )
%     set( gca, 'xscale', 'log' )        
%     xlabel( 'strand length (um)' )
% 
%     subplot( 2, 4, 2 )
% %         histogram( network_statistics.strand_ave_radii, 'Normalization', 'cdf' )
%     [ ~, edges ] = histcounts( log10( network_statistics.strand_ave_radii ));
%     histogram( network_statistics.strand_ave_radii, 10 .^ edges, 'Normalization', 'cdf' )
%     set( gca, 'xscale', 'log' )
%     ylim([ 0, 1 ])
% 
%     subplot( 2, 4, 6 )
% %         histogram( network_statistics.strand_ave_radii, 'Normalization', 'pdf' )    
%     [ ~, edges ] = histcounts( log10( network_statistics.strand_ave_radii ));
%     histogram( network_statistics.strand_ave_radii, 10 .^ edges, 'Normalization', 'pdf' )
%     set( gca, 'xscale', 'log' )        
%     xlabel( 'strand radius (um)' )
% 
%     subplot( 2, 4, 3 )
% %         histogram( network_statistics.strand_areas, 'Normalization', 'cdf' )
%     [ ~, edges ] = histcounts( log10( network_statistics.strand_areas ));
%     histogram( network_statistics.strand_areas, 10 .^ edges, 'Normalization', 'cdf' )
%     set( gca, 'xscale', 'log' )        
%     ylim([ 0, 1 ])
% 
%     subplot( 2, 4, 7 )
% %         histogram( network_statistics.strand_areas, 'Normalization', 'pdf' )
%     [ ~, edges ] = histcounts( log10( network_statistics.strand_areas ));
%     histogram( network_statistics.strand_areas, 10 .^ edges, 'Normalization', 'pdf' )
%     set( gca, 'xscale', 'log' )        
%     xlabel( 'strand area (um^2)' )
% 
%     subplot( 2, 4, 4 )
% %         histogram( network_statistics.strand_volumes, 'Normalization', 'cdf' )
%     [ ~, edges ] = histcounts( log10( network_statistics.strand_volumes ));
%     histogram( network_statistics.strand_volumes, 10 .^ edges, 'Normalization', 'cdf' )
%     set( gca, 'xscale', 'log' )        
%     ylim([ 0, 1 ])
% 
%     subplot( 2, 4, 8 )
% %         histogram( network_statistics.strand_volumes, 'Normalization', 'pdf' )
%     [ ~, edges ] = histcounts( log10( network_statistics.strand_volumes ));
%     histogram( network_statistics.strand_volumes, 10 .^ edges, 'Normalization', 'pdf' )
%     set( gca, 'xscale', 'log' )        
%     xlabel( 'strand volume (um^3)' )

end % FUNCTION
