function network_histogram_plotter( network_statistics, varargin ) %#ok<INUSL>

histogram_styles = { 'cdf', 'pdf' };

if ~ isempty( varargin )

    fields = varargin ;

else % ELSE use default

    fields = { 'strand_lengths', 'strand_ave_radii', 'strand_areas', 'strand_volumes' };

end % IF varargin is nonempty

number_of_statistics = length( fields );

statistic_indices = 1 : number_of_statistics ;

figure

for statistic_index = statistic_indices

    network_statistic = eval([ 'network_statistics.', fields{ statistic_index }, ';' ]);    
            
    for histogram_style_index = 1 : 2

        subplot( 2, number_of_statistics, ( histogram_style_index - 1 ) * number_of_statistics + statistic_index )
    
        switch fields{ statistic_index }

            case { 'strand_z_direction' } % statistics that should have linear scale on x axis

                histogram( network_statistic, 'Normalization', histogram_styles{ histogram_style_index })

            otherwise % use a log scale on x axis

                [ ~, edges ] = histcounts( log10( network_statistic ));
                histogram( network_statistic, 10 .^ edges, 'Normalization', histogram_styles{ histogram_style_index })
                set( gca, 'xscale', 'log' )

        end % SWITCH field by which ones should have a log-distributed x axis

        switch histogram_styles{ histogram_style_index }
            
            case 'cdf'
        
                ylim([ 0, 1 ])
                ylabel( 'CDF' )
                
            case 'pdf'
                
                ylabel( 'PDF' )

                switch fields{ statistic_index }
                    
                    case 'strand_lengths',      xlabel( 'strand length (um)' )
                    case 'strand_ave_radii',    xlabel( 'strand radius (um)' )
                    case 'strand_areas',        xlabel( 'strand area (um^2)' )
                    case 'strand_volumes',      xlabel( 'strand volume (um^3)' )
                    case 'strand_z_direction',  xlabel( 'strand direction z-component (component of unit vector)' )
                        
                end % SWITCH field
        end % SWITCH histo style
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

