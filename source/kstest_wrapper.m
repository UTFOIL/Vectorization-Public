%% kstest wrapper
% function kstest_wrapper( file )
global datasetName is_comparing_stroke_intermouse

% datasetName = 'brett' ;
% datasetName = 'eddy' ;
% datasetName = 'intermouse' ;
% % datasetName = 'depth-resolved' ;
datasetName = 'pre-/post-stroke' ;

is_comparing_stroke_intermouse = false ;

file = { };
%         file  =  cell(1,numFiles);

switch datasetName 

    case 'brett'
                
%         file{ 1 } = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 1\truncated\batch_220701-164800\vectors\network_230308-123034_Brett_Fused_Raw_w1__r_cropped_decomposed_stats.mat' ;
        file{ 1 } = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 1\truncated\batch_230318-150419\vectors\network_230318-150419_Brett_Fused_Raw_w1__r_cropped_decomposed_stats.mat' ;
%         file{ 2 } = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 2\New folder\Truncated\batch_230109-123958\vectors\network_230201-060057_Brett_Fused_Raw_w2__r_cropped_decomposed_stats.mat' ;
        file{ 2 } = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 2\New folder\Truncated\batch_230109-123958\vectors\network_230316-221628_Brett_Fused_Raw_w2__r_cropped_decomposed_stats.mat' ;
        file{ 3 } = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 3\truncated\batch_221106-224630\vectors\network_230125-100328_Brett_Fused_Raw_w3__r_cropped_decomposed_stats.mat' ;
        file{ 4 } = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 4\batch_221018-164034\vectors\network_230222-024534_Brett_Fused_Raw_w4__r_cropped_decomposed_stats.mat' ;
%         file{ 5 } = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 5\batch_220614-230343\vectors\network_230308-150236_Brett_Fused_Raw_w5__r_cropped_decomposed_stats.mat' ;
        file{ 5 } = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 5\batch_230318-084151\vectors\network_230318-145228_Brett_Fused_Raw_w5__r_cropped_decomposed_stats.mat' ;
        
        name{ 1 } = 'week 0' ;
        name{ 2 } = 'week 1' ;
        name{ 3 } = 'week 2' ;
        name{ 4 } = 'week 3' ;
        name{ 5 } = 'week 4' ;

    case 'eddy'

        file{ 1 } = 'E:\2P imaging\2021_Chronic_Imaging\Eddy week 1\batch_210604-203819\vectors\network_220531-103816_Eddy_Fused_Raw_w1_log_stable__r_cropped_decomposed_stats.mat' ;
        file{ 2 } = 'E:\2P imaging\2021_Chronic_Imaging\Eddy week 3\batch_210616-194614\vectors\network_210617-184739_Eddy_Fused_Raw_w3__r_cropped_decomposed_stats.mat' ;
        file{ 3 } = 'E:\2P imaging\2021_Chronic_Imaging\Eddy week 5\batch_210618-130838\vectors\network_221120-211737_Eddy_Fused_Raw_w5__r_cropped_decomposed_stats.mat' ;

        name{ 1 } = 'week 0' ;
        name{ 2 } = 'week 2' ;
        name{ 3 } = 'week 4' ;
        
    case 'intermouse'

%         file{ 1 } = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 1\truncated\batch_220701-164800\vectors\network_230308-123034_Brett_Fused_Raw_w1__r_cropped_decomposed_stats.mat' ;
        file{ 1 } = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 1\truncated\batch_230318-150419\vectors\network_230318-150419_Brett_Fused_Raw_w1__r_cropped_decomposed_stats.mat' ;
        file{ 2 } = 'E:\2P imaging\2021_Chronic_Imaging\Eddy week 1\batch_210604-203819\vectors\network_220531-103816_Eddy_Fused_Raw_w1_log_stable__r_cropped_decomposed_stats.mat' ;
        file{ 3 } = 'E:\2P imaging\2021_Chronic_Imaging\Eddy week 3\batch_210616-194614\vectors\network_210617-184739_Eddy_Fused_Raw_w3__r_cropped_decomposed_stats.mat' ;
        file{ 4 } = 'E:\2P imaging\2021_Chronic_Imaging\Eddy week 5\batch_210618-130838\vectors\network_221120-211737_Eddy_Fused_Raw_w5__r_cropped_decomposed_stats.mat' ;

%         name{ 1 } = 'mouse A' ;
%         name{ 2 } = 'mouse B' ;
%         name{ 3 } = 'mouse C' ;

        name{ 1 } = 'stroke week 0' ;
        name{ 2 } = 'control week 0' ;
        name{ 3 } = 'control week 2' ;
        name{ 4 } = 'control week 4' ;

    case 'depth-resolved'

        file{ 1 } = '' ;
        file{ 2 } = '' ;
        file{ 3 } = '' ;
        
        name{ 1 } = 'stroke week 0' ;
        name{ 1 } = 'stroke week 0' ;
        name{ 1 } = 'stroke week 0' ;

    case  'pre-/post-stroke'

        file{ 1 } = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 1\truncated\batch_230318-150419\vectors\network_230318-150419_Brett_Fused_Raw_w1__r_cropped_decomposed_stats.mat' ;
        file{ 2 } = 'E:\2P imaging\2021_Chronic_Imaging\Brett week 5\batch_230318-084151\vectors\network_230318-145228_Brett_Fused_Raw_w5__r_cropped_decomposed_stats.mat' ;
        file{ 3 } = 'E:\2P imaging\2021_Chronic_Imaging\Dwight week 1\batch_230521-012504\vectors\network_230521-131318_Dwight_Fused_Raw_w1__r_cropped_decomposed_stats.mat' ;
        file{ 4 } = 'E:\2P imaging\2021_Chronic_Imaging\Dwight week 5\batch_230520-113855\vectors\network_230520-234351_Dwight_Fused_Raw_w5__r_cropped_decomposed_stats.mat' ;
        file{ 5 } = 'E:\2P imaging\2021_Chronic_Imaging\Ant week 1\batch_230523-211303\vectors\network_230525-132949_Ant_Fused_Raw_w1__r_cropped_decomposed_stats.mat' ;
        file{ 6 } = 'E:\2P imaging\2021_Chronic_Imaging\Ant week 5\New folder\batch_230203-030500\vectors\network_230522-121519_Ant_Fused_Raw_w5__r_cropped_decomposed_stats.mat' ;
        
        name{ 1 } = 'Mouse A week 0' ;
        name{ 2 } = 'Mouse A week 4' ;
        name{ 3 } = 'Mouse C week 0' ;
        name{ 4 } = 'Mouse C week 4' ;
        name{ 5 } = 'Mouse D week 0' ;
        name{ 6 } = 'Mouse D week 4' ;

        file = file([1,3,5,2,4,6]);
        name = name([1,3,5,2,4,6]);

        is_comparing_stroke_intermouse = false ;

        if is_comparing_stroke_intermouse
            % % first set
%             file = file([1,3,5]);
%             name = name([1,3,5]);
            % % second set
            file = file([2,4,6]);
            name = name([2,4,6]);

            % run first set, then save variables, 
            % 
            % 
% xdT = {xdT,xd}
% xoT = {xoT,xo}
% cdfdDT = {cdfdDT,cdfdD}
% cdfoDT = {cdfoDT,cdfoD}            
            % 
            % run second set, append variables
% xdT = [xdT,xd]
% xoT = [xoT,xo]
% cdfdDT = [cdfdDT,cdfdD]
% cdfoDT = [cdfoDT,cdfoD]

% xd = xdT
% cdfdD = cdfdDT
% xo = xoT
% cdfoD = cdfoDT
            % comment out everything in this script but the last lines needed to generate the figure and plot,
            % highjacking the normal notation. Set 
% numFiles = 6

        end
end

%% extract data from files
numFiles = length( file );

    distance_cell = cell(numFiles,1);
%      radius_cell = cell(1,numFiles); % !!! radius was misnamed in *_decomposed_stats.mat
 orientation_cell = cell(numFiles,1);
delta_length_cell = cell(numFiles,1);

  x_cell = cell(numFiles,1);
pdf_cell = cell(numFiles,1);
cdf_cell = cell(numFiles,1);

length_array = zeros( numFiles, 1 );

NumStrands_array = zeros( numFiles, 1 );

idx_range = 1 : numFiles ;

is_strand_analysis   = false ;
is_strand_analysis_B = false ;
is_strand_analysis_C = false ;

% % for length-based calcs
    KStest_resolution_distance    = 0.05; % mm
    KStest_resolution_orientation = 0.05 ; % component

if is_strand_analysis
   
% % for strands-based calcs
    KStest_resolution_distance    = 0.05 ; % mm
    KStest_resolution_orientation = 0.1 ; % component

    if is_strand_analysis_B
    % % for B strands-based calcs
    KStest_resolution_distance    = 0.1 ; % length [decibel of um]
%     KStest_resolution_orientation = 0.125 ; % inverse tort
    KStest_resolution_orientation = 1/15 ; % asec * 2 / pi
% %     KStest_resolution_orientation = 0.1 ; % atan * 2 / pi

    end
end

for idx = idx_range

    load( file{ idx }(1:end-21), 'network_statistics')

    load( file{ idx })

%     % histograms
%     figure( 3 * idx     ), histogram(    distance )
% %     figure( 3 * idx + 1 ), histogram(      radius )
%     figure( 3 * idx + 2 ), histogram( orientation )
    
%         distance_cell{ idx } =   [ distance / 1000; 1; 1 ]; % [mm]
        distance_cell{ idx } =   [ distance / 1000; 1; 1 ] .^ 2 ; % [mm2]
%          radius_cell{ idx } =      radius ;
     orientation_cell{ idx } = [ orientation; 0; 1 ];

    delta_length_cell{ idx } = [ delta_length; 0; 0 ];

    length_array( idx ) = sum( delta_length );

    NumStrands_array( idx ) = NumStrands ;

    if is_strand_analysis

%            distance_cell{ idx } = [           network_statistics.strand_ave_d      / 1000 ;1;1]; % [mm2]
           distance_cell{ idx } = [           network_statistics.strand_ave_d      / 1000 ;1;1]; % [mm]
        orientation_cell{ idx } = [           network_statistics.strand_ave_r_component   ;0;1];
       delta_length_cell{ idx } = [ones( size(network_statistics.strand_ave_r_component ))/1000;0;0];

        if is_strand_analysis_B
%             angle_samples = 0 : 0.1 : 90 ;
%             x_coordinate = tand( angle_samples );
%             parabola_arc_length      =  x_coordinate.*(1+4*x_coordinate.^2).^0.5...
%                                  +log(2*x_coordinate +(1+4*x_coordinate.^2).^0.5)/2 ;
%             tortuosity_LUT_by_angle = parabola_arc_length / 2 ./ x_coordinate ;
% 
%             tortuosity_LUT_by_angle(1)=0 ;
% 
%             tortuosity_LUT_by_angle(end)=tortuosity_LUT_by_angle(end-1)*2 ;
% 
%             figure, plot( tortuosity_LUT_by_angle,angle_samples), hold on
%                         plot(secd(angle_samples),angle_samples)

           distance_cell{ idx } = [          log( network_statistics.strand_lengths ) / log( 10 );1;1]; % [decibel of um]
%         orientation_cell{ idx } = [     min( 1 ./ network_statistics.strand_tortuosity, 1 )      ;0;1]; % (inverse) unitless
         orientation_cell{ idx } = [     acsc(max(network_statistics.strand_tortuosity,1))*2/pi       ;0;1]; % arc-secant [degrees].,1))*2/pi       ;0;1]; % arc-secant [degrees]
%         orientation_cell{ idx } = [     max(log(network_statistics.strand_tortuosity),0).^0.25;0;0]; % log()^0.25
%         orientation_cell{ idx } = [
%         interp1(tortuosity_LUT_by_angle,angle_samples,max(network_statistics.strand_tortuosity,1))/90 ;0;1]; % arc-secant [degrees].,1))*2/pi       ;0;1]; % this one is the angle of a parabola with the given tortuosity, which is in close agreement with the acsc() calculation

            if is_strand_analysis_C

           distance_cell{ idx } =            network_statistics.strand_ave_radii ; % [um]

            end
        end
    end

    x_cell{   idx } = x   ;
    pdf_cell{ idx } = pdf ;
    cdf_cell{ idx } = cdf ;

end

%% total length/NumStrands time-course
figure
plot( length_array     /     length_array( 1 ) * 100, 'k-x' ) 
ylabel('Total Length [% of Initial]')
% xlabel('Time')
xticks(1:numFiles)
xticklabels(name)
display('Total Lengths: '), length_array
display(['Length rel std dev: ', num2str(std( length_array ) / mean( length_array ))])
hold on
plot( [1,5],100 - 100 * 0.0805 * [1, 1],'k--')

figure
plot( NumStrands_array / NumStrands_array( 1 ) * 100, 'k-x' )
ylabel('Total Num. Strands [% of Initial]')
xticks(1:numFiles)
xticklabels(name)
display('Total Num. Strands: '), NumStrands_array
display(['Num. Strands rel std dev: ', num2str(std( NumStrands_array ) / mean( NumStrands_array ))])
% hold on
% plot( [1,5],100 - 100 * 1 * [1, 1],'k--')

%% box plots

% !!!!! box plots are unweighted !!!!

idx_range_cell = mat2cell( idx_range', ones( size(idx_range)), 1 );

x1 = cell2mat(    distance_cell );
x2 = cell2mat( orientation_cell );

g_cell = cellfun( @( x ) repmat(name(x),numel(distance_cell{x}),1), idx_range_cell, 'UniformOutput', false );

g = { };

for idx = idx_range

    g = [ g ; g_cell{ idx }];

end

figure
boxplot(x1, g, 'Symbol','' ) % make outliers invisible
% title('Distance from Infarct\nCenter [mm]')
ylabel('Distance [mm]')
if is_strand_analysis_B, ylabel('Length [mm]'), end
set(gcf, 'Position',  [1136         912         269         189])

figure
boxplot(x2, g, 'Symbol','' ) % make outliers invisible
% title('Orientation to Infarct Center')
ylabel('Orientation [Unitless]')
if is_strand_analysis_B, ylabel('Tortuosity [Degrees]'), end
set(gcf, 'Position',  [1136         912         269         189])


kd = [];
ko = [];

kdU = [];
koU = [];

xd = {};
xo = {};

xd1 = {};
xo1 = {};

pdfd = {};
pdfo = {};

pdfdU = {};
pdfoU = {};

cdfdU = {};
cdfoU = {};

pdfdD  = {};
cdfdD  = {};
pdfdDU = {};
cdfdDU = {};

pdfoD  = {};
cdfoD  = {};
pdfoDU = {};
cdfoDU = {};

%% KS test B: all compared to inital
for idx = idx_range

%     [ ~, pd(idx), kd(idx)]=kstest2(    distance_cell{ 1 },    distance_cell{ idx });
%     [ ~, po(idx), ko(idx)]=kstest2( orientation_cell{ 1 }, orientation_cell{ idx });

        [ kd( idx), pdfdD{ idx}, cdfdD{ idx}, ...       
          kdU(idx), pdfdDU{idx}, cdfdDU{idx}, xd{idx}] = weighted_KStest2(    distance_cell{  1  }, delta_length_cell{  1  }, ...
                                                                              distance_cell{ idx }, delta_length_cell{ idx },       KStest_resolution_distance,        0    );
        [ ko( idx), pdfoD{ idx}, cdfoD{ idx}, ...
          koU(idx), pdfoDU{idx}, cdfoDU{idx}, xo{idx}] = weighted_KStest2( orientation_cell{  1   }, delta_length_cell{  1   }, ...
                                                                           orientation_cell{ idx  }, delta_length_cell{ idx  }, KStest_resolution_orientation, "extrap" );

    [ xd1{idx}, pdfdU{idx}, cdfdU{idx}] = smooth_hist(    distance_cell{ idx }, KStest_resolution_distance   , delta_length_cell{ idx });
    [ xo1{idx}, pdfoU{idx}, cdfoU{idx}] = smooth_hist( orientation_cell{ idx }, KStest_resolution_orientation, delta_length_cell{ idx });

    pdfd{idx} = pdfdU{idx} / max(cdfdU{idx});
    pdfo{idx} = pdfoU{idx} / max(cdfoU{idx});

    cdfd{idx} = cdfdU{idx} / max(cdfdU{idx});
    cdfo{idx} = cdfoU{idx} / max(cdfoU{idx});

%     disp([ 'compare ', name{ idx  }, ...
%                ' to ', name{  1   }, '; distance p value: ', num2str( pd( idx )), ...
%                                   ', orientation p value: ', num2str( po( idx ))])

% end
% 
% for idx = idx_range

    disp([ 'compare ', name{ idx  }, ...
               ' to ', name{  1   }, '; distance KS stat [%]: ', num2str( 100 * kd( idx )), ...
                                  ', orientation KS stat [%]: ', num2str( 100 * ko( idx ))])

end

baseline     = mean( 100 * kd( 2 : end ));
baseline_std =  std( 100 * kd( 2 : end ));
display(['KS stat    distances: baseline: ', num2str( baseline ), '   St.Dev.: ', num2str(baseline_std)])

baseline     = mean( 100 * ko( 2 : end ));
baseline_std =  std( 100 * ko( 2 : end ));
display(['KS stat orientations: baseline: ', num2str( baseline ), '   St.Dev.: ', num2str(baseline_std)])

hd = figure;
plot( 1 : numFiles - 1, 100 * kd( 2 : end ), '-x', 'Color', "#D95319" )
hold on

ho = figure;
plot( 1 : numFiles - 1, 100 * ko( 2 : end ), '-x', 'Color', "#D95319" )
hold on

% basis density
figure, for idx = 1 : numFiles,     distribution_plotter(xd1{idx},    pdfdU{idx},idx-1, true), end, set(gca,'fontname','times'),  set(gcf, 'Position',  [879   601   122   109])
figure, for idx = 1 : numFiles,     distribution_plotter(xo1{idx},    pdfoU{idx},idx-1, true), end, set(gca,'fontname','times'),  set(gcf, 'Position',  [879   601   122   109])
% PDF
figure, for idx = 1 : numFiles,     distribution_plotter(xd1{idx},100*pdfd{ idx},idx-1, true), end, set(gca,'fontname','times'),  set(gcf, 'Position',  [879   601   122   109])
figure, for idx = 1 : numFiles,     distribution_plotter(xo1{idx},100*pdfo{ idx},idx-1, true), end, set(gca,'fontname','times'),  set(gcf, 'Position',  [879   601   122   109])
% Cumulative basis
figure, for idx = 1 : numFiles,     distribution_plotter(xd1{idx},    cdfdU{idx},idx-1,false), end, set(gca,'fontname','times'),  set(gcf, 'Position',  [879   601   122   109])
figure, for idx = 1 : numFiles,     distribution_plotter(xo1{idx},    cdfoU{idx},idx-1,false), end, set(gca,'fontname','times'),  set(gcf, 'Position',  [879   601   122   109])
% CDF
figure, for idx = 1 : numFiles,     distribution_plotter(xd1{idx},100*cdfd{ idx},idx-1,false), end, set(gca,'fontname','times'),  set(gcf, 'Position',  [879   601   122   109])
figure, for idx = 1 : numFiles,     distribution_plotter(xo1{idx},100*cdfo{ idx},idx-1,false), end, set(gca,'fontname','times'),  set(gcf, 'Position',  [879   601   122   109])
% CDF difference
figure, for idx = 1 : numFiles, distribution_plotter(xd{ idx},100*cdfdD{idx},idx-1,false), end, set(gca,'fontname','times'),  set(gcf, 'Position',  [879   601   122   109])
figure, for idx = 1 : numFiles, distribution_plotter(xo{ idx},100*cdfoD{idx},idx-1,false), end, set(gca,'fontname','times'),  set(gcf, 'Position',  [879   601   122   109])

%% KS test A: Current Compared to previous

for idx = idx_range

    idx2 = mod( idx - 2, numFiles ) + 1 ;

%     [ ~, pd(idx), kd(idx)]=kstest2(    distance_cell{ idx2 },    ...
%                                        distance_cell{ idx  });
%     [ ~, po(idx), ko(idx)]=kstest2( orientation_cell{ idx2 },     ...
%                                     orientation_cell{ idx  });

    [ kd( idx), pdfdD{ idx}, cdfdD{ idx}, ...       
      kdU(idx), pdfdDU{idx}, cdfdDU{idx}, xd{idx} ] = weighted_KStest2(    distance_cell{ idx2 }, delta_length_cell{ idx2 }, ...
                                                                           distance_cell{ idx  }, delta_length_cell{ idx  }, KStest_resolution_distance,        0    );
    [ ko( idx), pdfoD{ idx}, cdfoD{ idx}, ...
      koU(idx), pdfoDU{idx}, cdfoDU{idx}, xo{idx} ] = weighted_KStest2( orientation_cell{ idx2 }, delta_length_cell{ idx2 }, ...
                                                                        orientation_cell{ idx  }, delta_length_cell{ idx  }, KStest_resolution_orientation, "extrap" );


%                     [ x, pdf, cdf ] = smooth_hist( log( network_statistic ) / log( 10 ), ( log( max_stat ) - log( min_stat )) / log( 10 ) / number_of_bins );
%                     
%                     if histogram_style_index == 1, plot( 10 .^ x, cdf,                 'Color', plot_color ), hold on
%                     else,                          plot( 10 .^ x, pdf,                 'Color', plot_color ), hold on
%                                                H = area( 10 .^ x, pdf ); set(H(1), 'FaceColor', plot_color, 'EdgeColor', 'none', 'FaceAlpha', 0.1 );
%                     end
%                         
%                     edges = x ;


    
%     disp([ 'compare ', name{ idx2 }, ...
%               ' and ', name{ idx  }, '; distance p value: ', num2str( pd( idx )), ...
%                                   ', orientation p value: ', num2str( po( idx ))])

%     idx2 = mod( idx - 2, numFiles ) + 1 ;

    disp([ 'compare ', name{ idx2 }, ...
              ' and ', name{ idx  }, '; distance KS stat [%]: ', num2str( 100 * kd( idx )), ...
                                  ', orientation KS stat [%]: ', num2str( 100 * ko( idx ))])

end


plot( hd.CurrentAxes, 1 : numFiles - 1, 100 * kd( 1 : end - 1), 'x-', 'Color', "#0072BD" )
plot( ho.CurrentAxes, 1 : numFiles - 1, 100 * ko( 1 : end - 1), 'x-', 'Color', "#0072BD" )


% % !!!! also plot pdf's not (differenced) here while we are at it, because it is easier than using the vectorize() fxn on repeat
% figure
% for idx = 1 : numFiles
% 
%     [ x, pdf, cdf ] = smooth_hist( distance_cell{ idx  }, KStest_resolution_distance, delta_length_cell{ idx  });
% 
% 
% end

% comment out everything up to here to highjack notation and make intermouse comparison at both pre-
% and post- time-points (numFiles == 6 not 3)

% CDF difference
figure, for idx = 1 : numFiles, distribution_plotter(xd{idx},100*cdfdD{idx},idx-1,false), end, set(gca,'fontname','times'), set(gcf, 'Position',  [879   601   122   109])
figure, for idx = 1 : numFiles, distribution_plotter(xo{idx},100*cdfoD{idx},idx-1,false), end, set(gca,'fontname','times'), set(gcf, 'Position',  [879   601   122   109])

baseline     = mean( 100 * kd( 1 : end ));
baseline_std =  std( 100 * kd( 1 : end ));
display(['KS stat    distances: baseline: ', num2str( baseline ), '   StdDev: ', num2str(baseline_std)])

baseline     = mean( 100 * ko( 1 : end ));
baseline_std =  std( 100 * ko( 1 : end ));
display(['KS stat orientations: baseline: ', num2str( baseline ), '   StdDev: ', num2str(baseline_std)])

%% KS test Baseline (!!!! Hardcoded... to softcode: fxn'alize the kstest wrapper and call it with the different datset inputs here to extract baseline values
is_plotting_baseline = false ;

if is_plotting_baseline

%     plot(hd.CurrentAxes,[0,4],2.7*[1,1]-1.2,'k--')
    plot(hd.CurrentAxes,[0,4],2.26*[1,1]+1.26,'k--')
%     H1 = area( hd.CurrentAxes,[0,4], 2.9*[1,1]+1.1, 2.9-1.1 ); set(H1, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.1 );    

    plot(hd.CurrentAxes,[0,2,4],[2.97,3.01,0.81],'ko')

%     plot(hd.CurrentAxes,[0,4],5.4*[1,1]-1.8,'m--')
    plot(hd.CurrentAxes,[0,4],5.23*[1,1]+1.60,'m--')
%     H2 = area( hd.CurrentAxes,[0,4], 5.5*[1,1]+1.7, 5.5-1.7 ); set(H2, 'FaceColor', 'm', 'EdgeColor', 'none', 'FaceAlpha', 0.1 );    

%     plot(hd.CurrentAxes,[0,2,4],[3.5,6.3,6.8],'mo')
    plot(hd.CurrentAxes,[0,2,4],[3.37,6.16,6.15],'mo')

%     rectangle('Position',[0,5.5-1.7,4,1.7*2],'FaceColor','m')
    
end
if is_plotting_baseline
%     plot(ho.CurrentAxes,[0,4],0.95*[1,1]-0.11,'k--')
    plot(ho.CurrentAxes,[0,4],0.81*[1,1]+0.43,'k--')
%     H3 = area( ho.CurrentAxes,[0,4], 0.88*[1,1]+0.45, 0.88-0.45 ); set(H3, 'FaceColor', 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.1 );    

    plot(ho.CurrentAxes,[0,2,4],[0.33,1.14,0.96],'ko')

%     plot(ho.CurrentAxes,[0,4],2.42*[1,1]-0.55,'m--')
    plot(ho.CurrentAxes,[0,4],1.67*[1,1]+0.55,'m--')
%     H4 = area( ho.CurrentAxes,[0,4], 2.01*[1,1]+0.59, 2.01-0.59 ); set(H4, 'FaceColor', 'm', 'EdgeColor', 'none', 'FaceAlpha', 0.1 );    

%     plot(ho.CurrentAxes,[0,2,4],[1.84,2.66,1.51],'mo')
    plot(ho.CurrentAxes,[0,2,4],[1.54,2.28,1.21],'mo')
end

set(hd.CurrentAxes,'XLim',[0,4])
set(ho.CurrentAxes,'XLim',[0,4])

y_limits = get(hd.CurrentAxes,'YLim');
           set(hd.CurrentAxes,'YLim',[0,y_limits(2)])
y_limits = get(ho.CurrentAxes,'YLim');
           set(ho.CurrentAxes,'YLim',[0,y_limits(2)])

set(hd, 'Position', [1002         584         184          98])
set(ho, 'Position', [1002         584         184          98])

function distribution_plotter(x,pdf,idx,area_flag)

    global datasetName is_comparing_stroke_intermouse
    
    switch idx 

        case 0, plot_color = [ 0, 0, 0 ];
        case 1, plot_color = [ 1, 0, 0 ];
        case 2, plot_color = [ 0, 1, 0 ];
        case 3, plot_color = [ 0, 0, 1 ];
        case 4, plot_color = [ 1, 1, 1 ] * 2/3 ;
    end  
    if strcmp(datasetName,'brett')

    switch idx 

        case 0, plot_color = [  0,   0 ,  0  ];
        case 1, plot_color = [  1,   0 ,  0  ];
        case 2, plot_color = [ 4/5, 2/9, 2/9 ];
        case 3, plot_color = [ 3/4, 4/9, 4/9 ];
        case 4, plot_color = [  1,   1 ,  1  ] * 2/3 ;
    end  

    end

    if strcmp(datasetName,'eddy')

        switch idx 
    
            case 0, plot_color = [ 0, 0.25, 0.75 ];
            case 1, plot_color = [ 0, 0.5,  0.5  ];
            case 2, plot_color = [ 0, 0.75, 0.25 ];

        end  
    end

    if strcmp(datasetName,'intermouse')

        switch idx 
    
            case 0, plot_color = [ 0, 0,    0    ];
            case 1, plot_color = [ 0, 0.25, 0.75 ];
            case 2, plot_color = [ 0, 0.5 , 0.5  ];
            case 3, plot_color = [ 0, 0.75, 0.25 ];

        end  
    end
         
    if strcmp(datasetName,'pre-/post-stroke')

% %         if is_comparing_intermouse, idx = idx*2 ; end
%         if is_comparing_intermouse, idx = idx*2 + 1; end


% % %             case 0, plot_color = [ 0, 0, 0 ];
% % %             case 1, plot_color = [ 1, 1, 1 ] * 2/3 ;
% % %             case 2, plot_color = [ 1, 0, 0 ] * 2/3 ;
% % %             case 3, plot_color = [ 0, 1, 1 ] * 2/3 ;
% % %             case 4, plot_color = [ 0, 0, 1 ] * 2/3 ;
% % %             case 5, plot_color = [ 1, 1, 0 ] * 2/3 ;

% %             case 0, plot_color = [ 0, 0, 0 ];
% %             case 1, plot_color = [ 1, 0, 0 ] * 2/3 ;
% %             case 2, plot_color = [ 0, 0, 1 ] * 2/3 ;
% %             case 3, plot_color = [ 1, 1, 1 ] * 2/3 ;
% %             case 4, plot_color = [ 0, 1, 1 ] * 2/3 ;                
% %             case 5, plot_color = [ 1, 1, 0 ] * 2/3 ;

%             case {0,1,2}, plot_color = [ 0, 0, 0 ];
%             case {3,4,5}, plot_color = [ 1, 1, 1 ] * 2/3 ;


% % % % % % % % % % % % % %             case 0, plot_color = [ 0, 0, 0 ];
% % % % % % % % % % % % % %             case 1, plot_color = [ 1, 0, 0 ] * 1/3 ;
% % % % % % % % % % % % % %             case 2, plot_color = [ 0, 0, 1 ] * 1/3 ;
% % % % % % % % % % % % % %             case 3, plot_color = [ 1, 1, 1 ] * 2/3 ;
% % % % % % % % % % % % % %             case 4, plot_color = [ 0.5, 1, 1 ] * 2/3 ;                
% % % % % % % % % % % % % %             case 5, plot_color = [ 1, 1, 0.5 ] * 2/3 ;


        if ~is_comparing_stroke_intermouse

            switch idx 
    
                case 0, plot_color = [ 0,  0 ,  0  ]       ;
                case 2, plot_color = [ 0, 1/8, 1/4 ] * 2/3 ;
                case 4, plot_color = [ 0, 1/4, 1/2 ] * 2/3 ;
                case 1, plot_color = [ 1,  1 ,  1  ] * 2/3 ;
                case 3, plot_color = [ 1, 7/8, 3/4 ] * 2/3 ;
                case 5, plot_color = [ 1, 3/4, 1/2 ] * 2/3 ;    
            end


        else

            switch idx

                case 0, plot_color = [ 0,  0 ,  0  ]       ;
                case 1, plot_color = [ 0, 1/8, 1/4 ] * 2/3 ;
                case 2, plot_color = [ 0, 1/4, 1/2 ] * 2/3 ;
                case 3, plot_color = [ 1,  1 ,  1  ] * 2/3 ;
                case 4, plot_color = [ 1, 7/8, 3/4 ] * 2/3 ;
                case 5, plot_color = [ 1, 3/4, 1/2 ] * 2/3 ;                
            end
        end
    end            
    
                      plot( x, pdf,                'Color', plot_color ), hold on
    if area_flag, H = area( x, pdf ); set(H(1),'FaceColor', plot_color, 'EdgeColor', 'none', 'FaceAlpha', 0.1 ); end

end