function [ KS_test_stat,    difference_distribution,    cumulative_difference_distribution, ...
           KS_test_stat_U,  difference_distribution_U,  cumulative_difference_distribution_U, x ] = weighted_KStest2(stats1, weights1, ...
                                                                                                                     stats2, weights2, resolution, extrap_option )
%% weighted_KStest2 - SAM 4/5/23

%% !!!!! remove outliers 
% ??? remove datapoints until the change is within the resolution (sd) to some factor

%% smooth the data at the requested resolution, downsample, and regular interpolate

[x1, pdf1, cdf1 ] = smooth_hist(stats1, resolution, weights1 );
[x2, pdf2, cdf2 ] = smooth_hist(stats2, resolution, weights2 );

% put PDF/CDF on common support which is the union of the two supports
x = unique([ x1 ; x2 ]);

% resample support to have approximate regular spacing at resolution
x = linspace( min(x), max(x), round((max(x)-min(x))/resolution));

pdf1 = interp1(x1,pdf1,x,"linear",extrap_option);
pdf2 = interp1(x2,pdf2,x,"linear",extrap_option);

cdf1 = interp1(x1,cdf1,x,"linear",NaN); cdf1(x<min(x1)) = 0 ; cdf1(x>max(x1)) = max(cdf1) ;
cdf2 = interp1(x2,cdf2,x,"linear",NaN); cdf2(x<min(x2)) = 0 ; cdf2(x>max(x2)) = max(cdf2) ;

%% retain unnormalized distributions for curiosity (method == 1)
for method_idx = 1 : 2
    
    if method_idx == 2
        %% normalize the distributions
        pdf1 = pdf1 / max(cdf1);
        pdf2 = pdf2 / max(cdf2);
        
        cdf1 = cdf1 / max(cdf1);
        cdf2 = cdf2 / max(cdf2);
    end
    
    %% compute difference distribution
    pdfD = pdf2 - pdf1 ;
    cdfD = cdf2 - cdf1 ;
    
    % %% integrate the difference distribution ??? complete ???? 
    
    %% find max abs value of integrated difference distribution            
    KSts = max(abs(cdfD));

    switch method_idx

        case 1

                       difference_distribution_U = pdfD ;
            cumulative_difference_distribution_U = cdfD ;
                                  KS_test_stat_U = KSts ;

        case 2

                       difference_distribution   = pdfD ;
            cumulative_difference_distribution   = cdfD ;
                                  KS_test_stat   = KSts ;

    end
end

end