function [x, pdf, cdf ] = smooth_hist(stats, sd, varargin )
%% Annie Zhou, March 2022
% 
% inputs: 
%   stats: stats vector of interest
%   sd: standard deviation of Gaussians to be calculated
%
% outputs: % SAM 4/22/22 
%   X: a regular sampling spanning the values attained by variable STATS at a rate >= 1/SD
%   CDF: smoothed cdf of STATS. Smoothed using a kernel (Normal(STATS, SD)).
%   PDF: (symmetric) difference quotient of the CDF vs X plot

%% 

if isempty( varargin ), weights = ones( size( stats )); else, weights = varargin{1}; end

% % set x values for calculating Gauss curves to be the unique values in stats
% % x = unique(stats)'; 
% regular_sampling = min( stats ): sd : ...
%                    max( stats )       ;
% x = unique([stats;regular_sampling']); % SAM 3/18/22

num_samples = ceil((   max( stats ) ...
                      - min( stats )) / sd );

x = linspace( min( stats ), max( stats ), num_samples )';

% % alternate method for setting x values to be equally spaced
% res = length(unique(stats)); % "resolution" of each calculated Gaussian curve: number of points used to calculate curve
% x = linspace(min(stats), max(stats), res)';

% sum Gaussian curves
s_hist = zeros(size(x));
for i = 1:length(stats) % !!! try to only loop through length(x) for speedup % SAM 3/11/22
    mu =   stats(i);
    w  = weights(i);
    
    s_hist = s_hist + w * gauss(x, mu, sd);
end

% % plot(s)
% figure; plot(x,s_hist)
% % figure; histogram(stats) 

cdf = cumsum( s_hist );

bin_widths = [  x( 2       ) - x( 1            )     ; ...
              ( x( 3 : end ) - x( 1 : end - 2 )) / 2 ;
                x(     end ) - x(     end - 1 )        ];

pdf = s_hist ./ bin_widths ;

function [gd] = gauss(x, mu, sd) % Calculate Gaussian curve
%    gd = 1/(2*pi*sd)*exp(-(x-mu).^2/(2*sd^2));
    gd = exp(-(x-mu).^2/(2*sd^2));
    gd = gd / sum( gd ); % normalized to have total weighting of 1 % SAM 3/11/22
end
end