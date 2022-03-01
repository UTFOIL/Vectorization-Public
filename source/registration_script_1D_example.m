%% 1D vector registration example

tic

% before and after images 1D
NumPoints = 10 ;

translation = 10 ;

epsilon = translation / 10 ; % translation units

domain = 100 ;

NumDeleted = round( 0.2 * NumPoints );

NumAdded = NumDeleted ;

PriorSigma = domain / 3 ;

PriorMean = 0 ; % initial guess for the translation parameter

noise = domain / NumPoints / 100 ;

before = domain * rand( 1, NumPoints );

after = [ before( 1 : end - NumDeleted ), domain * rand( 1, NumAdded )] ...
      + translation + noise * randn( 1, NumPoints );

NumDeltas = 100 ;

metric = zeros( 1, NumDeltas );

IsOptimizing = true ;

NumIterations = 0 ;

NumAfter  = numel( after  );
NumBefore = numel( before );

NumPairs = min( NumAfter, NumBefore );

after2before = spalloc(numel(after),1,NumPairs);

pairs = 1 : NumPairs ;

MetricsBefore2after = zeros( 1, NumBefore );
            
before2after        = zeros( 1, NumBefore );

while IsOptimizing

    NumIterations = NumIterations + 1 ;
    
    perturbation_idx = 0 ;
    
    ParamSet = PriorMean + PriorSigma * randn( 1, NumDeltas );

    for param = ParamSet

        perturbation_idx = perturbation_idx + 1 ;
        
        % gaussian of pointwise distances, a pairwise metric
        metric_matrix = exp( - ( before + param - after' ) .^ 2 / 2 / PriorSigma ^ 2 );
        
        % go from best pair to worst pair, forcing an invertible transform
        for pair = pairs
    
            [ PairMetric, idx ] = max( metric_matrix( : ));
            
            BeforeIdx = floor(( idx - 1 ) / NumAfter ) + 1 ;
            
            AfterIdx = idx - ( BeforeIdx - 1 ) * NumAfter ;
            
            MetricsBefore2after( BeforeIdx ) = PairMetric ;
            
            before2after( BeforeIdx ) = AfterIdx ;
            
            metric_matrix( AfterIdx,  : ) = 0 ;
            metric_matrix( :, BeforeIdx ) = 0 ;
            
        end
            
        metric( perturbation_idx ) = mean( MetricsBefore2after );
            
%         [ MetricsAfter2before, after2before ] = max( metric_matrix, [ ], 2 );
%         [ MetricsBefore2after, before2after ] = max( metric_matrix, [ ], 1 );
%         
%         % forwards and backwards directions pick closest neighbors as registration partners
%         metric( perturbation_idx ) = median( MetricsAfter2before ) ...
%                                    * median( MetricsBefore2after );
    
    end
    
    % first and second moments of the metric-weighted distribution of translations
    PosteriorMean  =   sum( metric .*   ParamSet                        ) / sum( metric )        ;
    PosteriorSigma = ( sum( metric .* ( ParamSet - PosteriorMean ) .^ 2 ) / sum( metric )) ^ 0.5 ;
        
%     MeanUpdate = PosteriorMean  - PriorMean ; % !!! dont let the post sigma fall below this value
    
    % is_optimizing if the mean isn't converging or the search is broadening
    IsOptimizing  = abs( PosteriorMean  - PriorMean  ) > epsilon ... 
                 ||    ( PosteriorSigma > PriorSigma )           ...
                 ||      PosteriorSigma                > epsilon ;
              
	% update transformation and metric formula ;
    PriorMean  = PosteriorMean  ;
    PriorSigma = PosteriorSigma ;
        
end

NumIterations
PriorMean
PriorSigma

% recalculate best registration
param = PriorMean ;

metric_matrix = exp( - ( before + param - after' ) .^ 2 / 2 / PriorSigma ^ 2 );

% go from best pair to worst pair, forcing an invertible transform
for pair = pairs

    [ PairMetric, idx ] = max( metric_matrix( : ));

    BeforeIdx = floor(( idx - 1 ) / NumAfter ) + 1 ;

    AfterIdx = idx - ( BeforeIdx - 1 ) * NumAfter ;

    MetricsBefore2after( BeforeIdx ) = PairMetric ;

    before2after( BeforeIdx ) = AfterIdx ;

    metric_matrix( AfterIdx,  : ) = 0 ;
    metric_matrix( :, BeforeIdx ) = 0 ;

end


figure
scatter( before,   ones( size( before )))
hold on
scatter(  after, 2*ones( size( after  )))
ylim([0,3])
set(gca,'YTickLabel',[ ]);
ylabel('Before                 After')

line( [ before ; after( before2after )], [ 1 ; 2 ] .* ones( 2, NumBefore ), 'Color','black','LineStyle','--' )


[ before, idx ] = sort( before );

plot( before, MetricsBefore2after(idx) )

% for BeforeIdx = 1 : NumBefore
% 
%     line(before2after( BeforeIdx ) = AfterIdx ;
% 
%     metric_matrix( AfterIdx,  : ) = 0 ;
%     metric_matrix( :, BeforeIdx ) = 0 ;
% 
% end



toc