function [ edge_metrics ] = get_edge_metric( edge_energies )
%% GET_EDGE_METRIC SAM + WAS 4/26/19
% calculates a quality metric (probability of existence) for the edge objects

% metric = 'mean' ;
metric = 'max' ; % 5/15/19 11:54
% metric = 'max_rel_min' ; % 11/8/19 1100
% metric = 'max_rel_max_vertex' ; % 11/11/19 1345 SAM

% metric = 'mean' ; % 5/16/19 1917

switch metric
    
    case 'mean'

        % % summarizing each trajectory by its mean energy 
        edge_metrics = cellfun( @mean,  edge_energies );
        
    case 'max'

        edge_metrics = cellfun( @max,  edge_energies );

    case 'max_rel_min'

        % % summarizing each trajectory by its relative change from min to max
        edge_metrics = cellfun( @( x ) 1000 * (( min( x ) - max( x )) / min( x ) - 1 ),  edge_energies );
        % mean_edge_energies = cellfun( @( x ) ( min( x ) - max( x )) / max( x ),  edge_energies_temp );
        
    case 'max_rel_max_vertex'
        
        edge_metrics = cellfun( @( x ) 1000 * ( - max( x ) / max( x( 1 ), x( end ))),  edge_energies );
        
        
% summarizing each trajectory by its relative standard deviation (a.k.a. coefficient of variation)
% mean_edge_energies    = cellfun( @( x ) mean( x ) / std( x ),  edge_energies );    

% % subtract off a linear fit to the energy vs edge ordinate.  This line is fully constrained by the energies of the terminal vertices.  SAM 2/21/19
% mean_edge_energies    = cellfun( @( x ) - std( x' - linspace( x( 1 ), x( end ), length( x ))) / mean( x ),  edge_energies );    


end % swtich

% for metrics in [ -1000, 0 ]
edge_metrics( isnan( edge_metrics )) = - 1000 ;

end % FUNCTION

