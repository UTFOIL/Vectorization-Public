function subsampled_strand_subscripts = subsample_vectors( strand_subscripts, samples_per_strand )

% expecting strand subscripts in units of microns for positions and microns of radius for size

% vessel directions may or may not be appended as the last three columns

strand_delta_lengths    = cellfun( @( x ) sum(((   x( 2 : end    , 1 : 3 )                                           ...
                                                 - x( 1 : end - 1, 1 : 3 )) ) .^ 2, 2 ) .^ 0.5, ...
                                   strand_subscripts, 'UniformOutput', false                                         );


% !!! exerpted from resample_vectors

% using the L-2 norm to get even spacing in real-space
edge_cumulative_lengths                                                                                ...
                = cellfun( @( delta_l ) [ 0; cumsum( delta_l )], ...
                           strand_delta_lengths, 'UniformOutput', false                                     );
                       
edge_sample_lengths = cellfun( @( v ) linspace( 0, v( end ), samples_per_strand )',  ...
                               edge_cumulative_lengths, 'UniformOutput', false         );
                               
% interpolate the vector positions and size indices along the edge index, round to closest
% subscripts, then remove redundant vectors with unique function
number_of_edges = length( strand_subscripts );

for edge_index = 1 : number_of_edges
        
    strand_subscripts{ edge_index }       = interp1( edge_cumulative_lengths{ edge_index }, ...
                                                           strand_subscripts{ edge_index }, ...
                                                         edge_sample_lengths{ edge_index }  );
                                                     
	% !!!! renormalize vessel_directions ( columns 5 - 7 ) to be unit vectors


end % FOR

% % !!!!! edit the resample function to handle the L2 norm instead of just L-infinity, also handle
% % variable spacing to get N samples per strand instead of a sample every d microns.
% resample_vectors( [], [], [], use_case, samples_per_strand );

% decompose to list form (forget strand structure)
strand_subscripts = cell2mat( strand_subscripts );

end %FUNCTION subsample