function [ network_statistics ] = calculate_network_statistics( strand_subscripts, bifurcation_vertices, lumen_radius_in_microns_range, microns_per_voxel, size_of_image )
%% calculate_network_statistics SAM 5/25/19
% This function calculates network statistics for the input vectorized network. Some statistics are
% network totals, some are densities (totals over total volume), while others are on a per strand
% basis. A strand is a vessel segment from branchpoint/endpoint to branchpoint/endpoint.  All
% statistics are in dimensions of microns or some power of microns.

%% unaveraged strand quantities

% length by strand
strand_delta_lengths    = cellfun( @( x ) sum(((   x( 2 : end    , 1 : 3 )                                           ...
                                                 - x( 1 : end - 1, 1 : 3 )) .* microns_per_voxel ) .^ 2, 2 ) .^ 0.5, ...
                                   strand_subscripts, 'UniformOutput', false                                         );
                            
strand_shortest_length  = cellfun( @( x ) sum(((   x(         end, 1 : 3 )                                           ...
                                                 - x(         1  , 1 : 3 )) .* microns_per_voxel ) .^ 2    ) .^ 0.5, ...
                                   strand_subscripts                                                                 );

strand_radii            = cellfun( @( x ) exp( interp1( log( lumen_radius_in_microns_range ), x( :, 4 ))), ...
                                   strand_subscripts, 'UniformOutput', false                               );

%% average strand quantities 

network_statistics.strand_lengths           = cellfun( @sum, strand_delta_lengths );

network_statistics.strand_tortuosity        = network_statistics.strand_lengths ./ strand_shortest_length ;

% diameter by strand (length weighted average of radius across each strand)
network_statistics.strand_ave_radii         = cellfun( @( delta_l, r ) sum( delta_l .* ( r( 2 : end ) + r( 1 : end - 1 )) / 2 ) / sum( delta_l ), strand_delta_lengths, strand_radii );

% direction by strand (length weighted average of direction across each strand)
% network_statistics.strand_ave_directions                                                                                                                                                                                ...
%                                   = cell2mat( cellfun( @( delta_l, v ) sum( delta_l .* ( v( 2 : end, : ) + v( 1 : end - 1, : )) / 2, 1 ) / sum( delta_l ), strand_delta_lengths, vessel_directions, 'UniformOutput', false ));

network_statistics.strand_ave_directions                                                                                                                                                                                ...
                                  = cell2mat( cellfun( @( x, l_s ) ( x( end, 1 : 3 ) - x( 1, 1 : 3 )) .* microns_per_voxel / l_s, strand_subscripts, num2cell( strand_shortest_length ), 'UniformOutput', false ));                              
                              
% surface area by strand (trapezoidal method integration of circumference (2*pi*r) with respect to length)
network_statistics.strand_areas             = cellfun( @( delta_l, r ) pi * sum( delta_l .* ( r( 2 : end )      + r( 1 : end - 1 ))),             strand_delta_lengths, strand_radii );

% volume by strand (trapezoidal method integration of cross-sectional area (pi*r^2) respect to length)
network_statistics.strand_volumes           = cellfun( @( delta_l, r ) pi * sum( delta_l .* ( r( 2 : end ) .^ 2 + r( 1 : end - 1 ) .^ 2 ) / 2 ),  strand_delta_lengths, strand_radii );

%% network total statistics

% total length
network_statistics.length           = sum( network_statistics.strand_lengths );

% total surface area
network_statistics.area             = sum( network_statistics.strand_areas   );

% total volume
network_statistics.volume           = sum( network_statistics.strand_volumes );

% total number of bifurcations
network_statistics.num_bifurcations = length( bifurcation_vertices );

%% density statistics

total_volume = prod( size_of_image .* microns_per_voxel );

network_statistics.length_density        = network_statistics.length           / total_volume ;

network_statistics.area_density          = network_statistics.area             / total_volume ;

network_statistics.volume_density        = network_statistics.volume           / total_volume ;

network_statistics.bifurcation_density   = network_statistics.num_bifurcations / total_volume ;

%% derived statistics

network_statistics.strand_z_direction = abs( network_statistics.strand_ave_directions( :, 3 ));

end % FUNCTION
