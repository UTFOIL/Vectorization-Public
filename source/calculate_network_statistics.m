function [ network_statistics ] = calculate_network_statistics( strand_subscripts, bifurcation_vertices, lumen_radius_in_microns_range, microns_per_voxel, size_of_image, origin_location )
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

% strand depth (z component)
network_statistics.strand_depths            = cellfun( @( x ) mean( x( :, 3 )) * microns_per_voxel( 3 ), strand_subscripts );                              

if exist( 'origin_location', 'var' ), if ~ isempty( origin_location ) %#ok<ALIGN> 
    
%     origin_location = varargin{ 1 };

    if abs( origin_location( 3 )) == inf, dimensionality = 2 ;
    else,                                 dimensionality = 3 ;
    end

    strand_r = cellfun( @( x ) (                 x( :, 1 :  dimensionality ) ...
                                 - origin_location( :, 1 :  dimensionality ))...
                             .*  microns_per_voxel( :, 1 :  dimensionality ) , strand_subscripts, 'UniformOutput', false );

    network_statistics.strand_r = strand_r ;

	strand_d = cellfun( @( x ) sum( x .^ 2, 2 ) .^ 0.5 , strand_r, 'UniformOutput', false );

	network_statistics.strand_ave_d          = cellfun( @mean, strand_d );
    
%     % !!!!! this line appears incorrect. Take the mean of the projection, not the projection of the mean !!!!!! SAM 10/27/22
%     network_statistics.strand_ave_r_component = cellfun( @( r, d, v ) abs( sum( mean( r ./ d, 1 ) .* v )), strand_r, ...
%                                                                                                            strand_d, ...
%                                                          mat2cell(                      network_statistics.strand_ave_directions( :, 1 : dimensionality ),           ...
%                                                                             ones( size( network_statistics.strand_ave_directions( :, 1 : dimensionality ), 1 ), 1 ), ...
%                                                                                                                                          dimensionality              ));

    network_statistics.strand_ave_r_component = cellfun( @( r, d, v ) mean( abs( sum( r ./ d .* v, 2 ))),  strand_r, ...
                                                                                                           strand_d, ...
                                                         mat2cell(                      network_statistics.strand_ave_directions( :, 1 : dimensionality ),           ...
                                                                            ones( size( network_statistics.strand_ave_directions( :, 1 : dimensionality ), 1 ), 1 ), ...
                                                                                                                                         dimensionality              ));

    network_statistics.origin_dimensionality = dimensionality ;
    
end, end

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

network_statistics.strand_depths      = cellfun( @( x ) mean( x( :, 3 )) .* microns_per_voxel( 3 ), strand_subscripts );

end % FUNCTION