function center_of_area = calculate_center_of_area( strand_subscripts )

% !!! exerpted from calculate_depth_statistics

strand_delta_lengths    = cellfun( @( x ) sum(((   x( 2 : end    , 1 : 3 )                                           ...
                                                 - x( 1 : end - 1, 1 : 3 ))  ) .^ 2, 2 ) .^ 0.5, ...
                                   strand_subscripts, 'UniformOutput', false                                         );
                                                           
% strand inner-position radii
strand_inner_pos_radii      = cellfun( @( x ) (   x( 2 : end     , 4 )          ...
                                                + x( 1 : end - 1 , 4 )) / 2,    ...
                                       strand_subscripts, 'UniformOutput', false );

% strand inner-position y-position
strand_inner_pos_y          = cellfun( @( x ) (   x( 2 : end    , 1 )                                 ...
                                                + x( 1 : end - 1, 1 )) / 2 , ...
                                       strand_subscripts, 'UniformOutput', false                      );

% strand inner-position z-position
strand_inner_pos_x          = cellfun( @( x ) (   x( 2 : end    , 2 )                                 ...
                                                + x( 1 : end - 1, 2 )) / 2 , ...
                                       strand_subscripts, 'UniformOutput', false                      );

% strand inner-position z-position
strand_inner_pos_z          = cellfun( @( x ) (   x( 2 : end    , 3 )                                 ...
                                                + x( 1 : end - 1, 3 )) / 2 , ...
                                       strand_subscripts, 'UniformOutput', false                      );

delta_lengths           = cell2mat( strand_delta_lengths           );
inner_pos_radii         = cell2mat( strand_inner_pos_radii         );
inner_pos_y             = cell2mat( strand_inner_pos_y             );
inner_pos_x             = cell2mat( strand_inner_pos_x             );
inner_pos_z             = cell2mat( strand_inner_pos_z             );

inner_area = delta_lengths .* inner_pos_radii ; % off by constant factor, 2*pi

total_inner_area = sum( inner_area );

center_of_area = sum([ inner_pos_y, inner_pos_x, inner_pos_z ] .* inner_area ) / total_inner_area ;

end %FUCNTION