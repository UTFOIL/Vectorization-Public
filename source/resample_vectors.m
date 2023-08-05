function [ size_of_image, lumen_radius_in_pixels_range, edge_subscripts, vessel_directions ] ...
                                            = resample_vectors( lumen_radius_in_pixels_range, resolution_factors, edge_subscripts, size_of_image, varargin )
%% exerpted from flow_field_subrouting SAM 4/13/19                                        
                                        
%% interpolation of size look up table

% resolution_factors = resolution_factor * [ 1, 1, z_per_xy_length_of_pxl_ratio ];

number_of_scales = size( lumen_radius_in_pixels_range, 1 );

increased_number_of_scales = ceil( resolution_factors( 4 ) * ( number_of_scales - 1 )) / resolution_factors( 4 );

scale_sampling_range = 0 : 1 / resolution_factors( 4 ) : increased_number_of_scales ;

lumen_radius_in_pixels_range = exp( interp1(( 0 : number_of_scales - 1 )',       ...
                                            log( lumen_radius_in_pixels_range ), ...
                                            scale_sampling_range                 ));

lumen_radius_in_pixels_range = lumen_radius_in_pixels_range .* resolution_factors( 1 : 3 );    

% Extension of image
size_of_image = ceil(( size_of_image - 1 ) .* resolution_factors( 1 : 3 )) + 1 ;
                                        
%% space and scale interpolation of vector subscripts
% % precalculate cumulative distances covered by the vectors in each edge for later interpolation (units of the new voxel lenghts) 

edge_subscripts = cellfun( @( v ) resolution_factors  .* ( v - 1 ) + 1, edge_subscripts, 'UniformOutput', false );                                
% edge_subscripts = cellfun( @( v ) unique( resolution_factors  .* ( v - 1 ) + 1, 'rows', 'stable' ), edge_subscripts, 'UniformOutput', false );                                

if ~ isempty( varargin ), vessel_directions = varargin{ 1 }; end

try 
    
    normed_space = varargin{ 2 }; 
%     normed_space = 'L^2' ;

catch

    normed_space = 'L^inf' ;

end

switch normed_space

    case 'L^inf'
        % using the L-infinity norm to get 1 voxel length spacing instead of the L-2 norm for real-space
        % distances
        edge_cumulative_lengths                                                                                ...
                        = cellfun( @( v ) [ 0; cumsum( max( abs((   v( 1 + 1 : end    , 1 : 3 )                ...
                                                                  - v( 1     : end - 1, 1 : 3 ))), [ ], 2 ))], ...
                                   edge_subscripts, 'UniformOutput', false                                     );
        % edge_cumulative_lengths                                                                                ...
        %                 = cellfun( @( v ) unique([ 0; cumsum( max( abs((   v( 1 + 1 : end    , 1 : 3 )                ...
        %                                                           - v( 1     : end - 1, 1 : 3 ))), [ ], 2 ))], 'stable' ), ...
        %                            edge_subscripts, 'UniformOutput', false                                     );
    
    case 'L^2'
        % using the L-2 norm to get real-space distances
        edge_cumulative_lengths                                                                                  ...
                        = cellfun( @( v ) [ 0; cumsum( sum((   v( 1 + 1 : end    , 1 : 3 )                       ...
                                                             - v( 1     : end - 1, 1 : 3 )) .^ 2, 2 ) .^ 0.5 )], ...
                                   edge_subscripts, 'UniformOutput', false                                     );
        
end

edge_sample_lengths = cellfun( @( v ) linspace( 0, v( end ), ceil( v( end )) + 1 )',  ...
                               edge_cumulative_lengths, 'UniformOutput', false         );
                               
% interpolate the vector positions and size indices along the edge index, round to closest
% subscripts, then remove redundant vectors with unique function
number_of_edges = length( edge_subscripts );

for edge_index = 1 : number_of_edges
        
%     [ ~, position_indices_1 ] = unique( edge_subscripts{ edge_index }( :, 1 : 3 ), 'rows', 'stable' );
%     
%             edge_subscripts{ edge_index } =         edge_subscripts{ edge_index }( position_indices_1, : );
%     edge_cumulative_lengths{ edge_index } = edge_cumulative_lengths{ edge_index }( position_indices_1    );

    if edge_cumulative_lengths{ edge_index }( end ) == 0

        switch normed_space
        
            case 'L_inf'
                edge_subscripts{ edge_index } = round( edge_subscripts{ edge_index });

            case 'L^2'

                % nothing
        end
    else
        switch normed_space

            case 'L_inf'
                          [ ~, unique_indices ]                                                                 ...
                                           = unique( round( interp1( edge_cumulative_lengths{ edge_index },     ...
                                                                             edge_subscripts{ edge_index },     ...
                                                                         edge_sample_lengths{ edge_index }  )), ...
                                                     'rows', 'stable'                                           );
            % !!!!! 'stable' call in unique funtion above is extremely important for smooth overwrite rendering
            % of flow field or graded visual output !!!!!
        end

                edge_subscripts{ edge_index }         = interp1( edge_cumulative_lengths{ edge_index }, ...
                                                                         edge_subscripts{ edge_index }, ...
                                                                     edge_sample_lengths{ edge_index }  );
    
        switch normed_space
        
            case 'L-inf'
            
                edge_subscripts{ edge_index } = edge_subscripts{ edge_index }( unique_indices, : );
            
        end


        if ~ isempty( varargin )
    
    %         vessel_directions{ edge_index } = vessel_directions{ edge_index }( position_indices_1, : );
    
                    vessel_directions{ edge_index }       = interp1( edge_cumulative_lengths{ edge_index }, ...
                                                                           vessel_directions{ edge_index }, ...
                                                                         edge_sample_lengths{ edge_index }  );
    
            switch normed_space
            
                case 'L-inf'
            
                    vessel_directions{ edge_index } = vessel_directions{ edge_index }( unique_indices, : );

            end
        end
    end
end % FOR
                       
end

