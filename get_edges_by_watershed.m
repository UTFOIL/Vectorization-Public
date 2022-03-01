function [ edge_locations,    ...
           edge_energies ,    ...
           edges2vertices,    ...
                  energy_map, ...
            vertex_index_map, ...
                 pointer_map, ...
                d_over_r_map, ...
            branch_order_map   ]                                 ...
                = get_edges_by_watershed(   energy_tolerance, edge_number_tolerance,      ...
                                          distance_tolerance, local_distance_tolerance_range, size_tolerance,      ...
                                          strel_linear_LUT_range, local_subscripts_range, strel_r_over_R_LUT_range, strel_unit_vectors_LUT_range,...
                                          vertex_locations, energy_map, size_map, microns_per_voxel )
%% initialization:
% vertex_index_map has all the vertices marked with unique indices in order from best to worst
% energy. pointer_map has all the strels surrounding the vertices filled in with pointers pointed
% toward their centers. Avaliable energy map has just the vertices and their strels filled in.

% method = 'four' ;

% number_of_edges_at_vertex_max = 3 ; % # of edges connecting at origin vertex
%              branch_order_max = 2 ; % # 1 means just parent traces, 2 means a generation of children

%     branch_order_max = 3  ; % SAM 11/4/21
% max_moves_per_vertex = 10 ; % SAM 11/4/21
% max_moves_per_vertex = 20 ; % SAM 11/4/21

% distance_map = [ ];
% branch_order_map = [ ];

% is_binarized = false ; % hard-coded input

% is_binarized_1 = true ;

% direction_tolerance = 1 / 3 ; % exponent of unitless ratio in [ 0, 1 ]
% direction_tolerance = 2 ; % exponent of unitless ratio in [ 0, 1 ]

% strel_length = numel( strel_distance_LUT_range ) .^ ( 1 / 3 );
% 
% strel_distance_LUT_range = reshape( strel_distance_LUT_range, strel_length, strel_length, strel_length );
% 
% microns_per_voxel = [ strel_distance_LUT_range(( end + 1 ) / 2 + 1, ( end + 1 ) / 2    , ( end + 1 ) / 2     ),...
%                       strel_distance_LUT_range(( end + 1 ) / 2    , ( end + 1 ) / 2 + 1, ( end + 1 ) / 2     ),...
%                       strel_distance_LUT_range(( end + 1 ) / 2    , ( end + 1 ) / 2    , ( end + 1 ) / 2 + 1 ) ];
%                   
% strel_displacement_LUT = strel_distance_LUT_range ;     strel_displacement_LUT( 1 : ( end + 1 ) / 2 )...
%                                               = - strel_displacement_LUT( 1 : ( end + 1 ) / 2 );
%                   
% strel_unit_vectors_LUT_range ...
%     = cat( 4, strel_displacement_LUT( :, ( end + 1 ) / 2,    ( end + 1 ) / 2    ) .* ones( size( strel_distance_LUT_range )),...
%               strel_displacement_LUT(    ( end + 1 ) / 2, :, ( end + 1 ) / 2    ) .* ones( size( strel_distance_LUT_range )),...
%               strel_displacement_LUT(    ( end + 1 ) / 2,    ( end + 1 ) / 2, : ) .* ones( size( strel_distance_LUT_range )) );
%           
% strel_unit_vectors_LUT_range =          strel_unit_vectors_LUT_range    ./   strel_distance_LUT_range       ;
% 
% strel_unit_vectors_LUT_range = reshape( strel_unit_vectors_LUT_range, numel( strel_distance_LUT_range ), 3 );
% 
% strel_unit_vectors_LUT_range(( end + 1 ) / 2, : ) = 0 ; % replace undefined 0/0 with 0
% 
% strel_distance_LUT_range = strel_distance_LUT_range( : );

% pointer_map is indices into the linear_strel. pointers show how to get back to the origin vertex,
% so the relative pointer indices are the negative (or reverse order) of the relative strel indices
% pointer_template = transpose( uint8( numel( linear_strel )) : - 1 : 1 );
%      strel_pointers_LUT_range = transpose( numel( strel_linear_LUT_range ) : - 1 : 1 );
% strel_pointers_LUT_range = cellfun( @( x ) transpose( numel( x ) : - 1 : 1 ), strel_linear_LUT_range, 'UniformOutput', false );
strel_pointers_LUT_range = cellfun( @( x ) transpose( 1 : numel( x )), strel_linear_LUT_range, 'UniformOutput', false );
    
%     distances_template = ;

% dot multiply to get cos(theta) with the other direction

number_of_vertices = length( vertex_locations );
    
% % flip the vertex indices and add one just for this function, so that higher index is lower energy
% backwards_indices = number_of_vertices : -1 : 1 ;
% 
% vertex_locations = int64( vertex_locations( backwards_indices ));
% 
% original_vertex_indices = [ number_of_vertices + 1; original_vertex_indices( backwards_indices )];  % background code concatenated on

vertex_locations = int64( vertex_locations );


% border_index = uint32( 1 );

% numel_in_image = numel( energy_map );

size_of_image = int64( size( energy_map ));

cum_prod_image_dims = cumprod( size_of_image );

% strel_apothem = max( round( strel_linear_LUT_range / cum_prod_image_dims( 2 )));
strel_apothem = 1 ;

% Linear indexing of the borders of the image:
border_locations_Y = [ 0 : strel_apothem                            - 1, cum_prod_image_dims( 1 ) - strel_apothem                            : cum_prod_image_dims( 1 ) - 1 ]' + ( 1 : cum_prod_image_dims( 1 ) : cum_prod_image_dims( 3 ));
border_locations_X = [ 0 : strel_apothem * cum_prod_image_dims( 1 ) - 1, cum_prod_image_dims( 2 ) - strel_apothem * cum_prod_image_dims( 1 ) : cum_prod_image_dims( 2 ) - 1 ]' + ( 1 : cum_prod_image_dims( 2 ) : cum_prod_image_dims( 3 ));
border_locations_Z = [ 0 : strel_apothem * cum_prod_image_dims( 2 ) - 1, cum_prod_image_dims( 3 ) - strel_apothem * cum_prod_image_dims( 2 ) : cum_prod_image_dims( 3 ) - 1 ]' + ( 1 : cum_prod_image_dims( 3 ) : cum_prod_image_dims( 3 ));

border_locations = unique([ border_locations_Y( : ); border_locations_X( : ); border_locations_Z( : )]);
                    
% available_locations_map = zeros( size_of_image, 'logical' ); available_locations_map( vertex_locations ) = 1                          ;
branch_order_map = zeros( size_of_image, 'uint8'   );
vertex_index_map = zeros( size_of_image, 'uint32'  ); vertex_index_map( border_locations ) = 1 + number_of_vertices ; ...
                                                      vertex_index_map( vertex_locations ) = 1 : number_of_vertices ; % last index is for the border N+first "original" vertex
     pointer_map = zeros( size_of_image, 'uint64'  ); % pointers point to the center of the strel so that when watersheds join, we can trace back to their origin vertices
    d_over_r_map = zeros( size_of_image            ); % distance traversed normalized by the radius of origin vertex
    
energy_map_backup = energy_map ;

% moves_per_vertex = zeros([ number_of_vertices + 1 , 1 ], 'uint8' );
        
% number_of_vertices_found = 0 ;
% index_of_last_vertex_found = number_of_vertices + 2 ; % one for background and one to subtract later

number_of_edges = uint32( 0 );

% number_of_edges_per_vertex = zeros( number_of_vertices + 1, 1 );

edges2vertices  = zeros( 0, 2, 'uint32' );
edge_locations  =      {                }; % empty cell array

% edge_locations_to_reset = {[ ], [ ]};

vertex_adjacency_matrix ...
           = logical( spalloc( number_of_vertices + 1, ... % the image border adds 1 "vertex" index
                               number_of_vertices + 1, ... 
                             ( number_of_vertices + 1 ) * round( edge_number_tolerance ))); 
vertex_adjacency_matrix( 1 : ( number_of_vertices + 2 )       ...
                           : ( number_of_vertices + 1 ) ^ 2 ) = 1 ; % mark all vertices as self-adjacent: sparse identity matrix results

%% main WHILE: tracing to lowest available energy location
% available_locations = vertex_locations ;
available_locations = vertex_locations( end : -1 : 1 ); % switch the order of the vertex locations so the best locations are at the bottom of the list, and we don't have to look far to get to them in the FIND( ___, 'last' ) call later

while true
    %% choose next index
    
    % requiring an available move
    if isempty( available_locations ), break, end        
    
    current_location           = available_locations( end ); % pre-sorted % SAM 1/27/22

%     [ min_available_energy, current_available_index ] = min( energy_map( available_locations ));
%     [ min_available_energy, current_available_index ] = min( energy_map_backup( available_locations ));
    min_available_energy = energy_map_backup( current_location ); % pre-sorted % SAM 1/27/22
    
%     % requiring an available move
%     if isempty( current_available_index ), break, end    

%     current_location           = available_locations( current_available_index );        
    current_strel_vertex_index =    vertex_index_map( current_location );
    is_new_trace_from_current_location = pointer_map( current_location ) == 0 ;        
    
%     if ~ is_new_trace_from_current_location % only check this if it is not a vertex location to save time
        % requiring negative edge energy to explore    
        if min_available_energy >= 0, break, end
%         % requiring edge energy above the tolerance to continue exploring !!!!!! not fair to verticies with worse energy
%         if min_available_energy >= energy_map( vertex_locations( current_strel_vertex_index )) ...
%                                  * ( 1 - energy_tolerance )                                  , break, end

%     end

    
%     available_locations_map( current_location ) = false ;

%     available_locations( available_locations == current_location ) = 1 ; % code for destroyed
%     available_locations( available_locations == current_location ) = [ ]; % slower than above but saves time when searching for next available minimum

% 	      moves_per_vertex( current_vertex_index )            = moves_per_vertex( current_vertex_index ) + 1 ;    
%   locations_with_indices{ current_vertex_index }( end + 1 ) =                   current_location ;        

    current_strel = strel_linear_LUT_range{ round( size_map( current_location ))} + current_location ;

    is_current_strel_in_map = current_strel >= 1                      ...
                            & current_strel <= cum_prod_image_dims( 3 );

    current_strel( ~ is_current_strel_in_map ) = current_location ; % to be removed later

   [ current_subscript_y,                             ...
     current_subscript_x,                             ...
     current_subscript_z ] = ind2sub(  size_of_image, ...
                                     current_location ); current_subscripts = [ current_subscript_y, ...
                                                                                current_subscript_x, ...
                                                                                current_subscript_z ];

    absolute_subscripts_range = local_subscripts_range{ round( size_map( current_location ))} + current_subscripts ;

    is_current_strel_in_map = logical( prod(   absolute_subscripts_range >= [ 1, 1, 1 ]             ...
                                             & absolute_subscripts_range <= size_of_image, 2 ));

    current_strel( ~ is_current_strel_in_map ) = [ ];

    current_strel_unit_vectors       = strel_unit_vectors_LUT_range{ round( size_map( current_location ))};
    current_strel_pointers_LUT_range =     strel_pointers_LUT_range{ round( size_map( current_location ))};
    current_strel_r_over_R_LUT_range =     strel_r_over_R_LUT_range{ round( size_map( current_location ))};

    current_strel_unit_vectors(        ~ is_current_strel_in_map, : ) = [ ];
    current_strel_pointers_LUT_range(  ~ is_current_strel_in_map    ) = [ ];
    current_strel_r_over_R_LUT_range(  ~ is_current_strel_in_map    ) = [ ];

    vertices_of_current_strel = vertex_index_map( current_strel );

    %% redefine the energy for the current neighborhood (strel), accounting for...

    current_strel_energies = energy_map_backup( current_strel );    
%     current_strel_energies = energy_map( current_strel ); % SAM 1/27/22 (and earlier?)

%         %% energy
% %     % redefine energy in terms of similarity to vertex eneergy to prevent low energy vertices from
% %     moving up
%     current_strel_energies =          ...    energy_map( vertex_locations( current_vertex_index )) * ...
%                   - 1e5 * exp( - ((    current_strel_energies ...
%                                      - energy_map( vertex_locations( current_vertex_index ))) ...
%                                *   3 / energy_map( vertex_locations( current_vertex_index ))   ) .^ 2 / 2 );
%     % redefine energy in terms of origin vertex
%         current_strel_energies = current_strel_energies ...
%                             / - energy_map( vertex_locations( current_vertex_index ));
%         current_strel_energies = abs(   current_strel_energies ...
%                                    / energy_map( vertex_locations( current_vertex_index )) ...
%                                    - 1 ) - 1 ;
%         current_strel_energies = max(     current_strel_energies ...
%                                    / - energy_map( vertex_locations( current_vertex_index )), ...
%                                    - 1 );
%         current_strel_energies = energy_map( vertex_locations( current_vertex_index )) ...
%                             * ( 1 - (   current_strel_energies ...
%                                       / energy_map( vertex_locations( current_vertex_index )) ...
%                                       - 1                                                         ) .^ 2 );
%         current_strel_energies = energy_map( vertex_locations( current_vertex_index )) ...
%                             * exp( - ((    current_strel_energies ...
%                                          - energy_map( vertex_locations( current_vertex_index ))) ...
%                                      * 6 / energy_map( vertex_locations( current_vertex_index ))   ) .^ 2 / 2 );
%         current_strel_energies = current_strel_energies ...
%                            .* exp( - ((    current_strel_energies ...
%                                          - energy_map( vertex_locations( current_vertex_index ))) ...
%                                      * 3 / energy_map( vertex_locations( current_vertex_index ))   ) .^ 2 / 2 );

        %% size

    size_index_differences = size_map( current_strel ) - size_map( vertex_locations( current_strel_vertex_index ));

    size_index_energy_adjustment_factors = exp( - ( size_index_differences / size_tolerance ) .^ 2 / 2 ); % gaussian with max of 1 and index_tolerance as standard deviation

    current_strel_energies = current_strel_energies .* size_index_energy_adjustment_factors ;

%         if exist( 'distance_tolerance_range', 'var' )
        %% distance (local)
%         if exist( 'local_distance_tolerance_range', 'var' )

% %         distance_energy_adjustment_factors = -1 + 2 .^ (  strel_distance_LUT_range{              size_map( current_location )} ...
% %                                                         / local_distance_tolerance_range( round( size_map(  vertex_locations( current_vertex_index )))));
%         distance_energy_adjustment_factors = ( 1 - max( exp( - (   current_strel_distance_LUT_range ...
%                                                                  /   local_distance_tolerance_range( round( size_map(  vertex_locations( current_vertex_index ))))) .^ 2 / 2 ), ...
%                                                         exp( - 1 / 2 )                                                                                                           )) ...
%                                            / ( 1 -      exp( - 1 / 2 ));

%     distance_energy_adjustment_factors = tanh((      current_strel_distance_LUT_range ...
%                                                  *   1.5                              ...
%                                                  /   local_distance_tolerance_range( round( size_map(  vertex_locations( current_vertex_index ))))) .^ 2 );

%     distance_energy_adjustment_factors = ( 1 - cos( pi * current_strel_distance_LUT_range ...
%                                                      / local_distance_tolerance_range( round( size_map(  vertex_locations( current_vertex_index )))))) ...
%                                        / 2 ; % SAM 1/20/22

%     distance_energy_adjustment_factors = ( 1 - sin( 2 * pi * max( min(   current_strel_distance_LUT_range ...
%                                                                        / local_distance_tolerance_range( round( size_map(  vertex_locations( current_vertex_index )))), ...
%                                                                        3/4 ), ...
%                                                                   1/4 )       )) / 2 ; % SAM 1/21/22

%     distance_energy_adjustment_factors = exp( - 8 * max( 0, 3 / 4 -   current_strel_distance_LUT_range ...
%                                                                     / local_distance_tolerance_range( round( size_map( vertex_locations( current_strel_vertex_index ))))) .^ 2 ); % SAM 1/21/22

%     distance_energy_adjustment_factor = sin( min( pi / 2, ...
%                                                    pi * current_strel_distance_LUT_range ...
%                                                       / 2    /  local_distance_tolerance_range( round( size_map( vertex_locations( current_strel_vertex_index )))))); % SAM 1/25/22

%     distance_energy_adjustment_factor = sin( pi / 2 * min( 1, current_strel_r_over_R_LUT_range )); % SAM 1/28/22
    distance_energy_adjustment_factor = ( 1 - cos( pi * min( 1, 4/3 * current_strel_r_over_R_LUT_range ))) / 2 ; % x and dx are 0 at r = 0, and x = 1 and dx is 0 at r = R. r/R is scaled to use more of the strel. cos is for a smooth transition at the donut whole into the donut  % SAM 2/4/22

    current_strel_energies = current_strel_energies .* distance_energy_adjustment_factor ;
        %% distance (total)
% 
%     distance_energy_adjustment_factors ...
%         = exp( - ((                 distance_map(                              current_location ) ...
%                     + current_strel_distance_LUT_range )...
%                   /                 distance_tolerance_range( round( size_map(  vertex_locations( current_strel_vertex_index ))))) .^ 2 / 2 / 2 ^ 2 ); % put the tolerance limit at a z score of 2 SAM 1/11/22
%     distance_energy_adjustment_factor = exp( - distance_map( current_location ) ^ 2 / 2 ); % put the tolerance limit at a z score of 1 SAM 1/11/22
    distance_energy_adjustment_factor = exp( - ( 3 * d_over_r_map( current_location ) / distance_tolerance ) ^ 2 / 2 ); % put the tolerance limit at a z score of 3 SAM 2/4/22

    current_strel_energies = current_strel_energies * distance_energy_adjustment_factor ;
% 
%         end % IF distance max constraint (for speed/memory)


%         end % IF distance max constraint (for speed/memory)    

%     if exist( 'edge_number_tolerance', 'var' )
% 
%         %% per vertex level
% %                 edge_number_energy_adjustment_factors ...
% %                        = exp( - abs(        number_of_edges_per_vertex( current_vertex_index )                                     ...
% %                                      / edge_number_tolerance                                   ) ^ ( 2 + is_binarized_1 * 10000 ) / 2 );                
% %                 edge_number_energy_adjustment_factors ...
% %                        = number_of_edges_per_vertex( current_vertex_index ) <= number_of_edges_at_vertex_max ;                
% % 
% %                 current_strel_energies = current_strel_energies * edge_number_energy_adjustment_factors ;
% 
%         %% current branch of this vertex level
% 
% %                 edge_number_energy_adjustment_factors ...
% %                        = exp( - abs( double( branch_order_map( current_location ))                                     ...
% %                                      / edge_number_tolerance                ) ^ ( 2 + is_binarized * 10000 ) / 2 );
%         edge_number_energy_adjustment_factors ...
%                = branch_order_map( current_location ) <= branch_order_max ;                
% 
%         current_strel_energies = current_strel_energies .* edge_number_energy_adjustment_factors ;
% 
%     end % IF edges per vertex constraint (connection complexity)
% 
        %% direction (local)

    if pointer_map( current_location ) % IF current location is not a vertex or edge
        direction_forwards = strel_unit_vectors_LUT_range{ round( size_map( current_location ))}( pointer_map( current_location ), : ); % ????? somethin wrong with the naming?

        % dot multiply to get cos(theta) with the other direction
%                 strel_units_along_direction_forwards = - sum(    current_strel_unit_vectors   ...
%                                                                .* direction_backwards,   2 );
        strel_units_along_direction_forwards = sum(    current_strel_unit_vectors   ...
                                                    .* direction_forwards,        2 );

        % not allowed to turn around (must go forward)
        strel_units_along_direction_forwards( strel_units_along_direction_forwards < 0 ) = 0 ;

%         current_strel_energies = current_strel_energies .*          strel_units_along_direction_forwards .^ ( direction_tolerance / 2 );
        current_strel_energies = current_strel_energies .*          strel_units_along_direction_forwards  ; % SAM 1/21/22
%         current_strel_energies = current_strel_energies .*          strel_units_along_direction_forwards .^ 0.5 ; % SAM 1/24/22
%         current_strel_energies = current_strel_energies .* logical( strel_units_along_direction_forwards ); % must move forward, cannot move against direction of travel % SAM 1/21/22
%                 current_strel_energies = current_strel_energies .* ( strel_units_along_direction_forwards < 1/2 ); % 60 degree cone permitted in the directino of travel

    end
%         %% direction (total)
% %         if pointer_map( current_location )
%     % %     [ direction_home ] = get_direction_home( current_location );
%     %         
%     %                 % direct away from home location (vertex or parent strand)
%     %     
%     %                 if branch_order_map( current_location ) == 1 % home is a vertex
%     %     
%         home_location = vertex_locations( current_strel_vertex_index ); 
% 
%         direction_home = index2position(    home_location, cum_prod_image_dims ) ...
%                        - index2position( current_location, cum_prod_image_dims );
% 
%         direction_home = double( direction_home' ) .* microns_per_voxel ; % conversion from voxels to microns 
% 
%         direction_home = direction_home ./ sum( direction_home .^ 2 ) .^ 0.5 ;
%     %     
%     %                 else % home is a parent edge
%     %     
%     %                     tracing_location = current_location ;
%     %     
%     %                     direction_home = [ 0, 0, 0 ];
%     %     
%     %                     % trace back until hitting a zero in the pointer index map.
%     %                     while true
%     %     
%     %                         if pointer_map( tracing_location ) == 0, break, end
%     %     
%     %                         direction_home   =   direction_home + strel_unit_vectors_LUT_range{ size_map( tracing_location )}( pointer_map( tracing_location ), : );
%     %     
%     %                         tracing_location = tracing_location +       strel_linear_LUT_range{ size_map( tracing_location )}( pointer_map( tracing_location )    );
%     %     
%     %                     end % WHILE tracing
%     %     
%     %                     home_location = tracing_location ;
%     %     
%     %                 end % IF tracing will end at vertex
%     %     
%         % dot multiply to get cos(theta) with the other direction
%         strel_units_away_from_home = - sum(    current_strel_unit_vectors   ...
%                                             .* direction_home,        2 );
% 
%         % not allowed to turn around (must go forward)
%         strel_units_away_from_home( strel_units_away_from_home < 0 ) = 0 ;
% 
% %             current_strel_energies = current_strel_energies .* strel_units_away_from_home .^ ( direction_tolerance / 2 );                    
%         current_strel_energies = current_strel_energies .* strel_units_away_from_home .^ 0.5 ;                    
% %     
% %             else % ELSE current location is a parent edge or a vertex)
% %     
% %                 home_location = current_location ;
% %     
%     end % IF pointer map is nonzero at current location                    
% % 
% % %                 if is_constrained_by_parent_edge_direction 
% % 

    %% write the newly revealed neighborhood to the energy, vertex, pointer, distance, and size maps

%     is_not_self_vertex_in_strel   = vertices_of_current_strel ~= current_vertex_index ;    
%     is_not_self_vertex_in_strel = vertices_of_current_strel == 0 ;     % unnassigned to a vertex is not self
%     is_not_self_vertex_in_strel(  vertices_of_current_strel ~= 0 ) ... % and where it is nonzero ...
%         = ~ vertex_adjacency_matrix( current_vertex_index, nonzeros( vertices_of_current_strel )); % check the (NOT-)adjacency matrix for "(ANTI-)self-ness"
    
    is_without_vertex_in_strel = vertex_index_map( current_strel ) == 0 ;
    
% %     is_new_in_strel =         pointer_map( current_strel ) == 0 ...
% %                     &    vertex_index_map( current_strel ) == 0 ;...    % locations that have not been explored by any vertex    
% %                 ...        & current_strel_energies                   < 0 ;
% %             ...            & is_energy_tolerated_in_strel ;
%     is_new_in_strel =     is_unclaimed_to_any_vertex_in_strel ...
%                     | (   is_not_self_vertex_in_strel ...
%                         & is_energy_lower_in_strel    );
% 
%                       current_strel_new   = current_strel(                    is_new_in_strel );
%                       
%     vertex_index_map( current_strel_new ) = current_vertex_index ;
%           energy_map( current_strel_new ) = current_strel_energies(           is_new_in_strel );
%          pointer_map( current_strel_new ) = current_strel_pointers_LUT_range( is_new_in_strel );
%         distance_map( current_strel_new ) = current_strel_distance_LUT_range( is_new_in_strel ) ...
%                                           + distance_map( current_location  );
                      current_strel_unclaimed   = current_strel(                    is_without_vertex_in_strel );
    vertex_index_map( current_strel_unclaimed ) = current_strel_vertex_index ;
          energy_map( current_strel_unclaimed ) = current_strel_energies(           is_without_vertex_in_strel );    
         pointer_map( current_strel_unclaimed ) = current_strel_pointers_LUT_range( is_without_vertex_in_strel );
        d_over_r_map( current_strel_unclaimed ) = current_strel_r_over_R_LUT_range( is_without_vertex_in_strel ) ...
                                                + d_over_r_map( current_location  );
            size_map( current_strel_unclaimed ) = size_map(     current_location ); % overwriting size map so that all strels from traces belonging to the same vertex are the same size
%     branch_order_map( current_strel_new ) =     branch_order_map( current_location )...
%                                           + uint8( ~ pointer_map( current_location ));        % !!!!!! new definition 11/4/21                            

%     % overwrite the current strel energies if a previous vertex watershed has explored a part of it with better energy
%     is_energy_lower_from_Vertex_B_in_strel =       energy_map( current_strel ) <  current_strel_energies     ...
%                                            & vertex_index_map( current_strel ) ~= current_strel_vertex_index ...
%                                            & ~ is_without_vertex_in_strel ;
%     
% 	current_strel_energies(          is_energy_lower_from_Vertex_B_in_strel ) ...
%         = energy_map( current_strel( is_energy_lower_from_Vertex_B_in_strel )); 

    %% reveal the next best available move(s), combining watersheds if one (or many) are present
%     is_energy_tolerated_in_strel = current_strel_energies                                       ...
%                                  < energy_map( vertex_locations( current_vertex_index )) * ...
%                            ... - 1e5 * ...
%                                    ( 1 - energy_tolerance )                                  ;

%     % pull a new strel of energies from the energy map which will still include the newly added
%     % entries as well as any older entries that were had better energy and were not overwritten.
%     current_strel_energies = energy_map( current_strel );

%     is_energy_tolerated_in_strel = current_strel_energies ...
%                                  < energy_map( vertex_locations( current_vertex_index )) ...
%                                  * ( 1 - energy_tolerance )                                  ;

    if is_new_trace_from_current_location % IF moves are being revealed from a source (vertex or parent edge) origin

        seed_index_range = 1 : edge_number_tolerance     ; % each vertex can go in edge_number_tolerance number of directions

    else % The current location is along an exploratory edge trajectory, not at a source location
        % only keep the best candidate inside the strel for possible moves. Greatly improves speed,
        % and should not decrease the degrees of freedom of the 1D traces

%             seed_index_range = 1 : edge_number_tolerance    ; % SAM 12/18/21
%         seed_index_range = 1 : edge_number_tolerance - 1 ; % SAM 11/2/21
        seed_index_range = 1 ; 

    end

	is_current_location_clear = false ;    
    
    % sequentially select the best next positions (i.e. directions), from the available
    % moves in the neighborhood, and remove locations in the same direction as previously
    % chosen locations
    for seed_idx = seed_index_range % number of traces to eminate from the current trace location (1 for a continued trajectory)

%         is_energy_above_watershed_in_strel = energy_map_backup( current_strel ) ...
%                                            > energy_map_backup( vertex_locations( current_strel_vertex_index ));
%                                        
%         current_strel_energies( is_energy_above_watershed_in_strel ) = 0 ;

        is_energy_tolerated_in_strel = current_strel_energies ...
                                     < energy_map_backup( vertex_locations( current_strel_vertex_index )) ...
                                     * ( 1 - energy_tolerance )                                  ;        

% %         is_energy_tolerated_in_strel = is_energy_tolerated_in_strel .... % also enforce that vertices can only explore the watershed in the "uphill" direction (edge traces will be monotonic increasing for the first half, and monotonic decreasing for the second half, in energy. (if the next best move is downstream towards another watershed, then it is illegal (because it is not within the current vertex's watershed).
% %                                      & current_strel_energies ...
% %                                      > min_available_energy ;        
% %         is_energy_tolerated_in_strel = is_energy_tolerated_in_strel .... % also enforce that vertices can only explore the watershed in the "uphill" direction (edge traces will be monotonic increasing for the first half, and monotonic decreasing for the second half, in energy. (if the next best move is downstream towards another watershed, then it is illegal (because it is not within the current vertex's watershed).
% %                                      & current_strel_energies ...
% %                                      > energy_map_backup( vertex_locations( current_strel_vertex_index ));        
% %         is_energy_tolerated_in_strel = is_energy_tolerated_in_strel .... % also enforce that vertices can only explore the watershed in the "uphill" direction (edge traces will be monotonic increasing for the first half, and monotonic decreasing for the second half, in energy. (if the next best move is downstream towards another watershed, then it is illegal (because it is not within the current vertex's watershed).
% %                                      & energy_map_backup( current_strel ) ...
% %                                      > energy_map_backup( current_location );        
%         is_energy_tolerated_in_strel = is_energy_tolerated_in_strel .... % also enforce that vertices can only explore the watershed in the "uphill" direction (edge traces will be monotonic increasing for the first half, and monotonic decreasing for the second half, in energy. (if the next best move is downstream towards another watershed, then it is illegal (because it is not within the current vertex's watershed).
%                                      & energy_map_backup( current_strel ) ...
%                                      > energy_map_backup( vertex_locations( current_strel_vertex_index ));        

%         if any( is_not_self_vertex_in_strel ) % IF any possible moves exist in this strel

        % pick the best location from the locations that have not been assigned to the current
        % vertex yet. If this location falls in uncharted territory, then record the location as
        % a possible move, otherwise, this is the end of this trace.
%         [ ~, strel_idx ] = min(( current_strel_energies - 1 ) .* is_not_self_vertex_in_strel ); !!!! this was a mistake because, we ... dido as below 2 lines
        [ ~, strel_idx ] = min( current_strel_energies );
%         [ ~, strel_idx ] = min(( energy_map( current_strel ) - 1 ) .* is_not_self_vertex_in_strel ); % !!!! this was a mistake because, we want the next best move to be chosen, and then evaluate whether we want to keep it, not the other way around

        next_location     =             current_strel( strel_idx );
        next_vertex_index = vertices_of_current_strel( strel_idx );
        
%             % !!!!!!!!!!!!!!!!! % reset the energy for this location % SAM 11/3/21
%             energy_map( next_location ) = current_strel_energies( strel_idx );
        if ~ is_energy_tolerated_in_strel( strel_idx )
%             if ~ is_current_location_clear, available_locations( available_locations == current_location ) = [ ];  is_current_location_clear = true ; end
            if ~ is_current_location_clear, is_current_location_clear = true ; available_locations( end ) = [ ]; end
        else
            if next_vertex_index == 0 % IF the next location has not yet been explored by another vertex before

                branch_order_map( next_location ) = branch_order_map( current_location ) + seed_idx - 1 ;

                if branch_order_map( next_location ) < edge_number_tolerance

%                             available_locations_map( next_location ) = true ;

%                             available_locations( end + 1 ) = next_location ; %#ok<AGROW> % slower than above but saves time when searching for next available minimum 
%                     if ~ is_current_location_cleared, available_locations( available_locations == current_location ) = next_location ;  is_current_location_cleared = true ; 
%                     else                              available_locations(                 end + 1                 ) = next_location ; %#ok<AGROW> 
%                     end

                    
                    if seed_idx == 1 % for speed, assume that a primary seed will be closer to the bottom (better) of the list of possible energies, and a secondary seed will be closer to the top (worse) of the list
                         if                                 energy_map_backup( available_locations(  1  )) ...
                                                         <= energy_map_backup(      next_location        )
                                                     
                                location_idx = 1 ;
                         else
                             
                                location_idx = 1 +  find(    energy_map_backup( available_locations       )            ...
                                                          >  energy_map_backup(      next_location        ), 1, 'last'  );
                         end
                    else,if                                  energy_map_backup( available_locations( end )) ...
                                                          >= energy_map_backup(      next_location        )
                                                      
                            if is_current_location_clear                     
                                location_idx = 1 + numel(                       available_locations                     );
                            else
                                location_idx = 	   numel(                       available_locations                     );
                            end
                        else
                                location_idx =      find(    energy_map_backup( available_locations       )             ...
                                                        <  energy_map_backup(      next_location        ), 1, 'first' );
                         end                            
                    end
%                     if isempty( location_idx )
%                         warning('here')
%                     end
                    if  ~ is_current_location_clear 
                          is_current_location_clear = true ;
                          available_locations = [ available_locations( 1 : location_idx - 1 ); next_location; available_locations( location_idx : end - 1 )]; % pre-sorted % SAM 1/27/22  
                    else, available_locations = [ available_locations( 1 : location_idx - 1 ); next_location; available_locations( location_idx : end     )]; % pre-sorted % SAM 1/27/22  
                    end                    
                end
            else % ELSE next move already belongs to a vertex

                %% join watersheds if it is not self

    %             % IF the next location is on the list of available locations, remove it
    %             available_locations_map( next_location ) = false ;
    % 
    %             available_locations( available_locations == next_location ) = [ ]; % remove the

                % IF locations belonging to the next vertex in the current strel are on the list of
                % available locations, remove them in a FOR
                is_next_vertex_in_strel = vertices_of_current_strel == next_vertex_index ;

%                 locations_to_reset = intersect( available_locations, current_strel( is_next_vertex_in_strel ))';
% 
%     %                 available_locations_map( locations_to_reset ) = false ;
% 
%                 for location = locations_to_reset, available_locations( available_locations == location ) = [ ]; end

                locations_to_reset = intersect( available_locations, current_strel( is_next_vertex_in_strel ))' ;

%                 if ~ is_current_location_clear,  locations_to_reset = [ locations_to_reset,  current_location ]; is_current_location_clear = true ; end %#ok<AGROW>
                if ~ is_current_location_clear, is_current_location_clear = true ; ...
                        locations_to_reset( locations_to_reset == available_locations( end )) = [ ]; ...
                                                                  available_locations( end )  = [ ]; ...
                end 
                                   
                numel_locations_to_reset = numel( locations_to_reset );

                locations_to_reset_reindexed = zeros( numel_locations_to_reset, 1 );
                
                idx_range = 1 : numel_locations_to_reset ;
                
                for idx = idx_range , locations_to_reset_reindexed( idx ) = find( available_locations == locations_to_reset( idx ), 1, 'first' ); end
                    
                available_locations( locations_to_reset_reindexed ) = [ ]; % SAM 1/27/22 condensed all array downsizing to one line (down from two seperate lines (one of which was inside the FOR loop))

                if ~ vertex_adjacency_matrix( next_vertex_index, current_strel_vertex_index ) % not( all( empty )) is false

                    number_of_edges = number_of_edges + 1 ;                

                    % Record the identities of the vertices involved.
                    edges2vertices( number_of_edges, [ 1, 2 ]) = [ current_strel_vertex_index, next_vertex_index ];

                    % label these vertices as adjacent on the adjacency matrix
                    vertex_adjacency_matrix( current_strel_vertex_index,          next_vertex_index ) = true ;
                    vertex_adjacency_matrix(          next_vertex_index, current_strel_vertex_index ) = true ;

    %                 is_location_from_vertex_B_in_strel = vertices_of_current_strel == next_vertex_index;            

    % %                 vertices_of_current_strel( vertices_of_current_strel == vertices_at_edge( 2 )) = current_vertex_index ; % reset these markers to continue seraching for more neighbor connections
    % % 
    % %                 is_not_self_vertex_in_strel   = vertices_of_current_strel ~= current_vertex_index ;    
    %                 is_not_self_vertex_in_strel = vertices_of_current_strel == 0 ;     % unnassigned to a vertex is not self
    %                 is_not_self_vertex_in_strel(  vertices_of_current_strel ~= 0 ) ... % and where it is nonzero ...
    %                     = ~ vertex_adjacency_matrix( current_vertex_index, nonzeros( vertices_of_current_strel )); % check the (NOT-)adjacency matrix for "(ANTI-)self-ness"

                    for edge_half = 1 : 2

                        if edge_half == 1 % half #1 is the first watershed to get there (lower energy route), Vertex A

                            tracing_location = current_location ; 

                        else % edge_half == 2, Vertex B
    
                            if next_vertex_index ~= 1 + number_of_vertices
                                
                                is_vertex_B_origin = pointer_map(      current_strel ) ==         0         ...
                                                   & vertex_index_map( current_strel ) == next_vertex_index ;

                                % If the vertex origin is in the current strel, then we should choose this point as the first point on the "other" half of the edge trace here    
                                if any( is_vertex_B_origin )

                                    next_location = current_strel( is_vertex_B_origin );                        

                                end
                            end
                            
                            tracing_location = next_location ;
                        
                        end

                        tracing_ordinate = 0 ;

                        % trace back until hitting a zero in the pointer index map, recording the linear indices
                        % into the image and setting the pointer index map to zero along the way.  If this
                        % process hits a zero index that doesn't belong to a vertex then record this for
                        % later (we will need to add a new vertex and break up the existing edge into two).
                        while true

                            % trace back from the point of watershed contact to the index origin (where
                            % there is no pointer, because it is a vertex, parent edge, or image boundary
                            tracing_ordinate = tracing_ordinate + 1 ;

                            edge_locations{ number_of_edges, edge_half }( tracing_ordinate ) = tracing_location ; %#ok<AGROW>

                            if pointer_map( tracing_location ) == 0, break, end 

            %                     new_location = tracing_location + strel_linear_LUT_range{ round( size_map( tracing_location ))}( pointer_map( tracing_location ));
                            new_location = tracing_location - strel_linear_LUT_range{ round( size_map( tracing_location ))}( pointer_map( tracing_location ));
            %                                                ------

            %                 % commented out 10/21/21 SAM #no fancy stuff !!!!!!!!! eliminating
            %                 children edges. only vertices % can source new trajectories (edges
            %                 cannot).
    %                         % Erase the pointers for the found trajectory, so that this path
    %         %                 will not contribute to % any other edge from this vertex.
    %         %     %                 clear_locations( tracing_location )
    %         %                     pointer_map( tracing_location ) = 0 ;
    %         %                 % END: commented out 10/21/21 SAM #no fancy stuff

                            tracing_location = new_location ;

                        end % WHILE tracing  

            %                 % commented out 10/21/21 SAM #no fancy stuff !!!!!!!!!
            %                 if branch_order_map( tracing_location ) == 0 && vertices_at_edge( edge_half ) ~= border_index % trace ended at vertex with one or zero edges
            % 
            %                     number_of_edges_per_vertex( vertices_at_edge( edge_half ))     ...
            %                   = number_of_edges_per_vertex( vertices_at_edge( edge_half )) + 1 ;
            % 
            %                 end
            % 
            %                 % !!!!!!!! this should depend on where the trace ended, and nothing to do with the
            %                 % number of edges vfor this vertex. Use branch_order_map instead?
            %                 if number_of_edges_per_vertex( vertices_at_edge( edge_half )) == 2
            % 
            %     %                 clear_locations( tracing_location )
            %                     pointer_map( tracing_location ) = 0 ;
            % 
            %                     branch_order_map( tracing_location ) = 1 ; % set the vertex branch order to 1 to make a continuous parent (order==1) backbone
            % 
            %                     number_of_edges_per_vertex( vertices_at_edge( edge_half )) = 3 ; % increment again so it is never chosen again (impossible to know if it is a vertex location easily after previous few linse.                
            % 
            %                     edge_locations_to_reset{ edge_half }  = edge_locations{ number_of_edges, edge_half }               ;
            % 
            %                 else
            % 
            %                     edge_locations_to_reset{ edge_half }  = edge_locations{ number_of_edges, edge_half }( 1 : end - 1 );
            % 
            %                 end % IF second edge for this vertex
            % 
            %                 vertex_index_map( edge_locations_to_reset{ edge_half }) = vertices_at_edge( edge_half );
            % 
                    end % FOR edge_half
            % 
            %             edge_locations_to_reset_mat = [ edge_locations_to_reset{ 1 },...
            %                                             edge_locations_to_reset{ 2 } ];

    %                 edge_locations_to_check = [ edge_locations{ number_of_edges, 1 },...
    %                                             edge_locations{ number_of_edges, 2 } ];
    % 
    %                 if numel( edge_locations_to_check ) ~= numel( unique( edge_locations_to_check ))
    % 
    %                     warning('repeated edge location')
    % 
    %                 end

            % 
            %     % commented out this conditional 10/12/21 SAM
            %     %         if vertices_at_edge( 2 ) ~= border_index 
            %     % 
            %     %             % allow children edge to come from these locations (if the second vertex is not the image boundary).
            %     %             available_locations_map( edge_locations_to_reset_mat ) = true ;
            %     % 
            %     %             available_locations( end + 1 : end + numel( edge_locations_to_reset_mat )) = edge_locations_to_reset_mat ;                    
            %     % 
            %     %         end
            %     % % commented in  this conditional 10/20/21 SAM
            % 
            %             energy_map( edge_locations_to_reset_mat ) = energy_map_backup( edge_locations_to_reset_mat );
            % 
            %                 % END: commented out 10/21/21 SAM #no fancy stuff

                end % WHILE joining watershed                     
            end % IF next best move is free or ELSE if it belongs to vertex already
        end
%                     end
% 
%                 end

        % remove from possibility the next positions that are in a similar direction to this
        % position
        cosine_of_current_min_to_strel ...
        = sum(    current_strel_unit_vectors                            ...
               .* current_strel_unit_vectors( strel_idx, : ), 2 );

%                    % graded weighting, not binary mask
%             current_strel_energies( is_new_in_strel ) = current_strel_energies( is_new_in_strel )         ...
%                                                   .* ( 1 - cosine_of_current_min_to_strel_new ) / 2 ;

%                 current_strel_energies( is_new_in_strel ) = current_strel_energies( is_new_in_strel )         ...
%                                                       .* ( cosine_of_current_min_to_strel_new ...
%                                                            <  2 ^ 0.5 / 2 ); % clear out a conic window at 45 degrees from the chosen direction % !!!!! make this window a plane instaed of cone

%         current_strel_energies = current_strel_energies         ...
%                            .* ( cosine_of_current_min_to_strel ...
%                             <  1 / 2 ); % clear out a conic window at 60 degrees from the chosen direction % !!!!! make this window a plane instaed of cone
        
        current_strel_energies = current_strel_energies         ... % SAM 1/25/22
                           .* ( 1 - cosine_of_current_min_to_strel ) / 2; % clear out a conic window using a continuous spatial filter that makes the chosen direction have zero energy, the opposite direction have best energy (factor of 1), and the orthogonal direction a factor of 1/2 as good of energy

%             current_strel_energies( is_new_in_strel ) = current_strel_energies( is_new_in_strel )         ...
%                                                   .* ( cosine_of_current_min_to_strel_new ...
%                                                        <  0 ); % clear out a hemisphere SAM 10/29/21

%             current_strel_energies( is_new_in_strel ) = current_strel_energies( is_new_in_strel )         ...
%                                                   .* ( cosine_of_current_min_to_strel_new ...
%                                                        <  - 1 / 2 ); % clear out an invertex cone of 60 degrees SAM 10/29/21 17:30

%         end % IF any moves exist: free space in the current neighborhood  
    end % FOR seed point per vertex
%         else % ELSE current neighborhood (strel) is centered on a origin belonging to an edge trace
%             
% 
%             [ ~, strel_idx ] = min(( current_strel_energies - 1 ) .* is_not_self_vertex_in_strel );
%     
%                       next_location =             current_strel( strel_idx );
%             vertex_of_next_location = vertices_of_current_strel( strel_idx );
% 
%             if vertex_of_next_location == 0 % IF the next location has not yet been explored by another vertex before
%             
%                 available_locations_map( next_location ) = true ;
% 
%                 available_locations( end + 1 ) = next_location ;
% 
%                 if pointer_map( next_location ) == 0
% 
%                     warning('next location has no pointer')
% 
%                 end
%             end
%         end
%     end % IF any moves exist: free space in the current neighborhood    
end % WHILE searching for edges

%% organizing output

% erase edge pairs that connect to boundary % !!!!!!!!!!!!!!!!! keep these insetad SAM 12/7/21 % but only keep the ones that touch the outer image boundary (not just the chunk boundary)% SAM 1/21/22
is_edge_connected_to_boundary = edges2vertices( :, 2 ) == 1 + number_of_vertices ;

edge_locations( is_edge_connected_to_boundary, : ) = [ ];
edges2vertices( is_edge_connected_to_boundary, : ) = [ ];

% % visulize maps
% mat2tif(               double( vertex_index_map ),                                               'vertex_index_map.tif' )
% mat2tif(                double( available_locations_map ),                                               'available_locations_map.tif' )
% mat2tif(                     double( pointer_map ),                                                    'pointer_map.tif' )
% mat2tif(                    distance_map,                                                   'distance_map.tif' )

% flip all of the first halves of edges so that they start with the vertex and end at the pass
edge_locations = cellfun( @( h1, h2 ) [ fliplr( h1 ), h2 ]', edge_locations( :, 1 ), edge_locations( :, 2 ), 'UniformOutput', false );
edge_energies  = cellfun( @(    x   )       energy_map( x ),                         edge_locations        , 'UniformOutput', false );

% scale the pointer map
numel_strel_linear_LUT_range = cellfun( @numel, strel_linear_LUT_range );

% pointer_map( pointer_map > 0 ) ...
%     = 1000 ./ numel_strel_linear_LUT_range( round( size_map( pointer_map > 0 ))) ...
%    .* double( pointer_map( pointer_map > 0 )); ...; % size_map( pointer_map > 0 );( pointer_map > 0 )
pointer_map( pointer_map > 0 ) ...
    = 1000 ./ numel_strel_linear_LUT_range( round( size_map( pointer_map > 0 ))) ...
   .* double( nonzeros( pointer_map )); ...; % size_map( pointer_map > 0 );( pointer_map > 0 )

% function clear_locations( location )
% 
%     strel = location + strel_linear_LUT_range{ round( size_map( location ))} ;
% 
% %     is_strel_in_map = strel >= 1 & strel <= cum_prod_image_dims( 3 );
%     
%     is_strel_in_map = strel >= 1                      ...
%                     & strel <= cum_prod_image_dims( 3 );
%                   
% 	strel( ~ is_strel_in_map ) = location ; % a safe choice no errors
% 
%           [ subscript_y, subscript_x, subscript_z ] ...
%                                                          = ind2sub( size_of_image, location );
%     
%     subscripts ...
%         = [ subscript_y, subscript_x, subscript_z ];
%     
%     absolute_subscripts = local_subscripts_range{ round( size_map( location ))} + subscripts ;
%     
%     is_strel_in_map  = logical( prod(   absolute_subscripts >= [ 1, 1, 1 ]      ...
%                                       & absolute_subscripts <= size_of_image, 2 ));
%     
%     strel( ~ is_strel_in_map  ) = [ ];
% 
%     are_neighbors_pointing_at_location ...
%         = strel_pointers_LUT_range{ round( size_map( location ))}( is_strel_in_map  ) ...
%         ==      pointer_map( strel );
%     
% %     are_neighbors_of_this_vertex        ...
% %         =  vertex_index_map( strel    ) ...
% %         == vertex_index_map( location );
% % 
% %     are_neighbors_of_no_vertex               ...
% %         =  vertex_index_map( strel    ) == 0 ;
% % 
% %     locations_to_clear = strel(       are_neighbors_pointing_at_location ...
% %                                 & (   are_neighbors_of_this_vertex       ...
% %                                     | are_neighbors_of_no_vertex          ))' ;
% %     locations_to_clear = strel(   are_neighbors_pointing_at_location )' ;
%     are_neighbors_the_same_size_as_this_vertex        ...
%         =  size_map( strel    ) ...
%         == size_map( location );
% 
%     locations_to_clear = strel(   are_neighbors_pointing_at_location            ...
%                                 & are_neighbors_the_same_size_as_this_vertex )' ; % pointer maps change based on the size of the strel
%   
%     vertex_index_map( location ) = 0 ;
%     available_locations_map( location ) = false ;
%          pointer_map( location ) = 0 ; 
% %           energy_map( location ) = 0 ;
%           
%           available_locations( available_locations == location ) = [ ];
% 
%     for location = locations_to_clear
% 
%         clear_locations( location )
% 
%     end % FOR location to clear
% 
% end % FUNCTION clear_locations

% % function [ direction_home, parent_edge_location ] = get_direction_home( location )
% function [ direction_home ] = get_direction_home( location )
% 
% %     if branch_order_map( location ) == 1 % home is a vertex
%                        
%         direction_home = index2position( vertex_locations( vertex_index_map( location ) - 1 ), cum_prod_image_dims )...
%                        - index2position(                                     location        , cum_prod_image_dims );
%                    
%         direction_home = double( direction_home' ) .* microns_per_voxel ; % conversion from voxels to microns 
%         
%         direction_home = direction_home ./ sum( direction_home .^ 2 ) .^ 0.5 ;
%         
% %         parent_edge_location = 0 ; % code for no parent edge
% %         
% %     else % home is a parent edge
% %     
% %         direction_home = [ 0, 0, 0 ];
% % 
% %         % trace back until hitting a zero in the pointer index map, recording the linear indices
% %         % into the image and setting the pointer index map to zero along the way.  If this
% %         % process hits a zero index that doesn't belong to a vertex then record this for
% %         % later (we will need to add a new vertex and break up the existing edge into two).
% %         while true
% % 
% %             if pointer_map( location ) == 0, break, end
% % 
% %             direction_home = direction_home + strel_unit_vectors_LUT_range( pointer_map( location ), : );
% % 
% %             location       =       location + strel_linear_LUT_range(       pointer_map( location )    );
% % 
% %         end % WHILE tracing
% %         
% %         parent_edge_location = location ;
% %         
% %     end % IF tracing will end at vertex
% end % FUNCTION find_direction_home

%     function edge_direction = get_parent_edge_direction( parent_edge_location )
%         
%         parent_unit_vectors_in_strel                                                                                      ...
%             = strel_unit_vectors_LUT(        pointer_map( parent_edge_location + strel_linear_LUT( : )) == 0              ...
%                                       & vertex_index_map( parent_edge_location + strel_linear_LUT( : )) > border_index, : );
%         
%         covariance_matrix = cov( parent_unit_vectors_in_strel ); %%% ??? scale by microns_per_voxel ???
%         
%         [ largest_principal_component, ~ ] = eigs( covariance_matrix, 1 );
%         
%         edge_direction = largest_principal_component' ;
%         
%     end % FUNCTION get_edge_direction

% max( moves_per_vertex )

end % FUNCTION get_edges_by_watershed


