function [ edges2vertices, lengths_of_edges, edge_space_subscripts, edge_scale_subscripts, edge_energies ] ...
                       = get_edges_V203( pixels_per_sigma_range, vertex_scale_subscripts,           ...
                                         vertex_space_subscripts, strel_apothem,     ...
                                         max_edge_length_per_radius, number_of_edges_per_vertex,    ...
                                         edge_walk_temperature, data_directory, energy_handle )
%% SAM 2/10/18
%
% This function takes the vertex positions, scales, and associated directions and looks for
% neighboring vertices to connect to and make edges.  If no vertices exist in the local
% neighborhood, it will place more vertices when possible and continue searching in a greater
% neighborhood.
%
% Edges will be created when the direction of one vertex is found to be within 1.5 sigma of a
% neighboring vertex. The creation of an edge will prevent this direction from placing new
% vertices. The method of placing new vertices is conceptually equivalent to the method of finding
% directions from phase 2, except that k will be forced to be 1, so that only one direction will be
% found, and the sampling element to find that direction will be a truncated cone (originating at
% the direction location with angle of thirty degree off the axis from owner vertex to direction
% location, and bounded in radial direction by 1 and 3 sigma of the owner vertex.
%
% V02 in which the inputs are changed to be ammenable to the vectorize_V02 script 3/5/18 SAM
%
% V10 in which the solution strategy is changed from resembling Phase II to more resembling Phase I.
% The directions serve as the starting points for walks toward the minimum laplacian.  The previous
% 26 neighbor strel is zero'd at every scale at the start of each new step. %SAM 4/5/18
%
% V14 in which the direction information from phase II is not used, and instead a probabilistic
% approach is taken to explore the most likely paths from starting vertex to terminal vertex by
% converting the laplacian image into probabilities using the Boltzmann distribution with a
% characteristic energy that makes the starting vertex voxel a characteristic likelihood when
% compared to its 26 neighbors (a min projection is taken across scale space).  A number of tokens
% are then placed at the starting vertex with a uniform distribution of percentiles between 0 and
% 100 assigned to them.  The possible next voxels are then ordered from most to least likely,
% forming bins in which the tokens will be placed according to their assigned percentiles.  If a
% token's percentile falls into the bin of a possible next voxel, then the token will move to that
% voxel.  The unconditional probability of landing on a certain voxel on a certain trajectory can
% then be calculated as the product of all of the bin widths of the options that led to that voxel.
% Voxels whose unconditional probabilities fall below 1 / number_of_tokens are removed from the pool
% of possible next voxels at each step.  If all bins are removed this way, then the token moves to
% the voxel with the least laplacian.  Token percentiles are adjusted to fall within 0 and 100 after
% each move by mapping the range of percentiles that would have led to that voxel onto (0,100). SAM
% 4/22/18
%
% note: we should decide the temperature for the whole field of view, then decide how many edges to
% give to each vertex to capture what is around it. Ask how closely do the edge percentiles need to
% be spaced to give at least one look at the, say, 20th percentile first option out of the vertex
%
% V130 in which the characteristic energy is the energy of the origin vertex times some scaling
% factor divided by the radius of the vertex in pixels (geomean of the x, y, and z pixel radii in
% anisotropic case). Also changing the word "laplacian" to "energy" in the variable names. SAM
% 5/3/18  Also the termination of trajectories is now probability based instead of number of voxels
% based. SAM 5/4/18
%
% V131 in which the characteristic energy is not divided by the radius of the vertex, instead that
% factor is applied to each minimum_probability_threshold on a per vertex basis. SAM 5/4/18
%
% V150 in which the enrgy file is min projected and the second index of the 4th dimension is the
% image of scale indices associated with each minimum.  % 5/7/18 SAM
%
% V160 in which the subscripts are returned instead of the "positions" SAM 5/14/18
%
% V190 significant improvements to speed and clarity were made SAM 8/15/18  Output is one step away
% from outputting space and scale subscrtips separately (to save memory as scale subscripts can be
% uint8 and space can be uint16
%
% V200 in which the energy h5 file is expecting the scale and energy images in the reverse order of
% the previous version.  Also, expecting cell vector edge_subscripts to be split into a cell vector
% of scale subscripts (uint8) and one of space subscripts (uint16).  SAM + WAS 12/5/18
%
% V201: replaced directories input with data_directory SAM 12/7/18
%
% V202: Included the edge_walk_temperature input parameter for use with the noise study.  SAM
% 12/17/18
%
% V203: removed the hard switch to deterministic minimum-searching.  Now the search is always
% probabilistic, up to rounding error.  SAM 1/7/19

%% get variable ranges

path_to_energy_data = [ data_directory, energy_handle ];

energy_file_info    = h5info( path_to_energy_data );

image_dims          = energy_file_info.Datasets.Dataspace.Size ;

strel_apothem = uint16( strel_apothem );

strel_width = 2 * strel_apothem + 1 ;

% % strel_dims  = uint8( strel_width * [ 1, 1, 1 ]);
% strel_dims  = [ strel_width, strel_width, strel_width ];

numel_of_strel = uint8( strel_width ^ 3 );

middle_of_strel = ( numel_of_strel + 1 ) / 2 ;

local_subscripts_range = 0 : strel_width - 1 ;

[ local_subscripts_mesh_y,                                     ...
  local_subscripts_mesh_x,                                     ...
  local_subscripts_mesh_z  ] = ndgrid( local_subscripts_range, ...
                                       local_subscripts_range, ...
                                       local_subscripts_range  );
                                  
local_subscripts_template = cat( 4, local_subscripts_mesh_y, ...
                                    local_subscripts_mesh_x, ...
                                    local_subscripts_mesh_z  );
                                
local_subscripts_template = permute( local_subscripts_template, [ 4, 1, 2, 3 ]);

local_subscripts_template = reshape( local_subscripts_template, 3, numel_of_strel );

local_subscripts_template = local_subscripts_template' ;

local_subscripts_range    =    local_subscripts_range' ;

average_radius_in_pixels_range = geomean( pixels_per_sigma_range, 2 );

max_edge_length_in_scales = min( max_edge_length_per_radius .* average_radius_in_pixels_range, 100 );

number_of_vertices  = length( vertex_scale_subscripts ); 

number_of_edges     = number_of_vertices * number_of_edges_per_vertex ;

lengths_of_edges      = zeros( number_of_edges, 1, 'uint16' ); % number of positions visited (units of voxels)
  edges2vertices      = zeros( number_of_edges, 2, 'uint32' ); % the indices of the vertices connected by each edge

edge_space_subscripts =  cell( number_of_edges, 1 );

edge_index_range      =    1 : number_of_edges ;

past_subscript_memory_limit = 100 ; % number of entries in the list of previous locations

% deterministic_probability_threshold = 1 / number_of_edges_per_vertex ;

% load the energy data and associated scale index data (two 3D images concatenated into 4D array )
energy_and_index_data = h52mat( path_to_energy_data );

energy_data = energy_and_index_data( :, :, :, 2 );

index_data  = energy_and_index_data( :, :, :, 1 );

clear( 'energy_and_index_data' )

% set characeteristic energy globally by looking at the energy fluctuations in the noise of the energy image
% characteristic_energy = mean( energy_data( energy_data( : ) < 0 )) * characteristic_energy_fraction ;

% infinity out the borders of the energy image so that the strel will never reach the end of the
% image
energy_data(       [ 1 : strel_apothem, image_dims( 1 ) - strel_apothem + 1 : image_dims( 1 )], :, : ) = Inf ;
energy_data(    :, [ 1 : strel_apothem, image_dims( 2 ) - strel_apothem + 1 : image_dims( 2 )], :    ) = Inf ;
energy_data( :, :, [ 1 : strel_apothem, image_dims( 3 ) - strel_apothem + 1 : image_dims( 3 )]       ) = Inf ;

% strel_apothem = uint16( strel_apothem ); SAM + WAS 12/5/18


%% main PARFOR: loop through potential edges
parfor edge_index = edge_index_range
        
    vertex_index        =  floor(( edge_index - 1 ) / number_of_edges_per_vertex ) + 1 ;    
    
    percentile_index    = edge_index - ( vertex_index - 1 ) * number_of_edges_per_vertex ;
    
    edge_percentile     = percentile_index / ( number_of_edges_per_vertex + 1 );
                         
    % initialize the neighbor vertex at zero.  zero in the output means no terminal vertex was found
    terminal_vertex_index = uint32( 0 );
    
    edge_is_below_length_maximum   = true ;
    searching_for_terminal_vertex  = true ;
    
    not_first_edge_index = false ;    
    
    current_space_subscripts = vertex_space_subscripts( vertex_index, : );
%     current_scale_subscripts = vertex_scale_subscripts( vertex_index    );
                
    current_subscript_ordinate = uint8( 1 );
    
%     characteristic_energy = energy_values( vertex_index ) * characteristic_energy_fraction ;

    past_subscripts = zeros( past_subscript_memory_limit + numel_of_strel, 3, 'uint16' );
    
    % put the current vertex location on the past_subscripts list so it won't be chosen again
    past_subscripts( 1, : ) = current_space_subscripts ;       
        
%     current_probability              = 1 ;
%     current_probability_renormalized = 1 ;    
    
    max_edge_length_at_vertex = max_edge_length_in_scales( round( vertex_scale_subscripts( vertex_index )));
    
%     % adjust the probability minimum by the number of pixels per sigma of the origin vertex
%     probability_minimum_at_vertex = probability_minima_in_scales( current_scale_subscripts );
    
    % note: larger vertices should get more trajectories because they intrinsically have more
    % options, because they are at a higher resolution. SAM 5/4/18
            
    %% growing the edges WHILE
    while edge_is_below_length_maximum && searching_for_terminal_vertex
        
        current_local_subscripts_range = current_space_subscripts + local_subscripts_range - strel_apothem ;
        
        energy_data_chunk = energy_data( current_local_subscripts_range( :, 1 ), ...
                                         current_local_subscripts_range( :, 2 ), ...
                                         current_local_subscripts_range( :, 3 )  );
                                                            
        past_subscripts_local_temp = past_subscripts + strel_apothem + 1 - current_space_subscripts ;

        past_subscripts_local                                                                   ...
            = past_subscripts_local_temp(   past_subscripts_local_temp( :, 1 ) >   0            ...
                                          & past_subscripts_local_temp( :, 2 ) >   0            ...
                                          & past_subscripts_local_temp( :, 3 ) >   0            ...
                                          & past_subscripts_local_temp( :, 1 ) <=  strel_width  ...
                                          & past_subscripts_local_temp( :, 2 ) <=  strel_width  ...
                                          & past_subscripts_local_temp( :, 3 ) <=  strel_width, ...
                                                                    :                           );

        current_energy = energy_data_chunk( middle_of_strel );
                                                                
        % Infinity out previous positions
        energy_data_chunk(     past_subscripts_local( :, 1 )                                 ...
                           + ( past_subscripts_local( :, 2 ) - 1 ) * strel_width             ...
                           + ( past_subscripts_local( :, 3 ) - 1 ) * strel_width ^ 2 ) = Inf ;

        valid_local_indices = find( energy_data_chunk( : ) < 0 );
        
        % IF there exist candidate next positions (requiring negative energy here)
        if ~ isempty( valid_local_indices )
            

            % converting energy to probability by Boltzmann distribution (user can change temperature),
            % characteristic energy variable carries the negative sign that is missing in below
            
            exponential_argument = double( lengths_of_edges( edge_index ) + 1 ) * energy_data_chunk( : ) / ( current_energy * edge_walk_temperature );
            
            if max( exponential_argument ) == Inf % IF edge_walk_temperature is zero (up to rounding error)
                
                [ ~, next_local_index ] = min( energy_data_chunk( : ));
                
            else
            
                exponential_argument = exponential_argument - max( exponential_argument );

                probability_chunk = exp( exponential_argument );

%                 probability_chunk =      probability_chunk ...
%                                   / sum( probability_chunk );       
                              
                probabilities_valid_renormalized =      probability_chunk ...
                                                 / sum( probability_chunk );

                cumulative_probabilies = [ 0; cumsum( probabilities_valid_renormalized )];

                next_local_index = find( edge_percentile < cumulative_probabilies( 2 : end ), 1 );

% %                 % calculate the current probability (raw and renormalized to only include the above threshold options)
% %                 current_probability              = current_probability                 ...
% %                                                  * probability_chunk( next_local_index );
% 
%                 current_probability_renormalized = current_probability_renormalized                 ...
%                                                  * probabilities_valid_renormalized( next_local_index );

                % renormalize the edge_percentile into the (0,1) range (in this variable we record
                % the conditional probability of future moves given the previous moves)
                edge_percentile                                                                                      ...
                    = ( edge_percentile                                - cumulative_probabilies( next_local_index )) ...
                    / ( cumulative_probabilies( next_local_index + 1 ) - cumulative_probabilies( next_local_index )) ;                 

            end
            
%             invalid_probability_indices =  current_probability_renormalized   ...
%                                         *  probability_chunk                  ...
%                                         < deterministic_probability_threshold ;
% 
%             % zero out probabilities below the threshold
%             probability_chunk( invalid_probability_indices ) = 0 ;

%             % if there exists at least one probabilistic option
%             if any( probability_chunk )

                % renormalize without the options of below threshold probabilities

%             else % there are no probabilistic options


%             end % if there exists at least one probabilistic option
%             end  % IF no options are probable enough to warrant non deterministic action
                
            next_space_subscripts = current_space_subscripts + local_subscripts_template( next_local_index, : ) - strel_apothem ;
            
        else % ELSE there is no candidate for the next position

            next_space_subscripts = [ ];

            searching_for_terminal_vertex = false ;

        end % IF there exist candidate next positions

        % IF not the first edge index
        if not_first_edge_index

            % check to see if the current step is a vertex or not (in a cascade of checks)
            terminal_vertex_index_candidates_1 ...
                = find( vertex_space_subscripts( :, 1                                  ) == current_space_subscripts( 1 )    );

            % IF vertices found
            if ~ isempty( terminal_vertex_index_candidates_1 )

                terminal_vertex_index_candidates_2 ...
                    = find( vertex_space_subscripts( terminal_vertex_index_candidates_1, 2 ) == current_space_subscripts( 2 )    );

                % IF vertices found                
                if ~ isempty( terminal_vertex_index_candidates_2 )

                    terminal_vertex_index_candidates_3 ...
                        = find( vertex_space_subscripts( terminal_vertex_index_candidates_1(                                    ...
                                               terminal_vertex_index_candidates_2  ), 3 ) == current_space_subscripts( 3 ), 1 );

                    % IF vertex found
                    if ~ isempty( terminal_vertex_index_candidates_3 )

                        % record the terminal vertex                        
                        candidate_terminal_vertex_index                                 ...
                                = uint32( terminal_vertex_index_candidates_1(           ...
                                              terminal_vertex_index_candidates_2(       ...
                                                  terminal_vertex_index_candidates_3 )));   

                        terminal_vertex_index = candidate_terminal_vertex_index ;

                        searching_for_terminal_vertex = false ;

                    end % IF vertex   found
                end     % IF vertices found
            end         % IF vertices found         
        end % IF not the first edge index
        
        not_first_edge_index = true ;
        
        % recording the 5D ( 3-spatial, size, laplacian ) position of this edge at this growth step
        edge_space_subscripts{ edge_index }          ...
            = [ edge_space_subscripts{ edge_index }; ...
                current_space_subscripts;            ];
                                
        % end of the loop housekeeping
        lengths_of_edges( edge_index )           ...
            = lengths_of_edges( edge_index ) + 1 ;
                             
        % check if current probability is below the hard minimum
%         edge_is_above_probability_minimum = current_probability > probability_minimum ;
        
        edge_is_below_length_maximum = lengths_of_edges( edge_index ) < max_edge_length_at_vertex ;

        additional_past_subscripts = local_subscripts_template( valid_local_indices, : ) ...
                                   + current_space_subscripts                            ...
                                   - strel_apothem ;

        numel_of_valid_local_voxels = numel( valid_local_indices );
        
        next_subscript_ordinate = current_subscript_ordinate + numel_of_valid_local_voxels ;
        
        past_subscripts(   current_subscript_ordinate + 1 : next_subscript_ordinate, : )            ...
                                                                       = additional_past_subscripts ;

        current_subscript_ordinate = next_subscript_ordinate ;
        % If we go over the past subscript memory limit, then overwrite the previous subscripts that
        % are farthest away from the current (at the start of the list)
        if current_subscript_ordinate + numel_of_valid_local_voxels > past_subscript_memory_limit
            
            current_subscript_ordinate = 0 ;
            
        end

        current_space_subscripts  = next_space_subscripts ;
%         current_energy            = next_energy_value ;

    end % WHILE growing    
    
    % mark the owner vertex as the start of the edge and the neighbor vertex as the end
    edges2vertices( edge_index, : ) = uint32([ vertex_index, terminal_vertex_index ]);
        
end % PARFOR edge_index

edge_scale_subscripts                                                                                                           ...
    = cellfun( @( space_subscripts )                                                                                            ...
                                 index_data(   double(  space_subscripts( :, 1 ))                                               ...
                                             + double(( space_subscripts( :, 2 ) - 1 )) * image_dims( 1 )                       ...
                                             + double(( space_subscripts( :, 3 ) - 1 )) * image_dims( 1 ) * image_dims( 2 )),   ...
               edge_space_subscripts, 'UniformOutput', false                                                                 ); ...
               
edge_energies                                                                                                                   ...
    = cellfun( @( space_subscripts )                                                                                            ...
                                energy_data(   double(  space_subscripts( :, 1 ))                                               ...
                                             + double(( space_subscripts( :, 2 ) - 1 )) * image_dims( 1 )                       ...
                                             + double(( space_subscripts( :, 3 ) - 1 )) * image_dims( 1 ) * image_dims( 2 )),   ...
               edge_space_subscripts, 'UniformOutput', false                                                                 ); ...

end % FUNCTION