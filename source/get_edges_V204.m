function [ edges2vertices, edge_lengths, edge_space_subscripts, edge_scale_subscripts, edge_energies ]  ...
                     = get_edges_V204(                lumen_radius_in_microns_range, microns_per_voxel, ...
                               length_dilation_ratio, vertex_scale_subscripts, vertex_space_subscripts, ...
                                                      strel_apothem, max_edge_length_per_origin_radius, ...
                                                            number_of_edges_per_vertex, data_directory, ...
                                                                                         energy_handle  )
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
%
% V204: abandoning the random walk for a deterministic search for the paths with the least mean energy
% from the origin vertex. SAM 4/18/19

%% get variable ranges

path_to_energy_data = [ data_directory, energy_handle ];

energy_file_info    = h5info( path_to_energy_data );

image_pair_dims     = energy_file_info.Datasets.Dataspace.Size ;

image_dims          = image_pair_dims( 1 : 3 );

cum_prod_image_dims = cumprod( image_dims );

% numel_image         = prod( image_dims( 1 : 3 ));
% numel_image_yx      = prod( image_dims( 1 : 2 ));

% strel_apothem = uint16( strel_apothem );

strel_width   = 2 * strel_apothem + 1 ;

% estimate of number of voxels per edge location:
% numel_of_strel = strel_width ^ 3 ;
numel_of_strel = round( strel_width ^ 3 / 2.5 );

% % strel_dims  = uint8( strel_width * [ 1, 1, 1 ]);
% strel_dims  = [ strel_width, strel_width, strel_width ];
% middle_of_strel = ( numel_of_strel + 1 ) / 2 ;

% % underestimate of number of voxels per edge location:
% numel_of_strel_cross_section = strel_width ^ 2 ;

local_subscripts_range = - strel_apothem : strel_apothem ;

linear_strel = local_subscripts_range                                                ;
linear_strel = local_subscripts_range * cum_prod_image_dims( 1 ) + linear_strel( : ) ;
linear_strel = local_subscripts_range * cum_prod_image_dims( 2 ) + linear_strel( : ) ;
linear_strel =                                              int64( linear_strel( : ));

[ strel_y,                                     ...
  strel_x,                                     ...
  strel_z  ] = ndgrid( local_subscripts_range, ...
                       local_subscripts_range, ...
                       local_subscripts_range  );

strel_yxz = [ strel_y( : ), strel_x( : ), strel_z( : )];

strel_distance_yxz = strel_yxz .* microns_per_voxel ;

strel_distance_LUT = sum( strel_distance_yxz .^ 2, 2 ) .^ 0.5 ;

% average_radius_in_pixels_range = geomean( lumen_radius_in_pixels_range, 2 );
% 
% max_edge_length_in_scales = max_edge_length_per_radius .* average_radius_in_pixels_range ;

max_edge_length_in_microns_range = max_edge_length_per_origin_radius .* lumen_radius_in_microns_range ;

number_of_vertices  = length( vertex_scale_subscripts ); 

% number_of_edges     = number_of_vertices * number_of_edges_per_vertex ;

% lengths_of_edges      = cell( number_of_vertices, 1 );
edges2vertices        = cell( number_of_vertices, 1 );
edge_indices_temp     = cell( number_of_vertices, 1 );

%      lengths_of_edges{ : } = zeros( number_of_edges_per_vertex, 1, 'uint16' ); % number of positions visited (units of voxels)
%        edges2vertices{ : } = zeros( number_of_edges_per_vertex, 2, 'uint32' ); % the indices of the vertices connected by each edge
% edge_indices{ : } = zeros( number_of_edges_per_vertex, vertex_length_maximum, 'uint16' );

% edge_indices{ : } = cell( number_of_edges_per_vertex, 1 );

vertex_unique_range = uint32( 1 : number_of_vertices );

% load the energy data and associated scale index data (two 3D images concatenated into 4D array )
energy_and_index_data = h52mat( path_to_energy_data );

energy_data = energy_and_index_data( :, :, :, 2 );

 index_data = energy_and_index_data( :, :, :, 1 );

clear( 'energy_and_index_data' )

% set characeteristic energy globally by looking at the energy fluctuations in the noise of the energy image
% characteristic_energy = mean( energy_data( energy_data( : ) < 0 )) * characteristic_energy_fraction ;

% zero out the borders of the energy image so that the strel will never reach the end of the image
energy_data(       [ 1 : strel_apothem, image_dims( 1 ) - strel_apothem + 1 : image_dims( 1 )], :, : ) = 0 ;
energy_data(    :, [ 1 : strel_apothem, image_dims( 2 ) - strel_apothem + 1 : image_dims( 2 )], :    ) = 0 ;
energy_data( :, :, [ 1 : strel_apothem, image_dims( 3 ) - strel_apothem + 1 : image_dims( 3 )]       ) = 0 ;

% strel_apothem = uint16( strel_apothem ); SAM + WAS 12/5/18

vertex_space_subscripts = int64( vertex_space_subscripts );

vertex_locations =   vertex_space_subscripts( :, 1 )                                  ...
                 + ( vertex_space_subscripts( :, 2 ) - 1 ) * cum_prod_image_dims( 1 ) ...
                 + ( vertex_space_subscripts( :, 3 ) - 1 ) * cum_prod_image_dims( 2 );
                  
% clear( 'vertex_space_subscripts' )

vertex_space_subscripts = double( vertex_space_subscripts );

vertex_image = sparse( double( vertex_locations ), 1, double( vertex_unique_range ), cum_prod_image_dims( 3 ), 1, number_of_vertices );

% % spherical structuring element (strel) templates
% number_of_scales = size( lumen_radius_in_microns_range, 1 );
% 
% lumen_radius_in_voxels_range = lumen_radius_in_microns_range ./ microns_per_voxel ;
% 
% dilated_lumen_radius_in_voxels_range = lumen_radius_in_voxels_range * length_dilation_ratio ;
% 
% strel_templates = cell( number_of_scales, 1 );
% 
% for scale_index = 1 : number_of_scales
% 
%     % find all pixel locations within the ellipsoid radii from the vertex position
%     
%     strel_templates{ scale_index }                                                                  ...
%         = construct_structuring_element_V190( dilated_lumen_radius_in_voxels_range( scale_index, : ), ...
%                                               image_dims                                              );
%         
% end % scale FOR


% % find all pixel locations within the ellipsoid radii from the vertex position
% 
% % only the largest subscript template is saved
% [ strel_template, subscript_template ]                                                                 ...
%     = construct_structuring_element_V190( dilated_lumen_radius_in_voxels_range( number_of_scales, : ), ...
%                                           image_dims                                                   );
        

% numel_of_origin_vertex_strel = round( 4 / 3 * pi * lumen_radius_in_microns_range( end ) .^ 3 / prod( microns_per_voxel ));

% vertex_index_list = [ 1353 ];

%% main PARFOR: loop through the vertices, searching for nearby vertices to connect to
% changed from PARFOR to FOR 6/8/20
for vertex_index = vertex_unique_range
               
%     if any( vertex_index == vertex_index_list )
%        
%         disp('here')
%         
%     end
    
    % initialize the terminal vertex at zero so that the original vertex doesn't end the search
    terminal_vertex_index = 0 ;
    
    current_scale_subscript      = vertex_scale_subscripts( vertex_index );    
    
	max_edge_length_in_microns   = exp( interp1( log( max_edge_length_in_microns_range ), current_scale_subscript ));
    
    max_edge_length_in_voxels    = round( max( max_edge_length_in_microns ./ microns_per_voxel )) + 1 ;
        
%     numel_of_origin_vertex_strel = numel( strel_templates( round( current_scale_subscript )));
    
%     numel_of_origin_vertex_strel = round( 4 / 3 * pi * lumen_radius_in_microns_range( round( current_scale_subscript )) .^ 3 / prod( microns_per_voxel ));
        
    max_number_of_indices        = max_edge_length_in_voxels * number_of_edges_per_vertex ;
        
    previous_indices_visited     = zeros( max_number_of_indices, 1, 'int64' );
    
    max_nnz_of_sparse_arrays     = max_number_of_indices * numel_of_strel ;    
    
    edge_indices_temp{ vertex_index } = zeros( number_of_edges_per_vertex, max_edge_length_in_voxels,  'int64' );
       edges2vertices{ vertex_index } = zeros( number_of_edges_per_vertex, 2,                         'uint32' ); % the indices of the vertices connected by each edge
                   
       displacement_vectors               = zeros( 3, number_of_edges_per_vertex );
%        displacement_vector_sqr_magnitudes = zeros( 1, number_of_edges_per_vertex );
                            
       pointer_index_map         = spalloc( cum_prod_image_dims( 3 ), 1, max_nnz_of_sparse_arrays );
    available_energy_map         = spalloc( cum_prod_image_dims( 3 ), 1, max_nnz_of_sparse_arrays );
%             distance_map         = spalloc( cum_prod_image_dims( 3 ), 1, max_nnz_of_sparse_arrays + 4 * numel_of_origin_vertex_strel );        
            distance_map         = spalloc( cum_prod_image_dims( 3 ), 1, max_nnz_of_sparse_arrays );        
    
    number_of_indices = 0 ;
    
    number_of_edges_found = 0 ;
    
    current_location = double( vertex_locations( vertex_index ));

    pos_xy  = rem( current_location - 1, cum_prod_image_dims( 2 )) + 1 ;

    pos_z   = ( current_location - pos_xy ) / cum_prod_image_dims( 2 ) + 1 ;

    pos_y   = rem( pos_xy - 1, cum_prod_image_dims( 1 )) + 1 ;

    pos_x   = ( pos_xy - pos_y ) / cum_prod_image_dims( 1 ) + 1 ;            

    origin_position = [ pos_y; pos_x; pos_z ] .* microns_per_voxel' ;
    
    distance_map( current_location ) = 1 ;
    
    is_below_edge_max              = true ;
    is_vertex_below_length_maximum = true ;
    there_exists_possible_move     = true ;    
    
    % don't enter the loop if the current vertex is in a location that can't be accessed by an edge
    if ~ energy_data( current_location ), there_exists_possible_move = false ; end
        
    % WHILE searching around in new locations
    while is_vertex_below_length_maximum && there_exists_possible_move && is_below_edge_max
        
        % record the newly visited location and add a neighborhood surrounding (size/shape of the
        % strel) to the available energies to choose from.  Also keep track of the distances of the
        % newly revealed locations and mark them with a pointer that points those locations back to
        % this one.  Also need to check if the added location belongs to a vertex (distinct from the
        % origin vertex).
        number_of_indices               = number_of_indices + 1 ;
        
        is_vertex_below_length_maximum  = number_of_indices < max_number_of_indices ;

        previous_indices_visited( number_of_indices )   = current_location ;

        current_linear_strel                            = current_location + linear_strel ;
        
        is_new_index_in_strel   = ~ distance_map( current_linear_strel );

        new_indices_considered  = current_linear_strel( is_new_index_in_strel );
                        
        distance_map( new_indices_considered )                                                      ...
                        = distance_map( current_location ) + strel_distance_LUT( is_new_index_in_strel );
        
        pointer_index_map( new_indices_considered ) = number_of_indices ;
   
        new_indices_considered = new_indices_considered(   distance_map( new_indices_considered ) ...
                                                         < max_edge_length_in_microns             );
                                                     
        % IF any edges have been found
        if any( edges2vertices{ vertex_index }( :, 1 ))
            
            length_of_new_indices_considered = length( new_indices_considered );
            
            if length_of_new_indices_considered

                indices_beyond_found_vertex = zeros( 1, length_of_new_indices_considered );

                % convert the linear indexing to spatial subscripts
                pos_xy  = rem( double( new_indices_considered' ) - 1, cum_prod_image_dims( 2 )) + 1 ;
                
                pos_z   = ( double( new_indices_considered' ) - pos_xy ) / cum_prod_image_dims( 2 ) + 1 ;
                
                pos_y   = rem( pos_xy - 1, cum_prod_image_dims( 1 )) + 1 ;
                
                pos_x   = ( pos_xy - pos_y ) / cum_prod_image_dims( 1 ) + 1 ;

                vectors_from_origin_to_indices = [ pos_y; pos_x; pos_z ] .* microns_per_voxel' - origin_position ;

                for terminal_index = 1 : number_of_edges_found
                    
                    terminal_vertex_displacements = displacement_vectors( :, terminal_index );

                    indices_beyond_this_found_vertex = sum( terminal_vertex_displacements .* vectors_from_origin_to_indices, 1 ) >= 1 ;

                    indices_beyond_found_vertex = indices_beyond_found_vertex | indices_beyond_this_found_vertex ;

                end

                new_indices_considered( indices_beyond_found_vertex ) = [ ];
                
            end

        end

        available_energy_map( new_indices_considered ) = energy_data( new_indices_considered );
                        
        % IF terminal vertex found
        if terminal_vertex_index
                        
            number_of_edges_found = number_of_edges_found + 1 ;
            
%             % zero out a sphere representing this vertex in the available energy map
%             terminal_vertex_linear_strel = min( max( current_location + strel_templates{ round( vertex_scale_subscripts( terminal_vertex_index ))}, 1 ), cum_prod_image_dims( 3 ));
            
%             % zero out a hemisphere located at this vertex and facing the origin vertex in the available energy map
%             r_from_A_to_B = microns_per_voxel .* ( vertex_space_subscripts( terminal_vertex_index, : ) - vertex_space_subscripts( vertex_index, : ));
%             
%             distance_from_A_to_B = sum( r_from_A_to_B .^ 2 ) ^ 0.5 ;
%             
%             r_unit_from_A_to_B = r_from_A_to_B / distance_from_A_to_B ;
%             
%             distance_from_A_to_B_mesh = sum( microns_per_voxel .* subscript_template .* r_unit_from_A_to_B, 2 );
%             
%             hemisphere_template = strel_template( distance_from_A_to_B_mesh >= 0 & distance_from_A_to_B_mesh <= 2 * sum(( r_unit_from_A_to_B .* microns_per_voxel ) .^ 2 ) .^ 0.5 );
%             
%             terminal_vertex_linear_strel = min( max( current_location + hemisphere_template, 1 ), cum_prod_image_dims( 3 ));

%                     distance_map( terminal_vertex_linear_strel ) = Inf ;
% 
%             available_energy_map( terminal_vertex_linear_strel ) =   0 ;
            
                    distance_map( current_location ) = Inf ;

            available_energy_map( current_location ) =   0 ;

%             % Uncomment to only block out the current strel, not a spherical strel for the vertex:
%                     distance_map( current_linaer_strel ) = Inf ;
%             
%             available_energy_map( current_linaer_strel ) = 0 ;
                        
            tracing_index = current_location ;
            
            tracing_ordinate = 0 ;
            
            % trace back until hitting a zero in the pointer index map, recording the linear indices
            % into the image and setting the pointer index map to zero along the way.  If this
            % process hits a zero location that doesn't belong to a vertex then record this for
            % later (we will need to add a new vertex and break up the existing edge into two).
            while pointer_index_map( tracing_index ) > 0
                
                tracing_ordinate = tracing_ordinate + 1 ;
                
                edge_indices_temp{ vertex_index }( number_of_edges_found, tracing_ordinate ) = tracing_index ;
            
                tracing_index = previous_indices_visited( pointer_index_map( tracing_index ));
                
            end % WHILE tracing

            % Erase the pointers for the found trajectory, so that this path will not contribute to
            % any other edge from this vertex. Replace the pointers with the negative of the parent
            % edge identity that owns this path.
            pointer_index_map( nonzeros( edge_indices_temp{ vertex_index }( number_of_edges_found, : ))) = - number_of_edges_found ;

            % include the terminal location as the last entry in the list
            edge_indices_temp{ vertex_index }( number_of_edges_found, tracing_ordinate + 1 ) = tracing_index ;
                                    
            % Checking if there exists a parent edge to this one with worse energy (flagging it for
            % removal, because we can't support children with better energy than their parents.  (A
            % child requires a parent to make a physical connection between two vertices.))
            parent_index = - pointer_index_map( tracing_index );
            
            % IF edge has a parent
            if parent_index
            
                % check if this parent edge is itself a child
                parent_pointers = unique( - pointer_index_map( nonzeros( edge_indices_temp{ vertex_index }( parent_index, : ))));
                
                parent_pointers( parent_pointers == 0 | parent_pointers == parent_index ) = 0 ;
                
                % parent is also child
                if any( parent_pointers )
                    
%                     is_edge_valid = false ;

                    terminal_vertex_index = 0 ; origin_vertex_index = 0 ;                    
                    
                else % ELSE edge is valid
                
                    parent_indices = nonzeros( edge_indices_temp{ vertex_index }(          parent_index, : ));
                     child_indices = nonzeros( edge_indices_temp{ vertex_index }( number_of_edges_found, : ));
                    
                    parent_energies = energy_data( parent_indices );
                     child_energies = energy_data(  child_indices );

%                     parent_energy = get_edge_metric({ parent_energies });
%                      child_energy = get_edge_metric({  child_energies });
                     
                    parent_energy = max( parent_energies );
                     child_energy = max(  child_energies );

%                     % IF child has higher (worse) energy than parent, AND parent was valid.
%                     is_edge_valid = child_energy > parent_energy                    ...
%                                   & edges2vertices{ vertex_index }( parent_index, 1 );

                    % IF child has higher (worse) energy than parent, AND parent was valid.
                    if child_energy > parent_energy && edges2vertices{ vertex_index }( parent_index, 1 )

%                         % set the origin vertex to the terminal vertex of the parent
%                         origin_vertex_index = edges2vertices{ vertex_index }( parent_index, 1 )

%                         origin_vertex_index = vertex_index ;

                        % set the origin vertex to the vertex that has the lowest energy connection
                        % to the bifurcation point where the child split off
                        bifurcation_index = find( parent_indices == child_indices( end ));

                        parent_1_energies = energy_data( parent_indices( 1 : bifurcation_index - 1       ));
                        parent_2_energies = energy_data( parent_indices(     bifurcation_index + 1 : end ));
                        
                        % in case the child stems directly from the terminal vertex of the parent
                        if isempty( parent_1_energies )
                            
                            parent_1_energy = -Inf ;
                            
                        else
                            
%                             parent_1_energy = get_edge_metric({ parent_1_energies });
                            
                            parent_1_energy = max( parent_1_energies );                            
                            
                        end
                        
                        parent_2_energy = max( parent_2_energies );                            

                        [ ~, better_half_of_parent ] = min([ parent_1_energy, parent_2_energy ]);
                        
%                         % sort by length
%                         parent_1_energy = length( parent_1_energies ) - 0.5 * ( 1 == better_half_of_parent );
%                         parent_2_energy = length( parent_2_energies ) - 0.5 * ( 2 == better_half_of_parent );
%                         
%                         [ ~, better_half_of_parent ] = min([ parent_1_energy, parent_2_energy ]);                        

                        origin_vertex_index = edges2vertices{ vertex_index }( parent_index, better_half_of_parent );

                    else % ELSE edge not valid
                                        
                        terminal_vertex_index = 0 ; origin_vertex_index = 0 ;                    
                    
                    end % ELSE edge not valid
                end % ELSE edge is valid
                 
            else % ELSE edge has no parent
                
%                 is_edge_valid = true ;

                origin_vertex_index = vertex_index ;

            end % IF edge has parent
             
%             if is_edge_valid

            % Record the identity of the origin vertex in the second position in edges2vertices
            % (corresponds to last entry in the list of positions).  Also record the terminal
            % vertex, first position.
            edges2vertices{ vertex_index }( number_of_edges_found, : ) = [ terminal_vertex_index, origin_vertex_index ];
            
            % convert the linear indexing to spatial subscripts                        
            pos_xy  = rem( current_location - 1, cum_prod_image_dims( 2 )) + 1 ;

            pos_z   = ( current_location - pos_xy ) / cum_prod_image_dims( 2 ) + 1 ;

            pos_y   = rem( pos_xy - 1, cum_prod_image_dims( 1 )) + 1 ;

            pos_x   = ( pos_xy - pos_y ) / cum_prod_image_dims( 1 ) + 1 ;            
            
            displacement_vectors( :, number_of_edges_found ) = [ pos_y; pos_x; pos_z ] .* microns_per_voxel' - origin_position ;
            
            % divide by its magnitude squared so that later projections can be compared to unity
            displacement_vectors( :, number_of_edges_found ) = displacement_vectors( :, number_of_edges_found ) / sum( displacement_vectors( :, number_of_edges_found ) .^ 2 );

            % If we find this edge not valid by the above rule, we will place a zero in the
            % edges2vertices variable for this edge's entry.  Thus, a zero in the
            % edges2vertices variable is code for "no edge was found"            
                
%             end % IF valid edge
            
        else % ELSE no vertex exists at the current location
            
            available_energy_map( current_location ) = 0 ;

        end % IF terminal vertex found
        
        is_below_edge_max = number_of_edges_found < number_of_edges_per_vertex ;

        % choose next location as the next available with lowest energy
        [ min_available_energy, current_location ] = min( available_energy_map );
                
        terminal_vertex_index = full( vertex_image( current_location ));
                
        % requiring negative edge energy
        there_exists_possible_move = min_available_energy < 0 ;
                            
    end % WHILE searching for edges
    
end % PARFOR vertex_index

% % looping back through in series looking to add new vertices and edges where needed
% for vertex_index = vertex_index_range
%     
% end % FOR vertex

%% Organizing outputs:

% extracting outputs from the edges2vertices and edge_indices cell arrays (whose format was
% temporarilly different to allow for the parrallel FOR to run):
edges2vertices = cell2mat( edges2vertices );

% looking for nonzero in the first position (means an edge was found for that place holder)
found_edge_index_range = find( edges2vertices( :, 1 ))';

% removing placeholders that don't hold edges
edges2vertices = edges2vertices( found_edge_index_range, : );

number_of_edges = numel( found_edge_index_range );

edge_indices = cell( number_of_edges, 1 );

vertex_indices = floor(( found_edge_index_range - 1 ) / number_of_edges_per_vertex ) + 1 ;

edge_indices_at_vertices = found_edge_index_range - ( vertex_indices - 1 ) * number_of_edges_per_vertex ;

edge_index_range = 1 : number_of_edges ;

for edge_index = edge_index_range
    
    edge_indices{ edge_index } = double( nonzeros( edge_indices_temp{ vertex_indices( edge_index )}( edge_indices_at_vertices( edge_index ), : )));
    
end % FOR found edge

clear( 'edge_indices_temp' )

edge_lengths = cellfun( @( x ) uint16( numel( x )), edge_indices ); % number of positions visited (units of voxels)

% convert the linear indexing to spatial subscripts

edge_space_subscripts_xy  = cellfun( @( t )  rem( t - 1 ,    cum_prod_image_dims( 2 )) + 1, edge_indices,                                      'UniformOutput', false );

edge_space_subscripts_z   = cellfun( @( t, xy ) ( t - xy ) / cum_prod_image_dims( 2 )  + 1, edge_indices, edge_space_subscripts_xy,            'UniformOutput', false );

edge_space_subscripts_y   = cellfun( @( xy ) rem( xy - 1,    cum_prod_image_dims( 1 )) + 1, edge_space_subscripts_xy,                          'UniformOutput', false );

edge_space_subscripts_x   = cellfun( @( xy, y ) ( xy - y ) / cum_prod_image_dims( 1 )  + 1, edge_space_subscripts_xy, edge_space_subscripts_y, 'UniformOutput', false );

clear( 'edge_space_subscripts_xy' )

edge_space_subscripts     = cellfun( @( y, x, z ) uint16([ y, x, z ]), edge_space_subscripts_y, edge_space_subscripts_x, edge_space_subscripts_z, 'UniformOutput', false );

clear( 'edge_space_subscripts_y', 'edge_space_subscripts_x', 'edge_space_subscripts_z' )

edge_scale_subscripts     = cellfun( @( x )  index_data( x ), edge_indices, 'UniformOutput', false );
               
edge_energies             = cellfun( @( x ) energy_data( x ), edge_indices, 'UniformOutput', false );

end % FUNCTION