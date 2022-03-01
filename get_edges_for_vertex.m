function [ edge_indices_temp, edges2vertices] ...
                = get_edges_for_vertex( strel_linear_indexing_templates, vertex_volume_image, energy_image, size_image, microns_per_voxel, strel_apothem, ...
                                      number_of_edges_per_vertex, use_case, current_index, ...
                       vertex_image, vertex_indices, reading_box_apothem, linear_strel, ...
                                        max_edge_length_in_microns_range, numel_of_strel, strel_distance_LUT, ...
                                 cum_prod_image_dims, size_of_image, path_to_energy_data  )

%     if any( vertex_index == vertex_index_list )
%        
%         disp('here')
%         
%     end

vertex_index = full( vertex_image( current_index )); % vertex_index == 0 if the use_case is add_vertex_to_edges

number_of_vertices = length( vertex_indices );
% number_of_vertices = max( vertex_image );

% initialize the terminal vertex at zero so that the original vertex doesn't end the search
terminal_vertex_index = 0 ;

current_vertex_position = double( index2position( current_index, cum_prod_image_dims ));

% reference h5 for current scale
current_scale_subscript = h52mat( path_to_energy_data, [ current_vertex_position; 1 ], [ 1; 1; 1; 1 ]);

% current_scale_subscript      = vertex_scale_subscripts( vertex_index );    

% if max_edge_length_in_microns == Inf
%     
%     max_edge_length_in_microns   = exp( interp1( log( max_edge_length_in_microns_range ), current_scale_subscript ));
%     
% end

max_edge_length_in_microns   = exp( interp1( log( max_edge_length_in_microns_range ), current_scale_subscript ));

max_edge_length_in_voxels    = round( max( max_edge_length_in_microns ./ microns_per_voxel )) + 1 ;

%     numel_of_origin_vertex_strel = numel( structuring_element_linear_indexing_templates( round( current_scale_subscript )));

%     numel_of_origin_vertex_strel = round( 4 / 3 * pi * lumen_radius_in_microns_range( round( current_scale_subscript )) .^ 3 / prod( microns_per_voxel ));

max_number_of_indices        = max_edge_length_in_voxels * number_of_edges_per_vertex ;

previous_indices_visited     = zeros( max_number_of_indices, 1, 'int64' );

max_nnz_of_sparse_arrays     = double( max_number_of_indices * numel_of_strel );    

edge_indices_temp = zeros( number_of_edges_per_vertex, max_edge_length_in_voxels,  'int64' );
   edges2vertices = zeros( number_of_edges_per_vertex, 2,                         'uint32' ); % the indices of the vertices connected by each edge

   displacement_vectors               = zeros( 3, number_of_edges_per_vertex );
%        displacement_vector_sqr_magnitudes = zeros( 1, number_of_edges_per_vertex );

cum_prod_image_dims_double = double( cum_prod_image_dims );

   pointer_index_map         = spalloc( cum_prod_image_dims_double( 3 ), 1, max_nnz_of_sparse_arrays );
available_energy_map         = spalloc( cum_prod_image_dims_double( 3 ), 1, max_nnz_of_sparse_arrays );
%             distance_map         = spalloc( cum_prod_image_dims_double( 3 ), 1, max_nnz_of_sparse_arrays + 4 * numel_of_origin_vertex_strel );        
        distance_map         = spalloc( cum_prod_image_dims_double( 3 ), 1, max_nnz_of_sparse_arrays );
  pointer_energy_map         = spalloc( cum_prod_image_dims_double( 3 ), 1, max_nnz_of_sparse_arrays );
        
switch use_case
    
    case 'get_edges'
        
                  energy_map =          spalloc( cum_prod_image_dims_double( 3 ), 1, max_nnz_of_sparse_arrays ) ;
%                     size_map =          spalloc( cum_prod_image_dims_double( 3 ), 1, max_nnz_of_sparse_arrays ) ;                  
         previously_read_map = logical( spalloc( cum_prod_image_dims_double( 3 ), 1, max_nnz_of_sparse_arrays ));
         
    case { 'add_vertex_to_edge', 'extend_dead_end_edge', 'add_manual_edge' }
        
        energy_map =      energy_image ; % sparse masked image 
          size_map = round( size_image );
          
%     case 'extend_dead_end_edge'
         
end

number_of_indices = 0 ;

number_of_edges_found = 0 ;

cached_variable_1 = int64([ - 1,  1 ]);
cached_variable_2 = strel_apothem * int64([ 1; 1; 1 ]);
                
% current_index = double( vertex_indices( vertex_index ));

origin_position = double( index2position( current_index, cum_prod_image_dims )) .* microns_per_voxel' ;

distance_map( current_index ) = 1 ;

is_below_edge_max              = true ;
is_vertex_below_length_maximum = true ;
%     there_exists_possible_move     = true ;    

%     % don't enter the loop if the current vertex is in a index that can't be accessed by an edge
%     if ~ energy_map( current_index ), there_exists_possible_move = false ; end    
current_voxel_position = index2position( current_index, cum_prod_image_dims );

% restrict image size so that the strel never oversteps the image boundary
size_of_image = size_of_image' - strel_apothem ;

there_exists_possible_move = all( current_voxel_position > strel_apothem ) && all( current_voxel_position <= size_of_image );

% WHILE searching around in new indices
while is_vertex_below_length_maximum && there_exists_possible_move && is_below_edge_max

    switch use_case

        case 'get_edges'
                
            is_outside_box = previously_read_map( current_index ) == false ;

            if is_outside_box % then load a new chunk of the energy image

                current_voxel_position = index2position( current_index, cum_prod_image_dims );

                reading_box_limits = current_voxel_position + reading_box_apothem' .* cached_variable_1 ;
                
                % check where image boundary is to adjust box limits.            
                reading_box_limits_int64 = min( max( reading_box_limits, strel_apothem + 1 ), size_of_image );
                
                reading_box_limits = double( reading_box_limits_int64 );

                reading_box_starts  = reading_box_limits( :, 1 );

                reading_box_lengths = reading_box_limits( :, 2 ) - reading_box_limits( :, 1 ) + 1 ;

                reading_indices = cell( 1, 3 );
                
                reading_indices{ 1 } = reading_box_limits_int64( 1, 1 ) : reading_box_limits_int64( 1, 2 );
                reading_indices{ 2 } = reading_box_limits_int64( 2, 1 ) : reading_box_limits_int64( 2, 2 );
                reading_indices{ 3 } = reading_box_limits_int64( 3, 1 ) : reading_box_limits_int64( 3, 2 );
                
                % mark the box marked as read in on the previously_read_map slightly smaller than
                % the box of energy values read in (unless the reading box is at the image boundary,
                % then the entries will coincide) to avoid edge of box issues. Linear indexing the
                % box inside the image:
                reading_indices_linear_mesh                                                                                                                                                            ...
                    = reshape(  reading_indices{ 1 }                                 , length( reading_indices{ 1 }) * [ 1, 0, 0 ] + [ 0, 1, 1 ]) ...
                    + reshape(( reading_indices{ 2 } - 1 ) * cum_prod_image_dims( 1 ), length( reading_indices{ 2 }) * [ 0, 1, 0 ] + [ 1, 0, 1 ]) ...
                    + reshape(( reading_indices{ 3 } - 1 ) * cum_prod_image_dims( 2 ), length( reading_indices{ 3 }) * [ 0, 0, 1 ] + [ 1, 1, 0 ]) ;
                
                is_not_at_image_boundary = int64( reading_box_limits_int64 ~= [ cached_variable_2, size_of_image ]);

                previously_read_map( reading_indices_linear_mesh(     1 + strel_apothem * is_not_at_image_boundary( 1, 1 )  ...
                                                                  : end - strel_apothem * is_not_at_image_boundary( 1, 2 ), ...
                                                                      1 + strel_apothem * is_not_at_image_boundary( 2, 1 )  ...
                                                                  : end - strel_apothem * is_not_at_image_boundary( 2, 2 ), ...
                                                                      1 + strel_apothem * is_not_at_image_boundary( 3, 1 )  ...
                                                                  : end - strel_apothem * is_not_at_image_boundary( 3, 2 )  ))  ...
                                                                                                                                = true ;
                
                % load the energy data and associated scale index data (two 3D images concatenated into 4D array )
                energy_map( reading_indices_linear_mesh )                                                               ...
                                                            = h52mat( path_to_energy_data, [ reading_box_starts  ; 2 ], ...
                                                                                           [ reading_box_lengths ; 1 ]) ;
                                                                                       
%                 switch use_case
%                     
%                     case { 'add_vertex_to_edge', 'extend_dead_end_edge', 'add_manual_edge' }
%                         
%                         size_map( reading_indices_linear_mesh )                                     ...
%                                         = h52mat( path_to_energy_data, [ reading_box_starts  ; 1 ], ...
%                                                                        [ reading_box_lengths ; 1 ]) ;
% 
%                 end

%                 % zero out the borders of the energy image so that the strel will never reach the
%                 % end of the image. Linear indexing of the borders of the energy map image are set
%                 % to zero:
%                 energy_map([ 0 : strel_apothem                            - 1, cum_prod_image_dims( 1 ) - strel_apothem                            : cum_prod_image_dims( 1 ) - 1 ]' + ( 1 : cum_prod_image_dims( 1 ) : cum_prod_image_dims( 3 ))) = 0 ;
%                 energy_map([ 0 : strel_apothem * cum_prod_image_dims( 1 ) - 1, cum_prod_image_dims( 2 ) - strel_apothem * cum_prod_image_dims( 1 ) : cum_prod_image_dims( 2 ) - 1 ]' + ( 1 : cum_prod_image_dims( 2 ) : cum_prod_image_dims( 3 ))) = 0 ;
%                 energy_map([ 0 : strel_apothem * cum_prod_image_dims( 2 ) - 1, cum_prod_image_dims( 3 ) - strel_apothem * cum_prod_image_dims( 2 ) : cum_prod_image_dims( 3 ) - 1 ]' + ( 1 : cum_prod_image_dims( 3 ) : cum_prod_image_dims( 3 ))) = 0 ;
                    
            end % IF loading new chunk of energy image
            
        case { 'add_vertex_to_edge', 'extend_dead_end_edge', 'add_manual_edge' }
            
            % energy and size data are already in workspace as sparse arrays. They are nonzero only
            % at edge centerline voxels

    end % SIWTCH use_case
    
    % record the newly visited index and add a neighborhood surrounding (size/shape of the
    % strel) to the available energies to choose from.  Also keep track of the distances of the
    % newly revealed indices and mark them with a pointer that points those indices back to
    % this one.  Also need to check if the added index belongs to a vertex (distinct from the
    % origin vertex).
    number_of_indices               = number_of_indices + 1 ;

    is_vertex_below_length_maximum  = number_of_indices < max_number_of_indices ;

    previous_indices_visited( number_of_indices )   = current_index ;

    current_linear_strel                            = current_index + linear_strel ;

%     is_index_new   = ~ distance_map( current_linear_strel );

    % location will not be considered a new index again
    pointer_energy_map( current_index ) = - Inf ;

    is_index_new   = pointer_energy_map( current_linear_strel ) > energy_map( current_index );    
    
    new_indices_considered  = current_linear_strel( is_index_new );

    % distance and pointer_index map will have a few more entries than the available_energy map,
    % because distance and pointer are filled out here. available energy is populated after
    % restricting the list of new_indices_considered based on the distance limit.
    distance_map( new_indices_considered )                                                      ...
                    = distance_map( current_index ) + strel_distance_LUT( is_index_new );
                
    pointer_index_map( new_indices_considered ) = number_of_indices ;
    
    pointer_energy_map( new_indices_considered ) = energy_map( current_index  );

    new_indices_considered = new_indices_considered(   distance_map( new_indices_considered ) ...
                                                     < max_edge_length_in_microns             );

    % IF any edges have been found
    if any( edges2vertices( :, 1 ))

        length_of_new_indices_considered = length( new_indices_considered );

        if length_of_new_indices_considered

            indices_beyond_found_vertex = zeros( 1, length_of_new_indices_considered );

            % convert the linear indexing to spatial subscripts
%             pos_xy  = rem( double( new_indices_considered' ) - 1, cum_prod_image_dims( 2 )) + 1 ;
% 
%             pos_z   = ( double( new_indices_considered' ) - pos_xy ) / cum_prod_image_dims( 2 ) + 1 ;
            pos_xy  = rem( new_indices_considered' - 1, cum_prod_image_dims( 2 )) + 1 ;

            pos_z   = ( new_indices_considered' - pos_xy ) / cum_prod_image_dims( 2 ) + 1 ;

            pos_y   = rem( pos_xy - 1, cum_prod_image_dims( 1 )) + 1 ;

            pos_x   = ( pos_xy - pos_y ) / cum_prod_image_dims( 1 ) + 1 ;

            vectors_from_origin_to_indices = double([ pos_y; pos_x; pos_z ]) .* microns_per_voxel' - origin_position ;

            for terminal_index = 1 : number_of_edges_found

                terminal_vertex_displacements = displacement_vectors( :, terminal_index );

                indices_beyond_this_found_vertex = sum( terminal_vertex_displacements .* vectors_from_origin_to_indices, 1 ) > 1 ;

                indices_beyond_found_vertex = indices_beyond_found_vertex | indices_beyond_this_found_vertex ;

            end

            new_indices_considered( indices_beyond_found_vertex ) = [ ];

        end

    end

    % IF terminal vertex found
    if terminal_vertex_index

        number_of_edges_found = number_of_edges_found + 1 ;

%                 distance_map( new_indices_considered ) = Inf ;
%         available_energy_map( new_indices_considered ) =   0 ;
        
%                 distance_map( current_index ) = Inf ;
%         available_energy_map( current_index ) =   0 ;

        tracing_index = current_index ;

        tracing_ordinate = 0 ;

        % trace back until hitting a zero in the pointer index map, recording the linear indices
        % into the image and setting the pointer index map to zero along the way.  If this
        % process hits a zero index that doesn't belong to a vertex then record this for
        % later (we will need to add a new vertex and break up the existing edge into two).
        while pointer_index_map( tracing_index ) > 0

            tracing_ordinate = tracing_ordinate + 1 ;

            edge_indices_temp( number_of_edges_found, tracing_ordinate ) = tracing_index ;

            tracing_index = previous_indices_visited( pointer_index_map( tracing_index ));

        end % WHILE tracing

        % Erase the pointers for the found trajectory, so that this path will not contribute to
        % any other edge from this vertex. Replace the pointers with the negative of the parent
        % edge identity that owns this path.
        pointer_index_map( nonzeros( edge_indices_temp( number_of_edges_found, : ))) = - number_of_edges_found ;

        % include the terminal index as the last entry in the list
        edge_indices_temp( number_of_edges_found, tracing_ordinate + 1 ) = tracing_index ;

        % Checking if there exists a parent edge to this one with worse energy (flagging it for
        % removal, because we can't support children with better energy than their parents.  (A
        % child requires a parent to make a physical connection between two vertices.))
        parent_index = - pointer_index_map( tracing_index );

%         % !!!! this section seems superfluous because it should be impossible to get a parent with
%         worse energy thatn the current edge, because the current edge has the worst energy
%         available, no? I see for other use_case values like add_vertex_to_edge, it might be
%         possible to get this situation. There should be a swtich case here then it would seem, so
%         that it doesn't run during the get_edges_V204 function. SAM 8/14/20
        
        % IF edge has a parent
        if parent_index

            % check if this parent edge is itself a child
            parent_pointers = unique( - pointer_index_map( nonzeros( edge_indices_temp( parent_index, : ))));

            parent_pointers( parent_pointers == 0 | parent_pointers == parent_index ) = 0 ;

            % parent is also child
            if any( parent_pointers )

%                     is_edge_valid = false ;

                terminal_vertex_index = 0 ; origin_vertex_index = 0 ;                    

            else % ELSE edge is valid

                parent_indices = nonzeros( edge_indices_temp(          parent_index, : ));
                 child_indices = nonzeros( edge_indices_temp( number_of_edges_found, : ));

                parent_energies = energy_map( parent_indices );
                 child_energies = energy_map(  child_indices );

%                     parent_energy = get_edge_metric({ parent_energies });
%                      child_energy = get_edge_metric({  child_energies });

                parent_energy = max( parent_energies );
                 child_energy = max(  child_energies );

%                     % IF child has higher (worse) energy than parent, AND parent was valid.
%                     is_edge_valid = child_energy > parent_energy                    ...
%                                   & edges2vertices( parent_index, 1 );

                % IF child has higher (worse) energy than parent, AND parent was valid.
                if child_energy > parent_energy && edges2vertices( parent_index, 1 )

%                         % set the origin vertex to the terminal vertex of the parent
%                         origin_vertex_index = edges2vertices( parent_index, 1 )

%                         origin_vertex_index = vertex_index ;

                    % set the origin vertex to the vertex that has the lowest energy connection
                    % to the bifurcation point where the child split off
                    bifurcation_index = find( parent_indices == child_indices( end ));

                    parent_1_energies = energy_map( parent_indices( 1 : bifurcation_index - 1       ));
                    parent_2_energies = energy_map( parent_indices(     bifurcation_index + 1 : end ));

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

                    origin_vertex_index = edges2vertices( parent_index, better_half_of_parent );

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
        
        if terminal_vertex_index < Inf % Inf index is code for new vertex needs to be created
            
            edges2vertices( number_of_edges_found, : ) = [ terminal_vertex_index, origin_vertex_index ];          
        
        else % Inf index is code for new vertex
        
            edges2vertices( number_of_edges_found, : ) = [ number_of_vertices + 1, origin_vertex_index ]; 
            
        end
        
        current_position = double( index2position( current_index, cum_prod_image_dims )) .* microns_per_voxel' ;

        displacement_vectors( :, number_of_edges_found ) = current_position - origin_position ;

        % divide by its magnitude squared so that later projections can be compared to unity
        displacement_vectors( :, number_of_edges_found ) = displacement_vectors( :, number_of_edges_found ) / sum( displacement_vectors( :, number_of_edges_found ) .^ 2 );

        % If we find this edge not valid by the above rule, we will place a zero in the
        % edges2vertices variable for this edge's entry.  Thus, a zero in the
        % edges2vertices variable is code for "no edge was found"            

%             end % IF valid edge

        available_energy_map( new_indices_considered ) = 0 ;

    else % ELSE no vertex exists at the current index

        available_energy_map( new_indices_considered ) = energy_map( new_indices_considered );
        
%         available_energy_map( current_index ) = 0 ;

    end % IF terminal vertex found
    
    % location will not be explored again
    available_energy_map( current_index ) = 0 ;  

    is_below_edge_max = number_of_edges_found < number_of_edges_per_vertex ;

    % choose next index as the next available with lowest energy
    [ min_available_energy, current_index ] = min( available_energy_map );
    
    current_index = int64( current_index );

    terminal_vertex_index = full( vertex_image( current_index ));
    
    switch use_case

        case { 'add_vertex_to_edge', 'extend_dead_end_edge', 'add_manual_edge' }
                        
            if ~ terminal_vertex_index % Then the current index is not the center of a neighboring vertex
                
                % So try to paint a vertex here and check volume conflicts with existing vertices
                if all( ~ vertex_volume_image( current_index + strel_linear_indexing_templates{ size_map( current_index )}))
%                      min( max( current_index + strel_linear_indexing_templates{ size_map( current_index )}, 1 ), cum_prod_image_dims( 3 )) % image boundaries
                    
                    terminal_vertex_index = Inf ; % Inf is code for needing to make a new vertex
                    
                end
            end
    end

    % requiring negative edge energy
    there_exists_possible_move = min_available_energy < 0 ;

end % WHILE searching for edges

end % FUNCTION get_edges_for_vertex