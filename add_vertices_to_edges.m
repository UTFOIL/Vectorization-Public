function [ edge_space_subscripts, edge_scale_subscripts, edge_energies, edges2vertices, vertex_space_subscripts, vertex_scale_subscripts, vertex_energies, bridge_edges ] ...
             = add_vertices_to_edges( edge_space_subscripts, edge_scale_subscripts, edge_energies, edges2vertices, vertex_space_subscripts, vertex_scale_subscripts, vertex_energies, size_of_image, microns_per_voxel, lumen_radius_in_microns_range, data_directory, energy_handle, strel_apothem, max_edge_length_per_origin_radius )
%% Add vertices where children edges meet their parents.  
% all of these set calculations must be done stably so that the new_vertex_index FOR loop is done
% from best energy to worst energy

% number_of_original_edges = length( edges2vertices );
% 
% original_edge_indices = 1 : number_of_original_edges ;

% edge_space_subscripts_matrix = double( cell2mat( edge_space_subscripts ));
% 
% vertex_space_subscripts_double = double( vertex_space_subscripts );

edge_space_subscripts_matrix = cell2mat( edge_space_subscripts );

% vertex_space_subscripts = vertex_space_subscripts ;

edge_indices   = int64( sub2ind( size_of_image, edge_space_subscripts_matrix( :, 1 ), edge_space_subscripts_matrix( :, 2 ), edge_space_subscripts_matrix( :, 3 )));
vertex_indices = int64( sub2ind( size_of_image,      vertex_space_subscripts( :, 1 ),      vertex_space_subscripts( :, 2 ),      vertex_space_subscripts( :, 3 )));

edge_lengths = cellfun( @length, edge_scale_subscripts );

edge_indices_cell = mat2cell( edge_indices, edge_lengths, 1 );

edge_indices = cell2mat( edge_indices_cell );

% edge_lengths   = cellfun( @length, edge_indices_cell );

% this unique need not be stable because all of these indices will be erased, so their order doesn't
% matter.  Throw out the last occurence so that the parent edge remains the representative of its
% parent/child intersection index.
[ ~, unique_edge_index_indices ] = unique( edge_indices, 'last' );

% throw out all of the unique ones and save the repeats
edge_index_repeats = edge_indices ; edge_index_repeats( unique_edge_index_indices ) = [ ];

% the vertex indices to be added are the edge indices that are repeated at the end point of at
% least one edge, and at the midpoint of at least one edge, and where there is not already a vertex.
exterior_edge_indices = cell2mat( cellfun( @( x ) x([ 1,    end   ]), edge_indices_cell, 'UniformOutput', false ));
interior_edge_indices = cell2mat( cellfun( @( x ) x(  2 : end - 1  ), edge_indices_cell, 'UniformOutput', false ));

% ?? computationally efficiency could be increased by subtracting out the vertex indices earlier??
non_vertex_edge_index_repeats          =   setdiff(            edge_index_repeats,                 vertex_indices, 'stable' ) ;
non_vertex_edge_exterior_index_repeats = intersect( non_vertex_edge_index_repeats,          exterior_edge_indices, 'stable' ) ;
new_vertex_indices                     = intersect( non_vertex_edge_exterior_index_repeats, interior_edge_indices, 'stable' )';

max_number_of_new_vertices = length( new_vertex_indices );

strel_linear_indexing_templates = construct_structuring_elements( lumen_radius_in_microns_range, microns_per_voxel, size_of_image );

% paint a (binary) reference image with the vertex volumes filled in
[ vertex_volume_image ] = logical( paint_vertex_image( vertex_space_subscripts, vertex_scale_subscripts, vertex_energies, ...
                                                       lumen_radius_in_microns_range ./ microns_per_voxel, size_of_image                             ));

% generate a sparsely indexed image of vertex indices at their centerpoints
% create a reading box element for reading in local cubes of energy image during exploration
% create conversion keys between subscripts and linear indices

% number_of_edges_per_vertex = 1 ;

[ vertex_center_image, vertex_indices, reading_box_apothem, linear_strel,     ...
                          max_edge_length_in_microns_range, numel_of_strel,   ...
                   cum_prod_image_dims, size_of_image, path_to_energy_data, strel_distance_LUT ] ...
  = generate_reference_image( lumen_radius_in_microns_range, microns_per_voxel, ...
                              vertex_scale_subscripts, vertex_space_subscripts, ...
                              strel_apothem, max_edge_length_per_origin_radius, ...
                                                 data_directory, energy_handle  );
                                      
                                      
                                           
% [ vertex_center_image, vertex_locations, reading_box_apothem, linear_strel,    ...
%                    max_edge_length_in_microns_range, numel_of_strel,             ...
%             cum_prod_image_dims, size_of_image, path_to_energy_data, strel_distance_LUT ]            ...
%   = generate_reference_image(  lumen_radius_in_microns_range, microns_per_voxel, ...
%                                vertex_scale_subscripts, vertex_space_subscripts, ...
%                                strel_apothem, max_edge_length_per_origin_radius, ...
%                                                   data_directory, energy_handle  )
                                                             
if max_number_of_new_vertices
                      
    new_vertex_space_subscripts = zeros( max_number_of_new_vertices, 3, 'uint16' );
    new_vertex_scale_subscripts = zeros( max_number_of_new_vertices, 1           );
    new_vertex_energies         = zeros( max_number_of_new_vertices, 1           );

    number_of_vertices       = length( vertex_indices );

    number_of_new_vertices = 0 ;

    % remove edges that will not be part of this child/parent reconfiguration FOR loop

    is_edge_parent_or_child = cellfun( @( x ) logical( length( find( x == new_vertex_indices, 1 ))), edge_indices_cell );

    % temporarilly remove the unnecessary edges for the following FOR loop
    % old_edge_indices_cell   =       edge_indices_cell( ~ is_edge_parent_or_child    );
    old_edges2vertices        =            edges2vertices( ~ is_edge_parent_or_child, : );
    old_edge_space_subscripts =     edge_space_subscripts( ~ is_edge_parent_or_child    );
    old_edge_scale_subscripts =     edge_scale_subscripts( ~ is_edge_parent_or_child    );
    old_edge_energies         =             edge_energies( ~ is_edge_parent_or_child    );

        edge_indices_cell     =         edge_indices_cell(   is_edge_parent_or_child    );
        edges2vertices        =            edges2vertices(   is_edge_parent_or_child, : );
        edge_space_subscripts =     edge_space_subscripts(   is_edge_parent_or_child    );
        edge_scale_subscripts =     edge_scale_subscripts(   is_edge_parent_or_child    );
        edge_energies         =             edge_energies(   is_edge_parent_or_child    );

    edge_indices = cell2mat(         edge_indices_cell )  ;
    edge_lengths = cellfun( @length, edge_indices_cell )' ;
    
    use_case = 'add_vertex_to_edge' ;
    
    % load the energy data and associated scale index data (two 3D images concatenated into 4D array )
    energy_and_index_data = h52mat( path_to_energy_data );

    energy_image = energy_and_index_data( :, :, :, 2 );
      size_image = energy_and_index_data( :, :, :, 1 );

    clear( 'energy_and_index_data' )

    edge_indices_unique = unique( edge_indices );
    
    % mask in the energy image where the edges are    
    energy_image = sparse( double( edge_indices_unique ), 1, energy_image( edge_indices_unique ), double( cum_prod_image_dims( 3 )), 1 );
      size_image = sparse( double( edge_indices_unique ), 1,   size_image( edge_indices_unique ), double( cum_prod_image_dims( 3 )), 1 );
      
        edge_indices_temp = cell( max_number_of_new_vertices, 1 );
      edges2vertices_temp = cell( max_number_of_new_vertices, 1 );
      
    number_of_bridge_edges = 0 ;

    % look for the new edges to add to each added vertex
    for new_vertex_index = new_vertex_indices

        number_of_bridge_edges = number_of_bridge_edges + 1 ;
        
        cumulative_edge_lengths = cumsum( edge_lengths );
        
        % Locate all the jndices into the list of edge indices (into the image) that coincide
        % with this vertex index (into the image)
        jndices = find( edge_indices == new_vertex_index )';

        number_of_jndices = length( jndices );
                
        edges_at_vertex = zeros( 1, number_of_jndices );

        for kndex = 1 : number_of_jndices

            edges_at_vertex( kndex ) = find( jndices( kndex ) <= cumulative_edge_lengths, 1 );
            
        end

        inter_edges_at_vertex = jndices - cumulative_edge_lengths( edges_at_vertex ) + edge_lengths( edges_at_vertex );

        % --- Determine whether this index is an interior or exterior edge index, based on
        % edge lengths:

        % IF new vertex is on the last index of the edge, then it is a child edge          
        is_child = inter_edges_at_vertex == edge_lengths( edges_at_vertex );

%         % IF the a child edge is the most energetically favorable (first on the list), then that
%         % means that it's parent was erased.  So erase children edges until a parent edge reveals itself
%         while is_child( 1 )
%             
%             % !!!! move this outside and into the clean orphans functinos for clear semantics
%             error
%             
%         end
        
        energy_image_temp = energy_image ;
        
        child_indices_to_be_deleted = cell2mat( cellfun( @( v ) v( 1 : end - 1 ), edge_indices_cell( edges_at_vertex( is_child )), 'UniformOutput', false ));
        
        energy_image_temp( child_indices_to_be_deleted ) = 0 ;
        
%         energy_image_temp( edge_indices( cumulative_edge_lengths( edges_at_vertex( kndex )))) = 0 ;

        % get the edge segment from the new_vertex_index to the new_vertex_index_updated which
        % is guaranteed to be a nonoverlapping vertex.
        [ edge_indices_temp{ number_of_bridge_edges }, edges2vertices_temp{ number_of_bridge_edges }]                                                ...
                           = get_edges_for_vertex( strel_linear_indexing_templates, vertex_volume_image, energy_image_temp, size_image, microns_per_voxel, strel_apothem, ...
                                      1, use_case, new_vertex_index, ...
                       vertex_center_image, vertex_indices, reading_box_apothem, linear_strel, ...
                                        max_edge_length_in_microns_range, numel_of_strel, strel_distance_LUT, ...
                                 cum_prod_image_dims, size_of_image, path_to_energy_data  );
                                     

                                     
        edge_indices_temp2 = nonzeros( edge_indices_temp{ number_of_bridge_edges } );
        
        new_vertex_index_updated = edge_indices_temp2( 1 );
                
%         is_extending_edge = new_vertex_index_updated ~= new_vertex_index ;
        
        is_creating_new_vertex = edges2vertices_temp{ number_of_bridge_edges }( 1 ) == number_of_vertices + 1 ;
        
%         % add the new edge segment to any children edges for this new_vertex_index
%         if is_extending_edge 
%             
%             
%             
%         end % IF `is_extending_edge

        % add the new vertex and break up edges accordingly in the vertex creation step
        if is_creating_new_vertex

            vertex_indices( end + 1 ) = new_vertex_index_updated ;
                        
            number_of_vertices     = number_of_vertices     + 1 ;
            number_of_new_vertices = number_of_new_vertices + 1 ;
            
            % record the new vertex properties from this index in the edge
            new_vertex_space_subscripts( number_of_new_vertices, : ) = index2position( new_vertex_index_updated, cum_prod_image_dims )';
            new_vertex_scale_subscripts( number_of_new_vertices    ) =     size_image( new_vertex_index_updated ) ;
                    new_vertex_energies( number_of_new_vertices    ) =   energy_image( new_vertex_index_updated ) ;
                    
            % Locate all the jndices into the list of edge indices (into the image) that coincide
            % with this vertex index (into the image)
            jndices = find( edge_indices == new_vertex_index )';
            
            number_of_jndices = length( jndices );
                       
            % rename the ending vertex for all the children edges
            edges2vertices( edges_at_vertex( is_child ), 2 ) = number_of_vertices ;

            % split parents at the new vertex into two edges whose concatenation equals the original parent
            parents_at_vertex = find( ~ is_child );
            
            for kndex = parents_at_vertex

                % record the new edge-vertex connectivity for the two added edges
                new_edges2vertices{ 1, 1 } = [ edges2vertices( edges_at_vertex( kndex ), 1 ); number_of_vertices ];
                new_edges2vertices{ 2, 1 } = [ number_of_vertices;  edges2vertices( edges_at_vertex( kndex ), 2 )];

                % break the edge into two pieces (the part before the new vertex and the part after it)
                new_edge_index_indices_into_old = { 1 : inter_edges_at_vertex( kndex )                                           ; ...
                                                        inter_edges_at_vertex( kndex ) : edge_lengths( edges_at_vertex( kndex ))};

                for pair_index = 1 : 2

                    % extract the information about the added edge from half of the old one
                             new_edge_indices{ pair_index, 1 } =     edge_indices_cell{ edges_at_vertex( kndex )}( new_edge_index_indices_into_old{ pair_index }    );
                    new_edge_space_subscripts{ pair_index, 1 } = edge_space_subscripts{ edges_at_vertex( kndex )}( new_edge_index_indices_into_old{ pair_index }, : );
                    new_edge_scale_subscripts{ pair_index, 1 } = edge_scale_subscripts{ edges_at_vertex( kndex )}( new_edge_index_indices_into_old{ pair_index }    );
                            new_edge_energies{ pair_index, 1 } =         edge_energies{ edges_at_vertex( kndex )}( new_edge_index_indices_into_old{ pair_index }    );

                            if length( new_edge_energies{ pair_index, 1 })==1
                                
                                'here'
                                
                            end
                            
                end % FOR added edge_pair to be

                % add the pair of edges to the edge lists
                edges2vertices          = [ edges2vertices        ;      [ new_edges2vertices{ : }]'];    
                edge_indices_cell       = [ edge_indices_cell     ;          new_edge_indices( : )  ];
                edge_space_subscripts   = [ edge_space_subscripts ; new_edge_space_subscripts( : )  ];
                edge_scale_subscripts   = [ edge_scale_subscripts ; new_edge_scale_subscripts( : )  ];
                edge_energies           = [ edge_energies         ;         new_edge_energies( : )  ];

                edge_lengths            = [ edge_lengths          , cellfun( @length, new_edge_indices( : )' )];
                edge_indices            = [ edge_indices          ;     cat(       1, new_edge_indices{ : }  )];                 

                %                 end % IF is_child_established

                %                 end % ELSE vertex is interior to edge

            end % FOR kindices into jindices

            edges_to_delete = parents_at_vertex ;

            jndices_to_delete = [ ];

            % delete the edge(s) that was replaced by the pair of edges
            for edge_at_vertex = edges_to_delete

                jndices_to_delete = [ jndices_to_delete, cumulative_edge_lengths( edge_at_vertex ) - edge_lengths( edge_at_vertex ) + 1 : cumulative_edge_lengths( edge_at_vertex )];

            end % FOR edge_indices_to_delete

            edge_indices( jndices_to_delete ) = [ ]; 

                   edges2vertices( edges_to_delete, : ) = [ ];    
                edge_indices_cell( edges_to_delete    ) = [ ];
            edge_space_subscripts( edges_to_delete    ) = [ ];      
            edge_scale_subscripts( edges_to_delete    ) = [ ];
                    edge_energies( edges_to_delete    ) = [ ];
                     edge_lengths( edges_to_delete    ) = [ ];

        else % not creating new vertex (links to existing vertex)

            edges2vertices( edges_at_vertex( is_child ), 2 ) = edges2vertices_temp{ number_of_bridge_edges }( 1 );
                     
        end % new vertex created

    end % FOR new_vertex_indices

    % add the new vertex(-ices)
%     vertex_indices           = [ vertex_indices        ; new_vertex_indices'       ]; % !!! this list is misalligned from the others
    vertex_space_subscripts    = [ vertex_space_subscripts ; new_vertex_space_subscripts ];
    vertex_scale_subscripts    = [ vertex_scale_subscripts ; new_vertex_scale_subscripts ];
    vertex_energies            = [ vertex_energies         ; new_vertex_energies         ];
    
    % remove vertex zero-entries (that were not placed)
    unplaced_vertices = vertex_scale_subscripts == 0 ;
    
    vertex_space_subscripts( unplaced_vertices, : ) = [ ];
    vertex_scale_subscripts( unplaced_vertices    ) = [ ];
    vertex_energies(         unplaced_vertices    ) = [ ];

    % put the temporarilly removed edges back in
    % edge_indices_cell     = [ old_edge_indices_cell   ;       edge_indices_cell ];
    edges2vertices          = [ old_edges2vertices        ;            edges2vertices ];    
    edge_space_subscripts   = [ old_edge_space_subscripts ;     edge_space_subscripts ];
    edge_scale_subscripts   = [ old_edge_scale_subscripts ;     edge_scale_subscripts ];
    edge_energies           = [ old_edge_energies         ;             edge_energies ];

    %% clean the edge pairs again

    [ edges2vertices, edge_indices_to_keep ] = clean_edge_pairs( edges2vertices, edge_energies, false );

    edge_space_subscripts  =  edge_space_subscripts( edge_indices_to_keep );
    edge_scale_subscripts  =  edge_scale_subscripts( edge_indices_to_keep );
    edge_energies          =          edge_energies( edge_indices_to_keep );

end % IF any vertices to add

[ edges2vertices_bridges, edge_lengths_bridges, edge_space_subscripts_bridges, edge_scale_subscripts_bridges, edge_energies_bridges ] ...
       = get_edge_vectors( 1, edge_indices_temp, edges2vertices_temp, path_to_energy_data, cum_prod_image_dims );
   
   edge_subscripts_bridges = cellfun( @( x, y ) [ x, y ], edge_space_subscripts_bridges, edge_scale_subscripts_bridges, 'UniformOutput', false );
   
   mean_edge_energies = cellfun( @mean, edge_energies_bridges );
   
bridge_edges = struct( 'edges2vertices',        edges2vertices_bridges,        ...
                       'edge_space_subscripts', edge_space_subscripts_bridges, ...
                       'edge_scale_subscripts', edge_scale_subscripts_bridges, ...
                       'edge_energies',         edge_energies_bridges,         ...
                       'mean_edge_energies',    mean_edge_energies             );
                       
% %   visualize the bridge edges that were created in this function
%                    
%                    
%         visualize_edges_V180( edge_subscripts_bridges, mean_edge_energies,    ...
%                               lumen_radius_in_microns_range ./ microns_per_voxel,               ...
%                               double( size_of_image ), 'bridge_edges_spheres.tif', ...
%                               'bridge_edges_centerlines.tif'              )

end % function

% function conflicting_indices = test_paint_new_vertex( to_be_painted_index )
% 
%    space_subscript_int64 = int64( vertex_space_subscripts( to_be_painted_index, : ));
% 
%    vertex_position_linear_index =  space_subscript_int64( :, 1 )                                               ...
%                                + ( space_subscript_int64( :, 2 ) - 1 ) * size_of_image( 1 )                    ...
%                                + ( space_subscript_int64( :, 3 ) - 1 ) * size_of_image( 1 ) * size_of_image( 2 );
% 
%     vertex_structure_positions_linear_indexing{ to_be_painted_index }                                    ...
%         = structuring_element_linear_indexing_templates{ round( vertex_scale_subscripts( to_be_painted_index ))} ...
%                                                                           + vertex_position_linear_index ;
% 
%     % look at the volume to be painted and see if there are vertices in the way
%     conflicting_indices = unique( index_image( vertex_structure_positions_linear_indexing{ to_be_painted_index }))';
% 
%     % don't try to delete backgrouund
%     conflicting_indices( conflicting_indices == 0 ) = [ ];
% 
% end % add_vertex