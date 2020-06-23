function [ edge_space_subscripts, edge_scale_subscripts, edge_energies, edges2vertices, vertex_space_subscripts, vertex_scale_subscripts, vertex_energies ] ...
             = add_vertices_to_edges( edge_space_subscripts, edge_scale_subscripts, edge_energies, edges2vertices, vertex_space_subscripts, vertex_scale_subscripts, vertex_energies, size_of_image )
%% Add vertices where children edges meet their parents.  
% all of these set calculations must be done stably so that the new_vertex_location FOR loop is done
% from best energy to worst energy

% number_of_original_edges = length( edges2vertices );
% 
% original_edge_indices = 1 : number_of_original_edges ;

edge_space_subscripts_matrix = double( cell2mat( edge_space_subscripts ));

vertex_space_subscripts_double = double( vertex_space_subscripts );

edge_locations   = sub2ind( size_of_image,   edge_space_subscripts_matrix( :, 1 ),   edge_space_subscripts_matrix( :, 2 ),   edge_space_subscripts_matrix( :, 3 ));
vertex_locations = sub2ind( size_of_image, vertex_space_subscripts_double( :, 1 ), vertex_space_subscripts_double( :, 2 ), vertex_space_subscripts_double( :, 3 ));

edge_lengths = cellfun( @length, edge_scale_subscripts );

edge_locations_cell = mat2cell( edge_locations, edge_lengths, 1 );

edge_locations = cell2mat( edge_locations_cell );

% edge_lengths   = cellfun( @length, edge_locations_cell );

% this unique need not be stable because all of these indices will be erased, so their order doesn't
% matter.  Throw out the last occurence so that the parent edge remains the representative of its
% parent/child intersection location.
[ ~, unique_edge_location_indices ] = unique( edge_locations, 'last' );

% throw out all of the unique ones and save the repeats
edge_location_repeats = edge_locations ; edge_location_repeats( unique_edge_location_indices ) = [ ];

% the vertex locations to be added are the edge locations that are repeated at the end point of at
% least one edge, and at the midpoint of at least one edge, and where there is not already a vertex.
exterior_edge_locations = cell2mat( cellfun( @( x ) x([ 1,    end   ]), edge_locations_cell, 'UniformOutput', false ));
interior_edge_locations = cell2mat( cellfun( @( x ) x(  2 : end - 1  ), edge_locations_cell, 'UniformOutput', false ));

non_vertex_edge_location_repeats          =   setdiff(            edge_location_repeats,                 vertex_locations, 'stable' ) ;

non_vertex_edge_exterior_location_repeats = intersect( non_vertex_edge_location_repeats,          exterior_edge_locations, 'stable' ) ;

new_vertex_locations                      = intersect( non_vertex_edge_exterior_location_repeats, interior_edge_locations, 'stable' )';

number_of_new_vertices = length( new_vertex_locations );

if number_of_new_vertices
                      
    new_vertex_space_subscripts = zeros( number_of_new_vertices, 3, 'uint16' );
    new_vertex_scale_subscripts = zeros( number_of_new_vertices, 1           );
    new_vertex_energies         = zeros( number_of_new_vertices, 1           );

    number_of_vertices       = length( vertex_locations );

    number_of_vertices_added = 0 ;

    % remove edges that will not be part of this child/parent reconfiguration FOR loop

    is_edge_parent_or_child = cellfun( @( x ) logical( length( find( x == new_vertex_locations, 1 ))), edge_locations_cell );

    % temporarilly remove the unnecessary edges for the following FOR loop
    % old_edge_locations_cell   =       edge_locations_cell( ~ is_edge_parent_or_child    );
    old_edges2vertices        =            edges2vertices( ~ is_edge_parent_or_child, : );
    old_edge_space_subscripts =     edge_space_subscripts( ~ is_edge_parent_or_child    );
    old_edge_scale_subscripts =     edge_scale_subscripts( ~ is_edge_parent_or_child    );
    old_edge_energies         =             edge_energies( ~ is_edge_parent_or_child    );

        edge_locations_cell   =       edge_locations_cell(   is_edge_parent_or_child    );
        edges2vertices        =            edges2vertices(   is_edge_parent_or_child, : );
        edge_space_subscripts =     edge_space_subscripts(   is_edge_parent_or_child    );
        edge_scale_subscripts =     edge_scale_subscripts(   is_edge_parent_or_child    );
        edge_energies         =             edge_energies(   is_edge_parent_or_child    );

    edge_locations = cell2mat( edge_locations_cell );

    edge_lengths   = cellfun( @length, edge_locations_cell );

    % look for the new edges to add to each added vertex

    for new_vertex_location = new_vertex_locations

%         is_parent_established = false ;

        cumulative_edge_lengths = cumsum( edge_lengths );

        edge_indices_to_delete  = [ ];        

%         % We will again force children to be worse energy than their parents.  In this context, this
%         % means that we can't start placing children (edges where the vertex is on an exterior location)
%         % until we get a parent at this location.  If the parent is the last edge in the list, then we
%         % also should not break it up.  Remember that the location list is sorted by the edge energy
%         % that contains each location.
%         is_child_established = false ;        

        % Locate all the locations on the list of edge locations that coincide with this vertex location
        % in order to remove vertices from the list that won't have a valid child on the basis of energy
        edge_location_indices_at_new_vertex = find( edge_locations == new_vertex_location )';

%         for edge_location_index = edge_location_indices_at_new_vertex
% 
%             edge_index = find( edge_location_index <= cumulative_edge_lengths, 1 );
% 
%             relative_edge_index_at_new_vertex = edge_location_index - cumulative_edge_lengths( edge_index ) + edge_lengths( edge_index );
% 
%             % IF vertex is on the first or last location of the edge
%             if relative_edge_index_at_new_vertex == edge_lengths( edge_index )
% 
%                 if is_parent_established 
%                                 
%                     is_child_established = true ;
% 
%                 end
% 
%             else % ELSE vertex is on an interior location of the edge
% 
%                 is_parent_established = true ;
% 
%             end % ELSE vertex is on an interior location of the edge
%         end % FOR edge_location_indices_at_new_vertex
% 
%         is_parent_established = false ;
% 
%         % IF there will be a parent-child relationship and thus a reason to place a new vertex

        for edge_location_index = edge_location_indices_at_new_vertex
            
            edge_index = find( edge_location_index <= cumulative_edge_lengths, 1 );

            relative_edge_index_at_new_vertex = edge_location_index - cumulative_edge_lengths( edge_index ) + edge_lengths( edge_index );

            % if first edge having a vertex added to it at this location
            if edge_location_index == edge_location_indices_at_new_vertex( 1 )
                
                number_of_vertices       = number_of_vertices       + 1 ;
                number_of_vertices_added = number_of_vertices_added + 1 ;

                % record the new vertex properties from this location in the edge
                new_vertex_space_subscripts( number_of_vertices_added, : ) = edge_space_subscripts{ edge_index }( relative_edge_index_at_new_vertex, : );
                new_vertex_scale_subscripts( number_of_vertices_added    ) = edge_scale_subscripts{ edge_index }( relative_edge_index_at_new_vertex    );        
                        new_vertex_energies( number_of_vertices_added    ) =         edge_energies{ edge_index }( relative_edge_index_at_new_vertex    );                    

            end            
            
            % determine whether this location is an interior or exterior edge location, based on edge lengths
            
            % IF vertex is on the last location of the edge            
            if relative_edge_index_at_new_vertex == edge_lengths( edge_index )

%                 if is_parent_established 

                    edges2vertices( edge_index, end ) = number_of_vertices ;
                    
%                 else % this edge must be deleted
% 
%                     % !!!! this line was not working as intended.  Not including it.  SAM 7/23/19
%                     edge_indices_to_delete  = [ edge_indices_to_delete,  edge_index ]
% 
%                 end

            else % ELSE vertex is on an interior location of the edge

%                 if is_child_established
                    
% % % % % commented out and partially moved to just before this FOR loop 1/13/20
% % % % 
% % % %                     % IF first time reaching a parent edge coinciding with this vertex                
% % % %                     if ~ is_parent_established
% % % % 
% % % %                         % add the vertex to the list of vertices
% % % %                         is_parent_established = true ;
% % % % 
% % % %                         number_of_vertices       = number_of_vertices       + 1 ;
% % % %                         number_of_vertices_added = number_of_vertices_added + 1 ;
% % % % 
% % % %                         % record the new vertex properties from this location in the edge
% % % %                         new_vertex_space_subscripts( number_of_vertices_added, : ) = edge_space_subscripts{ edge_index }( relative_edge_index_at_new_vertex, : );
% % % %                         new_vertex_scale_subscripts( number_of_vertices_added    ) = edge_scale_subscripts{ edge_index }( relative_edge_index_at_new_vertex    );        
% % % %                                 new_vertex_energies( number_of_vertices_added    ) =         edge_energies{ edge_index }( relative_edge_index_at_new_vertex    );                    
% % % % 
% % % %                     end % IF first time reaching a parent edge coinciding with this vertex

                    % break this parent edge into two

                    % interior edge locations need to have their edges deleted and split into two new ones
                    edge_indices_to_delete  = [ edge_indices_to_delete,  edge_index ];

                    % record the new edge-vertex connectivity for the two added edges
                    new_edges2vertices{ 1, 1 } = [ edges2vertices( edge_index, 1 ); number_of_vertices ];
                    new_edges2vertices{ 2, 1 } = [ number_of_vertices;  edges2vertices( edge_index, 2 )];
                    
                    % break the edge into two pieces (the part before the new vertex and the part after it)
                    new_edge_location_indices_into_old = { 1 : relative_edge_index_at_new_vertex                             ; ...
                                                               relative_edge_index_at_new_vertex : edge_lengths( edge_index )  };

                    for pair_index = 1 : 2

                        % extract the information about the added edge from half of the old one
                               new_edge_locations{ pair_index, 1 } =   edge_locations_cell{ edge_index }( new_edge_location_indices_into_old{ pair_index }    );
                        new_edge_space_subscripts{ pair_index, 1 } = edge_space_subscripts{ edge_index }( new_edge_location_indices_into_old{ pair_index }, : );
                        new_edge_scale_subscripts{ pair_index, 1 } = edge_scale_subscripts{ edge_index }( new_edge_location_indices_into_old{ pair_index }    );
                                new_edge_energies{ pair_index, 1 } =         edge_energies{ edge_index }( new_edge_location_indices_into_old{ pair_index }    );

                    end % FOR added edge_pair to be

                    % add the pair of edges to the edge lists
                    edges2vertices          = [ edges2vertices        ;      [ new_edges2vertices{ : }]'];    
                    edge_locations_cell     = [ edge_locations_cell   ;        new_edge_locations( : ) ];
                    edge_space_subscripts   = [ edge_space_subscripts ; new_edge_space_subscripts( : )  ];
                    edge_scale_subscripts   = [ edge_scale_subscripts ; new_edge_scale_subscripts( : )  ];
                    edge_energies           = [ edge_energies         ;         new_edge_energies( : )  ];

                    edge_lengths            = [ edge_lengths          ; cellfun( @length, new_edge_locations( : ))];

                    edge_locations          = [ edge_locations        ; cat( 1, new_edge_locations{ : })];                 

%                 end % IF is_child_established
                
            end % ELSE vertex is interior to edge

        end % FOR edge_location_indices_at_new_vertex

        edge_location_indices_to_delete = [ ];

        % delete the edge(s) that was replaced by the pair of edges
        for edge_index = edge_indices_to_delete

            edge_location_indices_to_delete = [ edge_location_indices_to_delete, cumulative_edge_lengths( edge_index ) - edge_lengths( edge_index ) + 1 : cumulative_edge_lengths( edge_index )];


        end % FOR edge_indices_to_delete

        edge_locations( edge_location_indices_to_delete ) = [ ]; 

               edges2vertices( edge_indices_to_delete, : ) = [ ];    
          edge_locations_cell( edge_indices_to_delete    ) = [ ];
        edge_space_subscripts( edge_indices_to_delete    ) = [ ];      
        edge_scale_subscripts( edge_indices_to_delete    ) = [ ];
                edge_energies( edge_indices_to_delete    ) = [ ];
                 edge_lengths( edge_indices_to_delete    ) = [ ];

    end % FOR new_vertex_locations

    % add the new vertex(-ices)
    % vertex_locations           = [ vertex_locations        ; new_vertex_locations'       ]; % !!! this list is misalligned from the others
    vertex_space_subscripts    = [ vertex_space_subscripts ; new_vertex_space_subscripts ];
    vertex_scale_subscripts    = [ vertex_scale_subscripts ; new_vertex_scale_subscripts ];
    vertex_energies            = [ vertex_energies         ; new_vertex_energies         ];
    
    unplaced_vertices = vertex_scale_subscripts == 0 ;
    
    vertex_space_subscripts( unplaced_vertices, : ) = [ ];
    vertex_scale_subscripts( unplaced_vertices    ) = [ ];
    vertex_energies(         unplaced_vertices    ) = [ ];

    % put the temporarilly removed edges back in
    % edge_locations_cell     = [ old_edge_locations_cell   ;       edge_locations_cell ];
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

end % function
