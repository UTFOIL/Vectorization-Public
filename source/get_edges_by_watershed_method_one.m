function [ edge_locations, edges2vertices]                                                            ...
                = get_edges_by_watershed( linear_strel, first_max_edges_per_vertex, second_max_edges_per_vertex, max_edge_energy,             ...
                                                                       vertex_locations, energy_map )
%% initialize maps
% vertex_index_map has all the vertices marked with unique indices in order from best to worst
% energy. pointer_map has all the strels surrounding the vertices filled in with pointers pointed
% toward their centers. Avaliable energy map has just the vertices and their strels filled in.
number_of_vertices = length( vertex_locations );

if ~ issorted( energy_map( vertex_locations ))

    [ vertex_locations, original_vertex_indices ] = sort( energy_map( vertex_locations ));
    
else % issorted( energy_map( vertex_locations ))
    
    original_vertex_indices = 1 : number_of_vertices ;
    
end

% flip the vertex indices and add one just for this function, so that higher index is lower energy
backwards_indices = number_of_vertices : -1 : 1 ;

vertex_locations = vertex_locations( backwards_indices );

original_vertex_indices = [ number_of_vertices + 1, backwards_indices( original_vertex_indices )];  % background code concatenated on

border_index = 1 ;

numel_in_image = numel( energy_map );

size_of_image = int64( size( energy_map ));

cum_prod_image_dims = cumprod( size_of_image );

strel_apothem = max( round( linear_strel / cum_prod_image_dims( 2 )));

% %     vertex_index_map = zeros(  size_of_image,    'uint32'  );
%     vertex_index_map = zeros(  size_of_image,    'logical'  );
% %          pointer_map = zeros( numel_in_image, 1, 'uint8'    );
%          
% % available_energy_map = zeros( numel_in_image, 1, 'uint64'  );
% %     is_available_map = sparse( zeros( numel_in_image, 1, 'logical' )); %!!!!! doesnt pre-allocate memory
%     
% %     is_available_map = logical( spalloc( numel_in_image, 1, round( numel_in_image / 100 ))); %!!!!! doesnt pre-allocate memory
% 
% % sparse logical with true values at the vertex locations for starters and enough mem to hold one
% % thousandth of the image numel
% %     is_available_map = sparse( double( vertex_locations ), ones( size( vertex_locations )), ones( size( vertex_locations ), 'logical' ), numel_in_image, 1, round( numel_in_image / 10 ));
% 
% % available_energy_map = spalloc( numel_in_image, 1, numel_in_image / 2 );
% 
% % mark the outside of the vertex_index_map with a border idx to stop edges there.
% border_index  = number_of_vertices + 1 ;
% %  inert_border_index  = number_of_vertices + 2 ;
% 
% % !!!!!!!!!!!!!!!!!! this could be done more elegantly with modular arithmetic !!!!!!!!!!!!!!!!!!
% vertex_index_map([ 1, end ],  :, : ) = true ;
% vertex_index_map( :, [ 1, end ], : ) = true ;
% vertex_index_map( :, :,  [ 1, end ]) = true ;
% 
% border_locations = find( vertex_index_map );

% Linear indexing of the borders of the image:
border_locations_Y = [ 0 : strel_apothem                            - 1, cum_prod_image_dims( 1 ) - strel_apothem                            : cum_prod_image_dims( 1 ) - 1 ]' + ( 1 : cum_prod_image_dims( 1 ) : cum_prod_image_dims( 3 ));
border_locations_X = [ 0 : strel_apothem * cum_prod_image_dims( 1 ) - 1, cum_prod_image_dims( 2 ) - strel_apothem * cum_prod_image_dims( 1 ) : cum_prod_image_dims( 2 ) - 1 ]' + ( 1 : cum_prod_image_dims( 2 ) : cum_prod_image_dims( 3 ));
border_locations_Z = [ 0 : strel_apothem * cum_prod_image_dims( 2 ) - 1, cum_prod_image_dims( 3 ) - strel_apothem * cum_prod_image_dims( 2 ) : cum_prod_image_dims( 3 ) - 1 ]' + ( 1 : cum_prod_image_dims( 3 ) : cum_prod_image_dims( 3 ));

border_locations = unique([ border_locations_Y( : ); border_locations_X( : ); border_locations_Z( : )]);


% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% !!!! this may have memory issues to have the whole energy_map on the workspace.  Load in the
% energy map piecewise (on-demand in small windows) from the h5 file into a sparse energy_map.
% Erase the energy_map values when the available_map is turned false. The pre-sorting to facilitate
% the minimum energy search may need to be done on the fly at the time of loading from h5, or it may
% need to abandoned. (the minimum available energy value can be found on the fly)

% !!!! this should be done at time of reading in h52mat
% apply upper energy threshold
energy_map( energy_map > max_edge_energy ) = 0 ;
      
% linearize to match the other maps
% vertex_index_map = vertex_index_map( : );
      energy_map =       energy_map( : );
      
[ ~, sorting_indices ] = sort( energy_map );

sorting_indices = int64( sorting_indices );

unsorting_indices = ( 1 : uint64( numel_in_image ))'; unsorting_indices( sorting_indices ) = unsorting_indices ;

% % energy_map = uint64( sorting_indices );

% pointer_map is indices into the linear_strel. pointers show how to get back to the origin vertex,
% so the relative pointer indices are the negative (or reverse order) of the relative strel indices
% pointer_template = transpose( uint8( numel( linear_strel )) : - 1 : 1 );
new_pointers_template = transpose( numel( linear_strel ) : - 1 : 1 );

       vertex_index_map = sparse( double( border_locations ), 1, border_index,             numel_in_image, 1, round( numel_in_image / 10  ));

is_available_map_sorted = sparse( double( unsorting_indices( vertex_locations )), 1, true, numel_in_image, 1, round( numel_in_image / 20 ));
            pointer_map = sparse([ ], [ ], [ ],                                            numel_in_image, 1, round( numel_in_image / 20  ));
            
%% main WHILE
% number_of_vertices_found = 0 ;
index_of_last_vertex_found = number_of_vertices + 2 ; % one for background and one to subtract later

number_of_edges = uint32( 0 );

edges2vertices  = uint32([ ]);
edge_locations    =        { } ; % empty cell array

vertex_adjacency_matrix = logical( spalloc( number_of_vertices + 1, number_of_vertices + 1, ( number_of_vertices + 1 ) * second_max_edges_per_vertex ));

% has_second_max_edges = zeros( number_of_vertices + 1, 1, 'logical' );
has__first_max_edges = zeros( number_of_vertices + 1, 1, 'logical' );

strel_locations_bordering_index = cell( number_of_vertices + 1, 1 );

% WHILE searching for lowest available energy
while true

    % choose next index as the next available with lowest energy

%     % for unsorted energy map:
%     locations_available = find( is_available_map );
% 
%     [ min_available_energy, current_available_index ] = min( energy_map( locations_available ));
%     [ min_available_energy, current_available_index ] = min( energy_map( is_available_map_i ));
%     
%     current_location = locations_available( current_available_index );
    
    current_location = sorting_indices( find( is_available_map_sorted , 1 ));
    
    % requiring an available move
    if isempty( current_location ), break, end
    
    min_available_energy = energy_map( current_location );
    
    % requiring negative edge energy to explore    
    if min_available_energy >= 0, break, end
    
    current_location =  int64( current_location );    

    if pointer_map( current_location ) == 0
                
%         number_of_vertices_found = number_of_vertices_found + 1 ;
% 
%         current_vertex_index = number_of_vertices_found ;

        index_of_last_vertex_found = index_of_last_vertex_found - 1 ;
        
        current_vertex_index = index_of_last_vertex_found ;

    else % pointer_map( current_location ) == 0
        
        % record which vertex index the current location points to
        current_vertex_index = full( vertex_index_map( current_location + linear_strel( pointer_map( current_location ))));

    end
            
    % location will not be explored again
%     available_energy_map( current_location ) = 0 ;
    is_available_map_sorted( unsorting_indices( current_location  )) = false ;
           vertex_index_map(                    current_location )   = current_vertex_index ; % won't make any more edges and will make an inert shell around the spent waterhsed
        
%     sparse_entry_to_remove = find( is_available_map_i == current_location );
%            
%     is_available_map_i( sparse_entry_to_remove ) = [ ];
%     
%     vertex_index_map_i = [ vertex_index_map_i; current_location     ];
%     vertex_index_map_v = [ vertex_index_map_v; current_vertex_index ]; % won't make any more edges and will make an inert shell around the spent waterhsed
           
% %     is_available_map_i( is_available_map_i == current_location ) = [ ];
% %     is_available_map_c = is_available_map_c - 1 ;
% %     
% %     vertex_index_map_c = vertex_index_map_c + 1 ;
% %     vertex_index_map_i( vertex_index_map_c ) = current_location ;
% %     vertex_index_map_v( vertex_index_map_c ) = current_vertex_index ; % won't make any more edges and will make an inert shell around the spent waterhsed
% %     


	current_strel = linear_strel + current_location ;
    
    % look for adjacent watersheds
    vertices_in_strel = unique( nonzeros( vertex_index_map( current_strel )))';
    
%     [ ~, current_strel_at_vertices ] = intersect( vertex_index_map_i, current_strel );
%     
%     vertices_in_strel = vertex_index_map_v( current_strel_at_vertices );    

    if ~has__first_max_edges( current_vertex_index ) % vertex has not reached first max number of edges 
            
        % look for the strel/neighborhood around the current vertex that has not yet been explored
        is_new_in_strel = full(      pointer_map( current_strel )) == 0 ...
                        & full( vertex_index_map( current_strel )) == 0 ;
                    
%         current_strel_new = setdiff( current_strel, [ current_strel_at_vertices; pointer_map_i ]);
% 
%         is_new_in_strel = any( current_strel == current_strel_new' );

        current_strel_new = current_strel( is_new_in_strel );

        pointer_map( current_strel_new ) = new_pointers_template( is_new_in_strel );

%         pointer_map_i = [ pointer_map_i;         current_strel( is_new_in_strel )];
%         pointer_map_v = [ pointer_map_v; new_pointers_template( is_new_in_strel )];
        
        % fill in the strell with the available energy
    %     available_energy_map( current_strel_new ) = energy_map( current_strel_new );
        is_available_map_sorted( unsorting_indices( current_strel_new )) = true  ;

%         is_available_map_i = [ is_available_map_i;             current_strel( is_new_in_strel )              ];
        
% %         has_second_max_edges( vertex_index ) =    number_of_edges__at_vertex ...
% %                                                >= second_max_edges_per_vertex ;                

    else % current vertex has reached first max number of edges

        % look to see if this watershed or an adjacent one is done growing so we can clear up memory        
        
        vertices_in_strel_at_edge_max = vertices_in_strel( has__first_max_edges( vertices_in_strel ));
        
        for vertex_index = vertices_in_strel_at_edge_max

            % IF no locations bordering this index are still available (to this index or another)
            if ~ any( is_available_map_sorted( unsorting_indices( strel_locations_bordering_index{ vertex_index })))

                % find all locations sourrounded by this index

                locations_with_index = int64( find( vertex_index_map == vertex_index ));

          %       locations_with_index = vertex_index_map_i( vertex_index_map_v == current_vertex_index );        

                locations_with_index_and_pointer = locations_with_index( logical( pointer_map( locations_with_index )));

    %             locations_with_index_and_pointer = intersect( locations_with_index, pointer_map_i );

                if ~isempty( locations_with_index_and_pointer )

                    strels_with_index_and_pointer = locations_with_index_and_pointer' + linear_strel ;

                    vertex_index_strels = full( vertex_index_map( strels_with_index_and_pointer ));

                    vertices_at_current_vertex = permute( unique( edges2vertices( any( edges2vertices == vertex_index, 2 ), : )), [ 3, 2, 1 ]);

                    % clear all locations with pointers (not the found edges) whose strels contain only this index, or the other index in the edge, or border index
                    locations_with_pointer_surrounded_by_index                                           ...
                        = locations_with_index_and_pointer( all( any(    vertex_index_strels             ...
                                                                      == vertices_at_current_vertex, 3 )));

                    vertex_index_map( locations_with_pointer_surrounded_by_index ) = 0 ; % free up memory

                end % IF there are any locations with this index and a pointer
            end % IF vertex watershed has no available, adjacent energy
        end % FOR vertices in current strel
    end % IF vertex does not have max edges
    
    other_vertices_in_strel = vertices_in_strel( vertices_in_strel ~= current_vertex_index );

%     vertices_at_current_vertex = unique( edges2vertices( any( edges2vertices == current_vertex_index, 2 ), : );
%     
%     other_vertices_at_current_vertex = vertices_at_current_vertex( vertices_at_current_vertex ~= current_vertex_index );        
        
    % IF one or more other watersheds meet the current watershed for the first timw
    while ~all( vertex_adjacency_matrix( other_vertices_in_strel, current_vertex_index )) % not( all( empty )) is false
        
%     while ~isempty( setdiff( vertices_in_strel, vertices_at_current_vertex ))

        number_of_edges = number_of_edges + 1 ;

        % in case there are more than one waterhseds meeting, select the first on the list
%         best_vertex_in_others = find( ~ vertex_adjacency_matrix( current_vertex_index, other_vertices_in_strel ), 1 );
        best_vertex_in_others = find( ~ vertex_adjacency_matrix( current_vertex_index, other_vertices_in_strel ), 1, 'last' );
        
        other_vertex_index = other_vertices_in_strel( best_vertex_in_others );        

        for edge_half = 1 : 2

            % half #1 is the first watershed to get there (lower energy route)
            if edge_half == 1 

                tracing_location = current_location ; 

            else % edge_half == 2

                location_with_other_vertex_in_strel = current_strel( vertex_index_map( current_strel ) == other_vertex_index );

                % in case of multiple locations to choose, pick the one with the lowest energy
                [ ~, tracing_location_in_strel ] = min( energy_map( location_with_other_vertex_in_strel ));
                
                tracing_location = location_with_other_vertex_in_strel( tracing_location_in_strel );

            end

            tracing_ordinate = 0 ;

            % trace back until hitting a zero in the pointer index map, recording the linear indices
            % into the image and setting the pointer index map to zero along the way.  If this
            % process hits a zero index that doesn't belong to a vertex then record this for
            % later (we will need to add a new vertex and break up the existing edge into two).
            while true

                tracing_ordinate = tracing_ordinate + 1 ;

                edge_locations{ number_of_edges, edge_half }( tracing_ordinate ) = tracing_location ;

                if pointer_map( tracing_location ) == 0, break, end

                tracing_location = tracing_location + linear_strel( pointer_map( tracing_location ));

            end % WHILE tracing            
        end % FOR edge_half

        % Erase the pointers for the found trajectory, so that this path will not contribute to
        % any other edge from this vertex. Replace the pointers with the negative of the parent
        % edge identity that owns this path.
        pointer_map( edge_locations{ number_of_edges, 1 }) = 0 ;
        pointer_map( edge_locations{ number_of_edges, 2 }) = 0 ;
 
%         [ ~, sparse_entry_to_erase ] = intersect( pointer_map_i, [ edge_locations{ number_of_edges, 1 };   ...
%                                                                    edge_locations{ number_of_edges, 2 } ]); ...
%                                                                  
%         pointer_map_i( sparse_entry_to_erase ) = [ ];
%         pointer_map_v( sparse_entry_to_erase ) = [ ];
                                                                              
        % Record the identity of the origin vertex in the second position in edges2vertices
        % (corresponds to last entry in the list of positions).  Also record the terminal
        % vertex, first position.
        edges2vertices( number_of_edges, 1 ) = current_vertex_index ;
        edges2vertices( number_of_edges, 2 ) =   other_vertex_index ;       

        % label these vertices as adjacent on the adjacency matrix
        vertex_adjacency_matrix( current_vertex_index,   other_vertex_index ) = true ;
        vertex_adjacency_matrix(   other_vertex_index, current_vertex_index ) = true ;

        for vertex_index = [ other_vertex_index, current_vertex_index ] % favor current vertex, so remove other_vertex_index first

            if vertex_index ~= border_index % the image boundary has no limit to the number of edges it can create
                
%                 had_second_max_edges = has_second_max_edges( vertex_index );
%                 had__first_max_edges = has__first_max_edges( vertex_index );

%                 % IF this vertex never had the second maximum number of edges
%                 if ~had_second_max_edges
                    % calculate max's and see if you need to clean up the watershed to save memory
                    
                if ~has__first_max_edges( vertex_index )

                    number_of_edges__at_vertex = full( sum( vertex_adjacency_matrix( :, vertex_index )));                    

%                     number_of_edges__at_vertex = sum( edges2vertices( : ) == vertex_index );                    

                    has__first_max_edges( vertex_index ) =    number_of_edges__at_vertex ...
                                                           >= first_max_edges_per_vertex ;
                                                       
                    if has__first_max_edges( vertex_index )

                        locations_with_index = int64( find( vertex_index_map == vertex_index ));

                %       locations_with_index = vertex_index_map_i( vertex_index_map_v == vertex_index );        

                        strel_locations_bordering_index{ vertex_index } = setdiff( locations_with_index' + linear_strel( linear_strel ~= 0 ), locations_with_index );

                    end % IF vertex has first max now
                end % IF vertex did not have first max before
%                 end % IF vertex never had second max
            end % IF vertex isn't code for border
        end % FOR vertex_index
    end % WHILE joining watershed 
end % WHILE searching for edges

%% organizing output

% reassign vertex indices to their input, flipped and possibly unsorted values    
edges2vertices = original_vertex_indices( edges2vertices ); 

% vertex_index_map( vertex_index_map> 0 ) = vertex_indices( vertex_index_map( vertex_index_map > 0 )); 

% % visulize maps
mat2tif( reshape( full( vertex_index_map ), size_of_image ), 'vertex_index_map.tif' )
mat2tif( reshape( full( is_available_map_sorted ), size_of_image ), 'is_available_map_sorted.tif' )
mat2tif( reshape( full( is_available_map_sorted( unsorting_indices )), size_of_image ), 'is_available_map.tif' )
mat2tif( reshape( full( pointer_map ), size_of_image ), 'pointer_map.tif' )

% flip all of the first halves of edges so that they start with the vertex and end at the pass
edge_locations = cellfun( @( h1, h2 ) [ fliplr( h1 ), h2 ]', edge_locations( :, 1 ), edge_locations( :, 2 ), 'UniformOutput', false );

end % FUNCTION get_edges_by_watershed