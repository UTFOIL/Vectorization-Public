function [ original_edge_indices ] = clean_edges_orphans( edge_space_subscripts, size_of_image, vertex_space_subscripts )
%% clean_edge_cycles
% SAM 7/1/19

number_of_edges = size( edge_space_subscripts, 1 );

original_edge_indices = 1 : number_of_edges ;

%% Remove child (grandchild, (etc.)) edges who do not connect to any other edge
% This section function should be called before the cycle-breaking section, if edges were erased
% before that section.  In that case we would want to eliminate any orphaned connecitons from the
% adjacency matrix before removing cycles and possibly removing non-orphaned edges before an
% orphaned, destined to be removed edge.
% 
% Note: In the add vertices section we again force children to have higher (worse) energy than their
% parents.

edge_space_subscripts_matrix = double( cell2mat( edge_space_subscripts ));

vertex_space_subscripts_double = double( vertex_space_subscripts );

edge_locations   = sub2ind( size_of_image,   edge_space_subscripts_matrix( :, 1 ),   edge_space_subscripts_matrix( :, 2 ),   edge_space_subscripts_matrix( :, 3 ));
vertex_locations = sub2ind( size_of_image, vertex_space_subscripts_double( :, 1 ), vertex_space_subscripts_double( :, 2 ), vertex_space_subscripts_double( :, 3 ));

edge_lengths = cellfun( @( x ) size( x, 1 ), edge_space_subscripts );

edge_locations_cell = mat2cell( edge_locations, edge_lengths, 1 );

searching_for_orphans = true ;

while searching_for_orphans
    
    edge_locations = cell2mat( edge_locations_cell );
    
%     % unstable unique, because all these get deleted anyway
%     [ ~, unique_edge_location_indices ] = unique( edge_locations );
% 
%     % throw out all of the unique ones and save the repeats
%     edge_location_repeats = edge_locations ; edge_location_repeats( unique_edge_location_indices ) = [ ];
    
    % redefine the edge look up table because the number of edge list numbering has changed
    number_of_edges       = length( edge_locations_cell );    
    
    edge_index_range_cell = num2cell( 1 : number_of_edges )';
    
    edge_index_LUT = cellfun( @( x, y ) x * ones( size( y )), edge_index_range_cell,      edge_locations_cell, 'UniformOutput', false  );    
    
    interior_edge_locations                 = cell2mat( cellfun( @( x ) x( 2 : end - 1 ), edge_locations_cell, 'UniformOutput', false ));    
    exterior_edge_locations                 = cell2mat( cellfun( @( x ) x([ 1, end ]   ), edge_locations_cell, 'UniformOutput', false ));
    exterior_edge_location_index2edge_index = cell2mat( cellfun( @( x ) x([ 1, end ]   ), edge_index_LUT     , 'UniformOutput', false ));  
    
    
%     % only the last location of an edge can be disjoint from any vertex location
%     exterior_edge_locations                 = cell2mat( cellfun( @( x ) x( end ), edge_locations_cell, 'UniformOutput', false ));
%     exterior_edge_location_index2edge_index = cell2mat( cellfun( @( x ) x( end ), edge_index_LUT     , 'UniformOutput', false ));    
%     % no longer true, since get_edges_by_watershed function. SAM 9/23/20

    % look for edge terminal locations that don't coincide with a vertex or any other edge location.
    % (no repeats returnded by setdiff, but no repeats sought at this line (repeated locations are
    % subtracted out by setdiff ))
    
%         [ ~, orphan_terminal_indices ] = setdiff( exterior_edge_locations( : ), union( edge_location_repeats, vertex_locations ));
    [ ~, orphan_terminal_indices ] = setdiff( exterior_edge_locations( : ), union( interior_edge_locations, vertex_locations ));
    
    edge_indices_to_remove = exterior_edge_location_index2edge_index( orphan_terminal_indices );
        
    % erase these edges
    original_edge_indices( edge_indices_to_remove    ) = [ ];
      edge_locations_cell( edge_indices_to_remove    ) = [ ];
%            edges2vertices( edge_indices_to_remove, : ) = [ ];
%     edge_space_subscripts( edge_indices_to_remove    ) = [ ];      
%     edge_scale_subscripts( edge_indices_to_remove    ) = [ ];
%             edge_energies( edge_indices_to_remove    ) = [ ];
            
    searching_for_orphans = logical( length( edge_indices_to_remove ));
    
end % WHILE searching for orphans (no parent)

% %% remove child edges who only connect to other children or who have better energy than their possible parents
% 
% % !!!!!!!!!!!!!!!!!!!! WHILE WHILE WHILE
% 
% % look for edge external locations that don't coincide with a vertex (no repeats returned by
% % setdiff)
% child_external_edge_locations = setdiff( exterior_edge_locations( : ), vertex_locations )';
% 
% % find the repeats for the repeats for these child termini (that are disjoing from vertex locations)
% child_external_edge_indices =  [ ];
% 
% for child_external_edge_location = child_external_edge_locations
%     
%     child_external_edge_indices = [ child_external_edge_indices, find( exterior_edge_locations( : ) == child_external_edge_location )' ];
%     
% end
% 
% edge_indices_for_children = exterior_edge_location_index2edge_index( child_external_edge_indices )';
% 
% is_child_edge_orphan = zeros( size( original_edge_indices ), 'logical' );
% 
% for child_edge_index = edge_indices_for_children
%     
%     child_parent_meeting_location = edge_locations_cell{ child_edge_index }( end );
%     
% %     if child_parent_meeting_location == 26028996
% %         
% %         pause
% %         
% %     end
%     
%     edges_coinciding_with_child = find( cellfun( @( x ) any( child_parent_meeting_location == x ), edge_locations_cell ));
%     
%     parents_coinciding_with_child = setdiff( edges_coinciding_with_child, edge_indices_for_children );
%     
%     if isempty( parents_coinciding_with_child )
%                 
%         is_child_edge_orphan( child_edge_index ) = true ;
%         
%     else
%     
%         best_parent_coinciding_with_child = parents_coinciding_with_child( 1 );
% 
%         is_child_edge_orphan( child_edge_index ) = child_edge_index < best_parent_coinciding_with_child ;
%         
%     end
%     
% end % FOR searching for orphans (parent of higher energy)
% 
% original_edge_indices( is_child_edge_orphan ) = [ ];

end % FUNCTION clean_edges_orphans
