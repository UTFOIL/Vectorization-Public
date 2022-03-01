function [ edge_space_subscripts, edge_scale_subscripts, edge_energies, edges2vertices ] ...
             = clean_edges( edge_space_subscripts, edge_scale_subscripts, edge_energies, edges2vertices )
%% choose_edges_V160
                                   
% this function selects the "best" edges of the possible edges found in the get_edges probabalistic
% walk through the 4D (projected to 3D acrosss scale) energy field created in the get_energy
% function.  For a given vertex pair A to B, the best trajectory connecting A to B is called the one
% whose maximal energy attained along the trajectory is the least among all the trajectories from A
% to B.  SAM 5/14/18
%
% V181: the previous version is combined with a function that resembles choose_vertices_V184 in that
% it paints edges and looks at the painted image to determine conflicts.  Edges conflict with other
% edges whenever they overlap away from the vertices that they connect.  All edges paint the same
% Inf value. For each edge, after it paints itself it also paints the index of its start vertex
% around that vertex and the index of its end vertex around its end.  To avoid conflict between
% adjacent edges targeting the same node, the vertex influence should be larger than the edge
% influence. SAM 5/30/18
%
% V182: edges don't paint, all that is painted is the vertices they start and stop on. The edges
% still have an influence volume for searching for violations of vertices, but this volume could be
% the same as that used for the vertices. SAM 5/31/18
%
% V190: mean edge energy is used instead of max edge energy as the summary statistic.  8/2/18 SAM
%
% function name changed from choose_edges_V190 to crop_edges_V200.  Keeping only the part of the
% choose edges function that does the cropping and the part that chooses best unique trajectories
% between vertices A and B mutually. eliminating references to sigma; instead referring to lumen
% radius.  Splitting edge subscripts into space (uint16) and scale (uint8) cell vectors.  SAM
% 11/14/18

%% Choosing best unique trajectories between vertices A and B mutually.

[ edges2vertices, original_edge_indices ] = clean_edge_pairs( edges2vertices, edge_energies );

edge_space_subscripts  =  edge_space_subscripts( original_edge_indices );
edge_scale_subscripts  =  edge_scale_subscripts( original_edge_indices );
edge_energies          =          edge_energies( original_edge_indices );

%% scratch

% %% Remove edges whose immediate parent no longer exist. 
% % Remove any edge whose endpoints according to the vertices (edges2vertices matrix) disagree with
% % its endpoints from the edge objects.  (This removes any edge who failed to go from a child to
% % having a vertex of its own, which means its immediate parent was lost along the way).
% 
% end_points_from_vertices = vertex_locations( edges2vertices );
% 
% end_points_from_edges = cell2mat( cellfun( @( x ) x([ 1, end ]), edge_locations_cell, 'UniformOutput', false ));
% 
% % edge is deleted if one or both endpoints disagree
% is_edge_deleted = logical( sum( end_points_from_vertices ~= end_points_from_edges, 2 ));
% 
%        edges2vertices( is_edge_deleted, : ) = [ ];
% edge_space_subscripts( is_edge_deleted    ) = [ ];
% edge_scale_subscripts( is_edge_deleted    ) = [ ];
%         edge_energies( is_edge_deleted    ) = [ ];
% edge_space_subscripts( is_edge_deleted    ) = [ ];
% 
% % vertex_space_subscripts_double = uint16( vertex_space_subscripts_double );

end % function
