function [ edges2vertices, edge_lengths, edge_space_subscripts, edge_scale_subscripts, edge_energies ]  ...
                     = get_edges_V204(  lumen_radius_in_microns_range, microns_per_voxel, ...
                                        vertex_scale_subscripts, vertex_space_subscripts, ...
                                        strel_apothem, max_edge_length_per_origin_radius, ...
                                              number_of_edges_per_vertex, data_directory, ...
                                                                 energy_handle, use_case  )
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

[ vertex_center_image, vertex_locations, reading_box_apothem, linear_strel,  ...
                   max_edge_length_in_microns_range, numel_of_strel,  ...
            cum_prod_image_dims, size_of_image, path_to_energy_data, strel_distance_LUT ] ...
  = generate_reference_image(  lumen_radius_in_microns_range, microns_per_voxel, ...
                               vertex_scale_subscripts, vertex_space_subscripts, ...
                               strel_apothem, max_edge_length_per_origin_radius, ...
                                                  data_directory, energy_handle  );
                                                                                                        
number_of_vertices  = length( vertex_scale_subscripts ); 

% number_of_edges     = number_of_vertices * number_of_edges_per_vertex ;

% lengths_of_edges      = cell( number_of_vertices, 1 );
edges2vertices        = cell( number_of_vertices, 1 );
edge_indices_temp     = cell( number_of_vertices, 1 );

%      lengths_of_edges{ : } = zeros( number_of_edges_per_vertex, 1, 'uint16' ); % number of positions visited (units of voxels)
%        edges2vertices{ : } = zeros( number_of_edges_per_vertex, 2, 'uint32' ); % the indices of the vertices connected by each edge
% edge_indices{ : } = zeros( number_of_edges_per_vertex, vertex_length_maximum, 'uint16' );

% edge_indices{ : } = cell( number_of_edges_per_vertex, 1 );

% switch use_case
%     
%     case 'standard'

vertex_unique_range = uint32( 1 : number_of_vertices );
        
%     case { 'add_vertex_to_edge', 'extend_dead_end_edge' }
%         
%         vertex_unique_range = varargin{ 1 };

use_case = 'get_edges' ;

%% main PARFOR: loop through the vertices, searching for nearby vertices to connect to
for vertex_index = vertex_unique_range

    [ edge_indices_temp{ vertex_index }, edges2vertices{ vertex_index }]                    ...
             = get_edges_for_vertex([ ], [ ], [ ], [ ], microns_per_voxel, strel_apothem,   ...
                                     number_of_edges_per_vertex, use_case,                  ...
                                                     vertex_locations( vertex_index ),      ...
                  vertex_center_image, vertex_locations, reading_box_apothem, linear_strel, ...
                                   max_edge_length_in_microns_range, numel_of_strel, strel_distance_LUT,    ...
                            cum_prod_image_dims, size_of_image, path_to_energy_data         );
                        
%                         any ( edge_indices_temp{ vertex_index }( : ) == 20667639 )
%   !!!  any( any ( edge_indices_temp{ vertex_index }( : ) == [ 7728, 25346, 27653 ] ))
%     if any( any ( edge_indices_temp{ vertex_index }( : ) == [  1863425, 2125058 ] ))
    if any( any ( edge_indices_temp{ vertex_index }( : ) ==     560105 ))
%         %
%         
    end

end % PARFOR vertex_index

%% Organizing outputs:

[ edges2vertices, edge_lengths, edge_space_subscripts, edge_scale_subscripts, edge_energies ]       ...
    = get_edge_vectors( number_of_edges_per_vertex, edge_indices_temp, edges2vertices, path_to_energy_data, cum_prod_image_dims );

end % FUNCTION

