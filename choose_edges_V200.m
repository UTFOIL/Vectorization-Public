function [ mean_edge_energies, original_edge_indices, edges2vertices, painted_image ]         ...
                  = choose_edges_V200( edges2vertices, edge_energies, edge_space_subscripts, edge_scale_subscripts, ...
                                       pixels_per_sigma_range,                               ...
                                       sigma_per_influence_vertices,                         ...
                                       sigma_per_influence_edges, size_of_image             )
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
% V200, Splitting edge subscripts into space (uint16) and scale (uint8) cell vectors.

%remember the indices of the original vertices for later reference

edge_subscripts = cellfun(      @(      space_subscripts,      scale_subscripts  )  ...
                           round([      space_subscripts,      scale_subscripts ]), ...
                                   edge_space_subscripts, edge_scale_subscripts,    ...
                                         'UniformOutput', false                     );

number_of_original_edges = length( edges2vertices );

original_edge_indices = 1 : number_of_original_edges ;

% remove trajectories that terminate on their original vertex
indices_of_non_self_referential_edges = find( edges2vertices( :, 1 ) ~= edges2vertices( :, 2 ));

edges2vertices          =          edges2vertices( indices_of_non_self_referential_edges, : );
% edge_sizes_temp       =       edge_subscripts( indices_of_non_self_referential_edges, 4 );
edge_energies_temp    =         edge_energies( indices_of_non_self_referential_edges    );
original_edge_indices = original_edge_indices( indices_of_non_self_referential_edges    );

% remove trajectories that don't land on any vertex
indices_of_terminal_edges = find( edges2vertices( :, 2 ) > 0 );

edges2vertices          =          edges2vertices( indices_of_terminal_edges, : );
% edge_sizes_temp       =       edge_sizes_temp( indices_of_terminal_edges    );
edge_energies_temp    =    edge_energies_temp( indices_of_terminal_edges    );
original_edge_indices = original_edge_indices( indices_of_terminal_edges    );

% remove edges that ever attained a negative energy on the trajectory
max_edge_energies = cellfun( @max,  edge_energies_temp );

indices_of_negative_energy_edges = find( max_edge_energies < 0 );

% summarizing each trajectory by its mean energy 
% mean_edge_energies = cellfun( @mean,  edge_energies_temp );

[ mean_edge_energies ] = get_edge_metric( edge_energies_temp );

% 10/4/18 Possible Improvement: weight the edge means by a factor that is low when a certain
% variance in the size coordinate of each trajectory is large.  This variance is first adjusted so
% that any variance explained by a linear fit is removed:
%
% characteristic_size_stdev = scales_per_octave ; % needs to be input to this function
%
% % subtract off the constant term
% edge_sizes_temp = cellfun( @(v) v - mean(v), edge_sizes_temp, 'UniformOutput', false )
%
% subtract off the (orthogonal) linear fits through the origin 
% edge_sizes_temp = cellfun( @(v) v -       ((1:length(v))-(length(v)-1)/2)'    ...
%                                     *     ((1:length(v))-(length(v)-1)/2)     ...
%                                     * v/( ((1:length(v))-(length(v)-1)/2)     ...
%                                          *((1:length(v))-(length(v)-1)/2)' ), ...
%                             edge_sizes_temp, 'UniformOutput', false           )
% 
% % calculate the variances of the sizes, take negative exponent, and multiply this factor by the
% % mean energy
% mean_edge_energies = cellfun( @mean,  edge_energies_temp ) .* exp( - cellfun( @stdev, edge_sizes_temp ) / characteristic_size_stdev, edge_sizes_temp ));

mean_edge_energies    =     mean_edge_energies( indices_of_negative_energy_edges    );
original_edge_indices = original_edge_indices(  indices_of_negative_energy_edges    );
edges2vertices          =   edges2vertices(         indices_of_negative_energy_edges, : );

% sorting the trajectories by lowest activation energy before choosing unique edges     
% SORT sorts in ascending order by default, so lowest energies at the top of the list.
[ mean_edge_energies, indices_sorted_by_mean_energy ] = sort( mean_edge_energies );

original_edge_indices = original_edge_indices(  indices_sorted_by_mean_energy    );
edges2vertices          =   edges2vertices(         indices_sorted_by_mean_energy, : );

% choosing only the best trajectory from vertex A to B for all A, B in the set of all vertices.
% UNIQUE chooses the first of the multiple instances, so top of the list/lowest energy is preferred.
[ edges2vertices, indices_of_unique_edges ] = unique( edges2vertices, 'rows', 'stable' );

original_edge_indices = original_edge_indices(  indices_of_unique_edges );
mean_edge_energies    =    mean_edge_energies( indices_of_unique_edges );

% choosing the better of the two directions of edges in the cases where we have pairs of mutual
% trajectories (from A to B and B to A)
[ ~, mutual_edges_sorted_indices, mutual_edges_indices_antiparallel_to_sorted ] ...
    = intersect([ edges2vertices( :, 1 ), edges2vertices( :, 2 )],  ...
                [ edges2vertices( :, 2 ), edges2vertices( :, 1 )],  ...
                'rows', 'stable'                                );

logical_of_worse_of_mutual_edge_pairs                                           ...
    = mutual_edges_sorted_indices > mutual_edges_indices_antiparallel_to_sorted ;  

original_edge_indices( mutual_edges_sorted_indices( logical_of_worse_of_mutual_edge_pairs     )) = [ ];
   mean_edge_energies( mutual_edges_sorted_indices( logical_of_worse_of_mutual_edge_pairs     )) = [ ];
       edges2vertices( mutual_edges_sorted_indices( logical_of_worse_of_mutual_edge_pairs ), : ) = [ ];      
         
% end % function
% 

edge_energies    =    edge_energies( original_edge_indices );
edge_subscripts  =  edge_subscripts( original_edge_indices );
% degrees_of_edges( original_edge_indices );
degrees_of_edges = cellfun( @( x ) size( x, 1 ), edge_subscripts );

% function [ original_vertex_indices, painted_image ]                                            ...
%             = choose_vertices_V184( space_subscripts, scale_subscripts, energy_values,         ...
%                                     pixels_per_sigma_range, sigma_per_influence, size_of_image )
%% Choose vertices_V180
% % the purpose of this function is to clean up the output of the get_vertices function. It does this
% % by sorting the vertices from best to worst (highest to lowest) energy and then looping through and
% % eliminating vertices whose volumes overlap with vertices of lower energy.
% %
% % SAM 5/29/18
% %
% % V181 in which the elimination is one sided: volume overlap isn't sufficient to eliminate vertices.
% % We look to see if the current vertex's center (not volume) is inside the volume of any previously
% % accepted vertex. SAM 5/29/18
% %
% % V182 is somewhere between V180 and V181 SAM 5/29/18, the inner sigma volumes need to overlap (not
% % the sizes as in the case of V180.
% %
% % V183 in which the size of the volumes needed to overlap to exclude a vertex is input as a factor
% % of sigma. SAM 5/29
% %
% % V184 where instead of searching for previously identified vertices that will interfere with the
% % current object, the method is to look at a painted image of previously identified vertices where
% % previous vertices have been painted. SAM 5/30/18

%% construct the relative sphere elements at all the scales


number_of_scales = size( pixels_per_sigma_range, 1 );

scale_subscript_range = 1 : number_of_scales ;

% ---------------------------------------- edge objects -------------------------------------------

edge_element_subscripts      = cell( number_of_scales, 1 );
edge_element_linear_indexing = cell( number_of_scales, 1 );

for scale_subscript = scale_subscript_range
    
    radii_in_pixels_at_scale = sigma_per_influence_edges                  ...
                             * pixels_per_sigma_range( scale_subscript, : );
    
    [ edge_element_linear_indexing{ scale_subscript },                                 ...
      edge_element_subscripts{      scale_subscript }, ]                               ...
        = construct_structuring_element( radii_in_pixels_at_scale, size_of_image  );
 
end % constructing relative elements FOR scale

% --------------------------------------- vertex objects ------------------------------------------

vertex_element_subscripts      = cell( number_of_scales, 1 );
vertex_element_linear_indexing = cell( number_of_scales, 1 );

for scale_subscript = scale_subscript_range
    
    radii_in_pixels_at_scale = sigma_per_influence_vertices               ...
                             * pixels_per_sigma_range( scale_subscript, : );
    
    [ vertex_element_linear_indexing{ scale_subscript },                               ...
      vertex_element_subscripts{      scale_subscript }, ]                             ...
            = construct_structuring_element( radii_in_pixels_at_scale, size_of_image  );
 
end % constructing relative elements FOR scale

%% begin the painting procedure

% enumerate the original edges for purposes of outputting the indices of the chosen edges
number_of_edges = length( edge_energies );

intermediary_edge_indices = 1 : number_of_edges ;

% sort the edges from best to worst in terms of energy value (lower energy is better)
[ ~, sorted_edge_indices ] = sort( mean_edge_energies );

intermediary_edge_indices_sorted = intermediary_edge_indices( sorted_edge_indices );

% initialize the vector to assign logical 1 or 0 to whehter each vertex is chosen
chosen_edge_logical = zeros( number_of_edges, 1, 'logical' );

% intialize the painted_image with a blank canvas
painted_image = zeros( size_of_image, 'uint32' );

% predict which templates will need to be trimmed so as to not wrap around or go over the last index

%   ones_array =     ones( number_of_edges, 1 );
  twos_array = 2 * ones( number_of_edges, 1 );
% threes_array = 3 * ones( number_of_edges, 1 );
%  fours_array = 4 * ones( number_of_edges, 1 );

% ---------------------------------------- edge objects -------------------------------------------

% pull the edge subscripts out of the cell array and concatenate into matrix
edge_subscripts_matrix = cell2mat( edge_subscripts );

radii_in_pixels                                                            ...
    = round((   sigma_per_influence_edges                                  ...
              * pixels_per_sigma_range( edge_subscripts_matrix( :, 4 ), : )));

subscript_maxs = edge_subscripts_matrix( :, 1 : 3 ) + radii_in_pixels ;
subscript_mins = edge_subscripts_matrix( :, 1 : 3 ) - radii_in_pixels ;

y_is_over  = subscript_maxs( :, 1 ) > size_of_image( 1 );
x_is_over  = subscript_maxs( :, 2 ) > size_of_image( 2 );
z_is_over  = subscript_maxs( :, 3 ) > size_of_image( 3 );

y_is_under = subscript_mins( :, 1 ) <         1         ;
x_is_under = subscript_mins( :, 2 ) <         1         ;
z_is_under = subscript_mins( :, 3 ) <         1         ;

% unforseen bug results from passing uint16 degrees of edges to mat2cell because it accumulates in
% the same data type and maxes out before reaching the sum of the degrees of edges

degrees_of_edges_uint_32 = uint32( degrees_of_edges );

y_is_over  = mat2cell( y_is_over,  degrees_of_edges_uint_32, 1 );
x_is_over  = mat2cell( x_is_over,  degrees_of_edges_uint_32, 1 );
z_is_over  = mat2cell( z_is_over,  degrees_of_edges_uint_32, 1 );

y_is_under = mat2cell( y_is_under, degrees_of_edges_uint_32, 1 );
x_is_under = mat2cell( x_is_under, degrees_of_edges_uint_32, 1 );
z_is_under = mat2cell( z_is_under, degrees_of_edges_uint_32, 1 );

% index_code_for_edge = max( chosen_edges( : )) + 1 ;

% --------------------------------------- vertex objects ------------------------------------------

% get the first and last positions from the edge since these are the start and end vertices
% corresponding to the first and second columns (respectively) of the chosen_edges variable
vertex_subscripts = cellfun( @( v ) v([ 1, end ], : ), edge_subscripts, 'UniformOutput', false );

vertex_subscripts_matrix = cell2mat( vertex_subscripts );

radii_in_pixels                                                              ...
    = round((   sigma_per_influence_vertices                                 ...
              * pixels_per_sigma_range( vertex_subscripts_matrix( :, 4 ), : )));

subscript_maxs = vertex_subscripts_matrix( :, 1 : 3 ) + radii_in_pixels ;
subscript_mins = vertex_subscripts_matrix( :, 1 : 3 ) - radii_in_pixels ;

y_is_over_vertices  = subscript_maxs( :, 1 ) > size_of_image( 1 );
x_is_over_vertices  = subscript_maxs( :, 2 ) > size_of_image( 2 );
z_is_over_vertices  = subscript_maxs( :, 3 ) > size_of_image( 3 );

y_is_under_vertices = subscript_mins( :, 1 ) <         1         ;
x_is_under_vertices = subscript_mins( :, 2 ) <         1         ;
z_is_under_vertices = subscript_mins( :, 3 ) <         1         ;

y_is_over_vertices  = mat2cell( y_is_over_vertices,  twos_array, 1 );
x_is_over_vertices  = mat2cell( x_is_over_vertices,  twos_array, 1 );
z_is_over_vertices  = mat2cell( z_is_over_vertices,  twos_array, 1 );

y_is_under_vertices = mat2cell( y_is_under_vertices, twos_array, 1 );
x_is_under_vertices = mat2cell( x_is_under_vertices, twos_array, 1 );
z_is_under_vertices = mat2cell( z_is_under_vertices, twos_array, 1 );

% vertex_subscripts   = mat2cell( vertex_subscripts_matrix, twos_array, 4 );

vertex_index_range = 1 : 2 ;

%% loop through the edges from best to worst (lowest to highest energy).
for edge_index = intermediary_edge_indices_sorted

    % zero out the spheres of influence of the two terminal vertices
    for vertex_index = vertex_index_range

        voxel_index = sub2ind( size_of_image, vertex_subscripts{ edge_index }( vertex_index, 1 ), ...
                                              vertex_subscripts{ edge_index }( vertex_index, 2 ), ...
                                              vertex_subscripts{ edge_index }( vertex_index, 3 )  );

        % shave off any overhang of the sphere element over the edge of the image
        vertex_linear_indexing{ vertex_index }                                                     ...
            = voxel_index                                                                          ...
            + vertex_element_linear_indexing{ vertex_subscripts{ edge_index }( vertex_index,   4   )};

        vertex_position_subscripts                                                                   ...
            =                                 vertex_subscripts{ edge_index }( vertex_index, 1 : 3 ) ...
            + vertex_element_subscripts{      vertex_subscripts{ edge_index }( vertex_index,   4   )};

        if y_is_over_vertices{  edge_index }( vertex_index ), vertex_linear_indexing{ vertex_index }( vertex_position_subscripts( :, 1 ) > size_of_image( 1 )) = voxel_index ; end
        if x_is_over_vertices{  edge_index }( vertex_index ), vertex_linear_indexing{ vertex_index }( vertex_position_subscripts( :, 2 ) > size_of_image( 2 )) = voxel_index ; end
        if z_is_over_vertices{  edge_index }( vertex_index ), vertex_linear_indexing{ vertex_index }( vertex_position_subscripts( :, 3 ) > size_of_image( 3 )) = voxel_index ; end

        if y_is_under_vertices{ edge_index }( vertex_index ), vertex_linear_indexing{ vertex_index }( vertex_position_subscripts( :, 1 ) <         1         ) = voxel_index ; end
        if x_is_under_vertices{ edge_index }( vertex_index ), vertex_linear_indexing{ vertex_index }( vertex_position_subscripts( :, 2 ) <         1         ) = voxel_index ; end
        if z_is_under_vertices{ edge_index }( vertex_index ), vertex_linear_indexing{ vertex_index }( vertex_position_subscripts( :, 3 ) <         1         ) = voxel_index ; end

    end % start and stop vertices FOR
    
	image_snapshot_at_vertices = painted_image( cell2mat( vertex_linear_indexing( : )));

	painted_image( cell2mat( vertex_linear_indexing( : ))) = 0 ;    
    
    edge_position_index_range = uint16( randperm( degrees_of_edges( edge_index )));    
    
%     edge_position_linear_indexing_cell = cell( degrees_of_edges( edge_index ), 1 );
    
    % assume the edge is good unless shown otherwise
    chosen_edge_logical( edge_index ) = true ;

    % loop through the positions along each edge
    for edge_position_index = edge_position_index_range

        voxel_index = sub2ind( size_of_image, edge_subscripts{ edge_index }( edge_position_index, 1 ), ...
                                              edge_subscripts{ edge_index }( edge_position_index, 2 ), ...
                                              edge_subscripts{ edge_index }( edge_position_index, 3 )  );

        % shave off any overhang of the sphere element over the edge of the image
        edge_position_linear_indexing                                                                 ...
            = voxel_index                                                                             ...
            + edge_element_linear_indexing{ edge_subscripts{ edge_index }( edge_position_index,   4   )};

        edge_position_subscripts                                                                        ...
            =                               edge_subscripts{ edge_index }( edge_position_index, 1 : 3 ) ...
            + edge_element_subscripts{      edge_subscripts{ edge_index }( edge_position_index,   4   )};

        if y_is_over{  edge_index }( edge_position_index ), edge_position_linear_indexing( edge_position_subscripts( :, 1 ) > size_of_image( 1 )) = voxel_index ; end
        if x_is_over{  edge_index }( edge_position_index ), edge_position_linear_indexing( edge_position_subscripts( :, 2 ) > size_of_image( 2 )) = voxel_index ; end
        if z_is_over{  edge_index }( edge_position_index ), edge_position_linear_indexing( edge_position_subscripts( :, 3 ) > size_of_image( 3 )) = voxel_index ; end

        if y_is_under{ edge_index }( edge_position_index ), edge_position_linear_indexing( edge_position_subscripts( :, 1 ) <         1         ) = voxel_index ; end
        if x_is_under{ edge_index }( edge_position_index ), edge_position_linear_indexing( edge_position_subscripts( :, 2 ) <         1         ) = voxel_index ; end
        if z_is_under{ edge_index }( edge_position_index ), edge_position_linear_indexing( edge_position_subscripts( :, 3 ) <         1         ) = voxel_index ; end

%         edge_position_linear_indexing_cell{ edge_position_index } = edge_position_linear_indexing ;
        
        % check if this area has already been painted (should only have zeros and potentially the start
        % and end vertices on the canvas here, otherwise this real estate is taken).
        conflicting_object_indices                                     ...
            = setdiff( painted_image( edge_position_linear_indexing ), ...
                       edges2vertices( edge_index, : )                   );

        % if there are any nonzero values in the conflicting objects set then we don't chose this
        % edge, we don't need to continue the FOR loop
        if any( conflicting_object_indices )
            
            chosen_edge_logical( edge_index ) = false ;
            
            break
            
        end        
    end % edge position FOR

    % If the edge is a keeper, paint the image in its volume of influence.
    if chosen_edge_logical( edge_index )
        
%         edge_entire_linear_indexing = cell2mat( edge_position_linear_indexing_cell );
%         
%         painted_image( edge_entire_linear_indexing ) = index_code_for_edge ;
        
        % also paint the vertices at either end with their indices        
        for vertex_index = vertex_index_range
        
            painted_image( vertex_linear_indexing{ vertex_index }) = edges2vertices( edge_index, vertex_index );

        end % start and stop vertices FOR
    else % edge is not keeper

        % put the image back how it was        
        painted_image( cell2mat( vertex_linear_indexing( : ))) = image_snapshot_at_vertices ;
        
    end % IF edge is chosen
end % sorted edge FOR

original_edge_indices = original_edge_indices( chosen_edge_logical    );
mean_edge_energies    =    mean_edge_energies( chosen_edge_logical    );
edges2vertices        =        edges2vertices( chosen_edge_logical, : );

end % FUNCTION
