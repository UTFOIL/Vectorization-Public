function [ edges2vertices, edge_numels, edge_space_subscripts, edge_scale_subscripts, edge_energies ]  ...
                     = get_edges_V300(  microns_per_voxel, lumen_radius_in_microns_range, vertex_space_subscripts, vertex_scale_subscripts, ...
                                        step_length_per_radius, max_edges_per_vertex, max_edge_energy, data_directory, visual_data_directory, energy_handle, edges_handle  )
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
%
% V300: no more walking from each vertex in parallael, but now it resembles a watershed algorithm
% and the entire image is operated on at once. SAM 8/14/20

%% hard-coded inputs:
% % all tolerances except energy are soft-thresholds: standard deviations of gaussian-weighted energy down-scaling factors
% % energy tolerance is a hard threshold after factoring in the soft-thresholds

max_voxels_per_node = 1e8 ;
% max_voxels_per_node = 1e6 ;
% max_voxels_per_node = 1e5 ;

% distance_tolerance_per_origin_radius = 5 ; % 10 % but this is half the distance because of symmetry
% distance_tolerance_per_origin_radius = 10 ; % 211020
% distance_tolerance_per_origin_radius = 1 ; % SAM 10/29/21
% distance_tolerance_per_origin_radius = 1.5 ; % SAM 10/29/21 16:40
% distance_tolerance_per_origin_radius = 2 ; % SAM 11/2/21, 11/7/21
distance_tolerance_per_origin_radius = 3 ; % SAM 1/11/22
% edge_number_tolerance = 2.01 ; % number of edges per vertex !!!!!!!!! step function, not gaussian
% edge_number_tolerance = 1.5 ; % number of edges per vertex 
% edge_number_tolerance = 4 ; % number of edges per vertex SAM before 10/29/21
% edge_number_tolerance = 2 ; % number of edges per vertex SAM 10/29/21 % SAM 12/18/21
edge_number_tolerance = 2 ; % number of edges per vertex SAM 10/29/21 % SAM 12/18/21
% edge_number_tolerance = 3 ; % number of edges per vertex  SAM 10/29/21 17:55 % SAM 12/24/21
% second_max_edges_per_vertex = max_edges_per_vertex ;
% radius_tolerance =  0.0801 ; % fractional change of radius       from origin vertex radius
%    radius_tolerance =  0.2  ; % fractional change of radius       from origin vertex radius % SAM 4/12/22
   radius_tolerance =  0.5  ; % fractional change of radius       from origin vertex radius % SAM 10/27/21
%     radius_tolerance = 1 ;  % SAM 12/24/21
%     radius_tolerance = 1.5 ;  % SAM 1/11/22
%     radius_tolerance = 2/3 ; % 11/7/21 % SAM 1/21/22
%     radius_tolerance = inf ; % 12/9/21, this did not help % SAM 12/18/21 increased total length and bifurcations marginally. Worth the move to have one less parameter
%    radius_tolerance =  3  ; % fractional change of radius       from origin vertex radius % SAM 10/27/21 %%%%% !!!!!!!!!!!!!!!!!! make a local radius tolerance
%    energy_tolerance =  0.95 ; % fractional change (positive only) from origin vertex energy
%    energy_tolerance =  0.999 ; % fractional change (positive only) from origin vertex energy
%    energy_tolerance =  1 - 1e-3 ; % fractional change (positive only) from origin vertex energy
%    energy_tolerance =  1 - 1e-6 ; % SAM 11/7/21
%    energy_tolerance =  1 - 1e-5 ; % SAM 1/24/22
%    energy_tolerance =  1 - 1e-3 ; % SAM 1/21/22
%    energy_tolerance =  1 - 1e-4 ; % SAM 11/11/21
   energy_tolerance =  1 ; % SAM 12/8/21 % SAM 1/14/21 this is working great and not taking forever to run
direction_tolerance = 1 ; % 90 degree turn (w.r.t. current trajectory) prohibited 

%% dependent inputs:
size_ratio_per_index = lumen_radius_in_microns_range( 2 ) ...
                     / lumen_radius_in_microns_range( 1 );
                 
size_tolerance = log( 1 + radius_tolerance ) / log( size_ratio_per_index );

% method = 'four' ;

% number_of_images_output = 5 ; % images being written to the h5 files for visualization purposes.
% number_of_images_output = 3 ;
number_of_images_output = 4 ;

% make chunk directory for the edges segmentation visual data
chunk_directory = [ data_directory, edges_handle, '_chunks', filesep ];

mkdir( chunk_directory );

%% get variable ranges

% %%%% removed SAM 10/1921
% 
% [ ~, ~, ~, ~, distance_tolerance_range, ~, ~, ~, path_to_energy_data, strel_distance_LUT ] ...
%   = generate_reference_image(  lumen_radius_in_microns_range, microns_per_voxel, ...
%                                vertex_space_subscripts, ...
%                                strel_apothem, distance_tolerance_per_origin_radius, ...
%                                                   data_directory, energy_handle );
% 
% 
%                                               
% %%%% END: removed SAM 10/1921


%%%% added SAM 10/1921

path_to_energy_data = [ data_directory, energy_handle ];

% distance_tolerance_range = distance_tolerance_per_origin_radius * lumen_radius_in_microns_range ;

% [ strel_linear_LUT_range, numel_of_strel_range, cum_prod_image_dims, local_subscripts_range, strel_distance_LUT_range ] ...
%     = calculate_linear_strel_range( size_of_image, microns_per_voxel, lumen_radius_in_microns_range );

%%%% END: added SAM 10/1921

vertex_space_subscripts = double( vertex_space_subscripts );
                                              
                                              
energy_file_info = h5info( path_to_energy_data );

size_of_image = energy_file_info.Datasets.Dataspace.Size ;

size_of_image = size_of_image( 1 : 3 );

strel_size_in_pixels = lumen_radius_in_microns_range( 1 ) ./ microns_per_voxel ;

[ chunk_lattice_dimensions, number_of_chunks ]                                                      ...
              = get_chunking_lattice_V190( strel_size_in_pixels, max_voxels_per_node, size_of_image );   
     
% chunk_index_range = randperm( number_of_chunks ); % not allowed by MATLAB's PARFOR (workaround?)
chunk_index_range = 1 : number_of_chunks ;     
     
% output variable declarations / initializations:
 edges2vertices(        1 : number_of_chunks, 1 ) = {[ ]};
 edge_numels(           1 : number_of_chunks, 1 ) = {[ ]};

 edge_space_subscripts( 1 : number_of_chunks, 1 ) = {{ }};
 edge_scale_subscripts( 1 : number_of_chunks, 1 ) = {{ }};
 edge_energies(         1 : number_of_chunks, 1 ) = {{ }};
%  edge_locations(        1 : number_of_chunks, 1 ) = {{ }};

% % load energy data chunk-wise in physical space and for all scales
% [ y_reading_starts, x_reading_starts, z_reading_starts,                                             ...
%   y_reading_counts, x_reading_counts, z_reading_counts  ]                                           ...
%       = get_starts_and_counts_V200( chunk_lattice_dimensions, space_strel_apothem * ones( 1, 3, 'int16' ), size_of_image, [ 1, 1, 1 ]);

[ ~, ~, ~,      ...
  ~, ~, ~,      ...
  y_writing_starts, x_writing_starts, z_writing_starts,      ...           
  y_writing_counts, x_writing_counts, z_writing_counts, ]    ...
          = get_starts_and_counts_V200( chunk_lattice_dimensions, [ 0, 0, 0], size_of_image, [ 1, 1, 1 ]);
                    
% inclusive limits
y_writing_limits = [ y_writing_starts; y_writing_starts + y_writing_counts - 1 ];
x_writing_limits = [ x_writing_starts; x_writing_starts + x_writing_counts - 1 ];
z_writing_limits = [ z_writing_starts; z_writing_starts + z_writing_counts - 1 ];
                
% numel_strel = ( 2 * space_strel_apothem + 1 ) ^ 3 ;
% 
% quantile_of_second_best = 3 / 2 / numel_strel + eps ; % 1 / 2 / numel_strel is first best
% quantile_of_third_best = 5 / 2 / numel_strel + eps ;

vertex_radii_in_microns = exp( interp1( log( lumen_radius_in_microns_range ), vertex_scale_subscripts )); %in microns 

% % distance tolerance ( 3 standard deviation of the Gaussian weighted length tolerance) of the largest
% % radius vertex in the writing volume of interest. Factor of 2 for vertex inclusion considerations:    
% vertex_distance_tolerance_in_pixels = 3 * 2 * uint16(   distance_tolerance_per_origin_radius         ...
%                                                       * vertex_radii_in_microns ./ microns_per_voxel );
vertex_distance_tolerance_in_pixels = uint16( 2 * distance_tolerance_per_origin_radius         ...
                                                * vertex_radii_in_microns ./ microns_per_voxel );

%% main PARFOR  
for chunk_index = chunk_index_range
        
    [ y_chunk_index, x_chunk_index, z_chunk_index ] = ind2sub( chunk_lattice_dimensions, chunk_index );
            
    is_vertex_in_writing_chunk                                                                  ...
                       = vertex_space_subscripts( :, 1 ) >= y_writing_limits( 1, y_chunk_index )...
                       & vertex_space_subscripts( :, 2 ) >= x_writing_limits( 1, x_chunk_index )...
                       & vertex_space_subscripts( :, 3 ) >= z_writing_limits( 1, z_chunk_index )...
                       & vertex_space_subscripts( :, 1 ) <= y_writing_limits( 2, y_chunk_index )...
                       & vertex_space_subscripts( :, 2 ) <= x_writing_limits( 2, x_chunk_index )...
                       & vertex_space_subscripts( :, 3 ) <= z_writing_limits( 2, z_chunk_index );
                       
%     if sum( is_vertex_in_writing_chunk ) == 0
%         
%         chunk_overlap_in_pixels = [ 0, 0, 0 ]; 
%     
%     else
        
        % factor of 2 needs to be divided out of distance tolerances
%         chunk_overlap_in_pixels = max( vertex_distance_tolerance_in_pixels( is_vertex_in_writing_chunk, : ) / 2 );
%         chunk_overlap_in_pixels = max( vertex_distance_tolerance_in_pixels( is_vertex_in_writing_chunk, : ));
%         chunk_overlap_in_pixels = max( vertex_distance_tolerance_in_pixels / 2 ); % !!!!! this needs to be symmetric pairwise between chunks, so make it the same for all. SAM 11/2/21
        chunk_overlap_in_pixels = 3 * max( vertex_distance_tolerance_in_pixels ); 
        
%     end
    
    [ y_reading_starts, x_reading_starts, z_reading_starts,      ...
      y_reading_counts, x_reading_counts, z_reading_counts,      ...
                                               ~, ~, ~, ~, ~, ~, ...
      y_offsets,        x_offsets,        z_offsets,        ]    ...
              = get_starts_and_counts_V200( chunk_lattice_dimensions, chunk_overlap_in_pixels, size_of_image, [ 1, 1, 1 ]);
          
    % exclusive limits for vertices because image boundary is not valid
    y_reading_limits = [ y_reading_starts; y_reading_starts + y_reading_counts - 1 ];
    x_reading_limits = [ x_reading_starts; x_reading_starts + x_reading_counts - 1 ];
    z_reading_limits = [ z_reading_starts; z_reading_starts + z_reading_counts - 1 ];
    
    is_vertex_in_reading_chunk                                                                  ...
                       = vertex_space_subscripts( :, 1 ) >  y_reading_limits( 1, y_chunk_index )...
                       & vertex_space_subscripts( :, 2 ) >  x_reading_limits( 1, x_chunk_index )...
                       & vertex_space_subscripts( :, 3 ) >  z_reading_limits( 1, z_chunk_index )...
                       & vertex_space_subscripts( :, 1 ) <  y_reading_limits( 2, y_chunk_index )...
                       & vertex_space_subscripts( :, 2 ) <  x_reading_limits( 2, x_chunk_index )...
                       & vertex_space_subscripts( :, 3 ) <  z_reading_limits( 2, z_chunk_index );

%     is_vertex_in_writing_chunk                                                                                                                ...
%                        = vertex_space_subscripts( :, 1 ) >  y_writing_limits( 1, y_chunk_index ) - vertex_distance_tolerance_in_pixels( :, 1 )...
%                        & vertex_space_subscripts( :, 2 ) >  x_writing_limits( 1, x_chunk_index ) - vertex_distance_tolerance_in_pixels( :, 2 )...
%                        & vertex_space_subscripts( :, 3 ) >  z_writing_limits( 1, z_chunk_index ) - vertex_distance_tolerance_in_pixels( :, 3 )...
%                        & vertex_space_subscripts( :, 1 ) <  y_writing_limits( 2, y_chunk_index ) + vertex_distance_tolerance_in_pixels( :, 1 )...
%                        & vertex_space_subscripts( :, 2 ) <  x_writing_limits( 2, x_chunk_index ) + vertex_distance_tolerance_in_pixels( :, 2 )...
%                        & vertex_space_subscripts( :, 3 ) <  z_writing_limits( 2, z_chunk_index ) + vertex_distance_tolerance_in_pixels( :, 3 );
                   
% 	vertices_in_chunk = find( is_vertex_in_reading_chunk & is_vertex_in_writing_chunk ); %
% 	!!!!!!!!!!!!!!!!!! this led to weird bug with long turtous strands
	vertices_in_reading_chunk = find( is_vertex_in_reading_chunk  );
    
    edges_chunk_name = [ int2str( chunk_index ), ' of ', int2str( number_of_chunks )];

    edges_chunk_file = [ chunk_directory, edges_chunk_name ];
    
    size_of_writing_chunk = [ y_writing_counts( y_chunk_index ), ...
                              x_writing_counts( x_chunk_index ), ...
                              z_writing_counts( z_chunk_index ), ...
                              number_of_images_output            ];    
    
    if isempty( vertices_in_reading_chunk )
        
        h5create( edges_chunk_file, '/d', size_of_writing_chunk )
        
    else % IF vertices
    
        size_of_reading_chunk = [ y_reading_counts( y_chunk_index ), ...
                                  x_reading_counts( x_chunk_index ), ...
                                  z_reading_counts( z_chunk_index )  ];

        vertex_space_subscripts_at_reading_chunk ...
                  = [  vertex_space_subscripts( vertices_in_reading_chunk, 1 ) - y_reading_starts( y_chunk_index ) + 1 ,  ...
                       vertex_space_subscripts( vertices_in_reading_chunk, 2 ) - x_reading_starts( x_chunk_index ) + 1 ,  ...
                       vertex_space_subscripts( vertices_in_reading_chunk, 3 ) - z_reading_starts( z_chunk_index ) + 1    ];

        vertex_locations_at_reading_chunk ...
            = sub2ind( size_of_reading_chunk, ...
                       vertex_space_subscripts_at_reading_chunk( :, 1 ), ...
                       vertex_space_subscripts_at_reading_chunk( :, 2 ), ...
                       vertex_space_subscripts_at_reading_chunk( :, 3 )  );

        size_of_index_and_energy_chunk = [ size_of_reading_chunk, 2 ];

        energy_and_index_data                                          ...
                     = h52mat(   path_to_energy_data,                  ...
                               [ y_reading_starts( y_chunk_index ),    ...
                                 x_reading_starts( x_chunk_index ),    ...
                                 z_reading_starts( z_chunk_index ),    ...
                                                 1                  ], ...
                                 size_of_index_and_energy_chunk        );

        %% main function: loop through the vertices, searching for nearby vertices to connect to
%         [ strel_linear_LUT, ~, cum_prod_image_dims ] = calculate_linear_strel( size_of_reading_chunk, step_length_per_radius );

        [ strel_linear_LUT_range, ~, cum_prod_image_dims, local_subscripts_range, strel_r_over_R_LUT_range, strel_unit_vectors_LUT_range ] ...
            = calculate_linear_strel_range( size_of_reading_chunk, microns_per_voxel, lumen_radius_in_microns_range * step_length_per_radius );

%         local_distance_tolerance_range = lumen_radius_in_microns_range * step_length_per_radius ;
%               distance_tolerance_range = lumen_radius_in_microns_range * distance_tolerance_per_origin_radius ;
        local_distance_tolerance_range = lumen_radius_in_microns_range ;
%               distance_tolerance_range = lumen_radius_in_microns_range ;
%               distance_tolerance_range = lumen_radius_in_microns_range * inf ;
        distance_tolerance = distance_tolerance_per_origin_radius ;

        strel_r_over_R_LUT_range = cellfun( @( x, y ) x * y, strel_r_over_R_LUT_range, num2cell( 1 ./ lumen_radius_in_microns_range ), 'UniformOutput', false );

        [ edge_locations               , ...
          edge_energies{  chunk_index }, ...
          edges2vertices{ chunk_index }, ...
                  energy_map, ... !!! should not be input to next function.
            vertex_index_map, ...
                 pointer_map, ...
                distance_map, ...
            branch_order_map  ]                                                   ...
                                  = get_edges_by_watershed(   energy_tolerance, edge_number_tolerance,      ...
                                                            distance_tolerance, local_distance_tolerance_range, size_tolerance,      ...
                                                            strel_linear_LUT_range, local_subscripts_range, strel_r_over_R_LUT_range, strel_unit_vectors_LUT_range, ...
                                                            vertex_locations_at_reading_chunk, ...
                                                            energy_and_index_data( :, :, :, 2 ),            ...
                                                            energy_and_index_data( :, :, :, 1 ), microns_per_voxel );
                                                        
        % remove edges that don't have at least one vertex in the writing chunk  % SAM 12/7/21
        is_edge_in_writing_chunk = is_vertex_in_writing_chunk( vertices_in_reading_chunk( edges2vertices{ chunk_index }( :, 1 ))) ...
                                 | is_vertex_in_writing_chunk( vertices_in_reading_chunk( edges2vertices{ chunk_index }( :, 2 )));
                             
        edge_locations                = edge_locations(                is_edge_in_writing_chunk    );
        edge_energies{  chunk_index } = edge_energies{  chunk_index }( is_edge_in_writing_chunk    );
        edges2vertices{ chunk_index } = edges2vertices{ chunk_index }( is_edge_in_writing_chunk, : );
        % END: SAM 12/7/21
        
% !!!!!!!! this is doen inside the get_edges_by_watershed() function
%         % erase edges that connect to the boundary
%         boundary_idx = numel( vertices_in_chunk ) + 1 ;
% 
%         included_edges_logical = edges2vertices{ chunk_index }( :, 2 ) ~= boundary_idx ;
% 
%         edge_locations                = edge_locations(                included_edges_logical    );
%     %     edge_energies{  chunk_index } = edge_energies{  chunk_index }( included_edges_logical, : );
%         edges2vertices{ chunk_index } = edges2vertices{ chunk_index }( included_edges_logical, : );

        [ edge_numels{           chunk_index }, ... !!!! waste of memory and overhead to record numels
          edge_space_subscripts{ chunk_index }, ...
          edge_scale_subscripts{ chunk_index }, ...
          edge_energies{         chunk_index } ]... % not overwriting the (adjusted) energy values output from the get_edges_by_watershed fxn
               = get_edge_vectors_V300( edge_locations, energy_map, energy_and_index_data( :, :, :, 1 ), cum_prod_image_dims );

%         if numel( find( cellfun( @( x ) numel( unique( x, 'rows' )) ~= numel( x ), edge_space_subscripts{ chunk_index } )))
% 
%             error
% 
%         end

% % there should be no orphan edges from the new get edges function % SAM 12/7/21
%         chosen_edge_indices = clean_edges_orphans( edge_space_subscripts{ chunk_index }, size_of_reading_chunk, vertex_space_subscripts_at_reading_chunk );
%         
%         if numel( chosen_edge_indices ) ~= numel( edge_scale_subscripts{ chunk_index })
%             warning('orphan edges erased')
%         end
% 
%         edges2vertices{ chunk_index }        = edges2vertices{ chunk_index }( chosen_edge_indices, : );
% 
%         edge_numels{           chunk_index } = edge_numels{           chunk_index }( chosen_edge_indices );
%         edge_space_subscripts{ chunk_index } = edge_space_subscripts{ chunk_index }( chosen_edge_indices );
%         edge_scale_subscripts{ chunk_index } = edge_scale_subscripts{ chunk_index }( chosen_edge_indices );
%         edge_energies{         chunk_index } = edge_energies{         chunk_index }( chosen_edge_indices );

        if any( is_edge_in_writing_chunk )

            %% map chunk-specific values/indices back to global values
            edges2vertices{ chunk_index } = vertices_in_reading_chunk( edges2vertices{ chunk_index });
            
            if size( edges2vertices{ chunk_index }, 2 ) ~= 2 

                edges2vertices{ chunk_index } = reshape( edges2vertices{ chunk_index }, [ 1, 2 ]);

            end
            
            edge_space_subscripts{ chunk_index } = cell2mat( edge_space_subscripts{ chunk_index });

            reading_starts = uint16([ y_reading_starts( y_chunk_index ), ...
                                      x_reading_starts( x_chunk_index ), ...
                                      z_reading_starts( z_chunk_index )  ]);

            edge_space_subscripts{ chunk_index } = edge_space_subscripts{ chunk_index } ...
                                                 + reading_starts - 1                   ;

            edge_space_subscripts{ chunk_index } = mat2cell( edge_space_subscripts{ chunk_index }, edge_numels{ chunk_index }, 3 );

            vertex_index_map( vertex_index_map - 1 > 0 ) = 1 + numel( vertex_scale_subscripts ) - vertices_in_reading_chunk( vertex_index_map( vertex_index_map - 1 > 0 ) - 1 );

        end
        
        %% writing h5 visual file

                             offsets = [ y_offsets( y_chunk_index ), ...
                                         x_offsets( x_chunk_index ), ...
                                         z_offsets( z_chunk_index )  ];    

        edges_chunk = cat( 4, double(       energy_map( offsets( 1 ) + ( 1 : size_of_writing_chunk( 1 )),    ...
                                                        offsets( 2 ) + ( 1 : size_of_writing_chunk( 2 )),    ...
                                                        offsets( 3 ) + ( 1 : size_of_writing_chunk( 3 )) )), ...
                              double( vertex_index_map( offsets( 1 ) + ( 1 : size_of_writing_chunk( 1 )),    ...
                                                        offsets( 2 ) + ( 1 : size_of_writing_chunk( 2 )),    ...
                                                        offsets( 3 ) + ( 1 : size_of_writing_chunk( 3 )) )), ...
                              double(      pointer_map( offsets( 1 ) + ( 1 : size_of_writing_chunk( 1 )),    ...
                                                        offsets( 2 ) + ( 1 : size_of_writing_chunk( 2 )),    ...
                                                        offsets( 3 ) + ( 1 : size_of_writing_chunk( 3 )) )),...); ..., ...
                              double( branch_order_map( offsets( 1 ) + ( 1 : size_of_writing_chunk( 1 )),    ...
                                                        offsets( 2 ) + ( 1 : size_of_writing_chunk( 2 )),    ...
                                                        offsets( 3 ) + ( 1 : size_of_writing_chunk( 3 )) ))  );
     ...                         double(     distance_map( offsets( 1 ) + ( 1 : size_of_writing_chunk( 1 )),    ...
     ...                                                   offsets( 2 ) + ( 1 : size_of_writing_chunk( 2 )),    ...
     ...                                                   offsets( 3 ) + ( 1 : size_of_writing_chunk( 3 )) )), ...
     
        h5create( edges_chunk_file, '/d', size_of_writing_chunk )

        mat2h5(   edges_chunk_file,                 edges_chunk );    

    end % IF exist vertices in writing chunk

end % PARFOR image chunk

%% Organizing outputs:

is_cell_not_empty = ~ cellfun( @isempty, edge_space_subscripts );

edges2vertices         = cell2mat( edges2vertices( is_cell_not_empty ));
edge_numels            = cell2mat( edge_numels(    is_cell_not_empty ));

edge_space_subscripts = edge_space_subscripts( is_cell_not_empty );
edge_scale_subscripts = edge_scale_subscripts( is_cell_not_empty );
edge_energies         = edge_energies(         is_cell_not_empty );

edge_space_subscripts  = mat2cell( cell2mat( cellfun( @cell2mat, edge_space_subscripts, 'UniformOutput', false )), edge_numels, 3 );
edge_scale_subscripts  = mat2cell( cell2mat( cellfun( @cell2mat, edge_scale_subscripts, 'UniformOutput', false )), edge_numels, 1 );
edge_energies          = mat2cell( cell2mat( cellfun( @cell2mat, edge_energies        , 'UniformOutput', false )), edge_numels, 1 );

sorted_edge_indices = sort_edges( edge_energies );

edge_numels           = edge_numels(           sorted_edge_indices );
edge_space_subscripts = edge_space_subscripts( sorted_edge_indices );
edge_scale_subscripts = edge_scale_subscripts( sorted_edge_indices );
edge_energies         = edge_energies(         sorted_edge_indices );

edges2vertices = edges2vertices( sorted_edge_indices, : );

% is_visualizing = false ;
is_visualizing = true ;

if is_visualizing

    %% make tiled master files
    % combine the chunk files into a master file outside of a parfor to avoid simulataneous writing
    % issues
    size_of_edges_image = [ size_of_image, number_of_images_output ];

    edges_file = [ data_directory, edges_handle ];

    h5create( edges_file, '/d', size_of_edges_image )

    for chunk_index = chunk_index_range

        [ y_chunk_index, ...
          x_chunk_index, ...
          z_chunk_index ] = ind2sub( chunk_lattice_dimensions, chunk_index );

        edges_chunk_name = [ int2str( chunk_index ), ' of ', int2str( number_of_chunks )];

        writing_starts = [ y_writing_starts( y_chunk_index ), ...
                           x_writing_starts( x_chunk_index ), ...
                           z_writing_starts( z_chunk_index ), ...
                           1                                  ];

        writing_counts = [ y_writing_counts( y_chunk_index ), ...
                           x_writing_counts( x_chunk_index ), ...
                           z_writing_counts( z_chunk_index ), ...
                           number_of_images_output            ];                        

        edges_chunk_file = [ chunk_directory, edges_chunk_name ];

        mat2h5( edges_file, h52mat( edges_chunk_file ), writing_starts, writing_counts );

    end % chunk FOR
end

%% Remove Chunk File directories that have been tiled and will not be referenced again

try rmdir( chunk_directory, 's' ), catch, end


%% IF visualizing: output final tifs from the edge exploration maps
if is_visualizing

    %% visualize intermediate output images from get_edges
    edges_image = h52mat( edges_file );

    visual_edges_data_handle = [ visual_data_directory, edges_handle ];

    mat2tif(         edges_image( :, :, :, 1 ) , [ visual_edges_data_handle,       '_energy_map.tif' ])
    mat2tif( uint16( edges_image( :, :, :, 2 )), [ visual_edges_data_handle, '_vertex_index_map.tif' ])
    mat2tif( uint16( edges_image( :, :, :, 3 )), [ visual_edges_data_handle,      '_pointer_map.tif' ])
    mat2tif( uint16( edges_image( :, :, :, 4 )), [ visual_edges_data_handle, '_branch_order_map.tif' ])
%     mat2tif(         edges_image( :, :, :, 5 ),  [ visual_edges_data_handle,     '_distance_map.tif' ])
    
    try delete( edges_file ); end % SAM 12/7/21
    
end % IF visualizing

end % FUNCTION

