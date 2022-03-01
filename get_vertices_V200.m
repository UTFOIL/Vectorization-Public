function [ space_subscripts, scale_subscripts, energy_values, positions ]                           ...
                             = get_vertices_V200( lumen_radius_in_microns_range, microns_per_pixel, ...
                                                          space_strel_apothem, max_voxels_per_node, ...
                                                                   energy_upper_bound, energy_file  )
                            
%% find all the candidate vertices: find minima in vesselness images
% SAM 1/8/18
%
% V02 in which the inputs are changed to be ammenable to the vectorize_V02 script 3/5/18 SAM
%
% V10 in which the octave structure is eliminated SAM 3/29/18
%
% V11 in which the maximized value is not the negative laplacian but the vesselness feature SAM
% 4/16/18
%
% V130 in which the number_of_scales bug was fixed to have one more scale than in V11. Also the
% vertex position is assigned to the integer scale and not the half integer scale as was the case
% with the difference of Gaussian approximation to the Laplacian). Also, The word "vesselness" is
% changed to "energy" in the function names.
%
% Also eliminated the need to pass the interpolate handle.  Size of image is now retrieved from the
% energy file. SAM 4/25/18
% 
% note: either remove subpixel localization  or extend it to finding the edges subpixelly too. SAM
% 4/25/18
%
% V131 in which subpixel localization is dropped SAM 5/5/18
%
% V150 in which the energy file is min projected and the second index of the 4th dimension is the
% image of scale indices associated with each minimum.  % 5/7/18 SAM
%
% V190 in which the location of the energy file to be loaded has changed.  SAM 8/26/18
%
% V191 location of energy file is changed back to where it was at V150. V190 came with many more
% updates than those stated.  SAM 11/2/18
%
% V200 in which the input energy file comes as a single string variable instead of a directory and a
% handle.

% number_of_scales = size( pixels_per_sigma_range, 1 );

scales_per_octave = log( 2 )                                              ...
                             / log(   lumen_radius_in_microns_range( 2 )  ...
                                    / lumen_radius_in_microns_range( 1 )) ...
                             / 3 ; % divide by three for the volume interperetation of octave
                       

space_strel_range = - space_strel_apothem : space_strel_apothem ;

energy_file_info = h5info( energy_file );

size_of_image = energy_file_info.Datasets.Dataspace.Size ;

size_of_image = size_of_image( 1 : 3 );

strel_size_in_pixels = lumen_radius_in_microns_range( 1 ) ./ microns_per_pixel ;

[ chunk_lattice_dimensions, number_of_chunks ]                                                      ...
              = get_chunking_lattice_V190( strel_size_in_pixels, max_voxels_per_node, size_of_image );   
     
chunk_index_range = 1 : number_of_chunks ;     
         
% output variable declarations / initializations:
 space_subscripts( 1 : number_of_chunks, 1 ) = { uint16([ ])};
 scale_subscripts( 1 : number_of_chunks, 1 ) = {        [ ] };
    energy_values( 1 : number_of_chunks, 1 ) = {        [ ] };
        positions( 1 : number_of_chunks, 1 ) = {        [ ] };

% load energy data chunk-wise in physical space and for all scales
[ y_reading_starts, x_reading_starts, z_reading_starts,                                             ...
  y_reading_counts, x_reading_counts, z_reading_counts  ]                                           ...
      = get_starts_and_counts_V200( chunk_lattice_dimensions, space_strel_apothem * ones( 1, 3, 'int16' ), size_of_image, [ 1, 1, 1 ]);
  
% numel_strel = ( 2 * space_strel_apothem + 1 ) ^ 3 ;
% 
% quantile_of_second_best = 3 / 2 / numel_strel + eps ; % 1 / 2 / numel_strel is first best
% quantile_of_third_best = 5 / 2 / numel_strel + eps ;

%% main PARFOR  
parfor chunk_index = chunk_index_range

    [ y_chunk_index, x_chunk_index, z_chunk_index ] = ind2sub( chunk_lattice_dimensions, chunk_index );

    number_of_vertices_found = 0 ;
    
    size_of_energy_chunk = [ y_reading_counts( y_chunk_index ), ...
                             x_reading_counts( x_chunk_index ), ...
                             z_reading_counts( z_chunk_index )  ];
                         
    size_of_index_and_energy_chunk = [ size_of_energy_chunk, 2 ];

    index_and_energy_chunk                                         ...
                 = h52mat(   energy_file,                          ...
                           [ y_reading_starts( y_chunk_index ),    ...
                             x_reading_starts( x_chunk_index ),    ...
                             z_reading_starts( z_chunk_index ),    ...
                                             1                  ], ...
                             size_of_index_and_energy_chunk        );
                         
%     number_of_elements_in_energy_chunk = prod( size_of_energy_chunk );

    % Make two copies of the chunk, one for reference and the other for marking places that
    % have been searched. 
    energy_chunk_temp = zeros( size_of_energy_chunk );

    % Mark the border as searched before we begin so we don't find false positives there.
    energy_chunk_temp(                                              ...
            1 + space_strel_apothem : end - space_strel_apothem,    ...
            1 + space_strel_apothem : end - space_strel_apothem,    ...
            1 + space_strel_apothem : end - space_strel_apothem   ) ...
                                                                    ...
        = index_and_energy_chunk(                                   ...
            1 + space_strel_apothem : end - space_strel_apothem,    ...
            1 + space_strel_apothem : end - space_strel_apothem,    ...
            1 + space_strel_apothem : end - space_strel_apothem,    ...
                                    2                             );

    energy_chunk_temp( energy_chunk_temp > energy_upper_bound ) = 0 ;
    
    [ energy_min_value, energy_min_index ] = min( energy_chunk_temp( : ));    
    
    while energy_min_value < 0

        % get physical and scale coordinates of the global min
        [ y_subscript, x_subscript, z_subscript ] ...
            = ind2sub( size_of_energy_chunk, energy_min_index );

        % pull out the local window using the space strel
        energy_window                                                  ...
            = index_and_energy_chunk( y_subscript + space_strel_range, ...
                                      x_subscript + space_strel_range, ...
                                      z_subscript + space_strel_range, ...
                                      2                                );

        % test if this point is a min of the window (not just a min of the temporary energy
        % chunk)
        is_vertex_found                    ...
            =  min( energy_window( : )) ...
            == energy_min_value;
        
%         is_vertex_found = energy_min_value <= quantile( energy_window( : ), quantile_of_third_best );

        %% vertex found
        if is_vertex_found

            number_of_vertices_found = number_of_vertices_found + 1 ;

            space_subscripts{ chunk_index }                            ...
                = [ space_subscripts{  chunk_index };                  ...
                      [ uint16( y_subscript ),                         ...
                        uint16( x_subscript ),                         ...
                        uint16( z_subscript )  ]                       ...
                    + [ uint16( y_reading_starts( y_chunk_index )),    ...
                        uint16( x_reading_starts( x_chunk_index )),    ...
                        uint16( z_reading_starts( z_chunk_index ))  ]  ...
                    - 1                                                ];

              s_subscript = index_and_energy_chunk( energy_min_index );
                
              scale_subscripts{ chunk_index }           ...
                  = [ scale_subscripts{  chunk_index }; ...
                      s_subscript                       ];

              energy_values{ chunk_index } = [ energy_values{ chunk_index }; ...
                                               energy_min_value              ];

            % the fourth column of positions summarizes the size.  Take that factor and multiply it
            % by the first row of pixels_per_sigma_range to get the 3 radii.
            positions{ chunk_index }                                                            ...
                = [ positions{ chunk_index };                                                   ...
                    ( double( space_subscripts{ chunk_index }( number_of_vertices_found, : ))), ...
                                                2 ^ (( s_subscript - 1 ) / scales_per_octave )  ];

        end % IF vertex found

        energy_chunk_temp(                                           ...
                              y_subscript + space_strel_range,       ...
                              x_subscript + space_strel_range,       ...
                              z_subscript + space_strel_range,       ...
                                              :                ) = 0 ;


        [ energy_min_value, energy_min_index ] = min( energy_chunk_temp( : ));

    end % WHILE searching for vertices               
    
end % PARFOR chunks

%% combining outputs from the parallel threads into single vectors/matrices

space_subscripts  = cell2mat( space_subscripts  ); 
scale_subscripts  = cell2mat( scale_subscripts  ); 
energy_values     = cell2mat( energy_values     ); 
positions         = cell2mat( positions         ); 

end % FUNCTION