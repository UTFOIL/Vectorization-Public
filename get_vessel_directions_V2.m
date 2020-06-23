function [ vessel_directions, edge_subscripts_smoothed, tissue_types ]                                ...
     = get_vessel_directions_V2( edge_subscripts, edge_indices_in_strands, edge_backwards_in_strands, ...
           smoothing_kernel_sigma_to_lumen_radius_ratio, tissue_type_cutoffs, lumen_radius_in_pixels_range  )
%% get_vessel_directions SAM for Chakameh 7/12/18
% The purpose of this function is to loop through the strand/junction objects putting the edges
% together into strands or junctions as encoded by the inputs edge_indices_in_strands and
% edge_backwards_in_strands.  Then for each strand doing some spatial and radial smoothing and then
% computing the spatial derivative with respect the axial direction of the strand.  The output is
% then put back in the organization of the edge_subscripts input (a cell array of edges).  Note that
% the length units of the different components of the vessel direction vectors may not be the same
% due to anisotropic sampling in y, x, and z.
%
% V2, the sigma of the gaussian kernel for smoothing is an input % SAM 7/16/18
%
% V2 also, this version was merged with V4 which includes the tissue type cutoffs input and
% tissue_types output. SAM 3/5/19

% concatenate infinity onto the tissue type cutoffs to largest tissue assignment work later.
tissue_type_cutoffs( end + 1 ) = Inf ;

number_of_strands = length( edge_indices_in_strands );

strand_index_range = 1 : number_of_strands ;

% count the number of edges in each strand, output into vector with length of the number of strands
numbers_of_edges_in_strands                                                                         ...
            = cellfun( @( indices_at_strand ) size( indices_at_strand, 1 ), edge_indices_in_strands );

vessel_directions        = cell( size( edge_subscripts ));
edge_subscripts_smoothed = cell( size( edge_subscripts ));

tissue_types            = zeros( size( edge_subscripts ));

for strand_index = strand_index_range

    %% Assembly into Strands from Edges                                                             
    % look up the edge_subscripts associated with the edges in this strand. The variable
    % edge_subscripts_at_strand_in_edges is a cell array of vector lists, one list for each edge in
    % the current strand. The edges are called in their order of appearance in the strand, but
    % within each edge the list is as likely to be backwards of the master list as it is to be in
    % the same order. So we will have to flip about half of them so that there is a preserved order
    % in listing between edges. The edges to be flipped are encoded in the edge_backwards_in_strands
    % input.
    edge_subscripts_at_strand_in_edges = edge_subscripts( edge_indices_in_strands{ strand_index });
    
    edge_index_in_strand_range = 1 : numbers_of_edges_in_strands( strand_index );
    
    % loop through the edges in the curent strand
    for edge_index_in_strand = edge_index_in_strand_range
        
        % flip this edge if needed
        if edge_backwards_in_strands{ strand_index }( edge_index_in_strand )
        
                      edge_subscripts_at_strand_in_edges{ edge_index_in_strand } ...
            = flipud( edge_subscripts_at_strand_in_edges{ edge_index_in_strand });
            
        end % backwards edge IF
    
        % trim the last vector of each edge except on the last edge (to avoid double listing of the
        % vertices
        if edge_index_in_strand < numbers_of_edges_in_strands( strand_index )
            
              edge_subscripts_at_strand_in_edges{ edge_index_in_strand }                 ...
            = edge_subscripts_at_strand_in_edges{ edge_index_in_strand }( 1 : end - 1, : );
        
        end % not last edge IF
    end % edge in strand FOR                                   
    
    % combine all the subscripts from all the edges into a single ordered list representing the
    % current strand
    edge_subscripts_at_strand = cell2mat( edge_subscripts_at_strand_in_edges );
    
    %% Smoothing along strands                                                                          
    
    % get the average size of the strand
    average_scale_index_at_strand = round( mean( edge_subscripts_at_strand( :, 4 )));

    % take geo mean of the voxel side lengths to get a characteristic voxel length.  Then convert the
    % average scale index to a lumen radius in units of voxels.
    average_lumen_radius_in_pixels_at_strand = prod( lumen_radius_in_pixels_range( average_scale_index_at_strand, : ), 2 ) .^ ( 1 / 3 );    
    
    % Smooth the subscripts in each x, y, z, and r_subscript dimension along the vessel axis. Use a
    % Gaussian 1D kernel (st. dev. = 1 [units of path index (physical space between indices may
    % vary)]).
    %
    % NOTE to programmer: in a future version, you should include the edge energies variable as an
    % input to this function and manipulate that cell array as we did the edge_subscripts input and
    % then use the edge energy as a weighting function for the following smoothing operation.  (It
    % will still include the gaussian weighting too (perhaps applied to the energies to avoid
    % numerical instability on the application of.)  Consider converting edge energy to weighting
    % using a sliding window operation that does an exponential transform of the energies inside
    % with characteristic energy given by an input fraction of the central energy value at each
    % window.  (This would be similar to the transformation done when originally looking for edges
    % on random walkds.) SAM 7/13/18 (Pad weighting vector by zeros on either end to give no weight
    % to the extrapolated values created to make the convolution possible near the ends.
    sigma_smoothing_kernel                             ...
        = smoothing_kernel_sigma_to_lumen_radius_ratio ...
        * average_lumen_radius_in_pixels_at_strand ;
    
    radius_of_smoothing_kernel = round( 3 * sigma_smoothing_kernel ); 
    
    kernel_index_range = - radius_of_smoothing_kernel : radius_of_smoothing_kernel ; 
    
    gaussian_kernel = exp( - kernel_index_range .^ 2 / ( 2 * sigma_smoothing_kernel ^ 2 ));
    
    gaussian_kernel = gaussian_kernel / sum( gaussian_kernel );    
    
    % Pad the strand by repeating the first and last entries before convolution. Pad with twice the
    % radius of the smoothing kernel, so that at 1 radius from the start, the convolution is valid
    % and 99.7 % guaranteed to coincide with the first vertex.
    number_of_vectors_at_strand = size( edge_subscripts_at_strand, 1 );
    
    edge_subscripts_at_strand_padded = zeros( 4 * radius_of_smoothing_kernel + number_of_vectors_at_strand, ...
                                              4                                                             );
    
    edge_subscripts_at_strand_padded(     1 : 2 * radius_of_smoothing_kernel,     : ) =  edge_subscripts_at_strand(  1, :  )     ...
                                                                                      .* ones( 2 * radius_of_smoothing_kernel, 4 );
                                                                                  
    edge_subscripts_at_strand_padded(         2 * radius_of_smoothing_kernel + 1                                   ...
                                      : end - 2 * radius_of_smoothing_kernel    , : ) =  edge_subscripts_at_strand ;

    edge_subscripts_at_strand_padded(   end - 2 * radius_of_smoothing_kernel + 1                                                 ...
                                      : end                                     , : ) =  edge_subscripts_at_strand( end, : )     ...
                                                                                      .* ones( 2 * radius_of_smoothing_kernel, 4 );
    
    edge_subscripts_at_strand_smoothed = zeros([ 2 * radius_of_smoothing_kernel, 0 ] + size( edge_subscripts_at_strand ));
    
    % apply smoothing convolution
    for subscript_ordinate = 1 : 4 
    
        edge_subscripts_at_strand_smoothed( :, subscript_ordinate )                                 ...
                                 = conv( edge_subscripts_at_strand_padded( :, subscript_ordinate ), ...
                                         gaussian_kernel', 'valid'                                  );
        
    end

    % guarantee that the first location coincides with the first vertex
    edge_subscripts_at_strand_smoothed(  1 , : ) = edge_subscripts_at_strand(  1 , : );
    edge_subscripts_at_strand_smoothed( end, : ) = edge_subscripts_at_strand( end, : );
    
    %% Interpolate back to the original number of vectors as the original strand
    
    % precalculate cumulative distances covered by the vectors in each edge for later interpolation
    edge_cumulative_length                                                                                        ...
        = [ 0; cumsum( sum((   edge_subscripts_at_strand_smoothed( 1 + 1 : end    , 1 : 3 )                       ...
                             - edge_subscripts_at_strand_smoothed( 1     : end - 1, 1 : 3 )) .^ 2, 2 ) .^ 0.5 )]; ...                               
        
    edge_subscripts_at_strand_smoothed                           ...
        = interp1( edge_cumulative_length,                       ...
                   edge_subscripts_at_strand_smoothed,           ...
                   linspace( 0, edge_cumulative_length( end ),   ...
                                number_of_vectors_at_strand   )' );
                            
    %% Spatial Derivative Approximation
                            
    % convert to a common real distance unit in each spatial dimension
    edge_subscripts_at_strand_smoothed( :, 1 : 3 ) = edge_subscripts_at_strand_smoothed( :, 1 : 3 ) ...
                                                  ./ lumen_radius_in_pixels_range( 1, : );    
    
    % approximate the spatial derivative with respect to the axial direction of the vessel using a
    % symmetric difference.
    vessel_directions_at_strand_cropped                                                             ...
                                     = edge_subscripts_at_strand_smoothed( 1 + 2 : end    , 1 : 3 ) ...
                                     - edge_subscripts_at_strand_smoothed( 1     : end - 2, 1 : 3 ) ;
                            
    % Pad the ends of the cropped output with repeats of the first and last entries
    vessel_directions_at_strand = [ vessel_directions_at_strand_cropped( 1,   : ); ...
                                    vessel_directions_at_strand_cropped          ; ...
                                    vessel_directions_at_strand_cropped( end, : )  ];
    
    % normalize the spatial derivatives with respect to the pythagorean sum (make them unit vectors)
    vessel_directions_at_strand_unit =         vessel_directions_at_strand                   ...
                                     ./ ( sum( vessel_directions_at_strand .^ 2, 2 )) .^ 0.5 ;
                                 
    % convert back to units of voxels
    edge_subscripts_at_strand_smoothed( :, 1 : 3 ) = edge_subscripts_at_strand_smoothed( :, 1 : 3 ) ...
                                                  .* lumen_radius_in_pixels_range( 1, : );                                 

    %% Dissassembly from strands back into edges                                                    
    % loop back through the strands and this time use the edge lookup table not to read edge
    % positions but to write the new new informtaion (smoothed coordinates and spatial derivatives)
    % from strand form into the edge form.
    
    % count the number of vectors inside of the current strand for each of its edges
    numbers_of_vectors_in_edges = cellfun( @( subscripts_at_edge ) size( subscripts_at_edge, 1 ),   ...
                                                   edge_subscripts_at_strand_in_edges               );
                                       
    cumulative_number_of_vectors_at_strand = [ 0; cumsum( numbers_of_vectors_in_edges )];
        
    % loop through the edges in the current strand
    for edge_index_in_strand = edge_index_in_strand_range
        
        % put the last vector of each edge (except last one) back in (every edge includes both vertices)
        if edge_index_in_strand < numbers_of_edges_in_strands( strand_index )
        
            vector_indices_of_edge_in_strand                                             ...
                = cumulative_number_of_vectors_at_strand( edge_index_in_strand     ) + 1 ...
                : cumulative_number_of_vectors_at_strand( edge_index_in_strand + 1 ) + 1 ; % note the " + 1 "
                        
        else % last edge ELSE
        
            vector_indices_of_edge_in_strand                                             ...
                = cumulative_number_of_vectors_at_strand( edge_index_in_strand     ) + 1 ...
                : cumulative_number_of_vectors_at_strand( edge_index_in_strand + 1 ) 	 ; % compare to note ^ in this IF
            
        end % not last edge IF
        
        % flip this edge if needed
        if edge_backwards_in_strands{ strand_index }( edge_index_in_strand )

            vector_indices_of_edge_in_strand = fliplr( vector_indices_of_edge_in_strand );

        end % backwards edge IF

        edge_index_in_edge_list = edge_indices_in_strands{ strand_index }( edge_index_in_strand );
        
        vessel_directions{ edge_index_in_edge_list }                                                ...
                            = vessel_directions_at_strand_unit( vector_indices_of_edge_in_strand, : );
                        
        edge_subscripts_smoothed{ edge_index_in_edge_list }                                         ...
                          = edge_subscripts_at_strand_smoothed( vector_indices_of_edge_in_strand, : );
                      
        tissue_types( edge_index_in_edge_list )                                                                ...
            = find(    tissue_type_cutoffs                                                                     ...
                    >= median( edge_subscripts_at_strand_smoothed( vector_indices_of_edge_in_strand, 4 )), 1 ) ...
            + 1 ; %add 1 so that capillary will be tissue type 2 and vessel tissue type 3

    end % edge in strand FOR
end % strand FOR
end % FUNCTION