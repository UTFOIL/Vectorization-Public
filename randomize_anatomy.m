function [ strand_space_subscripts, strand_scale_subscripts ] ...
                = randomize_anatomy( strand_space_subscripts, strand_scale_subscripts, microns_per_voxel, lumen_radius_in_microns_range, size_of_image, network_statistics, radii_threshold )
%% randomize_anatomy
% SAM 7/2/19
% 
% randomize the anatomy of the strand objects by uniform randomly translating each object below a
% certain radii_threshold to a new location in the image.

strand_subscripts = cellfun(      @(      space_subscripts,      scale_subscripts  )  ...
                                    [      space_subscripts,      scale_subscripts ], ...
                                   strand_space_subscripts, strand_scale_subscripts,    ...
                                         'UniformOutput', false                     );

degrees_of_strands = cellfun( @( x ) size( x, 1 ), strand_subscripts );

%% construct the relative sphere elements at all the scales
pixels_per_sigma_range = lumen_radius_in_microns_range ./ microns_per_voxel ;

number_of_scales = size( pixels_per_sigma_range, 1 );

scale_subscript_range = 1 : number_of_scales ;

% ---------------------------------------- strand objects -------------------------------------------

strand_element_subscripts      = cell( number_of_scales, 1 );
strand_element_linear_indexing = cell( number_of_scales, 1 );

for scale_subscript = scale_subscript_range
    
    [ strand_element_linear_indexing{ scale_subscript },                                 ...
      strand_element_subscripts{      scale_subscript }, ]                               ...
        = construct_structuring_element_V190( pixels_per_sigma_range( scale_subscript, : ), size_of_image  );
 
end % constructing relative elements FOR scale

%% measure sizes of strands for repositioning constraints

% intialize the painted_image with a blank canvas
painted_image = zeros( size_of_image, 'logical' );

strand_max_y = cellfun( @( x ) max( x( :, 1 ) + pixels_per_sigma_range( round( x( :, 4 )), 1 )), strand_subscripts );
strand_max_x = cellfun( @( x ) max( x( :, 2 ) + pixels_per_sigma_range( round( x( :, 4 )), 2 )), strand_subscripts );
strand_max_z = cellfun( @( x ) max( x( :, 3 ) + pixels_per_sigma_range( round( x( :, 4 )), 3 )), strand_subscripts );

strand_min_y = cellfun( @( x ) min( x( :, 1 ) - pixels_per_sigma_range( round( x( :, 4 )), 1 )), strand_subscripts );
strand_min_x = cellfun( @( x ) min( x( :, 2 ) - pixels_per_sigma_range( round( x( :, 4 )), 2 )), strand_subscripts );
strand_min_z = cellfun( @( x ) min( x( :, 3 ) - pixels_per_sigma_range( round( x( :, 4 )), 3 )), strand_subscripts );

strand_maxs = [ strand_max_y, strand_max_x, strand_max_z ];
strand_mins = [ strand_min_y, strand_min_x, strand_min_z ];

strand_lengths = ceil( strand_maxs - strand_mins );

% Fix vessels over a certain size threshold to not be randomize
if exist( 'radii_threshold', 'var' )

    is_strand_fixed = network_statistics.strand_ave_radii > radii_threshold ;
    
else
    
    is_strand_fixed = zeros( size( strand_space_subscripts ), 'logical' );
    
end

z_translation_factor = 0 ; % 0 or 1.  1 allows random translation to any image location in z, 0 is no translation in z

%% loop through the strands from largest to smallest to paint

[ ~, strand_indices_sorted_by_size ] = sort( network_statistics.strand_volumes, 'descend' );

strand_indices_sorted_by_size = strand_indices_sorted_by_size';

for strand_index = strand_indices_sorted_by_size
  
    placing_strand = true ;

    while placing_strand

        if is_strand_fixed( strand_index ) % paint without looking nor randomizing
            
            strand_space_subscripts_temp = strand_subscripts{ strand_index }( :, 1 : 3 );
            
            strand_space_subscripts_int = round( strand_space_subscripts_temp );

            strand_position_index_range = uint16( 1 : degrees_of_strands( strand_index ));

        else % ELSE attempt to randomize

            % choose a new random uniform location so that this object fits entirely in the image            
            strand_start_subscripts = 1 + ( size_of_image - 1 - strand_lengths( strand_index, : )) .* rand( 1, 3 );

            strand_space_subscripts_temp = strand_subscripts{ strand_index }( :, 1 : 3 ) + ( strand_start_subscripts - strand_mins( strand_index, : )) .* [ 1, 1, z_translation_factor ];

            strand_space_subscripts_int = round( strand_space_subscripts_temp );

            strand_position_index_range = uint16( randperm( degrees_of_strands( strand_index )));

        end % ELSE attempt to randomize


        strand_position_linear_indexing_cell = cell( degrees_of_strands( strand_index ), 1 );

        % assume the strand is good unless or until shown otherwise
        is_strand_placed = true ;

        % loop through the positions along each strand
        for strand_position_index = strand_position_index_range

            voxel_index = sub2ind( size_of_image, strand_space_subscripts_int( strand_position_index, 1 ), ...
                                                  strand_space_subscripts_int( strand_position_index, 2 ), ...
                                                  strand_space_subscripts_int( strand_position_index, 3 )  );

            strand_position_linear_indexing                                                                 ...
                = voxel_index                                                                             ...
                + strand_element_linear_indexing{ round( strand_subscripts{ strand_index }( strand_position_index,   4   ))};

            strand_position_linear_indexing_cell{ strand_position_index } = strand_position_linear_indexing ;

            if ~ is_strand_fixed( strand_index ) 

                % check if this area has already been painted (should only have zeros and potentially the start
                % and end vertices on the canvas here, otherwise this real estate is taken).
                conflicting_object_indices = painted_image( strand_position_linear_indexing );

                % if there are any nonzero values in the conflicting objects set then we don't chose this
                % strand, we don't need to continue the FOR loop
                if any( conflicting_object_indices )

                    is_strand_placed = false ;

                    break

                end % IF painting conflict
            end % IF randomized strand
        end % strand position FOR

        % If the strand is a keeper, paint the image in its volume of influence.
        if is_strand_placed

            placing_strand = false ;

            strand_entire_linear_indexing = cell2mat( strand_position_linear_indexing_cell );

    %         painted_image( strand_entire_linear_indexing ) = index_code_for_strand ;
            painted_image( strand_entire_linear_indexing ) = true ;

            strand_subscripts{ strand_index }( :, 1 : 3 ) = strand_space_subscripts_temp ;

        end % IF strand is placed
    end % WHILE placing strand
end % sorted strand FOR

strand_space_subscripts = cellfun( @( x ) x( :, 1 : 3 ), strand_subscripts, 'UniformOutput', false );
strand_scale_subscripts = cellfun( @( x ) x( :,   4   ), strand_subscripts, 'UniformOutput', false );

end % FUNCTION randomize_anatmoy