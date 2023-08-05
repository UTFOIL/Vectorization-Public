function [ chosen_vertex_original_indices, painted_image ]                                          ...
                        = choose_vertices_V200( space_subscripts, scale_subscripts, energy_values,  ...
                                                lumen_radius_in_microns_range, microns_per_pixel,   ...
                                                size_of_image, length_dilation_ratio                )
%% Choose vertices_V180
% the purpose of this function is to clean up the output of the get_vertices function. It does this
% by sorting the vertices from best to worst (highest to lowest) energy and then looping through and
% eliminating vertices whose volumes overlap with vertices of lower energy.
%
% SAM 5/29/18
%
% V181 in which the elimination is one sided: volume overlap isn't sufficient to eliminate vertices.
% We look to see if the current vertex's center (not volume) is inside the volume of any previously
% accepted vertex. SAM 5/29/18
%
% V182 is somewhere between V180 and V181 SAM 5/29/18, the inner sigma volumes need to overlap (not
% the sizes as in the case of V180.
%
% V183 in which the size of the volumes needed to overlap to exclude a vertex is input as a factor
% of sigma. SAM 5/29
%
% V184 where instead of searching for previously identified vertices that will interfere with the
% current object, the method is to look at a painted image of previously identified vertices where
% previous vertices have been painted. SAM 5/30/18
%
% V190 removing vertices that overlap the image boundary before entering the for loop. SAM 8/16/18
%
% V191 cropping and sorting by energy are done outside this function in the crop_vertices_V190
% function. eliminating references to sigma; instead referring to lumen radius.  SAM 8/23/18
%
% V200 added the length_dilation_ratio input. SAM 1/21/19

% create the relative sphere elements at all the scales
scale_subscripts = round( scale_subscripts );

unique_scale_subscripts = unique( scale_subscripts )';

largest_scale_subscript = max( unique_scale_subscripts );

structuring_element_linear_indexing_templates = cell( largest_scale_subscript, 1 );

radii_in_pixels_range = length_dilation_ratio * lumen_radius_in_microns_range ./ microns_per_pixel ;

for scale_subscript = unique_scale_subscripts
      
    % find all pixel locations within the ellipsoid radii from the vertex position    
    structuring_element_linear_indexing_templates{ scale_subscript }                                            ...
        = int32( construct_structuring_element( radii_in_pixels_range( scale_subscript, : ), size_of_image ));
     
end % constructing relative elements FOR scale

% begin the painting procedure

% enumerate the original vertices for purposes of outputting the indices of the chosen vertices
number_of_vertices = length( energy_values );

vertex_index_range = 1 : uint32( number_of_vertices );
     
 space_subscripts = int32( space_subscripts );

% % initialize the vector to assign logical 1 or 0 to whehter each vertex is chosen
% chosen_vertex_logical = zeros( number_of_vertices, 1, 'logical' );

% initialize the vector to hold the original indices of the chosen vertices
number_of_vertices = length( energy_values );

chosen_vertex_original_indices = zeros( number_of_vertices, 1, 'uint32' );

% intialize the painted_image with a blank canvas
painted_image = zeros( size_of_image, 'logical' );
  scale_image = zeros( size_of_image, 'uint8' ); % !!!!! number of scales limited to 256

num_voxels = prod( size_of_image );
  
% initialize chosen vertex counter
number_of_chosen_vertices = uint32( 0 );

% precalculate linear indices from space subscripts in the image
vertex_position_linear_indices =   space_subscripts( :, 1 )                                               ...
                               + ( space_subscripts( :, 2 ) - 1 ) * size_of_image( 1 )                    ...
                               + ( space_subscripts( :, 3 ) - 1 ) * size_of_image( 1 ) * size_of_image( 2 );

% vertex_position_linear_indices_cell = num2cell( vertex_position_linear_indices );
% 
% structuring_element_linear_indexing = structuring_element_linear_indexing_templates( scale_subscripts );
%                            
% vertex_structure_positions_linear_indexing = cellfun( @plus, structuring_element_linear_indexing,   ...
%                                                                vertex_position_linear_indices_cell, ...
%                                                                              'UniformOutput', false );
              
% loop through the vertices from best to worst (lowest to highest energy). Should be sorted before
% this function.
for vertex_index = vertex_index_range 
                
    vertex_structure_positions_linear_indexing = translate_strel_template( vertex_position_linear_indices(                              vertex_index ), ...
                                                                       structuring_element_linear_indexing_templates{ scale_subscripts( vertex_index )});
    
	% check if this area has already been painted
%     if all( ~ painted_image( vertex_structure_positions_linear_indexing{ vertex_index }))
    if all( ~ painted_image( vertex_structure_positions_linear_indexing ))
       
        % accumulate the current vertex into the list of the chosen ones
        number_of_chosen_vertices = number_of_chosen_vertices + 1 ;
        
        chosen_vertex_original_indices( number_of_chosen_vertices ) = vertex_index ;
        
%         scale_in_vertex = max( scale_image( vertex_structure_positions_linear_indexing ));
%         
%         % adopt previously defined scale if present
%         if scale_in_vertex, scale_subscripts( vertex_index ) = scale_in_vertex ; end
        
        vertex_structure_positions_linear_indexing = translate_strel_template( vertex_position_linear_indices(                              vertex_index ), ...
                                                                           structuring_element_linear_indexing_templates{ scale_subscripts( vertex_index )});
        
                    vertex_structure_positions_linear_indexing                  ...
        = min( max( vertex_structure_positions_linear_indexing, 1 ), num_voxels );
                                                                       
        % paint the image to mark this newly chosen vertex        
%         painted_image( vertex_structure_positions_linear_indexing{ vertex_index }) = true ;
        painted_image( vertex_structure_positions_linear_indexing ) = true ;
        
    else
        
        % record the scale of this vertex within its volume
                 scale_image( vertex_structure_positions_linear_indexing ) ...
          = max( scale_subscripts( vertex_index ), ...
                 scale_image( vertex_structure_positions_linear_indexing ));
        
    end % IF area is blank canvas (unpainted with any previous vertices)
                                           
end % vertex FOR

% original_vertex_indices = original_vertex_indices( chosen_vertex_logical );

chosen_vertex_original_indices = chosen_vertex_original_indices( 1 : number_of_chosen_vertices );

end % FUNCTION

function vertex_structure_positions_linear_indexing = translate_strel_template( vertex_position_linear_index, structuring_element_linear_indexing_template )

    vertex_structure_positions_linear_indexing ...
         = structuring_element_linear_indexing_template ...
         +     vertex_position_linear_index                            ;

     
end

