function visualize_strands( strand_subscripts, strand_fill_values, ellipsoid_radii_in_voxels, ...
                               size_of_image, spheres_visual_file, visual_contrast )
%% SAM 12/12/17 
% adapted from the function with the same name in the folder AA
%
% V2, in which subpixel scales and locations are passed and the ellipsoids are
% sized and placed appropriately in the image.
%
% V02 in which the inputs are changed to be ammenable to the vectorize_V02 script 3/5/18 SAM
%
% V10 in which the negative laplacian is used (instead of just the value 1) as the fill value for
% each vertex sphere SAM 4/4/18
%
% function name changed from visualize_vertices to visualize_edges
% 
% V11 which is exactly the same function as visualize_vertices_V10
%
% V12 in which the total mean_edge_energies are passed instead of the edge energies at every
% position along the trajectory.  Trajectories are sorted by their total_mean_edge_energies from
% worst to best, and the better trajectories overwrite the worse trajectories in the resulting
% image, essentially plotting the maximum of the total_mean_edge_energies at every point in the
% space that belongs to the volume of a trajectory.  Also the positions is now input as a cell array
% with entries corresponding to trajectories.  SAM 5/5/18
%
% V160 in which the max_edge_energies are passed instead of the "total mean_edge_energies" (inspired
% by thinking about the trajectories as reaction pathways and the max being like the activation
% energy).  Also combing the centerlines into this same function. % SAM 5/14/18
%
% V161 in which the position elements are constructed once into a look up table.  SAM 5/14/18
%
% V180 in which the subscript inputs are allowed to be non integer
%
% This function renders the edge objects by filling ellipsoidal objects in an image around each
% edge centerline location (whose size and shape is determined by the radius of the edge at that
% location and the resolution of the image).  SAM 7/21/19

strands_image   = zeros( size_of_image );
weighting_image = zeros( size_of_image ) + eps ;

% number_of_image_voxels = prod( size_of_image );

number_of_edges = length( strand_fill_values );

% round the edge subscripts to integers 
strand_subscripts = cellfun( @( v ) round( v ), strand_subscripts, 'UniformOutput', false );

% pre-calculating all of the ellipsoidal structuring elements to be used to paint in the image
number_of_scales = size( ellipsoid_radii_in_voxels, 1 );

strel_indices = cell( number_of_scales, 1 );

for scale_index = 1 : number_of_scales

    % find all pixel locations within the ellipsoid radii from the vertex position
    strel_indices{ scale_index }                                         ...
        = construct_structuring_element( ellipsoid_radii_in_voxels( scale_index, : ), size_of_image );
        
end % scale FOR

% loop through the edges
for edge_index = 1 : number_of_edges
    
    subscripts_at_edge = strand_subscripts{ edge_index };
    
    number_of_edge_positions = length( subscripts_at_edge( :, 1 ));
    
    % loop through the positions in each edge
    for edge_position_index = 1 : number_of_edge_positions
    
        % find the linear index of the center pixel of the sphere that defines the edge at this
        % position
        edge_position_linear_index = sub2ind( size_of_image,                                ...
                                              subscripts_at_edge( edge_position_index, 1 ), ...
                                              subscripts_at_edge( edge_position_index, 2 ), ...
                                              subscripts_at_edge( edge_position_index, 3 )  );        
        
        % adding up and averaging later
          weighting_image(   edge_position_linear_index                                ...
                           + strel_indices{ subscripts_at_edge( edge_position_index, 4 )}) ...
        = weighting_image(   edge_position_linear_index                                ...
                           + strel_indices{ subscripts_at_edge( edge_position_index, 4 )}) ...
        + 1 ;

            strands_image(   edge_position_linear_index                                ...
                           + strel_indices{ subscripts_at_edge( edge_position_index, 4 )}) ...
        =   strands_image(   edge_position_linear_index                                ...
                           + strel_indices{ subscripts_at_edge( edge_position_index, 4 )}) ...
        + strand_fill_values{ edge_index }( edge_position_index );

        
    end % for position
end % for trajectory

strands_image = strands_image ./ weighting_image ;

mat2tif( uint16( visual_contrast * strands_image ),    spheres_visual_file )

end % function