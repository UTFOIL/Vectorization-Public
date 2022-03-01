function [ flow_field, tissue_type_image ]                                                             ...
                        = render_flow_field_V3( edge_subscripts, vessel_directions, mean_edge_energies, ...
                                                tissue_types, lumen_radius_in_pixels_range,            ...
                                                size_of_image, tissue_type_visual_file,                ...
                                                centerline_visual_file, flow_visual_files              )
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
% function name changed from visualize_edges to render_flow_field SAM 7/13/18. This function is part
% of the effort simulate blood flow direction based on 2P images.  This function differs from the
% previous version in that it outputs a 4D image.  The first three dimensions are y, x, and z
% spatial coordinates in the image volume and te fourth dimension has indices of 1, 2, and 3
% standing for the y, x, and z components of the flow field at that location. SAM 7/13/18 
%
% Note that this function could be improved by inputting the edge_energies variable (instead of the
% max_edge_energies variable) and sorting each individual vector instead of just sorting at the edge
% level. SAM 7/13/18
%
% Note that this function could be improved by using the vessel direction to change the fill shape
% from a sphere to a disk.  This would alleviate the overwriting artifacts that we are sure to
% observe in the rendering attempted below. SAM 7/13/18
%
% V2 in which the tissue type is an input and writes to a TIFF image instead of the
% max_edge_energies, and it is also output as a matrix. SAM 7/25/18
%
% V3 in which the distance of each point in the strel from the centerline is calculated about every
% centerline voxel, and this is used to weight the rendered flow-field. (Weighting is parabolic with
% respect to distance from the centerline.  Weighting has average value of 1 over the
% cross-sectional area of vessel normal to the centerline axis, where the average is a mean over an
% estimated number of pixels.  (This noramlization calculation currently requires isotropic sampling
% in the 3D image.  To do the averaging correctly for different orientations of vessels in an
% anisotropic voxel lattice would require more attention.) SAM 9/19/18
%
% V3 also, removed the sigma_per_size factor and renamed pixels_per_sigma_range to
% lumen_radius_in_pixels_range

% !!! consider writing the flow field not in overwriting of cumulative fasion but in a cumulative
% average fashion. (flow field cumulative and a weighting factor cumulative (which has a disk or
% ellipsoidal shape containing the cross-section of interest). see the depth and direction
% vissualizations for example). Keep in mind that two adjacent vessels with flows 180 degrees to
% each other will average to no flow !!!

% these dimensions will be permuted later to be consistent with the version comment above
flow_field = zeros([ 3, size_of_image ]); % ??? where does it get permuted???

tissue_type_image = zeros( size_of_image );
centerlines_image = zeros( size_of_image );

number_of_image_voxels = prod( size_of_image );

number_of_edges = length( mean_edge_energies );

% sorting trajectories by mean energy in descending order (most negative at end) and then sorting
% trajectories by mean size in ascending order (largest at end)
mean_edge_sizes = cellfun( @( x ) mean( x( :, 4 )), edge_subscripts );

[ ~, indices_sorted_by_energy ] = sort( mean_edge_energies,                          'descend' );
[ ~, indices_sorted_by_size   ] = sort( mean_edge_sizes( indices_sorted_by_energy ),  'ascend' );

sorted_indices = indices_sorted_by_energy( indices_sorted_by_size );

edge_subscripts    =    edge_subscripts( sorted_indices );
mean_edge_energies = mean_edge_energies( sorted_indices );
vessel_directions  =  vessel_directions( sorted_indices );
tissue_types       =       tissue_types( sorted_indices );

% pre-calculating all of the ellipsoidal structuring elements to be used to paint in the image
number_of_scales = size( lumen_radius_in_pixels_range, 1 );

structuring_element_linear_indexing         = cell( number_of_scales, 1 );
structuring_element_directions_from_center  = cell( number_of_scales, 1 );
structuring_element_distances_from_center   = cell( number_of_scales, 1 );

normalizing_coefficient = zeros( number_of_scales, 1 );
middle_index            = zeros( number_of_scales, 1 );

for scale_index = 1 : number_of_scales

    % find all pixel locations within the ellipsoid radii from the vertex position
    [ structuring_element_linear_indexing{ scale_index }, structuring_element_subscripts ]                      ...
            = construct_structuring_element( lumen_radius_in_pixels_range( scale_index, : ), size_of_image );
            
    structuring_element_distances_from_center{ scale_index } = sum( structuring_element_subscripts .^ 2, 2 ) .^ 0.5 ;

    structuring_element_directions_from_center{ scale_index } = structuring_element_subscripts ./ structuring_element_distances_from_center{ scale_index } ;

    middle_index( scale_index ) = ( length( structuring_element_distances_from_center{ scale_index }) + 1 ) / 2 ;

    % normalizing coefficient (for normalizing mean speed across the cross-sectional area to be
    % one) as written here only works in isotropic voxel lattice:
    normalizing_coefficient( scale_index ) = 2 / lumen_radius_in_pixels_range( scale_index, 1 ) ^ 2 ;
        
end % scale FOR

% loop through the edges
for edge_index = 1 : number_of_edges
    
    subscripts_at_edge = edge_subscripts{ edge_index };
    
    number_of_edge_positions = length( subscripts_at_edge( :, 1 ));
    
    % loop through the positions in each edge
    for edge_position_index = 1 : number_of_edge_positions
    
        % find the linear index of the center pixel of the sphere that defines the edge at this
        % position
        edge_position_linear_index = sub2ind( size_of_image,                                ...
                                              subscripts_at_edge( edge_position_index, 1 ), ...
                                              subscripts_at_edge( edge_position_index, 2 ), ...
                                              subscripts_at_edge( edge_position_index, 3 )  );        
        
        % label the centerline in a cumulative fashion so that it will get brighter if it is doubly
        % labeled
        centerlines_image( edge_position_linear_index )         ...
            =   centerlines_image( edge_position_linear_index ) ...
              - mean_edge_energies( edge_index );
                          
        % label the sphere in an overwriting fashion so that only the lowest energy edges will shine
        % through in a multiply labeled region.  (Remember that we sorted by max energy attained
        % before this FOR loop so the later edges to be written are lower in energy). 
        tissue_type_image( max( min(   edge_position_linear_index                                                                ...
                                     + structuring_element_linear_indexing{ subscripts_at_edge( edge_position_index, 4 )},       ...
                                     number_of_image_voxels                                                                ),    ...
                                1                                                                                             )) ...
            = tissue_types( edge_index );
                
        % label flow field in overwriting fashion.  Fill everywhere that would be filled in the
        % spheres image above, but instead of writing a single value in, fill with the 3D vector
        % representing the vessel directions.
        
        cosine_of_angle_from_vessel_axis = sum(      vessel_directions{ edge_index }( edge_position_index, : )                                     ...
                                                 .* structuring_element_directions_from_center{ subscripts_at_edge( edge_position_index, 4 )}, 2 ) ;
                                             
        cosine_of_angle_from_vessel_axis( middle_index( subscripts_at_edge( edge_position_index, 4 ))) = 0 ;
                                             
        abs_sine_of_angle_from_vessel_axis = ( 1 - cosine_of_angle_from_vessel_axis .^ 2 ) .^ 0.5 ;
        
        distance_from_vessel_axis = structuring_element_distances_from_center{ subscripts_at_edge( edge_position_index, 4 )} ...
                                  .* abs_sine_of_angle_from_vessel_axis ;
        
%         % restricting the strel to the elements that are not in the imaginary cone inside of 30
%         % degrees from the vessel axis.
%         elements_of_strel_outside_coaxial_cone = find( abs_sine_of_angle_from_vessel_axis >= ( 1 / 2 ));
%         
%         % restricting the strel to the elements that are not in the imaginary cone inside of 45
%         % degrees from the vessel axis.
%         elements_of_strel_outside_coaxial_cone = find( abs_sine_of_angle_from_vessel_axis >= ( 1 / 2 ) ^ 2 );
        
        % no restriction applied to the elements used to paing the flow field image.  Each
        % centerline will paint in a sphere of radius equal to the vessel radius with a parabolic
        % velocity profile and a rotational symmetry about the vessel axis
        elements_of_strel_outside_coaxial_cone = find( abs_sine_of_angle_from_vessel_axis >= 0 );
        
        distance_from_vessel_axis = distance_from_vessel_axis( elements_of_strel_outside_coaxial_cone );
        
        % converting distance from vessel to a velocity of flow magnitude by assuming a parabeloid
        % flow profile across cross-sections of vessel (consistent with laminar flow assumption).
        % This profile will be normalized to have total veloticy of 1 across the cross-section
        % normal to the vessel direction.     
        flow_magnitude = normalizing_coefficient( subscripts_at_edge( edge_position_index, 4 )    )     ...
                       * (       lumen_radius_in_pixels_range( subscripts_at_edge( edge_position_index, 4 ), 1 ) ^ 2 ...
                            -  distance_from_vessel_axis .^ 2                                           )';
        
                
%         flow_field_outside_axial_cone = 
        
        number_of_pixels_in_strel = size( distance_from_vessel_axis, 1 );        
        
        flow_field( :, max( min( edge_position_linear_index                                                                                                                        ...
                                    + structuring_element_linear_indexing{ subscripts_at_edge( edge_position_index, 4 )}( elements_of_strel_outside_coaxial_cone ),                ...
                                 number_of_image_voxels                                                                                                             ),             ...
                            1                                                                                                                                          ))          ...
            =   vessel_directions{ edge_index }( edge_position_index, : )' .* ones( 3, number_of_pixels_in_strel ) .* flow_magnitude ;
        
    end % for position
end % for edge

mat2tif( uint16( tissue_type_image ), tissue_type_visual_file )

mat2tif( uint16( centerlines_image ),  centerline_visual_file )

clear centerlines_image

for component_ordinate = 1 : 3
    
    % publish a tif for each component of the flow field separately.  Scale the flow field values to
    % have radius 2000 and center 1000 for rendering in a uint16 scale.
    mat2tif( 1000 * squeeze( flow_field( component_ordinate, :, :, : )), ...
                       flow_visual_files{ component_ordinate }           )
    
end % component FOR

end % function