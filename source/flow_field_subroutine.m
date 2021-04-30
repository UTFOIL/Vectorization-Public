function [ flow_field, tissue_type_image ] = flow_field_subroutine( path_to_flow_field_export )

%% for Chakameh, 4/25/18, rendering vascular vectors, from Sam
% this script takes the positions list and puts a sphere on a matrix whose uniform value is equal to
% the corresponding laplacian value of that position.  The output of the 'visualize' functions
% below are TIFF files to be loaded into FIJI or imageJ and can then be overlayed over the original.
%
% for your purposes you may want to eliminate the laplacian as this is a contrast metric for telling
% us how "good" that vector is.  Also, you might want to adapt the visual functions to output the
% 3D image matrix before it is converted into a TIFF.
%
% tangent vector calculations should be performed on the positions in the trajectory form, as this
% is the form that assures that the index within the contents of each position in the cell array
% corresponds to a trajectory index.
%
% Enjoy! 
% SAM 4/25/18

%% V200
% In the second version, we put the contrast thresholding in this script (in V1 it was done before
% you load the data).  We also switched from a list view of the vector objects to a grouped cell
% array view where each cell array bin represents a single edge (a vessel segment between two
% vertices), and inside the bin is the list of physical positions in physical order from one vertex
% to another.
%
% The shortcomings of grouping the vector objects into edges were two-fold: 
% 
% 1) we couldn't easily tell which vertices were bifurcations 
%
% 2) we couldn't easily calculate spatial derivatives with respect to vessel axial direction (our
% approximation for local blood flow direction).
%
% SAM circa 7/5/18 (retro-dated on 7/12/18)

%% V300
% In the third version, we organized the edges into objects called "strands" and "junctions" which
% are both collections of edges arranged head to tail so that following the object down its list of
% vectors will take you along a path from one vertex to another across (possibly many) other
% vertices and edges.  What constitutes inclusion in a strand is that all edges in a strand will be
% either the first or second best edge for all the vertices in that strand and that those vertices
% will be strongly (i.e. mutually) connected by the edges.
%
% After strands are found, all edges involved in strands are erased and then all vertices included
% in strands (excluding the end vertices who have only one neighbor) are granted their next best
% edge.  Any new connections formed in this way that take us from one strand to a different strand
% are called junctions.
%
% Combining the strand and junction networks and then identifying vertices with three edges gives us
% the location of vertices that might represent bifurcations.  The contrast metric that would be
% appropriate for determining if the bifurcation is authentic or artifactual would be the max_energy
% for the edge(s) in the junction object that connects to the bifurcation vertex.
%
% The strand and junction organization of edges is encoded in the variables edge_indices_in_strands
% and edge_indices_in_junctions in the chakameh_V300.mat file.  Likewise there are variables that
% encode the organization of the vertices into strands and junctions.  These encoding variables
% index the edges or vertices used to construct each strand or junction objects.  The list of
% vectors to which these encoding variables index are the edge_subscripts variable and the
% vertex_subscript variable.  If when an edge is called in the formation of a strand its order needs
% to be flipped before correctly fitting into the list, this is will be encoded in the variable
% edge_backwards_in_strands (or edge_backwards_in_junctions). 
%
% V400: specific to 170315 dataset: max energy variable name changed to mean energy.  contrast limit
% set to zero (no thresholding).  Fixed flow-field and tissue type rendering resolution at 1 um^3
% voxel size.  Made the interpolation step more robust by adding objects when the interpolation
% creats gaps between objects.  Also added a manual edge for a large surface vessel before this
% script but will not do the smoothing operation on that edge in this script.  SAM 1/2/19
%
% V500: Removed the special treatment for manually added edges, as those methods have been improved
% and now the manually added edges have the same data format as the automatically generated ones.
% Therefore this script should be compatible with any input dataset, provided the .mat file is
% correctly located and loaded. SAM 2/27/19

%% MAIN SCRIPT
% loading, naming, and resolution selection

% load('flow_field_export_network_190303-210838_Fused_T123_T456_Raw_pad30.mat')

load( path_to_flow_field_export )

% root_path = 'change this to your desired root' ;
root_path = [ pwd, filesep ];

% time_stamp = char( datetime('now', 'TimeZone', 'local', 'Format', 'yyMMdd-HHmmss' )); 
% 
% file_base_name = 'edges_' ;
% 
% input: loaded from the .mat file
% file_name = [ file_base_name, time_stamp ];

path_to_tissue_type_visual = [ root_path, file_name, '_tissue_types.tif' ];
path_to_centerlines_visual = [ root_path, file_name,  '_centerlines.tif' ];   

flow_visual_files = cell( 3, 1 );

flow_visual_files{ 1 }      = [ root_path, file_name,      '_flow_y.tif' ];
flow_visual_files{ 2 }      = [ root_path, file_name,      '_flow_x.tif' ];
flow_visual_files{ 3 }      = [ root_path, file_name,      '_flow_z.tif' ];

desired_cubic_voxel_length_in_microns = 2 ;

% digital resolution enhancement: fixing voxel size at 1 um ^ 3
resolution_factor = microns_per_pixel_xy / desired_cubic_voxel_length_in_microns ; % double this factor to double the rendering resolution in each spatial dimension

% note to programmer (consider allowing for rendering resolution enhancement of radius as well as of
% space and possibly cut down on the number of scales required in processing) SAM 7/13/18

% %% smoothing and approximating directions
% % Arrange edges into strands to do smoothing and differentiation along the strand axis.
% %
% % calclation of vessel directions should be done before threseholding to minimize edge effects on
% % derivative calculation
% 
% % input: loaded from the .mat file
% % smoothing_kernel_sigma_to_lumen_radius_ratio = .5 ;
% 
% % note: make the tissue type cutoff decision by picking the largest radius in microns that will be
% % allowed into each tissue category.  To view the radii in microns and select the  cutoff indices,
% % look at the following look_up_table
% 
% % input: loaded from the .mat file
% % lumen_radius_in_microns_range = lumen_radius_in_pixels_range( :, 1 ) * microns_per_pixel_xy ;
% 
% % input: loaded from the .mat file
% % % % capillaries vs larger vessels
% % tissue_type_cutoffs_in_microns = [ 5.5 ]; % microns of the radius of the vessel SAM 8/7/18
% 
% % tissue_type_cutoffs_in_microns = lumen_radius_in_pixels_range( :, 1 ) * microns_per_pixel_xy ; % microns of the radius of the vessel SAM 8/7/18
% 
% % below tissue_type_cutoffs_in_microns( 1 ) is tissue type 1, between
% % tissue_type_cutoffs_in_microns( 1 ) and tissue_type_cutoffs_in_microns( 2 ) is tissue type 2, ...
% % bigger than tissue_type_cutoffs_in_microns( end ) is the last tissue type. There is no limit to
% % number of tissue types. Tissue type zero is background (extravascular).
% 
% [ strand_subscripts, strand_energies ]                                                                           ...
%         = get_strand_objects( edge_subscripts, edge_energies, edge_indices_in_strands, edge_backwards_in_strands );
%     
% strand_space_subscripts = cellfun( @( x ) x( :, 1 : 3 ), strand_subscripts, 'UniformOutput', false );
% strand_scale_subscripts = cellfun( @( x ) x( :,   4   ), strand_subscripts, 'UniformOutput', false );
% 
% [ strand_space_subscripts, strand_scale_subscripts, strand_energies ]  ...
%                              = smooth_edges_V2( strand_space_subscripts, strand_scale_subscripts, ...
%                                                 strand_energies,                                  ...
%                                                 smoothing_kernel_sigma_to_lumen_radius_ratio,     ...
%                                                 lumen_radius_in_microns_range, microns_per_voxel  );
%     
% [ vessel_directions ] = get_vessel_directions_V3( strand_space_subscripts, microns_per_voxel );

%% Assigning tissue types to strands:

% converting the tissue type cutoff from microns to an index in the lumen_radius_in_pixels_range LUT
number_of_cutoffs = length( tissue_type_cutoffs_in_microns );

tissue_type_cutoffs = zeros( 1, number_of_cutoffs );

for tissue_type_cutoff_index = 1 : number_of_cutoffs
    
    tissue_type_cutoffs( tissue_type_cutoff_index )                                                 ...
                               = find(   tissue_type_cutoffs_in_microns( tissue_type_cutoff_index ) ...
                                       / microns_per_pixel_xy                                       ...
                                       > lumen_radius_in_pixels_range( :, 1 ),                      ...
                                       1, 'last'                                                    );
    
end % tissue type FOR

% concatenate infinity onto the tissue type cutoffs for largest tissue assignment work later.
tissue_type_cutoffs( end + 1 ) = Inf ;

tissue_types = cellfun( @( x ) find( tissue_type_cutoffs >= median( x ), 1 ) + 1, strand_scale_subscripts );
% add 1 to the FIND call so that capillary will be tissue type 2 and vessel tissue type 3

%% interpolation of vectors

% guarantee isotropic interpolation
resolution_factors = resolution_factor * [ 1, 1, z_per_xy_length_of_pxl_ratio ];

strand_subscripts = cellfun( @( x, y ) [ x, y ], strand_space_subscripts, strand_scale_subscripts, 'UniformOutput', false );

[ size_of_image, lumen_radius_in_pixels_range, strand_subscripts, vessel_directions ] ...
                                    = resample_vectors( lumen_radius_in_pixels_range, [ resolution_factors, 2 ], strand_subscripts, size_of_image, vessel_directions );

strand_subscripts = cellfun( @round, strand_subscripts, 'UniformOutput', false );

%% visualization and flow field output

[ mean_strand_energies ] = get_edge_metric( strand_energies );
                         
[ flow_field, tissue_type_image ]                                                                   ...
           = render_flow_field_V3( strand_subscripts, vessel_directions, mean_strand_energies, ...
                                   tissue_types, lumen_radius_in_pixels_range,                      ...
                                   size_of_image, path_to_tissue_type_visual,                       ...
                                   path_to_centerlines_visual, flow_visual_files                    );
                            
% %% locations and sizes of bifurcation vertices
% 
% bifurcation_subscripts = vertex_subscripts( bifurcation_vertices, : );