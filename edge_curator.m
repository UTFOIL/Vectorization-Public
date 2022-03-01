function [ mean_edge_energies, edges2vertices, edge_energies, edge_space_subscripts,                ...
           edge_scale_subscripts, index_image, true_edges, displayed_edges, deleted_edges,...
           vertex_space_subscripts, vertex_scale_subscripts                               ]         ...
            = edge_curator( edge_energies, edge_space_subscripts, edge_scale_subscripts,         ...
                              edges2vertices, vertex_space_subscripts,                              ...
                              vertex_scale_subscripts, lumen_radius_in_microns_range,               ...
                              length_dilation_ratio_vertices, length_dilation_ratio_edges,          ...
                              microns_per_pixel, path_to_original_data, path_to_saved_curation,     ...
                              path_to_energy_data, intensity_limits, energy_limits, threshold_percent_overlap )
%% vertex_curator SAM 7/19/18
% This function creates and runs the graphical user interface for the manual curation of vertices.
% It shows a projection of the vertices and (inverted) raw data overlay.  User has control over the
% depth and thickness of the volume to be projected.
%
% note: make the display overlay occur by having two different handles to the two images to avoid
% unnecessary rendering.  Also consider making alpha value proportional to - vertex_energies
%
% V2: coming back to this GUI after a few weeks and needed to make a new version. New changes will
% include a clear negatives button and a repopulate button and a place vertex button.  Minor changes
% will be to the color scheme and minimap.  This version will be blue for positive and red for
% negative.  Minimap will be a 3D rectangular prism.
%
% V3: Many edits made by William Andrew Sikora (WAS) in October, 2018. 
% (10/20/18)
% -Made histogram figures invisible
% -Removed all uicontrol elements from main window but toggle, sliders, and slider labels
% -Locked buttons and sliders into their final locations with �normalized� �Position� fields
% -Adjusted display_axes �Position� field to give more room for buttons
% -Added toggling variable under global list (0 = false, 1 = true) and if statement to toggle_callback to make the toggle button toggleable (using the added variable, of course)
% -Implemented what I *think* displays all the circles as logical, rather than scaled values. But I will need to verify this with you
% 
% (10/24/18)
% -removed �toggling� global variable and all changes to toggle button made on 10/20/18
% -implemented an ideal toggling system by simply disabling the toggle button when the recursive toggle_callback is in operation
% -added minimap figure
% -made minimap resolution 100 by 100 to save on computation time (very slow at full res)
% -got it to render as we discussed
% -made minimap axes locked square just like display_axes
% -added color picker functionality for minimap (control+f �color selection tool�)
% 
% (10/24/18) Part 2
% -added ticks to the axes of minimap
% -guaranteed that data on minimap is displayed at equal scale on all axes (i think. I used daspet command so search for that to double check)  
%
% Edits made by SAM on 11/7/18:
%
% 1. rearranged buttons
% 2. made all text font sizes normalized
% 3. automatically saves upon exiting, and the save function has more variables being saved now
% 4. energy_limits is an input and vertices are displayed in many shades
% 5. fixed bugs in saving, loading, undo and redo that prevented proper restoring of the histograms
% and slider positions upon loading a stored curation.
% 6. added red and cyan boxes to energy histogram to visualize the energy threshold
% 7. Extended backup functionality to threshold drawing and toggling.  Now thresholding and toggling 
% can be undone with the undo button.
%
% V4: Depth information was added to the reference volume figure. SAM 11/8/18
%
% V5: Edits made by WAS on 11/15/18:
%
% 1. Make minimap reflect currect state of thresholding
%   a. Added update_minimap_threshold function to track thresholding
%   b. Changed update_minimap function to implement changes
% 2. Changed 'William Drew Sikora (WDS)' to 'William Andrew Sikora (WAS)'
% 3. Implemented undo/redo and save/load for changes made in part 1 of these edits
%
% Function name was changed from vertex_curator to edge_curator.  Inputs changed from vertices to
% edges.  SAM+WAS 11/16/18
%
% V2, in which the vertices are also displayed all in the same shade of yellow.  SAM + WAS 12/6/18
%
% V3: Edits made by WAS on 12/13/18:
%
% 1. Reflected thresholding on histogram. 
%   -Note: Histogram only paints thresholds in current view, but I commented 
%          code that can make it do all thresholds presently applied to index image
%   -MORE IMPORTANT NOTE: because threshold_mat_3D cannon't be built up
%          from previously curated data, this makes vertex_curator_V6 incompatible
%          with data partially curated (or at least thresholded) data.
% 2. Made colormap matrix and color picking variables global to use with 
%    the histogram as well as the minimap
% 3. Added load_colormap function and called it just before first render
% 4. Saved 3D representation of all thresholds at reolution of
%    'size_of_minimap'. Named it 'threshold_mat_3D'. SEE 'MORE IMPORTANT
%    NOTE' IN EDIT 1!
% 5. Changed 'update_minimap_threshold' to 'update_threshold_visualization_matricies'
%
% V4: paint_index_image function updated to include edge-edge volume exclusion WAS 2/11/19
%
% V5: converted mean_edge_energies to a z-score.  Edited V4 changes (removed the multiple edge
% conflict test) SAM 2/13/19


%% Initializations         
   
initial_load_variables = {  'edge_energies', ...
                    'edge_space_subscripts', ...
                    'edge_scale_subscripts', ...
                           'edges2vertices', ...
                'max_number_of_added_edges', ...
                         'intensity_limits', ...
                            'energy_limits'  };

max_number_of_added_edges = 1000 ;

listing_for_path_to_saved_curation = dir([ path_to_saved_curation, '.mat' ]);

if ~ isempty( listing_for_path_to_saved_curation )

    answer = questdlg( 'Load the previously saved curation?', 'Edge Curator' );

    if strcmp( answer, 'Yes' )
            
        load( path_to_saved_curation, initial_load_variables{ : });
        
        % these are placeholder variables, will be rewritten later, added objects removed for
        % accounting purposes only
        
        % cell vectors
        edge_energies         =         edge_energies( max_number_of_added_edges + 2 : end );
        edge_space_subscripts = edge_space_subscripts( max_number_of_added_edges + 2 : end );      
        edge_scale_subscripts = edge_scale_subscripts( max_number_of_added_edges + 2 : end );  
        
        % numerical vector and matrix
        edges2vertices        =        edges2vertices( max_number_of_added_edges + 2 : end, : ); 
%         degrees_of_edges      =      degrees_of_edges( max_number_of_added_edges + 2 : end    ); 
        
        load_upon_startup = true ;
    
    elseif isempty( answer ) || strcmp( answer, 'Cancel' ) % The "X" or 'Cancel'

              edge_energies = [ ]; % sam 1/27/22
        
%         mean_edge_energies = get_edge_metric( edge_energies );
              mean_edge_energies = [ ];
          
                index_image = [ ];
        
                 true_edges = [ ];
            displayed_edges = [ ];
              deleted_edges = [ ];
              
        return
        
    else % 'No'
        
        load_upon_startup = false ;
        
    end
else
    
    load_upon_startup = false ;

end

vertex_color =  ones( 1, 1, 3 );
% vertex_color( 3 ) = 0 ; % yellow/gold
% vertex_color( 1 ) = 0 ; % green+blue=aqua
% vertex_color( 1 ) = 0.5 ; % green+blue=aqua
% % gray
% vertex_color([ 1, 3 ]) = 0.5 ; % bright green
% vertex_color( 2 ) = 0 ; % purple
% vertex_color( : ) =  [ 0.5, 0.5, 1 ]; % bb blue
% vertex_color( : ) =  [ 1, 0, 1 ]; % bb purple
% vertex_color( : ) =  [ 0, 0, 1 ]; % blue
% vertex_color( : ) =  [ 0.3, 0.7, 1 ]; % edge color (true)
vertex_color( : ) =  [ 0.7, 0.3, 1 ]; % hot blue

centerline_color = zeros( 1, 1, 3 ); 
centerline_color( 2 ) = 1 ; % green
% centerline_color( : ) = [ 0.7, 1, 0.3 ]; % hot green
centerline_color( : ) = [ 0.3, 1, 0.7 ]; % cool green

    %% Figures                                  
origin = [ 0.1, 0.4 ];
    
scaling = 0.7 ;
    
bottom_buffer = 0.05 ;

left_buffer = bottom_buffer / 2 ;

spacer = 0.005 ;

histo_widths = 0.19 ;

histo_heights = 0.45 ;

half_screen = 1 - left_buffer - 2 * histo_widths - spacer ;

filesep_indices = regexp( path_to_saved_curation, filesep );

edge_set_name =  path_to_saved_curation( filesep_indices( end ) + 1 : end );

digit_indices = regexp( edge_set_name, '\d' ); % \d searches for digits

dataset_name =  edge_set_name( digit_indices( 12 ) + 2 : end ); % 12 digits in the date-time stamp

display_figure = figure( 'Name'     ,[ dataset_name, ': Edge Curator: Volume Display' ], ...
                         'Units'    , 'normalized',                     ...
                         'SizeChangedFcn',@main_figure_size_update,     ...
                         'Visible', 'off',                              ...
                         'OuterPosition' , scaling * [ origin + [ left_buffer, bottom_buffer ], [ half_screen - spacer - left_buffer , 2 * histo_heights + spacer ]]);
    
display_figure.Units = 'pixels' ; 

% % force square dimensions
% display_figure.Position([ 3, 4 ]) = min( display_figure.Position([ 3, 4 ]));

display_axes = axes( 'OuterPosition', [0.075, 0.15, 0.85, 0.85]);
display_axes.Units = 'pixels';
    
h = pan( display_figure );
h.ActionPostCallback = @update_xy_zoom ;

h = zoom( display_figure );
h.ActionPostCallback = @update_xy_zoom ;
h.Enable             =            'on' ;

% intensity_histo figure and axes initializations
intensity_histo_figure = figure( 'Name'     , [ dataset_name, ': Edge Curator: Intensity Histogram' ], ...
                                 'Units'    , 'normalized',                          ...
                                 'OuterPosition' , scaling * [ origin + [ half_screen, bottom_buffer ], [ histo_widths, histo_heights ]]);

% intensity_histo_figure.Units = 'pixels' ; 
% 
% % force square dimensions
% intensity_histo_figure.Position([ 3, 4 ]) = min( intensity_histo_figure.Position([ 3, 4 ]));
% 
% intensity_histo_figure.Visible = 'off' ;

intensity_histo_axes = axes( 'OuterPosition', [ 0.1, 0.3, 0.8, 0.6 ]);

zoom( intensity_histo_figure, 'yon' )

% intensity_histo figure and axes initializations
energy_histo_figure = figure( 'Name'     , [ dataset_name, ': Edge Curator: Energy Histogram' ], ...
                              'Units'    , 'normalized',                       ...
                              'OuterPosition' , scaling * [ origin + [ half_screen + histo_widths + spacer, bottom_buffer ], [ histo_widths, histo_heights ]]);

% energy_histo_figure.Units = 'pixels' ; 
% 
% % force square dimensions
% energy_histo_figure.Position([ 3, 4 ]) = min( energy_histo_figure.Position([ 3, 4 ]));
% 
% energy_histo_figure.Visible = 'off' ;

energy_histo_axes = axes( 'OuterPosition', [ 0.1, 0.3, 0.8, 0.6 ]);

zoom( energy_histo_figure, 'yon' )

% % mini_map figure and axes initializations
minimap_figure = figure( 'Name'     , [ dataset_name, ': Edge Curator: Volume Map' ], ...
                         'Units'    , 'normalized',                       ...
                    'SizeChangedFcn',@minimap_figure_size_update,         ...
                           'Visible','off',                               ...
                         'OuterPosition' , scaling * [ origin + [ half_screen, bottom_buffer + histo_heights + spacer ], [ 2 * histo_widths + spacer, histo_heights ]]);
                           
minimap_figure.Units = 'pixels' ; 
 
% minimap_figure.Position([ 3, 4 ]) = min( minimap_figure.Position([ 3, 4 ]));

minimap_axes = axes(minimap_figure, 'OuterPosition', [0.2, 0.2, 0.7, 0.7]);
minimap_axes.Units = 'pixels';
view(minimap_axes,3);
    %% Variables                                
% depricated variables:
y_limits = [ ];
x_limits = [ ];
z_limits = [ ];

z_thickness = [ ];
z_depth     = [ ];
z_range     = [ ];

depth_dimension = [ ];

% initialize global variables
depth_range                 = [ ];
FOV_limits                  = [ ];
edge_centers                = [ ];
index_image_2D              =  1 ;
vertex_index_image_2D       =  1 ;
centerlines_image_2D        = [ ];
index2intensity             = [ ];
index_image_crop            = [ ];
centerlines_image_crop      = [ ];
truth_image_uint8           = [ ];
in_view_edges_z             = [ ];
% gaussian_weight_in_z        = [ ];
original_image_crop_2D      = [ ];
adjusted_original_image_2D      = [ ]; 
are_intensity_limits_binarized  = true ;
edge_display_method = 1 ;

is_intensity_inverted = true ;

counter_toggle              = 0 ;
counter_threshold           = 0 ;
counter_add_vertex          = 0 ;
counter_add_edge            = 0 ;
elapsed_time                = 0 ;

tstart = tic ;

projection_dimension = 3 ;

dimension_order_matrix ...
    = [ 2, 3, 1 ; ...% y-projection
        3, 1, 2 ; ...% x-projection
        1, 2, 3   ]; % dimension permutation to project in different dimensions

dimension_letters = { 'Y', 'X', 'Z' };

dimension_order = dimension_order_matrix( projection_dimension, : );

dimension_order = [ 1, 2, 3 ]; % dimension permutation to project in different dimensions

vertex_structure_positions_linear_indexing = [ ];
  edge_structure_positions_linear_indexing = [ ];
  
  edge_positions_linear_indices = [ ];
   
size_of_minimap             = 100 ;
energy_threshold            = 0   ;

threshold_mat_3D            = energy_threshold.*ones(size_of_minimap,size_of_minimap,size_of_minimap);
xy_threshold_mat            = energy_threshold.*ones(size_of_minimap);
xz_threshold_mat            = energy_threshold.*ones(size_of_minimap);
yz_threshold_mat            = energy_threshold.*ones(size_of_minimap);
colormap_mat                = [];
threshold_color_range       = [0.5 1];
out_of_bounds_color         = 0.0;
border_color                = 0.25;
in_view_color               = 0.49;

% last_vertex_index = max( edges2vertices( : ));
% last_vertex_index = length( vertex_space_subscripts );

% add a row for the background when indexing by the index image (that image has ones on background)
number_of_edges = length( edge_energies );

number_of_added_edges = 0 ;

[ mean_edge_energies ] = get_edge_metric( edge_energies );
    
mean_edge_energies    = [ zeros( max_number_of_added_edges + 1, 1 ); mean_edge_energies ];

edges2vertices        = [ zeros( max_number_of_added_edges + 1, 2, 'uint32' ); edges2vertices ];

% degrees_of_edges      = [  ones( max_number_of_added_edges + 1, 1, 'uint16' ); degrees_of_edges ];

edge_space_subscripts = [ cell( max_number_of_added_edges + 1, 1 ); edge_space_subscripts ];
edge_scale_subscripts = [ cell( max_number_of_added_edges + 1, 1 ); edge_scale_subscripts ];                      
edge_energies         = [ cell( max_number_of_added_edges + 1, 1 ); edge_energies         ];

edge_space_subscripts( 1 : max_number_of_added_edges + 1 ) = { zeros( 1, 3 )};
edge_scale_subscripts( 1 : max_number_of_added_edges + 1 ) = { 1 };
        edge_energies( 1 : max_number_of_added_edges + 1 ) = { 0 };
                            
length_of_edge_lists = number_of_edges + max_number_of_added_edges + 1 ;

% vertex_code_in_index_image_2D = last_vertex_index + length_of_edge_lists + 1 ;

% intialize logical vectors of class assignment
   true_edges          = ones(  length_of_edge_lists, 1, 'logical' );
in_view_edges_xy       = ones(  length_of_edge_lists, 1, 'logical' );
in_view_edges          = ones(  length_of_edge_lists, 1, 'logical' );
deleted_edges          = zeros( length_of_edge_lists, 1, 'logical' );
displayed_edges        = zeros( length_of_edge_lists, 1, 'logical' );

% mark the background and place holders as deleted and not displayed or in view
deleted_edges(     1 : max_number_of_added_edges + 1 ) = true  ;
in_view_edges(     1 : max_number_of_added_edges + 1 ) = false ;
in_view_edges_xy(  1 : max_number_of_added_edges + 1 ) = false ;

% % make analogy to the following statement to replace the prior statement
% deleted_vertices(    vertex_space_subscripts( :, 1 ) <= 0 ) = true  ;
% in_view_vertices(    vertex_space_subscripts( :, 1 ) <= 0 ) = false ;
% in_view_vertices_xy( vertex_space_subscripts( :, 1 ) <= 0 ) = false ;

maximum_number_of_backups = 20 ;

backup_ordinate = 0 ;

% pull original from h5 files       
original_image = h52mat( path_to_original_data );

size_of_image = size( original_image );

% % main display image initializtion, background is coded by last_vertex_index plus one
% index_image = 1 + last_vertex_index * ones( size_of_image, 'uint32' );

% main display image initializtion, background is coded by zero
      index_image = zeros( size_of_image, 'int32' );

centerlines_image = zeros( size_of_image, 'logical' );

number_of_pixels_2D = size_of_image( 1 ) * size_of_image( 2 ) ;

% initialize axis limits
FOV_limits( 1, 1 : 2 ) = [ 1, size_of_image( 1 )];
FOV_limits( 2, 1 : 2 ) = [ 1, size_of_image( 2 )];

display_axes.XLim = FOV_limits( 2, : );
display_axes.YLim = FOV_limits( 1, : );

% initialize depth_of_view and view_thickness for initial viewing and slider limits
depth_of_view      = round( size_of_image( 3 ) / 2 );
view_thickness_max =        size_of_image( 3 );

% view_thickness = min( 5, view_thickness_max );     
view_thickness = 0 ;

autosave_variables = { 'true_edges', ...
                   'deleted_edges', ...
                 'displayed_edges', ...
                      'FOV_limits', ...
                         'depth_of_view', ...
                     'view_thickness', ...
                  'edges2vertices', ...
           'edge_space_subscripts', ...
           'edge_scale_subscripts', ...
           'vertex_space_subscripts', ...
           'vertex_scale_subscripts', ...
           'edge_energies'        , ...
      'mean_edge_energies'        , ...
                'intensity_limits', ...
                   'energy_limits', ...
           'number_of_added_edges', ...
                'energy_threshold', ...
                'xy_threshold_mat', ...
                'yz_threshold_mat', ...
                'xz_threshold_mat', ...
                'threshold_mat_3D', ...
              'counter_add_edge'  , ...
              'counter_add_vertex', ...
              'counter_toggle'    , ...
              'counter_threshold' , ...
           'is_intensity_inverted', ...
                 'dimension_order', ...
      'projection_dimension'      , ...
             'number_of_pixels_2D'  };
          
save_variables = [ autosave_variables{ : }, ...
              {          'backup_ordinate', ...
               'max_number_of_added_edges', ...
                            'elapsed_time', }];

    %% Templates and Images                     

% precalculating the volumes of influence of the vertices as linear indices into the 3D image for
% quick ease of use later during live performance.

number_of_scales = length( lumen_radius_in_microns_range );

scale_subscript_range = 1 : number_of_scales ;

vertex_element_linear_indexing_templates = cell( number_of_scales, 1 );
  edge_element_linear_indexing_templates = cell( number_of_scales, 1 );

vertex_radii_in_pixels_range =  length_dilation_ratio_vertices ...
                             .* lumen_radius_in_microns_range  ...
                             ./ microns_per_pixel              ;
                  
  edge_radii_in_pixels_range =  length_dilation_ratio_edges   ...
                             .* lumen_radius_in_microns_range ...
                             ./ microns_per_pixel             ;
                  
for scale_subscript = scale_subscript_range
      
    % find all pixel locations within the ellipsoid radii from the element position
    vertex_element_linear_indexing_templates{ scale_subscript }                                                        ...
        = int64( construct_structuring_element( vertex_radii_in_pixels_range( scale_subscript, : ), size_of_image ));
     
      edge_element_linear_indexing_templates{ scale_subscript }                                                        ...
        = int64( construct_structuring_element(   edge_radii_in_pixels_range( scale_subscript, : ), size_of_image ));
     
end % constructing edge and vertex element templates FOR scale

% initialize_structuring_elements

%% User interfaces                          

% main display controls

% positions normalized to lock them in place upon re-scaling WDS oct. '18
save_button           = uicontrol( display_figure, 'Style'   , 'pushbutton',                    ...
                                                   'Units'   , 'normalized',                    ...
                                                   'Position', [ 0.045, 0.075, 0.125, 0.05],    ...
                                                   'String'  , 'Save',                          ...
                                                   'FontName', 'Times New Roman',               ...
                                                   'FontUnits', 'normalized',                   ...
                                                   'FontSize', 0.5,                             ...                                                   
                                                   'Callback', @store_curation                  );
                                                        
load_button           = uicontrol( display_figure, 'Style'   , 'pushbutton',                    ...
                                                   'Units'   , 'normalized',                    ...
                                                   'Position', [ 0.045, 0.0125, 0.125, 0.05],   ...
                                                   'String'  , 'Load',                          ...
                                                   'FontName', 'Times New Roman',               ...
                                                   'FontUnits', 'normalized',                   ...
                                                   'FontSize', 0.5,                             ...                                                                                                      
                                                   'Callback', @load_curation                   );
                                               
undo_button           = uicontrol( display_figure, 'Style'   , 'pushbutton',                    ...
                                                   'String'  , 'Undo',                          ...
                                                   'Units'   , 'normalized',                    ...
                                                   'FontName', 'Times New Roman',               ...
                                                   'FontUnits', 'normalized',                   ...
                                                   'FontSize', 0.5,                             ...                                                                                                      
                                                   'Position', [ 0.175, 0.075, 0.125, 0.05],    ...
                                                   'Callback', @undo_curation                   );
                                               
redo_button           = uicontrol( display_figure, 'Style'   , 'pushbutton',                    ...
                                                   'String'  , 'Redo',                          ...
                                                   'FontName', 'Times New Roman',               ...
                                                   'FontUnits', 'normalized',                   ...
                                                   'FontSize', 0.5,                             ...
                                                   'Units'   , 'normalized',                    ...
                                                   'Position', [ 0.175, 0.0125, 0.125, 0.05],   ...
                                                   'Callback', @redo_curation                   );      
                                                                                          
depth_of_view_label         = uicontrol( display_figure, 'Style'   , 'pushbutton',                          ...
                                                   'Units'   , 'normalized',                    ...
                                                   'Position', [ 0.4125, 0.105, 0.165, 0.03 ],  ...
                                                   'String'  , 'Z-Depth',                       ...
                                                   'FontName', 'Times New Roman',               ...
                                                   'FontUnits', 'normalized',                   ...
                                                   'FontSize', 0.5,                             ...
                                                   'Callback', @projection_dimension_callback   );                                               
                                               
depth_of_view_slider        = uicontrol( display_figure, 'Style'   , 'slider',                        ...
                                                   'Min'     , 1,                               ...
                                                   'Max'     , size_of_image( 3 ),              ...
                                                 'SliderStep', [        1   / size_of_image( 3 ),  ...
                                                 ( 2 * view_thickness + 1 ) / size_of_image( 3 ) ],... % slider moves the distance of the slice thickness ...
                                                   'Value'   , 1,                               ...
                                                   'Units'   , 'normalized',                    ...
                                                   'Position', [ 0.305, 0.075, 0.385, 0.04],    ...
                                                   'Callback', @change_depth_of_view                  ); 

view_thickness_slider    = uicontrol( display_figure, 'Style'   , 'slider',                        ...
                                                   'Min'     , 0,                               ...
                                                   'Max'     ,   log((      size_of_image( dimension_order( 3 )) - 1 ) * 3 / 2 + 1 ),  ...
                                                 'SliderStep', [ log( 1 / ( size_of_image( dimension_order( 3 )) - 1 ) * 3 / 2 + 1 ),  ...
                                                                 log( 5 / ( size_of_image( dimension_order( 3 )) - 1 ) * 3 / 2 + 1 ) ],...                                                   
                                                   'Value'   , 0,                               ...
                                                   'Units'   , 'normalized',                    ...
                                                   'Position', [ 0.305, 0.0325, 0.385, 0.03],   ...
                                                   'Callback', @change_view_thickness              );

view_thickness_label     = uicontrol( display_figure, 'Style'   , 'pushbutton',                          ...
                                                   'Units'   , 'normalized',                    ...
                                                   'Position', [ 0.4125, 0, 0.165, 0.03 ],      ...
                                                   'String'  , 'Z-Thickness',                   ...
                                                   'FontName', 'Times New Roman',               ...
                                                   'FontUnits', 'normalized',                   ...
                                                   'FontSize', 0.5,                             ...
                                                   'Callback', @projection_dimension_callback   );                                           
					                                     
paint_button          = uicontrol( display_figure, 'Style'   , 'pushbutton',                    ...
                                                   'Units'   , 'normalized',                    ...
                                                   'Position', [ 0.695, 0.075, 0.125, 0.05],    ...
                                                   'String'  , 'Paint',                         ...
                                                   'FontName', 'Times New Roman',               ...
                                                   'FontUnits', 'normalized',                   ...
                                                   'FontSize', 0.5,                             ...                                                   
                                                   'Callback', @paint_callback                  );                                               
                                               
sweep_button          = uicontrol( display_figure, 'Style'   , 'pushbutton',                    ...
                                                   'Units'   , 'normalized',                    ...
                                                   'Position', [ 0.695, 0.0125, 0.125, 0.05],   ...
                                                   'String'  , 'Sweep',                         ...
                                                   'FontName', 'Times New Roman',               ...
                                                   'FontUnits', 'normalized',                   ...
                                                   'FontSize', 0.5,                             ...                                                                                                      
                                                   'Callback', @sweep_callback                  );
                                               
add_button            = uicontrol( display_figure, 'Style'   , 'pushbutton',                    ...
                                                   'Units'   , 'normalized',                    ...
                                                   'Position', [ 0.825, 0.09166666, 0.125, 0.0333333],...
                                                   'String'  , 'Add Edge',                      ...
                                                   'FontName', 'Times New Roman',               ...
                                                   'FontUnits', 'normalized',                   ...
                                                   'FontSize', 0.5,                             ...                                                                                                      
                                                   'Callback', @add_callback                    );
                                               
add_vertex_button            = uicontrol( display_figure, 'Style'   , 'pushbutton',             ...
                                                   'Units'   , 'normalized',                    ...
                                                   'Position', [ 0.825, 0.05208333, 0.125, 0.0333333],...
                                                   'String'  , 'Add Vertex',                    ...
                                                   'FontName', 'Times New Roman',               ...
                                                   'FontUnits', 'normalized',                   ...
                                                   'FontSize', 0.5,                             ...                                                                                                      
                                                   'Callback', @add_vertex_callback             );

toggle_button         = uicontrol( display_figure, 'Style'   , 'pushbutton',                    ...
                                                   'String'  , 'Toggle',                        ...
                                                   'Units'   , 'normalized',                    ...
                                                   'FontName', 'Times New Roman',               ...
                                                   'FontUnits', 'normalized',                   ...
                                                   'FontSize', 0.5,                             ...                                                                                                      
                                                   'Position', [ 0.825, 0.0125, 0.125, 0.0333333],...
                                                   'Callback', @toggle_callback                 );  
                                                                                                                                                                                                                                                         
%% intensity histogram controls                                               
intensity_min_text_box = uicontrol( intensity_histo_figure,                                     ...
                                                    'Style'   , 'edit',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.1, 0.1, 0.3, 0.075 ],       ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 0.5,                            ...                                                   
                                                    'Callback', @change_intensity_min,          ...
                                                    'String'  , num2str( intensity_limits( 1 )));
                                                        
intensity_max_text_box = uicontrol( intensity_histo_figure,                                     ...
                                                    'Style'   , 'edit',                         ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 0.5,                            ...                                                   
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.6, 0.1, 0.3, 0.075 ],       ...    
                                                    'Callback', @change_intensity_max,          ...
                                                    'String'  , num2str( intensity_limits( 2 )));
                                                        
intensity_min_label    = uicontrol( intensity_histo_figure,                                     ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.1, 0.05, 0.3, 0.05 ],       ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 0.5,                            ...                                                   
                                                    'String'  , 'Min'                           );
                                                        
intensity_max_label    = uicontrol( intensity_histo_figure,                                     ...
                                                    'Style'   , 'text',                         ...
                                                    'Position', [ 210, 40,  80, 20 ],           ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.6, 0.05, 0.3, 0.05 ],       ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 0.5,                            ...                                                    
                                                    'String'  , 'Max'                           );
                                                

invert_button           = uicontrol( intensity_histo_figure,                                      ...
                                                    'Style'   , 'pushbutton',                    ...
                                                    'Units'   , 'normalized',                    ...
                                                    'Position', [ 0.4, 0.93, 0.2, 0.05 ],        ...
                                                    'String'  , 'Inverted',                      ...
                                                    'FontName', 'Times New Roman',               ...
                                                    'FontUnits', 'normalized',                   ...
                                                    'FontSize', 0.5,                             ...                                                                                                      
                                                    'Callback', @invert_callback               );
                                                        
%% energy histogram controls
energy_min_text_box       = uicontrol( energy_histo_figure,                                     ...
                                                    'Style'   , 'edit',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.1, 0.1, 0.2, 0.075 ],       ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 0.5,                            ...
                                                    'Callback', @change_energy_min,             ...
                                                    'String'  , num2str( energy_limits( 1 ))    );
                                                        
energy_max_text_box       = uicontrol( energy_histo_figure,                                     ...
                                                    'Style'   , 'edit',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.7, 0.1, 0.2, 0.075 ],       ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 0.5,                            ...                                                   
                                                    'Callback', @change_energy_max,             ...
                                                    'String'  , num2str( energy_limits( 2 ))    );

energy_threshold_text_box = uicontrol( energy_histo_figure,                                     ...
                                                    'Style'   , 'edit',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.4, 0.1, 0.2, 0.075 ],       ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 0.5,                            ...                                                   
                                                    'String'  , '0',                            ...
                                                    'Callback', @apply_energy_threshold         );
                                                        
energy_min_label          = uicontrol( energy_histo_figure, ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.1, 0.05, 0.2, 0.05 ],       ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 0.5,                            ...                                                    
                                                    'String'  , 'Min'                           );
                                                        
energy_max_label          = uicontrol( energy_histo_figure,                                     ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.7, 0.05, 0.2, 0.05 ],       ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 0.5,                            ...                                                    
                                                    'String'  , 'Max'                           );

                                                           
energy_threshold_label    = uicontrol( energy_histo_figure,                                     ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.4, 0.05, 0.2, 0.05 ],       ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 0.5,                            ...                                                    
                                                    'String'  , 'Threshold'                     );
                                                
true_label                = uicontrol( energy_histo_figure,                                     ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.1, 0.93, 0.2, 0.05 ],       ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 0, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , 'True'                          );

binarize_button           = uicontrol( energy_histo_figure, ...
                                                    'Style'   , 'pushbutton',                    ...
                                                    'Units'   , 'normalized',                    ...
                                                    'Position', [ 0.4, 0.93, 0.2, 0.05 ],        ...
                                                    'String'  , 'Binary',                        ...
                                                    'FontName', 'Times New Roman',               ...
                                                    'FontUnits', 'normalized',                   ...
                                                    'FontSize', 0.5,                             ...                                                                                                      
                                                    'Callback', @binarize_callback               );
                                                                                                                                                
false_label               = uicontrol( energy_histo_figure,                                     ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.7, 0.93, 0.2, 0.05 ],       ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 0, 0 ],                    ...                                                                                                        
                                                    'String'  , 'False'                         );
                                                
%% Minimap controls
micron_label              = uicontrol( minimap_figure,                                          ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.8, 0.93, 0.18, 0.05 ],       ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ .3, .2, 1 ],                  ...  
                                             'ForegroundColor', [ 1, 1, 0.1 ],                  ...                                                                                                                                                     
                                                    'String'  , 'Microns'                       );

voxel_label               = uicontrol( minimap_figure,                                          ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.01, 0.14, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , 'Voxels'                        );                                                

x_label                  = uicontrol( minimap_figure,                                           ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.11, 0.14, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , 'X'                             );

y_label                  = uicontrol( minimap_figure,                                           ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.21, 0.14, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , 'Y'                             );



z_label                  = uicontrol( minimap_figure,                                           ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.31, 0.14, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , 'Z'                             );
                  
start_label               = uicontrol( minimap_figure,                                          ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.01, 0.10, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , 'Start'                         );                                                
                                                
x_start_label             = uicontrol( minimap_figure,                                          ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.11, 0.10, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , ''                              );

y_start_label            = uicontrol( minimap_figure,                                           ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.21, 0.10, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , ''                              );



z_start_label            = uicontrol( minimap_figure,                                           ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.31, 0.10, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , ''                              );
                                                
end_label                = uicontrol( minimap_figure,                                           ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.01, 0.06, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , 'End'                           );                                                
                                                
                                                
x_end_label             = uicontrol( minimap_figure,                                            ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.11, 0.06, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , ''                              );

y_end_label            = uicontrol( minimap_figure,                                             ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.21, 0.06, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , ''                              );



z_end_label            = uicontrol( minimap_figure,                                             ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.31, 0.06, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , ''                              );
                                             
total_label            = uicontrol( minimap_figure,                                             ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.01, 0.02, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , 'Total'                         );                                                
                                                
                                                
x_total_label            = uicontrol( minimap_figure,                                           ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.11, 0.02, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , ''                              );

y_total_label           = uicontrol( minimap_figure,                                            ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.21, 0.02, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , ''                              );



z_total_label           = uicontrol( minimap_figure,                                            ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.31, 0.02, 0.09, 0.03 ],     ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 1, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , ''                              );
                                             
%% initialize displays
load_colormap

if load_upon_startup, load_curation, else initialize_structuring_elements, end %#ok<SEPEX>

paint_index_image, update_edge_intensities, update_z_zoom
update_z_sliders, update_minimap, update_xy_zoom
                     
% turn all figures visible     

        display_figure.Visible = 'on';
        minimap_figure.Visible = 'on';
intensity_histo_figure.Visible = 'on';
   energy_histo_figure.Visible = 'on';
    
    %% Visual Updating Functions        
    
    %% initialize structuring element               
    function initialize_structuring_elements

        % precalculate linear indices from space subscripts in the image
%         degrees_of_edges_uint_32 = uint32([ ones( max_number_of_added_edges + 1, 1 ); degrees_of_edges ]);
%         degrees_of_edges_uint_32 = uint32( degrees_of_edges );
        degrees_of_edges_uint_32 = uint32( cellfun( @length, edge_scale_subscripts ));

        edge_space_subscripts_mat =         cell2mat( edge_space_subscripts );
        edge_scale_subscripts_mat = uint8(  cell2mat( edge_scale_subscripts ));        

        space_subscripts_int64 = int64( edge_space_subscripts_mat );

        position_linear_indices =   space_subscripts_int64( :, 1 )                                               ...
                                + ( space_subscripts_int64( :, 2 ) - 1 ) * size_of_image( 1 )                    ...
                                + ( space_subscripts_int64( :, 3 ) - 1 ) * size_of_image( 1 ) * size_of_image( 2 );

        sphere_position_linear_indices_cell = num2cell( position_linear_indices );
        
        edge_positions_linear_indices        = mat2cell( position_linear_indices, degrees_of_edges_uint_32 );

        structuring_element_linear_indexing = edge_element_linear_indexing_templates( edge_scale_subscripts_mat );

        sphere_structure_positions_linear_indexing = cellfun( @plus, structuring_element_linear_indexing, ...
                                                                     sphere_position_linear_indices_cell, ...
                                                                                  'UniformOutput', false  ); 

%         vertex_structure_positions_linear_indexing = structuring_element_linear_indexing_templates( edge_scale_subscripts_mat );

        % !!!!!!!!!!!!! mat2cell may have bug for type non-double second inputs
        sphere_structure_positions_linear_indexing_edge_cell = mat2cell( sphere_structure_positions_linear_indexing,  degrees_of_edges_uint_32, 1 );

%         end_vertices_structure_positions_linear_indexing_1 = cellfun( @( x ) x{  1  }, sphere_structure_positions_linear_indexing_edge_cell, 'UniformOutput', false );
%         end_vertices_structure_positions_linear_indexing_2 = cellfun( @( x ) x{ end }, sphere_structure_positions_linear_indexing_edge_cell, 'UniformOutput', false );

        space_subscripts_int64 = int64( vertex_space_subscripts );

        position_linear_indices =   space_subscripts_int64( :, 1 )                                               ...
                                + ( space_subscripts_int64( :, 2 ) - 1 ) * size_of_image( 1 )                    ...
                                + ( space_subscripts_int64( :, 3 ) - 1 ) * size_of_image( 1 ) * size_of_image( 2 );

        vertex_position_linear_indices_cell = num2cell( position_linear_indices );

        structuring_element_linear_indexing = vertex_element_linear_indexing_templates( round( vertex_scale_subscripts ));

        vertex_structure_positions_linear_indexing = cellfun( @plus, structuring_element_linear_indexing,   ...
                                                                       vertex_position_linear_indices_cell, ...
                                                                                     'UniformOutput', false ); 

%         vertex_structure_positions_linear_indexing = [ end_vertices_structure_positions_linear_indexing_1, ...
%                                                        end_vertices_structure_positions_linear_indexing_2  ];

        edge_structure_positions_linear_indexing = cellfun( @cell2mat, sphere_structure_positions_linear_indexing_edge_cell, 'UniformOutput', false );

        edge_structure_positions_linear_indexing = cellfun( @( x ) unique( x, 'rows' ), edge_structure_positions_linear_indexing, 'UniformOutput', false );

        edge_structure_positions_subscript_xy = cellfun( @( x )      1 + mod( x - 1,  prod( size_of_image([ 1, 2 ]))), edge_structure_positions_linear_indexing, 'UniformOutput', false );
        
        edge_structure_positions_subscript_z  = cellfun( @( x ) floor( double( x ) ./ prod( size_of_image([ 1, 2 ]))), edge_structure_positions_linear_indexing, 'UniformOutput', false );
        
        edge_structure_positions_subscript_y  = cellfun( @( x )      1 + mod( x - 1,        size_of_image( 1 )),       edge_structure_positions_subscript_xy,    'UniformOutput', false );
        
        edge_structure_positions_subscript_x  = cellfun( @( x ) floor( double( x ) ./       size_of_image( 1 )),       edge_structure_positions_subscript_xy,    'UniformOutput', false );        
                
        edge_structure_positions_subscripts = cellfun( @( y, x, z ) [ y, x, z ], edge_structure_positions_subscript_y,                        ...
                                                                                 edge_structure_positions_subscript_x,                        ...
                                                                                 edge_structure_positions_subscript_z, 'UniformOutput', false );
                                            
        edge_centers = cellfun( @( x ) round( mean( x, 1 )), edge_structure_positions_subscripts, 'UniformOutput', false );

        edge_centers = cell2mat( edge_centers );

    end
        %% overlay                                  
    function update_overlay

    %     % Make all truth images render at same intensity
    %     %note: this does not cancel background vector which should help
    %     %false positives shine though
    %     truth_image_uint8 = uint8(logical(truth_image_uint8).*128);

        % combine the raw and truth images into overlay
        overlay_image = truth_image_uint8 + adjusted_original_image_2D ;

        % draw the image in the display axes
        image( display_axes, overlay_image )

        % reset axes properties and zoom position
        display_axes.XTick      = [ ];
        display_axes.YTick      = [ ];

        display_axes.XTickLabel = '' ;
        display_axes.YTickLabel = '' ;

        display_axes.XLim = FOV_limits( dimension_order( 2 ), : );
        display_axes.YLim = FOV_limits( dimension_order( 1 ), : );

    %     display_axes.ButtonDownFcn = @vertex_toggler2 ;

    end % update_overlay
        %% vector                                   
    function update_index_image_crop

        % create 2D index image from 3D index image by taking first non-background (background has index
        % 1) entry in third dimension note: permute matrix before to project in different dimension
        
        % % for a y-z projection
        % permute( index_image, [ 3, 1, 2 ])
        % 
        % % for a x-z projection
        % permute( index_image, [ 2, 3, 1 ])
        
              index_image_crop = permute(       index_image,        dimension_order );
              index_image_crop =                index_image_crop( :, :, depth_range );
        
        centerlines_image_crop = permute( centerlines_image,        dimension_order );
        centerlines_image_crop =          centerlines_image_crop( :, :, depth_range );
        
        update_index_image_2D_crop

    end % update_index_image_2D_crop

    function update_index_image_2D_crop
        
        centerlines_image_2D = max( centerlines_image_crop, [ ], 3 ); % opaque ray tracing in depth
        
        for index = 1 : 2 % 1. vertices 2. edges
        
            object_index_image_crop = index_image_crop ;

            if index == 1, object_index_image_crop = - object_index_image_crop ; end    

            object_index_image_crop( object_index_image_crop < 0 ) = 0 ; % remove 1. edges xor 2. vertices
            
%             object_index_image_2D = max( object_index_image_crop ~= 0, [ ], 3 ); % opaque ray tracing to first object in depth dimension
            
            [ ~ , argmax_vertex_image_2D ] = max( object_index_image_crop ~= 0, [ ], 3 ); % opaque ray tracing in depth
            
            linear_index_of_max_image_2D = zeros( size( argmax_vertex_image_2D ));

            linear_index_of_max_image_2D( : ) = ( argmax_vertex_image_2D( : ) - 1 ) * number_of_pixels_2D + ( 1 : number_of_pixels_2D )';
            
            if index == 1 % vertices

                vertex_index_image_2D = object_index_image_crop( linear_index_of_max_image_2D );

            else % edges

                       index_image_2D = object_index_image_crop( linear_index_of_max_image_2D ) ;
                       index_image_2D( index_image_2D  < 0 ) = 0 ; % remove vertices
                       index_image_2D( index_image_2D == 0 ) = 1 ; % set background from 0 to 1

            end    
        end % FOR vertices and edges
            
%         [ ~, argmax_vertex_image_2D ] = max( index_image_crop ~= last_vertex_index + 1 , [ ], 3 );
%         [ ~ , argmax_vertex_image_2D ] = max( index_image_crop ~= 0 || index_image_crop ~= 1, [ ], 3 );
%         [ ~ , argmax_vertex_image_2D ] = max( object_index_image_crop ~= 0, [ ], 3 );

%         [ ~ , argmax_vertex_image_2D ] = min( index_image_crop ~= 0, [ ], 3 );
%         argmax_vertex_image_2D = max(argmax_vertex_image_2D, argmax_index_image_2D);

        % number_of_pixels_2D = numel( argmax );

%         linear_index_of_max_image_2D = zeros( size( argmax_vertex_image_2D ));
% 
%         linear_index_of_max_image_2D( : ) = ( argmax_vertex_image_2D( : ) - 1 ) * number_of_pixels_2D + ( 1 : number_of_pixels_2D )';

        % index2intensity has the first index dedicated to background and no indices dedicated to
        % vertices (most indies are dedicated to edges).  Subtracting off vertices and background
        % (unsigned integer type will round those to zero), then adding one so that the first real edge
        % points to the second entry in index2intensity.                        
%         vertex_index_image_2D = index_image_crop( linear_index_of_max_image_2D ); 
%         vertex_index_image_2D = - index_image_crop( linear_index_of_max_image_2D ); 

        % background coded as zero for vertex index image
%         vertex_index_image_2D( vertex_index_image_2D >= last_vertex_index + 1 ) = 0 ;
%         vertex_index_image_2D( vertex_index_image_2D < 0 ) = 0 ;
                
%         vertex_index_image_2D = vertex_index_image_2D + 1 ;

%         centerlines_image_2D = max( centerlines_image_crop, [ ], 3 );
                
% %         index_image_2D = index_image_crop( linear_index_of_max_image_2D ) - ( last_vertex_index + 1 ) + 1 ;
% %        index_image_2D = index_image_crop ;
%         index_image_2D = index_image_crop( linear_index_of_max_image_2D ) ;
%         
%         index_image_2D( index_image_2D < 0 ) = 0 ;
%         
% %        index_image_2D = index_image_2D + 1 ;
%         index_image_2D(index_image_2D == 0) = 1 ;
                
        update_edge_intensities

    end % update_index_image

    function update_edge_intensities

%         % !!!! change this to an update edge color and make the color code the direction, or assign
%         % random colors. expand IF to SWITCH with at least 4 cases !!!!!!!!!!!        
%         if ( are_intensity_limits_binarized  )
% 
%             index2intensity = 128 * ones( size( mean_edge_energies ));
%             
%         else
%             
%             index2intensity = 128 * min( 1, (     energy_limits( 2 ) - mean_edge_energies )   ...
%                                               / ( energy_limits( 2 ) - energy_limits( 1 )));
% 
%         end % IF binarized intensity
        
        switch edge_display_method
            
            case 1 % binary
                
                index2intensity = 128 * ones( size( mean_edge_energies ));
                
            case 2 % energy weighted
                
                index2intensity = 128 * min( 1, (     energy_limits( 2 ) - mean_edge_energies )   ...
                                                  / ( energy_limits( 2 ) - energy_limits( 1 )));
                
            case 3 % random unique
                
                index2intensity = 128 * randperm( numel( mean_edge_energies )) / numel( mean_edge_energies );
                
            case 4
                
                
        end

        update_truth_image

    end % update_truth_image_intensity

    function update_truth_image
        
        truth_image_uint8 =  uint8( index2intensity( index_image_2D )                     ...
                                    .* double( displayed_edges( index_image_2D ) & ~ deleted_edges( index_image_2D ))       ...
                                    .* cat( 3,         double( ~ true_edges( index_image_2D ))  ...
                                               + 0.3 * double(   true_edges( index_image_2D )), ...
                                                 0.3 * double( ~ true_edges( index_image_2D ))  ...
                                               + 0.7 * double(   true_edges( index_image_2D )), ...
                                                 0.3 * double( ~ true_edges( index_image_2D ))  ...
                                               +       double(   true_edges( index_image_2D ))) ...
                                    .* ~ centerlines_image_2D                                  ) ... % show edge wherever there is no centerline
                          +  uint8( 128 *     vertex_color .* logical( vertex_index_image_2D ) .* ( index_image_2D == 1 ) .* ~ centerlines_image_2D )      ... % show vertex wherever there is no edge or centerline
                          +  uint8( 128 * centerline_color .*           centerlines_image_2D  );
                                                            

    end % update_truth_image

%     function update_centerlines_image
%         
%         centerlines_image = 
%         
%     end % update_centerlines_image
        %% raw                                      
    function update_original_image

        % crop 3D volume
        original_image_crop = permute( original_image, dimension_order );

        original_image_crop = original_image_crop( :, :, depth_range );

        % take a cube of the 3D image, max project in the third dimension and if requested, 
        [ original_image_crop_2D, argmax_original_image_2D ] = max( original_image_crop, [ ], 3 );

    %     if gaussian_depth

%         % do Guassian weighting in the third dimension based on the z location of the argument of
%         % the max projection
%         minimum_intensity_in_crop = double( min( original_image_crop_2D( : )));
% 
%         original_image_crop_2D = double(   original_image_crop_2D                  ...
%                                          - minimum_intensity_in_crop )             ...
%                                .* gaussian_weight_in_z( argmax_original_image_2D ) ...
%                                + minimum_intensity_in_crop ;

        original_image_crop_2D = double( original_image_crop_2D );

    %     end % IF gaussian_depth

        update_original_image_intensity

    end % update_original_image

    function update_original_image_intensity

        % adjust (and possibly invert the contrast for the 2D original image and restrict intensity
        % values to the integers 0 : 128
        adjusted_original_image_2D                                                 ...
                    = uint8( is_intensity_inverted * 128           ...
                 + ( - 1 ) ^ is_intensity_inverted                 ...
                 * min( 1,   ( original_image_crop_2D - intensity_limits( 1 )) ...
                           / ( intensity_limits( 2 )  - intensity_limits( 1 )))...
                 * 128                                                         );  
        intensity_histo_axes.XLim = intensity_limits([ 1, 2 ]);                 

    end % update_original_image_contrast
        %% histogram                                
    function update_intesnity_histogram

        y_range =   max(               1                     , round( FOV_limits( dimension_order( 1 ), 1 ))) ...
                  : min( size_of_image( dimension_order( 1 )), round( FOV_limits( dimension_order( 1 ), 2 ))) ;

        x_range =   max(               1                     , round( FOV_limits( dimension_order( 2 ), 1 ))) ...
                  : min( size_of_image( dimension_order( 2 )), round( FOV_limits( dimension_order( 2 ), 2 ))) ;
              
        histogram( original_image_crop_2D( y_range, x_range ), 'Parent', intensity_histo_axes, 'FaceColor', [ 0, 0, 0 ], 'EdgeColor', 'none' );

        try
            intensity_histo_axes.XLim = intensity_limits ;    
        catch
            warning('The intensity min should be below the max');
        end

        xlabel(intensity_histo_axes, 'Pixel Intensity')

        zoom( intensity_histo_axes, 'yon' )

    end % update_intesnity_histogram   

    function update_energy_histogram

        histogram( mean_edge_energies( in_view_edges & displayed_edges & ~ deleted_edges ), 'Parent', energy_histo_axes, 'FaceColor', [0, 0, 0], 'EdgeColor', 'white' );
        
        colormap( energy_histo_axes, colormap_mat );
        
        %COMMENT OUT TO SHOW ALL THRESHOLDS ON HISTOGRAM
        adj_y_lim = 1 + [round((size_of_minimap-1)*FOV_limits( 1, 1)/size_of_image(1)) round((size_of_minimap-1)*FOV_limits( 1, 2)/size_of_image(1))];    
        adj_x_lim = 1 + [round((size_of_minimap-1)*FOV_limits( 2, 1)/size_of_image(2)) round((size_of_minimap-1)*FOV_limits( 2, 2)/size_of_image(2))];
        adj_z_lim = 1 + [round((size_of_minimap-1)*FOV_limits( 3, 1)/size_of_image(3)) round((size_of_minimap-1)*FOV_limits( 3, 2)/size_of_image(3))];

        adj_z_lim( 1 ) = max( adj_z_lim( 1 ) - 1, 1   );
        adj_z_lim( 2 ) = min( adj_z_lim( 2 ) + 1, size_of_minimap );
        
%         unique_thresholds =   intersect(intersect(xy_threshold_mat(adj_y_lim(1):adj_y_lim(2), adj_x_lim(1):adj_x_lim(2)),  ...
%                                                   xz_threshold_mat(adj_z_lim(1):adj_z_lim(2), adj_x_lim(1):adj_x_lim(2))), ...
%                                                   yz_threshold_mat(adj_z_lim(1):adj_z_lim(2), adj_y_lim(1):adj_y_lim(2)))  ;
        unique_thresholds = unique(threshold_mat_3D(adj_y_lim(1):adj_y_lim(2), ...
                                                    adj_x_lim(1):adj_x_lim(2), ...
                                                    adj_z_lim(1):adj_z_lim(2)));
        %END COMMENT OUT 
                        
        %UNCOMMENT TO SHOW ALL THRESHOLDS
%        unique_thresholds = [xy_threshold_mat, xz_threshold_mat, yz_threshold_mat];
%        unique_thresholds = unique(unique_thresholds(:,:));
        %END UNCOMMENT

        unique_thresholds = sort(unique_thresholds);
        unique_thresholds(unique_thresholds == 0) = NaN;
        
        color_factors = unique_thresholds;
        color_factors( color_factors == 0 ) = NaN;
        threshold_mat_3D(threshold_mat_3D == 0) = NaN;
        color_factors = color_factors - min(threshold_mat_3D(:));
        if (max(threshold_mat_3D(:)-min(threshold_mat_3D(:)))) ~= 0
            color_factors = color_factors ./ (max(threshold_mat_3D(:)-min(threshold_mat_3D(:))));
        end
        color_factors = 1 - color_factors;
        color_factors = color_factors .* (threshold_color_range(2) - threshold_color_range(1)) + threshold_color_range(1);
        color_factors( isnan( color_factors ) ) = out_of_bounds_color ;
        threshold_mat_3D(isnan(threshold_mat_3D)) = 0;

        color_factors( isnan(color_factors) ) = out_of_bounds_color;    
        unique_thresholds( isnan(unique_thresholds)) = 0;
        
        if min(unique_thresholds) > energy_limits(1)
            unique_thresholds = [energy_limits(1); unique_thresholds];
            color_factors = [out_of_bounds_color; color_factors];
        end
        if max(unique_thresholds) < energy_limits(2)
            unique_thresholds = [unique_thresholds; energy_limits(2)];
            color_factors = [color_factors; color_factors(end)];
        end
        
        color_factors = 1 + round(255 .* color_factors);
        for threshold_index = 1:length(unique_thresholds)-1  
            try rectangle( energy_histo_axes, 'Position', [ unique_thresholds(threshold_index),   0, unique_thresholds(threshold_index+1) - unique_thresholds(threshold_index), energy_histo_axes.YLim( 2 ) ], 'FaceColor', [colormap_mat(color_factors(threshold_index),1) colormap_mat(color_factors(threshold_index),2) colormap_mat(color_factors(threshold_index),3) 0.66], 'EdgeColor', 'none'), end
        end
        
%         try rectangle( energy_histo_axes, 'Position', [ energy_threshold,   0, energy_limits( 2 ) - energy_threshold, energy_histo_axes.YLim( 2 ) ], 'FaceColor', [1 0 0 0.5]), end

%         try rectangle( energy_histo_axes, 'Position', [ energy_limits( 1 ), 0, energy_threshold - energy_limits( 1 ), energy_histo_axes.YLim( 2 ) ], 'FaceColor', [0 1 1 0.5]), end

        try energy_histo_axes.XLim = energy_limits ; end

%        histogram( vertex_energies( in_view_vertices & displayed_vertices ), 'Parent', energy_histo_axes, 'FaceColor', [0, 0, 0], 'EdgeColor', 'white' );
%        energy_histo_axes.Color = [1 1 1 0];
%        h.Color = 'none';
        xlabel(energy_histo_axes, 'Vertex Energy')    

        zoom( energy_histo_axes, 'yon' )

    end % update_energy_histogram

%     function update_energy_histogram
% 
%         histogram( mean_edge_energies( in_view_edges & displayed_edges ), 'Parent', energy_histo_axes, 'FaceColor', [ 0, 0, 0 ], 'EdgeColor', 'none' );
% 
%         try rectangle( energy_histo_axes, 'Position', [ energy_threshold,   0, energy_limits( 2 ) - energy_threshold, energy_histo_axes.YLim( 2 ) ], 'FaceColor', [1 0 0 0.5]), end
% 
%         try rectangle( energy_histo_axes, 'Position', [ energy_limits( 1 ), 0, energy_threshold - energy_limits( 1 ), energy_histo_axes.YLim( 2 ) ], 'FaceColor', [0 1 1 0.5]), end
% 
%         try energy_histo_axes.XLim = energy_limits ; end
% 
%         xlabel(energy_histo_axes, 'Vertex Energy')    
% 
%         zoom( energy_histo_axes, 'yon' )
% 
%     end % update_energy_histogram

        %% mini-map visual updates                  

    function update_threshold_visualization_matricies

        adj_y_lim = 1 + [round(99*FOV_limits( 1, 1)/size_of_image(1)) round(99*FOV_limits( 1, 2)/size_of_image(1))];    
        adj_x_lim = 1 + [round(99*FOV_limits( 2, 1)/size_of_image(2)) round(99*FOV_limits( 2, 2)/size_of_image(2))];
        adj_z_lim = 1 + [round(99*FOV_limits( 3, 1)/size_of_image(3)) round(99*FOV_limits( 3, 2)/size_of_image(3))];

        adj_z_lim( 1 ) = max( adj_z_lim( 1 ) - 1, 1   );
        adj_z_lim( 2 ) = min( adj_z_lim( 2 ) + 1, size_of_minimap );

        xy_threshold_mat(adj_y_lim(1):adj_y_lim(2), adj_x_lim(1):adj_x_lim(2)) = energy_threshold;
        xz_threshold_mat(adj_z_lim(1):adj_z_lim(2), adj_x_lim(1):adj_x_lim(2)) = energy_threshold;
        yz_threshold_mat(adj_z_lim(1):adj_z_lim(2), adj_y_lim(1):adj_y_lim(2)) = energy_threshold;
        
        threshold_mat_3D(adj_y_lim(1):adj_y_lim(2),                 ...
                         adj_x_lim(1):adj_x_lim(2),                 ...
                         adj_z_lim(1):adj_z_lim(2)) = energy_threshold;

    end %update_minimap_threshold

    function update_minimap

        % speed saver by locking resolution to size_of_minimap computationally
        adj_y_lim = 1 + [round((size_of_minimap-1)*FOV_limits( 1, 1)/size_of_image(1)) round((size_of_minimap-1)*FOV_limits( 1, 2)/size_of_image(1))];    
        adj_x_lim = 1 + [round((size_of_minimap-1)*FOV_limits( 2, 1)/size_of_image(2)) round((size_of_minimap-1)*FOV_limits( 2, 2)/size_of_image(2))];
        adj_z_lim = 1 + [round((size_of_minimap-1)*FOV_limits( 3, 1)/size_of_image(3)) round((size_of_minimap-1)*FOV_limits( 3, 2)/size_of_image(3))];

        adj_z_lim( 1 ) = max( adj_z_lim( 1 ) - 1, 1   );
        adj_z_lim( 2 ) = min( adj_z_lim( 2 ) + 1, size_of_minimap );

        colormap( minimap_figure, colormap_mat );

        in_view_opacity = 0.6;        
        
        %xy
        [x_mini,y_mini] = meshgrid(linspace(0,size_of_image(2),size_of_minimap),linspace(0,-size_of_image(1),size_of_minimap));
    %     z_mini = 0.*x_mini + 0.*y_mini;
        z_mini = zeros( size_of_minimap );
        %EDIT: changed for threshold mapping
        view_mini = NaN .* zeros(size_of_minimap);        
        c_mini = xy_threshold_mat;
        c_mini( c_mini == 0 ) = NaN;
        c_mini = c_mini - min(c_mini(:));
        if max(c_mini(:)) ~= 0
            c_mini = c_mini ./ max(c_mini(:));
        end
        c_mini = 1 - c_mini;
        c_mini = c_mini .* (threshold_color_range(2) - threshold_color_range(1)) + threshold_color_range(1);
        c_mini( isnan( c_mini ) ) = out_of_bounds_color ;
        view_mini(max(1,adj_y_lim(1)):min(size_of_minimap,adj_y_lim(2)), max(1,adj_x_lim(1)):min(size_of_minimap,adj_x_lim(2))) = in_view_color;
        c_mini(1,:) = border_color; c_mini(end,:) = border_color; c_mini(:,1) = border_color; c_mini(:,end) = border_color;
        hold(minimap_axes,'off'); %Makes sure minimap is reset by first surf call of function
        surf(minimap_axes, x_mini.*microns_per_pixel(2), y_mini.*microns_per_pixel(1), z_mini.*microns_per_pixel(3), c_mini,'LineStyle','none'); 
        hold(minimap_axes,'on');
        surf(minimap_axes, x_mini.*microns_per_pixel(2), y_mini.*microns_per_pixel(1), z_mini.*microns_per_pixel(3), view_mini,'LineStyle','none','FaceAlpha',in_view_opacity);
        z_mini = z_mini - size_of_image(3);
        surf(minimap_axes, x_mini.*microns_per_pixel(2), y_mini.*microns_per_pixel(1), z_mini.*microns_per_pixel(3), c_mini,'LineStyle','none');
        surf(minimap_axes, x_mini.*microns_per_pixel(2), y_mini.*microns_per_pixel(1), z_mini.*microns_per_pixel(3), view_mini,'LineStyle','none','FaceAlpha',in_view_opacity);
        
        %yz
        [y_mini,z_mini] = meshgrid(linspace(0,-size_of_image(1),size_of_minimap),linspace(0,-size_of_image(3),size_of_minimap));
    %     x_mini = 0.*y_mini + 0.*z_mini;
        x_mini = zeros( size_of_minimap );
    %    c_mini = out_of_bounds_color.*ones(size_of_minimap);
        %EDIT: changed for threshold mapping
        view_mini = NaN .* zeros(size_of_minimap);        
        c_mini = yz_threshold_mat;
        c_mini( c_mini == 0 ) = NaN;
        c_mini = c_mini - min(c_mini(:));
        if max(c_mini(:)) ~= 0
            c_mini = c_mini ./ max(c_mini(:));
        end
        c_mini = 1 - c_mini;
        c_mini = c_mini .* (threshold_color_range(2) - threshold_color_range(1)) + threshold_color_range(1);
        c_mini( isnan(c_mini)) = out_of_bounds_color ;
        view_mini(max(1,adj_z_lim(1)):min(size_of_minimap,adj_z_lim(2)),max(1,adj_y_lim(1)):min(size_of_minimap,adj_y_lim(2))) = in_view_color;
        c_mini(1,:) = border_color; c_mini(end,:) = border_color; c_mini(:,1) = border_color; c_mini(:,end) = border_color;
        surf(minimap_axes, x_mini.*microns_per_pixel(2), y_mini.*microns_per_pixel(1), z_mini.*microns_per_pixel(3), c_mini,'LineStyle','none','CDataMode','manual');
        surf(minimap_axes, x_mini.*microns_per_pixel(2), y_mini.*microns_per_pixel(1), z_mini.*microns_per_pixel(3), view_mini,'LineStyle','none','CDataMode','manual','FaceAlpha',in_view_opacity);
        x_mini = x_mini + size_of_image(2);
        surf(minimap_axes, x_mini.*microns_per_pixel(2), y_mini.*microns_per_pixel(1), z_mini.*microns_per_pixel(3), c_mini,'LineStyle','none');
        surf(minimap_axes, x_mini.*microns_per_pixel(2), y_mini.*microns_per_pixel(1), z_mini.*microns_per_pixel(3), view_mini,'LineStyle','none','CDataMode','manual','FaceAlpha',in_view_opacity);

        %xz
        [x_mini,z_mini] = meshgrid(linspace(0,size_of_image(2),size_of_minimap),linspace(0,-size_of_image(3),size_of_minimap));
    %     y_mini = 0.*x_mini + 0.*z_mini;
        y_mini = zeros( size_of_minimap );
        %c_mini = out_of_bounds_color.*ones(size_of_minimap);
        %EDIT: changed for threshold mapping
        view_mini = NaN * zeros(size_of_minimap);        
        c_mini = xz_threshold_mat;
        c_mini( c_mini == 0 ) = NaN;
        c_mini = c_mini - min(c_mini(:));
        if max(c_mini(:)) ~= 0
            c_mini = c_mini ./ max(c_mini(:));
        end
        c_mini = 1 - c_mini;
        c_mini = c_mini .* (threshold_color_range(2) - threshold_color_range(1)) + threshold_color_range(1);
        c_mini( isnan(c_mini)) = out_of_bounds_color ;
        view_mini(max(1,adj_z_lim(1)):min(size_of_minimap,adj_z_lim(2)),max(1,adj_x_lim(1)):min(size_of_minimap,adj_x_lim(2))) = in_view_color;        c_mini(1,:) = border_color; c_mini(end,:) = border_color; c_mini(:,1) = border_color; c_mini(:,end) = border_color;
        surf(minimap_axes, x_mini.*microns_per_pixel(2), y_mini.*microns_per_pixel(1), z_mini.*microns_per_pixel(3), c_mini,'LineStyle','none');
        surf(minimap_axes, x_mini.*microns_per_pixel(2), y_mini.*microns_per_pixel(1), z_mini.*microns_per_pixel(3), view_mini,'LineStyle','none','FaceAlpha',in_view_opacity);
        y_mini = y_mini - size_of_image(1);
        surf(minimap_axes, x_mini.*microns_per_pixel(2), y_mini.*microns_per_pixel(1), z_mini.*microns_per_pixel(3), c_mini,'LineStyle','none');
        surf(minimap_axes, x_mini.*microns_per_pixel(2), y_mini.*microns_per_pixel(1), z_mini.*microns_per_pixel(3), view_mini,'LineStyle','none','FaceAlpha',in_view_opacity);

        caxis(minimap_axes,[0 1.0]); %constrained color axis
        axis(minimap_axes,[0 size_of_image(2)*microns_per_pixel(2)   ...
                            -size_of_image(1)*microns_per_pixel(1) 0 ...
                            -size_of_image(3)*microns_per_pixel(3) 0 ]);
        xlabel(minimap_axes,'X (um)');
        ylabel(minimap_axes,'Y (um)');
        zlabel(minimap_axes,'Z (um)');
        daspect(minimap_axes,[1 1 1]); 
        xticks(minimap_axes,0:round( size_of_image(2)*microns_per_pixel(2)/4) :round(size_of_image(2)*microns_per_pixel(2)     ));
        yticks(minimap_axes,  round(-size_of_image(1)*microns_per_pixel(1))   :round(size_of_image(1)*microns_per_pixel(1)/4):0 );
        zticks(minimap_axes, fliplr( 0   :-round(size_of_image(3)*microns_per_pixel(3)/4):round(-size_of_image(3)*microns_per_pixel(3))));

        set( x_start_label, 'String', num2str( round( FOV_limits( 2, 1 ))));
        set( y_start_label, 'String', num2str( round( FOV_limits( 1, 1 ))));
        set( z_start_label, 'String', num2str(   max( FOV_limits( 3, 1 ), 1 )));

        set( x_end_label, 'String', num2str( round( FOV_limits( 2, 2 ))));
        set( y_end_label, 'String', num2str( round( FOV_limits( 1, 2 ))));
        set( z_end_label, 'String', num2str(   min( FOV_limits( 3, 2 ), size_of_image( 3 ))));

        set( x_total_label, 'String', num2str( round( FOV_limits( 2, 2 )) - round( FOV_limits( 2, 1 )) + 1 ));
        set( y_total_label, 'String', num2str( round( FOV_limits( 1, 2 )) - round( FOV_limits( 1, 1 )) + 1 ));
        set( z_total_label, 'String', num2str(   min( FOV_limits( 3, 2 ), size_of_image( 3 )) ...
                                               - max( FOV_limits( 3, 1 ), 1 ) + 1             ));

    end % update_minimap

        %% z navigation                             
    function update_z_zoom
        % update main display

        % compute limits from the z depth and thickness settings
        FOV_limits( dimension_order( 3 ), 1 : 2 ) = [ depth_of_view - view_thickness, ...
                                                      depth_of_view + view_thickness  ];    

        depth_range = max(               1                     , FOV_limits( dimension_order( 3 ), 1 )) ...
                    : min( size_of_image( dimension_order( 3 )), FOV_limits( dimension_order( 3 ), 2 )) ;
            
        FOV_limits( dimension_order( 3 ), 1 : 2 ) = [ depth_range( 1 ), depth_range( end )];
        
%         view_thickness_slider_thickness = max( round(( view_thickness ) / 4 ), 1 );        
        set(  depth_of_view_slider, 'SliderStep', [ 1 /        size_of_image( dimension_order( 3 )), ...
                                  ( 2 * view_thickness + 1 ) / size_of_image( dimension_order( 3 )) ]);
        set( view_thickness_slider, 'SliderStep', [ log( 1 / ( size_of_image( dimension_order( 3 ) ) - 1 ) * 3 / 2 + 1 ), ...
                                                    log( 5 / ( size_of_image( dimension_order( 3 ) ) - 1 ) * 3 / 2 + 1 ) ]);

        set(  depth_of_view_slider, 'Max',       size_of_image( dimension_order( 3 ))) ;
        set( view_thickness_slider, 'Max', log(( size_of_image( dimension_order( 3 )) - 1 ) * 3 / 2 + 1 ));

        % find vertices in the view in z
        in_view_edges_z = edge_centers( :, dimension_order( 3 )) >= depth_range(  1  ) ... 
                        & edge_centers( :, dimension_order( 3 )) <= depth_range( end ) ;

        in_view_edges = in_view_edges_xy & in_view_edges_z ;    

%         % create gaussian for psuedo focus weighting    
%         gaussian_weight_in_z = exp( - ( depth_range - depth_of_view ) .^ 2 / view_thickness ^ 2 / 4 );
% 
%         if view_thickness == 0, gaussian_weight_in_z( 1 ) = 1 ; end

        update_index_image_crop, update_original_image, update_overlay

        update_intesnity_histogram, update_energy_histogram, 

        update_minimap

    end % update z zoom

    function update_z_sliders

        set( depth_of_view_slider,  'Value',      1 ) % to refresh the slider, for some reason it doesn't refresh        
        
        set(  depth_of_view_slider, 'Value',      depth_of_view       )
        set( view_thickness_slider, 'Value', log( view_thickness + 1 ))

    end % update_sliders

    %% Callback Functions
        %% volume navigation                        
    function update_xy_zoom( ~, ~ )
        
        % record the x-y limits from display axes
        
%         x_length = display_axes.XLim(2) - display_axes.XLim(1);
%         y_length = display_axes.YLim(2) - display_axes.YLim(1);
% 
%         display_length = min( x_length, y_length );
% 
%         FOV_limits( dimension_order( 2 ),  = [display_axes.XLim(1) display_axes.XLim(1) + display_length];
%         FOV_limits( dimension_order( 1 ),  = [display_axes.YLim(1) display_axes.YLim(1) + display_length];
% 
%         display_axes.XLim(2) = FOV_limits( dimension_order( 2 ), (2);
%         display_axes.YLim(2) = FOV_limits( dimension_order( 1 ), (2);

        FOV_limits( dimension_order( 2 ), 1 : 2 ) = [ display_axes.XLim( 1 ), display_axes.XLim( 2 )];
        FOV_limits( dimension_order( 1 ), 1 : 2 ) = [ display_axes.YLim( 1 ), display_axes.YLim( 2 )];
        
        ylim([ 1, size_of_image( dimension_order( 1 ))])
        xlim([ 1, size_of_image( dimension_order( 2 ))])
        
        zoom( 'reset' )
        
        ylim( FOV_limits( dimension_order( 1 ), : ))
        xlim( FOV_limits( dimension_order( 2 ), : ))        

        in_view_edges_xy = edge_centers( :, dimension_order( 1 )) >= FOV_limits( dimension_order( 1 ), 1 ) ...
                         & edge_centers( :, dimension_order( 1 )) <= FOV_limits( dimension_order( 1 ), 2 ) ... 
                         & edge_centers( :, dimension_order( 2 )) >= FOV_limits( dimension_order( 2 ), 1 ) ... 
                         & edge_centers( :, dimension_order( 2 )) <= FOV_limits( dimension_order( 2 ), 2 ) ;

        in_view_edges = in_view_edges_xy & in_view_edges_z ;

        % propogate the new display to the histograms and mini-map
        update_intesnity_histogram, update_energy_histogram, update_minimap, main_figure_size_update

    end % update_xy_zoom

    function change_depth_of_view( source, ~ )

        % update the z depth from the slide value, then update the crop
        depth_of_view = round( source.Value ); 

        update_z_zoom

    end % change_depth_of_view

    function change_view_thickness( source, ~ )

        % update the z thickness from the slide values
        view_thickness = round( exp( source.Value )) - 1 ; 

        update_z_zoom, update_z_sliders

    end % change_view_thickness

    function projection_dimension_callback( ~, ~ )

        view_thickness = FOV_limits( dimension_order( 3 ), 2 ) ...
                       - FOV_limits( dimension_order( 3 ), 1 ) ;  

        if view_thickness == 0

            FOV_limits( dimension_order( 3 ), 1 ) = max(               1,                      FOV_limits( dimension_order( 3 ), 1 ) - 1 );
            FOV_limits( dimension_order( 3 ), 2 ) = min( size_of_image( dimension_order( 3 )), FOV_limits( dimension_order( 3 ), 2 ) + 1 );

        end

        projection_dimension = mod( projection_dimension, 3 ) + 1 ;

%         dimension_letters = { 'Y', 'X', 'Z' };

%         dimension_order_matrix ...
%             = [ 2, 3, 1 ; ...% y-projection
%                 3, 1, 2 ; ...% x-projection
%                 1, 2, 3   ]; % dimension permutation to project in different dimensions

        dimension_letter = dimension_letters{ projection_dimension };

        set( view_thickness_label, 'string', [ dimension_letter, '-Thickness' ]);
        set(  depth_of_view_label, 'string', [ dimension_letter, '-Depth' ]);

        dimension_order = dimension_order_matrix( projection_dimension, : );

        number_of_pixels_2D = prod( size_of_image( dimension_order( 1 : 2 )));

        view_thickness = round((   FOV_limits( dimension_order( 3 ), 2 ) ...
                                 - FOV_limits( dimension_order( 3 ), 1 ) - 1 ) / 2 );

        depth_of_view = round((   FOV_limits( dimension_order( 3 ), 2 )      ...
                                + FOV_limits( dimension_order( 3 ), 1 )) / 2 );

        update_overlay, update_z_zoom, update_xy_zoom, update_z_sliders

    end

        %% utility functions                        
    %EDIT: added this function to globalize colormap
    function load_colormap
        n = 64;  
                                %make white spot at end of first quarter
        r1 = zeros(1,n);        r1(end-3:end) = 1;
        g1 = linspace(0,0.6,n); g1(end-3:end) = 1;
        b1 = linspace(0.6,0,n); b1(end-3:end) = 1;
%         a1 = linspace(1,1,n);

        r2 = linspace(0,1,n);
        g2 = linspace(0.6,1,n);
        b2 = linspace(0,0,n);
%         a2 = linspace(1,1,n); a2(end-3:end) = 0.25;

        r3 = linspace(1,0.9,n);
        g3 = linspace(0.75,0.25,n);
        b3 = linspace(0.75,0.25,n);
%         a3 = linspace(1,1,n);
        
        r4 = linspace(0.9,0.5,n);
        g4 = linspace(0.25,0,n);
        b4 = linspace(0.25,0,n);
%         a4 = linspace(1,1,n);
        
        colormap_mat = [ r1(:), g1(:), b1(:); ...
                         r2(:), g2(:), b2(:); ...
                         r3(:), g3(:), b3(:); ...
                         r4(:), g4(:), b4(:)  ];
         
        %color selection tool
        threshold_color_range = [0.5 1];
        out_of_bounds_color = 0.0;
        border_color = 0.245;
        in_view_color = 0.245;
                     
    end
        
    function main_figure_size_update( ~, ~) 
        
        display_figure_position = get(display_figure,'Position');
        
%         display_axes.Position = [display_figure_position(3)*0.075, display_figure_position(4)*0.15, min(display_figure_position(3)*0.85,display_figure_position(4)*0.85), min(display_figure_position(3)*0.85,display_figure_position(4)*0.85)];        
        
        image_aspect_x_over_y = microns_per_pixel( dimension_order( 2 )) * ( FOV_limits( dimension_order( 2 ), 2 ) - FOV_limits( dimension_order( 2 ), 1 ) + 1 ) ...
                              / microns_per_pixel( dimension_order( 1 )) / ( FOV_limits( dimension_order( 1 ), 2 ) - FOV_limits( dimension_order( 1 ), 1 ) + 1 );

        display_figure_position_min = min( display_figure_position( 3 ) / image_aspect_x_over_y, ...
                                           display_figure_position( 4 ));
        
        display_axes.Position = [ display_figure_position( 3 ) * 0.075, ...
                                  display_figure_position( 4 ) * 0.15,  ...
                                  display_figure_position_min * 0.85 * image_aspect_x_over_y, ...
                                  display_figure_position_min * 0.85 ];
                              
    end

    function minimap_figure_size_update( ~, ~) 
        pos = get(minimap_figure,'Position');
        minimap_axes.Position = [pos(3)*0.2, pos(4)*0.2, min(pos(3)*0.7,pos(4)*0.7), min(pos(3)*0.7,pos(4)*0.7)];
    end

        %% automatic backup operation               
    function backup_curation

        save([ path_to_saved_curation, '_backup_', num2str( backup_ordinate )], autosave_variables{ : });

        backup_ordinate = mod( backup_ordinate + 1, maximum_number_of_backups ); 
        
        is_active = true ;

    end % backup_curation
        %% undo and redo                            
    function undo_curation( ~, ~ )

        backup_curation

        backup_ordinate = mod( backup_ordinate - 2, maximum_number_of_backups );

        try

            load([ path_to_saved_curation, '_backup_', num2str( backup_ordinate )]);

        catch

            warning('The requested earlier backup curation was not found')

            backup_ordinate = mod( backup_ordinate + 1, maximum_number_of_backups );

        end

        refresh_and_paint_index_image
        
    end % undo_curation

    function redo_curation( ~, ~ )

        backup_ordinate = mod( backup_ordinate + 1, maximum_number_of_backups );    

        try

            load([ path_to_saved_curation, '_backup_', num2str( backup_ordinate )]);

        catch

            warning('The requested latter backup curation was not found')

            backup_ordinate = mod( backup_ordinate - 1, maximum_number_of_backups );        

        end

        refresh_and_paint_index_image

    end % redo_curation
        %% save and load                            

    function store_curation( ~, ~ )

        listing_for_path_to_saved_curation = dir([ path_to_saved_curation, '.mat' ]);

        if ~ isempty( listing_for_path_to_saved_curation )

            answer = questdlg( 'Overwrite the previously saved curation?', 'Edge Curator' );

            switch answer, case { 'no', 'cancel', '' }, return, end

        end
        
        elapsed_time = elapsed_time + is_active * toc( tstart );
        
        tstart = tic ;

        save( path_to_saved_curation, save_variables{ : });

    end % store_curation

    function load_curation( ~, ~ )

        backup_curation

        load( path_to_saved_curation )

        display_axes.XLim = FOV_limits( dimension_order( 2 ), : );
        display_axes.YLim = FOV_limits( dimension_order( 1 ), : );

        initialize_structuring_elements

        refresh_and_paint_index_image
        
        tstart = tic ;        

    end % load_curation

        %% sweeping and painting operations         
    function sweep_index_image

        % move any in view edges, displayed false, and not yet deleted to the deleted category
        to_be_deleted_edges = ~ true_edges & in_view_edges & displayed_edges & ~ deleted_edges ;

        % If ever once displayed, the edge will remain in the displayed vector
%             displayed_edges = ~ to_be_deleted_edges & displayed_edges ; 
              deleted_edges =   to_be_deleted_edges |   deleted_edges ;
        
        refresh_and_paint_index_image
        
%         % background is coded by last_vertex_index plus one, set it to not deleted to avoid searching
%         % and deleting the background (and whatever bug that might create) in the following
%         to_be_deleted_edges( last_vertex_index + 1 ) = false ;        
%         
%         to_be_deleted_indices = find( to_be_deleted_edges )' ;
%         
%         for edge_index = to_be_deleted_indices
% 
%             % set all the indices belonging to deleted edges to background (background is last_vertex_index plus one)
%             index_image( edge_structure_positions_linear_indexing{ edge_index }) = last_vertex_index + 1 ;
%             
%             % repaint any vertices that still are connected to a not deleted edge
%             vertex_1_valid = any(   edges2vertices( :, 1 ) == edges2vertices( edge_index, 1 ) ...
%                                   | edges2vertices( :, 2 ) == edges2vertices( edge_index, 1 ) ...
%                                   & ~ deleted_edges                                           );
%             
%             vertex_2_valid = any(   edges2vertices( :, 1 ) == edges2vertices( edge_index, 2 ) ...
%                                   | edges2vertices( :, 2 ) == edges2vertices( edge_index, 2 ) ...
%                                   & ~ deleted_edges                                           );                              
%                               
%             if vertex_1_valid
%                 
%                 index_image( end_vertices_structure_positions_linear_indexing{ edge_index, 1 }) = edges2vertices( edge_index, 1 );
%                 
%             end
%             
%             if vertex_2_valid
%                 
%                 index_image( end_vertices_structure_positions_linear_indexing{ edge_index, 2 }) = edges2vertices( edge_index, 2 );
%             
%             end
%         end % FOR deleted vertex indices
    end % sweep_index_image

    function paint_index_image
        
        % attempt to paint any vertices in the current view that aren't displayed already or deleted
        to_be_painted_edges = find( ~ ( deleted_edges | displayed_edges ) & in_view_edges )';

        render_edges_in_index_image( to_be_painted_edges )
        
    end % FUNCTION paint_index_image

    function refresh_and_paint_index_image

        % select the vertices to attempt to paint in the following FOR loop
        to_be_painted_edges = find( displayed_edges & ~ deleted_edges )';

        % % main display image initializtion, background is coded by last_vertex_index plus one
        % index_image = 1 + last_vertex_index * ones( size_of_image, 'uint32' );

        % main display image initializtion, background is coded by zero
              index_image = zeros( size_of_image, 'int32'   );
        
        centerlines_image = zeros( size_of_image, 'logical' );

        render_edges_in_index_image( to_be_painted_edges )

        update_z_zoom, update_xy_zoom, update_z_sliders
        update_edge_intensities, update_energy_histogram, update_minimap

    end % FUNCTION refresh_and_paint_index_image

    function refresh_index_image_all_vertices

        % select the vertices to attempt to paint in the following FOR loop
        to_be_painted_edges = find( displayed_edges & ~ deleted_edges )';

        % % main display image initializtion, background is coded by last_vertex_index plus one
        % index_image = 1 + last_vertex_index * ones( size_of_image, 'uint32' );

        % main display image initializtion, background is coded by zero
        index_image = zeros( size_of_image, 'int32' );

        % loop through the vertices from best to worst (lowest to highest energy). Should be sorted before
        % this function.
        for edge_index = to_be_painted_edges 
            
%             index_image( edge_structure_positions_linear_indexing{ edge_index }) = last_vertex_index + edge_index ;
            index_image( edge_structure_positions_linear_indexing{ edge_index }) = edge_index ;
            
        end % edge FOR
        
        to_be_painted_vertices = 1 : length( vertex_scale_subscripts ) ;
        
        for vertex_index = to_be_painted_vertices
            
%             index_image( vertex_structure_positions_linear_indexing{ vertex_index }) = vertex_index ;
            index_image( vertex_structure_positions_linear_indexing{ vertex_index }) = - vertex_index ;

        end % vertex FOR

        update_z_zoom, update_xy_zoom, update_z_sliders
        update_edge_intensities, update_energy_histogram, update_minimap

    end % FUNCTION paint_from_scratch

    function render_edges_in_index_image( to_be_painted_edges )
        
        % loop through the vertices from best to worst (lowest to highest energy). Should be sorted
        % before this function.        
        for edge_A_index = to_be_painted_edges 
            
% SAM 7/22/19, removing volume exclusion            
%             index_image_snapshot = index_image( cell2mat( vertex_structure_positions_linear_indexing( edges2vertices( edge_A_index, : )')));
%             
%             % this is slightly inconsistent with the choose_edges_V200 function where the vertices
%             % for each edge are defined as the first and last positions of that edge
% %             index_image( vertex_structure_positions_linear_indexing{ edges2vertices( edge_A_index, 1 )}) = edges2vertices( edge_A_index, 1 );
% %             index_image( vertex_structure_positions_linear_indexing{ edges2vertices( edge_A_index, 2 )}) = edges2vertices( edge_A_index, 2 );            
%             
%             index_image( vertex_structure_positions_linear_indexing{ edges2vertices( edge_A_index, 1 )}) = - edges2vertices( edge_A_index, 1 );
%             index_image( vertex_structure_positions_linear_indexing{ edges2vertices( edge_A_index, 2 )}) = - edges2vertices( edge_A_index, 2 );            
% 
%             index_image_at_edge = index_image( edge_structure_positions_linear_indexing{ edge_A_index });
%                         
% %             number_of_voxels_at_edge = numel( index_image_at_edge );
% 
%             % check if this volume has any previously painted vertices that don't belong to this edge
% %             if all(   index_image_at_edge >  last_vertex_index                 ...
% %                     | index_image_at_edge == edges2vertices( edge_A_index, 1 ) ...
% %                     | index_image_at_edge == edges2vertices( edge_A_index, 2 ) )
% 
% % SAM 7/22/19, removing volume exclusion
% %             if all(   index_image_at_edge >  0                                   ...
% %                     | index_image_at_edge == - edges2vertices( edge_A_index, 1 ) ...
% %                     | index_image_at_edge == - edges2vertices( edge_A_index, 2 ) )
%                 
%                 % identify the edge objects in conflict with this edge           
%                 unique_indices_at_edge = unique( index_image_at_edge );
%                 
%                 % exclude the background and vertices
% %                 unique_indices_at_edge ...
% %                                 = unique_indices_at_edge( unique_indices_at_edge > last_vertex_index + 1 );
%                 unique_indices_at_edge ...
%                                 = unique_indices_at_edge( unique_indices_at_edge > 0 );
%                                                            
%                 % translate the index image convention to the edge listing convention
% %                 edge_indices_at_edge = unique_indices_at_edge - last_vertex_index ;
%                 edge_indices_at_edge = unique_indices_at_edge ;
%                 
%                 % isolate the edge_structures for those edges in conflict with this edge
%                 edge_structures_at_edge = edge_structure_positions_linear_indexing( edge_indices_at_edge );
%                 
%                 % measure the number of voxels in each edge structure in conflict                
% %                 edge_numels_at_edge = cellfun( @numel, edge_structures_at_edge );
% 
% %                 % total number of edge voxels considered is those that are in the edge but not
% %                 % either of its vertices.                
% %                 edge_numels_at_edge = cellfun( @( x, y, z ) setdiff( x, union( y, z )), edge_structures_at_edge );
%                 
%                 number_of_edge_conflicts = length( edge_indices_at_edge );
%                 
%                 edge_conflict_ordinates = 1 : number_of_edge_conflicts ;                
%                 
% %                 edge_conflict_percentages = zeros( number_of_edge_conflicts, 1 );
%                 edge_conflict_percentages    = zeros( number_of_edge_conflicts + 1, 1 );
%                 
%                 number_of_conflicting_voxels = zeros( number_of_edge_conflicts + 1, 1 );
% 
%                 for edge_conflict_ordinate = edge_conflict_ordinates
%                     
%                     % total number of edge voxels considered is those that are in the edge but not
%                     % either of its vertices.                
%                     edge_numels_at_edge = numel( setdiff(                                                                    edge_structures_at_edge{ edge_conflict_ordinate },       ...
%                                                           cell2mat( vertex_structure_positions_linear_indexing( edges2vertices( edge_indices_at_edge( edge_conflict_ordinate ), : )'))));                    
%                     
%                     % count the voxels at edge A that are in conflict with edge B
%                     number_of_conflicting_voxels( edge_conflict_ordinate, 1 ) = sum( index_image_at_edge == unique_indices_at_edge( edge_conflict_ordinate ));
%                     
%                     % normalized by edge B
%                     edge_conflict_percentages( edge_conflict_ordinate, 1 ) = number_of_conflicting_voxels( edge_conflict_ordinate, 1 ) / edge_numels_at_edge ;
%                     
% %                     % normalized by edge A
% %                     edge_conflict_percentages( edge_conflict_ordinate, 2 ) = number_of_conflicting_voxels / number_of_voxels_at_edge ;
%                     
%                 end % edge B FOR
%                 
%                 edge_numels_at_edge_A = numel( setdiff(             edge_structure_positions_linear_indexing{                 edge_A_index },    ...
%                                                         cell2mat( vertex_structure_positions_linear_indexing( edges2vertices( edge_A_index, : )'))));  
%                 
% %                 edge_conflict_percentages( edge_conflict_ordinate + 1 ) = sum( number_of_conflicting_voxels ) / number_of_voxels_at_edge ;
%                 edge_conflict_percentages( edge_conflict_ordinate + 1 ) = sum( number_of_conflicting_voxels ) / edge_numels_at_edge_A ;
% 
%                 if max( edge_conflict_percentages( : )) > threshold_percent_overlap
%                     
%                      displayed_edges( edge_A_index ) = false ;
%                     
%                 else %  volume has no significant conflicts with other edges
                    
                    % The reason for the IF/ELSE (instead of just IF):  This way, this line is also
                    % reached when edge_conflict_percentages is empty (no edge conflicts)
                    
                    displayed_edges( edge_A_index ) = true ;

                    % paint the edge volume first with its unique index in the index image to mark and then
                    % overwrite the beginning and ending vertices with their unique indices
%                     index_image( edge_structure_positions_linear_indexing{ edge_A_index }) = last_vertex_index + edge_A_index ;
                    index_image( edge_structure_positions_linear_indexing{ edge_A_index }) = edge_A_index ;

        %             index_image( vertex_structure_positions_linear_indexing{ edges2vertices( edge_A_index, 1 )}) = edges2vertices( edge_A_index, 1 );
        %             index_image( vertex_structure_positions_linear_indexing{ edges2vertices( edge_A_index, 2 )}) = edges2vertices( edge_A_index, 2 );            

                    index_image( vertex_structure_positions_linear_indexing{ edges2vertices( edge_A_index, 1 )}) = - int32(edges2vertices( edge_A_index, 1 ));
                    index_image( vertex_structure_positions_linear_indexing{ edges2vertices( edge_A_index, 2 )}) = - int32(edges2vertices( edge_A_index, 2 ));            
                    
                    centerlines_image( edge_positions_linear_indices{ edge_A_index }) = true ;
                    
%                 end % IF volume has significant conflicts with other edges                
%             else % volume conflicts with at least one vertex of a higher energy edge
%                 
%                 displayed_edges( edge_A_index ) = false ;
%                 
%                 % put the image back how it was                
%                 index_image( cell2mat( vertex_structure_positions_linear_indexing( edges2vertices( edge_A_index, : )'))) = index_image_snapshot ;
%                 
%             end % IF volume is blank canvas (no conflicts with vertices of a higher energy edge)
        end % edge A FOR
    end % FUNCTION render_edges_in_index_image

    function paint_new_edge( to_be_painted_edge_index )
        
        space_subscript_int64 = int64( edge_space_subscripts{ to_be_painted_edge_index });

        edge_positions_linear_index =   space_subscript_int64( :, 1 )                                               ...
                                    + ( space_subscript_int64( :, 2 ) - 1 ) * size_of_image( 1 )                    ...
                                    + ( space_subscript_int64( :, 3 ) - 1 ) * size_of_image( 1 ) * size_of_image( 2 );

        edge_positions_linear_index_cell = num2cell( edge_positions_linear_index );

        edge_positions_linear_indices{ to_be_painted_edge_index } = edge_positions_linear_index ;
        
        edge_structuring_element_linear_indexing = edge_element_linear_indexing_templates( uint8( edge_scale_subscripts{ to_be_painted_edge_index }));

        sphere_structure_positions_linear_indexing_for_edge = cellfun( @plus, edge_structuring_element_linear_indexing, ...
                                                                                      edge_positions_linear_index_cell, ...
                                                                                                'UniformOutput', false  ); 

%         sphere_structure_positions_linear_indexing( to_be_painted_edge_index, : ) = sphere_structure_positions_linear_indexing_for_edge([ 1, end ]);

        edge_structure_positions_linear_indexing{ to_be_painted_edge_index } = cell2mat( sphere_structure_positions_linear_indexing_for_edge );

        edge_structure_positions_linear_indexing{ to_be_painted_edge_index } = unique( edge_structure_positions_linear_indexing{ to_be_painted_edge_index }, 'rows' );

        edge_centers( to_be_painted_edge_index, : ) = round( mean( edge_structure_positions_linear_indexing{ to_be_painted_edge_index }));
        
        % volume exclusion for edge body: need to add the percentage based volume exclusion for the
        % edge bodies (not just their vertices).
        
        % zero out the vertices previous to checking the edge volume
%         index_image( vertex_structure_positions_linear_indexing{ edges2vertices( to_be_painted_edge_index, 1 )}) = last_vertex_index + 1 ;
%         index_image( vertex_structure_positions_linear_indexing{ edges2vertices( to_be_painted_edge_index, 2 )}) = last_vertex_index + 1  ;        
        index_image( vertex_structure_positions_linear_indexing{ edges2vertices( to_be_painted_edge_index, 1 )}) = 0 ;
        index_image( vertex_structure_positions_linear_indexing{ edges2vertices( to_be_painted_edge_index, 2 )}) = 0 ;        
        
        % look at the edge volume to be painted and see if there are vertices in the way to delete
%         indices_in_edge_volume = unique( index_image( edge_structure_positions_linear_indexing{ to_be_painted_edge_index }))';
        indices_in_edge_volume = - unique( index_image( edge_structure_positions_linear_indexing{ to_be_painted_edge_index }))';

%         to_be_deleted_vertex_indices = indices_in_edge_volume(    indices_in_edge_volume <= last_vertex_index                             ...
%                                                                &  indices_in_edge_volume ~= edges2vertices( to_be_painted_edge_index, 1 ) ...
%                                                                &  indices_in_edge_volume ~= edges2vertices( to_be_painted_edge_index, 2 ));
        to_be_deleted_vertex_indices = indices_in_edge_volume(    indices_in_edge_volume >  0                                               ...
                                                               &  indices_in_edge_volume ~= edges2vertices( to_be_painted_edge_index, 1 ) ...
                                                               &  indices_in_edge_volume ~= edges2vertices( to_be_painted_edge_index, 2 ));
                                                           
        to_be_deleted_vertex_indices = [ ]; % SAM 7/22/19, removing volume exclusion

        % once vertices are detected in the edge volume (that do not belong to this edge), erase any
        % edges corresponding to those vertices
        for vertex_index = to_be_deleted_vertex_indices

            to_be_deleted_edge_indices = mod( find( edges2vertices == vertex_index ), length_of_edge_lists )';

            % set all the indices belonging to deleted edges to background (background is
            % last_vertex_index plus one)
            for edge_index = to_be_deleted_edge_indices

%                 index_image( edge_structure_positions_linear_indexing{ edge_index }) = last_vertex_index + 1 ;
                      index_image( edge_structure_positions_linear_indexing{ edge_index }) = 0 ;
                
                centerlines_image(             edge_positions_linear_indices{ edge_index }) = false ;
                
                  deleted_edges( edge_index ) = true  ;
                displayed_edges( edge_index ) = false ;

            end % FOR deleted edge indices
        end % FOR deleted vertex indices

        % paint the edge volume first with its unique index in the index image to mark and then
        % overwrite the beginning and ending vertices with their unique indices
%         index_image( edge_structure_positions_linear_indexing{ to_be_painted_edge_index }) = last_vertex_index + to_be_painted_edge_index ;
%         
%         index_image( vertex_structure_positions_linear_indexing{ edges2vertices( to_be_painted_edge_index, 1 )}) = edges2vertices( to_be_painted_edge_index, 1 );
%         index_image( vertex_structure_positions_linear_indexing{ edges2vertices( to_be_painted_edge_index, 2 )}) = edges2vertices( to_be_painted_edge_index, 2 );
        index_image( edge_structure_positions_linear_indexing{ to_be_painted_edge_index }) = to_be_painted_edge_index ;
        
        index_image( vertex_structure_positions_linear_indexing{ edges2vertices( to_be_painted_edge_index, 1 )}) = - edges2vertices( to_be_painted_edge_index, 1 );
        index_image( vertex_structure_positions_linear_indexing{ edges2vertices( to_be_painted_edge_index, 2 )}) = - edges2vertices( to_be_painted_edge_index, 2 );

        centerlines_image( edge_positions_linear_indices{ to_be_painted_edge_index }) = true ;
        
    end % paint_new_edge
        %% sweeping and painting                    
    function sweep_callback( ~, ~ )

        backup_curation

        % Find the vertices that are false, delete them, don't display them, and turn them into
        % background in the 3D index image.
        sweep_index_image

        update_index_image_crop, update_overlay, update_energy_histogram

        % re-initialize energy histogram and energy_limits for visualization and thresholding

    end % sweep_falses

    function paint_callback( ~, ~ )

        backup_curation

        % repopulate the current field of view with any vertices that have not yet been deleted and are
        % not currently displayed.  Start with the best contrast objects and go down the list, trying to
        % place new objects where the volumes don't conflict
        paint_index_image

        update_index_image_crop, update_overlay, update_energy_histogram

    %     add_repop_to_mini_map

    end
        %% toggling and adding                      
    function toggle_callback( ~, ~ ) 

        toggle_button.Enable = 'off';

        % prompt user to click the image
        [ x_click, y_click ] = ginput( 1 );

        % Get the edge index from the position of the mouse click. 
        if    y_click > FOV_limits( dimension_order( 1 ), 1 ) && y_click < FOV_limits( dimension_order( 1 ), 2 ) ...
           && x_click > FOV_limits( dimension_order( 2 ), 1 ) && x_click < FOV_limits( dimension_order( 2 ), 2 )

            % This will fail if user clicks outside the image. Exit the recursive function in that case.
            try edge_index = index_image_2D( round( y_click ), round( x_click )); catch, return, end

        else % IF click 

            toggle_button.Enable = 'on';
            return % Exit the recursive function

        end % IF click inside current field of view

        backup_curation

        true_edges( edge_index ) = ~ true_edges( edge_index );
        
%         % section added SAM 9/4/20
%         edge_space_subscripts_temp = edge_space_subscripts ;
%         
%         edge_space_subscripts_temp{ ~ true_edges } = [ ];
%                                 
%         % !!! will empty edges pass thru nicely?
%         chosen_edge_indices = clean_edges_orphans( edge_space_subscripts_temp, size_of_image, vertex_space_subscripts );
% 
%         rejected_edge_indices = setdif( 1 : length( edge_space_subscripts_temp ), chosen_edge_indices );
%         
%         true_edges( rejected_edge_indices ) = false ;
%         
%         % end: section added SAM 9/4/20

        counter_toggle = counter_toggle + 1 ;        
        
        update_truth_image, update_overlay

        % recall this function to select anew
        toggle_callback;

    end % vertex_toggler

%     function add_callback( ~, ~ )
% 
%         % !!!!!! this function needs to be updated to create a temporary vertex acting as seed point to
%         % a monte carlo random walk through the energy image.  the top two associated extant vertices
%         % associated with this edge will be the beginning and ending vertices for the new edge to be
%         % created.  Also, update the degrees_of_edges global variable for the added index.
% 
%     %     % prompt user to click the image
%     %     [ x_click, y_click ] = ginput( 1 );
%     % 
%     %     % Get the vertex index from the position of the mouse click. 
%     %     if    y_click > FOV_limits( dimension_order( 1 ), ( 1 ) && y_click < FOV_limits( dimension_order( 1 ), ( 2 ) ...
%     %        && x_click > FOV_limits( dimension_order( 2 ), ( 1 ) && x_click < FOV_limits( dimension_order( 2 ), ( 2 )
%     %    
%     %         backup_curation
%     %    
%     %         z_start = max( FOV_limits( dimension_order( 3 ), ( 1 ), 1 ) ;
%     %         z_end   = min( FOV_limits( dimension_order( 3 ), ( 2 ), size_of_image( dimension_order( 3 ) ));
%     %         z_count = z_end - z_start + 1 ;
%     %    
%     %         added_vertex_xy_subscripts = [ round( y_click ), round( x_click )];
%     %         
%     %         vertex_energy_and_scale_index = h52mat(                        path_to_energy_data, ...
%     %                                         	    [ added_vertex_xy_subscripts, z_start, 1 ], ...
%     %                                                                       [ 1, 1, z_count, 2 ]  ); 
%     %                                                                           
%     %         [ added_vertex_energy, added_vertex_z_relative ]                                            ...
%     %                                                   = min( vertex_energy_and_scale_index( z_count + 1 : end ));
%     %         
%     %         added_vertex_scale_index = vertex_energy_and_scale_index( added_vertex_z_relative );
%     %         
%     %         added_vertex_z_subscript = added_vertex_z_relative + z_start - 1 ;
%     %                                                                           
%     %         to_be_painted_index = number_of_added_edges + 2 ;
%     %         
%     %         edge_energies(         to_be_painted_index    ) = added_vertex_energy ;
%     %         edge_scale_subscripts( to_be_painted_index    ) = added_vertex_scale_index ;
%     %         edge_space_subscripts( to_be_painted_index, : ) = [ added_vertex_xy_subscripts, ...
%     %                                                               added_vertex_z_subscript    ]; 
%     %                                                           
%     %         %% need to add warning like vertex curator example SAM+WAS 11/19/18
%     %                                                           
%     %         paint_new_edge( to_be_painted_index )
%     %                                         
%     %           deleted_edges( to_be_painted_index ) = false ;
%     %         displayed_edges( to_be_painted_index ) = true ;
%     %         
%     %         update_vertex_intensities
%     %         
%     %         number_of_added_edges = number_of_added_edges + 1 ;          
%     %                                      
%     %         update_index_image_2D_crop, update_overlay, update_energy_histogram
%     % 
%     %         % recall this function to select anew
%     %         add_callback        
%     %                    
%     %     end % IF click inside current field of view
%     end % edge_placer

        %% edge_placer (simple)
    function add_callback( ~, ~ ) 
        
        add_button.Enable = 'off';        
        
        % update display to show all possible vertices, (not just the ones with edges attached)
        refresh_index_image_all_vertices
        
        % prompt user to click the image
        [ x_click, y_click ] = ginput( 1 );

        % Get the vertex index from the position of the mouse click. 
        if    y_click > FOV_limits( dimension_order( 1 ), 1 ) && y_click < FOV_limits( dimension_order( 1 ), 2 ) ...
           && x_click > FOV_limits( dimension_order( 2 ), 1 ) && x_click < FOV_limits( dimension_order( 2 ), 2 )

            % This will fail if user clicks outside the image. Exit the recursive function in that case.
            try vertex_index_1 = vertex_index_image_2D( round( y_click ), round( x_click )); catch, return, end
            
            % !!!! put the Stopper here to navigate in the volume between vertex selections
            
            % prompt user to click the image
            [ x_click, y_click ] = ginput( 1 );         
            
            % Get the vertex index from the position of the mouse click. 
            if    y_click > FOV_limits( dimension_order( 1 ), 1 ) && y_click < FOV_limits( dimension_order( 1 ), 2 ) ...
               && x_click > FOV_limits( dimension_order( 2 ), 1 ) && x_click < FOV_limits( dimension_order( 2 ), 2 )

                % This will fail if user clicks outside the image. Exit the recursive function in that case.
                try vertex_index_2 = vertex_index_image_2D( round( y_click ), round( x_click )); catch, return, end
                
            else % IF click 

                refresh_and_paint_index_image
                add_button.Enable = 'on';
                return % Exit the recursive function

            end % IF click inside current field of view                
        else % IF click 

            refresh_and_paint_index_image
            add_button.Enable = 'on';
            return % Exit the recursive function

        end % IF click inside current field of view

        % check if user selected two vertices
        if vertex_index_1 > 0 && vertex_index_2 > 0
        
            backup_curation

            new_edge_subscripts = [ double( vertex_space_subscripts( vertex_index_1, 1 : 3 )), vertex_scale_subscripts( vertex_index_1 ); ...
                                    double( vertex_space_subscripts( vertex_index_2, 1 : 3 )), vertex_scale_subscripts( vertex_index_2 )  ];

            cumulative_lengths = [ 0; sum( max( abs((   new_edge_subscripts( 2, 1 : 3 )            ...
                                                      - new_edge_subscripts( 1, 1 : 3 ))), [ ], 2 ))];

            % interpolate the vector positions and size indices along the edge index, with even spacing
            % along units of voxel lengths (anisotropically)
            new_edge_subscripts_interpolated                               ...
                    = interp1( cumulative_lengths,                         ...
                               new_edge_subscripts,                        ...
                               linspace( 0, cumulative_lengths( 2 ),       ...
                                            cumulative_lengths( 2 ) + 1 )' );

            to_be_painted_edge = number_of_added_edges + 2 ;

            last_vertex_index = max( edges2vertices( : ));

            % erase any edge that already connects these two vertices (must check A->B and B->A)
            old_adjacency_indices = double( edges2vertices( :, 1 )) + double(last_vertex_index * ( edges2vertices( :, 2 ) - 1 ));

            new_adjacency_index   =  [  double(  vertex_index_1  )  + last_vertex_index * ( double(    vertex_index_2     - 1 )), ...
                                        double(  vertex_index_2  )  + last_vertex_index * ( double(    vertex_index_1     - 1 ))  ];

            is_edge_redundant = old_adjacency_indices == new_adjacency_index( 1 ) | old_adjacency_indices == new_adjacency_index( 2 );

              deleted_edges( is_edge_redundant ) = true  ;
            displayed_edges( is_edge_redundant ) = false ; % this edge likely was unnoticed by curator

            % record the new edge in the edge lists
                edges2vertices( to_be_painted_edge, : ) = [ vertex_index_1, vertex_index_2 ];

            mean_edge_energies( to_be_painted_edge, : ) = -Inf ;

    % % for the energy convention, - stdev / mean, the best response is 0.
    %         mean_edge_energies( to_be_painted_edge, : ) = 0 ;       

         edge_space_subscripts{ to_be_painted_edge } = new_edge_subscripts_interpolated( :, 1 : 3 );
         edge_scale_subscripts{ to_be_painted_edge } = new_edge_subscripts_interpolated( :,   4   );

                 edge_energies{ to_be_painted_edge } = - Inf * ones( size( new_edge_subscripts_interpolated( :, 1 )));
    %           degrees_of_edges( to_be_painted_edge ) =       uint16( size( new_edge_subscripts_interpolated, 1     ));

            number_of_added_edges = number_of_added_edges + 1 ;

            paint_new_edge( to_be_painted_edge );

              deleted_edges( to_be_painted_edge ) = false ;
            displayed_edges( to_be_painted_edge ) = true  ;

            counter_add_edge = counter_add_edge + 1 ;
            
        end

        refresh_and_paint_index_image

        add_button.Enable = 'on';         

    end % edge_placer (simple)

        %% vertex placer
    function add_vertex_callback( ~, ~ )
        
        refresh_index_image_all_vertices

        % prompt user to click the image
        [ x_click, y_click ] = ginput( 1 );

        % Get the vertex index from the position of the mouse click. 
        if    y_click > FOV_limits( dimension_order( 1 ), 1 ) && y_click < FOV_limits( dimension_order( 1 ), 2 ) ...
           && x_click > FOV_limits( dimension_order( 2 ), 1 ) && x_click < FOV_limits( dimension_order( 2 ), 2 )
       
            z_start = max( FOV_limits( dimension_order( 3 ), 1 ), 1 ) ;
            z_end   = min( FOV_limits( dimension_order( 3 ), 2 ), size_of_image( dimension_order( 3 ) ));
            z_count = z_end - z_start + 1 ;

            added_vertex_xy_subscripts = [ round( y_click ), round( x_click )];

            starts([ dimension_order, 4 ]) = [ added_vertex_xy_subscripts, z_start, 1 ];
            counts([ dimension_order, 4 ]) = [ 1, 1, z_count, 2 ];
                        
            vertex_energy_and_scale_index = h52mat( path_to_energy_data, starts, counts ); 

            [ added_vertex_energy, added_vertex_z_relative ]                                            ...
                                                      = min( vertex_energy_and_scale_index( z_count + 1 : end ));

            added_vertex_scale_index = vertex_energy_and_scale_index( added_vertex_z_relative );

            added_vertex_z_subscript = added_vertex_z_relative + z_start - 1 ;
            
            to_be_painted_subscripts( dimension_order ) = uint16([ added_vertex_xy_subscripts, ...
                                                                   added_vertex_z_subscript    ]);

%            to_be_painted_index = number_of_added_vertices + 2 ;
            to_be_painted_index = numel(vertex_scale_subscripts) + 1;

            radius_in_pixels                                                                              ...
                    = uint16( round(    lumen_radius_in_microns_range( round( added_vertex_scale_index )) ...
                                     ./ microns_per_pixel                                                 ));

            % predict which vertices will overlap the outside of the image
            subscript_max = to_be_painted_subscripts + radius_in_pixels ;
            subscript_min = to_be_painted_subscripts - radius_in_pixels ;

            y_is_over  = subscript_max( :, 1 ) > size_of_image( 1 );
            x_is_over  = subscript_max( :, 2 ) > size_of_image( 2 );
            z_is_over  = subscript_max( :, 3 ) > size_of_image( 3 );

            y_is_under = subscript_min( :, 1 ) <         1         ;
            x_is_under = subscript_min( :, 2 ) <         1         ;
            z_is_under = subscript_min( :, 3 ) <         1         ;

            % don't include vertices that touch the image boundary
            excluded_vertex_flag = y_is_over | x_is_over | z_is_over | y_is_under | x_is_under | z_is_under ;

            if excluded_vertex_flag

                warning('Attempted to add object across the image boundary.')

            else % vertex will not cross the image boundary

                backup_curation

%                vertex_energies(         to_be_painted_index    ) = added_vertex_energy      ;
                vertex_scale_subscripts( to_be_painted_index    ) = added_vertex_scale_index ;
                vertex_space_subscripts( to_be_painted_index, : ) = to_be_painted_subscripts ;                 
                
%                paint_new_vertex( to_be_painted_index )

                space_subscript_int64 = int64( vertex_space_subscripts( to_be_painted_index, : ));

                position_linear_index =   space_subscript_int64( :, 1 )                                               ...
                                      + ( space_subscript_int64( :, 2 ) - 1 ) * size_of_image( 1 )                    ...
                                      + ( space_subscript_int64( :, 3 ) - 1 ) * size_of_image( 1 ) * size_of_image( 2 );

                structuring_element_linear_indexing = vertex_element_linear_indexing_templates{ round( vertex_scale_subscripts( to_be_painted_index ))};
                
                vertex_structure_positions_linear_indexing{ to_be_painted_index } = structuring_element_linear_indexing + position_linear_index ;
                                                                                         
                index_image( vertex_structure_positions_linear_indexing{ to_be_painted_index }) = - to_be_painted_index ;

                counter_add_vertex = counter_add_vertex + 1 ;                
                
                update_z_zoom, update_xy_zoom, update_z_sliders
                update_edge_intensities, update_energy_histogram, update_minimap
                
%                  deleted_vertices( to_be_painted_index ) = false ;
%                displayed_vertices( to_be_painted_index ) = true ;

%                update_vertex_intensities

%                number_of_added_vertices = number_of_added_vertices + 1 ;          

%                 update_index_image_2D_crop, update_overlay
                %update_energy_histogram
%                 refresh_index_image_all_vertices
                
                % recall this function to select anew
                add_vertex_callback                        

            end % vertex would cross the image boundary
        end % IF click inside current field of view
        refresh_and_paint_index_image
    end % vertex_placer

        %% intensity histogram                      

    function invert_callback( ~, ~) 

        is_intensity_inverted = ~ is_intensity_inverted;

        if is_intensity_inverted

            invert_button.String = 'Inverted' ;

        else

            invert_button.String = 'Original' ;

        end

        update_original_image_intensity, update_overlay

    end %invert_callback
    
    function change_intensity_min( source, ~ )

        intensity_limits( 1 ) = str2num( source.String );

        update_original_image_intensity, update_overlay

    end % change_intensity_limits

    function change_intensity_max( source, ~ )

        intensity_limits( 2 ) = str2num( source.String );

        update_original_image_intensity, update_overlay

    end % change_intensity_limits
        %% energy histogram      
    function binarize_callback( ~, ~) 
        
%         are_intensity_limits_binarized = ~ are_intensity_limits_binarized;
%         
%         if are_intensity_limits_binarized
%             
%             binarize_button.String = 'Binary' ;
%             
%         else
%             
%             binarize_button.String = 'Graded' ;
%             
%         end

        edge_display_method = mod( edge_display_method, 3 ) + 1 ;
        
        switch edge_display_method
            
            case 1
                
                binarize_button.String = 'Binary' ;
                
                
            case 2
                
                binarize_button.String = 'Graded' ;
                
            case 3

                binarize_button.String = 'Unique' ;
                
            case 4
                
        end
         
        update_edge_intensities, update_overlay
        
    end %binarize_callback	
	
    function change_energy_min( source, ~ )

        energy_limits( 1 ) = str2num( source.String );

        if energy_limits( 1 ) >= energy_limits( 2 ), warning('The energy min should be below the max'), end    

        update_edge_intensities, update_overlay, update_energy_histogram

    end % change_intensity_limits

    function change_energy_max( source, ~ )

        energy_limits( 2 ) = str2num( source.String );

        if energy_limits( 1 ) >= energy_limits( 2 ), warning('The energy max should be above the min'), end        

        update_edge_intensities, update_overlay, update_energy_histogram

    end % change_intensity_limits

    function apply_energy_threshold( source, ~ )

        backup_curation        

        energy_threshold = max( str2num( source.String ), min( mean_edge_energies ));

        true_edges( in_view_edges ) = mean_edge_energies( in_view_edges ) ...
                                    < energy_threshold                    ;
                    
%         % section added SAM 9/4/20
%         edge_space_subscripts_temp = edge_space_subscripts ;
%         
%         edge_space_subscripts_temp{ ~ true_edges } = [ ];
%                                 
%         % !!! will empty edges pass thru nicely?
%         chosen_edge_indices = clean_edges_orphans( edge_space_subscripts_temp, size_of_image, vertex_space_subscripts );
% 
%         rejected_edge_indices = setdif( 1 : length( edge_space_subscripts_temp ), chosen_edge_indices );
%         
%         true_edges( rejected_edge_indices ) = false ;
%         
%         % end: section added SAM 9/4/20
        
        update_threshold_visualization_matricies, update_minimap
        update_truth_image, update_overlay, update_energy_histogram

        counter_threshold = counter_threshold + 1 ;        
        
        if energy_threshold > 0 || imag( energy_threshold ) ~= 0, warning( 'energy threshold should be non-positive real number' ), end

        if energy_limits( 2 ) <= energy_threshold, warning('The energy threshold should be below the max'), end

        if energy_limits( 1 ) >= energy_threshold, warning('The energy threshold should be above the min'), end

    end % apply_energy_threshold

%% Finalizations                                

% wait for user to play on the curator
time_to_check_activity = 60 ; % seconds

is_curating = true ;

while is_curating

    is_active = false ;
    
    uiwait( display_figure, time_to_check_activity )
    
    elapsed_time = elapsed_time + is_active * toc( tstart );
    
    tstart = tic ;
    
    is_curating = ishandle(display_figure) && strcmp(get(display_figure, 'type'), 'figure');

end

is_active = false ; % time spent after closing curator is not active time (timer was just reset)

store_curation

chosen_edges = true_edges & displayed_edges & ~ deleted_edges ;

% background should not be chosen
chosen_edges( 1 ) = 0 ;

% performing the selections proposed by either choose_vertices or vertex_curator
edge_space_subscripts = edge_space_subscripts( chosen_edges    );
edge_scale_subscripts = edge_scale_subscripts( chosen_edges    );
edge_energies         =         edge_energies( chosen_edges    );    
mean_edge_energies    =    mean_edge_energies( chosen_edges    );  
edges2vertices        =        edges2vertices( chosen_edges, : );

close([ minimap_figure energy_histo_figure intensity_histo_figure ]);

end % main FUNCTION