function [ vertex_energies, vertex_space_subscripts, vertex_scale_subscripts, index_image, true_vertices, displayed_vertices, deleted_vertices ]         ...
            = vertex_curator_V6( vertex_energies, vertex_space_subscripts, vertex_scale_subscripts, ...
                                                  lumen_radius_in_microns_range, microns_per_pixel, ...
                                path_to_original_data, path_to_saved_curation, path_to_energy_data, ...
                                                                    intensity_limits, energy_limits )                                                     
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
% -Locked buttons and sliders into their final locations with ‘normalized’ ‘Position’ fields
% -Adjusted display_axes ‘Position’ field to give more room for buttons
% -Added toggling variable under global list (0 = false, 1 = true) and if statement to toggle_callback to make the toggle button toggleable (using the added variable, of course)
% -Implemented what I *think* displays all the circles as logical, rather than scaled values. But I will need to verify this with you
% 
% (10/24/18)
% -removed “toggling” global variable and all changes to toggle button made on 10/20/18
% -implemented an ideal toggling system by simply disabling the toggle button when the recursive toggle_callback is in operation
% -added minimap figure
% -made minimap resolution 100 by 100 to save on computation time (very slow at full res)
% -got it to render as we discussed
% -made minimap axes locked square just like display_axes
% -added color picker functionality for minimap (control+f “color selection tool”)
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
% V6: Edits made by WAS on 12/13/18:
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

%% Initializations   

max_number_of_added_vertices = 1000 ;

listing_for_path_to_saved_curation = dir([ path_to_saved_curation, '.mat' ]);

if ~ isempty( listing_for_path_to_saved_curation )

    answer = questdlg( 'Load the previously saved curation?', 'Vertex Curator' );

    if strcmp( answer, 'Yes' )
            
        load( path_to_saved_curation, 'vertex_energies', 'vertex_space_subscripts', 'vertex_scale_subscripts', 'max_number_of_added_vertices', 'intensity_limits', 'energy_limits' );
            
        vertex_energies         =         vertex_energies( max_number_of_added_vertices + 2 : end    );
        vertex_scale_subscripts = vertex_scale_subscripts( max_number_of_added_vertices + 2 : end    );        
        vertex_space_subscripts = vertex_space_subscripts( max_number_of_added_vertices + 2 : end, : );   
        
        load_upon_startup = true ;
    
    else
        
        load_upon_startup = false ;
        
    end
else
    
    load_upon_startup = false ;

end

    %% Figures                
bottom_buffer = 0.05 ;

left_buffer = bottom_buffer / 2 ;

spacer = 0.005 ;

histo_widths = 0.19 ;

histo_heights = 0.45 ;

half_screen = 1 - left_buffer - 2 * histo_widths - spacer ;

display_figure = figure( 'Name'     , 'Vertex Curator: Volume Display', ...
                         'Units'    , 'normalized',                     ...
                         'SizeChangedFcn',@main_figure_size_update,     ...
                         'Visible', 'off',                              ...
                         'OuterPosition' , [ left_buffer, bottom_buffer, half_screen - spacer - left_buffer, 2 * histo_heights + spacer ]);
    
display_figure.Units = 'pixels' ; 
% 
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
intensity_histo_figure = figure( 'Name'     , 'Vertex Curator: Intensity Histogram', ...
                                 'Units'    , 'normalized',                          ...
                                 'OuterPosition' , [ half_screen, bottom_buffer, histo_widths, histo_heights ]);

% intensity_histo_figure.Units = 'pixels' ; 
% 
% % force square dimensions
% intensity_histo_figure.Position([ 3, 4 ]) = min( intensity_histo_figure.Position([ 3, 4 ]));
% 
% intensity_histo_figure.Visible = 'off' ;

intensity_histo_axes = axes( 'OuterPosition', [ 0.1, 0.3, 0.8, 0.6 ]);

zoom( intensity_histo_figure, 'yon' )

% intensity_histo figure and axes initializations
energy_histo_figure = figure( 'Name'     , 'Vertex Curator: Energy Histogram', ...
                              'Units'    , 'normalized',                       ...
                              'OuterPosition' , [ half_screen + histo_widths + spacer, bottom_buffer, histo_widths, histo_heights ]);

% energy_histo_figure.Units = 'pixels' ; 
% 
% % force square dimensions
% energy_histo_figure.Position([ 3, 4 ]) = min( energy_histo_figure.Position([ 3, 4 ]));

energy_histo_figure.Visible = 'off' ;

energy_histo_axes = axes( 'OuterPosition', [ 0.1, 0.3, 0.8, 0.6 ]);

zoom( energy_histo_figure, 'yon' )

% % mini_map figure and axes initializations
minimap_figure = figure( 'Name'     , 'Vertex Curator: Reference Volume', ...
                         'Units'    , 'normalized',                       ...
                    'SizeChangedFcn',@minimap_figure_size_update,         ...
                           'Visible','off',                               ...
                         'OuterPosition' , [ half_screen, bottom_buffer + histo_heights + spacer, 2 * histo_widths + spacer, histo_heights ]);
                           
minimap_figure.Units = 'pixels' ; 
%  
% minimap_figure.Position([ 3, 4 ]) = min( minimap_figure.Position([ 3, 4 ]));

minimap_axes = axes(minimap_figure, 'OuterPosition', [0.2, 0.2, 0.7, 0.7]);
minimap_axes.Units = 'pixels';
view(minimap_axes,3);
    %% Variables                                
    
% initialize global variables
z_range                     = [ ];
z_limits                    = [ ];
index_image_2D              = [ ];
index2intensity             = [ ];
index_image_crop            = [ ];
truth_image_uint8           = [ ];
in_view_vertices_z          = [ ];
% gaussian_weight_in_z        = [ ];
argmax_vertex_image_2D      = [ ];
original_image_crop_2D      = [ ];
inverted_original_image_2D  = [ ];
intensity_limits_binarized  = true ;
counter_toggle              = 0 ;
counter_swipe_toggle        = [ ] ;
counter_threshold           = 0 ;
counter_add_vertex          = 0 ;
elapsed_time                = 0 ;

tstart = tic ;

vertex_structure_positions_linear_indexing = [ ];
   
size_of_minimap             = 100;
energy_threshold            = 0 ;

threshold_mat_3D            = energy_threshold.*ones(size_of_minimap,size_of_minimap,size_of_minimap);
xy_threshold_mat            = energy_threshold.*ones(size_of_minimap);
xz_threshold_mat            = energy_threshold.*ones(size_of_minimap);
yz_threshold_mat            = energy_threshold.*ones(size_of_minimap);
colormap_mat                = [];
threshold_color_range       = [0.5 1];
out_of_bounds_color         = 0.0;
border_color                = 0.25;
in_view_color               = 0.49;

% add a row for the background when indexing by the index image (that image has ones on background)
number_of_vertices = length( vertex_energies );   

    number_of_added_vertices =    0 ;

vertex_space_subscripts = [ zeros( max_number_of_added_vertices + 1, 3 ); vertex_space_subscripts ];
vertex_scale_subscripts = [  ones( max_number_of_added_vertices + 1, 1 ); vertex_scale_subscripts ];                      
vertex_energies         = [ zeros( max_number_of_added_vertices + 1, 1 ); vertex_energies         ];
                            
length_of_vertex_lists = number_of_vertices + max_number_of_added_vertices + 1 ;

% intialize logical vectors of class assignment
   true_vertices          = ones(  length_of_vertex_lists, 1, 'logical' );
in_view_vertices_xy       = ones(  length_of_vertex_lists, 1, 'logical' );
in_view_vertices          = ones(  length_of_vertex_lists, 1, 'logical' );
deleted_vertices          = zeros( length_of_vertex_lists, 1, 'logical' );
displayed_vertices        = zeros( length_of_vertex_lists, 1, 'logical' );

% mark the background and place holders as deleted and not displayed or in view
deleted_vertices(    vertex_space_subscripts( :, 1 ) <= 0 ) = true  ;
in_view_vertices(    vertex_space_subscripts( :, 1 ) <= 0 ) = false ;
in_view_vertices_xy( vertex_space_subscripts( :, 1 ) <= 0 ) = false ;

maximum_number_of_backups = 20 ;

backup_ordinate = 0 ;

% pull original from h5 files       
original_image = h52mat( path_to_original_data );

size_of_image = size( original_image );

% main display image initializtion
index_image    = ones( size_of_image, 'uint32' );

number_of_pixels_2D = size_of_image( 1 ) * size_of_image( 2 ) ;

% initialize axis limits
y_limits = [ 1, size_of_image( 1 )];
x_limits = [ 1, size_of_image( 2 )];

display_axes.XLim = x_limits ;
display_axes.YLim = y_limits ;

% initialize z_depth and z_thickness for initial viewing and slider limits
z_depth         = round( size_of_image( 3 ) / 2 );
z_thickness_max = size_of_image( 3 );

z_thickness = min( 5, z_thickness_max );         
    %% Templates and Images                     

% precalculating the volumes of influence of the vertices as linear indices into the 3D image for
% quick ease of use later during live performance.

number_of_scales = length( lumen_radius_in_microns_range );

scale_subscript_range = 1 : number_of_scales ;

structuring_element_linear_indexing_templates = cell( number_of_scales, 1 );

radii_in_pixels_range = lumen_radius_in_microns_range ./ microns_per_pixel ;

for scale_subscript = scale_subscript_range

    % find all pixel locations within the ellipsoid radii from the vertex position    
    structuring_element_linear_indexing_templates{ scale_subscript }                                            ...
        = int64( construct_structuring_element_V190( radii_in_pixels_range( scale_subscript, : ), size_of_image ));

end % constructing relative elements FOR scale

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
                                                                                          
z_depth_label         = uicontrol( display_figure, 'Style'   , 'text',                          ...
                                                   'Units'   , 'normalized',                    ...
                                                   'Position', [ 0.4125, 0.105, 0.165, 0.03 ],  ...
                                                   'String'  , 'Z-Depth',                       ...
                                                   'FontName', 'Times New Roman',               ...
                                                   'FontUnits', 'normalized',                   ...
                                                   'FontSize', 0.5                              );                                               
                                               
z_depth_slider        = uicontrol( display_figure, 'Style'   , 'slider',                        ...
                                                   'Min'     , 1,                               ...
                                                   'Max'     , size_of_image( 3 ),              ...
                                                 'SliderStep', [ 1 / size_of_image( 3 ),        ...
                                                                 5 / size_of_image( 3 )  ],     ...
                                                   'Value'   , 0,                               ...
                                                   'Units'   , 'normalized',                    ...
                                                   'Position', [ 0.305, 0.075, 0.385, 0.04],    ...
                                                   'Callback', @change_z_depth                  ); 

z_thickness_slider    = uicontrol( display_figure, 'Style'   , 'slider',                        ...
                                                   'Min'     , 0,                               ...
                                                   'Max'     , z_thickness_max,                 ...
                                                 'SliderStep', [ 1 / z_thickness_max,           ...
                                                                 5 / z_thickness_max  ],        ...                                                   
                                                   'Value'   , 0,                               ...
                                                   'Units'   , 'normalized',                    ...
                                                   'Position', [ 0.305, 0.0325, 0.385, 0.03],   ...
                                                   'Callback', @change_z_thickness              );

z_thickness_label     = uicontrol( display_figure, 'Style'   , 'text',                          ...
                                                   'Units'   , 'normalized',                    ...
                                                   'Position', [ 0.4125, 0, 0.165, 0.03 ],      ...
                                                   'String'  , 'Z-Thickness',                   ...
                                                   'FontName', 'Times New Roman',               ...
                                                   'FontUnits', 'normalized',                   ...
                                                   'FontSize', 0.5                              );                                               
					                                     
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
                                                   'Position', [ 0.825, 0.075, 0.125, 0.05],    ...
                                                   'String'  , 'Add',                           ...
                                                   'FontName', 'Times New Roman',               ...
                                                   'FontUnits', 'normalized',                   ...
                                                   'FontSize', 0.5,                             ...                                                                                                      
                                                   'Callback', @add_callback                    );

toggle_button         = uicontrol( display_figure, 'Style'   , 'pushbutton',                    ...
                                                   'String'  , 'Toggle',                        ...
                                                   'Units'   , 'normalized',                    ...
                                                   'FontName', 'Times New Roman',               ...
                                                   'FontUnits', 'normalized',                   ...
                                                   'FontSize', 0.5,                             ...                                                                                                      
                                                   'Position', [ 0.825, 0.0125, 0.125, 0.05],   ...
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
                                                
binarize_button           = uicontrol( energy_histo_figure,                                      ...
                                                    'Style'   , 'pushbutton',                    ...
                                                    'Units'   , 'normalized',                    ...
                                                    'Position', [ 0.4, 0.93, 0.2, 0.05 ],        ...
                                                    'String'  , 'Binary',                        ...
                                                    'FontName', 'Times New Roman',               ...
                                                    'FontUnits', 'normalized',                   ...
                                                    'FontSize', 0.5,                             ...                                                                                                      
                                                    'Callback', @binarize_callback               );
                                                
true_label                = uicontrol( energy_histo_figure,                                     ...
                                                    'Style'   , 'text',                         ...
                                                    'Units'   , 'normalized',                   ...
                                                    'Position', [ 0.1, 0.93, 0.2, 0.05 ],       ...
                                                    'FontName', 'Times New Roman',              ...
                                                    'FontUnits', 'normalized',                  ...
                                                    'FontSize', 1,                              ...         
                                             'BackgroundColor', [ 0, 1, 1 ],                    ...                                                                                                        
                                                    'String'  , 'True'                          );
                                                
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

if load_upon_startup, load_curation, 

else
    
    initialize_structuring_elements

    paint_index_image, update_vertex_intensities, update_z_zoom
    update_z_sliders, update_minimap, update_xy_zoom    

end



% turn all figures visible     

        display_figure.Visible = 'on';
        minimap_figure.Visible = 'on';
intensity_histo_figure.Visible = 'on';
   energy_histo_figure.Visible = 'on';
    
    %% Visual Updating Functions                    
    
   %% initialize structuring element                   
    function initialize_structuring_elements

        % precalculate linear indices from space subscripts in the image
        space_subscripts_int64 = int64( vertex_space_subscripts );

        vertex_position_linear_indices =   space_subscripts_int64( :, 1 )                                               ...
                                       + ( space_subscripts_int64( :, 2 ) - 1 ) * size_of_image( 1 )                    ...
                                       + ( space_subscripts_int64( :, 3 ) - 1 ) * size_of_image( 1 ) * size_of_image( 2 );

        vertex_position_linear_indices_cell = num2cell( vertex_position_linear_indices );

        structuring_element_linear_indexing = structuring_element_linear_indexing_templates( round( vertex_scale_subscripts ));

        vertex_structure_positions_linear_indexing = cellfun( @plus, structuring_element_linear_indexing,   ...
                                                                       vertex_position_linear_indices_cell, ...
                                                                                     'UniformOutput', false );   
          
    end

        %% overlay                                  
    function update_overlay

    %     % Make all truth images render at same intensity
    %     %note: this does not cancel background vector which should help
    %     %false positives shine though
    %     truth_image_uint8 = uint8(logical(truth_image_uint8).*128);

        % combine the raw and truth images into overlay
        overlay_image = truth_image_uint8 + inverted_original_image_2D ;

        % draw the image in the display axes
        image( display_axes, overlay_image )

        % reset axes properties and zoom position
        display_axes.XTick      = [ ];
        display_axes.YTick      = [ ];

        display_axes.XTickLabel = '' ;
        display_axes.YTickLabel = '' ;

        display_axes.XLim = x_limits ;
        display_axes.YLim = y_limits ;

    %     display_axes.ButtonDownFcn = @vertex_toggler2 ;

    end % update_overlay
        %% vector                                   
    function update_index_image_2D_crop

        index_image_crop = index_image( :, :, z_range );

        update_index_image_2D

    end % update_index_image_2D_crop

    function update_index_image_2D

        % create 2D index image from 3D index image by taking first non-background (background has index
        % 1) entry in third dimension note: permute matrix before to project in different dimension
        [ ~, argmax_vertex_image_2D ] = max( logical( index_image_crop - 1 ), [ ], 3 );

        % number_of_pixels_2D = numel( argmax );

        linear_index_of_max_image_2D = zeros( size( argmax_vertex_image_2D ));

        linear_index_of_max_image_2D( : ) = ( argmax_vertex_image_2D( : ) - 1 ) * number_of_pixels_2D + ( 1 : number_of_pixels_2D )';

        index_image_2D = index_image_crop( linear_index_of_max_image_2D );

        update_vertex_intensities

    end % update_index_image

    function update_vertex_intensities

        if ( intensity_limits_binarized )

            index2intensity = 128 * ones( size( vertex_energies ));
            
        else
            
            index2intensity = 128 * min( 1, (     energy_limits( 2 ) - vertex_energies   )  ...
                                              / ( energy_limits( 2 ) - energy_limits( 1 )));

        end % IF binarized intensity

        update_truth_image

    end % update_truth_image_intensity

    function update_truth_image

        truth_image_uint8 =  uint8( index2intensity( index_image_2D )                     ...
                                    .* double( displayed_vertices( index_image_2D ) & ~ deleted_vertices( index_image_2D )) ...
                                    .* cat( 3,         double( ~ true_vertices( index_image_2D ))  ...
                                               + 0.3 * double(   true_vertices( index_image_2D )), ...
                                                 0.3 * double( ~ true_vertices( index_image_2D ))  ...
                                               + 0.7 * double(   true_vertices( index_image_2D )), ...
                                                 0.3 * double( ~ true_vertices( index_image_2D ))  ...
                                               +       double(   true_vertices( index_image_2D ))  ));

    end % update_truth_image
        %% raw                                      
    function update_original_image

        % crop 3D volume
        original_image_crop = original_image( :, :, z_range );

        % take a cube of the 3D image, max project in the third dimension and if requested, 
        [ original_image_crop_2D, argmax_original_image_2D ] = max( original_image_crop, [ ], 3 );

%     %     if gaussian_depth
% 
%         % do Guassian weighting in the third dimension based on the z location of the argument of
%         % the max projection
%         minimum_intensity_in_crop = double( min( original_image_crop_2D( : )));
% 
%         original_image_crop_2D = double(   original_image_crop_2D                  ...
%                                          - minimum_intensity_in_crop )             ...
%                                .* gaussian_weight_in_z( argmax_original_image_2D ) ...
%                                + minimum_intensity_in_crop ;
% 
%     %     end % IF gaussian_depth

        original_image_crop_2D = double( original_image_crop_2D );

        update_original_image_intensity

    end % update_original_image

    function update_original_image_intensity

        % adjust contrast for the 2D original image and 
        % restrict intensity values to the integers 0 : 255
        inverted_original_image_2D                                                    ...
            = uint8( 128                                                              ...
                     - min( 1, (     original_image_crop_2D - intensity_limits( 1 ))  ...
                                 / ( intensity_limits( 2 )  - intensity_limits( 1 ))) ...
                     * 128                                                            );  

        intensity_histo_axes.XLim = intensity_limits([ 1, 2 ]);                 

    end % update_original_image_contrast
        %% histogram                                
    function update_intesnity_histogram

        y_range =   max(      1            , round( y_limits( 1 ))) ...
                  : min( size_of_image( 1 ), round( y_limits( 2 ))) ;

        x_range =   max(      1            , round( x_limits( 1 ))) ...
                  : min( size_of_image( 2 ), round( x_limits( 2 ))) ;

        intensity_histo_axes_limits([ 3, 4 ]) = intensity_histo_axes.YLim ;

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

        histogram( vertex_energies( in_view_vertices & displayed_vertices & ~ deleted_vertices ), 'Parent', energy_histo_axes, 'FaceColor', [0, 0, 0], 'EdgeColor', 'white' );
        
        colormap( energy_histo_axes, colormap_mat );
        
        %COMMENT OUT TO SHOW ALL THRESHOLDS ON HISTOGRAM
        adj_y_lim = 1 + [round((size_of_minimap-1)*y_limits(1)/size_of_image(1)) round((size_of_minimap-1)*y_limits(2)/size_of_image(1))];    
        adj_x_lim = 1 + [round((size_of_minimap-1)*x_limits(1)/size_of_image(2)) round((size_of_minimap-1)*x_limits(2)/size_of_image(2))];
        adj_z_lim = 1 + [round((size_of_minimap-1)*z_limits(1)/size_of_image(3)) round((size_of_minimap-1)*z_limits(2)/size_of_image(3))];

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

        %% mini-map visual updates   

    function update_threshold_visualization_matricies

        adj_y_lim = 1 + [round((size_of_minimap-1)*y_limits(1)/size_of_image(1)) round((size_of_minimap-1)*y_limits(2)/size_of_image(1))];    
        adj_x_lim = 1 + [round((size_of_minimap-1)*x_limits(1)/size_of_image(2)) round((size_of_minimap-1)*x_limits(2)/size_of_image(2))];
        adj_z_lim = 1 + [round((size_of_minimap-1)*z_limits(1)/size_of_image(3)) round((size_of_minimap-1)*z_limits(2)/size_of_image(3))];
        
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
        adj_y_lim = 1 + [round(99*y_limits(1)/size_of_image(1)) round(99*y_limits(2)/size_of_image(1))];    
        adj_x_lim = 1 + [round(99*x_limits(1)/size_of_image(2)) round(99*x_limits(2)/size_of_image(2))];
        adj_z_lim = 1 + [round(99*z_limits(1)/size_of_image(3)) round(99*z_limits(2)/size_of_image(3))];

        adj_z_lim( 1 ) = max( adj_z_lim( 1 ) - 1, 1   );
        adj_z_lim( 2 ) = min( adj_z_lim( 2 ) + 1, size_of_minimap );    

        %color selection tool
        %relocated to load_colormap function
        
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
%        view_mini(find(isnan(view_mini))) = [];
        c_mini(1,:) = border_color; c_mini(end-1,:) = border_color; c_mini(:,1) = border_color; c_mini(:,end-1) = border_color;
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
        c_mini(1,:) = border_color; c_mini(end-1,:) = border_color; c_mini(:,1) = border_color; c_mini(:,end-1) = border_color;
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
        view_mini(max(1,adj_z_lim(1)):min(size_of_minimap,adj_z_lim(2)),max(1,adj_x_lim(1)):min(size_of_minimap,adj_x_lim(2))) = in_view_color;
        c_mini(1,:) = border_color; c_mini(end-1,:) = border_color; c_mini(:,1) = border_color; c_mini(:,end-1) = border_color;
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

        set( x_start_label, 'String', num2str( round( x_limits( 1 ))));
        set( y_start_label, 'String', num2str( round( y_limits( 1 ))));
        set( z_start_label, 'String', num2str(   max( z_limits( 1 ), 1 )));

        set( x_end_label, 'String', num2str( round( x_limits( 2 ))));
        set( y_end_label, 'String', num2str( round( y_limits( 2 ))));
        set( z_end_label, 'String', num2str(   min( z_limits( 2 ), size_of_image( 3 ))));

        set( x_total_label, 'String', num2str( round( x_limits( 2 )) - round( x_limits( 1 )) + 1 ));
        set( y_total_label, 'String', num2str( round( y_limits( 2 )) - round( y_limits( 1 )) + 1 ));
        set( z_total_label, 'String', num2str(   min( z_limits( 2 ), size_of_image( 3 )) ...
                                               - max( z_limits( 1 ), 1 ) + 1             ));

    end % update_minimap

        %% z navigation                             
    function update_z_zoom
        % update main display

        % compute limits from the z depth and thickness settings
        z_limits = [ z_depth - z_thickness, ...
                     z_depth + z_thickness  ];    

        z_range = max(      1            , z_limits( 1 )) ...
                : min( size_of_image( 3 ), z_limits( 2 )) ;

        % find vertices in the view in z
        in_view_vertices_z = vertex_space_subscripts( :, 3 ) >= z_range(  1  ) ... 
                           & vertex_space_subscripts( :, 3 ) <= z_range( end ) ;

        in_view_vertices = in_view_vertices_xy & in_view_vertices_z ;    

%         % create gaussian for psuedo focus weighting    
%         gaussian_weight_in_z = exp( - ( z_range - z_depth ) .^ 2 / z_thickness ^ 2 / 4 );
% 
%         if z_thickness == 0, gaussian_weight_in_z( 1 ) = 1 ; end

        update_index_image_2D_crop, update_original_image, update_overlay

        update_intesnity_histogram, update_energy_histogram, 

        update_minimap

    end % update_displays

    function update_z_sliders

        set( z_depth_slider,     'Value', z_depth     )
        set( z_thickness_slider, 'Value', z_thickness )

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
%         x_limits = [display_axes.XLim(1) display_axes.XLim(1) + display_length];
%         y_limits = [display_axes.YLim(1) display_axes.YLim(1) + display_length];
% 
%         display_axes.XLim( 2 ) = x_limits( 2 );
%         display_axes.YLim( 2 ) = y_limits( 2 );

        x_limits = [ display_axes.XLim( 1 ), display_axes.XLim( 2 )];
        y_limits = [ display_axes.YLim( 1 ), display_axes.YLim( 2 )];
        
        ylim([ 1, size_of_image( 1 )])
        xlim([ 1, size_of_image( 2 )])
        
        zoom( 'reset' )
        
        ylim( y_limits )
        xlim( x_limits )

        in_view_vertices_xy = vertex_space_subscripts( :, 1 ) >= y_limits( 1 ) ...
                            & vertex_space_subscripts( :, 1 ) <= y_limits( 2 ) ... 
                            & vertex_space_subscripts( :, 2 ) >= x_limits( 1 ) ... 
                            & vertex_space_subscripts( :, 2 ) <= x_limits( 2 ) ;

        in_view_vertices = in_view_vertices_xy & in_view_vertices_z ;

        % propogate the new display to the histograms and mini-map
        update_intesnity_histogram, update_energy_histogram, update_minimap, main_figure_size_update

    %     move_view_box_in_mini_map

    end % update_xy_zoom

    function change_z_depth( source, ~ )

        % update the z depth from the slide value, then update the crop
        z_depth = round( source.Value ); 

        update_z_zoom

    end % change_z_depth

    function change_z_thickness( source, ~ )

        % update the z depth from the slide value
        z_thickness = round( source.Value ); 

        update_z_zoom

    end % change_z_thickness

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
        
        

        image_aspect_x_over_y = ( x_limits( 2 ) - x_limits( 1 ) + 1 ) ...
                              / ( y_limits( 2 ) - y_limits( 1 ) + 1 );

        display_figure_position_min = min( display_figure_position( 3 ) / image_aspect_x_over_y, ...
                                           display_figure_position( 4 ));
        
        display_axes.Position = [ display_figure_position( 3 ) * 0.075, ...
                                  display_figure_position( 4 ) * 0.15,  ...
                                  display_figure_position_min * 0.85 * image_aspect_x_over_y, ...
                                  display_figure_position_min * 0.85 ];

    end

    function minimap_figure_size_update( ~, ~) 
        pos = get(minimap_figure,'Position');
        minimap_axes.OuterPosition = [pos(3)*0.2, pos(4)*0.2, min(pos(3)*0.7,pos(4)*0.7), min(pos(3)*0.7,pos(4)*0.7)];
    end

        %% automatic backup operation               
    function backup_curation

        save([ path_to_saved_curation, '_backup_', num2str( backup_ordinate )],        'true_vertices', ...
                                                                                    'deleted_vertices', ...
                                                                                  'displayed_vertices', ...
                                                                                            'x_limits', ...
                                                                                            'y_limits', ...
                                                                                             'z_depth', ...
                                                                                         'z_thickness', ...
                                                                             'vertex_space_subscripts', ...
                                                                             'vertex_scale_subscripts', ...
                                                                             'vertex_energies'        , ...
                                                                                    'intensity_limits', ...
                                                                                       'energy_limits', ...
                                                                            'number_of_added_vertices', ...
                                                                                    'energy_threshold', ...
                                                                                    'xy_threshold_mat', ...
                                                                                    'yz_threshold_mat', ...
                                                                                    'xz_threshold_mat', ...
                                                                                    'threshold_mat_3D', ...
                                                                                'counter_add_vertex',   ...
                                                                                'counter_toggle'    ,   ...
                                                                                'counter_swipe_toggle', ...
                                                                                'counter_threshold'     );

        backup_ordinate = mod( backup_ordinate + 1, maximum_number_of_backups );                                                                                  
        
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

            answer = questdlg( 'Overwrite the previously saved curation?', 'Vertex Curator' );

            switch answer, case { 'No', 'Cancel', '' }, return, end

        end
        
        elapsed_time = elapsed_time + toc( tstart );
        
        tstart = tic ;

        save( path_to_saved_curation,      'true_vertices', ...
                                        'deleted_vertices', ...
                                      'displayed_vertices', ...
                                 'vertex_space_subscripts', ...
                                 'vertex_scale_subscripts', ...
                                 'vertex_energies'        , ...
                           ...          'path_to_energy_data', ...
                           ...        'path_to_original_data', ...
                                        'intensity_limits', ...
                                           'energy_limits', ...
                                'number_of_added_vertices', ...
                                         'backup_ordinate', ...
                                                'x_limits', ...
                                                'y_limits', ...
                                                 'z_depth', ...
                                             'z_thickness', ...
                                        'energy_threshold', ...
                                        'xy_threshold_mat', ...
                                        'yz_threshold_mat', ...
                                        'xz_threshold_mat', ...
                            'max_number_of_added_vertices', ...
                                        'threshold_mat_3D', ...
                                    'counter_add_vertex'  , ...
                                    'counter_toggle'      , ...
                                    'counter_swipe_toggle', ...
                                    'counter_threshold',    ...
                                    'elapsed_time'          );

    end % store_curation

    function load_curation( ~, ~ )

        backup_curation

        load( path_to_saved_curation )

        display_axes.XLim = x_limits ;
        display_axes.YLim = y_limits ;

        initialize_structuring_elements
        
        refresh_and_paint_index_image
        
        tstart = tic ;        

    end % load_curation

        %% sweeping and painting operations         
    function sweep_index_image

        % move any in view vertices, displayed false, and not yet deleted to the deleted category
        to_be_deleted_vertices = ~ true_vertices & in_view_vertices & displayed_vertices & ~ deleted_vertices ;

        % If ever once displayed, the vertex will remain in the displayed vector
%             displayed_vertices = ~ to_be_deleted_vertices & displayed_vertices ; 
              deleted_vertices =   to_be_deleted_vertices |   deleted_vertices ;

        % background is coded by vertex of index 1, set it to not deleted to avoid searching and
        % deleting the background (and whatever bug that might create) in the following
        to_be_deleted_vertices( 1 ) = false ;

        to_be_deleted_indices = find( to_be_deleted_vertices )' ;

        for vertex_index = to_be_deleted_indices

            % set all the indices belonging to deleted ones to background (background is one)
            index_image( vertex_structure_positions_linear_indexing{ vertex_index }) = 1 ;

        end % FOR deleted vertex indices
    end % sweep_index_image

    function paint_index_image

        % attempt to paint any vertices in the current view that aren't displayed already or deleted
        to_be_painted_vertices = find( ~ ( deleted_vertices | displayed_vertices ) & in_view_vertices )';

        % subtract off one so that background is coded by zero for the following FOR loop
        index_image = index_image - 1 ;

        % loop through the vertices from best to worst (lowest to highest energy). Should be sorted before
        % this function.
        for vertex_index = to_be_painted_vertices 

            % check if this area has already been painted
            if all( ~ index_image( vertex_structure_positions_linear_indexing{ vertex_index }))

                displayed_vertices( vertex_index ) = true ;

                % paint the image to mark this newly displayed vertex        
                index_image( vertex_structure_positions_linear_indexing{ vertex_index }) = vertex_index - 1 ;

            end % IF area is blank canvas (no conflicts with any higher energy vertices)
        end % vertex FOR

        index_image = index_image + 1 ;

    end % volume_exclude_vertices

    function refresh_and_paint_index_image

        % select the vertices to attempt to paint in the following FOR loop
        to_be_painted_vertices = find( displayed_vertices & ~ deleted_vertices )';

        % background is temporarily coded by zero for the following FOR loop
        index_image = zeros( size_of_image, 'uint32' );

        % loop through the vertices from best to worst (lowest to highest energy). Should be sorted before
        % this function.
        for vertex_index = to_be_painted_vertices 

            % check if this area has already been painted
            if all( ~ index_image( vertex_structure_positions_linear_indexing{ vertex_index }))

                % paint the image to mark this newly displayed vertex        
                index_image( vertex_structure_positions_linear_indexing{ vertex_index }) = vertex_index - 1 ;

            else % volume conflicts with at least one higher energy vertex
                
                displayed_vertices( vertex_index ) = false ;
                
            end % IF area is blank canvas (no conflicts with any higher energy vertices)
        end % vertex FOR

        index_image = index_image + 1 ;

        update_z_zoom, update_xy_zoom, update_z_sliders
        update_vertex_intensities, update_energy_histogram, update_minimap

    end % paint_from_scratch

    function paint_new_vertex( to_be_painted_index )

       space_subscript_int64 = int64( vertex_space_subscripts( to_be_painted_index, : ));

       vertex_position_linear_index =  space_subscript_int64( :, 1 )                                               ...
                                   + ( space_subscript_int64( :, 2 ) - 1 ) * size_of_image( 1 )                    ...
                                   + ( space_subscript_int64( :, 3 ) - 1 ) * size_of_image( 1 ) * size_of_image( 2 );

        vertex_structure_positions_linear_indexing{ to_be_painted_index }                                    ...
            = structuring_element_linear_indexing_templates{ round( vertex_scale_subscripts( to_be_painted_index ))} ...
                                                                              + vertex_position_linear_index ;

        % look at the volume to be painted and see if there are vertices in the way to delete
        to_be_deleted_indices = unique( index_image( vertex_structure_positions_linear_indexing{ to_be_painted_index }))';

        % don't try to delete backgrouund
        to_be_deleted_indices( to_be_deleted_indices == 1 ) = [ ];

        for vertex_index = to_be_deleted_indices

            % set all the indices belonging to deleted ones to background (background is one)
            index_image( vertex_structure_positions_linear_indexing{ vertex_index }) = 1 ;

        end % FOR deleted vertex indices

        index_image( vertex_structure_positions_linear_indexing{ to_be_painted_index }) = to_be_painted_index ;

    end % add_vertex
        %% sweeping and painting                    
    function sweep_callback( ~, ~ )

        backup_curation

        % Find the vertices that are false, delete them, don't display them, and turn them into
        % background in the 3D index image.
        sweep_index_image

        update_index_image_2D_crop, update_overlay, update_energy_histogram

        % re-initialize energy histogram and energy_limits for visualization and thresholding

    end % sweep_falses

    function paint_callback( ~, ~ )

        backup_curation

        % repopulate the current field of view with any vertices that have not yet been deleted and are
        % not currently displayed.  Start with the best contrast objects and go down the list, trying to
        % place new objects where the volumes don't conflict
        paint_index_image

        update_index_image_2D_crop, update_overlay, update_energy_histogram

    %     add_repop_to_mini_map

    end
        %% toggling and placing                     
    function toggle_callback( ~, ~ ) 

        toggle_button.Enable = 'off';

        % prompt user to click the image
%         [ x_click, y_click ] = ginput( 1 );
        rect_click = getrect( display_axes );
        
        x_click = rect_click( 1 );            
        y_click = rect_click( 2 );
        
        x_min = rect_click( 1 );
        y_min = rect_click( 2 );
        x_max = rect_click( 1 ) + rect_click( 3 );
        y_max = rect_click( 2 ) + rect_click( 4 );
        x_ave = rect_click( 1 ) + rect_click( 3 ) / 2 ;
        y_ave = rect_click( 2 ) + rect_click( 4 ) / 2 ;

        % IF center of rectangle is inside the field of view  
        if    y_ave > y_limits( 1 ) && y_ave < y_limits( 2 ) ...
           && x_ave > x_limits( 1 ) && x_ave < x_limits( 2 )
       
            backup_curation
            
            is_toggling_threshold = 2 ;       
                       
            is_toggling = rect_click( 3 ) ^ 2 + rect_click( 4 ) ^ 2 < is_toggling_threshold ^ 2 ;
       
            if is_toggling

%                 % This will fail if user clicks outside the image. Exit the recursive function in that case.
%                 try vertex_index = index_image_2D( round( y_click ), round( x_click )); catch, return, end

                vertex_indices = index_image_2D( round( y_click ), round( x_click ));

                counter_toggle = counter_toggle + 1 ;

            else % is_swipe_toggling

                in_swipe_vertices_xy = vertex_space_subscripts( :, 1 ) >= y_min ...
                                     & vertex_space_subscripts( :, 1 ) <= y_max ...
                                     & vertex_space_subscripts( :, 2 ) >= x_min ...
                                     & vertex_space_subscripts( :, 2 ) <= x_max ;

                % logical vector, not truly vertex_indices ( = find( vertex_indices ))
                vertex_indices = in_swipe_vertices_xy & in_view_vertices_z ;

                counter_swipe_toggle( end + 1 ) = sum( vertex_indices );

            end % is_toggling
            
            % Toggle the truth of the selected vertex.      
            true_vertices( vertex_indices ) = ~ true_vertices( vertex_indices );
            
            update_truth_image, update_overlay            
            
            % recall this function to select anew
            toggle_callback ;

        else % ELSE center of rectangle is outside the field of view

            toggle_button.Enable = 'on';
            return % Exit the recursive function

        end % IF center of rectangle is inside the field of view    

    end % vertex_toggler

    function add_callback( ~, ~ )

        % prompt user to click the image
        [ x_click, y_click ] = ginput( 1 );

        % Get the vertex index from the position of the mouse click. 
        if    y_click > y_limits( 1 ) && y_click < y_limits( 2 ) ...
           && x_click > x_limits( 1 ) && x_click < x_limits( 2 )

            z_start = max( z_limits( 1 ), 1 ) ;
            z_end   = min( z_limits( 2 ), size_of_image( 3 ));
            z_count = z_end - z_start + 1 ;

            added_vertex_xy_subscripts = [ round( y_click ), round( x_click )];

            vertex_energy_and_scale_index = h52mat(                        path_to_energy_data, ...
                                                    [ added_vertex_xy_subscripts, z_start, 1 ], ...
                                                                          [ 1, 1, z_count, 2 ]  ); 

            [ added_vertex_energy, added_vertex_z_relative ]                                            ...
                                                      = min( vertex_energy_and_scale_index( z_count + 1 : end ));

            added_vertex_scale_index = vertex_energy_and_scale_index( added_vertex_z_relative );

            added_vertex_z_subscript = added_vertex_z_relative + z_start - 1 ;
            
            to_be_painted_subscripts = uint16([ added_vertex_xy_subscripts, ...
                                                added_vertex_z_subscript    ]);

            to_be_painted_index = number_of_added_vertices + 2 ;

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

                vertex_energies(         to_be_painted_index    ) = added_vertex_energy      ;
                vertex_scale_subscripts( to_be_painted_index    ) = added_vertex_scale_index ;
                vertex_space_subscripts( to_be_painted_index, : ) = to_be_painted_subscripts ;                 
                
                paint_new_vertex( to_be_painted_index )

                  deleted_vertices( to_be_painted_index ) = false ;
                displayed_vertices( to_be_painted_index ) = true ;
                
                counter_add_vertex = counter_add_vertex + 1 ;

                update_vertex_intensities

                number_of_added_vertices = number_of_added_vertices + 1 ;          

                update_index_image_2D_crop, update_overlay, update_energy_histogram
                
                % recall this function to select anew
                add_callback                        

            end % vertex would cross the image boundary
        end % IF click inside current field of view
    end % vertex_placer
        %% intensity histogram                      

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
        
        intensity_limits_binarized = ~intensity_limits_binarized;
        
        if intensity_limits_binarized
            
            binarize_button.String = 'Binary' ;
            
        else
            
            binarize_button.String = 'Graded' ;
            
        end
        
        update_vertex_intensities, update_overlay
        
    end %binarize_callback	
        
    function change_energy_min( source, ~ )

        energy_limits( 1 ) = str2num( source.String );

        if energy_limits( 1 ) >= energy_limits( 2 ), warning('The energy min should be below the max'), end

        update_vertex_intensities, update_overlay, update_energy_histogram

    end % change_intensity_limits

    function change_energy_max( source, ~ )

        energy_limits( 2 ) = str2num( source.String );

        if energy_limits( 1 ) >= energy_limits( 2 ), warning('The energy max should be above the min'), end    

        update_vertex_intensities, update_overlay, update_energy_histogram

    end % change_intensity_limits

    function apply_energy_threshold( source, ~ )                     

        backup_curation

        energy_threshold = max( str2num( source.String ), min( vertex_energies ));

        true_vertices( in_view_vertices ) = vertex_energies( in_view_vertices ) ...
                                          < energy_threshold                    ;

        update_threshold_visualization_matricies, update_minimap
        update_truth_image, update_overlay, update_energy_histogram

        counter_threshold = counter_threshold + 1 ;

        if energy_threshold > 0 || imag( energy_threshold ) ~= 0, warning( 'energy threshold should be non-positive real number' ), end

        if energy_limits( 2 ) <= energy_threshold, warning('The energy threshold should be below the max'), end

        if energy_limits( 1 ) >= energy_threshold, warning('The energy threshold should be above the min'), end

    end % apply_energy_threshold

%% Finalizations                                

% wait for user to play on the curator
uiwait( display_figure )

store_curation

chosen_vertices = true_vertices & displayed_vertices & ~ deleted_vertices ;

% background should not be chosen
chosen_vertices( 1 ) = 0 ;

% performing the selections proposed by either choose_vertices or vertex_curator
vertex_space_subscripts = vertex_space_subscripts( chosen_vertices, : );
vertex_scale_subscripts = vertex_scale_subscripts( chosen_vertices    );
vertex_energies         =         vertex_energies( chosen_vertices    );

close([ minimap_figure energy_histo_figure intensity_histo_figure ]);

end % main FUNCTION