%% Vectorize                                                                                        
% SAM 12/12/17
%
% This Script is the master file of the vectorization algorithm.  The point of the vectorization
% algorithm is to convert a raw 3D image of neural vasculature to a vectorized model of the vascular
% components.  This vectorized model can be rendered as a 3-dimensional image at any requested
% resolution, or it could be analyzed for statistical properties such as volume fraction of vessel.
% This vectorzied model could then be registered with other models that were formed from a set of
% images in a tiled pattern to obtain a larger field of view.  Alternatively (or additionally) the
% other vectorzied models could come from the same field of view as the original but offest in time
% to attain a movie of the same field of view. Statistics from each snapshot could then be tracked
% over time.
%
%

%% ------------------------------------------  Inputs --------------------------------------------- 
%
% The main input is the location of the folders of raw 2D images of fluorscence intensities.  The
% fluorescnce signal should originate from a fluorophore that is localized to the blood vessel
% lumens.  The folders of data should be arranged into the following directory structure:
%
%   root_directory
%     ->
%       [ ROI_base_name, ROI_names{ ROI_index( 1 )}]
%       [ ROI_base_name, ROI_names{ ROI_index( 2 )}]
%                         ...
%       [ ROI_base_name, ROI_names{ ROI_index( end )}]
%
% Where ROOT_DIRECTORY is a string of the full path to the location of the folder that contains all
% of the folders from each region of interest from the imaging session.  Inside of each ROI folder
% should be the raw pages of data read at each z-location from the scanning mirror system.
%
% Other inputs are the parameter settings which can be adjusted in this script:
%
% Note: this input scheme is compatible with the data organization used in Andy Dunn's FOI
% Laboratory in August and September of 2017
%
%% ------------------------------------------ Outputs --------------------------------------------- 
%
% The output is the vectorized model of the vascular network depicted in the 3D input image.  The
% vectorized model is composed of a list of objects that we can think of as 5 dimensional vectors
% (3D position, 1D size, and a contrast metric). These objects when plotted as volume filled spheres
% are the 3-dimensional volume model of the vasculature depicted in the input image.  
% 
% 
%% ------------------------------------------ Methods --------------------------------------------- 
%
% Vectorization is accomplished in four steps:  (1) centerline enhancement, (2) vertex detection,
% (3) edge formation, and (4) strand and junction assignment.
% 
% 1) Multi-scale gradient and curvature information from the original 3D image is combined to form a
% 4-dimensional, multi-scale, centerline-enhanced, 3D image (4DI).
% 
% 2) Spherical objects (vertices) are detected as local maxima in the 4DI (with associated x, y, z,
% radius, and 4DI value, similar method to the first part of the SIFT algorithm (David Lowe,
% International Journal of Computer Vision, 2004)).
% 
% 3) Edges are formed between vertices by voxel to voxel random walks through the 4DI image that
% favor higher 4DI values along the path.
% 
% 4) Strands (single random color in the colored strands image) are defined as the strongly
% connected components in a directed graph that includes the best two edges from each vertex
% (ordered by 4DI value). Junctions (labeled white in the colored strands image) are defined as the
% additional components that are introduced when vertices are allowed three edges apiece.

%% version notes:                                                                                   
%
% V02: in which the directory structure is overhauled and the calling of functions, and the saving
% of files and settings are all standardized.
%
% V10: in which:
%
% 1) there is no downsampling (all calculations are done at the same resolution).  This allows easy
% laplacian comparisons across octaves, easy zeroing of spaces across all octaves, and easy
% traversal through the 4D when placing edges.  This comes at the cost of increased storage and
% time.  Some of this storage will be won back by not needing to overlap the octaves at their edges
% in scale space.  More storage can be won back by eliminating the need to interpolate.  This would
% require smart sampling of the gaussian kernels (perhaps in the Fourier space after analytic
% transformation) and proper anisotropic weighting of the direction (and symmetry) feature(s) based
% on their orientations.  This increase in storage and time is not a primary concern.
%
% 2) without the octave breaks, we can compare all the laplacians on even playing field and
% eliminate vertices that lie within a radius of another vertex whose laplacian is more negative.
%
% 3) the edges are placed in the following way:  for each vertex, and each direction, a neighborhood
% of laplacian image is loaded around that vertex. The sphere representing the vertex itself is
% zero'd out, except for the cone representing the current direction.  A 26 neighbor strel that
% spans all the scales is placed at the center of the vertex and the lowest laplacian valued pixel
% is selected as the first mini_vertex of the edge.  The sterel is then centered at that mini-vertex
% and the previous strel is zero'd out.  The next however many mini-vertices in the edge are
% selected by succesively choosing the lowest laplacian from the partially zerod strels.  If at any
% point, there are no negative laplacians to choose from, then this edge is terminated.  If at any
% point, this edge finds itself at the same pixel and scale as another vertex, then that vertex is
% the terminal vertex of this edge.
%
% 4) laying the infrastructure to input vessel sizes in microns instead of in pixels and to
% incorporate deconvolution of the PSF during the blurring stages.
%
% 5) fixing a bug that put the image size in the settings folder, which would be incompatible with
% the intended use of the settings folder to apply to multiple ROI's equally.  Instead, we can read
% the image size from the h5 file like this:
%
% original_file = [ directories{ 4, ROI_index }, interpolate_handle ]; 
% 
% original_info = h5info( original_file );
% 
% size_of_image = uint16( original_info.Datasets.Dataspace.Size );
%
% V12 in which phase II is removed.  also fixed a naming error in the blurring settings
% "size_per_sigma" -> "sigma_per_size" SAM 4/22/18
%
% V13 in which the Laplacian is not necessarilly calculated, but an enery field is created to lie
% over the 4D image whose exact formula is an input variable but generally depends on the curvature
% values which are fleetingly calculated in the get_energy function.  The output of the get energy
% function is the energy field (which could be the laplacian) and is fed to the get_vertices and
% get_edges functions. 
%
% Also removing interpolation as an workflow option.  So the image size is retrieved from the
% original file hereafter and the voxel aspect ratio is declared in the blurring step. 
%
% note, we should make the energy function two 3D images instead of a 4D image.  We should min
% project the energy function across the size dimension but remember the scale coordinate where the
% min occurs. Both the vertices and edges end up doing a min projection anyway, so doing it up front
% would save computation a bit and a lot of storage.  Will have to adapt the get_vertices and
% get_edges functions along with get_energy to enact this change. (not yet implemented 5/5/18 SAM).
% 
% SAM 4/25/18 
%
% V14 in wich the derivatives are computed in the fourier domain instead of by finite differencing
% in the voxel array space. Blurring and energy calculation are thus merged into one step.  SAM
% 5/5/18
%
% V15 in which the energy is min projected across the scale dimension before saving into h5 file SAM
% 5/7/18
%
% V16 in which the clean up phase is added.  We attempt to assign existence probabilities to the
% different model components so that we may easily threshold them later.  SAM 5/11/18
%
% V162 in which the transition matrix is just multiplied and there is no connection test or true
% eigenvalue solving SAM 5/16/18
%
% V161 in which the transitino matrix is decomposed into connected components and eigenvalues are
% found for each connected component SAM 5/16/18
%
% V170 in which only the top two edges are kept from each vertex in the transition matrix and then
% the strongly connected components of that transition matrix (viewed as an adjacency matrix) are
% called strands. SAM 5/22/18
%
% V180, Vertices are chosen such that no vertex center is inside of a vertex of lower energy. SAM
% 5/29/18  Also, The edges are chosen such that no edge will conflict with a beginning or ending
% vertex of a lower energy edge. SAM 5/31/18 (retro-dated from 7/11/18 following the choose edges
% version notes)
%
% V190: SAM 7/18/18
%
% Manual curation is added (to the vertices section for starters).
%
% vector positions and sizes are smoothed along the axis of the strand that owns them (amount of
% smoothing is an input).
%
% Bifurcation Vertices are added to the visual ouptuts of the network section
%
% A 3D vector field of blood flow direction is estimated by the spatial derivatives with respect to
% the strand index (this vector is painted inside of the volume of the associated spherical object
% and edges that are better contrast overwrite edges of lesser contrast).  
% 
% A bug was fixed in the function get_network (previously our definition of strand objects as being
% a series of vectors was not achieved by the code.  Some strands contained bifurcations before.)
%
% The input, output, and methods sections at the top of this script were updated.  (The methods
% section could be expanded and the outputs section is incomplete as it does not include edges,
% strands, junctions, bifurcations, or flow field direction.
%
% To-Do:  
%
% add manual curation to the edges and network sections as needed.  
%
% Edge curation could rate every edge by its rank in the following way: lower score is better: every
% edge is assigned a number of points equal to the ordinate of the choice of this edge with respect
% to either of its vertices added together.  The lowest/best score would thus be a 2, meaning that
% both vertices choose this edge first.  If the edge is chosen as first by one vertex and second by
% the other, it would get a score of three.  The user could manipulate the edges under this metric,
% for instance, they could set all edges with scores 4 or greater to false.
%
% Or, edges could be rated by ( max_energy - min_energy ), and lower values would be favored.
%
% For the creation of the flow field make two improvements: (1) input the contrast metrics on a per
% vector basis instead of a per edge basis to make the resulting image more accurate.  (2) use the
% vessel direction to make a cylindrical object (not a spherical object) for the filling of the flow
% field directions.
%
% SAM 7/18/18

%% 0: Requested workflows, utilizing previous data, data deletion, and visual outputs               

%% ------------------------------------- Workflow Requests ---------------------------------------- 
% Each workflow has its own inputs and those are stored in individual settings files.  The workflow
% file remembers which settings were used for which runs

ordinate_of_start_workflow = 5 ;
ordinate_of_final_workflow = 5 ;

% ------------------------------------ Workflow Definitions ---------------------------------------
% these are ordered.  Down the list is downstream in workflow

                                % ordinate
workflows = { 'process_raw' ... % 1
              'energy'      ... % 2 
              'vertices'    ... % 3 
              'edges'       ... % 4
              'network'     };  % 5
              
number_of_workflows = length( workflows );

% create the PRODUCTIVE variable to encode the scheduled workflows
productive = ( 1 : number_of_workflows )' >= ordinate_of_start_workflow ...
           & ( 1 : number_of_workflows )' <= ordinate_of_final_workflow ;

% for example, 
%
%                          % ordinate
% % productive = [ 0 ; ... % 1
% %                1 ; ... % ordinate_of_start_workflow
% %                1 ; ... % 3
% %                1 ; ... % ordinate_of_final_workflow
% %                0 ; ];  % 5

                       % ordinate
forgetful  = [ 0 ; ... % 1
               0 ; ... % 2
               0 ; ... % 3
               0 ; ... % 4
               0   ];  % 5
           
visual     = [ 0 ; ... % 1
               0 ; ... % 2
               0 ; ... % 3
               0 ; ... % 4
               1   ];  % 5 
           
curation   = { '' ;      ... % 1
               '' ;      ... % 2
              'none' ; ... % 3
              'none' ; ... % 4
               ''        };  % 5           
%% --------------------------------------- Data Handling ------------------------------------------ 

% this is the folder that keeps the ROI folders of the raw imaging data
root_directory = 'E:\2P imaging\20170315_MouseT326_session2_fused_crop\';

% ROI_base_name = 'Compiled' ;
% ROI_suffixes = { 'Tile' };

% 12/13/18
ROI_base_name = 'Fused_T123_T456' ;
ROI_suffixes = { '_Raw' };

% Batches are groups of ROI's that get fed through workflows as a group.  This organization is
% useful for comparisons between ROI's with identical workflow parameters.  However, a batch may can
% contain a single ROI if requested.
load_most_recent_batch    = true ;
load_most_recent_workflow = false ;

% If any of the above switches are set to false then input the previous time stamps below.
previous_batch_time_stamp    = '181203-174759' ;
% 180806-180208
previous_workflow_time_stamp = '190123-180304' ;
% 180827-041056
% 180827-032234
%% ---------------------------------- Downstream of Settings -------------------------------------- 

time_stamp = char( datetime('now', 'TimeZone', 'local', 'Format', 'yyMMdd-HHmmss' )); 

% IF new batch
if productive( 1 ) == 1
    
    number_of_ROIs = length( ROI_suffixes );   

    ROI_index_range = 1 : number_of_ROIs ;
           
    % This is the folder that keeps the vectorization progress of one or more vectorization attempts
    % of a set of batch-processed ROI's
    batch_handle    = [ 'batch_'      , time_stamp        ];
    batch_directory = [ root_directory, batch_handle, '\' ];

    [ ~, ~ ] = mkdir( batch_directory );
    
    % input directories
                    ROI_names = cell( 1, number_of_ROIs );    
              raw_directories = cell( 1, number_of_ROIs );
        processed_directories = cell( 1, number_of_ROIs );
         
    for ROI_index = ROI_index_range
    
        % input directories for the individual ROI's
%                raw_directories{ ROI_index } = [      root_directory, [ ROI_base_name, ROI_suffixes{ ROI_index }]];     
         processed_directories{ ROI_index } =        root_directory ;
             
        ROI_names{ ROI_index } = [ '_', ROI_base_name, ROI_suffixes{ ROI_index }] ;                                             

    end % FOR ROI
    
    % input directories
    input_directories = [       raw_directories; ...
                          processed_directories  ];
        
              ROI_directory = [ batch_directory, ''                ];
             data_directory = [ batch_directory, 'data\'           ];
      visual_data_directory = [ batch_directory, 'visual_data\'    ];
           vector_directory = [ batch_directory, 'vectors\'        ];
    visual_vector_directory = [ batch_directory, 'visual_vectors\' ];
         curation_directory = [ batch_directory, 'curations\'      ];
         
    % output directories
    output_directories = { ROI_directory; ...
                          data_directory; ...
                   visual_data_directory; ...
                        vector_directory; ...
                 visual_vector_directory; ...
                      curation_directory  };         

    % create output directories
    [ ~, ~ ] = mkdir(           ROI_directory );
    [ ~, ~ ] = mkdir(          data_directory );
    [ ~, ~ ] = mkdir(   visual_data_directory );    
    [ ~, ~ ] = mkdir(        vector_directory );
    [ ~, ~ ] = mkdir( visual_vector_directory );
    [ ~, ~ ] = mkdir(      curation_directory );
    
    % This is where all the settings for this batch will be saved.
    settings_directory   = [ batch_directory, 'settings\' ];
    
    [ ~, ~ ] = mkdir( settings_directory );
    
    % This file holds the directory tree that is specific to this batch.
    batch_settings = [ settings_directory, 'batch\' ];        
    
    save( batch_settings        , ...
               'ROI_names'      , ...
            'output_directories', ...
             'input_directories', ...
               'raw_directories', ...
         'processed_directories', ...
               'ROI_directory'  , ...
              'data_directory'  , ...
       'visual_data_directory'  , ...
            'vector_directory'  , ...
     'visual_vector_directory'  , ...
          'curation_directory'  , ...
          'settings_directory'  , ...          
          'ROI_index_range'       );

else % Referencing Old Batch and Workflow

    if load_most_recent_batch

        % find most recent batch time stamp in the root_directory
        batch_listing = dir([ root_directory, 'batch_*-*' ]); 

        there_exists_a_previous_batch = ~ isempty( batch_listing );

        if there_exists_a_previous_batch

            previous_batch_time_stamp = batch_listing( end ).name( 7 : 19 );

        else  

            error('No batches were found for the input root_directory. To create a new batch, set the ordinate_of_start_workflow input to 1.')

        end

    end % load_most_recent_batch

    if load_most_recent_workflow

        % find most recent workflow time stamp in the root_directory
        workflow_listing = dir([ root_directory, 'batch_', previous_batch_time_stamp, '\settings\workflow*.mat' ]); 

        there_exists_a_previous_workflow = ~ isempty( workflow_listing );

        open_most_recent_workflow = true ;

        if there_exists_a_previous_workflow

            previous_workflow_time_stamp = workflow_listing( end ).name( 10 : 22 );

        else  

            error('No workflowes were found for the input root_directory. To create a new workflow, set the ordinate_of_start_workflow input to 1.')

        end

    end % load_most_recent_workflow    
    
    % load an old batch to continue working from
    previous_batch_handle       = [ 'batch_'                   , previous_batch_time_stamp      ];
    previous_batch_directory    = [ root_directory             , previous_batch_handle    , '\' ];
    previous_settings_directory = [ previous_batch_directory   , 'settings\'                    ];
    previous_batch_settings     = [ previous_settings_directory, 'batch\'                       ];   
    
    load( previous_batch_settings )    
    
    % load an old workflow schedule to get the production times to handle the old data
    previous_workflow_handle     = [ 'workflow_'                , previous_workflow_time_stamp ];
    previous_workflow            = [ previous_settings_directory, previous_workflow_handle     ];
    
    load( previous_workflow, 'production_times' )    

end % IF new batch
       
% Overwrite old time stamps that won't be used in this reproduction
production_times( ordinate_of_start_workflow : ordinate_of_final_workflow ) = { time_stamp };

workflow_handle = [ 'workflow_', time_stamp ];

path_to_workflow_settings = [ settings_directory, workflow_handle ];        

save( path_to_workflow_settings, ...
      'workflows'              , ...
      'productive'             , ...
      'production_times'       , ...
      'visual'                 , ...
      'forgetful'                );

%% 1: Process raw data from the microscope and put it into a tif file in the input_directory        
% Note to Programmer: this step should be taken out of this function.  Most users won't want this.
% Instead, start with a folder full of "processed tiffs" and skip to the interpolation step

original_TIF_handle  = 'Fused_T123_T456_Raw.tif' ;
original_data_handle = 'original'                ;

%% ---------------------------------------- Production -------------------------------------------- 

if productive( 1 )
    
    for ROI_index = ROI_index_range
        tic
        
        path_to_original_TIF  = [ processed_directories{ ROI_index }, original_TIF_handle  ]; % TIF file path
        
        path_to_original_data = [ data_directory, original_data_handle, ROI_names{ ROI_index }]; % h5  file path
        
%         process_data( raw_directories{ ROI_index })

        original_data = tif2mat( path_to_original_TIF );
        
        %% !!!!!!!!!!!!!!!!!  IF mirroring the top part of the image
        original_data_temp = cat(3, original_data( :, :, 30 : -1 : 1 ), original_data );
        original_data = original_data_temp;        
        
        size_of_image = size( original_data );
        
        h5create( path_to_original_data, '/d', size_of_image );
        
        mat2h5( path_to_original_data, original_data );
    
        toc
    end % FOR ROI
end % IF productive
%% -------------------------------------- Visualization ------------------------------------------- 

if visual( 1 )
    
    for ROI_index = ROI_index_range
        
        path_to_original_data = [ data_directory, original_data_handle, ROI_names{ ROI_index }]; % h5  file path
        path_to_original_visual = [ visual_data_directory, original_data_handle, ROI_names{ ROI_index }, '.tif' ]; % TIF file path       
        
        mat2tif( h52mat( path_to_original_data ), path_to_original_visual );
        
    end % FOR ROI
end % IF visual

%% 2:  calculate the multi-scale energy field  from the original images                             

energy_handle = [ 'energy_', production_times{ 2 }];

path_to_energy_settings = [ settings_directory, energy_handle ];

%% ---------------------------------------- Production -------------------------------------------- 

if productive( 2 )

    % Energy Settings   
    
%     z_per_xy_length_of_pxl_ratio = 3; %um before 4/22/18 SAM
%     z_per_xy_length_of_pxl_ratio = 2.5; %um  4/22/18 SAM
%     z_per_xy_length_of_pxl_ratio = 2.75; %um  4/24/18 SAM
%     z_per_xy_length_of_pxl_ratio = 4.67; %um  8/6/18 SAM
    z_per_xy_length_of_pxl_ratio = 5.820; %um  11/30/18 SAM
    
    voxel_aspect_ratio = [ 1, 1, z_per_xy_length_of_pxl_ratio ];
        
%     microns_per_pixel_xy                 = 2.75 ; % microns per pxl
%     microns_per_pixel_xy                 = 2.25 ; % microns per pxl  4/22/18
%     microns_per_pixel_xy                 = 2.5  ; % microns per pxl  4/23/18
%     microns_per_pixel_xy                 = 2.25  ; % microns per pxl  4/24/18
%     microns_per_pixel_xy                 = 1.07  ; % microns per pxl  8/6/18
    microns_per_pixel_xy                 = 0.8591  ; % microns per pxl  11/30/18 SAM
    
    microns_per_pixel = voxel_aspect_ratio .* microns_per_pixel_xy ;
    
%     radius_of_smallest_vessel_in_microns = 5    ; % microns per size
%     radius_of_smallest_vessel_in_microns = 4    ; % microns per size 8/7/18
%     radius_of_smallest_vessel_in_microns = 3    ; % microns per size 8/10/18
%     radius_of_smallest_vessel_in_microns = 2    ; % microns per size 8/26/18
    
%     radius_of_smallest_vessel_in_microns = 2      ; % microns per size 10/4/18
    radius_of_smallest_vessel_in_microns = 1      ; % microns per size 10/5/18
%     radius_of_smallest_vessel_in_microns = 30      ; % microns per size 12/11/18

%     radius_of_largest_vessel_in_microns  = 100  ; % microns per size
%     radius_of_largest_vessel_in_microns  = 90   ; % microns per size 4/22/18
%     radius_of_largest_vessel_in_microns  = 40   ; % microns per size 8/7/18

%     radius_of_largest_vessel_in_microns  = 38   ; % microns per size 8/14/18 SAM
%     radius_of_largest_vessel_in_microns  = 42   ; % microns per size 8/26/18 SAM

%     radius_of_largest_vessel_in_microns  = 5   ; % microns per size 12/1/18 SAM
%     radius_of_largest_vessel_in_microns  = 25   ; % microns per size 12/8/18 SAM
%     radius_of_largest_vessel_in_microns  = 75 ; % microns per size 12/11/18 SAM
%     radius_of_largest_vessel_in_microns  = 60 ; % microns per size 12/11/18 SAM
%     radius_of_largest_vessel_in_microns  = 95 ; % microns per size 1/7/19 SAM
    radius_of_largest_vessel_in_microns  = 110 ; % microns per size 1/8/19 SAM

    
%     sigma_per_size                       = 1.7  ; % sigma per size
    sigma_per_size                       = 1 ; % essentially removing this factor from downstream workflows because a new variable, "sigma_er_radius" is taking it's place and it is now behind the scenes in only this workflow.  SAM 8/10/18
    
%     microns_per_sigma_microscope_PSF     = [ 1, 1, 4 ];
%     microns_per_sigma_microscope_PSF     = [ 0, 0, 0 ]; % SAM 4/1/18
%     requested_microns_per_pixel_PSF = 0.5 * [ 1, 1, 4 ];
%     
%     microns_per_sigma_microscope_PSF     = min( requested_microns_per_pixel_PSF, ...
%                                                  starting_microns_per_pixel_xy   ); % SAM 4/2/18 0600
                                             
%     microns_per_sigma_microscope_PSF = 0.5 * [ 1, 1, 2 ]; % SAM 4/2/18 0700
    
    sample_index_of_refraction = 1.33 ;
    
    numerical_aperture = 0.95 ;
    
    excitation_wavelength = 1.3 ; % microns

    microns_per_sigma_microscope_PSF                            ...
        = excitation_wavelength / 2 ^ 0.5                       ...
        * [ .325 / numerical_aperture ^ 0.91 * [ 1, 1 ],        ...
            .532 / ( sample_index_of_refraction                 ...
                     - ( sample_index_of_refraction ^ 2         ...
                         - numerical_aperture ^ 2       ) ^ 0.5 )];
                                         
    % from Figure 4 of Nonlinear magic: multiphoton microscopy in the biosciences by Warren R
    % Zipfel, Rebecca M Williams & Watt W Webb, with sigma of the intensity gaussian = omega for the
    % squared intensity gaussian, valid for NA bigger than 0.7  SAM 4/11/18
    
    % microscope_PSF_fudge_factor is for enlarging the microscope PSF without theoretical
    % justification for the end of smoothing the matching kernel

%     microscope_PSF_fudge_factor = 1 ; % SAM 8/26/18
%     microscope_PSF_fudge_factor = 2 ; % SAM 8/26/18
%     microscope_PSF_fudge_factor = 8 ; % SAM 8/26/18

    microscope_PSF_fudge_factor = 1 ; % SAM 9/4/18 1500

    microns_per_sigma_microscope_PSF = microns_per_sigma_microscope_PSF * microscope_PSF_fudge_factor ;
                              
    pixels_per_sigma_PSF = microns_per_sigma_microscope_PSF ./ microns_per_pixel ;                          
    
%     largest_per_smallest_scale_size_ratio = radius_of_largest_vessel_in_microns  ...
%                                           / radius_of_smallest_vessel_in_microns ;

    % volume interpretation, SAM 8/14/18
    largest_per_smallest_scale_size_ratio = (   radius_of_largest_vessel_in_microns        ...
                                              / radius_of_smallest_vessel_in_microns ) ^ 3 ;

                                      
%     scales_per_octave       = uint8( 8 )    ; % 01/13/18
%     scales_per_octave       = uint8( 4 )    ; % 03/23/18
%     scales_per_octave       = uint8( 6 )    ; % 03/26/18
%     scales_per_octave       = 6 ; % 03/28/18
%     scales_per_octave       = 10 ; % 04/22/18
%     scales_per_octave       = 7 ; % 5/7/18

%     scales_per_octave       = 5 ; % 8/14/18
%     scales_per_octave       = 2 ; % 8/14/18 1500
%     scales_per_octave       = 4 ; % 8/23/18
%     scales_per_octave       = 2 ; % 8/26/18
%     scales_per_octave       = 1 ; % 8/26/18
%     scales_per_octave       = 1/2 ; % 8/26/18

%     scales_per_octave       = 3 ; % 10/?/18
%     scales_per_octave       = 1 ; % 10/8/18
    scales_per_octave       = 2 ; % 1/8/19

%     scales_per_octave       = 1/3 ; % 12/11/18

    
    final_scale = round( log( largest_per_smallest_scale_size_ratio ) / log( 2 ) * scales_per_octave );
    
    % need to pad by 1 scale on the top and bottom because we need at least one above or below to
    % call it an extreme point later when finding vertices (local energy minima in physical and
    % scale spaces (4D locale)).
    scale_ordinates = ( - 1 : final_scale + 1 )' ;
    
%     scale_factor_range = 2 .^ ( scale_ordinates / scales_per_octave );
    
    % volume interpretation, SAM 8/14/18
    scale_factor_range = 2 .^ ( scale_ordinates / scales_per_octave / 3 );
    
    lumen_radius_in_microns_range = radius_of_smallest_vessel_in_microns * scale_factor_range ;
        
    pixels_per_sigma_range = lumen_radius_in_microns_range ./ microns_per_pixel ;
    
%     max_voxels_per_node = 10 ^ 5.5 ; % 8/14/18
    max_voxels_per_node = 10 ^ 6 ; % 1/7/19
        
%     get_energy_version = 'V130' ; % 5/3/18 SAM
%     get_energy_version = 'V131' ; % 5/4/18 SAM
%     get_energy_version = 'V133' ; % 5/4/18 1630 SAM 
%     get_energy_version = 'V131' ; % 5/4/18 2350 SAM 
%     get_energy_version    = 'V140' ; % 5/5/18 0340 SAM 
%     get_energy_version    = 'V151' ; % 5/7/18 0340 SAM 
%     get_energy_version    = 'V160' ; % 5/17/18 SAM 
%     get_energy_version    = 'V161' ; % 5/18/18 SAM 
%     get_energy_version    = 'V162' ; % 5/21/18 SAM 
%     get_energy_version    = 'V163' ; % 5/22/18 SAM 
%     get_energy_version    = 'V170' ; % 5/23/18 SAM 
%     get_energy_version    = 'V171' ; % 5/25/18 0500 SAM 
%     get_energy_version    = 'V172' ; % 5/25/18 1645 SAM 
%     get_energy_version    = 'V171' ; % 5/25/18 1735 SAM 
%     get_energy_version    = 'V162' ; % 5/28/18 SAM 
%     get_energy_version    = 'V193' ; % 8/28/18 SAM 
    get_energy_version    = 'V202' ; % 12/8/18 SAM 
    
%     energy_filter_version = 'V140' ; % 5/5/18 0340 SAM 
%     energy_filter_version = 'V153' ; % 5/6/18 0340 SAM 
%     energy_filter_version = 'V160' ; % 5/17/18 SAM 
%     energy_filter_version = 'V161' ; % 5/18/18 SAM 
%     energy_filter_version = 'V162' ; % 5/21/18 SAM     
%     energy_filter_version = 'V163' ; % 5/22/18 SAM     
%     energy_filter_version = 'V171' ; % 5/25/18 0500 SAM    
%     energy_filter_version = 'V172' ; % 5/25/18 1645 SAM  
%     energy_filter_version = 'V171' ; % 5/25/18 1735 SAM  
%     energy_filter_version = 'V162' ; % 5/28/18 SAM  
    energy_filter_version = 'V191' ; % 8/10/18 SAM 
  
%     symmetry_ratio_factor = 2 ; % SAM 9/4/18 1500
%     symmetry_ratio_factor = 0 ; % 10/4/18
%     symmetry_ratio_factor = 1 ; % SAM 1/8/19
    symmetry_ratio_factor = 1/2 ; % SAM 1/8/19
    
    % this is like sigma per size but it only applies to the gaussian matching kernel
%     sigma_per_radius = 1.7 ; % 8/10/18
    sigma_per_radius = 1 ; % 9/28/18  % this fudge factor has been removed.  In practice it doesn't help at the capillary scale
    
    vessel_wall_thickness_in_microns = 1 ; % SAM 8/13/18 estimate from a cursory literature survey by David Miller
    
%     matching_kernel_string = '3D gaussian' ; % SAM 8/10/18
%     matching_kernel_string = 'spherical pulse' ; % SAM 8/11/18 1300
%     matching_kernel_string = '3D gaussian' ; % SAM 8/13/18 % derivative normalization back to vectorize_V180 style where it is w.r.t. radius
%     matching_kernel_string = 'spherical pulse' ; % SAM 8/14/18 1300

%     matching_kernel_string = '3D gaussian' ; % SAM 8/23/18 1300
%     matching_kernel_string = '3D gaussian conv spherical pulse' ; 

    matching_kernel_string = '3D gaussian' ; % SAM 9/4/18 1500

    save(  path_to_energy_settings          , ...
          'pixels_per_sigma_range'          , ...
          'scales_per_octave'               , ...
          'sigma_per_size'                  , ...  
          'max_voxels_per_node'             , ...
          'pixels_per_sigma_PSF'            , ...
          'voxel_aspect_ratio'              , ...
          'get_energy_version'              , ...
          'energy_filter_version'           , ...
          'symmetry_ratio_factor'           , ...
          'microns_per_pixel_xy'            , ...
          'microns_per_pixel'               , ...
          'matching_kernel_string'          , ...
          'microscope_PSF_fudge_factor'     , ...
          'lumen_radius_in_microns_range'   , ...
          'vessel_wall_thickness_in_microns', ...
          'sigma_per_radius'                  )
      
    for ROI_index = ROI_index_range
        tic
        
        get_energy_V202( matching_kernel_string, lumen_radius_in_microns_range,                  ...
                         vessel_wall_thickness_in_microns, microns_per_pixel,  ...
                         pixels_per_sigma_PSF, max_voxels_per_node, output_directories{ 2 }, ...
                         [ original_data_handle, ROI_names{ ROI_index }], ...
                         [        energy_handle, ROI_names{ ROI_index }], symmetry_ratio_factor              )
        
        toc
    end % FOR ROI    
end % IF productive
%% -------------------------------------- Visualization ------------------------------------------- 

if visual( 2 )
    
    for ROI_index = ROI_index_range
        
        path_to_energy_data = [        data_directory, energy_handle, ROI_names{ ROI_index }];               
        energy_visual_file  = [ visual_data_directory, energy_handle, ROI_names{ ROI_index }];
                
        energy_file_info = h5info( path_to_energy_data );

        size_of_image = energy_file_info.Datasets.Dataspace.Size ;        
            
        mat2tif( h52mat( path_to_energy_data,           ...
                        [ 1, 1, 1, 2 ],                 ...
                        [ size_of_image( 1 : 3 ), 1 ]), ...
                 [ energy_visual_file, '.tif' ]         )
             
        mat2tif( h52mat( path_to_energy_data,           ...
                        [ 1, 1, 1, 1 ],                 ...
                        [ size_of_image( 1 : 3 ), 1 ]), ...
                 [ energy_visual_file, '_indices.tif' ] )
             
    end % FOR ROI
end % IF visual
%% ----------------------------------------- Deletion --------------------------------------------- 

if forgetful( 1 )
    for ROI_index = ROI_index_range
        
        delete([ data_directory, original_data_handle, ROI_names{ ROI_index }])
        
    end % FOR ROI
end % IF forgetul

%% 3: get the vertices from the local minima of the 4D energy field                                 

vertices_handle = [ 'vertices_', production_times{ 3 }];

path_to_vertices_settings = [ settings_directory, vertices_handle ];

exporting_vertex_curation = false ;

%% ----------------------------------------- Production ------------------------------------------- 
if productive( 3 )
    
    load( path_to_energy_settings )
    
    % vertex settings
    space_strel_apothem     = 1 ; % 4/22/18

%     energy_upper_bound      = -3 ; % 11/2/18
    energy_upper_bound      = -4 ; % 12/13/18

    max_voxels_per_node     = 6000 ;
    
    get_vertices_version  = 'V200' ; % 12/4/18

    crop_vertices_version = 'V200' ; % 12/4/18
    
%     sigma_per_influence = sigma_per_size ; % SAM 6/5/18
    sigma_per_influence = 2 ; % 1/22/19
%     sigma_per_influence = 1 ; % 1/22/19

    save( path_to_vertices_settings, ...
          'space_strel_apothem'    , ...
          'energy_upper_bound'  ,    ...
          'max_voxels_per_node' ,    ...
          'get_vertices_version',    ...
          'crop_vertices_version',   ...
          'sigma_per_influence'      );

    for ROI_index = ROI_index_range
        tic
        
        path_to_energy_data        = [     data_directory,        energy_handle, ROI_names{ ROI_index }]; % mat file path        
        path_to_original_data      = [     data_directory, original_data_handle, ROI_names{ ROI_index }]; % h5  file path
        path_to_saved_curation     = [ curation_directory,      vertices_handle, ROI_names{ ROI_index }]; % logicals path
        
        energy_file_info = h5info( path_to_energy_data );

        size_of_image = energy_file_info.Datasets.Dataspace.Size ;  
        
        size_of_image = size_of_image( 1 : 3 );

        [ vertex_space_subscripts, vertex_scale_subscripts, vertex_energies ]                       ...
            = get_vertices_V200(                  lumen_radius_in_microns_range, microns_per_pixel, ...
                                      space_strel_apothem, max_voxels_per_node, energy_upper_bound, ...
                                    [ output_directories{ 2 }, energy_handle, ROI_names{ ROI_index  }]);
                            
        path_to_vertices_data = [ vector_directory, vertices_handle, ROI_names{ ROI_index }];
        
        [ vertex_space_subscripts, vertex_scale_subscripts, vertex_energies ]                       ...
                         = crop_vertices_V200( vertex_space_subscripts, vertex_scale_subscripts,    ...
                                               vertex_energies, sigma_per_influence * lumen_radius_in_microns_range,      ...
                                               microns_per_pixel, size_of_image                     );
                                           
        [ vertex_energies, sorted_indices ] = sort( vertex_energies );
                
        vertex_space_subscripts = vertex_space_subscripts( sorted_indices, : );
        vertex_scale_subscripts = vertex_scale_subscripts( sorted_indices    );
        
        save(   path_to_vertices_data,     ...
                'vertex_space_subscripts', ...
                'vertex_scale_subscripts', ...
                'vertex_energies'          );
            
        if exporting_vertex_curation
            % run this line to export curation inputs to a third party
            save([ path_to_saved_curation, '_curator_inputs' ], 'vertex_energies', ...
                             'vertex_space_subscripts', 'vertex_scale_subscripts', ...
                             'lumen_radius_in_microns_range', 'microns_per_pixel', ...
                                'path_to_original_data', 'path_to_saved_curation', ...
                                                             'path_to_energy_data' )            

        end
        toc
    end % FOR ROI       
end % IF productive
%% ------------------------------------------ Curation -------------------------------------------- 

for ROI_index = ROI_index_range
    
    path_to_energy_data             = [     data_directory,               energy_handle, ROI_names{ ROI_index }]; % mat file path        
    path_to_original_data           = [     data_directory,        original_data_handle, ROI_names{ ROI_index }]; % h5  file path
    path_to_vertices_data           = [   vector_directory,             vertices_handle, ROI_names{ ROI_index }]; %  vectors path
    path_to_curated_vertices_data   = [   vector_directory, 'curated_', vertices_handle, ROI_names{ ROI_index }]; %  vectors path    
    path_to_saved_curation          = [ curation_directory,             vertices_handle, ROI_names{ ROI_index }]; % logicals path

    original_file_info = h5info( path_to_original_data );

    size_of_image = original_file_info.Datasets.Dataspace.Size ;  

    load( path_to_vertices_data   )
    load( path_to_energy_settings )

    tic    
    switch curation{ 3 }
        
        case 'manual'
                        
            [ vertex_energies, vertex_space_subscripts, vertex_scale_subscripts ]                   ...
                      = vertex_curator_V5(                vertex_energies, vertex_space_subscripts, ...
                                            vertex_scale_subscripts, sigma_per_influence * lumen_radius_in_microns_range, ...
                                                          microns_per_pixel, path_to_original_data, ...
                                                       path_to_saved_curation, path_to_energy_data, ...
                                                                            [ 0, 16 ], [ -16, 0 ]);

        case 'auto'

            % instead do a volume conflict test and select the best contrast object at conflicts
            [ chosen_vertex_indices ]                                                    ...
                           = choose_vertices_V191( vertex_space_subscripts, vertex_scale_subscripts,    ...
                                                   vertex_energies, sigma_per_influence * lumen_radius_in_microns_range,      ...
                                                   microns_per_pixel, size_of_image                     );

            % performing the selections proposed by either choose_vertices or vertex_curator
            vertex_space_subscripts = vertex_space_subscripts( chosen_vertex_indices, : );
            vertex_scale_subscripts = vertex_scale_subscripts( chosen_vertex_indices    );
            vertex_energies         =         vertex_energies( chosen_vertex_indices    );    
            
        otherwise % do nothing

    end % IF curation
    toc                                                       
      
    switch curation{ 3 }
        case { 'manual', 'auto' }
    
            save(   path_to_curated_vertices_data, ...
                    'vertex_space_subscripts', ...
                    'vertex_scale_subscripts', ...
                    'vertex_energies'          );   
                
        otherwise % do nothing

    end
end % ROI FOR    
%% --------------------------------------- Visualization ------------------------------------------ 

if visual( 3 )
    
    load( path_to_energy_settings )        
        
    for ROI_index = ROI_index_range
        
        path_to_vertices_data              = [        vector_directory,             vertices_handle, ROI_names{ ROI_index }         ];
        path_to_curated_vertices_data      = [        vector_directory, 'curated_', vertices_handle, ROI_names{ ROI_index }         ]; %  vectors path            
        path_to_vertices_visuals           = [ visual_vector_directory,             vertices_handle, ROI_names{ ROI_index }, '.tif' ];    
        path_to_curated_vertices_visuals   = [ visual_vector_directory, 'curated_', vertices_handle, ROI_names{ ROI_index }, '.tif' ]; %  vectors path            
        path_to_energy_data                = [          data_directory,             energy_handle,   ROI_names{ ROI_index }         ];       
                
        energy_file_info = h5info( path_to_energy_data );

        size_of_image = energy_file_info.Datasets.Dataspace.Size ;  
        
        size_of_image = size_of_image( 1 : 3 );
        
        load( path_to_vertices_data )
                                                     
        visualize_vertices_V200(  vertex_space_subscripts, ...
                                  vertex_scale_subscripts, ...
                                          vertex_energies, ...
                    pixels_per_sigma_range, size_of_image, ...
                                 path_to_vertices_visuals  )
                       
        load( path_to_curated_vertices_data   )
                                                     
        visualize_vertices_V200(  vertex_space_subscripts, ...
                                  vertex_scale_subscripts, ...
                                          vertex_energies, ...
                    pixels_per_sigma_range, size_of_image, ...
                         path_to_curated_vertices_visuals  )
                                                          
    end % FOR ROI
end % IF visual

%% 4: get the edges by following vessel leads probabilistically to other vertices                   

edges_handle = [ 'edges_', production_times{ 4 }];

path_to_edges_settings = [ settings_directory, edges_handle ]; 

%% ----------------------------------------- Production ------------------------------------------- 
if productive( 4 )

    load( path_to_energy_settings )
    load( path_to_vertices_settings )    
        
    pixels_per_sigma_aspect_ratio = 1 ./ voxel_aspect_ratio;
    
    %% ----------------------------------------- settings ------------------------------------------
    
%     maximum_edge_degree = 128; % 4/3/18 SAM
%     maximum_edge_degree = 256; % 4/22/18 SAM   
%     maximum_edge_degree = 128; % 4/22/18 SAM    
    % note, repalce the maximum edge degree to a minimum_current_probability function.  Don't forget
    % to make the edge degree a uint16 instead of uint8 to allow for more than 256 degrees
    
%     probability_minimum = 10 ^ -5 ; % 5/4/18 SAM
%     probability_minimum = 10 ^ -5 ; % 5/4/18 1700 SAM
%     probability_minimum = 10 ^ -10 ; % 5/4/18 1830 SAM
%     probability_minimum = 10 ^ -5 ; % 5/22/18 1830 SAM
%     probability_minimum = 10 ^ -50 ; % 5/30/18 0140 SAM
%     probability_minimum = 10 ^ -25 ; % 5/30/18 2300
%     probability_minimum = 10 ^ -100 ; % 5/31/18 0240
%     probability_minimum = 10 ^ -70 ; % 8/14/18 0240

    probability_minimum = 20 ; % 8/23/18 %this is an edge length (per radius length) not a probability

    space_strel_apothem = 1 ;  % 4/7/18 SAM
    
%     scale_strel_apothem = Inf ; % 4/10/18 SAM
%     scale_strel_apothem = 2  ; % 4/13/18 SAM    
%     scale_strel_apothem = 1 ; % 4/10/18 SAM
%     scale_strel_apothem = 2  ; % 4/22/18 SAM
%     scale_strel_apothem = 10 ; % 5/4/18 1745 SAM
%     scale_strel_apothem = 5 ; % 5/4/18 1830 SAM
    
    % note: put scale_strel_apothem in terms of percent size change per pixel move, then convert
    % that to a scale_strel_apothem using the scales per octave
    
%     number_of_edges_per_vertex = 100 ;
%     number_of_edges_per_vertex = 20 ; % 5/4/18
%     number_of_edges_per_vertex = 100 ; % 5/4/18 1830
%     number_of_edges_per_vertex = 10 ; % 5/22/18 1830
%     number_of_edges_per_vertex = 20 ; % 5/29/18 1830
%     number_of_edges_per_vertex = 3 ; % 5/30/18 2300
%     number_of_edges_per_vertex = 10 ; % 5/31/18 0324
%     number_of_edges_per_vertex = 20 ; % 5/31/18 0400
%     number_of_edges_per_vertex = 100 ; % 5/31/18 0450
%     number_of_edges_per_vertex = 3 ; % 5/31/18 1013
%     number_of_edges_per_vertex = 10 ; % 5/31/18 1040
%     number_of_edges_per_vertex = 20 ; % 5/31/18 1200
%     number_of_edges_per_vertex = 50 ; % 5/31/18 1600
%     number_of_edges_per_vertex = 10 ; % 8/14/18

    number_of_edges_per_vertex = 50 ; % 8/23/18


%     characteristic_energy_fraction = 0.5 ; % 4/22/18 1315 SAM
%     characteristic_energy_fraction = 5 ; % 4/22/18 1648 SAM   
%     characteristic_energy_fraction = 1 ; % 5/3/18 SAM   
%     characteristic_energy_fraction = 1/10 ; % 5/4/18 SAM   
%     characteristic_energy_fraction = 1/4 ; % 5/4/18 1700 SAM   
%     characteristic_energy_fraction = 1/10 ; % 5/4/18 2000 SAM   
%     characteristic_energy_fraction = 1/100 ; % 5/30/18 0000 SAM   
%     characteristic_energy_fraction = 1 ; % 5/30/18 0140 SAM   
%     characteristic_energy_fraction = 1/10 ; % 5/30/18 0230 SAM   
%     characteristic_energy_fraction = 1/2 ; % 5/31/18 0300 SAM   
%     characteristic_energy_fraction = 1/10 ; % 5/31/18 0300 SAM   

    characteristic_energy_fraction = 100 ; % 8/23/18 SAM new definition of energy fraction with get_edges V190

    get_edges_version   = 'V190' ; % 11/2/18 SAM
    
%     choose_edges_version = 'V161' ; % 5/27/18 SAM
%     choose_edges_version = 'V162' ; % 5/28/18 SAM
%     choose_edges_version = 'V160' ; % 5/28/18 SAM
%     choose_edges_version = 'V181' ; % 5/30/18 SAM
%     choose_edges_version = 'V182' ; % 5/31/18 SAM
    choose_edges_version = 'V190' ; % 8/2/18 SAM
    
%     sigma_per_influence_vertices = sigma_per_influence ;

    sigma_per_influence_vertices = 1 ;
    
%     sigma_per_influence_edges    = 1 ;
%     sigma_per_influence_edges    = sigma_per_influence_vertices / 2 ; % SAM 5/31/18 0450
%     sigma_per_influence_edges    = sigma_per_influence_vertices ; % SAM 5/31/18 1025
    sigma_per_influence_edges    = sigma_per_influence_vertices / 2 ; % SAM 5/31/18 1110
    
    save(  path_to_edges_settings,          ...
          'probability_minimum',            ...
          'space_strel_apothem',            ...
          'number_of_edges_per_vertex',     ...
          'characteristic_energy_fraction', ...
          'get_edges_version',              ...
          'choose_edges_version',           ...
          'sigma_per_influence_edges',      ...
          'sigma_per_influence_vertices'    );
          
    for ROI_index = ROI_index_range
        tic
        
        path_to_curated_vertices_data = [ vector_directory, 'curated_', vertices_handle, ROI_names{ ROI_index }]; %  vectors path            
        path_to_edges_data            = [ vector_directory,                edges_handle, ROI_names{ ROI_index }];  
        path_to_energy_data           = [   data_directory,               energy_handle, ROI_names{ ROI_index }];       
        
        energy_file_info = h5info( path_to_energy_data );

        size_of_image = energy_file_info.Datasets.Dataspace.Size ;  
        
        size_of_image = size_of_image( 1 : 3 );
        
        load( path_to_curated_vertices_data )
                
        [ edges2vertices, edge_lengths, edge_space_subscripts, edge_scale_subscripts, edge_energies ]                        ...
            = get_edges_V200( pixels_per_sigma_range, vertex_scale_subscripts,                      ...
                              vertex_space_subscripts, vertex_energies, space_strel_apothem,        ...
                              probability_minimum, number_of_edges_per_vertex,                      ...
                              characteristic_energy_fraction,                                       ...
                              [ input_directories( :, ROI_index ); output_directories ],            ...
                              [ energy_handle, ROI_names{ ROI_index }]                              );
                          
                          
        
                          
%         % clean up the output to just keep the trajectories that found a neighbor and to only keep
%         % the "best" trajectory from each pair of vertices. Choosing the better trajectory from A
%         to B and B to A.  !!! Any direct visualization of the trajectory variability % should be
%         done before this step !!!
        [ edge_space_subscripts, edge_scale_subscripts, edge_energies, edge_lengths, edges2vertices ]                   ...
                          = crop_edges_V200( edge_space_subscripts, edge_scale_subscripts, edge_energies, edge_lengths, ...
                                             edges2vertices, lumen_radius_in_microns_range, microns_per_pixel, ...
                                             size_of_image                                     );
                
        % SAM and DRM 1/10/19
        % manual addition of two large surface vessels
        %
        % If re-running this data set, place the 6 vertices corresponding to the 1 large vessel in a
        % row at the start of the curator.  And then add 2 vertices for the otehr large vessel that
        % is kind of out of view.  #kludgeLife
        
        % vessel A
        for vertex_index = 1 : 5
        
            edges2vertices( end + 1, 1 : 2 ) = uint32([ vertex_index, vertex_index + 1 ]);
            edge_space_subscripts{ end + 1 } = [ vertex_space_subscripts( vertex_index, : ); vertex_space_subscripts( vertex_index + 1, : )];
            edge_scale_subscripts{ end + 1 } = [ vertex_scale_subscripts( vertex_index, : ); vertex_scale_subscripts( vertex_index + 1, : )];
            edge_energies{ end + 1 } = [ -1000; -1000];      
            edge_lengths( end + 1 ) = 2 ;
        
        end
        
        % vessel B
        edges2vertices( end + 1, 1 : 2 ) = uint32([ 7, 8 ]);
        edge_space_subscripts{ end + 1 } = [ vertex_space_subscripts( 7, : ); vertex_space_subscripts( 8, : )];
        edge_scale_subscripts{ end + 1 } = [ vertex_scale_subscripts( 7, : ); vertex_scale_subscripts( 8, : )];
        edge_energies{ end + 1 } = [ -1000; -1000];      
        edge_lengths( end + 1 ) = 2 ;
                                         
        % saving edge outputs
        save(  path_to_edges_data, ...
              'edge_energies'    , ...
          'edge_space_subscripts', ... 
          'edge_scale_subscripts', ...
                   'edge_lengths', ...
                 'edges2vertices'  );
          
        toc
    end % FOR ROI
end % IF productive
              
path_to_visual_edges_settings = [ settings_directory, 'visual_', edges_handle ];
%% ------------------------------------------ Curation -------------------------------------------- 

for ROI_index = ROI_index_range
        
    path_to_energy_data             = [     data_directory,               energy_handle, ROI_names{ ROI_index }]; % mat file path        
    path_to_original_data           = [     data_directory,        original_data_handle, ROI_names{ ROI_index }]; % h5  file path
    path_to_edges_data              = [   vector_directory,                edges_handle, ROI_names{ ROI_index }];  
    path_to_curated_edges_data      = [   vector_directory, 'curated_',    edges_handle, ROI_names{ ROI_index }]; %  vectors path    
    path_to_saved_curation          = [ curation_directory,                edges_handle, ROI_names{ ROI_index }]; % logicals path

    original_file_info = h5info( path_to_original_data );

    size_of_image = original_file_info.Datasets.Dataspace.Size ;  

    load( path_to_edges_data      )
    load( path_to_energy_settings )

    tic
    switch curation{ 4 }
        
        case 'manual'
            
            [ mean_edge_energies, edges2vertices, edge_energies, edge_space_subscripts, edge_scale_subscripts ] ...
                       = edge_curator_V2( edge_energies, edge_space_subscripts, edge_scale_subscripts,             ...
                                          edges2vertices, edge_lengths, lumen_radius_in_microns_range,         ...
                                          microns_per_pixel, path_to_original_data, path_to_saved_curation,        ...
                                          path_to_energy_data, [ 0, 16 ], [ -16, 0 ]);   
        
        edge_lengths = []; % put this in the curator function
                                              
        case 'auto'

            % instead do a volume conflict test and select the best contrast object at conflicts                                           
            [ mean_edge_energies, chosen_edge_indices, edges2vertices ]                                 ...
                              = choose_edges_V200( edges2vertices, edge_energies, edge_space_subscripts, edge_scale_subscripts,      ...
                                                   edge_lengths, pixels_per_sigma_range,            ...
                                                   sigma_per_influence_vertices,                        ...
                                                   sigma_per_influence_edges, size_of_image             );

            % performing the selections proposed by either choose_edges or edge_curator
            edge_lengths          = edge_lengths(          chosen_edge_indices );              
            edge_energies         = edge_energies(         chosen_edge_indices );
            edge_space_subscripts = edge_space_subscripts( chosen_edge_indices );
            edge_scale_subscripts = edge_scale_subscripts( chosen_edge_indices );
    
        otherwise % do nothing

    end % IF curation
    toc                                                        
      
    switch curation{ 4 }
        case { 'manual', 'auto' }
            
            save( path_to_curated_edges_data, ...
                             'edge_energies', ...
                     'edge_space_subscripts', ...
                     'edge_scale_subscripts', ...
                            'edges2vertices', ...
                              'edge_lengths', ...
                        'mean_edge_energies'  );      
           
        otherwise % do nothing           
           
    end
end % ROI FOR  
%% --------------------------------------- Visualization ------------------------------------------ 

if visual( 4 )
    
    load( path_to_energy_settings ) 
          
    for ROI_index = ROI_index_range
        
        path_to_energy_data                    = [          data_directory,            energy_handle,                 ROI_names{ ROI_index }         ];
        
        path_to_edges_data                     = [        vector_directory,             edges_handle,                 ROI_names{ ROI_index }         ]; %  vectors path 
        path_to_curated_edges_data             = [        vector_directory, 'curated_', edges_handle,                 ROI_names{ ROI_index }         ]; %  vectors path            
        
        path_to_spheres_visual_file            = [ visual_vector_directory,             edges_handle, '_spheres',     ROI_names{ ROI_index }, '.tif' ];
        path_to_curated_spheres_visual_file    = [ visual_vector_directory, 'curated_', edges_handle, '_spheres',     ROI_names{ ROI_index }, '.tif' ];

        path_to_centerline_visual_file         = [ visual_vector_directory,             edges_handle, '_centerlines', ROI_names{ ROI_index }, '.tif' ];
        path_to_curated_centerline_visual_file = [ visual_vector_directory, 'curated_', edges_handle, '_centerlines', ROI_names{ ROI_index }, '.tif' ];       

        energy_file_info = h5info( path_to_energy_data );

        size_of_image = energy_file_info.Datasets.Dataspace.Size ;  
        
        size_of_image = size_of_image( 1 : 3 );
        
%         load( path_to_edges_data )        
%         
%         edge_subscripts = cellfun(      @(       space_subscripts,      scale_subscripts  )  ...
%                                    double([      space_subscripts,      scale_subscripts ]), ...
%                                             edge_space_subscripts, edge_scale_subscripts,    ...
%                                                   'UniformOutput', false                     );
% 
%         visualize_edges_V180( edge_subscripts, mean_edge_energies, pixels_per_sigma_range,  ...
%                               sigma_per_size, size_of_image, path_to_spheres_visual_file,   ...
%                               path_to_centerline_visual_file                                )
%                           
%         visualize_depth_via_color_V6( edge_subscripts, mean_edge_energies, pixels_per_sigma_range, ...
%                                                                sigma_per_size, output_directories, ...
%                                                   [ original_data_handle, ROI_names{ ROI_index }], ...
%                                       [ 1, 1505 ], [ 1, 980 ], [ 40, 80 ], 0, [ - 12, 160 ])                                  
                    
        load( path_to_curated_edges_data )        
        
        edge_subscripts = cellfun(      @(       space_subscripts,      scale_subscripts  )  ...
                                   double([      space_subscripts,      scale_subscripts ]), ...
                                            edge_space_subscripts, edge_scale_subscripts,    ...
                                                  'UniformOutput', false                     );

        visualize_edges_V180( edge_subscripts, mean_edge_energies, pixels_per_sigma_range,        ...
                              sigma_per_size, size_of_image, path_to_curated_spheres_visual_file, ...
                              path_to_curated_centerline_visual_file                              )
                          
        visualize_depth_via_color_V6( edge_subscripts, mean_edge_energies, pixels_per_sigma_range, ...
                                                               sigma_per_size, output_directories, ...
                                                  [ original_data_handle, ROI_names{ ROI_index }], ...
                                                [ 1, 1505 ], [ 1, 980 ], [ 40, 80 ], 0, [ - 12, 160 ])                                  
                                  
    end % FOR ROI
end % IF visual
%% ------------------------------------------ Deletion -------------------------------------------- 

% delete energy data
if forgetful( 3 )
    
    for ROI_index = ROI_index_range
        
        load( path_to_downsample_settings )
        
        path_to_energy_data = [ data_directory, energy_handle, ROI_names{ ROI_index }];  

        delete( path_to_energy_data )
    
    end % FOR ROI
end % IF forgetful

%% 5: assign strand edges and bifurcation vertices to the network                                   

network_handle = [ 'network_', production_times{ 5 }];

path_to_network_settings = [ settings_directory, network_handle ]; 

%% ----------------------------------------- Production ------------------------------------------- 

if productive( 5 )

%     load( path_to_energy_settings   )
%     load( path_to_vertices_settings )  
    
    % get the characteristic energy fraction from the edges settings
    
    load( path_to_edges_settings )
    % loading the chosen edges, their positions and max_energies, and the
    % characteristic_energy_fraction used to find them
    
    %% ----------------------------------------- settings ------------------------------------------
    
%     get_network_version = 'V180' ; % 6/1/18
%     get_network_version = 'V181' ; % 6/12/18
%      get_network_version = 'V182' ; % 7/11/18
%      get_network_version = 'V183' ; % 8/14/18
     get_network_version = 'V190' ; % 9/1/18
    
    sort_network_version = 'V180' ; % 7/11/18
    
%     sigma_strand_smoothing = 1 ; % 7/16/18
%     sigma_strand_smoothing = 3 ; % 7/16/18 1800
%     sigma_strand_smoothing = 6 ; % 7/18/18 1800

    sigma_strand_smoothing = 2 ; % units of vertex radius 7/18/18 1800
    
    save(  path_to_network_settings, ...
           'get_network_version'   , ...
          'sort_network_version'   , ...
          'sigma_strand_smoothing'   );
          
    for ROI_index = ROI_index_range
        tic
        
        path_to_curated_edges_data        = [ vector_directory, 'curated_',    edges_handle, ROI_names{ ROI_index }]; %  vectors path    
        path_to_network_data              = [ vector_directory,              network_handle, ROI_names{ ROI_index }];          

        load( path_to_curated_edges_data )        
                
        [ bifurcation_vertices, vertex_indices_in_strands,                                          ...
           edge_indices_in_strands, end_vertices_of_strands ]      = get_network_V190( edges2vertices );
                                            
        % sort the strand output
        [ vertex_indices_in_strands, edge_indices_in_strands, edge_backwards_in_strands ]           ...
                                      = sort_network_V180( edges2vertices, end_vertices_of_strands,   ...
                                                           edge_indices_in_strands                  );
                                            
        edge_subscripts = cellfun(      @(       space_subscripts,      scale_subscripts  )  ...
                                   double([      space_subscripts,      scale_subscripts ]), ...
                                            edge_space_subscripts, edge_scale_subscripts,    ...
                                                  'UniformOutput', false                     );
                          
        % save function 7/13/18
        vertex_subscripts = [ double( vertex_space_subscripts ), double( vertex_scale_subscripts )];
        
        load( path_to_energy_settings )    
        
        z_per_xy_length_of_pxl_ratio = 5.820 ; % kludge SAM 1/22/19
        
        save('chakameh_V403','vertex_indices_in_strands',    ...
                               'edge_indices_in_strands',    ...
                             'edge_backwards_in_strands',    ...
                             'bifurcation_vertices',         ...
                             'edge_subscripts',              ...
                             'mean_edge_energies',            ...
                             'pixels_per_sigma_range',       ...
                             'sigma_per_size',               ...
                             'size_of_image',                ...
                             'microns_per_pixel_xy',         ...
                             'z_per_xy_length_of_pxl_ratio', ...
                             'vertex_subscripts'             )                                              
                                              
        % smooth the edge subscripts along the strands
        [ vessel_directions_from_strands, edge_subscripts_from_strands ]                   ...
             = get_vessel_directions_V5( edge_subscripts, edge_indices_in_strands,         ...
                                         edge_backwards_in_strands, sigma_strand_smoothing );
                                                                          
        % overwrite the old edge_subscripts with the smoothed edge_subscripts
        number_of_edges = length( edge_subscripts );
        
        % manual addition of edges: first six objects should not be smoothed. SAM 1/23/19
        edge_subscripts_from_strands( end - 5 : end ) = edge_subscripts( end - 5 : end );
        
        edge_subscripts   = cell( number_of_edges, 1 );
        vessel_directions = cell( number_of_edges, 1 );
        
        edge_index_range = 1 : number_of_edges ;
        
        for edge_index = edge_index_range 
            
            % IF edge index re-defined by the strands smoothing
            if ~ isempty( edge_subscripts_from_strands{ edge_index })
            
                  edge_subscripts{ edge_index } =   edge_subscripts_from_strands{ edge_index };
                vessel_directions{ edge_index } = vessel_directions_from_strands{ edge_index };

            end % IF edge index re-defined by strand smoothing
            
        end % edges FOR
                                                                                 
        % saving network outputs
        save(  path_to_network_data        , ...
              'vertex_indices_in_strands'  , ...
                'edge_indices_in_strands'  , ...
                'end_vertices_of_strands'  , ...
              'edge_backwards_in_strands'  , ...
              'bifurcation_vertices'       , ...
              'edge_subscripts'            , ...
              'vessel_directions'            );
          
        toc
    end % FOR ROI
end % IF productive
              
path_to_visual_strands_settings = [ settings_directory, 'visual_', network_handle ];
%% --------------------------------------- Visualization ------------------------------------------ 

if visual( 5 )
    
    load( path_to_energy_settings )    

    for ROI_index = ROI_index_range
                        
        path_to_energy_data                = [   data_directory,               energy_handle, ROI_names{ ROI_index }];           
        path_to_curated_vertices_data      = [ vector_directory, 'curated_', vertices_handle, ROI_names{ ROI_index }];
        path_to_curated_edges_data         = [ vector_directory, 'curated_',    edges_handle, ROI_names{ ROI_index }]; %  vectors path    
        path_to_network_data               = [ vector_directory,              network_handle, ROI_names{ ROI_index }];  
                
        strands_visual_start_vertices_file = [ visual_vector_directory,  network_handle, '_strands', '_start_vertices', ROI_names{ ROI_index }, '.tif' ];
        strands_visual_end_vertices_file   = [ visual_vector_directory,  network_handle, '_strands',   '_end_vertices', ROI_names{ ROI_index }, '.tif' ]; 
        strands_visual_bifurcations_file   = [ visual_vector_directory,  network_handle, '_strands',   '_bifurcations', ROI_names{ ROI_index }, '.tif' ];             
        strands_visual_spheres_file        = [ visual_vector_directory,  network_handle, '_strands',        '_spheres', ROI_names{ ROI_index }, '.tif' ];
        strands_visual_centerline_file     = [ visual_vector_directory,  network_handle, '_strands',    '_centerlines', ROI_names{ ROI_index }, '.tif' ];      
        
        load( path_to_curated_edges_data    )
        load( path_to_curated_vertices_data )
        load( path_to_network_data          )
           
        energy_file_info = h5info( path_to_energy_data );

        size_of_image = energy_file_info.Datasets.Dataspace.Size ;  
        
        size_of_image = size_of_image( 1 : 3 );
        
%         logical_nonring_strands = ~ cellfun( 'isempty', end_vertices_of_strands );
% 
%         start_vertex_indices = cellfun( @( v ) v( 1 ), end_vertices_of_strands( logical_nonring_strands ));
%         end_vertex_indices   = cellfun( @( v ) v( 2 ), end_vertices_of_strands( logical_nonring_strands ));
% 
%         %% network start, end, and bifurcation vertices visualization                               
%         visualize_vertices_V190( vertex_space_subscripts( start_vertex_indices, : ), ...
%                                  vertex_scale_subscripts( start_vertex_indices    ), ...
%                                          vertex_energies( start_vertex_indices    ), ...
%                                  pixels_per_sigma_range, sigma_per_size,             ...
%                                  size_of_image, strands_visual_start_vertices_file   )
%                              
%         visualize_vertices_V190( vertex_space_subscripts(   end_vertex_indices, : ), ...
%                                  vertex_scale_subscripts(   end_vertex_indices    ), ...
%                                          vertex_energies(   end_vertex_indices    ), ...
%                                  pixels_per_sigma_range, sigma_per_size,             ...
%                                  size_of_image, strands_visual_end_vertices_file     )

        %% network bifurction vertices visualization                                          
%         visualize_vertices_V190( vertex_space_subscripts( bifurcation_vertices, : ), ...
%                                  vertex_scale_subscripts( bifurcation_vertices    ), ...
%                                          vertex_energies( bifurcation_vertices    ), ...
%                                  pixels_per_sigma_range, sigma_per_size,             ...
%                                  size_of_image, strands_visual_bifurcations_file    )                             
        %% network depth visualization                                          
        % visualize strands
        edge_strand_indices   = cell2mat( edge_indices_in_strands );
                        
%         visualize_edges_V180( edge_subscripts(   edge_strand_indices ),   ...
%                               mean_edge_energies( edge_strand_indices ),   ...
%                               pixels_per_sigma_range, sigma_per_size,     ...
%                               size_of_image, strands_visual_spheres_file, ...
%                               strands_visual_centerline_file              )
%                                                     
%         visualize_depth_via_color_V5( edge_subscripts(   edge_strand_indices ),                 ...
%                                       mean_edge_energies( edge_strand_indices ),                 ...
%                                       pixels_per_sigma_range, sigma_per_size,                   ...
%                                       [ input_directories( :, ROI_index ); output_directories ], [ original_data_handle, ROI_names{ ROI_index }],        ...
%                                       [ 1, 1505 ], [ 1, 980 ], [ 40, 80 ], -5, [ - 12, 160 ])
        %% 2D strand visualization (random color assignment)                                        
                                    
        number_of_strands = length( edge_indices_in_strands );
        
        strand_index_range = 1 : number_of_strands ;
        
        edge_strand_assignments_in_strands = edge_indices_in_strands ;        
        
        % label the edges of each strand by the strand index of that strand for later colormapping
        for strand_index = strand_index_range
            
            edge_strand_assignments_in_strands{ strand_index }( : ) = strand_index ;
            
        end % strand FOR
                                    
        edge_strand_assignments = cell2mat( edge_strand_assignments_in_strands );
                                                
        visualize_strands_via_color_V2(     edge_subscripts(   edge_strand_indices ),              ...
                                            mean_edge_energies( edge_strand_indices ),              ...
                                            edge_strand_assignments,                           ...
                                            pixels_per_sigma_range, sigma_per_size,            ...
                                            [ input_directories( :, ROI_index ); output_directories ], [ original_data_handle, ROI_names{ ROI_index }], ...
                                            [ 1, 1505 ], [ 1, 980 ], [ 40, 80 ], -4,         ...
                                            [ -12, 160 ]                                   )
        %% 3D strand visualization (random color assignment)                                        

%         max_edge_energies( edge_junction_indices ) = 0 ; % to make image without junctions 

         visualize_strands_via_color_3D_V2( edge_subscripts(   edge_strand_indices ),          ...
                                            mean_edge_energies( edge_strand_indices ),          ...
                                            edge_strand_assignments, microns_per_pixel_xy, ...
                                            pixels_per_sigma_range, sigma_per_size, 0.5,     ...
                                   ...         [ 100, 1475 ], [ 301, 500 ], [ 0, 180 ], -2    )
                                            [ 1, 1505 ], [ 1, 980 ], [ 0, 75 ], -2 )
                                    
    end % FOR ROI
end % IF visual

%% Note:  make the following changes when you get the chance:                                       
%
% SAM 2/9/18
%
% 3 ) we should have separate feature thresholds (e.g. energy_upper_bound) that don't affect data
% collection or storage, just for displaying the objects and passing them to subsequent workflows.
% This is so that we can retro-actively apply thresholds and display images without having to re-do
% the previous workflows. 
%
% 7) figure out how to make multi-channel images in matlab if that's possible. then add those
% possibilities to the visual output request list.  Run Image J inside of matlab so that we dont
% have to make the overlays by hand every time
%
% SAM 3/4/18
%
% SAM 7/18/18 updated 
