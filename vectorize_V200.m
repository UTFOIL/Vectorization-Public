function [ time_stamp, ROI_names ] = vectorize_V200( varargin )
%% Vectorize_V200 - Samuel Alexander Mihelic - Novemeber 8th, 2018                                  
% VECTORIZE( ) prompts the user at the command window for all required inputs.  It first asks 
%     whether to vectorize a new batch of images or continue with a previous batch.  A batch is a
%     group of images that the VECTORIZE function organizes together with two properties:
%
%         1) The same set of input parameters applies to every image in a batch. 
% 
%         2) The images in a batch are processed in parallel at each step of the vectorization. 
%            (see Methods below for descriptions of the four steps in the vectorization algorithm).
%
%       If the user continues with a previous batch, VECTORIZE prompts the user to select a previous
%       batch folder with data to recycle.
%
%       Alternatively, if the user starts a new batch, VECTORIZE prompts the user to select a folder
%       with some image file(s) to be vectorized.  It makes a new batch folder in a location
%       specified by the user.
% 
%     In either case, VECTORIZE prompts the user for a few logistical inputs: which vectorization
%     step(s) to execute, what previous data or settings (if any) to recycle, which visual(s) to
%     output (if any), and whether or not to open a graphical curator interface. It also prompts the
%     user for workflow-specific parameters: It displays imported parameters for review, and prompts
%     the user for any missing required parameters.  VECTORIZE writes any outputs to the batch
%     folder with a time stamp of the form YYMMDD_HHmmss.
% 
%   Conventions:  Greater values in the IMAGE_MATRIX correspond to greater vascular signal
%                 The IMAGE_MATRIX dimensions correspond to the physical dimensions y, x, and z
%                 (1,x,z) is the top  border of the x-y image at height z
%                 (y,1,z) is the left border of the x-y image at height z
%                 (y,x,1) is the x-y image nearest to the objective
%
%   Supported input image file types: .tif
%
% For in-line function calls that do not require manual interfacing (e.g. for writing wrapper
% functions or for keeping a concise record of VECTORIZE function calls in a script file), see the
% Optional Input, Logistical Parameters, and Workflow Specific Parameters Sections.
% 
% Note:  For more organizational/navigational control over this document in MATLAB:
%           1) open the Preferences Window                                       (HOME>>PREFERENCES)
%           2) enable Code Folding for Sections              (MATLAB>>Editor/Debugger>>Code Folding)
%           3) fold all of the sections in this document                                    (ctrl,+)
%
%% ------------------------------------------- Overview ------------------------------------------- 
%
% The purpose of the vectorization algorithm is to convert a grayscale, 3D image of vasculature (see
% Inputs section) to a vectorized model of the vascular components.  The output model consists of a
% list of many spherical objects with 3D-position (x,y,z), radius (r), and a contrast metric (c, see
% Methods section). These objects are vectors because each object is a 5-tuples of real numbers:
% [x,y,z,r,c].  The output vectors can then be rendered as a 2- or 3-dimensional image at any
% requested resolution, or it could be analyzed for statistical properties such as volume fraction
% or bifurcation density. With these objects in hand, many analyses are greatly simplified. Many
% such demonstrations and visualizations are automatically output (see Outputs section).
%
%% ---------------------------------------- Optional Input ---------------------------------------- 
%
% VECTORIZE( IMAGE_MATRIX ) vectorizes the numerical array IMAGE_MATRIX.  A batch folder is made 
%     in the OutputDirectory specified by the user. The user is prompted for the other logistical
%     and workflow specific parameters as in the VECTORIZE( ) call.
%
%   Supported IMAGE_MATRIX variable types: 3D array of doubles
% 
%   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   
%
% VECTORIZE( IMAGE_MATRICES ) vectorizes each IMAGE_MATRIX in the cell vector IMAGE_MATRICES.  The
%     outputs in the batch folder are numbered by the input order of the images in the cell vector.
% 
%   Supported IMAGE_MATRICES variable types: Cell vector of 3D array of doubles
%
%    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   
%
% VECTORIZE( FILE_NAME ) vectorizes the IMAGE_MATRIX(-CES) specified by the path(s) in FILE_NAME.
%
%   Supported FILE_NAME variable types: character vectors
%
%   FILE_NAME is an absolute or relative paths to current working folder. Wild card commands (i.e.
%       '*' or '**' ) in the FILE_NAME are also supported. For example:
%
%     VECTORIZE(    '*.tif'  ) vectorizes all .tif files in the current directory.
%
%     VECTORIZE([ '**', filesep, '*.tif' ]) vectorizes all .tif files in the current directory or any
%         subdirectory thereof.
% 
%  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
%
% VECTORIZE( FILE_NAMES ) vectorizes the IMAGE_MATRICES specified by the cell vector of FILE_NAMES.
%
%   Supported FILE_NAMES variable types: Cell vector of character vectors
%
%% -------------------------------------- Logistical Parameters ----------------------------------- 
%
% VECTORIZE( ..., NAME, VALUE ) 
%     uses the following NAME/VALUE pair assignments for logistical inputs:
%
% ------- NAME -------      -------------------------------- VALUE --------------------------------
%
% 'OutputDirectory'         Char or string specifying the folder that contains the output batch
%                           folder.  The default is to prompt the user at the command window.
% 
% 'NewBatch'                'prompt' - Prompts the user at the command window.
%                           'yes'    - Makes a new batch folder in the OutputDirectory folder. 
%                                      This is the only option if the user provides the Optional
%                                      Input.
%                           'no'     - Writes into an existing batch folder in the OutputDirectory.
%
% 'PreviousBatch'           'prompt' (default) - Prompts the user to input a previous batch folder.
%                           'none'             - Does not import any existing data or settings.  If
%                                                NewBatch is 'no', this is  not an  option.
%                           'yyMMdd-HHmmss'    - Imports from the batch_yyMMdd-HHmmss folder in the 
%                                                OutputDirectory.
%                           'recent'           - Imports from the most recent batch_* folder in the
%                                                OutputDirectory.
%
% 'PreviousWorkflow'        'prompt' (default) - Prompts the user to input a previous settings file.
%                           'none'             - Does not import any existing data or settings.  If
%                                                NewBatch is 'no', this is not an option.
%                           'yyMMdd-HHmmss'    - Imports from the workflow_yyMMdd-HHmmss folder in
%                                                the batch folder. 
%                           'recent'           - Imports from the most recent workflow_settings_* 
%                                                folder in the batch folder.
%
% 'StartWorkflow'           Note: These are listed in order:  Energy is the first process.
% 
%                           'prompt' (default) - Prompts the user at the command window.
%                           'none'             - Does not run any workflow steps.  The user may run
%                                                curatation or visual steps.
%                           'next'      - Starts the vectorization process at the   next   step for
%                                         the selected PreviousWorkflow.
%                           'energy'    - Starts the vectorization process at the  Energy  step.
%                                         This is the only option if NewBatch is 'yes'.
%                           'vertices'  - Starts the vectorization process at the Vertices step.
%                           'edges'     - Starts the vectorization process at the   Edges  step.
%                           'network'   - Starts the vectorization process at the  Network step.
%
% 'FinalWorkflow'           'prompt'(default) - Prompts the user at the command window.
%                           'none'            - Does not run any vectorization processes.  
%                                               This is the only option if StartWorkflow is 'none'
%                           'one'       - Ends the vectorization process at the StartWorkflow step.
%                           'energy'    - Ends the vectorization process at the    Energy     step.
%                           'vertices'  - Ends the vectorization process at the   Vertices    step.
%                           'edges'     - Ends the vectorization process at the     Edges     step.
%                           'network'   - Ends the vectorization process at the    Network    step.
%
% 'Visual'                  'none'                 - Does not write visual outputs.
%                           'original'             - Writes visuals for just the input images.
%                           'energy'               - Writes visuals for just the   Energy step.
%                           'vertices'             - Writes visuals for just the Vertices step.
%                           'edges'                - Writes visuals for just the    Edges step.
%                           'network'              - Writes visuals for just the  Network step.
%                           { ... }                - Writes visuals for just the    ...   steps.
%                           'productive' (default) - Writes visuals for just the workflow steps.
%                                                    being executed at this vectorization call.
%                           'all'                  - Writes visuals for all vectorization steps.
% 
% 'SpecialOutput'           'none'           - Does not create any special network outputs.
%                           'histograms'     - (defualt) shows strand, length statistic histograms
%                           'depth-stats'    - Shows depth-resolved statistics.
%                           'flow-field'     - Writes x, y, and z component .tif's of flow field.
%                           'depth'          - Shows vectors over raw with color-coded depth.               
%                           'strands'        - Shows vectors over raw with color-coded strands.
%                           'directions'     - Shows vectors over raw with color-coded direcions.
%                           '3D-strands'     - Shows 3D volume rendering with color-coded strands.
%                           'casX'           - Creates .casX equivalent representation of strands.
%                                              Format .casX is due to LPPD in Chicago.
%                           { ... }          - Creates ... special network outputs.
%                           'all'            - Creates all special network outputs.
% 
% 'VertexCuration'          'auto'             - All non-overlapping vertices are passed to edges. 
%                                                Chooses least energy vector upon volume conflict. 
%                           'manual' (default) - Prompts user with a graphical curation interface.
%                           'machine-manual'   - Applies neural network categorization and then
%                                                prompts user with a graphical interface.
%                           'machine-auto'     - Applies neural network categorization and then
%                                                all non-overlapping vertices are passed to edges.
%                                                Chooses best vector according to the neural network
%                                                upon volume conflict.
%
% 'EdgeCuration'            'auto'             - All edges are passed to network.
%                           'manual' (default) - Prompts user with a graphical curation interface.
%                           'machine-manual'   - Applies neural network categorization and then
%                                                prompts user with a graphical interface.
%                           'machine-auto'     - Applies neural network categorization.  All edges 
%                                                are passed to network.
% 
% 'NetworkPath'             'prompt
%                           'built-in' (default) - Built-in network ...
%                           'train'              - Trains new network from all curation files found
%                                                  in the training folder in the vectorization base 
%                                                  directory.
%                           'yyMMdd-HHmmss'      - Imports network trained at yyMMdd-HHmmss from the 
%                                                  network folder in the vectorization base
%                                                  directory.
%                           'recent'             - Imports network trained most recently from the 
%                                                  network folder in the vectorization base
%                                                  directory.
%
% 'Forgetful'               'none' (default)   - No intermediate data files will be deleted
%                           'original'         - Deletes intermediate copies of the input images
%                           'energy'           - Deletes intermediate Energy data files
%                           'both'             - Deletes both intermediate data files when done
%
% 'Presumptive'             false (default) - Prompts user for all required workflow-specific inputs
%                           true            - Assumes previous settings or default if no previous
%% ------------------------------- Workflow-Specific Parameters ----------------------------------- 
%
% VECTORIZE( ..., NAME, VALUE ) 
%     uses the NAME/VALUE pair assignments listed below to input workflow-specific parameters.  Each
%     parameter is listed below under the first workflow step that it modifies.
%   
%   Resolving conflicts between NAME/VALUE pairs and imported parameters from the PreviousWorkflow:
% 
%       If a NAME/VALUE pais is provided for a parameter that modifies a workflow step that is
%       upstream of the StartWorkflow, then that value is ignored and the value from the
%       PreviousWorkflow value will remain. A warning is produced.  This must occur because the
%       relevant workflow has already been executed and is not scheduled to run again, therefore the
%       parameters that modify it will not change.
%
%       If a NAME/VALUE pair is provided for a parameter that only modifies workflows that are equal
%       to or downstream of the StartWorkflow, then that value will overwrite any value that may
%       have been retrieved from the PreviousWorkflow.  The overwriting is displayed in the commaind
%       window.
%
%%  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - Energy -  -  -  -  -  -  -  -  -  -  -  -  -  
% ---------------- NAME ----------------  ------------------------- VALUE -------------------------
% 
% 'microns_per_voxel'                     Real, positive, three-element vector specifying the voxel
%                                         size in microns in y, x, and z dimensions, respectively.
%                                         Default: [ 1, 1, 1 ]
% 
% 'radius_of_smallest_vessel_in_microns'  Real, positive scalar specifying the radius of the
%                                         smallest vessel to be detected in microns.  Default: 1.5
% 
% 'radius_of_largest_vessel_in_microns'   Real, positive scalar specifying the radius of the largest 
%                                         vessel to be detected in microns.  Default: 50
% 
% 'approximating_PSF'                     Logical scalar specifying whether to approximate the PSF
%                                         using "Nonlinear Magic: Multiphoton Microscopy in the 
%                                         Biosciences" (Zipfel, W.R. et al.).  Default: true
% 
% 'sample_index_of_refraction'            Real, positive scalar specifying the index of refraction 
%                                         of the sample.  This parameter is only used if
%                                         approximating the PSF.  Default: 1.33
% 
% 'numerical_aperture'                    Real, positive scalar specifying the numerical aperture of
%                                         the microscope objective.  Default: 0.95
% 
% 'excitation_wavelength_in_microns'      Real, positive scalar specifying the excitation wavelength
%                                         of the laser in microns.  Default: 1.3
% 
% 'scales_per_octave'                     Real, positive scalar specifying the number of vessel 
%                                         sizes to detected per doubling of the radius cubed.  
%                                         Default: 1.5
% 
% 'max_voxels_per_node_energy'            Real, positive scalar specifying Default: 1e5
% 
% 'gaussian_to_ideal_ratio'               Real scalar between 0 and 1 inclusive specifying the
%                                         standard deviation of the Gaussian kernel per the total
%                                         object length for objects that are much larger than the
%                                         PSF. Default: 1
% 
% 'spherical_to_annular_ratio'            Real scalar between 0 and 1 inclusive specifying the
%                                         weighting factor of the spherical pulse over the combined
%                                         weights of spherical and annular pulses. This parameter is
%                                         only used if gaussian_to_ideal_ratio is strictly less than
%                                         unity.  Default: 1
%
%
%%   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   Edges   -  -  -  -  -  -  -  -  -  -  -  -  - 
% ---------------- NAME ----------------  ------------------------- VALUE -------------------------
% 
% 'max_edge_length_per_origin_radius'     Real, positive scalar specifying the maximum length of an
%                                         edge trace per the radius of the seed vertex. Default: 30
%
% 'number_of_edges_per_vertex'            Real, positive integer specifying the maximum number of
%                                         edge traces per seed vertex. Default: 4
%
%%  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   Network -  -  -  -  -  -  -  -  -  -  -  -  -  
% ---------------- NAME ----------------  ------------------------- VALUE -------------------------
% 
% 'sigma_strand_smoothing'                Real, non-negative integer specifying the standard
%                                         deviation of the Gaussian smoothing kernel per the radius
%                                         of the strand at every position along the strand vector.
%                                         Default: 1
%
%% ------------------------------------------- Methods -------------------------------------------- 
%
% Vectorization is accomplished in four steps:  (1) energy image formation, (2) vertex extraction,
% (3) edge extraction, and (4) network extraction.  The raw output is a superposition
% of possible models until it is segmented in some way.  The objects are assigned a contrast metric
% based on the values from the energy image, and thresholding them on this value provides direct
% control over the sensitivity and specificity of the vectorization.  Alternatively, the graphical
% curation interface provides a platform for manual segmentation such as local threshold selection
% or point-and-click object selection.
% 
% 1: Energy Image Formation
%     Multi-scale (at many pre-defined sizes) gradient and curvature information from the original
%     3D image is combined to form a 4-dimensional, multi-scale, centerline-enhanced, image, known
%     as the energy image.  Ideally, the voxels with the lowest energy value the energy image will
%     be the most likely to be centerline voxels for vessels in the pre-defined size range.  This
%     can be visually verified by inspecting the energy*.tif visual output in the visual output
%     directory.  The energy image should be very negative right at the vessel centerlines and close
%     to zero or positive elsewhere.  The value of the size image at a given voxel shows the index
%     of the pre-defined sizes that is most likely assuming a vessel is centered at that voxel.
% 
% 2: Vertex Extraction
%     Vertices are extracted as local minima in the 4D energy image (with associated x, y, z,
%     radius, and energy value).  This method was inspired by the first part of the SIFT algorithm
%     (David Lowe, International Journal of Computer Vision, 2004)).  Ideally, local minima of the
%     energy image correspond to voxels that are locally the most likely to be along a vessel
%     centerline. The size coordinate is also required to be at a local energy minimum.  In theory,
%     the vertices are ordered from most to least likely to exist by the energy values at their
%     locations.
% 
% 3: Edge Extraction
%     Edges are extracted as voxel to voxel random walks through the (min. projected to 3D) energy
%     image.  Therefore edges are lists of spherical objects like vertices.  Edge trajectories seek
%     lower energy values and are forced to connect exactly two vertices.  The trajectories between
%     vertices are in theory ordered from most to least likely to exist by their mean energy values.
% 
% 4: Network Extraction
%     Strands are defined as the sequences of non-branching edges (single random color in the
%     colored strands image).  Strands are found by counting the number of edges in the adjacency
%     matrix of the vertices.  Strands are the connected components of the adjacency matrix that
%     only includes vertices with two edges.  With the strands of the network in hand, we
%     equivalently know where the bifurcations in the network are.  Network information unlocks many
%     doors.  For instance, we can smooth the positions and sizes of the extracted vectors along
%     their strands and approximate local blood flow fields.
%
%% ------------------------------------------- Outputs -------------------------------------------- 
%
% Standard network ouptut format(s): .casx
%
% The vectorization output is the set of 3-space locations and radii of all vessels as well as 
% their connectivity. This information is stored in several matlab variables in the network output 
% file in the vector directory of the batch_* directory output of the Vectorize function. This 
% vector information is also available in the casx file format if the user selects the 
% 'SpecialOutputs'/'casx' name/value pair input. The .casx standard file is due to LPPD at 
% University of Illinois at Chicago, Department of Bioengineering (https://lppd.bioe.uic.edu/)
%
% TIME_STAMP = VECTORIZE( ... )
%     returns the TIME_STAMP that was assigned to any new data or settings output
%
% [ TIME_STAMP, ROI_NAMES ] = VECTORIZE( ... )
%     also returns the FILE_NAMES that were used as input or the names assigned to the IMAGE_MATRICES in
%     the order that they were passed to VECTORIZE.
%
%% -------------------------------------------  Notes  -------------------------------------------- 
% VECTORIZE( ) prompts the user at the command window for all required inputs.  It first asks 
%     whether to vectorize a new batch of images or continue with a previous batch.  A batch is a
%     group of images that the VECTORIZE function organizes together with two properties:
%
%         1) The same set of input parameters applies to every image in a batch. 
% 
%         2) The images in a batch are processed in parallel at each step of the vectorization. 
%            (see Methods below for descriptions of the four steps in the vectorization algorithm).
%
%       If the user continues with a previous batch, VECTORIZE prompts the user to select a previous
%       batch folder with data to recycle.
%
%       Alternatively, if the user starts a new batch, VECTORIZE prompts the user to select a folder
%       with some image file(s) to be vectorized.  It makes a new batch folder in a location
%       specified by the user.
% 
%     In either case, VECTORIZE prompts the user for a few logistical inputs: which vectorization
%     step(s) to execute, what previous data or settings (if any) to recycle, which visual(s) to
%     output (if any), and whether or not to open a graphical curator interface. It also prompts the
%     user for workflow-specific parameters: It displays imported parameters for review, and prompts
%     the user for any missing required parameters.  VECTORIZE writes any outputs to the batch
%     folder with a time stamp of the form YYMMDD_HHmmss.
% 
%   Conventions:  Greater values in the IMAGE_MATRIX correspond to greater vascular signal
%                 The IMAGE_MATRIX dimensions correspond to the physical dimensions y, x, and z
%                 (1,x,z) is the top  border of the x-y image at height z
%                 (y,1,z) is the left border of the x-y image at height z
%                 (y,x,1) is the x-y image nearest to the objective
%
%   Supported input image file types: .tif
%
% For in-line function calls that do not require manual interfacing (e.g. for writing wrapper
% functions or for keeping a concise record of VECTORIZE function calls in a script file), see the
% Optional Input, Logistical Parameters, and Workflow Specific Parameters Sections.
% 
% Note:  For more organizational/navigational control over this document in MATLAB:
%           1) open the Preferences Window                                       (HOME>>PREFERENCES)
%           2) enable Code Folding for Sections              (MATLAB>>Editor/Debugger>>Code Folding)
%           3) fold all of the sections in this document                                    (ctrl,+)
%
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
%
% V200:  Attempting to create an end-user version of the vectorization software for a general
% neuroscientist that requires minimal directory structure and is intuitive to use. Moving from a
% script to a function.  SAM 11/8/18

%% 0: Input Parser                                                                                  

% time stamp
time_stamp = char( datetime('now', 'TimeZone', 'local', 'Format', 'yyMMdd-HHmmss' )); 

q = inputParser ;

q.FunctionName = 'vectorize' ;

q.KeepUnmatched = true ;
%% ---------------------------------------- Optional Input ---------------------------------------- 
  
% IF optional input is provided (the total number of inputs will be odd).
if mod( length( varargin ), 2 ) % odd number of inputs
        
    p = inputParser ;

    p.FunctionName = 'vectorize' ;

    p.KeepUnmatched = true ;
    
    addRequired(  p, 'OptionalInput' )
    
    parse( p, varargin{ : });
            
    NewBatch_values  = { 'yes' };
    NewBatch_default =   'yes'  ;
    
    optional_input_provided = true  ;
                
else % ELSE optional input not provided
    
    p = inputParser ;

    p.FunctionName = 'vectorize' ;

    p.KeepUnmatched = true ;    
            
    NewBatch_values  = { 'prompt', 'yes', 'no' };
    NewBatch_default =   'prompt'               ;
    
    optional_input_provided = false ;    
    
end % IF no optional input provided

    function [ validation_flag, ROI_names ] = validate_the_OptionalInput( x )

        validation_flag = true;

        % IF cell array
        if iscell( x )

            input_index_range = 1 : numel( x );

            validation_type = cell( size( x ));

            for input_index = input_index_range

                validation_type{ input_index } = validate_one_OptionalInput( x{ input_index }, ' cell array entry ' );

            end % FOR input index
        else % ELSE not a cell

            validation_type = { validate_one_OptionalInput( x, '' )}; % cell scalar
            
            input_index_range = 1 ;
            
            x = { x };
            
        end % IF cell array
               
        number_of_inputs_total = 0 ; 
        
        number_of_matrices = 0 ;
                
        ROI_names         = { };
        
%         matrix_file_name = [ current_directory, filesep, 'ans_' ];
        
        for input_index = input_index_range
        
            switch validation_type{ input_index }

                case 'path'
                    
                    input_file_listing = dir( x{ input_index });
                    
                    number_of_input_files = length( input_file_listing );
                    
                    if number_of_input_files == 0
                        
                        error([ 'No files matched the input file: ', x{ input_index }])
                    
                    end
                    
                    input_file_cell  = struct2cell( input_file_listing );

                    input_file_names = input_file_cell( 1, : );    
                    
                    input_file_names_temp = input_file_names ;
                                                                                                    
                    for file_index = 1 : number_of_input_files
                        
                        number_of_inputs_total = number_of_inputs_total + 1 ;
                                                
                        input_image_path = [ input_file_cell{ 2, file_index }, filesep, input_file_names{ file_index }]; % TIF file path
                        
                        ROI_names{ number_of_inputs_total } = [ '_', input_file_names{ file_index }( 1 : end - 4 )];
                        
                        file_extension = input_file_names{ file_index }( end - 3 : end );     
                        
                        if ~ strcmp( file_extension, '.tif' ) 
                            
                            error([ 'File path (', input_image_path, '), in OptionalInput pointed to file without ".tif" extension. Input file type must be .tif' ])
    
                        end
                        
                        if any( strcmp( input_file_names_temp{ file_index }, input_file_names_temp([ 1 : file_index - 1, file_index + 1 : end ])))
                            
                            error( 'Multiple inputs with the same name were referenced in .tif files in the OptionalInput.  Input file names must be unique.' )
                            
                        end
                                                                        
                        original_data = tif2mat( input_image_path );

                        size_of_image = size( original_data );
                        
                        path_to_original_data = fullfile( root_directory, [ 'batch_', time_stamp ], 'data', [ 'original',  ROI_names{ number_of_inputs_total }]); % h5  file path
                        
                        original_data_type = class( original_data );

                        h5create( path_to_original_data, '/d', size_of_image, 'Datatype', original_data_type );

                        mat2h5( path_to_original_data, original_data );

                    end
                    
                case 'matrix'

                    number_of_inputs_total = number_of_inputs_total + 1 ;    
                    
                    number_of_matrices = number_of_matrices + 1 ;
                    
                    ROI_names{ number_of_inputs_total } = [ '_ans_', num2str( number_of_matrices )];                    
                    
                    path_to_original_data = fullfile( root_directory, [ 'batch_', time_stamp ], 'data', [ 'original',  ROI_names{ number_of_inputs_total }]); % h5  file path
                                        
                    h5create( path_to_original_data, '/d', size( x{ input_index }), 'Datatype', class( x{ input_index }));

                    mat2h5( path_to_original_data, x{ input_index });
                                        
            end
        end % FOR ROI index

        function [ validation_type ] = validate_one_OptionalInput( x, phrase )

            if ischar( x ) || isstring( x )
                
                validation_type = 'path' ;
                                
            else % is not char or string
                if isnumeric( x )
                    
                    if length( size( x )) == 3, validation_type = 'matrix' ;

                    else, error([ 'OptionalInput', phrase, ' was numeric, but was not three-dimensional' ])

                    end
                else, error([ 'OptionalInput', phrase, 'was neither type char, string, nor numeric' ])

                end
            end % is char or string
        end % FUNCTION VALIDATE_ONE_OPTIONALINPUT
    end % FUNCTION VALIDATE_THE_OPTIONALINPUT
%% -------------------------------------- Logistical Parameters -----------------------------------  
    %%   -  -  -  -  -  -  -  -  -  -  -  -  -  - output directory -  -  -  -  -  -  -  -  -  -  -  

addParameter( p, 'OutputDirectory', '', @( x ) isstring( x ) || ischar( x ));

parse( p, varargin{ : });

root_directory = p.Results.OutputDirectory ;

if isempty( root_directory )
    
    root_directory = input( 'Please enter an absolute path to an output directory [current directory]:', 's' );
    
    % default text response    
    if isempty( root_directory ), root_directory = pwd ; end
    
end

root_directory = [ root_directory, filesep ];

disp([ 'OutputDirectory: ', root_directory ]);
    %% -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   new batch -  -  -  -  -  -  -  -  -  -  -  -   

batch_listing = dir([ root_directory, 'batch_*-*' ]); 

there_is_no_previous_batch = isempty( batch_listing );
 
addParameter( p, 'NewBatch', NewBatch_default, @( x ) any( validatestring( x, NewBatch_values, 'vectorize', 'NewBatch' )));

% if ~ isempty( varargin

parse( p, varargin{ : });

new_batch = p.Results.NewBatch ;

if strcmp( validatestring( new_batch, NewBatch_values ), 'prompt' )
    
    if there_is_no_previous_batch
                
        disp( 'No previous batch folders found in the chosen output directory. Creating a new one...' )
        
        new_batch = 'yes' ;
        
    else % there is a previous batch

        disp( 'Previous batch folder(s) were found in the chosen output directory.' )
        
        disp( 'You may import one of these or create a new folder.' )
        
        new_batch = input( 'Would you like to create a new batch folder in the OutputDirectory [no]?: ', 's' );

    end
    
    % default text response        
    if isempty( new_batch ), new_batch = 'no' ; end
    
    NewBatch_values = { 'yes', 'no' };
    
end

batch_listing_cell        = struct2cell( batch_listing );

batch_listing_time_stamps = cellfun( @( x ) x( 7 : 19 ), batch_listing_cell( 1, : ), 'UniformOutput', false );        

new_batch = validatestring( new_batch, NewBatch_values );

if strcmp( new_batch, 'yes' )
            
    new_batch = true ; 
    
%     PreviousBatch_values  = [{ 'none', 'prompt', 'recent' }, batch_listing_time_stamps ];
%     PreviousBatch_default =            'prompt'                                             ;

    PreviousBatch_values  = { 'none' };
    PreviousBatch_default =   'none'  ;

    % This is the folder that keeps the vectorization progress of one or more vectorization attempts
    % of a set of batch-processed ROI's
    
    batch_handle = [ 'batch_', time_stamp ];
    
    [         batch_directory, ...
               data_directory, ...
        visual_data_directory, ...
             vector_directory, ...
      visual_vector_directory, ...
           curation_directory, ...
           settings_directory  ]    = get_directories( root_directory, batch_handle );
       
    disp([ 'Creating a new batch folder (for parallel image processing) in the OutputDirectory: ' batch_directory ])

    [ ~, ~ ] = mkdir(         batch_directory );         
    
    [ ~, ~ ] = mkdir(          data_directory );
    [ ~, ~ ] = mkdir(   visual_data_directory );    
    [ ~, ~ ] = mkdir(        vector_directory );
    [ ~, ~ ] = mkdir( visual_vector_directory );
    [ ~, ~ ] = mkdir(      curation_directory );         
         
    if ~ optional_input_provided
        
        [ file, path ] = uigetfile( '*.tif', 'Please select one or more input .tif file.', 'MultiSelect', 'on');
        
        if isnumeric( path ) || isnumeric( path ), error( 'No paths were collected by the folder selection interface.  Please retry selection.' ); end
        
        optional_input = fullfile( path, file );
            
    else % OptionalInput provided
        
        optional_input = p.Results.OptionalInput ;
                
    end % IF Optional Input Not Provided
    
    disp( 'Attempting to write an h5 files for each input into the "Data" subdirectory of the OutputDirectory... ' )
    
    [ ~, ROI_names ] = validate_the_OptionalInput( optional_input );   
    
    number_of_ROIs = length( ROI_names );   

    ROI_index_range = 1 : number_of_ROIs ;  
    
    if number_of_ROIs > 1
    
        names_list = '' ;
        
        for ROI_index = ROI_index_range( 1 : end - 1 )

            names_list = [ names_list, ROI_names{ ROI_index }, ', ' ];

        end
        
        names_list = [ names_list, 'and ', ROI_names{ end }( 2 : end ), '.' ];
        
        disp([ 'Input names were set to: ', names_list ])            
        
    else % one single ROI
        
        disp([ 'Input name was set to: ',  ROI_names{  1  }( 2 : end ), '.' ])
        
        if exist( 'path', 'var' ), disp([ 'Input file was located at path: ', path ]), end
        
    end
    
    % This is where all the settings for this batch will be saved.
    [ ~, ~ ] = mkdir( settings_directory );
    
    % This file holds the directory tree that is specific to this batch.
    batch_settings = [ settings_directory, 'batch' ]; 
    
    save( batch_settings      , ...
              'optional_input', ...
               'ROI_names'    , ...
  ...            'data_directory', ...
 ...      'visual_data_directory', ...
  ...          'vector_directory', ...
 ...    'visual_vector_directory', ...
  ...        'curation_directory', ...
  ...        'settings_directory', ...
          'ROI_index_range'    );
          
else % ELSE old batch
        
    there_exists_a_previous_batch = ~ isempty( batch_listing );
    
    if there_exists_a_previous_batch
        
        new_batch = false ; 
        
        PreviousBatch_values  = [{ 'prompt', 'recent' }, batch_listing_time_stamps ];
        PreviousBatch_default =    'prompt'                                             ;
        
        disp( 'Importing settings and data from a previous batch folder in the OutputDirectory.' )
        
    else % ELSE no previous batch exists

        error( 'No previous batches were found in the OutputDirectory. Previous batch folders need to be in the OutputDirectory to import data or settings.' )        

    end % IF previous batch exists
end % IF new batch
    %%  -  -  -  -  -  -  -  -  -  -  -  -  -   previous batch   -  -  -  -  -  -  -  -  -  -  -  - 

addParameter( p, 'PreviousBatch', PreviousBatch_default, @( x ) any( validatestring( x, PreviousBatch_values, 'vectorize', 'PreviousBatch' )));

parse( p, varargin{ : });    

previous_batch = validatestring( p.Results.PreviousBatch, PreviousBatch_values );

if strcmp( previous_batch, 'prompt' )
    
    number_of_previous_batches = length( batch_listing_time_stamps );
        
    disp( 'Previous batch time-stamps:' )
    
    for batch_index = 1 : number_of_previous_batches
       
        disp([ '          ', batch_listing_time_stamps{ batch_index }])
        
    end
            
    previous_batch = input( 'Please select a PreviousBatch from the OutputDirectory [recent]: ', 's' );
    
    % default text response        
    if isempty( previous_batch )
        
        previous_batch = 'recent' ; 
    
    else % validate response against possible choices
        
        previous_batch = validatestring( previous_batch, [{ 'recent' }, batch_listing_time_stamps ]);
        
    end
end

switch previous_batch
    
    case 'none'
        
        previous_batch_handle = 'none' ;
                        
    case 'recent'
        
        most_recent_previous_batch_time_stamp = batch_listing( end ).name( 7 : 19 );             
        
        previous_batch_handle = [ 'batch_', most_recent_previous_batch_time_stamp ];
                
    otherwise % specific time-stamp referenced by user

        previous_batch_handle = [ 'batch_', previous_batch ];

end

if strcmp( previous_batch_handle, 'none' )
    
    disp( 'Not loading any previous data or settings from a previous batch' )
    
    previous_batch_directory = '' ;
    
    PreviousWorkflow_values  = { 'none' };
    
    PreviousWorkflow_default =   'none'  ;

else
    
    % load an old batch to continue working from    
    disp([ 'Loading previous settings and data from    ', previous_batch_handle, ' in the OutputDirectory...' ])
    
    previous_batch_directory = [ root_directory, previous_batch_handle, filesep ];
    
    previous_settings_directory = [ previous_batch_directory, 'settings', filesep ];
    
    previous_batch_settings     = [ previous_settings_directory, 'batch' ];   

    load( previous_batch_settings )
    
    [                       ~, ...
               data_directory, ...
        visual_data_directory, ...
             vector_directory, ...
      visual_vector_directory, ...
           curation_directory, ...
           settings_directory  ]    = get_directories( root_directory, previous_batch_handle );
    
    % find most recent workflow time stamp in the most recent batch directory
    workflow_listing = dir([ previous_batch_directory, 'settings', filesep, 'workflow_*.mat' ]); 

    there_exists_a_previous_workflow = ~ isempty( workflow_listing );

    if there_exists_a_previous_workflow

        workflow_listing_cell        = struct2cell( workflow_listing );

        workflow_listing_time_stamps = cellfun( @( x ) x( 10 : 22 ), workflow_listing_cell( 1, : ), 'UniformOutput', false );        
        
        PreviousWorkflow_values  = [{ 'prompt', 'recent' }, workflow_listing_time_stamps ];
        
        PreviousWorkflow_default =    'prompt' ;
        
    else % ELSE no previous workflow exists

        error('No previous workflows were found in the PreviousBatch directory. Previous workflow .mat files need to be in the PreviousBatch directory to use previous data or settings.')        

    end % IF previous workflow exists
end
    %%   -  -  -  -  -  -  -  -  -  -  -  -  -  previous workflow  -  -  -  -  -  -  -  -  -  -  -  
    
addParameter( p, 'PreviousWorkflow', PreviousWorkflow_default, @( x ) any( validatestring( x, PreviousWorkflow_values, 'vectorize', 'PreviousWorkflow' )));

parse( p, varargin{ : });

previous_workflow = validatestring( p.Results.PreviousWorkflow, PreviousWorkflow_values );

if strcmp( previous_workflow, 'prompt' )
    
    number_of_previous_workflows = length( workflow_listing_time_stamps );
    
    disp( 'Previous workflow time-stamps:' )
    
    for workflow_index = 1 : number_of_previous_workflows
              
        disp([ '          ', workflow_listing_time_stamps{ workflow_index }])
        
    end
            
    previous_workflow = input( 'Please select a PreviousWorkflow .mat file from the PreviousBatch directory [recent]: ', 's' );
    
    % default text response        
    if isempty( previous_workflow )
        
        previous_workflow = 'recent' ; 
    
    else % validate response against possible choices
        
        previous_workflow = validatestring( previous_workflow, [{ 'recent' }, workflow_listing_time_stamps ]);
        
    end
end

switch previous_workflow
                    
    case 'recent'
        
        most_recent_previous_workflow_time_stamp = workflow_listing( end ).name( 10 : 22 );

        previous_workflow_handle = [ 'workflow_', most_recent_previous_workflow_time_stamp ];
        
    otherwise % specific time-stamp referenced by user

        previous_workflow_handle = [ 'workflow_', previous_workflow ];

end

% load or instantiate the production_times variable
if ~ strcmp( previous_workflow, 'none' )
        
    disp([ 'Loading previous settings and data from ', previous_workflow_handle, '...' ])
    
    % load an old workflow schedule to get the production times to handle the old data
    path_to_previous_workflow_settings = [ previous_settings_directory, previous_workflow_handle ];

    load( path_to_previous_workflow_settings, 'production_times', 'attempted_production_times' )  
            
else
    
    % initialize the production times variable
    production_times = cell( 1, 5 );
    
    production_times{ 1 } = time_stamp ;
    
    production_times( 2 : 5 ) = { 'yyMMdd-HHmmss' };
    
    attempted_production_times = production_times ;

end
    %% -  -  -  -  -  -  -  -  -  -  -  -  -  -  - start workflow -  -  -  -  -  -  -  -  -  -  -   

% workflow steps in order for reference:
workflows = { 'original', ... % 1
              'energy'  , ... % 2
              'vertices', ... % 3
              'edges'   , ... % 4
              'network' , };  % 5

% if new_batch
%     
%     StartWorkflow_values  = [ workflows( 2 ), { 'none' }];
%     StartWorkflow_default = 'prompt' ;
%     
%     time_index = 2 ;
% 
% else
    
% for time_index = 1 : 5
% 
%     if strcmp( production_times{ time_index }, 'yyMMdd-HHmmss' ), break, end
% 
% end

time_index = find( strcmp( 'yyMMdd-HHmmss', production_times ), 1 ); 

if isempty( time_index ), time_index = 5 ; end

StartWorkflow_values  = [{ 'prompt', 'next' }, workflows( 2 : time_index ), { 'none' }]; 
StartWorkflow_default =    'prompt' ;    

% end

% VertexCuration_values = { 'none', 'auto', 'manual' };
%   EdgeCuration_values = { 'none', 'auto', 'manual' };
%   
%      Forgetful_values = { 'none', 'original', workflows{ 1 }, 'both' };  
  
addParameter( p, 'StartWorkflow', StartWorkflow_default, @( x ) any( validatestring( x, StartWorkflow_values, 'vectorize', 'StartWorkflow' )))

% addParameter( p, 'Visual'          , 'productive', @validate_visual )
% 
% addParameter( p, 'Forgetful'       , 'none'   , @( x ) any( validatestring( x, Forgetful_values,      'vectorize', 'Forgetful'      )))
% addParameter( p, 'VertexCuration'  , 'auto'   , @( x ) any( validatestring( x, VertexCuration_values, 'vectorize', 'VertexCuration' )))
% addParameter( p, 'EdgeCuration'    , 'manual' , @( x ) any( validatestring( x, EdgeCuration_values,   'vectorize',   'EdgeCuration' )))

parse( p, varargin{ : });

start_workflow = validatestring( p.Results.StartWorkflow, StartWorkflow_values );

if strcmp( start_workflow, 'prompt' )
        
    disp( 'Possible StartWorkflow values:' )
        
    for StartWorkflow_value = StartWorkflow_values( 2 : end )

        disp([ '          ', StartWorkflow_value{ 1 }])
        
    end
    
    start_workflow = input( 'Please select a StartWorkflow value [next]: ', 's' );
    
    % default text response
    if isempty( start_workflow )
        
        start_workflow = 'next' ; 
    
    else % validate response against possible choices
        
        start_workflow = validatestring( start_workflow, StartWorkflow_values( 2 : end ));
        
    end
end

switch start_workflow
    
    case 'next'
        
        ordinate_of_start_workflow = time_index ;
        
        if ordinate_of_start_workflow == 5, start_workflow = workflows{ 5 }; end
        
    case workflows{ 2 }
        
        ordinate_of_start_workflow = 2 ;
        
    case workflows{ 3 }
        
        ordinate_of_start_workflow = 3 ;
        
    case workflows{ 4 }
        
        ordinate_of_start_workflow = 4 ;
        
    case workflows{ 5 }
        
        ordinate_of_start_workflow = 5 ;
        
end

switch start_workflow
    
    case 'none'
        
        ordinate_of_start_workflow = 1 ; % must set to something
        
        FinalWorkflow_values  = { 'none' };
        FinalWorkflow_default =   'none'  ;
        
        disp( 'no StartWorkflow selected, visual and curation processes are still possible.' )
        
    case 'network'
        
        ordinate_of_start_workflow = 5 ;
        
        FinalWorkflow_values  = { 'network', 'one' };
        FinalWorkflow_default =   'network' ;
        
        disp( 'StartWorkflow selected: network' )
        
    otherwise
        
        FinalWorkflow_values  = [{ 'prompt', 'one' }, workflows( ordinate_of_start_workflow : end )];
        FinalWorkflow_default =    'prompt' ;
        
        disp([ 'StartWorkflow selected: ', workflows{ ordinate_of_start_workflow }])
        
end
    %%  -  -  -  -  -  -  -  -  -  -  -  -  -  -   final workflow   -  -  -  -  -  -  -  -  -  -  - 

addParameter( p, 'FinalWorkflow', FinalWorkflow_default, @( x ) any( validatestring( x, FinalWorkflow_values, 'vectorize', 'FinalWorkflow' )))

parse( p, varargin{ : });

final_workflow = validatestring( p.Results.FinalWorkflow, FinalWorkflow_values );

if strcmp( final_workflow, 'prompt' )
        
    disp( 'Possible FinalWorkflow values:' )
        
    for FinalWorkflow_value = FinalWorkflow_values( 2 : end )

        disp([ '          ', FinalWorkflow_value{ 1 }])
        
    end
    
    final_workflow = input( 'Please select a FinalWorkflow value [network]: ', 's' );
    
    % default text response
    if isempty( final_workflow )
        
        final_workflow = workflows{ 5 }; 
    
    else % validate response against possible choices
        
        final_workflow = validatestring( final_workflow, FinalWorkflow_values( 2 : end ));
        
    end
end

switch final_workflow
    
    case 'none'
        
        ordinate_of_final_workflow = 0 ;
        
    case 'one'
        
        ordinate_of_final_workflow = ordinate_of_start_workflow ;
        
    case workflows{ 2 }
        
        ordinate_of_final_workflow = 2 ;
        
    case workflows{ 3 }
        
        ordinate_of_final_workflow = 3 ;
        
    case workflows{ 4 }
        
        ordinate_of_final_workflow = 4 ;
        
    case workflows{ 5 }
        
        ordinate_of_final_workflow = 5 ;
        
end

if ~ strcmp( final_workflow, 'none' )
        
    disp([ 'FinalWorkflow selected: ', workflows{ ordinate_of_final_workflow }])

end

 previous_production_times = attempted_production_times ;
attempted_production_times =           production_times ;

productive_indices = ordinate_of_start_workflow : ordinate_of_final_workflow ;

production_times( productive_indices ) = { 'yyMMdd-HHmmss' };

production_times{ 1 } = attempted_production_times{ 1 };

productive = zeros( 5, 1, 'logical' );  

productive( productive_indices ) = true ;

if any( productive ), production_times( productive_indices( 1 ) : end ) = { 'yyMMdd-HHmmss' }; end

attempted_production_times( productive_indices ) = { time_stamp };
    %% -  -  -  -  -  -  -  -  -  -  -  -  -   forgetful workflows   -  -  -  -  -  -  -  -  -  -   
Forgetful_values = { 'none', 'original', workflows{ 1 : 2 }, 'both' };  
    
addParameter( p, 'Forgetful', 'none', @( x ) any( validatestring( x, Forgetful_values, 'vectorize', 'Forgetful' )))

parse( p, varargin{ : });

switch p.Results.Forgetful
    
    case 'none'
        
        forgetful =                zeros( 5, 1, 'logical' );
        
    case workflows{ 1 }
        
        forgetful = [ true;        zeros( 4, 1, 'logical' )];
        
    case workflows{ 2 }
        
        forgetful = [ false; true; zeros( 3, 1, 'logical' )];
        
    case 'both'
        
        forgetful = [ true;  true; zeros( 3, 1, 'logical' )];        
        
end
    %%  -  -  -  -  -  -  -  -  -  -  -  -  -   curation workflows  -  -  -  -  -  -  -  -  -  -  - 
    
VertexCuration_values = { 'none', 'auto', 'manual', 'machine-manual', 'machine-auto' };
  EdgeCuration_values = { 'none', 'auto', 'manual', 'machine-manual', 'machine-auto' };

% 'none' is not an option if running the workflow to be curated.  'auto' is minimum workflow
if productive( 3 ), VertexCuration_default = 'manual'; VertexCuration_values( 1 ) = [ ]; else, VertexCuration_default = 'none'; end
if productive( 4 ),   EdgeCuration_default = 'manual';   EdgeCuration_values( 1 ) = [ ]; else,   EdgeCuration_default = 'none'; end

addParameter( p, 'VertexCuration', VertexCuration_default, @( x ) any( validatestring( x, VertexCuration_values, 'vectorize', 'VertexCuration' )))
addParameter( p,   'EdgeCuration',   EdgeCuration_default, @( x ) any( validatestring( x,   EdgeCuration_values, 'vectorize',   'EdgeCuration' )))

parse( p, varargin{ : }); 

vertex_curation = validatestring( p.Results.VertexCuration, VertexCuration_values );
  edge_curation = validatestring( p.Results.EdgeCuration  ,   EdgeCuration_values );
    %%   -  -  -  -  -  -  -  -  -  -  -  -  -  - visual workflows -  -  -  -  -  -  -  -  -  -  -  

addParameter( p, 'Visual', 'productive', @validate_visual )

parse( p, varargin{ : });    
    
[ ~, visual ] = validate_visual( p.Results.Visual );
    %%   -  -  -  -  -  -  -  -  -  -  -  -  -  -  presumptive  -  -  -  -  -  -  -  -  -  -  -  -  
    
addParameter( p, 'Presumptive', false, @islogical )

parse( p, varargin{ : });

presumptive = p.Results.Presumptive ;
    %% -  -  -  -  -  -  -  -  -  -  -  -  -  -  special outputs  -  -  -  -  -  -  -  -  -  -  -   

if visual( 5 ), special_output_default = 'histograms' ; else special_output_default = 'none' ; end
    
addParameter( p, 'SpecialOutput', special_output_default, @validate_special_output )

parse( p, varargin{ : });    
    
[ ~, special_output ] = validate_special_output( p.Results.SpecialOutput );
    %%  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  saving   -  -  -  -  -  -  -  -  -  -  -  -  - 
    
workflow_handle = [ 'workflow_', time_stamp ];

path_to_workflow_settings = [ settings_directory, workflow_handle ];    

save( path_to_workflow_settings, ...
      'workflows'              , ...
      'productive'             , ...
             'production_times', ...
   'attempted_production_times', ...
      'visual'                 , ...
      'special_output'         , ...
      'forgetful'              , ...
      'vertex_curation'        , ...
      'edge_curation'          , ...
      'presumptive'              );
  
%% ------------------------------- Workflow-Specific Parameters ----------------------------------- 

inputs_required     =     find( productive            )'         ;

if any( productive )

    inputs_inherited    = 1 : find( productive, 1         )  - 1     ;
    inputs_not_required =     find( productive, 1, 'last' )      : 5 ;

else
    
    inputs_inherited    = 1 : max([ find( visual, 1, 'last' );               ...
                                    3 * max(strcmp( vertex_curation, {'manual','machine-manual','machine-auto', 'auto' })); ...
                                    4 * max(strcmp(   edge_curation, {'manual','machine-manual','machine-auto', 'auto' })); ...
                                    5 * any( special_output )]);
    
    % no visual or curation specified                                
	if isempty( inputs_inherited )
        
        warning( 'Please specify a Visual, SpecialOutput, VertexCuration, or EdgeCuration name/value pair to perform a task.' )
        
%         return
        
    else
        
        inputs_not_required = inputs_inherited( end ) + 1 : 5 ;
    
    end
end

original_data_handle = 'original' ;

  previous_energy_handle = [  'energy_', previous_production_times{ 2 }];
previous_vertices_handle = ['vertices_', previous_production_times{ 3 }];
   previous_edges_handle = [   'edges_', previous_production_times{ 4 }];
 previous_network_handle = [ 'network_', previous_production_times{ 5 }];

           energy_handle = [  'energy_', attempted_production_times{ 2 }];
         vertices_handle = ['vertices_', attempted_production_times{ 3 }];
            edges_handle = [   'edges_', attempted_production_times{ 4 }];
          network_handle = [ 'network_', attempted_production_times{ 5 }];

previous_path_to_energy_settings   = [ settings_directory,   previous_energy_handle ];
previous_path_to_vertices_settings = [ settings_directory, previous_vertices_handle ];
previous_path_to_edges_settings    = [ settings_directory,    previous_edges_handle ];
previous_path_to_network_settings  = [ settings_directory,  previous_network_handle ];          
          
         path_to_energy_settings   = [ settings_directory,            energy_handle ];
         path_to_vertices_settings = [ settings_directory,          vertices_handle ];
         path_to_edges_settings    = [ settings_directory,             edges_handle ];
         path_to_network_settings  = [ settings_directory,           network_handle ];
    %% -  -  -  -  -  -  -  -  -  -  -  -  -   Inherited Inputs   -  -  -  -  -  -  -  -  -  -  -   

for workflow_index = inputs_inherited
    switch workflow_index
        case 2
        %%   -  -  -  -  -  -  -  -  -  -  -  - Energy -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
            excitation_wavelength = []; % To avoid 'attempted to add variable to a static workplace' error
            symmetry_ratio_factor = [];
        
            try load( path_to_energy_settings )
                                            
            catch
                                
                disp([ 'The previous energy settings file was missing from ', path_to_energy_settings ])
%                 error(  'Vectorization must be re-run from the energy step.' )
                            
            end
        case 3
        %%  -  -  -  -  -  -  -  -  -  -  -  - Vertices  -  -  -  -  -  -  -  -  -  -  -  -  -  -   

            try load( path_to_vertices_settings )
                                            
            catch
                                
                disp([ 'The previous vertices settings file was missing: ', path_to_vertices_settings ])
%                 error(  'Vectorization must be re-run from the vertex step.' )
                            
            end
        case 4
        %% -  -  -  -  -  -  -  -  -  -  -  -  - Edges  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

            try load( path_to_edges_settings )
                                            
            catch
                                
                disp([ 'The previous edges settings file was missing: ', path_to_edges_settings ])
%                 error(  'Vectorization must be re-run from the edges step.' )
                            
            end
        case 5
        %%   -  -  -  -  -  -  -  -  -  -  -   Network -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  

            try load( path_to_network_settings )
                                            
            catch
                                
                disp([ 'The previous network settings file was missing: ', path_to_network_settings ])
%                 error(  'Vectorization must be re-run from the network step.' )
                            
            end
    end
end % FOR INPUTS_REQUIRED         
    %%   -  -  -  -  -  -  -  -  -  -  -  -  -  Required Inputs -  -  -  -  -  -  -  -  -  -  -  -  

for workflow_index = inputs_required
    
    switch workflow_index
        
        case 2
        %%   -  -  -  -  -  -  -  -  -  -  -  - Energy -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  

            try load( previous_path_to_energy_settings )
                
                disp([  'Previous energy settings loaded from file at ', previous_path_to_energy_settings ])
                
                if ~ presumptive
                
                    disp(  '[Previous Settings] are provided when prompted for paramater values at the MATLAB command window.' )
            
                end
            catch
                                
                disp([ 'No previous energy settings file was found at ', previous_path_to_energy_settings ])
                
                if ~ presumptive
                
                    disp(   '[Default Settings] are provided when prompted for paramater values at the MATLAB command window.' )
                
                end
                
                % load defaults
                microns_per_voxel                    = [ 1, 1, 1 ]   ;
                radius_of_smallest_vessel_in_microns = 1.5           ;
                radius_of_largest_vessel_in_microns  = 50            ;
                sample_index_of_refraction           = 1.33          ;
                numerical_aperture                   = 0.95          ;
                excitation_wavelength_in_microns     = 1.3           ;
                scales_per_octave                    = 1.5           ;
                max_voxels_per_node_energy           = 1e5           ;
                gaussian_to_ideal_ratio              = 1             ;
                spherical_to_annular_ratio           = 1             ;
                vessel_wall_thickness_in_microns     = 0             ;
                matching_kernel_string               = '3D gaussian conv annular pulse' ;
                approximating_PSF                    = true          ;                        
            
            end
                        
            if presumptive
                
                addParameter( p, 'microns_per_voxel'                   ,                    microns_per_voxel )
                addParameter( p, 'radius_of_smallest_vessel_in_microns', radius_of_smallest_vessel_in_microns )
                addParameter( p, 'radius_of_largest_vessel_in_microns' ,  radius_of_largest_vessel_in_microns )
                addParameter( p, 'sample_index_of_refraction'          ,           sample_index_of_refraction )
                addParameter( p, 'numerical_aperture'                  ,                   numerical_aperture )
                addParameter( p, 'excitation_wavelength_in_microns'    ,     excitation_wavelength_in_microns )
                addParameter( p, 'scales_per_octave'                   ,                    scales_per_octave )
                addParameter( p, 'max_voxels_per_node_energy'          ,           max_voxels_per_node_energy )
                addParameter( p, 'gaussian_to_ideal_ratio'             ,              gaussian_to_ideal_ratio )
                addParameter( p, 'spherical_to_annular_ratio'          ,           spherical_to_annular_ratio )
                addParameter( p, 'approximating_PSF'                   ,                    approximating_PSF )                    

            else

                addParameter( p, 'microns_per_voxel'                   , 'prompt' )
                addParameter( p, 'radius_of_smallest_vessel_in_microns', 'prompt' )
                addParameter( p, 'radius_of_largest_vessel_in_microns' , 'prompt' )
                addParameter( p, 'sample_index_of_refraction'          , 'prompt' )
                addParameter( p, 'numerical_aperture'                  , 'prompt' )
                addParameter( p, 'excitation_wavelength_in_microns'    , 'prompt' )
                addParameter( p, 'scales_per_octave'                   , 'prompt' )
                addParameter( p, 'max_voxels_per_node_energy'          , 'prompt' )
                addParameter( p, 'gaussian_to_ideal_ratio'             , 'prompt' )
                addParameter( p, 'spherical_to_annular_ratio'          , 'prompt' )
                addParameter( p, 'approximating_PSF'                   , 'prompt' )
                

            end

            % always presumptive list:
            addParameter( p, 'matching_kernel_string'          ,           matching_kernel_string ) % !!!! depricated input
            addParameter( p, 'vessel_wall_thickness_in_microns', vessel_wall_thickness_in_microns )                            
            
            parse( p, varargin{ : });
            
            temp_microns_per_voxel                       = p.Results.microns_per_voxel                    ;
            temp_radius_of_smallest_vessel_in_microns    = p.Results.radius_of_smallest_vessel_in_microns ;
            temp_radius_of_largest_vessel_in_microns     = p.Results.radius_of_largest_vessel_in_microns  ;
            temp_sample_index_of_refraction              = p.Results.sample_index_of_refraction           ;
            temp_numerical_aperture                      = p.Results.numerical_aperture                   ;
            temp_excitation_wavelength_in_microns        = p.Results.excitation_wavelength_in_microns     ;
            temp_scales_per_octave                       = p.Results.scales_per_octave                    ;
            temp_max_voxels_per_node_energy              = p.Results.max_voxels_per_node_energy           ;
            temp_gaussian_to_ideal_ratio                 = p.Results.gaussian_to_ideal_ratio              ;
            temp_spherical_to_annular_ratio              = p.Results.spherical_to_annular_ratio           ;
            temp_vessel_wall_thickness_in_microns        = p.Results.vessel_wall_thickness_in_microns     ;
            temp_matching_kernel_string                  = p.Results.matching_kernel_string               ;
            temp_approximating_PSF                       = p.Results.approximating_PSF                    ; 
            %% --------------------------- microns_per_voxel                                        

            if strcmp( temp_microns_per_voxel, 'prompt' )

                temp_microns_per_voxel                                                                                                                  ...
                    = input([ 'Please input the y, x, and z voxel lengths in microns of the input image [[ ', num2str( microns_per_voxel( 1 )), ', ',   ...
                                                                                                              num2str( microns_per_voxel( 2 )), ', ',   ...
                                                                                                              num2str( microns_per_voxel( 3 )), ' ]]: ' ]);    

                if isempty( temp_microns_per_voxel ), temp_microns_per_voxel = microns_per_voxel ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_microns_per_voxel )      ...
               &&  isvector( temp_microns_per_voxel )      ...
               &&    length( temp_microns_per_voxel ) == 3

                if all( temp_microns_per_voxel > 0 )

                    microns_per_voxel = reshape( temp_microns_per_voxel, [ 1, 3 ]);

                else

                    error( 'Parameter microns_per_voxel must contain all positive entries' )

                end                
            else

                error( 'Parameter microns_per_voxel must be a vector with length 3' )

            end % IF valid parameter input
            
            disp([ 'microns_per_voxel                    = [ ', num2str( microns_per_voxel( 1 )), ', ', ...
                                                                num2str( microns_per_voxel( 2 )), ', ', ...
                                                                num2str( microns_per_voxel( 3 )), ' ]'  ])
            %% --------------------------- radius_of_smallest_vessel_in_microns                     

            if strcmp( temp_radius_of_smallest_vessel_in_microns, 'prompt' )

                temp_radius_of_smallest_vessel_in_microns ...
                    = input([ 'Please input the radius of the smallest vessel desired to be detected in microns [ ', num2str( radius_of_smallest_vessel_in_microns ), ' ]: ']);    

                if isempty( temp_radius_of_smallest_vessel_in_microns ), temp_radius_of_smallest_vessel_in_microns = radius_of_smallest_vessel_in_microns ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_radius_of_smallest_vessel_in_microns ) ...
               &&  isscalar( temp_radius_of_smallest_vessel_in_microns )

                if temp_radius_of_smallest_vessel_in_microns > 0 

                    radius_of_smallest_vessel_in_microns = temp_radius_of_smallest_vessel_in_microns ;

                else

                    error( 'Parameter radius_of_smallest_vessel_in_microns must be positive' )

                end                
            else

                error( 'Parameter radius_of_smallest_vessel_in_microns must be a numeric scalar' )

            end % IF valid parameter input

            disp([ 'radius_of_smallest_vessel_in_microns =   ', num2str( radius_of_smallest_vessel_in_microns )])
            %% --------------------------- radius_of_largest_vessel_in_microns                      

            if strcmp( temp_radius_of_largest_vessel_in_microns, 'prompt' )

                temp_radius_of_largest_vessel_in_microns ...
                    = input([ 'Please input the radius of the largest vessel desired to be detected in microns [ ', num2str( radius_of_largest_vessel_in_microns ), ' ]: ']);    

                if isempty( temp_radius_of_largest_vessel_in_microns ), temp_radius_of_largest_vessel_in_microns = radius_of_largest_vessel_in_microns ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_radius_of_largest_vessel_in_microns ) ...
               &&  isscalar( temp_radius_of_largest_vessel_in_microns )

                if temp_radius_of_largest_vessel_in_microns > 0 

                    radius_of_largest_vessel_in_microns = temp_radius_of_largest_vessel_in_microns ;

                else

                    error( 'Parameter radius_of_largest_vessel_in_microns must be positive' )

                end

            else

                error( 'Parameter radius_of_largest_vessel_in_microns must be a numeric scalar' )

            end % IF valid parameter input

            disp([ 'radius_of_largest_vessel_in_microns  =   ', num2str( radius_of_largest_vessel_in_microns )])
            %% --------------------------- approximating_PSF                                        

            if strcmp( temp_approximating_PSF, 'prompt' )

                disp( 'Enter "1" (or "true") if you would you like to approximate the' )
                disp( 'PSF of your two photon microscope according to' )                
                
                temp_approximating_PSF                                                                                                                              ...
                    = logical( input([ '"Nonlinear Magic: Multiphoton Microscopy in the Biosciences" (Zipfel, W.R. et al.) [ ', num2str( approximating_PSF ), ' ]: ']));    

                if isempty( temp_approximating_PSF ), temp_approximating_PSF = approximating_PSF ; end                                                                                                                            
                                
            end
            
            if    islogical( temp_approximating_PSF ) ...
               &&  isscalar( temp_approximating_PSF )

                approximating_PSF = temp_approximating_PSF ;

            else

                error( 'Parameter approximating_PSF must be a scalar logical' )

            end % IF valid parameter input
            
            disp([ 'approximating_PSF                    =   ', num2str( approximating_PSF )])
            
            if approximating_PSF
                %% --------------------------- sample_index_of_refraction                           

                if strcmp( temp_sample_index_of_refraction, 'prompt' )

                    temp_sample_index_of_refraction ...
                        = input([ 'Please input the index of refraction of the sample [ ', num2str( sample_index_of_refraction ), ' ]: ']);    

                    if isempty( temp_sample_index_of_refraction ), temp_sample_index_of_refraction = sample_index_of_refraction ; end                                                                                                                            

                end
                
                if    isnumeric( temp_sample_index_of_refraction ) ...
                   &&  isscalar( temp_sample_index_of_refraction )

                    if temp_sample_index_of_refraction > 0 

                        sample_index_of_refraction = temp_sample_index_of_refraction ;

                    else

                        error( 'Parameter sample_index_of_refraction must be positive' )

                    end
                else

                    error( 'Parameter sample_index_of_refraction must be a numeric scalar' )

                end % IF valid parameter input
                
                disp([ 'sample_index_of_refraction           =   ', num2str( sample_index_of_refraction )])                
                %% --------------------------- numerical_aperture                                   

                if strcmp( temp_numerical_aperture, 'prompt' )

                    temp_numerical_aperture ...
                        = input([ 'Please input the numerical aperture of the microscope objective [ ', num2str( numerical_aperture ), ' ]: ']);    

                    if isempty( temp_numerical_aperture ), temp_numerical_aperture = numerical_aperture ; end                                                                                                                            

                end
                
                if    isnumeric( temp_numerical_aperture ) ...
                   &&  isscalar( temp_numerical_aperture )

                    if temp_numerical_aperture > 0 

                        numerical_aperture = temp_numerical_aperture ;

                    else

                        error( 'Parameter numerical_aperture must be positive' )

                    end
                else

                    error( 'Parameter numerical_aperture must be a numeric scalar' )

                end % IF valid parameter input
            
                disp([ 'numerical_aperture                   =   ', num2str( numerical_aperture )])
                %% --------------------------- excitation_wavelength_in_microns                                

                if strcmp( temp_excitation_wavelength_in_microns, 'prompt' )

                    temp_excitation_wavelength_in_microns ...
                        = input([ 'Please input the wavelength of the LASER in microns (proportional to the size of the PSF) [ ', num2str( excitation_wavelength_in_microns ), ' ]: ']);    

                    if isempty( temp_excitation_wavelength_in_microns ), temp_excitation_wavelength_in_microns = excitation_wavelength_in_microns ; end                                                                                                                            

                end
                
                if    isnumeric( temp_excitation_wavelength_in_microns ) ...
                   &&  isscalar( temp_excitation_wavelength_in_microns )

                    if temp_excitation_wavelength_in_microns > 0 

                        excitation_wavelength_in_microns = temp_excitation_wavelength_in_microns ;

                    else

                        error( 'Parameter excitation_wavelength_in_microns must be positive' )

                    end
                else

                    error( 'Parameter excitation_wavelength_in_microns must be a numeric scalar' )

                end % IF valid parameter input

                disp([ 'excitation_wavelength_in_microns                =   ', num2str( excitation_wavelength_in_microns )])

            else % not approximating the PSF
                
                sample_index_of_refraction       = [ ];
                numerical_aperture               = [ ];
                excitation_wavelength_in_microns = [ ];
                
            end % IF approximating_PSF
            %% --------------------------- scales_per_octave                                        
            
            if strcmp( temp_scales_per_octave, 'prompt' )

                temp_scales_per_octave                                                                                                                               ...
                    = input([ 'Please input the desired number of vessel sizes to detected per doubling of the radius cubed [ ', num2str( scales_per_octave ), ' ]: ']);    

                if isempty( temp_scales_per_octave ), temp_scales_per_octave = scales_per_octave ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_scales_per_octave ) ...
               &&  isscalar( temp_scales_per_octave )

                if temp_scales_per_octave > 0 

                    scales_per_octave = temp_scales_per_octave ;

                else

                    error( 'Parameter scales_per_octave must be positive' )

                end
            else

                error( 'Parameter scales_per_octave must be a numeric scalar' )

            end % IF valid parameter input
        
            disp([ 'scales_per_octave                    =   ', num2str( scales_per_octave )])
            %% --------------------------- max_voxels_per_node_energy                               
            
            if strcmp( temp_max_voxels_per_node_energy, 'prompt' )

                temp_max_voxels_per_node_energy                                                                                                                           ...
                    = input([ 'Please input the desired number of voxels per computational node in parallel computation [ ', num2str( max_voxels_per_node_energy ), ' ]: ']);    

                if isempty( temp_max_voxels_per_node_energy ), temp_max_voxels_per_node_energy = max_voxels_per_node_energy ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_max_voxels_per_node_energy ) ...
               &&  isscalar( temp_max_voxels_per_node_energy )

                if temp_max_voxels_per_node_energy > 0 

                    max_voxels_per_node_energy = temp_max_voxels_per_node_energy ;

                else

                    error( 'Parameter max_voxels_per_node_energy must be positive' )

                end
            else

                error( 'Parameter max_voxels_per_node_energy must be a numeric scalar' )

            end % IF valid parameter input
            
            disp([ 'max_voxels_per_node_energy           =   ', num2str( max_voxels_per_node_energy )])
            %% --------------------------- gaussian_to_ideal_ratio                                  
            
            if strcmp( temp_gaussian_to_ideal_ratio, 'prompt' )

                temp_gaussian_to_ideal_ratio                                                                                                ...
                    = input([ 'Please input the desired amount of Gaussian pulse (vs spherical or annular) [ ', num2str( gaussian_to_ideal_ratio ), ' ]: ']);    

                if isempty( temp_gaussian_to_ideal_ratio ), temp_gaussian_to_ideal_ratio = gaussian_to_ideal_ratio ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_gaussian_to_ideal_ratio ) ...
               &&  isscalar( temp_gaussian_to_ideal_ratio )

                if temp_gaussian_to_ideal_ratio >= 0 && temp_gaussian_to_ideal_ratio <= 1

                    gaussian_to_ideal_ratio = temp_gaussian_to_ideal_ratio ;

                else

                    error( 'Parameter gaussian_to_ideal_ratio must be between 0 and 1 inclusive' )

                end                
            else

                error( 'Parameter gaussian_to_ideal_ratio must be a numeric scalar' )

            end % IF valid parameter input

            disp([ 'gaussian_to_ideal_ratio                =   ', num2str( gaussian_to_ideal_ratio )]) 
            %% --------------------------- spherical_to_annular_ratio                               
            
            if strcmp( temp_spherical_to_annular_ratio, 'prompt' )

                temp_spherical_to_annular_ratio                                                                                                ...
                    = input([ 'Please input the desired amount of Gaussian pulse (vs spherical or annular) [ ', num2str( spherical_to_annular_ratio ), ' ]: ']);    

                if isempty( temp_spherical_to_annular_ratio ), temp_spherical_to_annular_ratio = spherical_to_annular_ratio ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_spherical_to_annular_ratio ) ...
               &&  isscalar( temp_spherical_to_annular_ratio )

                if temp_spherical_to_annular_ratio >= 0 && temp_spherical_to_annular_ratio <= 1

                    spherical_to_annular_ratio = temp_spherical_to_annular_ratio ;

                else

                    error( 'Parameter spherical_to_annular_ratio must be between 0 and 1 inclusive' )

                end                
            else

                error( 'Parameter spherical_to_annular_ratio must be a numeric scalar' )

            end % IF valid parameter input

            disp([ 'spherical_to_annular_ratio                =   ', num2str( spherical_to_annular_ratio )])
            %% --------------------------- matching_kernel_string                                   
            
            matching_kernel_string_values = {    '3D gaussian',                      ...
                                                 '3D gaussian conv spherical pulse', ...
                                                 '3D gaussian conv annular pulse', ...
                                                  'spherical pulse',                 ...
                                                    'annular pulse',                 ...
                                                    'annular pulse V2',              ...
                                                 '3D annular gaussian',              ...
                                                 '3D annular gaussian V2',           ...
                                                     'radial gaussian'               };

            matching_kernels_of_vessel_wall = matching_kernel_string_values([ 3, 5 : end ]);
                        
            if strcmp( temp_matching_kernel_string, 'prompt' )

                disp( 'Matching kernel values:' )
                
                number_of_matching_kernel_string_values = length( matching_kernel_string_values );                

                for matching_kernel_index = 1 : number_of_matching_kernel_string_values

                    disp([ '          ', matching_kernel_string_values{ matching_kernel_index }])

                end

                temp_matching_kernel_string                                                                          ...
                    = input([ 'Please input the matching kernel to be used [ ', matching_kernel_string, ' ]: '], 's' );    

                if isempty( temp_matching_kernel_string ), temp_matching_kernel_string = matching_kernel_string ; end                                                                                                                            
                                
            end
            
            if      ischar( temp_matching_kernel_string ) ...
               || isstring( temp_matching_kernel_string )

                matching_kernel_string = validatestring( temp_matching_kernel_string, matching_kernel_string_values );

            else

                error( 'Parameter matching_kernel_string must be type char or string' )

            end % IF valid parameter input

            disp([ 'matching_kernel_string               =   ', matching_kernel_string ])
            
            switch matching_kernel_string
                
                case matching_kernels_of_vessel_wall
                    %% --------------------------- vessel_wall_thickness_in_microns                         

                    if strcmp( temp_vessel_wall_thickness_in_microns, 'prompt' )

                        temp_vessel_wall_thickness_in_microns ...
                            = input([ 'Please input the thickness of the vessel wall to be detected in microns [ ', num2str( vessel_wall_thickness_in_microns ), ' ]: ']);    

                        if isempty( temp_vessel_wall_thickness_in_microns ), temp_vessel_wall_thickness_in_microns = vessel_wall_thickness_in_microns ; end                                                                                                                            

                    end
                    
                    if    isnumeric( temp_vessel_wall_thickness_in_microns ) ...
                       &&  isscalar( temp_vessel_wall_thickness_in_microns )

                        if vessel_wall_thickness_in_microns >= 0 

                            vessel_wall_thickness_in_microns = temp_vessel_wall_thickness_in_microns ;

                        else

                            error( 'Parameter vessel_wall_thickness_in_microns must be nonnegative' )

                        end                
                    else

                        error( 'Parameter vessel_wall_thickness_in_microns must be a numeric scalar' )

                    end % IF valid parameter input

                    disp([ 'vessel_wall_thickness_in_microns  =   ', num2str( vessel_wall_thickness_in_microns )])
                
               otherwise % no vessel wall thickness required
                   
%                    vessel_wall_thickness_in_microns = 1 ;
            
            end
            %% ---------------------- dependent variables                                           
            
            if approximating_PSF
                
                % from Figure 4 of Nonlinear magic: multiphoton microscopy in the biosciences by Warren R
                % Zipfel, Rebecca M Williams & Watt W Webb, with sigma of the intensity gaussian = omega for the
                % squared intensity gaussian.  SAM 4/11/18            
                if numerical_aperture <= 0.7
                    
                    coefficient = 0.320 ; exponent = 1    ;

                else
                    
                    coefficient = 0.325 ; exponent = 0.91 ;
                                            
                end
                
                microns_per_sigma_microscope_PSF                                                     ...
                    = excitation_wavelength_in_microns / 2 ^ 0.5                                                ...
                    * [ coefficient / numerical_aperture ^ exponent * [ 1, 1 ],                      ...
                        .532 / (     sample_index_of_refraction                                      ...
                                 - ( sample_index_of_refraction ^ 2 - numerical_aperture ^ 2 ) ^ 0.5 )];
            else % ideal PSF assumed
               
                microns_per_sigma_microscope_PSF = [ 0, 0, 0 ];
                
            end
                        
            pixels_per_sigma_PSF = microns_per_sigma_microscope_PSF ./ microns_per_voxel ;                          

            largest_per_smallest_scale_volume_ratio = (   radius_of_largest_vessel_in_microns        ...
                                                        / radius_of_smallest_vessel_in_microns ) ^ 3 ;

            final_scale = round( log( largest_per_smallest_scale_volume_ratio ) / log( 2 ) * scales_per_octave );

            % need to pad by 1 scale on the top and bottom because we need at least one above or below to
            % call it an extreme point later when finding vertices (local energy minima in physical and
            % scale spaces (4D locale)).
            scale_ordinates = ( - 1 : final_scale + 1 )' ;

            scale_factor_range = 2 .^ ( scale_ordinates / scales_per_octave / 3 );

            lumen_radius_in_microns_range = radius_of_smallest_vessel_in_microns * scale_factor_range ;

            lumen_radius_in_pixels_range = lumen_radius_in_microns_range ./ microns_per_voxel ;
            
            disp([ 'microns_per_sigma_microscope_PSF     = [ ', num2str( microns_per_sigma_microscope_PSF( 1 )), ', ', ...
                                                                num2str( microns_per_sigma_microscope_PSF( 2 )), ', ', ...
                                                                num2str( microns_per_sigma_microscope_PSF( 3 )), ' ]'  ])
            
            phrase = num2str( lumen_radius_in_microns_range( 1 ));
            
            for lumen_radius = lumen_radius_in_microns_range( 2 : end )'
                
                phrase = [ phrase, ', ', num2str( lumen_radius )];
                
            end
                                                            
            disp([ 'lumen_radius_in_microns_range        = [ ', phrase, ' ]' ])                                                            
            %% ---------------------- saving                                                        

            disp([ 'Saving energy workflow settings file: ', path_to_energy_settings ])

            save(  path_to_energy_settings              , ...
                  'lumen_radius_in_microns_range'       , ...                
                  'lumen_radius_in_pixels_range'        , ...
                  'pixels_per_sigma_PSF'                , ...
                  'microns_per_voxel'                   , ...
                  'radius_of_largest_vessel_in_microns' , ...
                  'radius_of_smallest_vessel_in_microns', ...
                  'approximating_PSF'                   , ...                  
                  'sample_index_of_refraction'          , ...
                  'numerical_aperture'                  , ...
                  'excitation_wavelength_in_microns'    , ...
                  'scales_per_octave'                   , ...
                  'max_voxels_per_node_energy'          , ...                  
                  'gaussian_to_ideal_ratio'             , ...
                  'spherical_to_annular_ratio'          , ...
                  'matching_kernel_string'              , ...
                  'vessel_wall_thickness_in_microns'      )
        case 3
        %%  -  -  -  -  -  -  -  -  -  -  -  - Vertices  -  -  -  -  -  -  -  -  -  -  -  -  -  -   

            try load( previous_path_to_vertices_settings )
                
                disp([  'Previous vertices settings loaded from file at ', previous_path_to_vertices_settings ])
                
                if ~ presumptive
                
                    disp(  '[Previous Settings] are provided when prompted for paramater values at the MATLAB command window.' )
            
                end
            catch
                                
                disp([ 'No previous vertices settings file was found at ', previous_path_to_vertices_settings ])
                
                if ~ presumptive
                
                    disp(   '[Default Settings] are provided when prompted for paramater values at the MATLAB command window.' )

                end
                
                % load defaults                
                space_strel_apothem       =     1 ;
                energy_upper_bound        =     0 ;
                max_voxels_per_node       =  6000 ;
                length_dilation_ratio     =     1 ;
                exporting_vertex_curation = false ;
                
            end
            
%             if presumptive
% 
%                 addParameter( p, 'energy_upper_bound' , energy_upper_bound )
%                 
%             else
%                 
%                 addParameter( p, 'energy_upper_bound' , 'prompt' )
%                 
%             end

            % always presumptive list:
            addParameter( p, 'energy_upper_bound'       ,        energy_upper_bound )
            addParameter( p, 'space_strel_apothem'      ,       space_strel_apothem )
            addParameter( p, 'length_dilation_ratio'    ,     length_dilation_ratio )
            addParameter( p, 'exporting_vertex_curation', exporting_vertex_curation )     
            addParameter( p, 'max_voxels_per_node'      ,       max_voxels_per_node )            
            
            parse( p, varargin{ : });
            
            temp_space_strel_apothem        = p.Results.space_strel_apothem       ;
            temp_energy_upper_bound         = p.Results.energy_upper_bound        ;
            temp_max_voxels_per_node        = p.Results.max_voxels_per_node       ;
            temp_length_dilation_ratio      = p.Results.length_dilation_ratio     ;
            temp_exporting_vertex_curation  = p.Results.exporting_vertex_curation ;
            %% --------------------------- space_strel_apothem                                      

            if strcmp( temp_space_strel_apothem, 'prompt' )

                temp_space_strel_apothem                                                                                                                          ...
                    = input([ 'Please input the desired minimum number of voxels between vertex detections [ ', num2str( space_strel_apothem * 2 ), ' ]: ' ]) / 2 ;    

                if isempty( temp_space_strel_apothem ), temp_space_strel_apothem = space_strel_apothem ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_space_strel_apothem ) ...
               &&  isscalar( temp_space_strel_apothem )

                if temp_space_strel_apothem > 0

                    if ~ logical( temp_space_strel_apothem - round( temp_space_strel_apothem ))

                        space_strel_apothem = temp_space_strel_apothem ;

                    else 

                        error( 'Parameter space_strel_apothem must be an integer' )

                    end
                else

                    error( 'Parameter space_strel_apothem must be positive' )

                end
            else

                error( 'Parameter space_strel_apothem must be a numeric scalar' )

            end % IF valid parameter input
            
            disp([ 'space_strel_apothem       =   ', num2str( space_strel_apothem )])
            %% --------------------------- energy_upper_bound                                       

            if strcmp( temp_energy_upper_bound, 'prompt' )

                temp_energy_upper_bound                                                                                                    ...
                    = input([ 'Please input the desired maximum energy of the vertex detections [ ', num2str( energy_upper_bound ), ' ]: ' ]);    

                if isempty( temp_energy_upper_bound ), temp_energy_upper_bound = energy_upper_bound ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_energy_upper_bound ) ...
               &&  isscalar( temp_energy_upper_bound )

                if temp_energy_upper_bound <= 0

                    energy_upper_bound = temp_energy_upper_bound ;

                else

                    error( 'Parameter energy_upper_bound must be non-positive' )

                end
            else

                error( 'Parameter energy_upper_bound must be a numeric scalar' )

            end % IF valid parameter input
            
            disp([ 'energy_upper_bound        =   ', num2str( energy_upper_bound )])
            %% --------------------------- max_voxels_per_node                                      

            if strcmp( temp_max_voxels_per_node, 'prompt' )

                temp_max_voxels_per_node                                                                                                                            ...
                    = input([ 'Please input the desired number of voxels per computational node in parallel computation [ ', num2str( max_voxels_per_node ), ' ]: ' ]);    

                if isempty( temp_max_voxels_per_node ), temp_max_voxels_per_node = max_voxels_per_node ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_max_voxels_per_node ) ...
               &&  isscalar( temp_max_voxels_per_node )

                if temp_max_voxels_per_node > 0

                    max_voxels_per_node = temp_max_voxels_per_node ;

                else

                    error( 'Parameter max_voxels_per_node must be positive' )

                end
            else

                error( 'Parameter max_voxels_per_node must be a numeric scalar' )

            end % IF valid parameter input
            
            disp([ 'max_voxels_per_node    =   ', num2str( max_voxels_per_node )])
            %% --------------------------- length_dilation_ratio                                    

            if strcmp( temp_length_dilation_ratio, 'prompt' )

                temp_length_dilation_ratio                                                                                                                                               ...
                    = input([ 'Please input the ratio of rendering length to detection length for the vertex volume exclusion segmentation [ ', num2str( length_dilation_ratio ), ' ]: ' ]);    

                if isempty( temp_length_dilation_ratio ), temp_length_dilation_ratio = length_dilation_ratio ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_length_dilation_ratio ) ...
               &&  isscalar( temp_length_dilation_ratio )

                if temp_length_dilation_ratio > 0

                    length_dilation_ratio = temp_length_dilation_ratio ;

                else

                    error( 'Parameter length_dilation_ratio must be positive' )

                end
            else

                error( 'Parameter length_dilation_ratio must be a numeric scalar' )

            end % IF valid parameter input
            
            disp([ 'length_dilation_ratio       =   ', num2str( length_dilation_ratio )])
            %% --------------------------- exporting_vertex_curation                                

            if strcmp( temp_exporting_vertex_curation, 'prompt' )

                temp_exporting_vertex_curation                                                                                                          ...
                    = logical( input([ 'Please input 1 or true to save a curation file for exporting [ ', num2str( exporting_vertex_curation ), ' ]: ' ]));    

                if isempty( temp_exporting_vertex_curation ), temp_exporting_vertex_curation = exporting_vertex_curation ; end                                                                                                                            
                                
            end
            
            if    islogical( temp_exporting_vertex_curation ) ...
               &&  isscalar( temp_exporting_vertex_curation )

                exporting_vertex_curation = temp_exporting_vertex_curation ;

            else

                error( 'Parameter exporting_vertex_curation must be a scalar logical' )

            end % IF valid parameter input
            
            disp([ 'exporting_vertex_curation =   ', num2str( exporting_vertex_curation )])
            %% ---------------------- saving                                                        

            disp([ 'Saving vertices workflow settings file: ', path_to_vertices_settings ])        

            save(  path_to_vertices_settings, ...
                  'space_strel_apothem'     , ...
                  'energy_upper_bound'      , ...
                  'max_voxels_per_node'     , ...
                  'length_dilation_ratio'   , ...
                 'exporting_vertex_curation'  );
        case 4
        %% -  -  -  -  -  -  -  -  -  -  -  -  - Edges  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 

            try load( previous_path_to_edges_settings )
                
                disp([  'Previous edges settings loaded from file at ', previous_path_to_edges_settings ])
                
                if ~ presumptive
                    
                    disp(  '[Previous Settings] are provided when prompted for paramater values at the MATLAB command window.' )
            
                end
            catch
                                
                disp([ 'No previous edges settings file was found at ', previous_path_to_edges_settings ])
                
                if ~ presumptive
                
                    disp(   '[Default Settings] are provided when prompted for paramater values at the MATLAB command window.' )
                
                end
                
                % load defaults
                max_edge_length_per_origin_radius = 60  ;
                space_strel_apothem_edges         = 1   ;
                number_of_edges_per_vertex        = 4   ;
                sigma_edge_smoothing              = 0   ;
                length_dilation_ratio_vertices    = 2   ;
                length_dilation_ratio_edges       = 2/3 ;            

            end
            
            if presumptive
                
                addParameter( p, 'max_edge_length_per_origin_radius', max_edge_length_per_origin_radius )
                addParameter( p, 'number_of_edges_per_vertex'       ,        number_of_edges_per_vertex )
%                 addParameter( p, 'sigma_edge_smoothing'             ,              sigma_edge_smoothing )
                
            else
            
                addParameter( p, 'max_edge_length_per_origin_radius', 'prompt' )
                addParameter( p, 'number_of_edges_per_vertex'       , 'prompt' )
%                 addParameter( p, 'sigma_edge_smoothing'             , 'prompt' )                
                
            end
            
            % always presumptive list:
            addParameter( p, 'space_strel_apothem_edges'     ,      space_strel_apothem_edges )
            addParameter( p, 'sigma_edge_smoothing'          ,           sigma_edge_smoothing )                
            addParameter( p, 'length_dilation_ratio_vertices', length_dilation_ratio_vertices )
            addParameter( p, 'length_dilation_ratio_edges'   ,    length_dilation_ratio_edges )

            parse( p, varargin{ : });

            temp_max_edge_length_per_origin_radius = p.Results.max_edge_length_per_origin_radius ;
            temp_space_strel_apothem_edges         = p.Results.space_strel_apothem_edges         ;
            temp_number_of_edges_per_vertex        = p.Results.number_of_edges_per_vertex        ;
            temp_sigma_edge_smoothing              = p.Results.sigma_edge_smoothing              ;
            temp_length_dilation_ratio_vertices    = p.Results.length_dilation_ratio_vertices    ;
            temp_length_dilation_ratio_edges       = p.Results.length_dilation_ratio_edges       ;
            %% --------------------------- max_edge_length_per_origin_radius                        

            if strcmp( temp_max_edge_length_per_origin_radius, 'prompt' )

                temp_max_edge_length_per_origin_radius                                                                                                                                   ...
                    = input([ 'Please input the ratio of the desired maximum length of each edge to its original vertex radius [ ', num2str( max_edge_length_per_origin_radius ), ' ]: ' ]);    

                if isempty( temp_max_edge_length_per_origin_radius ), temp_max_edge_length_per_origin_radius = max_edge_length_per_origin_radius ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_max_edge_length_per_origin_radius ) ...
               &&  isscalar( temp_max_edge_length_per_origin_radius )

                if temp_max_edge_length_per_origin_radius > 0

                    max_edge_length_per_origin_radius = temp_max_edge_length_per_origin_radius ;

                else

                    error( 'Parameter max_edge_length_per_origin_radius must be positive' )

                end
            else

                error( 'Parameter max_edge_length_per_origin_radius must be a numeric scalar' )

            end % IF valid parameter input
            
            disp([ 'max_edge_length_per_origin_radius    =   ', num2str( max_edge_length_per_origin_radius )])
            %% --------------------------- space_strel_apothem_edges                                

            if strcmp( temp_space_strel_apothem_edges, 'prompt' )

                temp_space_strel_apothem_edges                                                                                                                                      ...
                    = input([ 'Please input the desired apothem of the structuring element for the energy neighborhood in units of voxels [ ', num2str( space_strel_apothem_edges ), ' ]: ' ]);    

                if isempty( temp_space_strel_apothem_edges ), temp_space_strel_apothem_edges = space_strel_apothem_edges ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_space_strel_apothem_edges ) ...
               &&  isscalar( temp_space_strel_apothem_edges )

                if temp_space_strel_apothem_edges > 0

                    if ~ logical( temp_space_strel_apothem_edges - round( temp_space_strel_apothem_edges ))

                        space_strel_apothem_edges = temp_space_strel_apothem_edges ;

                    else 

                        error( 'Parameter space_strel_apothem_edges must be an integer' )

                    end
                else

                    error( 'Parameter space_strel_apothem_edges must be positive' )

                end
            else

                error( 'Parameter space_strel_apothem_edges must be a numeric scalar' )

            end % IF valid parameter input
            
            disp([ 'space_strel_apothem_edges            =   ', num2str( space_strel_apothem_edges )])
            %% --------------------------- number_of_edges_per_vertex                               

            if strcmp( temp_number_of_edges_per_vertex, 'prompt' )

                temp_number_of_edges_per_vertex                                                                                        ...
                    = input([ 'Please input the max number of edges per seed vertex [ ', num2str( number_of_edges_per_vertex ), ' ]: ' ]);    

                if isempty( temp_number_of_edges_per_vertex ), temp_number_of_edges_per_vertex = number_of_edges_per_vertex ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_number_of_edges_per_vertex ) ...
               &&  isscalar( temp_number_of_edges_per_vertex )

                if temp_number_of_edges_per_vertex > 0

                    if ~ logical( temp_number_of_edges_per_vertex - round( temp_number_of_edges_per_vertex ))

                        number_of_edges_per_vertex = temp_number_of_edges_per_vertex ;

                    else 

                        error( 'Parameter number_of_edges_per_vertex must be an integer' )

                    end
                else

                    error( 'Parameter number_of_edges_per_vertex must be positive' )

                end
            else

                error( 'Parameter number_of_edges_per_vertex must be a numeric scalar' )

            end % IF valid parameter input
            
            disp([ 'number_of_edges_per_vertex           =   ', num2str( number_of_edges_per_vertex )])                        
            %% --------------------------- sigma_edge_smoothing                                    

            if strcmp( temp_sigma_edge_smoothing, 'prompt' )

                temp_sigma_edge_smoothing                                                                                                                                               ...
                    = input([ 'Please input the standard deviation of the Gaussian to be used to smooth along the edges as a factor of the vessel radius [ ', num2str( sigma_edge_smoothing ), ' ]: ' ]);    

                if isempty( temp_sigma_edge_smoothing ), temp_sigma_edge_smoothing = sigma_edge_smoothing ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_sigma_edge_smoothing ) ...
               &&  isscalar( temp_sigma_edge_smoothing )

                if temp_sigma_edge_smoothing >= 0

                    sigma_edge_smoothing = temp_sigma_edge_smoothing ;

                else

                    error( 'Parameter sigma_edge_smoothing must be nonnegative' )

                end
            else

                error( 'Parameter sigma_edge_smoothing must be a numeric scalar' )

            end % IF valid parameter input
            
            disp([ 'sigma_edge_smoothing                 =   ', num2str( sigma_edge_smoothing )])
            %% --------------------------- length_dilation_ratio_edges                              

            if strcmp( temp_length_dilation_ratio_edges, 'prompt' )

                temp_length_dilation_ratio_edges                                                                                                                                                  ...
                    = input([ 'Please input the ratio of edge rendering length to detection length for the edge volume exclusion segmentation [ ', num2str( length_dilation_ratio_edges ), ' ]: ' ]);    

                if isempty( temp_length_dilation_ratio_edges ), temp_length_dilation_ratio_edges = length_dilation_ratio_edges ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_length_dilation_ratio_edges ) ...
               &&  isscalar( temp_length_dilation_ratio_edges )

                if temp_length_dilation_ratio_edges > 0

                    length_dilation_ratio_edges = temp_length_dilation_ratio_edges ;

                else

                    error( 'Parameter length_dilation_ratio_edges must be positive' )

                end
            else

                error( 'Parameter length_dilation_ratio_edges must be a numeric scalar' )

            end % IF valid parameter input
            
            disp([ 'length_dilation_ratio_edges          =   ', num2str( length_dilation_ratio_edges )])
            %% --------------------------- length_dilation_ratio_vertices                           

            if strcmp( temp_length_dilation_ratio_vertices, 'prompt' )

                temp_length_dilation_ratio_vertices                                                                                                                                                    ...
                    = input([ 'Please input the ratio of vertex rendering length to detection length for the edge volume exclusion segmentation [ ', num2str( length_dilation_ratio_vertices ), ' ]: ' ]);    

                if isempty( temp_length_dilation_ratio_vertices ), temp_length_dilation_ratio_vertices = length_dilation_ratio_vertices ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_length_dilation_ratio_vertices ) ...
               &&  isscalar( temp_length_dilation_ratio_vertices )

                if temp_length_dilation_ratio_vertices > 0
                    
%                     if temp_length_dilation_ratio_vertices <= length_dilation_ratio

                        length_dilation_ratio_vertices = temp_length_dilation_ratio_vertices ;

%                     else
% 
%                         error( 'Parameter length_dilation_ratio_vertices must be less than or equal to the parameter, length_dilation_ratio' )
% 
%                     end
                else

                    error( 'Parameter length_dilation_ratio_vertices must be positive' )

                end
            else

                error( 'Parameter length_dilation_ratio_vertices must be a numeric scalar' )

            end % IF valid parameter input
            
            disp([ 'length_dilation_ratio_vertices       =   ', num2str( length_dilation_ratio_vertices )])            
            %% ---------------------- saving                                                        

            disp([ 'Saving edges workflow settings file: ', path_to_edges_settings ])        

            save(  path_to_edges_settings,             ...
                  'max_edge_length_per_origin_radius', ...
                  'space_strel_apothem_edges',         ...
                  'number_of_edges_per_vertex',        ...
                  'sigma_edge_smoothing',             ...
                  'length_dilation_ratio_edges',       ...
                  'length_dilation_ratio_vertices'     );
        case 5
        %%   -  -  -  -  -  -  -  -  -  -  -   Network -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  

            try load( previous_path_to_network_settings )
                                
                disp([  'Previous network settings loaded from file at ', previous_path_to_network_settings ])
                
                if ~ presumptive
                
                    disp(  '[Previous Settings] are provided when prompted for paramater values at the MATLAB command window.' )
    
                end
            catch
                                
                disp([ 'No previous network settings file was found at ', previous_path_to_network_settings ])
                
                if ~ presumptive
                
                    disp(   '[Default Settings] are provided when prompted for paramater values at the MATLAB command window.' )
    
                end
                
                % load defaults
                sigma_strand_smoothing = 1 ;                
                
            end
                        
            if presumptive
                
                addParameter( p, 'sigma_strand_smoothing', sigma_strand_smoothing )
            
            else
                
                addParameter( p, 'sigma_strand_smoothing', 'prompt' )
                
            end
            
            parse( p, varargin{ : });
            
            temp_sigma_strand_smoothing = p.Results.sigma_strand_smoothing ;
            %% --------------------------- sigma_strand_smoothing                                   

            if strcmp( temp_sigma_strand_smoothing, 'prompt' )

                temp_sigma_strand_smoothing                                                                                                                                               ...
                    = input([ 'Please input the standard deviation of the Gaussian to be used to smooth along the strands as a factor of the vessel radius [ ', num2str( sigma_strand_smoothing ), ' ]: ' ]);    

                if isempty( temp_sigma_strand_smoothing ), temp_sigma_strand_smoothing = sigma_strand_smoothing ; end                                                                                                                            
                                
            end
            
            if    isnumeric( temp_sigma_strand_smoothing ) ...
               &&  isscalar( temp_sigma_strand_smoothing )

                if temp_sigma_strand_smoothing > 0

                    sigma_strand_smoothing = temp_sigma_strand_smoothing ;

                else

                    error( 'Parameter sigma_strand_smoothing must be nonnegative' )

                end
            else

                error( 'Parameter sigma_strand_smoothing must be a numeric scalar' )

            end % IF valid parameter input
            
            disp([ 'sigma_strand_smoothing    =   ', num2str( sigma_strand_smoothing )])
            %% ---------------------- saving                                                        

            disp([ 'Saving network workflow settings file: ', path_to_network_settings ])        

            save(  path_to_network_settings, ...
                  'sigma_strand_smoothing'   );            

    end
end % FOR INPUTS_REQUIRED
    %%  -  -  -  -  -  -  -  -  -  -  -  -  -  -   Other Inputs  -  -  -  -  -  -  -  -  -  -  -  - 

for workflow_index = inputs_not_required
    switch workflow_index
        case 2
        %%   -  -  -  -  -  -  -  -  -  -  -  - Energy -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
        case 3
        %%  -  -  -  -  -  -  -  -  -  -  -  - Vertices  -  -  -  -  -  -  -  -  -  -  -  -  -  -   
        case 4
        %% -  -  -  -  -  -  -  -  -  -  -  -  - Edges  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - 
        case 5
        %%   -  -  -  -  -  -  -  -  -  -  -   Network -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  
    end
end % FOR INPUTS_NOT_REQUIRED

%% 1: Original Image
%% --------------------------------------- Visualization ------------------------------------------ 

if visual( 1 )
    
    for ROI_index = ROI_index_range
                
        path_to_original_visual = [ visual_data_directory, original_data_handle, ROI_names{ ROI_index }, '.tif' ]; % TIF file path       
        
        disp([ 'Writing .tif of the h5 original data to: ', path_to_original_visual ])        
        
        path_to_original_data = [ data_directory, original_data_handle, ROI_names{ ROI_index }]; % h5  file path
        
        mat2tif( h52mat( path_to_original_data ), path_to_original_visual );
        
    end % FOR ROI
end % IF visual

%% 2: Calculate the multi-scale energy field  from the original images
%% ----------------------------------------- Production ------------------------------------------- 

if productive( 2 )
        
    disp( 'Running Energy workflow...' )        
    
    % Suppress warnings saying that temporary variables will be removed on each iteration in PARFOR
    warning( 'off', 'MATLAB:mir_warning_maybe_uninitialized_temporary' )
    
    for ROI_index = ROI_index_range
        tic
                
        get_energy_V202(     matching_kernel_string, lumen_radius_in_microns_range, ...
                               vessel_wall_thickness_in_microns, microns_per_voxel, ...
                  pixels_per_sigma_PSF, max_voxels_per_node_energy, data_directory, ...
                                   [ original_data_handle, ROI_names{ ROI_index }], ...
                                   [        energy_handle, ROI_names{ ROI_index }], ...
                                                      gaussian_to_ideal_ratio, spherical_to_annular_ratio )
        
%         runtime_measurement_string = evalc( 'toc' );
%          
%         digit_indices = regexp( runtime_measurement_string, '\d' );        
%         
%         energy_runtime_in_seconds = sscanf( runtime_measurement_string( digit_indices( 1 ) : end ), '%f' );

        energy_runtime_in_seconds = toc ;
        
        disp([ 'Runtime for energy workflow for image ', ROI_names{ ROI_index }( 2 : end ), ' was ', num2str( round( energy_runtime_in_seconds )), ' seconds' ])

        
        path_to_energy_data = [ data_directory, energy_handle, ROI_names{ ROI_index }];
        
%         energy_file_info = h5info( path_to_energy_data );
% 
%         size_of_image = energy_file_info.Datasets.Dataspace.Size ;  
%         
%         size_of_image = size_of_image( 1 : 3 );
        
        path_to_original_data = [ data_directory, original_data_handle, ROI_names{ ROI_index }]; % h5  file path        
        
        original_image = h52mat( path_to_original_data );
        
        size_of_image = size( original_image );
        
        intensity_limits = quantile( double( original_image( : )), [ 0.01, 0.99 ]);
        
        clear original_image
        
        save( path_to_energy_data, ...
                                          'energy_runtime_in_seconds', ...
                                          'size_of_image',             ...
                                          'intensity_limits'           )
        
    end % FOR ROI
    
    % Restore warnings saying that temporary variables will be removed on each iteration in PARFOR
    warning( 'on', 'MATLAB:mir_warning_maybe_uninitialized_temporary' )    
        
    production_times{ 2 } = attempted_production_times{ 2 };      
    
    save( path_to_workflow_settings, ...
          'production_times'       , ...
          '-append'                  );
          
end % IF productive
%% --------------------------------------- Visualization ------------------------------------------ 

if visual( 2 )
    
    for ROI_index = ROI_index_range
        
        path_to_energy_data = [        data_directory, energy_handle, ROI_names{ ROI_index }];               
        energy_visual_file  = [ visual_data_directory, energy_handle, ROI_names{ ROI_index }];
                
        load([ path_to_energy_data, '.mat' ])
        
        mat2tif( h52mat( path_to_energy_data,           ...
                        [ 1, 1, 1, 1 ],                 ...
                        [ size_of_image( 1 : 3 ), 1 ]), ...
                   [ energy_visual_file, '_indices.tif' ])
             
        mat2tif( h52mat( path_to_energy_data,           ...
                        [ 1, 1, 1, 2 ],                 ...
                        [ size_of_image( 1 : 3 ), 1 ]), ...
                   [ energy_visual_file,         '.tif' ])
             
    end % FOR ROI
end % IF visual

%% 3: get the vertices from the local minima of the 4D energy field
%% ----------------------------------------- Production ------------------------------------------- 
if productive( 3 )
        
    disp( 'Running Vertices workflow...' )
    
    % Suppress warnings saying that temporary variables will be removed on each iteration in PARFOR
    warning( 'off', 'MATLAB:mir_warning_maybe_uninitialized_temporary' )

    for ROI_index = ROI_index_range
        tic
        
        path_to_energy_data        = [     data_directory,        energy_handle, ROI_names{ ROI_index }]; % mat file path        
        path_to_original_data      = [     data_directory, original_data_handle, ROI_names{ ROI_index }]; % h5  file path
        path_to_saved_curation     = [ curation_directory,      vertices_handle, ROI_names{ ROI_index }]; % logicals path
        
        load([ path_to_energy_data, '.mat' ])

        [ vertex_space_subscripts, vertex_scale_subscripts, vertex_energies ]                       ...
            = get_vertices_V200(                  lumen_radius_in_microns_range, microns_per_voxel, ...
                                      space_strel_apothem, max_voxels_per_node, energy_upper_bound, ...
                                    [ data_directory, energy_handle, ROI_names{ ROI_index  }]);
                            
        path_to_vertices      = [ vector_directory, vertices_handle, ROI_names{ ROI_index }];
        
        [ vertex_space_subscripts, vertex_scale_subscripts, vertex_energies ]                       ...
                            = crop_vertices_V200( vertex_space_subscripts, vertex_scale_subscripts, ...
                                                                                   vertex_energies, ... 
                                             length_dilation_ratio * lumen_radius_in_microns_range, ...
                                                                  microns_per_voxel, size_of_image  );
                                           
        [ vertex_energies, sorted_indices ] = sort( vertex_energies );
                
        vertex_space_subscripts = vertex_space_subscripts( sorted_indices, : );
        vertex_scale_subscripts = vertex_scale_subscripts( sorted_indices    );
        
%         runtime_measurement_string = evalc( 'toc' );
%          
%         digit_indices = regexp( runtime_measurement_string, '\d' );        
%         
%         vertices_runtime_in_seconds = sscanf( runtime_measurement_string( digit_indices( 1 ) : end ), '%f' );

        vertices_runtime_in_seconds = toc ;
        
        disp([ 'Runtime for vertices workflow for image ', ROI_names{ ROI_index }( 2 : end ), ' was ', num2str( round( vertices_runtime_in_seconds )), ' seconds' ])
                
        save(   path_to_vertices             , ...
                'vertex_space_subscripts'    , ...
                'vertex_scale_subscripts'    , ...
                'vertex_energies'            , ...
                'vertices_runtime_in_seconds'  );
            
        if exporting_vertex_curation
            % run this line to export curation inputs to a third party
            save([ path_to_saved_curation, '_curator_inputs' ], 'vertex_energies', ...
                             'vertex_space_subscripts', 'vertex_scale_subscripts', ...
                             'lumen_radius_in_microns_range', 'microns_per_voxel', ...
                                'path_to_original_data', 'path_to_saved_curation', ...
                                                             'path_to_energy_data' )            

        end         
    end % FOR ROI
    
    % Restore warnings saying that temporary variables will be removed on each iteration in PARFOR
    warning( 'on', 'MATLAB:mir_warning_maybe_uninitialized_temporary' )        
    
    production_times{ 3 } = attempted_production_times{ 3 };    
    
    save( path_to_workflow_settings, ...
          'production_times'       , ...
          '-append'                  );    
    
end % IF productive
%% ------------------------------------------ Curation -------------------------------------------- 

for ROI_index = ROI_index_range
    
    switch vertex_curation
        
        case { 'manual', 'auto', 'machine-manual', 'machine-auto' }
            
            tic
            
            path_to_energy_data             = [     data_directory,               energy_handle, ROI_names{ ROI_index }]; % mat file path        
            path_to_original_data           = [     data_directory,        original_data_handle, ROI_names{ ROI_index }]; % h5  file path
            path_to_vertices                = [   vector_directory,             vertices_handle, ROI_names{ ROI_index }]; %  vectors path
            path_to_curated_vertices        = [   vector_directory, 'curated_', vertices_handle, ROI_names{ ROI_index }]; %  vectors path    
            path_to_saved_curation          = [ curation_directory,             vertices_handle, ROI_names{ ROI_index }]; % logicals path

            original_file_info = h5info( path_to_original_data );

%             original_image = h52mat( path_to_original_data );
% 
%             intensity_limits = quantile( double( original_image( : )), [ 0.01, 0.99 ]);
% 
%             clear original_image                

            load([ path_to_energy_data, '.mat' ])
            load(  path_to_vertices             )
            load(  path_to_energy_settings      )
                
        otherwise % do nothing

    end

    isMachine = false;
    isManual = false;
    isAuto = false;
   
    switch vertex_curation
        
        case 'machine-auto'
            isMachine = true;
            isAuto = true;
            
        case 'machine-manual'
            isMachine = true;
            isManual = true;
            
        case 'manual'
            isManual = true;
            
        case 'auto'
            isAuto = true;
            
    end
        
    if isMachine
        
        disp([ 'Running machine vertex_curation for image ', ROI_names{ ROI_index }( 2 : end )])            

        path_to_vertex_features = [ curation_directory, vertices_handle, '_featureSet_', ROI_names{ ROI_index }, '.mat'];
        
        if exist(path_to_vertex_features,'file'), delete( path_to_vertex_features ); end                
        
%         if ~exist(path_to_vertex_features,'file')
            vertex_info_extractor({'uncurated'}, {'derivatives'}, {data_directory(1:end-5)}, path_to_vertex_features, ROI_index, time_stamp );
%         end
        
        load(path_to_vertex_features,'vertexFeaturePool')        
        
        vertex_energies = - 1000 * vertexCuratorNetwork_V3(simpleFeatureArray(vertexFeaturePool)');
        
    end %IF isMachine
    
    if isManual
        
        disp([ 'Loading vertex_curation interface for image ', ROI_names{ ROI_index }( 2 : end )])            

        [ vertex_energies, vertex_space_subscripts, vertex_scale_subscripts ]                   ...
                  = vertex_curator(                   vertex_energies, vertex_space_subscripts, ...
                                                                       vertex_scale_subscripts, ...
                                         length_dilation_ratio * lumen_radius_in_microns_range, ...
                                                      microns_per_voxel, path_to_original_data, ...
                                                   path_to_saved_curation, path_to_energy_data, ...
                                                  intensity_limits, [ min( vertex_energies ), 0 ]);

    end %IF isManual
    
    if isAuto
            
        disp([ 'Running automated size exclusion vertex_curation for image ', ROI_names{ ROI_index }( 2 : end )])            

        % instead do a volume conflict test and select the best contrast object at conflicts
        [ chosen_vertex_indices ]                                                               ...
                = choose_vertices_V200( vertex_space_subscripts, vertex_scale_subscripts,       ...
                                        vertex_energies, lumen_radius_in_microns_range,         ...
                                        microns_per_voxel, size_of_image, length_dilation_ratio );

        % performing the selections proposed by either choose_vertices or vertex_curator
        vertex_space_subscripts = vertex_space_subscripts( chosen_vertex_indices, : );
        vertex_scale_subscripts = vertex_scale_subscripts( chosen_vertex_indices    );
        vertex_energies         =         vertex_energies( chosen_vertex_indices    );    
        
    end % IF isAuto
    
    switch vertex_curation
        
        case { 'manual', 'auto', 'machine-manual', 'machine-auto' }
            
%             runtime_measurement_string = evalc( 'toc' );
% 
%             digit_indices = regexp( runtime_measurement_string, '\d' );        
% 
%             vertex_curation_runtime_in_seconds = sscanf( runtime_measurement_string( digit_indices( 1 ) : end ), '%f' );

            vertex_curation_runtime_in_seconds = toc ;

            disp([ 'Runtime for vertex_curation workflow for image ', ROI_names{ ROI_index }( 2 : end ), ' was ', num2str( round( vertex_curation_runtime_in_seconds )), ' seconds' ])            
    
            save(   path_to_curated_vertices            , ...
                    'vertex_space_subscripts'           , ...
                    'vertex_scale_subscripts'           , ...
                    'vertex_energies'                   , ...
                    'vertex_curation'                   , ...
                    'vertex_curation_runtime_in_seconds'  );
                
        otherwise % do nothing

    end
end % ROI FOR    
%% --------------------------------------- Visualization ------------------------------------------ 

if visual( 3 )
    
    load( path_to_energy_settings )        
        
    for ROI_index = ROI_index_range
        
        path_to_vertices                   = [        vector_directory,             vertices_handle, ROI_names{ ROI_index }         ];
        path_to_curated_vertices           = [        vector_directory, 'curated_', vertices_handle, ROI_names{ ROI_index }         ]; %  vectors path            
        path_to_vertices_visuals           = [ visual_vector_directory,             vertices_handle, ROI_names{ ROI_index }, '.tif' ];    
        path_to_curated_vertices_visuals   = [ visual_vector_directory, 'curated_', vertices_handle, ROI_names{ ROI_index }, '.tif' ]; %  vectors path            
        path_to_energy_data                = [          data_directory,             energy_handle,   ROI_names{ ROI_index }         ];       
                
        load([ path_to_energy_data, '.mat' ])        
        load(  path_to_vertices             )
        
        [ vertex_energies, sorted_indices ] = sort( vertex_energies, 'descend' );
                
        vertex_space_subscripts = vertex_space_subscripts( sorted_indices, : );
        vertex_scale_subscripts = vertex_scale_subscripts( sorted_indices    );        
                                                     
        visualize_vertices_V200(  vertex_space_subscripts, ...
                                  vertex_scale_subscripts, ...
                                          vertex_energies, ...
                    lumen_radius_in_pixels_range, size_of_image, ...
                                 path_to_vertices_visuals  )
                       
        try
                             
            load( path_to_curated_vertices        )

            visualize_vertices_V200(  vertex_space_subscripts, ...
                                      vertex_scale_subscripts, ...
                                              vertex_energies, ...
                        lumen_radius_in_pixels_range, size_of_image, ...
                             path_to_curated_vertices_visuals  )
                     
        end % TRY
    end % FOR ROI
end % IF visual

%% 4: get the edges by tracing from each vertex to some other vertices
%% ----------------------------------------- Production ------------------------------------------- 

if productive( 4 )

	disp( 'Running Edges workflow...' )                
    
    % Suppress warnings saying that temporary variables will be removed on each iteration in PARFOR
    warning( 'off', 'MATLAB:mir_warning_maybe_uninitialized_temporary' )    
                      
    for ROI_index = ROI_index_range
        tic
        
        path_to_vertices              = [ vector_directory,             vertices_handle, ROI_names{ ROI_index }]; %  vectors path                    
        path_to_curated_vertices      = [ vector_directory, 'curated_', vertices_handle, ROI_names{ ROI_index }]; %  vectors path            
        path_to_edges                 = [ vector_directory,                edges_handle, ROI_names{ ROI_index }];  
        path_to_energy_data           = [   data_directory,               energy_handle, ROI_names{ ROI_index }];       
        
        load([ path_to_energy_data, '.mat' ])
        
        try load( path_to_curated_vertices )
        
        catch
            
            disp([ 'No vertex curation was found at ', path_to_curated_vertices ])
            
            disp([ 'Extracting edges from uncurated vertices at ', path_to_vertices ])
            
            load( path_to_vertices )
        
        end
                
%         [ edges2vertices, edge_lengths, edge_space_subscripts, edge_scale_subscripts, edge_energies ]  ...
%             = get_edges_V203( lumen_radius_in_pixels_range, vertex_scale_subscripts,                   ...
%                               vertex_space_subscripts, space_strel_apothem_edges,                      ...
%                               max_edge_length_per_origin_radius, number_of_edges_per_vertex,           ...
%                               edge_walk_temperature,                                                   ...
%                               data_directory, [ energy_handle, ROI_names{ ROI_index }]                 );
                          
        [ edges2vertices, ~, edge_space_subscripts, edge_scale_subscripts, edge_energies ]             ...
                   = get_edges_V204(                 lumen_radius_in_microns_range, microns_per_voxel, ...
                              length_dilation_ratio_vertices, vertex_scale_subscripts, vertex_space_subscripts, ...
                                         space_strel_apothem_edges, max_edge_length_per_origin_radius, ...
                                                           number_of_edges_per_vertex, data_directory, ...
                                                             [ energy_handle, ROI_names{ ROI_index }]  );
                                                                                      
        % clean up the output to just keep the trajectories that found a neighbor and to only keep
        % the "best" trajectory from each pair of vertices. Choosing the better trajectory from A
        % to B and B to A.  !!! Any direct visualization of the trajectory variability % should be
        % done before this step !!!
        
%         is_keeping_only_mutual_edge_pairs = number_of_edges_per_vertex >= 4 ;
        is_keeping_only_mutual_edge_pairs = false ;
        
        [ edges2vertices, original_edge_indices ] = clean_edge_pairs( edges2vertices, edge_energies, is_keeping_only_mutual_edge_pairs );
        
        edge_space_subscripts  =  edge_space_subscripts( original_edge_indices );
        edge_scale_subscripts  =  edge_scale_subscripts( original_edge_indices );
        edge_energies          =          edge_energies( original_edge_indices );             

        % preserve the unsmoothed edges.  Only use the smoothed edges for the crop_edges function.
        % Need unsmoothed edges for many clean_edges functions such as for identifying children
        % parent relationships in clean_edge_cycles
        edge_space_subscripts_unsmoothed = edge_space_subscripts ;
        edge_scale_subscripts_unsmoothed = edge_scale_subscripts ;
                edge_energies_unsmoothed = edge_energies         ;
                     
%         sigma_edge_smoothing = 0.5 ;
%         sigma_edge_smoothing = 0 ;
                          
        if sigma_edge_smoothing
            
            [ edge_space_subscripts, edge_scale_subscripts, edge_energies ]                               ...
                                         = smooth_edges_V2( edge_space_subscripts, edge_scale_subscripts, ...
                                                            edge_energies,                                ...
                                                            sigma_edge_smoothing, ...
                                                            lumen_radius_in_microns_range, microns_per_voxel );         
                                                        
        end
                          
        [ excluded_edges_logical ]     ...
                          = crop_edges_V200( edge_space_subscripts, edge_scale_subscripts, edge_energies,       ...
                                             lumen_radius_in_microns_range, microns_per_voxel,  ...
                                             size_of_image                                                      );
                                         
                          edges2vertices( excluded_edges_logical, : ) = [ ];
                          
%                    edge_space_subscripts( excluded_edges_logical )    = [ ];
%                    edge_scale_subscripts( excluded_edges_logical )    = [ ];
%                            edge_energies( excluded_edges_logical )    = [ ];
       
        edge_space_subscripts_unsmoothed( excluded_edges_logical )    = [ ];
        edge_scale_subscripts_unsmoothed( excluded_edges_logical )    = [ ];
                edge_energies_unsmoothed( excluded_edges_logical )    = [ ];
                
        edge_space_subscripts = edge_space_subscripts_unsmoothed ;
        edge_scale_subscripts = edge_scale_subscripts_unsmoothed ;
        edge_energies         =         edge_energies_unsmoothed ;        
                                                                               
% %         % instead do a volume conflict test and select the best contrast object at conflicts                                           
% %         [ mean_edge_energies, chosen_edge_indices, edges2vertices ]                                   ...
% %                           = choose_edges_V200( edges2vertices, edge_energies, edge_space_subscripts, edge_scale_subscripts,      ...
% %                                                lumen_radius_in_pixels_range,            ...
% %                                                length_dilation_ratio_vertices,                        ...
% %                                                length_dilation_ratio_edges, size_of_image             );
% % 
% %         edge_space_subscripts            = edge_space_subscripts(            chosen_edge_indices );
% %         edge_scale_subscripts            = edge_scale_subscripts(            chosen_edge_indices );
% %         edge_energies                    = edge_energies(                    chosen_edge_indices );
% % 
% %         edge_space_subscripts_unsmoothed = edge_space_subscripts_unsmoothed( chosen_edge_indices );
% %         edge_scale_subscripts_unsmoothed = edge_scale_subscripts_unsmoothed( chosen_edge_indices );
% %         edge_energies_unsmoothed         = edge_energies_unsmoothed(         chosen_edge_indices );            
% % 

        chosen_edge_indices = clean_edges_cycles( edges2vertices );

        edges2vertices        = edges2vertices( chosen_edge_indices, : );

        edge_space_subscripts = edge_space_subscripts( chosen_edge_indices );
        edge_scale_subscripts = edge_scale_subscripts( chosen_edge_indices );
        edge_energies         = edge_energies(         chosen_edge_indices );

        chosen_edge_indices = clean_edges_orphans( edge_space_subscripts, size_of_image, vertex_space_subscripts );

        edges2vertices        = edges2vertices( chosen_edge_indices, : );

        edge_space_subscripts = edge_space_subscripts( chosen_edge_indices );
        edge_scale_subscripts = edge_scale_subscripts( chosen_edge_indices );
        edge_energies         = edge_energies(         chosen_edge_indices );

%         chosen_edge_indices = clean_edges_vertex_degree_excess( edges2vertices, 4 ); % 4 is max degree of junction allowed
% 
%         edges2vertices        = edges2vertices( chosen_edge_indices, : );
%         % 
%         edge_space_subscripts = edge_space_subscripts( chosen_edge_indices );
%         edge_scale_subscripts = edge_scale_subscripts( chosen_edge_indices );
%         edge_energies         = edge_energies(         chosen_edge_indices );

        % add vertices where children edges meet their parents
        [ edge_space_subscripts, edge_scale_subscripts, edge_energies, edges2vertices, vertex_space_subscripts, vertex_scale_subscripts, vertex_energies ] ...
                     = add_vertices_to_edges( edge_space_subscripts, edge_scale_subscripts, edge_energies, edges2vertices, vertex_space_subscripts, vertex_scale_subscripts, vertex_energies, size_of_image );
                 
        edge_space_subscripts = cellfun( @double, edge_space_subscripts, 'UniformOutput', false );

        if sigma_edge_smoothing

            [ edge_space_subscripts, edge_scale_subscripts, edge_energies ]                               ...
                                         = smooth_edges_V2( edge_space_subscripts, edge_scale_subscripts, ...
                                                            edge_energies,                                ...
                                                            sigma_edge_smoothing, ...
                                                            lumen_radius_in_microns_range, microns_per_voxel );
                                                        
        end

        mean_edge_energies = get_edge_metric( edge_energies );
                 
%         runtime_measurement_string = evalc( 'toc' );
% 
%         digit_indices = regexp( runtime_measurement_string, '\d' );        
% 
%         edges_runtime_in_seconds = sscanf( runtime_measurement_string( digit_indices( 1 ) : end ), '%f' );

        edges_runtime_in_seconds = toc ;

        disp([ 'Runtime for edges workflow for image ', ROI_names{ ROI_index }( 2 : end ), ' was ', num2str( round( edges_runtime_in_seconds )), ' seconds' ])
                                         
        % saving edge outputs
        save(  path_to_edges            , ...
              'edge_energies'           , ...
         'mean_edge_energies'           , ...
              'edge_space_subscripts'   , ... 
              'edge_scale_subscripts'   , ...
     ...         'edge_lengths'            , ...
              'edges2vertices'          , ...
            'vertex_space_subscripts'   , ...
            'vertex_scale_subscripts'   , ...
            'vertex_energies'           , ...
              'edges_runtime_in_seconds'  );
          
    end % FOR ROI
    
    % Restore warnings saying that temporary variables will be removed on each iteration in PARFOR
    warning( 'on', 'MATLAB:mir_warning_maybe_uninitialized_temporary' )    
    
    production_times{ 4 } = attempted_production_times{ 4 };
    
    save( path_to_workflow_settings, ...
          'production_times'       , ...
          '-append'                  );    
    
end % IF productive          
%% ------------------------------------------ Curation -------------------------------------------- 

for ROI_index = ROI_index_range

    switch edge_curation

        case { 'auto', 'machine-auto', 'machine-manual', 'manual' }
            
            tic

            path_to_energy_data             = [     data_directory,               energy_handle, ROI_names{ ROI_index }]; % h5 file path        
            path_to_original_data           = [     data_directory,        original_data_handle, ROI_names{ ROI_index }]; % h5 file path
            path_to_edges                   = [   vector_directory,                edges_handle, ROI_names{ ROI_index }];  
%             path_to_curated_vertices        = [   vector_directory, 'curated_', vertices_handle, ROI_names{ ROI_index }];          
            path_to_curated_edges           = [   vector_directory, 'curated_',    edges_handle, ROI_names{ ROI_index }]; %  vectors path    
            path_to_saved_curation          = [ curation_directory,                edges_handle, ROI_names{ ROI_index }]; % logicals path

            original_file_info = h5info( path_to_original_data );

%             original_image = h52mat( path_to_original_data );
% 
%             intensity_limits = quantile( double( original_image( : )), [ 0.01, 0.99 ]);
% 
%             clear original_image            

            load(  path_to_energy_settings      )
            load([ path_to_energy_data, '.mat' ])
%             load(  path_to_curated_vertices     )
            load(  path_to_edges                )            

        otherwise % do nothing

    end

    isMachine = false;
    isManual = false;
    switch edge_curation
        
        case 'machine-auto'
            isMachine = true;
            
        case 'machine-manual'
            isMachine = true;
            isManual = true;
            
        case 'manual'
            isManual = true;
    end %SWITCH edge_curation
    
    if isMachine
        
        disp([ 'Running machine edge_curation for image ', ROI_names{ ROI_index }( 2 : end )])            

        path_to_edge_features = [ curation_directory, edges_handle, '_featureSet_', ROI_names{ ROI_index }, '.mat'];
        
        if exist(path_to_edge_features,'file'), delete( path_to_edge_features ); end        
        
%         if ~exist(path_to_edge_features,'file')
            edge_info_extractor({'uncurated'}, {'all'}, {data_directory(1:end-5)}, path_to_edge_features, ROI_index, time_stamp );
%         end
               
        load(path_to_edge_features,'edgeFeaturePool') 
        
        mean_edge_energies = - 1000 * edgeCuratorNetwork_V4_20(simpleFeatureArray(edgeFeaturePool)') ;        
        
        edge_energies = cellfun(@(x,y) x*ones(size(y)), num2cell( mean_edge_energies ), edge_energies, 'UniformOutput', false);
        
%         edge_numels = cellfun( @numel, edge_energies );
%         
%         edge_energies = cellfun(@(x,y) x*ones(y,1), num2cell( mean_edge_energies ), num2cell( edge_numels ), 'UniformOutput', false);
        
        
    end %IF isMachine
        
    if isManual

        disp([ 'Loading edge_curation interface for image ', ROI_names{ ROI_index }( 2 : end )])            

%         if ~ exist( 'length_dilation_ratio',       'var' ), length_dilation_ratio       = 1;   end
        if ~ exist( 'length_dilation_ratio_edges', 'var' ), length_dilation_ratio_edges = 2/3; end
        
        [ mean_edge_energies, edges2vertices, edge_energies, edge_space_subscripts, edge_scale_subscripts, ~, ~, ~, ~, vertex_space_subscripts, vertex_scale_subscripts ] ...
                                  = edge_curator( edge_energies, edge_space_subscripts, edge_scale_subscripts,          ...
                                                  edges2vertices, vertex_space_subscripts,                              ...
                                                  vertex_scale_subscripts, lumen_radius_in_microns_range,               ...
                                                  1, length_dilation_ratio_edges,          ...
                                                  microns_per_voxel, path_to_original_data, path_to_saved_curation,     ...
                                                  path_to_energy_data, intensity_limits, [ min( mean_edge_energies ), 0 ], 1 );

        number_of_vertices_added_in_curator = numel(vertex_scale_subscripts) - numel(vertex_energies);
        vertex_energies = [vertex_energies; -Inf*ones(number_of_vertices_added_in_curator,1)];

%       edge_lengths = []; % put this in the curator function

    end %IF isManual

    switch edge_curation

        case { 'auto', 'machine-auto', 'machine-manual', 'manual' }

%             runtime_measurement_string = evalc( 'toc' );
% 
%             digit_indices = regexp( runtime_measurement_string, '\d' );        
% 
%             edge_curation_runtime_in_seconds = sscanf( runtime_measurement_string( digit_indices( 1 ) : end ), '%f' );
%             
            edge_curation_runtime_in_seconds = toc ;

            disp([ 'Runtime for edge_curation workflow for image ', ROI_names{ ROI_index }( 2 : end ), ' was ', num2str( round( edge_curation_runtime_in_seconds )), ' seconds' ])

            save( path_to_curated_edges         , ...
                                 'edge_energies', ...
                         'mean_edge_energies', ...
                         'edge_space_subscripts', ...
                         'edge_scale_subscripts', ...
                                'edges2vertices', ...
              ...                    'edge_lengths', ...
                       'vertex_space_subscripts', ...
                       'vertex_scale_subscripts', ...
                       'vertex_energies'        , ...
                                 'edge_curation', ...
              'edge_curation_runtime_in_seconds'  );      

        otherwise % do nothing

    end

end % ROI FOR  
%% --------------------------------------- Visualization ------------------------------------------ 
% path_to_visual_edges_settings = [ settings_directory, 'visual_', edges_handle ];

if visual( 4 )
    
    load( path_to_energy_settings ) 
          
    for ROI_index = ROI_index_range
        
        path_to_original_data                  = [          data_directory,     original_data_handle,                 ROI_names{ ROI_index }         ]; % h5  file path                
        path_to_energy_data                    = [          data_directory,            energy_handle,                 ROI_names{ ROI_index }         ];
        
        path_to_edges                          = [        vector_directory,             edges_handle,                 ROI_names{ ROI_index }         ]; %  vectors path 
        path_to_curated_edges                  = [        vector_directory, 'curated_', edges_handle,                 ROI_names{ ROI_index }         ]; %  vectors path  
        
	    path_to_decomposed_visual_file         = [ visual_vector_directory,             edges_handle, '_decomposed',  ROI_names{ ROI_index }, '.tif' ];
	    path_to_curated_decomposed_visual_file = [ visual_vector_directory, 'curated_', edges_handle, '_decomposed',  ROI_names{ ROI_index }, '.tif' ];

        path_to_spheres_visual_file            = [ visual_vector_directory,             edges_handle, '_spheres',     ROI_names{ ROI_index }, '.tif' ];
        path_to_curated_spheres_visual_file    = [ visual_vector_directory, 'curated_', edges_handle, '_spheres',     ROI_names{ ROI_index }, '.tif' ];

        path_to_centerline_visual_file         = [ visual_vector_directory,             edges_handle, '_centerlines', ROI_names{ ROI_index }, '.tif' ];
        path_to_curated_centerline_visual_file = [ visual_vector_directory, 'curated_', edges_handle, '_centerlines', ROI_names{ ROI_index }, '.tif' ];       
        
        load([ path_to_energy_data, '.mat' ])
        load(  path_to_edges                )         

	if ~ exist( 'intensity_limits', 'var' )

		original_image = h52mat( path_to_original_data );

		intensity_limits = quantile( double( original_image( : )), [ 0.01, 0.99 ]);

		clear original_image        

	end % if intensity limits doesn't exist

        edge_subscripts = cellfun(      @(       space_subscripts,      scale_subscripts  )  ...
                                   double([      space_subscripts,      scale_subscripts ]), ...
                                            edge_space_subscripts, edge_scale_subscripts,    ...
                                                  'UniformOutput', false                     );

        visualize_edges_V180( edge_subscripts, mean_edge_energies, lumen_radius_in_pixels_range,  ...
                              size_of_image, path_to_spheres_visual_file,                      ...
                              path_to_centerline_visual_file                                      )
                          
%         decomposed_edge_subscripts = cell2mat( edge_subscripts );
%         decomposed_edge_energies   = cell2mat( edge_energies   );
%         
%         [ decomposed_edge_energies, sorted_indices ] = sort( decomposed_edge_energies, 'descend' );
%         
%         decomposed_edge_subscripts = decomposed_edge_subscripts( sorted_indices, : );
%                           
%         visualize_vertices_V200(  round( decomposed_edge_subscripts( :, 1 : 3 )), ...
%                                          decomposed_edge_subscripts( :,   4   ), ...
%                                          decomposed_edge_energies,               ...
%                                     lumen_radius_in_pixels_range, size_of_image, ...
%                                         path_to_decomposed_visual_file           )                          
                          
                          
                          
%         visualize_depth_via_color_V200( edge_subscripts, mean_edge_energies, lumen_radius_in_pixels_range, ...
%                                                                      1, data_directory, ...
%                                                     [ original_data_handle, ROI_names{ ROI_index }], ...
%                                         [ 1, size_of_image( 1 )], [ 1, size_of_image( 2 )], [ 1, size_of_image( 3 )], max( mean_edge_energies ), intensity_limits )                                  
                    
        try
                                    
            load( path_to_curated_edges      )        

            edge_subscripts = cellfun(      @(       space_subscripts,      scale_subscripts  )  ...
                                       double([      space_subscripts,      scale_subscripts ]), ...
                                                edge_space_subscripts, edge_scale_subscripts,    ...
                                                      'UniformOutput', false                     );

            visualize_edges_V180( edge_subscripts, mean_edge_energies, lumen_radius_in_pixels_range, ...
                                  size_of_image, path_to_curated_spheres_visual_file, ...
                                  path_to_curated_centerline_visual_file                              )

%             visualize_depth_via_color_V200( edge_subscripts, mean_edge_energies, lumen_radius_in_pixels_range, ...
%                                                                          1, data_directory, ...
%                                                         [ original_data_handle, ROI_names{ ROI_index }], ...
%                                             [ 1, size_of_image( 1 )], [ 1, size_of_image( 2 )], [ 1, size_of_image( 3 )], max( mean_edge_energies ), intensity_limits   )
                                        
            decomposed_edge_subscripts = cell2mat( edge_subscripts );
            decomposed_edge_energies   = cell2mat( edge_energies   );

            [ decomposed_edge_energies, sorted_indices ] = sort( decomposed_edge_energies, 'descend' );

            decomposed_edge_subscripts = decomposed_edge_subscripts( sorted_indices, : );

            visualize_vertices_V200(  round( decomposed_edge_subscripts( :, 1 : 3 )), ...
                                             decomposed_edge_subscripts( :,   4   ), ...
                                             decomposed_edge_energies,               ...
                                        lumen_radius_in_pixels_range, size_of_image, ...
                                            path_to_curated_decomposed_visual_file   )
                                        
        end % TRY
    end % FOR ROI
end % IF visual
%% ------------------------------------------ Deletion -------------------------------------------- 

if forgetful( 1 )
    for ROI_index = ROI_index_range
        
        delete([ data_directory, original_data_handle, ROI_names{ ROI_index }])
        
    end % FOR ROI
end % IF forgetul
%% ------------------------------------------ Deletion -------------------------------------------- 

% delete energy data
if forgetful( 2 )
    
    for ROI_index = ROI_index_range
        
%         load( path_to_downsample_settings )
        
        path_to_energy_data = [ data_directory, energy_handle, ROI_names{ ROI_index }];  

        delete( path_to_energy_data )
    
    end % FOR ROI
end % IF forgetful

%% 5: assign strand edges and bifurcation vertices to the network                                   

% exporting_flow_field = special_output( 2 );

% exporting_casX_file = special_output( 7 );
%% ----------------------------------------- Production ------------------------------------------- 

if productive( 5 )

	disp( 'Running Network workflow...' )
    
    load( path_to_energy_settings )
    
    for ROI_index = ROI_index_range
        tic
        
%         path_to_curated_vertices          = [ vector_directory,           'curated_', vertices_handle, ROI_names{ ROI_index }];          
        path_to_curated_edges             = [ vector_directory,           'curated_',    edges_handle, ROI_names{ ROI_index }]; %  vectors path    
        path_to_network                   = [ vector_directory,                        network_handle, ROI_names{ ROI_index }];          
        path_to_energy_data               = [   data_directory,                         energy_handle, ROI_names{ ROI_index }]; % .mat file path
        
        load([ path_to_energy_data, '.mat' ])        
        % !!!! consider try catch statement in analogy to loading the vertices at the edge step
%         load(  path_to_curated_vertices     )
        load(  path_to_curated_edges        )
        
%         chosen_edge_indices = clean_edges_hairs( edges2vertices );
% 
%         edges2vertices        = edges2vertices( chosen_edge_indices, : );
% 
%         edge_space_subscripts    = edge_space_subscripts( chosen_edge_indices );
%         edge_scale_subscripts    = edge_scale_subscripts( chosen_edge_indices );
%         edge_energies            = edge_energies(         chosen_edge_indices );        
                
        [ bifurcation_vertices, vertex_indices_in_strands,                                            ...
           edge_indices_in_strands, end_vertices_of_strands ]      = get_network_V190( edges2vertices );
                                            
        % sort the strand output
        [ vertex_indices_in_strands, edge_indices_in_strands, edge_backwards_in_strands ]           ...
                                      = sort_network_V180( edges2vertices, end_vertices_of_strands, ...
                                                           edge_indices_in_strands                  );
                                                       
        strands2vertices = [ cellfun( @( x ) x(  1  ), vertex_indices_in_strands ), ...
                             cellfun( @( x ) x( end ), vertex_indices_in_strands )  ];
                                                       
        edge_subscripts = cellfun(      @(       space_subscripts,      scale_subscripts  )  ...
                                   double([      space_subscripts,      scale_subscripts ]), ...
                                            edge_space_subscripts, edge_scale_subscripts,    ...
                                                  'UniformOutput', false                     );

        [ strand_subscripts, strand_energies ]                                                                           ...
                = get_strand_objects( edge_subscripts, edge_energies, edge_indices_in_strands, edge_backwards_in_strands );

        strand_space_subscripts = cellfun( @( x ) x( :, 1 : 3 ), strand_subscripts, 'UniformOutput', false );
        strand_scale_subscripts = cellfun( @( x ) x( :,   4   ), strand_subscripts, 'UniformOutput', false );

        if sigma_strand_smoothing
        
            [ strand_space_subscripts, strand_scale_subscripts, strand_energies ]  ...
                                         = smooth_edges_V2( strand_space_subscripts, strand_scale_subscripts, ...
                                                            strand_energies,                                  ...
                                                            sigma_strand_smoothing,     ...
                                                            lumen_radius_in_microns_range, microns_per_voxel  );
                                                        
        end % IF smoothing
                                                    
        [ vessel_directions ] = get_vessel_directions_V3( strand_space_subscripts, microns_per_voxel );
                                                                                                                                                                                                                    
        strand_subscripts = cellfun( @( x, y ) [ x, y ], strand_space_subscripts, strand_scale_subscripts, 'UniformOutput', false );
        
        mean_strand_energies = get_edge_metric( strand_energies );        

%         runtime_measurement_string = evalc( 'toc' );
% 
%         digit_indices = regexp( runtime_measurement_string, '\d' );        
% 
%         network_runtime_in_seconds = sscanf( runtime_measurement_string( digit_indices( 1 ) : end ), '%f' );

        network_runtime_in_seconds = toc ;

        disp([ 'Runtime for network workflow for image ', ROI_names{ ROI_index }( 2 : end ), ' was ', num2str( round( network_runtime_in_seconds )), ' seconds' ])
                
        network_statistics = calculate_network_statistics( strand_subscripts, bifurcation_vertices, lumen_radius_in_microns_range, microns_per_voxel, size_of_image );

        % saving network outputs
        save(  path_to_network             , ...
              'bifurcation_vertices'       , ...
              'strand_subscripts'          , ...
              'strand_energies'            , ...
              'mean_strand_energies'       , ...
              'vessel_directions'          , ...
              'network_runtime_in_seconds' , ...
              'network_statistics'         , ...
              'strands2vertices'             );
          
    end % FOR ROI
    
    production_times{ 5 } = attempted_production_times{ 5 };
    
    save( path_to_workflow_settings, ...
          'production_times'       , ...
          '-append'                  );    
    
end % IF productive
              
path_to_visual_strands_settings = [ settings_directory, 'visual_', network_handle ];
%% --------------------------------------- Visualization ------------------------------------------ 

if visual( 5 )
    
    load( path_to_energy_settings )

    for ROI_index = ROI_index_range
                        
        path_to_original_data              = [   data_directory,    original_data_handle, ROI_names{ ROI_index }]; % h5  file path        
        path_to_energy_data                = [   data_directory,           energy_handle, ROI_names{ ROI_index }];           
%         path_to_curated_vertices           = [ vector_directory, 'curated_', vertices_handle, ROI_names{ ROI_index }];
        path_to_curated_edges              = [ vector_directory, 'curated_', edges_handle, ROI_names{ ROI_index }]; %  vectors path    
        path_to_network                    = [ vector_directory,           network_handle, ROI_names{ ROI_index }];  

        strands_visual_start_vertices_file = [ visual_vector_directory,  network_handle, '_strands', '_start_vertices', ROI_names{ ROI_index }, '.tif' ];
        strands_visual_end_vertices_file   = [ visual_vector_directory,  network_handle, '_strands',   '_end_vertices', ROI_names{ ROI_index }, '.tif' ]; 
        strands_visual_bifurcations_file   = [ visual_vector_directory,  network_handle, '_strands',   '_bifurcations', ROI_names{ ROI_index }, '.tif' ];             
        strands_visual_spheres_file        = [ visual_vector_directory,  network_handle, '_strands',        '_spheres', ROI_names{ ROI_index }, '.tif' ];
%         strands_visual_annuli_file         = [ visual_vector_directory,  network_handle, '_strands',         '_annuli', ROI_names{ ROI_index }, '.tif' ];
        strands_visual_centerline_file     = [ visual_vector_directory,  network_handle, '_strands',    '_centerlines', ROI_names{ ROI_index }, '.tif' ];      
        
        load([ path_to_energy_data, '.mat' ])
%         load(  path_to_curated_vertices     )        
        load(  path_to_curated_edges        )
        load(  path_to_network              ) 
	
	if ~ exist( 'intensity_limits', 'var' )

	%         % !!!!!! remove this (and other instances of) intensity_limits calculation in the future.
	%         % that variable will be saved in the [ path_to_energy_data, '.mat' ]	

	    original_image = h52mat( path_to_original_data );

	    intensity_limits = quantile( double( original_image( : )), [ 0.01, 0.99 ]);

	    clear original_image        

	end % if intensity limits doesn't exist	
        
            
%         logical_nonring_strands = ~ cellfun( 'isempty', end_vertices_of_strands );
% 
%         start_vertex_indices = cellfun( @( v ) v( 1 ), end_vertices_of_strands( logical_nonring_strands ));
%         end_vertex_indices   = cellfun( @( v ) v( 2 ), end_vertices_of_strands( logical_nonring_strands ));
% 
%         %% network start, end, and bifurcation vertices visualization                               
%         visualize_vertices_V190( vertex_space_subscripts( start_vertex_indices, : ), ...
%                                  vertex_scale_subscripts( start_vertex_indices    ), ...
%                                          vertex_energies( start_vertex_indices    ), ...
%                                  lumen_radius_in_pixels_range, sigma_per_size,             ...
%                                  size_of_image, strands_visual_start_vertices_file   )
%                              
%         visualize_vertices_V190( vertex_space_subscripts(   end_vertex_indices, : ), ...
%                                  vertex_scale_subscripts(   end_vertex_indices    ), ...
%                                          vertex_energies(   end_vertex_indices    ), ...
%                                  lumen_radius_in_pixels_range, sigma_per_size,             ...
%                                  size_of_image, strands_visual_end_vertices_file     )

        %% network bifurction vertices visualization                                                
        visualize_vertices_V200( vertex_space_subscripts( bifurcation_vertices, : ), ...
                                 vertex_scale_subscripts( bifurcation_vertices    ), ...
                                         vertex_energies( bifurcation_vertices    ), ...
             lumen_radius_in_pixels_range, size_of_image, strands_visual_bifurcations_file )                             
        %% visualize strands                                                                        
%         edge_strand_indices   = cell2mat( edge_indices_in_strands );
        visualize_edges_V180( strand_subscripts, mean_strand_energies,    ...
                              lumen_radius_in_pixels_range,               ...
                              size_of_image, strands_visual_spheres_file, ...
                              strands_visual_centerline_file              )

    end % FOR ROI
end % IF visual
%% -------------------------------------- Special Outputs ----------------------------------------- 

if any( special_output )
    
    load( path_to_energy_settings )

    for ROI_index = ROI_index_range
                        
        path_to_original_data              = [   data_directory,                 original_data_handle, ROI_names{ ROI_index }]; % h5  file path        
        path_to_energy_data                = [   data_directory,                        energy_handle, ROI_names{ ROI_index }];           
        path_to_flow_field_export          = [   data_directory, 'flow_field_export_', network_handle, ROI_names{ ROI_index }]; % .mat file path for exporting to flow field calculation        
        path_to_curated_edges              = [ vector_directory,           'curated_',   edges_handle, ROI_names{ ROI_index }]; %  vectors path    
        path_to_network                    = [ vector_directory,                       network_handle, ROI_names{ ROI_index }];  
        
        path_to_saved_curation          = [ curation_directory, edges_handle, ROI_names{ ROI_index }];        
        
        strands_visual_spheres_upsampled_file    = [ visual_vector_directory,  network_handle, '_strands_spheres_upsampled',     ROI_names{ ROI_index }, '.tif' ];
        strands_visual_centerline_upsampled_file = [ visual_vector_directory,  network_handle, '_strands_centerlines_upsampled', ROI_names{ ROI_index }, '.tif' ];              
        strands_visual_annuli_file               = [ visual_vector_directory,  network_handle, '_strands_annuli',                ROI_names{ ROI_index }, '.tif' ];
        
        load([ path_to_energy_data, '.mat' ])
        load(  path_to_curated_edges        )
        load(  path_to_network              )
        
        try load( path_to_saved_curation, 'intensity_limits' ), end % if this fails, the intensity limits will be automatically chosen from earlier
        
        %% network statistics (strand histograms )                                                  
        
        if special_output( 1 )
            
            network_histogram_plotter( network_statistics )
            
            number_of_bins = 45 ;
        
            area_histogram_plotter( strand_subscripts, lumen_radius_in_microns_range, microns_per_voxel, size_of_image, number_of_bins );
        
        end        
        %% flow field export                                                                        
        
        if special_output( 2 )
            
            % save function 7/13/18
            vertex_subscripts = [ double( vertex_space_subscripts ), double( vertex_scale_subscripts )];
            
%             edge_subscripts   = cellfun( @( x, y ) [ double( x ), double( y )], edge_space_subscripts, edge_scale_subscripts, 'UniformOutput', false );

            microns_per_pixel_xy = microns_per_voxel( 1 );

            z_per_xy_length_of_pxl_ratio = microns_per_voxel( 3 ) / microns_per_voxel( 1 );
            
%             sigma_per_size = 1 ;
            
%             pixels_per_sigma_range = lumen_radius_in_pixels_range ;
                        
%             sigma_edge_smoothing = sigma_strand_smoothing ;
            
            file_name = network_handle ;
                        
%             % a cutoff at every detected scale
%             tissue_type_cutoffs_in_microns = lumen_radius_in_pixels_range( :, 1 )' * microns_per_pixel_xy ;

            tissue_type_cutoffs_in_microns = [ 5.5 ];

            save( path_to_flow_field_export,                       ...
                      ...            'vertex_indices_in_strands',     ...
                      ...              'edge_indices_in_strands',     ...
                      ...            'edge_backwards_in_strands',     ...
                                  'bifurcation_vertices',          ...
                      ...            'edge_subscripts',               ...
                      ...            'edge_energies',                 ...
                      ...            'mean_edge_energies',            ...
                                       'strand_space_subscripts',  ...
                                       'strand_scale_subscripts',  ...
                                       'strand_energies',          ...
                                             'vessel_directions',  ...
                                  'lumen_radius_in_pixels_range',  ...
                                  'lumen_radius_in_microns_range', ...
     ...   'sigma_edge_smoothing',            ...                                  
                                  'size_of_image',                 ...
                                  'microns_per_pixel_xy',          ...
                                  'microns_per_voxel',             ...
                                  'z_per_xy_length_of_pxl_ratio',  ...
                                  'vertex_subscripts',             ...
                                  'file_name',                     ...
                                  'tissue_type_cutoffs_in_microns' )

            % !! pass in the directory structure to output these fields into the data output
            [ flow_field, tissue_type_image ] = flow_field_subroutine( path_to_flow_field_export );
                             
        end % IF exporting        
        %% 2D depth visualization (red to white gradient color)                                     
        
        is_inverted_original = false ;
        
        are_vectors_opaque = false ;
        
        number_of_slices = 1 ;
        
        size_dilation = 0.25 ;
        
%         y_crop_limits = [ 1, size_of_image( 1 )];
%         x_crop_limits = [ 1, size_of_image( 2 )];

        y_crop_limits = round([ size_of_image( 1 ) / 3, 2 *size_of_image( 1 ) / 3 ]);
        x_crop_limits = round([ size_of_image( 2 ) / 3, 2 *size_of_image( 2 ) / 3 ]);


        slice_index_range = ( 1 : number_of_slices )' ;
        
%         z_crop_limits = [ round( size_of_image( 3 ) / number_of_slices * ( slice_index_range - 1 )) + 1,   ...
%                           round( size_of_image( 3 ) / number_of_slices *   slice_index_range )          ];
                      
        z_crop_limits = 1 + [ round(     ( size_of_image( 3 ) - 1 ) / 3 ), ...
                              round( 2 * ( size_of_image( 3 ) - 1 ) / 3 )  ];                      
        
        if special_output( 3 )

%             visualize_depth_via_color_V200( strand_subscripts, mean_strand_energies,                         ...
%                                             lumen_radius_in_pixels_range, 1,                                 ...
%                                             data_directory, [ original_data_handle, ROI_names{ ROI_index }], ...
%                                             x_crop_limits, x_crop_limits,              ...
%                                             [ 1, round( size_of_image( 3 ) / 4 )], max( mean_strand_energies ), intensity_limits, is_inverted_original )

            for slice_index = 1 : number_of_slices 

                visualize_strands_via_color_V200(                                                    ...
                        strand_subscripts, vessel_directions, mean_strand_energies,                  ...
                        lumen_radius_in_pixels_range, size_dilation,                                             ...
                        data_directory, [ original_data_handle, ROI_names{ ROI_index }],             ...
                        y_crop_limits, x_crop_limits, z_crop_limits( slice_index, : ), ...
                        intensity_limits, is_inverted_original, are_vectors_opaque, 'depth', microns_per_voxel )

            end            
        end
        %% 2D strand visualization (random color assignment)                                        
        
        if special_output( 4 )
        %         number_of_strands = length( edge_indices_in_strands );
        %         
        %         strand_index_range = 1 : number_of_strands ;
        %         
        %         edge_strand_assignments_in_strands = edge_indices_in_strands ;        
        %         
        %         % label the edges of each strand by the strand index of that strand for later colormapping
        %         for strand_index = strand_index_range
        %             
        %             edge_strand_assignments_in_strands{ strand_index }( : ) = strand_index ;
        %             
        %         end % strand FOR
        %                                     
        %         edge_strand_assignments = cell2mat( edge_strand_assignments_in_strands );

            for slice_index = 1 : number_of_slices 

                visualize_strands_via_color_V200(                                                    ...
                        strand_subscripts, vessel_directions, mean_strand_energies,                  ...
                        lumen_radius_in_pixels_range, size_dilation,                                             ...
                        data_directory, [ original_data_handle, ROI_names{ ROI_index }],             ...
                        y_crop_limits, x_crop_limits, z_crop_limits( slice_index, : ), ...
                        intensity_limits, is_inverted_original, are_vectors_opaque, 'strands', microns_per_voxel )

            end
        end
        %% 2D direction visualization                                                               
        
        if special_output( 5 )
            
            for slice_index = 1 : number_of_slices 

                visualize_strands_via_color_V200(                                                    ...
                        strand_subscripts, vessel_directions, mean_strand_energies,                  ...
                        lumen_radius_in_pixels_range, size_dilation,                                             ...
                        data_directory, [ original_data_handle, ROI_names{ ROI_index }],             ...
                        y_crop_limits, x_crop_limits, z_crop_limits( slice_index, : ), ...
                        intensity_limits, is_inverted_original, are_vectors_opaque, 'directions', microns_per_voxel )

            end            
        end
        %% 3D strand visualization (random color assignment)                                        
        
        if special_output( 6 )

    %         max_edge_energies( edge_junction_indices ) = 0 ; % to make image without junctions 

            number_of_slices = 1 ;
            
            resolution_factor = 0.5 ;

            for slice_index = 1 : number_of_slices 

                visualize_strands_via_color_3D_V2(                                                     ...
                        strand_subscripts,                                                             ...
...                        mean_strand_energies - max( mean_strand_energies ), microns_per_voxel( 1 ),    ...
                        mean_strand_energies, microns_per_voxel( 1 ),    ...
                        lumen_radius_in_pixels_range, resolution_factor,                                               ...
                        [ round( size_of_image( 1 ) / number_of_slices ) * ( slice_index - 1 ) + 1,    ...
                          round( size_of_image( 1 ) / number_of_slices ) *   slice_index            ], ...
                        [ 1, size_of_image( 2 )], [ 1, size_of_image( 3 )], 0                          )

            end
        end      
        %% casX export                                                                              
        
        if special_output( 7 )
            
            [ point_coordinates, arc_connectivity, arc_diameters ] = strand2casx( vertex_space_subscripts, strands2vertices, strand_subscripts, microns_per_voxel, lumen_radius_in_microns_range );
        
            save([ path_to_network, '_casX' ], 'point_coordinates', 'arc_connectivity', 'arc_diameters' );

        end % IF exporting_casX_file
        %% upsampled rendering                                                                      
        if special_output( 8 )
%         For the noise study: upsample the vectors for rendering the ground truth image before blurring.

%             resolution_factors = [ 1, 1, 5, 3 ];
            resolution_factors = [ round( microns_per_voxel / min( microns_per_voxel )), 3 ];

% %             exerpt from noise_sensitivity_study script
%             resolution_factors = round( 60 * size_of_image_upsampled ./ size_of_image ) / 60 ;

            [ size_of_image_upsampled, lumen_radius_in_pixels_range_upsampled, strand_subscripts_upsampled ]                                           ...
                                                = resample_vectors( lumen_radius_in_pixels_range, resolution_factors, strand_subscripts, size_of_image );

            vessel_wall_thickness_in_voxels = vessel_wall_thickness_in_microns ./ microns_per_voxel .* resolution_factors( 1 : 3 );

            visualize_edges_V180( strand_subscripts_upsampled, mean_strand_energies,    ...
                                  lumen_radius_in_pixels_range_upsampled,               ...
                                  size_of_image_upsampled, strands_visual_spheres_upsampled_file, ...
                                  strands_visual_centerline_upsampled_file                        )
                              
            visualize_edges_annuli( strand_subscripts_upsampled, mean_strand_energies,  ...
                                    lumen_radius_in_pixels_range_upsampled,             ...
                                    vessel_wall_thickness_in_voxels,                    ...
                                    size_of_image_upsampled, strands_visual_annuli_file )

        end
        %% depth-statistics                                                                         
        if special_output( 9 )
            
            number_of_bins = 30 ;
            
            [ z_bin_ave_z, z_bin_length_densities, z_bin_SA_density, z_bin_vol_densities, z_bin_ave_radius ] = calculate_depth_statistics( strand_subscripts, lumen_radius_in_microns_range, microns_per_voxel, size_of_image, number_of_bins );
                                    
        end
        
    end % FOR ROI
end % IF visual

%% Auxiliary Functions                                                                              

    function [      batch_directory, ...
                     data_directory, ...
              visual_data_directory, ...
                   vector_directory, ...
            visual_vector_directory, ...
                 curation_directory, ...
                 settings_directory  ]    = get_directories( root_directory, batch_handle )

    batch_directory = [ root_directory, batch_handle, filesep ];           

             data_directory = [ batch_directory, 'data',           filesep ];
      visual_data_directory = [ batch_directory, 'visual_data',    filesep ];
           vector_directory = [ batch_directory, 'vectors',        filesep ];
    visual_vector_directory = [ batch_directory, 'visual_vectors', filesep ];
         curation_directory = [ batch_directory, 'curations',      filesep ];
         settings_directory = [ batch_directory, 'settings',       filesep ];

    end % FUNCTION get_directories

    function [ validation_flag, visual ] = validate_visual( x )
        
        curative = zeros( 5, 1, 'logical' );
        
        switch vertex_curation
        
            case { 'manual', 'auto', 'machine-manual', 'machine-auto' }

                curative( 3 ) = true ;
            
        end
        
        switch edge_curation
        
            case { 'manual', 'auto', 'machine-manual', 'machine-auto' }

                curative( 4 ) = true ;
            
        end
        
        productive_for_visual = productive | curative ;
        
        visual_values = [{ 'none' }, workflows, { 'productive', 'all' }];     

        validation_flag = true ;

        % IF cell array
        if iscell( x ) 
            
            % IF cell vector or scalar
            if length( size( x )) == 2 && min( size( x )) == 1

                vector_length = length( x );

                % IF cell scalar
                if vector_length == 1
                    
                    visual = validate_one_visual( x{ 1 }, productive_for_visual, visual_values );
                    
                else % cell vector
                    
                    visual = zeros( 5, 1, 'logical' );
                    
                    for input_index = 1 : vector_length

                        visual = visual + validate_one_visual( x{ input_index }, [ ], workflows );
                        
                    end % FOR input index  
                end % IF cell scalar
            else
                
                error( 'Visual parameter value was type cell array, but it was not a cell vector or scalar' )
                
            end
        else % not type cell
            
            if ischar( x ) || isstring( x )
                
                visual = validate_one_visual( x, productive_for_visual, visual_values );
                
            else
                
                error('Visual parameter value was not type cell, string, or char')
                
            end % IF type char or string
        end % IF type cell
        
        function [ visual ] = validate_one_visual( x, productive, visual_values )

            one_validation_result = validatestring( x, visual_values, 'vectorize', 'Visual' ); 

            switch one_validation_result

                case 'all'

                    visual =  ones( 5, 1, 'logical' );       
                    
                case 'productive'
                            
                    visual = productive ;  
                    
                    if new_batch, visual( 1 ) = true ; end

                otherwise % one workflow

                    visual = strcmp( one_validation_result, workflows )' ;
                    
%                     visual = zeros( 5, 1, 'logical' );                            
% 
%                     switch one_validation_result
% 
%                         case workflows{ 1 }, visual( 1 ) = true ;                                
%                         case workflows{ 2 }, visual( 2 ) = true ;
%                         case workflows{ 3 }, visual( 3 ) = true ;
%                         case workflows{ 4 }, visual( 4 ) = true ;
%                         case workflows{ 5 }, visual( 5 ) = true ;
%                             
%                     end % SWITCH
            end % SWITCH
        end % FUNCTION validate_one_visual      
    end % FUNCTION validate_visual

    function [ validation_flag, special_output ] = validate_special_output( x )
        
        SpecialOutput_core_values = { 'histograms', 'flow-field', 'depth', 'strands', 'directions', '3D-strands', 'casX', 'upsampled', 'depth-stats' };
        
        SpecialOutput_values = [{ 'none' }, SpecialOutput_core_values, { 'all' }];

        validation_flag = true ;

        % IF cell array
        if iscell( x ) 
            
            % IF cell vector or scalar
            if length( size( x )) == 2 && min( size( x )) == 1

                vector_length = length( x );

                % IF cell scalar
                if vector_length == 1
                    
                    special_output = validate_one_special_output( x{ 1 }, SpecialOutput_values );
                    
                else % cell vector
                    
                    special_output = zeros( length( SpecialOutput_core_values ), 1, 'logical' );
                    
                    for input_index = 1 : vector_length

                        special_output = special_output + validate_one_special_output( x{ input_index }, SpecialOutput_core_values );
                            
                    end % FOR input index  
                end % IF cell scalar
            else
                
                error( 'SpecialOutput parameter value was type cell array, but it was not a cell vector or scalar' )
                
            end
        else % not type cell
            
            if ischar( x ) || isstring( x )
                
                special_output = validate_one_special_output( x, SpecialOutput_values );
                
            else
                
                error('SpecialOutput parameter value was not type cell, string, or char')
                
            end % IF type char or string
        end % IF type cell
        
        function [ special_output ] = validate_one_special_output( x, special_output_values )

            one_validation_result = validatestring( x, special_output_values, 'vectorize', 'SpecialOutput' ); 

            switch one_validation_result

                case 'all'

                    special_output =  ones( length( SpecialOutput_core_values ), 1, 'logical' );       
                    
                otherwise % one workflow

                    special_output = strcmp( one_validation_result, SpecialOutput_core_values )' ;
                                        
            end % SWITCH
        end % FUNCTION validate_one_special_output      
    end % FUNCTION validate_special_output
end % FUNCTION vectorize

