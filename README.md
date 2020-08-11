# Vectorization-Public
Public (GNU GPL-3.0 license) repository for SLAVV software.

The main function is called 'vectorize_V200' (documentation at the end of this ReadMe document).

Example calls/uses are shown in the 'vectorization_script_....' files 

Enjoy, leave comments/suggestions, download, change, share (please include the 'licence' file), ....

please cite the methods paper: 

https://www.biorxiv.org/content/biorxiv/early/2020/06/16/2020.06.15.151076.full.pdf

bibtex citation:

@article{mihelic2020segmentation,
  title={Segmentation-less, automated vascular vectorization robustly extracts neurovascular network statistics from in vivo two-photon images},
  author={Mihelic, Samuel and Sikora, William and Hassan, Ahmed and Williamson, Michael and Jones, Theresa and Dunn, Andrew},
  journal={bioRxiv},
  year={2020},
  publisher={Cold Spring Harbor Laboratory}
}

Documentation for the main vectorization function, vectorize_V200

<br /> Vectorize_V200 - Samuel Alexander Mihelic - Novemeber 8th, 2018                                  
<br /> VECTORIZE( ) prompts the user at the command window for all required inputs.  It first asks 
<br />     whether to vectorize a new batch of images or continue with a previous batch.  A batch is a
<br />     group of images that the VECTORIZE function organizes together with two properties:
<br />
<br />         1) The same set of input parameters applies to every image in a batch. 
<br /> 
<br />         2) The images in a batch are processed in parallel at each step of the vectorization. 
<br />            (see Methods below for descriptions of the four steps in the vectorization algorithm).
<br />
<br />       If the user continues with a previous batch, VECTORIZE prompts the user to select a previous
<br />       batch folder with data to recycle.
<br />
<br />       Alternatively, if the user starts a new batch, VECTORIZE prompts the user to select a folder
<br />       with some image file(s) to be vectorized.  It makes a new batch folder in a location
<br />       specified by the user.
<br /> 
<br />     In either case, VECTORIZE prompts the user for a few logistical inputs: which vectorization
<br />     step(s) to execute, what previous data or settings (if any) to recycle, which visual(s) to
<br />     output (if any), and whether or not to open a graphical curator interface. It also prompts the
<br />     user for workflow-specific parameters: It displays imported parameters for review, and prompts
<br />     the user for any missing required parameters.  VECTORIZE writes any outputs to the batch
<br />     folder with a time stamp of the form YYMMDD_HHmmss.
<br /> 
<br />   Conventions:  Greater values in the IMAGE_MATRIX correspond to greater vascular signal
<br />                 The IMAGE_MATRIX dimensions correspond to the physical dimensions y, x, and z
<br />                 (1,x,z) is the top  border of the x-y image at height z
<br />                 (y,1,z) is the left border of the x-y image at height z
<br />                 (y,x,1) is the x-y image nearest to the objective
<br />
<br />   Supported input image file types: .tif
<br />
<br /> For in-line function calls that do not require manual interfacing (e.g. for writing wrapper
<br /> functions or for keeping a concise record of VECTORIZE function calls in a script file), see the
<br /> Optional Input, Logistical Parameters, and Workflow Specific Parameters Sections.
<br /> 
<br /> Note:  For more organizational/navigational control over this document in MATLAB:
<br />           1) open the Preferences Window                                       (HOME>>PREFERENCES)
<br />           2) enable Code Folding for Sections              (MATLAB>>Editor/Debugger>>Code Folding)
<br />           3) fold all of the sections in this document                                    (ctrl,+)
<br />
<br /><br /> ------------------------------------------- Overview ------------------------------------------- 
<br />
<br /> The purpose of the vectorization algorithm is to convert a grayscale, 3D image of vasculature (see
<br /> Inputs section) to a vectorized model of the vascular components.  The output model consists of a
<br /> list of many spherical objects with 3D-position (x,y,z), radius (r), and a contrast metric (c, see
<br /> Methods section). These objects are vectors because each object is a 5-tuples of real numbers:
<br /> [x,y,z,r,c].  The output vectors can then be rendered as a 2- or 3-dimensional image at any
<br /> requested resolution, or it could be analyzed for statistical properties such as volume fraction
<br /> or bifurcation density. With these objects in hand, many analyses are greatly simplified. Many
<br /> such demonstrations and visualizations are automatically output (see Outputs section).
<br />
<br /><br /> ---------------------------------------- Optional Input ---------------------------------------- 
<br />
<br /> VECTORIZE( IMAGE_MATRIX ) vectorizes the numerical array IMAGE_MATRIX.  A batch folder is made 
<br />     in the OutputDirectory specified by the user. The user is prompted for the other logistical
<br />     and workflow specific parameters as in the VECTORIZE( ) call.
<br />
<br />   Supported IMAGE_MATRIX variable types: 3D array of doubles
<br /> 
<br />   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   
<br />
<br /> VECTORIZE( IMAGE_MATRICES ) vectorizes each IMAGE_MATRIX in the cell vector IMAGE_MATRICES.  The
<br />     outputs in the batch folder are numbered by the input order of the images in the cell vector.
<br /> 
<br />   Supported IMAGE_MATRICES variable types: Cell vector of 3D array of doubles
<br />
<br />    -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   
<br />
<br /> VECTORIZE( FILE_NAME ) vectorizes the IMAGE_MATRIX(-CES) specified by the path(s) in FILE_NAME.
<br />
<br />   Supported FILE_NAME variable types: character vectors
<br />
<br />   FILE_NAME is an absolute or relative paths to current working folder. Wild card commands (i.e.
<br />       '*' or '**' ) in the FILE_NAME are also supported. For example:
<br />
<br />     VECTORIZE(    '*.tif'  ) vectorizes all .tif files in the current directory.
<br />
<br />     VECTORIZE([ '**', filesep, '*.tif' ]) vectorizes all .tif files in the current directory or any
<br />         subdirectory thereof.
<br /> 
<br />  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -
<br />
<br /> VECTORIZE( FILE_NAMES ) vectorizes the IMAGE_MATRICES specified by the cell vector of FILE_NAMES.
<br />
<br />   Supported FILE_NAMES variable types: Cell vector of character vectors
<br />
<br /><br /> -------------------------------------- Logistical Parameters ----------------------------------- 
<br />
<br /> VECTORIZE( ..., NAME, VALUE ) 
<br />     uses the following NAME/VALUE pair assignments for logistical inputs:
<br />
<br /> ------- NAME -------      -------------------------------- VALUE --------------------------------
<br />
<br /> 'OutputDirectory'         Char or string specifying the folder that contains the output batch
<br />                           folder.  The default is to prompt the user at the command window.
<br /> 
<br /> 'NewBatch'                'prompt' - Prompts the user at the command window.
<br />                           'yes'    - Makes a new batch folder in the OutputDirectory folder. 
<br />                                      This is the only option if the user provides the Optional
<br />                                      Input.
<br />                           'no'     - Writes into an existing batch folder in the OutputDirectory.
<br />
<br /> 'PreviousBatch'           'prompt' (default) - Prompts the user to input a previous batch folder.
<br />                           'none'             - Does not import any existing data or settings.  If
<br />                                                NewBatch is 'no', this is  not an  option.
<br />                           'yyMMdd-HHmmss'    - Imports from the batch_yyMMdd-HHmmss folder in the 
<br />                                                OutputDirectory.
<br />                           'recent'           - Imports from the most recent batch_* folder in the
<br />                                                OutputDirectory.
<br />
<br /> 'PreviousWorkflow'        'prompt' (default) - Prompts the user to input a previous settings file.
<br />                           'none'             - Does not import any existing data or settings.  If
<br />                                                NewBatch is 'no', this is not an option.
<br />                           'yyMMdd-HHmmss'    - Imports from the workflow_yyMMdd-HHmmss folder in
<br />                                                the batch folder. 
<br />                           'recent'           - Imports from the most recent workflow_settings_* 
<br />                                                folder in the batch folder.
<br />
<br /> 'StartWorkflow'           Note: These are listed in order:  Energy is the first process.
<br /> 
<br />                           'prompt' (default) - Prompts the user at the command window.
<br />                           'none'             - Does not run any workflow steps.  The user may run
<br />                                                curatation or visual steps.
<br />                           'next'      - Starts the vectorization process at the   next   step for
<br />                                         the selected PreviousWorkflow.
<br />                           'energy'    - Starts the vectorization process at the  Energy  step.
<br />                                         This is the only option if NewBatch is 'yes'.
<br />                           'vertices'  - Starts the vectorization process at the Vertices step.
<br />                           'edges'     - Starts the vectorization process at the   Edges  step.
<br />                           'network'   - Starts the vectorization process at the  Network step.
<br />
<br /> 'FinalWorkflow'           'prompt'(default) - Prompts the user at the command window.
<br />                           'none'            - Does not run any vectorization processes.  
<br />                                               This is the only option if StartWorkflow is 'none'
<br />                           'one'       - Ends the vectorization process at the StartWorkflow step.
<br />                           'energy'    - Ends the vectorization process at the    Energy     step.
<br />                           'vertices'  - Ends the vectorization process at the   Vertices    step.
<br />                           'edges'     - Ends the vectorization process at the     Edges     step.
<br />                           'network'   - Ends the vectorization process at the    Network    step.
<br />
<br /> 'Visual'                  'none'                 - Does not write visual outputs.
<br />                           'original'             - Writes visuals for just the input images.
<br />                           'energy'               - Writes visuals for just the   Energy step.
<br />                           'vertices'             - Writes visuals for just the Vertices step.
<br />                           'edges'                - Writes visuals for just the    Edges step.
<br />                           'network'              - Writes visuals for just the  Network step.
<br />                           { ... }                - Writes visuals for just the    ...   steps.
<br />                           'productive' (default) - Writes visuals for just the workflow steps.
<br />                                                    being executed at this vectorization call.
<br />                           'all'                  - Writes visuals for all vectorization steps.
<br /> 
<br /> 'SpecialOutput'           'none'           - Does not create any special network outputs.
<br />                           'histograms'     - (defualt) shows strand, length statistic histograms
<br />                           'depth-stats'    - Shows depth-resolved statistics.
<br />                           'flow-field'     - Writes x, y, and z component .tif's of flow field.
<br />                           'depth'          - Shows vectors over raw with color-coded depth.               
<br />                           'strands'        - Shows vectors over raw with color-coded strands.
<br />                           'directions'     - Shows vectors over raw with color-coded direcions.
<br />                           '3D-strands'     - Shows 3D volume rendering with color-coded strands.
<br />                           'casX'           - Creates .casX equivalent representation of strands.
<br />                                              Format .casX is due to LPPD in Chicago.
<br />                           { ... }          - Creates ... special network outputs.
<br />                           'all'            - Creates all special network outputs.
<br /> 
<br /> 'VertexCuration'          'auto'             - All non-overlapping vertices are passed to edges. 
<br />                                                Chooses least energy vector upon volume conflict. 
<br />                           'manual' (default) - Prompts user with a graphical curation interface.
<br />                           'machine-manual'   - Applies neural network categorization and then
<br />                                                prompts user with a graphical interface.
<br />                           'machine-auto'     - Applies neural network categorization and then
<br />                                                all non-overlapping vertices are passed to edges.
<br />                                                Chooses best vector according to the neural network
<br />                                                upon volume conflict.
<br />
<br /> 'EdgeCuration'            'auto'             - All edges are passed to network.
<br />                           'manual' (default) - Prompts user with a graphical curation interface.
<br />                           'machine-manual'   - Applies neural network categorization and then
<br />                                                prompts user with a graphical interface.
<br />                           'machine-auto'     - Applies neural network categorization.  All edges 
<br />                                                are passed to network.
<br /> 
<br /> 'NetworkPath'             'prompt
<br />                           'built-in' (default) - Built-in network ...
<br />                           'train'              - Trains new network from all curation files found
<br />                                                  in the training folder in the vectorization base 
<br />                                                  directory.
<br />                           'yyMMdd-HHmmss'      - Imports network trained at yyMMdd-HHmmss from the 
<br />                                                  network folder in the vectorization base
<br />                                                  directory.
<br />                           'recent'             - Imports network trained most recently from the 
<br />                                                  network folder in the vectorization base
<br />                                                  directory.
<br />
<br /> 'Forgetful'               'none' (default)   - No intermediate data files will be deleted
<br />                           'original'         - Deletes intermediate copies of the input images
<br />                           'energy'           - Deletes intermediate Energy data files
<br />                           'both'             - Deletes both intermediate data files when done
<br />
<br /> 'Presumptive'             false (default) - Prompts user for all required workflow-specific inputs
<br />                           true            - Assumes previous settings or default if no previous
<br /><br /> ------------------------------- Workflow-Specific Parameters ----------------------------------- 
<br />
<br /> VECTORIZE( ..., NAME, VALUE ) 
<br />     uses the NAME/VALUE pair assignments listed below to input workflow-specific parameters.  Each
<br />     parameter is listed below under the first workflow step that it modifies.
<br />   
<br />   Resolving conflicts between NAME/VALUE pairs and imported parameters from the PreviousWorkflow:
<br /> 
<br />       If a NAME/VALUE pais is provided for a parameter that modifies a workflow step that is
<br />       upstream of the StartWorkflow, then that value is ignored and the value from the
<br />       PreviousWorkflow value will remain. A warning is produced.  This must occur because the
<br />       relevant workflow has already been executed and is not scheduled to run again, therefore the
<br />       parameters that modify it will not change.
<br />
<br />       If a NAME/VALUE pair is provided for a parameter that only modifies workflows that are equal
<br />       to or downstream of the StartWorkflow, then that value will overwrite any value that may
<br />       have been retrieved from the PreviousWorkflow.  The overwriting is displayed in the commaind
<br />       window.
<br />
<br /><br />  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - Energy -  -  -  -  -  -  -  -  -  -  -  -  -  
<br /> ---------------- NAME ----------------  ------------------------- VALUE -------------------------
<br /> 
<br /> 'microns_per_voxel'                     Real, positive, three-element vector specifying the voxel
<br />                                         size in microns in y, x, and z dimensions, respectively.
<br />                                         Default: [ 1, 1, 1 ]
<br /> 
<br /> 'radius_of_smallest_vessel_in_microns'  Real, positive scalar specifying the radius of the
<br />                                         smallest vessel to be detected in microns.  Default: 1.5
<br /> 
<br /> 'radius_of_largest_vessel_in_microns'   Real, positive scalar specifying the radius of the largest 
<br />                                         vessel to be detected in microns.  Default: 50
<br /> 
<br /> 'approximating_PSF'                     Logical scalar specifying whether to approximate the PSF
<br />                                         using "Nonlinear Magic: Multiphoton Microscopy in the 
<br />                                         Biosciences" (Zipfel, W.R. et al.).  Default: true
<br /> 
<br /> 'sample_index_of_refraction'            Real, positive scalar specifying the index of refraction 
<br />                                         of the sample.  This parameter is only used if
<br />                                         approximating the PSF.  Default: 1.33
<br /> 
<br /> 'numerical_aperture'                    Real, positive scalar specifying the numerical aperture of
<br />                                         the microscope objective.  Default: 0.95
<br /> 
<br /> 'excitation_wavelength_in_microns'      Real, positive scalar specifying the excitation wavelength
<br />                                         of the laser in microns.  Default: 1.3
<br /> 
<br /> 'scales_per_octave'                     Real, positive scalar specifying the number of vessel 
<br />                                         sizes to detected per doubling of the radius cubed.  
<br />                                         Default: 1.5
<br /> 
<br /> 'max_voxels_per_node_energy'            Real, positive scalar specifying Default: 1e5
<br /> 
<br /> 'gaussian_to_ideal_ratio'               Real scalar between 0 and 1 inclusive specifying the
<br />                                         standard deviation of the Gaussian kernel per the total
<br />                                         object length for objects that are much larger than the
<br />                                         PSF. Default: 1
<br /> 
<br /> 'spherical_to_annular_ratio'            Real scalar between 0 and 1 inclusive specifying the
<br />                                         weighting factor of the spherical pulse over the combined
<br />                                         weights of spherical and annular pulses. This parameter is
<br />                                         only used if gaussian_to_ideal_ratio is strictly less than
<br />                                         unity.  Default: 1
<br />
<br />
<br /><br />   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   Edges   -  -  -  -  -  -  -  -  -  -  -  -  - 
<br /> ---------------- NAME ----------------  ------------------------- VALUE -------------------------
<br /> 
<br /> 'max_edge_length_per_origin_radius'     Real, positive scalar specifying the maximum length of an
<br />                                         edge trace per the radius of the seed vertex. Default: 30
<br />
<br /> 'number_of_edges_per_vertex'            Real, positive integer specifying the maximum number of
<br />                                         edge traces per seed vertex. Default: 4
<br />
<br /><br />  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   Network -  -  -  -  -  -  -  -  -  -  -  -  -  
<br /> ---------------- NAME ----------------  ------------------------- VALUE -------------------------
<br /> 
<br /> 'sigma_strand_smoothing'                Real, non-negative integer specifying the standard
<br />                                         deviation of the Gaussian smoothing kernel per the radius
<br />                                         of the strand at every position along the strand vector.
<br />                                         Default: 1
<br />
<br /><br /> ------------------------------------------- Methods -------------------------------------------- 
<br />
<br /> Vectorization is accomplished in four steps:  (1) energy image formation, (2) vertex extraction,
<br /> (3) edge extraction, and (4) network extraction.  The raw output is a superposition
<br /> of possible models until it is segmented in some way.  The objects are assigned a contrast metric
<br /> based on the values from the energy image, and thresholding them on this value provides direct
<br /> control over the sensitivity and specificity of the vectorization.  Alternatively, the graphical
<br /> curation interface provides a platform for manual segmentation such as local threshold selection
<br /> or point-and-click object selection.
<br /> 
<br /> 1: Energy Image Formation
<br />     Multi-scale (at many pre-defined sizes) gradient and curvature information from the original
<br />     3D image is combined to form a 4-dimensional, multi-scale, centerline-enhanced, image, known
<br />     as the energy image.  Ideally, the voxels with the lowest energy value the energy image will
<br />     be the most likely to be centerline voxels for vessels in the pre-defined size range.  This
<br />     can be visually verified by inspecting the energy*.tif visual output in the visual output
<br />     directory.  The energy image should be very negative right at the vessel centerlines and close
<br />     to zero or positive elsewhere.  The value of the size image at a given voxel shows the index
<br />     of the pre-defined sizes that is most likely assuming a vessel is centered at that voxel.
<br /> 
<br /> 2: Vertex Extraction
<br />     Vertices are extracted as local minima in the 4D energy image (with associated x, y, z,
<br />     radius, and energy value).  This method was inspired by the first part of the SIFT algorithm
<br />     (David Lowe, International Journal of Computer Vision, 2004)).  Ideally, local minima of the
<br />     energy image correspond to voxels that are locally the most likely to be along a vessel
<br />     centerline. The size coordinate is also required to be at a local energy minimum.  In theory,
<br />     the vertices are ordered from most to least likely to exist by the energy values at their
<br />     locations.
<br /> 
<br /> 3: Edge Extraction
<br />     Edges are extracted as voxel to voxel random walks through the (min. projected to 3D) energy
<br />     image.  Therefore edges are lists of spherical objects like vertices.  Edge trajectories seek
<br />     lower energy values and are forced to connect exactly two vertices.  The trajectories between
<br />     vertices are in theory ordered from most to least likely to exist by their mean energy values.
<br /> 
<br /> 4: Network Extraction
<br />     Strands are defined as the sequences of non-branching edges (single random color in the
<br />     colored strands image).  Strands are found by counting the number of edges in the adjacency
<br />     matrix of the vertices.  Strands are the connected components of the adjacency matrix that
<br />     only includes vertices with two edges.  With the strands of the network in hand, we
<br />     equivalently know where the bifurcations in the network are.  Network information unlocks many
<br />     doors.  For instance, we can smooth the positions and sizes of the extracted vectors along
<br />     their strands and approximate local blood flow fields.
