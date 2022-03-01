# **Repository ReadMe**
## Segmentation-Less, Automated, Vascular Vectorization (SLAVV)
#### Public ([GNU GPL-3.0 license](https://github.com/UTFOIL/Vectorization-Public/blob/master/LICENSE)) Code Repository 

For summarizing and analyzing vascular networks from grayscale, volumetric images. Use this software for locating, sizing, and mapping the connectivity of the depicted vascular networks.

### Table of Contents:
- [Software (MATLAB) Overview](#Software-Overview)
- [Curator Interface Tutorial](#tutorial)
- [Vectorized Example Datasets](#Example-Datasets)
- [Software (MATLAB) Documentation](#Documentation)
- [Methodology Manuscript](#Methodology-Manuscript)
- [Application Manuscripts](#Application-Manuscripts)

#### Enjoy, leave comments/suggestions, download, change, share :) Please include the [LICENSE file](https://github.com/UTFOIL/Vectorization-Public/blob/master/LICENSE), ....

#### Please cite the [SLAVV methodology manuscript](https://github.com/UTFOIL/Vectorization-Public#Methodology-Manuscript).

## Software Overview
The [main function](https://github.com/UTFOIL/Vectorization-Public/blob/master/vectorize_V200.m) is called ```vectorize_V200.m``` ([documentation](#Documentation)).

The [supporting functions](https://github.com/UTFOIL/Vectorization-Public/blob/master/source) are in the ```source``` folder.

To vectorize an image using the SLAVV software in MATLAB, simply call [```vectorize_V200.m```](https://github.com/UTFOIL/Vectorization-Public/blob/master/vectorize_V200.m) in the Command Window with no inputs to be prompted for all required inputs.

Accepted input image formats: 
- three-dimensional array of doubles (variable in MATLAB workspace)
- ```*.tif``` (file location)

Example inline calls can be found in the vectorization scripts, ```vectorization_script_*.m```. (e.g. [Example 1](https://github.com/UTFOIL/Vectorization-Public/blob/master/source/vectorization_script_2017MMDD_TxRed_chronic.m), [Example 2](https://github.com/UTFOIL/Vectorization-Public/blob/master/source/vectorization_script_michael.m), ...), and in the [performance sensitivity to image quality testing script](https://github.com/UTFOIL/Vectorization-Public/blob/master/source/noise_sensitivity_study_V2.m).

The resulting vectorization can be exported in the documented file formats: 
- ```*.vmv``` [VessMorphoVis plugin](https://github.com/BlueBrain/VessMorphoVis) to [Blender](https://www.blender.org/) due to [BBP](https://www.epfl.ch/research/domains/bluebrain/), described by [Abdellah et al. in Bioinformatics in 2020](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7355309/). The documentation is in the supplementary file for this publication. 
- ```*.casx``` Casx file format ([documentation](https://github.com/UTFOIL/Vectorization-Public/blob/master/references/casxFileDocumentation.pdf)) for storing vascular network data designed by G. Hartung and A. Linninger, UIC, 2016-2019.

<!--and [documentation](https://github.com/UTFOIL/Vectorization-Public/blob/master/references/btaa461_supplementary_data.pdf)) -->

## Tutorial

<!--[This folder](https://github.com/UTFOIL/Vectorization-Public/blob/master/tutorial/) contains screenshots-->
[Demonstrating](https://github.com/UTFOIL/Vectorization-Public/blob/master/SLAVV_tutorial_readme.md) the use of the graphical curator interface (built-into the SLAVV software) on three large (~1 mm<sup>3</sup>) images of living adult mouse brain vasculature. The three images are all from the same mouse, approximately the same field of view, and timed 2 weeks apart.

## Example Datasets

[This Google Drive folder](https://drive.google.com/drive/folders/1IO9PmpO0W_E1UaedvzdP5_PzHzfYGO_m?usp=sharing) and [this Texas Digital Library repository](https://doi.org/10.18738/T8/NA08NU) contain the three example images as they were processed for the [SLAVV methodology manuscript](https://github.com/UTFOIL/Vectorization-Public#Methodology-Manuscript) (_distinct_ from the three images used in the [Tutorial](https://github.com/UTFOIL/Vectorization-Public/blob/master/SLAVV_tutorial_readme.md)).  
  
Each of the three ```Image *``` directories contains:  
	- The original image under its original title as a ```*.tif```  
	- The wrapper script that was used to vectorize that image in a ```vectorization_script_*.m``` file  
	- All of the vectorization intermediates and outputs in a ```batch_*``` folder

## Documentation 
For the main vectorization function in MATLAB:

[```https://github.com/UTFOIL/Vectorization-Public/blob/master/vectorize_V200.m```](https://github.com/UTFOIL/Vectorization-Public/blob/master/vectorize_V200.m)

Copied here for reference:
```
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
% The purpose of the vectorization algorithm is to convert a grayscale, 3D image of vasculature to a
% vectorized model of the vascular components. With the vectorized objects in hand, many analyses
% are greatly simplified: The output vectors can be rendered as a 2- or 3-dimensional image at any
% requested resolution, or analyzed for statistical properties such as volume fraction or
% bifurcation density. These calculations as well as others can be specified using the SpecialOutput
% input parameter described in the Logistical Parameters section below.
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
%                           'histograms'     - Shows strand, length statistic histograms          (defualt) 
%                           'depth-stats'    - Shows depth-resolved statistics.                   (defualt) 
%                           'flow-field'     - Writes x, y, and z component .tif's of flow field.
%                           'depth'          - Shows vectors over raw with color-coded depth.     (defualt)           
%                           'strands'        - Shows vectors over raw with color-coded strands.   (defualt) 
%                           'directions'     - Shows vectors over raw with color-coded direcions. (defualt) 
%                           '3D-strands'     - Shows 3D volume rendering with color-coded strands.
%                           '3D-directions'  - Shows 3D volume rendering with color-coded strand 
%                                              directions.
%                           'original-stats' - Calculates signal, noise, background, etc, from
%                                              network binary. Saves in the network_*.mat vector file.
%                           'casX'           - Creates .casX equivalent representation of strands.
%                                              Format .casX is due to LPPD in Chicago, IL.
%                           'vmv'            - Creates .vmv equivalent representation of strands. (defualt)
%                                              Format .vmv is due to Blue Brain Project in
%                                              Switzerland.
%                           { ... }          - Creates ... (a list of) special network outputs.
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
% % % !!!!!!!!!!!!!!! not yet implemented
% % % 'NetworkPath'             'prompt
% % %                           'built-in' (default) - Built-in network ...
% % %                           'train'              - Trains new network from all curation files found
% % %                                                  in the training folder in the vectorization base 
% % %                                                  directory.
% % %                           'yyMMdd-HHmmss'      - Imports network trained at yyMMdd-HHmmss from the 
% % %                                                  network folder in the vectorization base
% % %                                                  directory.
% % %                           'recent'             - Imports network trained most recently from the 
% % %                                                  network folder in the vectorization base
% % %                                                  directory.
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
% %   Resolving conflicts between NAME/VALUE pairs and imported parameters from the PreviousWorkflow:
% % 
% %       If a NAME/VALUE pais is provided for a parameter that modifies a workflow step that is
% %       upstream of the StartWorkflow, then that value is ignored and the value from the
% %       PreviousWorkflow value will remain. A warning is produced.  This must occur because the
% %       relevant workflow has already been executed and is not scheduled to run again, therefore the
% %       parameters that modify it will not change.
% %
% %       If a NAME/VALUE pair is provided for a parameter that only modifies workflows that are equal
% %       to or downstream of the StartWorkflow, then that value will overwrite any value that may
% %       have been retrieved from the PreviousWorkflow.  The overwriting is displayed in the commaind
% %       window.
% %
%%  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  - Energy -  -  -  -  -  -  -  -  -  -  -  -  -  
% ---------------- NAME ----------------  ------------------------- VALUE -------------------------
% 
% 'microns_per_voxel'                     Real, positive, three-element vector specifying the voxel
%                                         size in microns in y (up/down), x (left/right), and z {out
%                                         of/into) dimensions, respectively. Default: [ 1, 1, 1 ]
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
% %%   -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   Edges   -  -  -  -  -  -  -  -  -  -  -  -  - 
% % ---------------- NAME ----------------  ------------------------- VALUE -------------------------
% % 
% % 'max_edge_length_per_origin_radius'     Real, positive scalar specifying the maximum length of an
% %                                         edge trace per the radius of the seed vertex. Default: 30
% %
% % 'number_of_edges_per_vertex'            Real, positive integer specifying the maximum number of
% %                                         edge traces per seed vertex. Default: 4
% %
% %%  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -   Network -  -  -  -  -  -  -  -  -  -  -  -  -  
% % ---------------- NAME ----------------  ------------------------- VALUE -------------------------
% % 
% % 'sigma_strand_smoothing'                Real, non-negative integer specifying the standard
% %                                         deviation of the Gaussian smoothing kernel per the radius
% %                                         of the strand at every position along the strand vector.
% %                                         Default: 1
% %
%% ------------------------------------------- Methods -------------------------------------------- 
%
% Vectorization is accomplished in four steps:  (1) energy image formation, (2) vertex extraction,
% (3) edge extraction, and (4) network extraction.  The raw output is a superposition of possible
% models until it is segmented in some way. The vertex and edge objects are assigned a contrast
% metric based on the values from the energy-filtered image, and thresholding them on this value
% provides direct control over the sensitivity and specificity of the vectorization.  Alternatively,
% the graphical curation interfaces (one for the edges, and one for the vertices) provides a
% platform for manual segmentation such as local threshold selection or point-and-click object
% selection.
% 
% 1: Energy Image Formation
%     The input image is linearly (matched-)filtered at many scales to extract curvature (and
%     gradient) information. The Hessian matrix is diagonalized to extract principle curvatures at
%     all voxels (and scales) where the Laplacian is negative (local bright spot in original image).
%     The energy function is an objective function to select for large negative principle curvatures
%     while each curvature is seperately weighted by a symmetry factor using the gradient. The 4D
%     multi-scale energy image is projected across the scale coordinate using a minimum projection
%     to select the most probable scale. The result is an enhancement of the vessel centerlines
%     while simultaneously selecting the scale coordinate (index into the size range LUT input).
% 
% 2: Vertex Extraction
%     Vertices are extracted as local minima in the 4D energy image (with associated x, y, z,
%     radius, and energy value). (This part of the method was inspired by the first part of the SIFT
%     algorithm (David Lowe, International Journal of Computer Vision, 2004)).  The vertex objects
%     are points of high contrast and symmetry (bright spots or bifurcations) along the vessel
%     segments. The labeling is sufficiently dense if there is at least one vertex per strand (see
%     network section for strand definition). The vertices are ordered by the energy values to rank
%     them from most likely to least likely (based on the energy filter model; without knowledge of
%     a ground truth).
% 
% 3: Edge Extraction
%     Edges are extracted as voxel walks through the (3D) energy image.  Therefore edges are
%     1-Dimensional objects such as a list or a trace, where each location along the trace is a
%     spherical object with an energy value (like a vertex).  Each edge walk starts at a vertex and
%     seeks the lowest energy values under the constraint that it must move away from its origin.
%     The trajectories between vertices are ordered by their maximum energy value attained, in order
%     to give an first estimate of which edges are likely to be true (based on the energy filter
%     model; without knowledge of a ground truth) and which edges are likely to be false.
% 
% 4: Network Extraction
%     The final network output is the minimal set of 1-Dimensional objects (strands) that connect
%     all of the bifurcations/endpoints according to the adjacency matrix exracted from the edges
%     and vertices. The strands are like the edges, but in general are longer, and composed of
%     multiple edges (at least 1). (Theerefore, each strand has at each location along its trace a
%     3-space position, radius, and an energy value). With the strands of the network in hand, we
%     equivalently know where the bifurcations in the network are.  Network information unlocks many
%     doors.  For instance, we can smooth the positions and sizes of the extracted vectors along
%     their strands and approximate local blood flow fields.
%
%% ------------------------------------------- Outputs -------------------------------------------- 
%
% TIME_STAMP = VECTORIZE( ... )
%     returns the TIME_STAMP that was assigned to any new data or settings output
%
% [ TIME_STAMP, ROI_NAMES ] = VECTORIZE( ... )
%     also returns the FILE_NAMES that were used as input or the names assigned to the IMAGE_MATRICES in
%     the order that they were passed to VECTORIZE.
%
% The final vectorized model of the network is the minimal set of 1-Dimensional objects (strands)
% that connect all of the bifurcations/endpoints to each other. Each strand has at each location
% along its trace a 3-space position, radius, and an energy value. This output (without the energy
% values) can be exported in .vmv or .casx formats (see the SpecialOutput input parameter).
% 
% Other intermediate outputs are output to the OutputDirectory in the batch folder like so:
% 
%  batch_YYMMDD_hhmmss
%       curations
%       data
%       settings
%       vectors
%       visual_data
%       visual_vectors
%       
% The curations folder contains the automatic and user-requested curation backups. 
% 
% The data folder has the hdf5 files of the original, energy, and intermediate edges. 
% 
% The settings folder holds the record of vectorization executions (time-stamp: YYMMDD_hhmmss), the
% workflows that were requested, and the values of the parameters used.
% 
% The vectors folder holds the vertex, edge, and strand objects. The final network vectors (strands)
% are saved in the internal SLAVV format, or one of the exported formats (.vmv, .casx; see
% SpecialOutput). A statistical summary of the network can be found in the network*.mat file under
% the struct variable network_statistics.
% 
% The visual_vectors folder holds .tif 3D images of vector renderings (with weightings corresponding
% to their energy values or thereof sorting index) in the same resolution as the original for
% overlaying and comparison purposes.
% 
% The visual_data folder holds the .tif 3D images of the original, the energy filtered image (along
% with the scale (index) image), and the intermediate edges. This is also where any figures
% specified for generation by the SpecialOutput input will be saved.
%
```

## Methodology Manuscript

Public Library of Science, Computational Biology (PLOS CompBio) publication:

[```https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009451```](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1009451)

```
@article{mihelic2021segmentation,
  title={Segmentation-Less, Automated, Vascular Vectorization},
  author={Mihelic, Samuel A and Sikora, William A and Hassan, Ahmed M and Williamson, Michael R and Jones, Theresa A and Dunn, Andrew K},
  journal={PLOS Computational Biology},
  volume={17},
  number={10},
  pages={e1009451},
  year={2021},
  publisher={Public Library of Science San Francisco, CA USA}
}
```


BioRxiv pre-print:

[```https://www.biorxiv.org/content/biorxiv/early/2020/06/16/2020.06.15.151076.full.pdf```](https://www.biorxiv.org/content/biorxiv/early/2020/06/16/2020.06.15.151076.full.pdf)

bibtex:
```
@article{mihelic2020segmentation,
  title={Segmentation-less, automated vascular vectorization robustly extracts neurovascular network statistics from in vivo two-photon images},
  author={Mihelic, Samuel and Sikora, William and Hassan, Ahmed and Williamson, Michael and Jones, Theresa and Dunn, Andrew},
  journal={bioRxiv},
  year={2020},
  publisher={Cold Spring Harbor Laboratory}
}
```

## Application Manuscripts

### Mathematical synthesis of the cortical circulation for the whole mouse brain—part II: Microcirculatory closure

[```https://neurophysics.ucsd.edu/publications/micc.12687.pdf```](https://neurophysics.ucsd.edu/publications/micc.12687.pdf)

bibtex:

```
@article{hartung2021mathematical,
  title={Mathematical synthesis of the cortical circulation for the whole mouse brain—part II: Microcirculatory closure},
  author={Hartung, Grant and Badr, Shoale and Mihelic, Samuel and Dunn, Andrew and Cheng, Xiaojun and Kura, Sreekanth and Boas, David A and Kleinfeld, David and Alaraj, Ali and Linninger, Andreas A},
  journal={Microcirculation},
  volume={28},
  number={5},
  pages={e12687},
  year={2021}
}
```

### Effect of vascular structure on laser speckle contrast imaging

[```https://opg.optica.org/boe/fulltext.cfm?uri=boe-11-10-5826```](https://opg.optica.org/boe/fulltext.cfm?uri=boe-11-10-5826)

bibtex:

```
@article{jafari2020effect,
  title={Effect of vascular structure on laser speckle contrast imaging},
  author={Jafari, Chakameh Z and Sullender, Colin T and Miller, David R and Mihelic, Samuel A and Dunn, Andrew K},
  journal={Biomedical Optics Express},
  volume={11},
  number={10},
  pages={5826--5841},
  year={2020},
  publisher={Optical Society of America}
}
```

### Evaluation of resonant scanning as a high-speed imaging technique for two-photon imaging of cortical vasculature

[```https://opg.optica.org/boe/fulltext.cfm?uri=boe-13-3-1374```](https://opg.optica.org/boe/fulltext.cfm?uri=boe-13-3-1374)

bibtex:

```
@article{zhou2022evaluation,
  title={Evaluation of resonant scanning as a high-speed imaging technique for two-photon imaging of cortical vasculature},
  author={Zhou, Annie and Engelmann, Shaun A and Mihelic, Samuel A and Tomar, Alankrit and Hassan, Ahmed M and Dunn, Andrew K},
  journal={Biomedical Optics Express},
  volume={13},
  number={3},
  pages={1374--1385},
  year={2022},
  publisher={Optica Publishing Group}
}
```
