%% Noise Sensitivity Study                                                                          
% The purpose of this study is to investigate the robustness of the vectorization algorithm under
% various amounts of contrast to noise and contrast to background conditions.  A ground truth image
% will be generated from raw data with the help of the curator.  Noise will be added onto this image
% on a pixel to pixel basis, and each pixel intensity will be drawn from a Gaussian distribution
% whose variance is proportional to its mean (approx. Poisson distributed for large mean)).  The
% proportionality constant will be varied on a logarithmic scale.  
%
% Vectorization will be performed algorithmically on many different input images of varying quality
% in order to assess its performance when presented with an input image of a given quality.  Key
% parameter values will be swept to determine the best settings for a given input image.  The energy
% transformation will be performed with many different values for the symmetry_sensitivity
% parameter, swept on a logarithmic scale. The detected vertices will be rendered at the same
% resolution as a ground truth image and compared voxel to voxel to determine their sensitivity and
% specificity.  The set of vertices used in the rendering will be restricted by thresholding the
% energy values. An ROC curve will be produced by sweeping the threshold value. Some of these
% thresholded sets will have edges extracted from them.  They will be chosen by constraining the
% desired specificity of the vertex detection step at a few different values (e.g. 0.999, 0.99,
% 0.9). The extracted edges will be rendered and compared to the ground truth in the same way as
% before, treating each location along the edge as a vertex.  The performances of the different
% parameter value combinations will be summarized by the area under the ROC curve of their Edge
% Extraction step.  Optimal parameters by this metric will be reported for a given input image
% quality.  Performance of the vectorization process as a segmentation strategy will be compared to
% an effective (non-vectorized) intensity thresholding approach.
%
% SAM 11/2/18
% 
% 
% 

%% buid ground truth from raw data using manual vectorization                                       
% 
% clear, close all, vectorize_V200( 'E:\2P imaging\20170802_TxRed_Chronic\Processed_Images_Stack01\PMT01_Red_Images\PMT01_Red_Raw_Data_16bit.tif', 'OutputDirectory', 'C:\Users\sam7343\Box Sync\AA\vectorization\vectorization_V200\' )
% 
% OutputDirectory: C:\Users\sam7343\Box Sync\AA\vectorization\vectorization_V200\\
% Importing settings and data from a previous batch folder in the OutputDirectory.
% Loading previous settings and data from batch_190329-111418 in the OutputDirectory...
% Loading previous settings and data from workflow_190405-125647...
% StartWorkflow selected: energy
% FinalWorkflow selected: network
% Previous energy settings loaded from file at C:\Users\sam7343\Box Sync\AA\vectorization\vectorization_V200\batch_190329-111418\settings\energy_190405-125647
% microns_per_voxel                    = [ 1.07, 1.07, 5 ]
% radius_of_smallest_vessel_in_microns =   1.5
% radius_of_largest_vessel_in_microns  =   8
% approximating_PSF                    =   1
% sample_index_of_refraction           =   1.33
% numerical_aperture                   =   0.95
% excitation_wavelength                =   1.3
% scales_per_octave                    =   3
% max_voxels_per_node_energy           =   500000
% symmetry_ratio_factor                =   0.1 !!! changed to gaussian_to_ideal_ratio with change of units too
% matching_kernel_string               =   3D gaussian
% microns_per_sigma_microscope_PSF     = [ 0.31303, 0.31303, 1.2251 ]
% lumen_radius_in_microns_range        = [ 1.3363, 1.5, 1.6837, 1.8899, 2.1213, 2.3811, 2.6727, 3, 3.3674, 3.7798, 4.2426, 4.7622, 5.3454, 6, 6.7348, 7.5595, 8.4853 ]
% Saving energy workflow settings file: C:\Users\sam7343\Box Sync\AA\vectorization\vectorization_V200\batch_190329-111418\settings\energy_190405-130731
% Previous vertices settings loaded from file at C:\Users\sam7343\Box Sync\AA\vectorization\vectorization_V200\batch_190329-111418\settings\vertices_190405-125647
% space_strel_apothem       =   1
% energy_upper_bound        =   -10
% max_voxels_per_node    =   6000
% length_dilation_ratio       =   2
% exporting_vertex_curation =   0
% Saving vertices workflow settings file: C:\Users\sam7343\Box Sync\AA\vectorization\vectorization_V200\batch_190329-111418\settings\vertices_190405-130731
% Previous edges settings loaded from file at C:\Users\sam7343\Box Sync\AA\vectorization\vectorization_V200\batch_190329-111418\settings\edges_190405-125647
% max_edge_length_per_origin_radius    =   20
% space_strel_apothem_edges            =   1
% number_of_edges_per_vertex           =   10
% edge_walk_temperature                =   0.1
% length_dilation_ratio_edges          =   0.5
% length_dilation_ratio_vertices       =   1
% Saving edges workflow settings file: C:\Users\sam7343\Box Sync\AA\vectorization\vectorization_V200\batch_190329-111418\settings\edges_190405-130731
% Previous network settings loaded from file at C:\Users\sam7343\Box Sync\AA\vectorization\vectorization_V200\batch_190329-111418\settings\network_190405-125647
% sigma_strand_smoothing    =   1
% Saving network workflow settings file: C:\Users\sam7343\Box Sync\AA\vectorization\vectorization_V200\batch_190329-111418\settings\network_190405-130731
% Running Energy workflow for image PMT01_Red_Raw_Data_16bit
% 
% 

%% blur ground truth image to build input images                                   

% shape of signal in simulated image:
signal_shape = 'spherical' ;
% signal_shape = 'annular' ;

cropping = false ;

downsampling = true ;

% is_testing_machine_curation = true ;
is_testing_machine_curation = false ;

% output_directory = 'E:\2P imaging\20170802_TxRed_Chronic\Processed_Images_Stack04\PMT01_Red_Images\' ;
% batch_time_stamp = '190531-164319' ;
% network_time_stamp = '190620-224054' ;

% output_directory = 'E:\Michael\' ;
% batch_time_stamp = '190912-113508' ;
% network_time_stamp = '190917-232930' ;

output_directory = 'E:\2P imaging\20170809_TxRed_Chronic\tiling 180628\' ;
batch_time_stamp = '191030-183626' ;
network_time_stamp = '191113-160636' ;

[ batch_directory, visual_data_directory, visual_vector_directory, vector_directory ] = get_directories( output_directory, batch_time_stamp );


% % !!! load the image name from the settings file instead of hardcoding the string: 
% ROI_name = 'PMT01_Red_Raw_Data_16bit' ;

load([ output_directory, 'batch_', batch_time_stamp, '\settings\batch' ], 'ROI_names' )

ROI_name = ROI_names{ 1 };

if strcmp( signal_shape, 'annular' ), strands_annuli_image = tif2mat([ visual_vector_directory, 'network_', network_time_stamp, '_strands_annuli',            ROI_name, '.tif' ]); end
strands_spheres_image                                      = tif2mat([ visual_vector_directory, 'network_', network_time_stamp, '_strands_spheres_upsampled', ROI_name, '.tif' ]);
strands_bifurcations_image                                 = tif2mat([ visual_vector_directory, 'network_', network_time_stamp, '_strands_bifurcations',      ROI_name, '.tif' ]);

load([ vector_directory, 'network_', network_time_stamp, ROI_name, '.mat' ])

original_network_statistics = network_statistics ;

size_of_image_upsampled     = size( strands_spheres_image      );
size_of_image               = size( strands_bifurcations_image );

% guess the resolution factors as the ratios of the different dimensions of the up and down-sampled
% versions of the image.  Round to make sure it is a recognizable ratio (one that would multiply 60
% evenly)
resolution_factors = round( 60 * ( size_of_image_upsampled - 1 ) ./ ( size_of_image - 1 )) / 60 ;

ground_truth_image_upsampled = logical( strands_spheres_image );

clear( 'strands_spheres_image' )

switch signal_shape

    case 'spherical', original_binary_upsampled_double = double( logical( ground_truth_image_upsampled ));
    case   'annular', original_binary_upsampled_double = double( logical( strands_annuli_image         ));

end

clear( 'strands_annuli_image' )

visual_contrast = 1e4 ;

mat2tif( visual_contrast * original_binary_upsampled_double, [ visual_data_directory, 'original_binary_upsampled.tif' ])

% microns_per_voxel                    = [ 1.07, 1.07, 5 ];
% microns_per_sigma_microscope_PSF     = [ 0.31303, 0.31303, 1.2251 ];

load([ output_directory, 'batch_', batch_time_stamp, '\settings\energy_', batch_time_stamp ], 'microns_per_voxel', 'pixels_per_sigma_PSF' )

microns_per_sigma_microscope_PSF = pixels_per_sigma_PSF .* microns_per_voxel ;

% PSF_fudge_factor                     = 5 ;
PSF_fudge_factor                     = 3 ;
% PSF_fudge_factor                     = 6 ;
% PSF_fudge_factor                     = 1 ;

pixels_per_sigma_gaussian = PSF_fudge_factor * microns_per_sigma_microscope_PSF ./ microns_per_voxel .* resolution_factors ;

max_voxels_per_node = 10 ^ 7 ;

% blurred_image_upsampled = abs( gaussian_blur( original_binary_upsampled_double, pixels_per_sigma_gaussian ));

blurred_image_upsampled = abs( gaussian_blur_in_chunks( original_binary_upsampled_double, pixels_per_sigma_gaussian, max_voxels_per_node ));

clear( 'original_binary_upsampled_double' )

mat2tif( visual_contrast * blurred_image_upsampled, [ visual_data_directory, 'original_blurred_upsampled.tif' ])

if downsampling

    % downsampling back to the original dimensions
    [ mesh_Y, mesh_X, mesh_Z ] = ndgrid( 1 + ( 0 : resolution_factors( 1 ) : resolution_factors( 1 ) * ( size_of_image( 1 ) - 1 )), ... 
                                         1 + ( 0 : resolution_factors( 2 ) : resolution_factors( 2 ) * ( size_of_image( 2 ) - 1 )), ...
                                         1 + ( 0 : resolution_factors( 3 ) : resolution_factors( 3 ) * ( size_of_image( 3 ) - 1 ))  ); 

    blurred_image      =                 interp3( blurred_image_upsampled,                mesh_X, mesh_Y, mesh_Z              );
%     ground_truth_image = logical( round( interp3( double( ground_truth_image_upsampled ), mesh_X, mesh_Y, mesh_Z, 'nearest' )));    
    ground_truth_image = logical( round( interp3( double( ground_truth_image_upsampled ), mesh_X, mesh_Y, mesh_Z            )));    
    
else
    
    blurred_image      =      blurred_image_upsampled ;
    
    ground_truth_image = ground_truth_image_upsampled ;
    
    microns_per_voxel = microns_per_voxel ./ resolution_factors ;
    
    size_of_image = size_of_image_upsampled ;
    
end

clear( 'mesh_X', 'mesh_Y', 'mesh_Z', 'blurred_image_upsampled', 'ground_truth_image_upsampled' );

if cropping
    
    cropping_factors = [ 1, 1, 3 ];
    
    size_of_image = round( size_of_image ./ cropping_factors );
    
    blurred_image       =      blurred_image( 1 : size_of_image( 1 ), 1 : size_of_image( 2 ), 1 : size_of_image( 3 ));
    ground_truth_image  = ground_truth_image( 1 : size_of_image( 1 ), 1 : size_of_image( 2 ), 1 : size_of_image( 3 ));

end % IF cropping

mat2tif( ground_truth_image,                   [ visual_data_directory, 'ground_truth_image.tif' ])
mat2tif( visual_contrast * blurred_image,      [ visual_data_directory,      'blurred_image.tif' ])

is_calculating_intensity_thesholded_statistics = true ;

if is_calculating_intensity_thesholded_statistics

% tic
%     original_area   = area_from_binary_image( ground_truth_image, microns_per_voxel );
%     toc

    original_volume = sum( ground_truth_image( : )) * prod( microns_per_voxel );

end

%% Add contrast, noise

[ batch_directory, visual_data_directory, visual_vector_directory, vector_directory ] = get_directories( output_directory, batch_time_stamp );

if ~exist('blurred_image','var')
    blurred_image = double( tif2mat([ visual_data_directory, 'blurred_image.tif' ])) / visual_contrast ;
end

switch signal_shape
    
    case 'spherical'
        
         contrast_to_noise_range = [ 1; 1.5; 2.5; 4 ]; % 1/13/20, new noise model

        % CNR is contrast to noise is C / sqrt( 2 * B + C ) given by the mean and
        % uncerainty of the operation I(foreground) - I(background) = (B + C) - B. Variance of background is
        % B and variance of foreground is B+C. Expected difference is contrast, C.


        %  background_base =  448 ; 
        %  contrast_base   =  128 ;

         background_base =  1920 ; 
         contrast_base   =  256 ;

        %      noise_base   = ( 2 * background_base + contrast_base ) ^ 0.5 ;
        %       CNR_base   =  contrast_base / noise_base   ; % = 4

%         sweep_direction = 'contrast' ;
        sweep_direction = 'background' ;

        switch sweep_direction

            case 'contrast'

                contrast_range = contrast_to_noise_range .^ 2 / 2 + contrast_to_noise_range .* ( contrast_to_noise_range .^ 2 / 4 + 2 * background_base ) .^ 0.5 ;

                   noise_base   = ( 2 * background_base + contrast_base ) ^ 0.5 ;        

                   noise_range = noise_base ;

            case 'background'

                  contrast_range = contrast_base ;

                background_range = (( contrast_base ./ contrast_to_noise_range ) .^ 2 - contrast_base ) / 2 ;

                     noise_range = ( 2 * background_range + contrast_base ) .^ 0.5 ;

        end

        % %         contrast_range = [ 0.2 ; 0.1 ];
        % %         contrast_range = [ 0.1 ];
        % %         contrast_range = [ 0.05 ];
        % %         contrast_range = [ 0.05 ; 0.1 ; 0.2 ];
        % %         contrast_range = [ 0.05 ; 0.2 ];
        % %         contrast_range = [ 0.025 ];
        % %         contrast_range = 0.2 ; % goes with noise sweep
        % %         contrast_range = ( 0.025 : 0.025 : 0.2 )'; % contrast sweep
        % %         contrast_range = ( 0.025 : 0.06666 : 0.225 )'; % contrast sweep large tiled geom
        %         contrast_range = 0.225 ; % goes with noise sweep large tiled geom
        % 
        % 
        % 
        % %         noise_range    = [ 10  ; 20  ];
        % %         noise_range    = 10 ; % goes with contrast sweep (and for large tiled geom)
        % %         noise_range    = 10 ./ ( 0.375 : 0.125 : 1 )'; % noise sweep
        % %         noise_range = 10 ./ ( 0.375 : 0.5 : 1 )';
        % %         noise_range = 10 ./ ( 0.375 : 0.625 : 1 )';
        % %         noise_range    = 10 ./ linspace( 0.375, 1, 3 )'; % noise sweep large tiled geom
        %         noise_range    = 10 / 0.375 .* ( 0.025 : 0.06666 : 0.225 )' / 0.025 ; % noise sweep large tiled geom V2: 1/10/20
        %         % noise sweep large tiled geom

    case 'annular'

        contrast_range = 1 ;

        noise_range    = [ 0.1 ];
        % noise_range    = [ 0.1; 0.5 ];
        
end

number_of_contrast_indices = numel( contrast_range );
number_of_noise_indices    = numel(    noise_range );

contrast_index_range = 1 : number_of_contrast_indices ;
noise_index_range    = 1 : number_of_noise_indices    ;

input_tif_names   = cell( number_of_contrast_indices, number_of_noise_indices );
input_tif_paths   = cell( number_of_contrast_indices, number_of_noise_indices );

% noisy_tif_names                                                                           ...
%     = cellfun( @( x ) [ 'noisy_', num2stripped( x, '%0.1e' ), '_image.tif' ],             ...
%                mat2cell( noise_range, ones( size( noise_range ))), 'UniformOutput', false );

for contrast_index = contrast_index_range
    
    contrast = contrast_range( contrast_index );
        
    % contrasted_image = contrast * indicator + background
 	contrasted_image = reshape( blurred_image * contrast, [ 1, size_of_image ]) + ( noise_range .^ 2 - contrast ) / 2 ;
    
%     % approximate poisson distribution for large mean input
%     noisy_images = reshape( contrasted_image,                                  [ 1, size_of_image ]) ...
%                  + reshape( contrasted_image .^ 0.5 .* randn( size_of_image ), [ 1, size_of_image ]) ...
%                  .* noise_range                                                                      ;

    % approximate poisson distribution for large mean input
    noisy_images = contrasted_image                                                               ...
                 + contrasted_image .^ 0.5 .* reshape( randn( size_of_image ), [ 1, size_of_image ]);

    % noisy_ground_truth_images ...
    %                      = reshape(                blurred_ground_truth_image,    [ 1, size_of_image ]) ...
    %                      + reshape( poissrnd( abs( blurred_ground_truth_image )), [ 1, size_of_image ]) ...
    %                      .* noise_range                                                                         ;

    for noise_index = noise_index_range

%             contrasted_image = max( blurred_image, 0 ) * contrast + noise_range( noise_index ) ^ 2 ;
% 
%         % exact poisson distribution. Quite slow.  A faster version is available on MATLAB file
%         % exchange, Poissrnd. 
% %         noisy_image = poissrnd( contrasted_image );
%         noisy_image = Poissrnd( contrasted_image );
%         
%         % Poissrnd(0) gets mapped to NaN for some reason
%         noisy_image( isnan( noisy_image )) = 0 ;
            
        input_tif_names{ contrast_index, noise_index } = [ 'contrast_', num2stripped( contrast, '%0.1e' ), '_noise_', num2stripped( noise_range( noise_index ), '%0.1e' ), '.tif' ];

        input_tif_paths{ contrast_index, noise_index } = [ visual_data_directory, input_tif_names{ contrast_index, noise_index }];
        
        mat2tif( squeeze( uint16( noisy_images( noise_index, :, :, : ))), input_tif_paths{ contrast_index, noise_index })
%         mat2tif( noisy_image, input_tif_paths{ contrast_index, noise_index })

    end % FOR noise
end % FOR contrast

clear( 'blurred_image', 'contrasted_image', 'noisy_images' )
%% Run Energy Image calculations (multi-scale derivative filtering)                                 

output_directory_study = 'E:\noise_study_out\' ;
% output_directory_study = 'E:\noise_study_out' ;


switch signal_shape
    
    case 'spherical'
    
        % % k_range     = [ 0, 10 .^ linspace( -1, 1, 3 )]';
        % % k_range     = 10 .^ linspace( -1, 1, 3 );

        % k_range        = [ 0 ; 0.1 ; 1 ];
        % k_range        = [ 0 ; 1 ];
        % k_range        = 0 ;
        % k_range        = [ 0.1 ; 1 ];
        % k_range        = [ 10 ; 100 ];
        % k_range        = 1 ;
        % k_range        = [ 1, 5 ];
        % k_range        = [ 0.5 ; 1 ];
        % k_range         = [ 0.2, 1, 5 ];
        % k_range        = 5 ;
        % k_range        = 2 ;
        % k_range        = [ 1 , 2 ];
        % k_range = [ 0, 0.5 ];
        % k_range        = [ 0, 2 ] ;
        % k_range = Inf ;
        % k_range = [ 2, Inf ];
        % k_range = [ 1, 2, Inf ];
        % k_range = [ 1, Inf ];

        % new def'n of k
%         k_range = [ 0.05, 0.1, 0.2, 0.5, 1 ];
%         k_range = [ 0.3, 0.4, 0.5, 0.6 ]; % .4 .5 amd /6 are the best. .4 was closest on stats. .6 had best voxel class. accuracy
%         k_range = [ .4, .5 ];
%         k_range = 0.5 ;
%         k_range = 1 ;
%         k_range = 0.75 ;
%         k_range = [ 0.4, 0.6, 0.8 ]; % for publication
        k_range = [ 0.6, 0.8, 1 ]; % for publication % 1/7/20
%         k_range = 0.6 ;
                
    case 'annular'
        
        % k_range = 0.5 ;
        % k_range = [ 0, 0.1, 0.5, 1, Inf ];
        % k_range = [ 0.1, 0.5, 1 ];
        
        % k_range = [ 0.25, 0.5, 0.75 ];
        % k_range = [ 0.1, 0.15, 0.2 ];
%         k_range = [ 0.025, 0.05, 0.075 ]; % 0.05 and 0.025 are the best so far for annular
        k_range = [ 0.025, 0.05 ];

end

if is_testing_machine_curation, k_range = [ 0.6, 0.6 ]; end

number_of_k_indices = numel( k_range );

k_index_range = 1 : number_of_k_indices ;

time_stamps_energy = cell( number_of_k_indices, 1 );

matching_kernel_string = '3D gaussian conv annular pulse' ;

switch signal_shape

%     case 'spherical', matching_kernel_string = '3D gaussian conv spherical pulse' ;
    case 'spherical', spherical_to_annular_ratio = 1 ;
    case   'annular', spherical_to_annular_ratio = 0 ;

end

if is_testing_machine_curation, k_index_range = 1 ; end

for k_index = k_index_range

    if isempty( time_stamps_energy{ 1 })

             NewBatch_value     = 'yes'  ;
        PreviousBatch_value     = 'none' ;
     PreviousWorkflow_value     = 'none' ;

    else

             NewBatch_value     = 'no'                    ;
        PreviousBatch_value     = time_stamps_energy{ 1 } ;
     PreviousWorkflow_value     = 'recent'                ;

    end
    
    name_value_pair_inputs = {  'OutputDirectory',                   output_directory_study,            ...
                                'PreviousBatch',                     PreviousBatch_value,               ...
                                'PreviousWorkflow',                  PreviousWorkflow_value,            ...
                                'StartWorkflow',                    'energy',                           ...
                                'FinalWorkflow',                    'one',                              ...
                                'Visual',                           'productive',                       ...
                                'NewBatch',                          NewBatch_value,                    ...
                                'Presumptive',                       true,                              ...
                                'matching_kernel_string',            matching_kernel_string,            ...
                                'gaussian_to_ideal_ratio',                  k_range( k_index ),         ...
                                'spherical_to_annular_ratio',               spherical_to_annular_ratio, ...
                                'microns_per_voxel',                        microns_per_voxel,          ...
                                'radius_of_smallest_vessel_in_microns',     1,                          ...
                                'radius_of_largest_vessel_in_microns',      40,                         ...
                                'approximating_PSF',                        true,                       ...
                                'excitation_wavelength',                    1.3 * PSF_fudge_factor,     ...
                                'scales_per_octave',                        3,                          ... 
                                'max_voxels_per_node_energy',               1000000                     };

    if isempty( time_stamps_energy{ 1 })

        time_stamps_energy{ k_index } = vectorize_V200( input_tif_paths, name_value_pair_inputs{ 1, : }); 
        
    else
        
        time_stamps_energy{ k_index } = vectorize_V200(                  name_value_pair_inputs{ 1, : });            
        
    end
end % FOR k range

if is_testing_machine_curation, k_index_range = 1 : 2 ; time_stamps_energy{ 2 } = time_stamps_energy{ 1 }; end

% the first time stamp is the time stamp for the batch (off which the other directories stem)
[ batch_directory, visual_data_directory, visual_vector_directory, vector_directory ] = get_directories( output_directory_study, time_stamps_energy{ 1 });
%% Extract vertices from the energy image                                                           

time_stamps_vertices = cell( number_of_k_indices, 1 );

for k_index = k_index_range
 
    if is_testing_machine_curation

        % k_index encodes machine learning vs auto curation
        if k_index == 1, VertexCuration_value = 'machine-auto' ; else, VertexCuration_value = 'auto' ; end

    else

        VertexCuration_value = 'auto' ;

    end    
    
    name_value_pair_inputs = {  'OutputDirectory',                   output_directory_study,        ...
                                'PreviousBatch',                     time_stamps_energy{ 1 },       ...
                                'PreviousWorkflow',                  time_stamps_energy{ k_index }, ...
                                'VertexCuration',                    VertexCuration_value,          ... 
                                'StartWorkflow',                    'vertices',                     ...
                                'FinalWorkflow',                    'one',                          ...
                                'Visual',                           'productive',                   ...
                                'NewBatch',                     	'no',                           ...
                                'Presumptive',                       true                           };

    time_stamps_vertices{ k_index } = vectorize_V200( name_value_pair_inputs{ 1, : });
        
end % FOR k range

is_first_time_running_vertices = true ;
%% Compare classification strength of rendered, weighted vertices and of the observed intensity     

number_of_thresholds = 101 ;

areas_under_ROC_curves_energy    = zeros( number_of_contrast_indices, number_of_noise_indices, number_of_k_indices );
areas_under_ROC_curves_intensity = zeros( number_of_contrast_indices, number_of_noise_indices,          1          );

best_intensity_thresholds        = zeros( number_of_contrast_indices, number_of_noise_indices );
max_accuracies_intensity         = zeros( number_of_contrast_indices, number_of_noise_indices );
max_accuracy_index               = zeros( number_of_contrast_indices, number_of_noise_indices );

% sensitivities_intensity          = zeros( number_of_contrast_indices, number_of_noise_indices,           1,         number_of_thresholds );
% specificities_intensity          = zeros( number_of_contrast_indices, number_of_noise_indices,           1,         number_of_thresholds );

sensitivities_energies           = zeros( number_of_k_indices, 1, number_of_thresholds );
specificities_energies           = zeros( number_of_k_indices, 1, number_of_thresholds );
   thresholds_energies           = zeros( number_of_k_indices, 1, number_of_thresholds );
   accuracies_energies           = zeros( number_of_k_indices, 1, number_of_thresholds );
   
 max_accuracies_energy           = zeros( number_of_contrast_indices, number_of_noise_indices, number_of_k_indices, 1 );   
   
% specificity_constraints = [ 0.9, 0.99, 0.999 ];

% specificity_constraints = [ 0.98, 0.99 ];
% specificity_constraints = [ 0.996 ];

specificity_constraints = NaN ; % this is an accuracy constraint now

% specificity_constraints = [ 0.9, 0.99 ];
% specificity_constraints = 0.9 ;
% specificity_constraints = 0.99 ;
% specificity_constraints = [ 0.99, 0.999 ];
% specificity_constraints = 0.999 ;

number_of_vertex_sets = numel( specificity_constraints );

vertex_set_index_range = 1 : number_of_vertex_sets ;

k_index_vertex_set_mesh = ( k_index_range )' + number_of_k_indices * ( vertex_set_index_range - 1 );

number_of_k_indices_by_vertex_sets = number_of_k_indices * number_of_vertex_sets ;

vertex_energy_thresholds = zeros( number_of_contrast_indices, number_of_noise_indices, number_of_k_indices, number_of_vertex_sets );

vertex_names            = cell( number_of_contrast_indices, number_of_noise_indices, number_of_k_indices );
path_to_vertices        = cell( number_of_contrast_indices, number_of_noise_indices, number_of_k_indices );
path_to_visual_vertices = cell( number_of_contrast_indices, number_of_noise_indices, number_of_k_indices );

if ~exist('ground_truth_image','var')
    ground_truth_image = tif2mat( [ output_directory, 'batch_', batch_time_stamp, '\visual_data\ground_truth_image.tif' ]);
end

for contrast_index = contrast_index_range

    for noise_index = noise_index_range
        
        input_image = double( tif2mat( input_tif_paths{ contrast_index, noise_index }));     
        
        [ sensitivities_intensity, specificities_intensity, thresholds_intensities, accuracies_intensities ] ...
                         = get_sens_and_spec( ground_truth_image, - input_image, number_of_thresholds );
                     
                     
        [ max_accuracies_intensity( contrast_index, noise_index ), max_accuracy_index( contrast_index, noise_index )] = max( accuracies_intensities );

        best_intensity_thresholds( contrast_index, noise_index ) = - thresholds_intensities( max_accuracy_index( contrast_index, noise_index ) );                     
        
        clear( 'input_image' );

        for k_index = k_index_range

            vertex_names{ contrast_index, noise_index, k_index } = [ 'curated_vertices_', time_stamps_vertices{ k_index }, '_', input_tif_names{ contrast_index, noise_index }( 1 : end - 4 )];
            
                   path_to_vertices{ contrast_index, noise_index, k_index } = [        vector_directory, vertex_names{ contrast_index, noise_index, k_index }, '.mat' ];              
            path_to_visual_vertices{ contrast_index, noise_index, k_index } = [ visual_vector_directory, vertex_names{ contrast_index, noise_index, k_index }, '.tif' ];  
            
            vertex_energy_image = - double( tif2mat( path_to_visual_vertices{ contrast_index, noise_index, k_index }));         

            [ sensitivities_energies( k_index, 1, : ),                                                 ...
              specificities_energies( k_index, 1, : ),                                                 ...
                 thresholds_energies( k_index, 1, : ),                                                 ...
                 accuracies_energies( k_index, 1, : )  ]                                               ...
                 = get_sens_and_spec( ground_truth_image, vertex_energy_image, number_of_thresholds );
             
        end % FOR k

        clear( 'vertex_energy_image' );
        
        [ max_accuracies_energy( contrast_index, noise_index, :, : ),                               ...
          max_accuracies_indices                                      ]                             ...
                                                                = max( accuracies_energies, [ ], 3 );
        
        % Compute ROC curves and the areas undereath them
        [ areas_under_ROC_curves_energy(    contrast_index, noise_index, : ),                       ...
          areas_under_ROC_curves_intensity( contrast_index, noise_index    ) ]                      ...
                                       = get_ROC( sensitivities_intensity, specificities_intensity, ...
                                                  sensitivities_energies,  specificities_energies,  ...
                                                  contrast_index, noise_index,                      ...
                                             max_accuracies_energy( contrast_index, noise_index, :, : ),   ...
                                                                                    max_accuracies_indices, ...
                                          max_accuracies_intensity( contrast_index, noise_index ),  ...
                                                max_accuracy_index( contrast_index, noise_index ), k_range   );
                                              
%         vertex_energy_thresholds( contrast_index, noise_index, :, : )                                             ...
%             = apply_specificity_constraints( specificities_energies, thresholds_energies, specificity_constraints, 'vertices' );

        vertex_energy_thresholds( contrast_index, noise_index, :, : ) = thresholds_energies( k_index_vertex_set_mesh + number_of_k_indices_by_vertex_sets * ( max_accuracies_indices - 1 ) );                                            

    end % FOR noise
end % FOR contrast
%% Extract edges from the different vertex sets                                                     

% Threshold the vertices with a few different specifities constrained (with the help of the ROC
% curve in the previous section).

vertex_variables_to_edit = { 'vertex_space_subscripts'    , ...
                             'vertex_scale_subscripts'    , ...
                             'vertex_energies'              };
                         
time_stamps_edges = cell( number_of_k_indices, number_of_vertex_sets );

for k_index = k_index_range

    for vertex_set_index = vertex_set_index_range
        
        for contrast_index = contrast_index_range

            for noise_index = noise_index_range

                if is_first_time_running_vertices
                    
                    load( path_to_vertices{ contrast_index, noise_index, k_index }, vertex_variables_to_edit{ 1, : })

                    if vertex_set_index == 1

                        save([ path_to_vertices{ contrast_index, noise_index, k_index }( 1 : end - 4 ), '_backup.mat' ], ...
                               vertex_variables_to_edit{ 1, : }                                                          )

                    end
                else
                    
                    load([ path_to_vertices{ contrast_index, noise_index, k_index }( 1 : end - 4 ), '_backup.mat' ], ...
                           vertex_variables_to_edit{ 1, : }                                                          )

                end
                
                % threshold the vertex objects based on their energy.  Each successive threshold in
                % this vertex set FOR loop is more and more restrictive.  So each vertex set is
                % contained in the previous and may therefore overwrite the previous.
                is_vertex_in_set = vertex_energies <= vertex_energy_thresholds( contrast_index, noise_index, k_index, vertex_set_index );

%                 number_of_vertices_excluded = sum( ~ is_vertex_in_set )

                vertex_energies         =         vertex_energies( is_vertex_in_set    );
                vertex_scale_subscripts = vertex_scale_subscripts( is_vertex_in_set    );
                vertex_space_subscripts = vertex_space_subscripts( is_vertex_in_set, : );

                save( path_to_vertices{ contrast_index, noise_index, k_index }, ...
                      vertex_variables_to_edit{ 1, : }, '-append'               );
                  
            end % FOR noise
        end % FOR contrast

        if is_testing_machine_curation
            
            if k_index == 1, EdgeCuration_value = 'machine-auto' ; else, EdgeCuration_value = 'auto' ; end
            
        else
            
            EdgeCuration_value = 'auto' ;
            
        end
        
        name_value_pair_inputs = {  'OutputDirectory',                   output_directory_study,            ...
                                    'PreviousBatch',                     time_stamps_energy{ 1 },           ...
                                    'PreviousWorkflow',                  time_stamps_vertices{ k_index },   ...
                                    'EdgeCuration',                      EdgeCuration_value,                ...
                                    'StartWorkflow',                    'edges',                            ...
                                    'FinalWorkflow',                    'one',                              ...
                                    'Visual',                           'productive',                       ...
                                    'NewBatch',                     	'no',                               ...
                                    'Presumptive',                       true,                              ...
                                    'max_edge_length_per_origin_radius', 100,                               ...
                                    'space_strel_apothem_edges',         1                                  };

        time_stamps_edges{ k_index, vertex_set_index } = vectorize_V200( name_value_pair_inputs{ 1, : });

    end % FOR vertex set
end % FOR k range

is_first_time_running_vertices = false ;
%% Compare classification strength of rendered, weighted edges and of the observed intensity        

% number_of_thresholds = 101 ;

areas_under_ROC_curves_energy    = zeros( number_of_contrast_indices, number_of_noise_indices, number_of_k_indices, number_of_vertex_sets );
areas_under_ROC_curves_intensity = zeros( number_of_contrast_indices, number_of_noise_indices,                     1                      );

edge_energy_thresholds           = zeros( number_of_contrast_indices, number_of_noise_indices, number_of_k_indices, number_of_vertex_sets );
 max_accuracies_energy           = zeros( number_of_contrast_indices, number_of_noise_indices, number_of_k_indices, number_of_vertex_sets );

% edge_names                       = cell(  number_of_contrast_indices, number_of_noise_indices, number_of_k_indices, number_of_vertex_sets );
path_to_edges                    = cell(  number_of_contrast_indices, number_of_noise_indices, number_of_k_indices, number_of_vertex_sets );
path_to_visual_edges             = cell(  number_of_contrast_indices, number_of_noise_indices, number_of_k_indices, number_of_vertex_sets );

sensitivities_energies           = zeros( number_of_k_indices, number_of_vertex_sets, number_of_thresholds );
specificities_energies           = zeros( number_of_k_indices, number_of_vertex_sets, number_of_thresholds );
   thresholds_energies           = zeros( number_of_k_indices, number_of_vertex_sets, number_of_thresholds );
   accuracies_energies           = zeros( number_of_k_indices, number_of_vertex_sets, number_of_thresholds );
   
if ~exist('ground_truth_image','var')
    ground_truth_image = tif2mat( [ output_directory, 'batch_', batch_time_stamp, '\visual_data\ground_truth_image.tif' ]);
end

for contrast_index = contrast_index_range

    for noise_index = noise_index_range
        
        input_image = double( tif2mat( input_tif_paths{ contrast_index, noise_index }));     
        
        [ sensitivities_intensity, specificities_intensity, thresholds_intensities, accuracies_intensities ] ...
                         = get_sens_and_spec( ground_truth_image, - input_image, number_of_thresholds );
                     
        [ max_accuracies_intensity( contrast_index, noise_index ), max_accuracy_index( contrast_index, noise_index )] = max( accuracies_intensities );
                     
        clear( 'input_image' );

        for k_index = k_index_range

            for vertex_set_index = vertex_set_index_range

%                 edge_name_visual = [ 'curated_edges_', time_stamps_edges{ k_index, vertex_set_index }, '_decomposed_', input_tif_names{ contrast_index, noise_index }( 1 : end - 4 )];
                edge_name_visual = [ 'curated_edges_', time_stamps_edges{ k_index, vertex_set_index }, '_spheres_', input_tif_names{ contrast_index, noise_index }( 1 : end - 4 )];
                edge_name_vector = [ 'curated_edges_', time_stamps_edges{ k_index, vertex_set_index },         '_', input_tif_names{ contrast_index, noise_index }( 1 : end - 4 )];

                       path_to_edges{ contrast_index, noise_index, k_index, vertex_set_index } = [        vector_directory, edge_name_vector, '.mat' ];              
                path_to_visual_edges{ contrast_index, noise_index, k_index, vertex_set_index } = [ visual_vector_directory, edge_name_visual, '.tif' ];  

                edge_energy_image = - double( tif2mat( path_to_visual_edges{ contrast_index, noise_index, k_index, vertex_set_index }));         

                [ sensitivities_energies( k_index, vertex_set_index, : ),                                      ...
                  specificities_energies( k_index, vertex_set_index, : ),                                      ...
                     thresholds_energies( k_index, vertex_set_index, : ),                                      ...
                     accuracies_energies( k_index, vertex_set_index, : )  ]                                    ...
                            = get_sens_and_spec( ground_truth_image, edge_energy_image, number_of_thresholds );

            end % FOR vertex set
        end % FOR k

        clear( 'edge_energy_image' );
        
        [ max_accuracies_energy( contrast_index, noise_index, :, : ),                               ...
          max_accuracies_indices                                      ]                             ...
                                                                = max( accuracies_energies, [ ], 3 );

%         [ max_accuracies_energy( contrast_index, noise_index, :, : ),                               ...
%           max_accuracies_indices                                      ]                             ...
%                                                                 = max( sensitivities_energies .^ 2 + specificities_energies .^ 2, [ ], 3 );

%         [ max_accuracies_energy( contrast_index, noise_index, :, : ),                               ...
%           max_accuracies_indices                                      ]                             ...
%                                                                 = min(( 1 - sensitivities_energies ) .^ 2 + ( 1 - specificities_energies ).^ 2, [ ], 3 );
% 
        edge_energy_thresholds( contrast_index, noise_index, :, : ) = thresholds_energies( k_index_vertex_set_mesh + number_of_k_indices_by_vertex_sets * ( max_accuracies_indices - 1 ) );

%         edge_energy_thresholds( contrast_index, noise_index, :, : )                                                        ...
%             = apply_specificity_constraints( specificities_energies, thresholds_energies, specificity_constraints, 'edges' );
            
        % Compute ROC curves and the areas undereath them
        [ areas_under_ROC_curves_energy(    contrast_index, noise_index, :, : ),                    ...
          areas_under_ROC_curves_intensity( contrast_index, noise_index       ) ]                   ...
                                       = get_ROC( sensitivities_intensity, specificities_intensity, ...
                                                  sensitivities_energies,  specificities_energies,  ...
                                                  contrast_index, noise_index,                      ...
                                             max_accuracies_energy( contrast_index, noise_index, :, : ),   ...
                                                                                    max_accuracies_indices, ...
                                          max_accuracies_intensity( contrast_index, noise_index ),  ...
                                                max_accuracy_index( contrast_index, noise_index ), k_range   );

    end % FOR noise
end % FOR contrast
%% Select most accurate rendering to observe the confusion images and compute network graph         

% Threshold the edges by choosing the most accurate threshold

edge_variables_to_edit = { 'edge_space_subscripts', ...
                           'edge_scale_subscripts', ...
                           'edge_energies'        , ...
                           'edges2vertices'       , ...
                           'mean_edge_energies'     };
                       
time_stamps_network = cell( number_of_k_indices, number_of_vertex_sets );                       
                         
for k_index = k_index_range

    for vertex_set_index = vertex_set_index_range
        
        for contrast_index = contrast_index_range

            for noise_index = noise_index_range

                % threshold the rendered edge image:                
                edge_energy_image = - double( tif2mat( path_to_visual_edges{ contrast_index, noise_index, k_index, vertex_set_index }));
                
                edge_binary_image = edge_energy_image <= edge_energy_thresholds( contrast_index, noise_index, k_index, vertex_set_index );
                
                % compute confusion images
                false_positive_image = ~ ground_truth_image &   edge_binary_image ;
                false_negative_image =   ground_truth_image & ~ edge_binary_image ;
                 true_positive_image =   ground_truth_image &   edge_binary_image ;
                 true_negative_image = ~ ground_truth_image & ~ edge_binary_image ;
                 
                mat2tif( false_positive_image, [ path_to_visual_edges{ contrast_index, noise_index, k_index, vertex_set_index }( 1 : end - 4 ), '_FP.tif' ])
                mat2tif( false_negative_image, [ path_to_visual_edges{ contrast_index, noise_index, k_index, vertex_set_index }( 1 : end - 4 ), '_FN.tif' ])
                mat2tif(  true_positive_image, [ path_to_visual_edges{ contrast_index, noise_index, k_index, vertex_set_index }( 1 : end - 4 ), '_TP.tif' ])
                mat2tif(  true_negative_image, [ path_to_visual_edges{ contrast_index, noise_index, k_index, vertex_set_index }( 1 : end - 4 ), '_TN.tif' ])
                 
                % threshold the edge objects:
                load( path_to_edges{ contrast_index, noise_index, k_index, vertex_set_index }, edge_variables_to_edit{ 1, : })

                save([ path_to_edges{ contrast_index, noise_index, k_index, vertex_set_index }( 1 : end - 4 ), '_backup.mat' ], ...
                       edge_variables_to_edit{ 1, : }                                                         )
                
                % threshold the edge objects based on their energy.
                is_edge_in_set = mean_edge_energies <= edge_energy_thresholds( contrast_index, noise_index, k_index, vertex_set_index );

%                 number_of_edges_excluded = sum( ~ is_edge_in_set )

                edge_energies         =         edge_energies( is_edge_in_set    );
                edge_scale_subscripts = edge_scale_subscripts( is_edge_in_set    );
                edge_space_subscripts = edge_space_subscripts( is_edge_in_set    );
                mean_edge_energies    =    mean_edge_energies( is_edge_in_set    );
                edges2vertices        =        edges2vertices( is_edge_in_set, : );

                save( path_to_edges{ contrast_index, noise_index, k_index, vertex_set_index }, ...
                      edge_variables_to_edit{ 1, : }, '-append'                                 );                 

            end % FOR noise
        end % FOR contrast
        
        clear( 'edge_binary_image', 'edge_energy_image', 'true_negative_image', 'false_negative_image', 'true_positive_image', 'false_positive_image' )
        
        % compute network statistics:
        name_value_pair_inputs = {  'OutputDirectory',                   output_directory_study,            ...
                                    'PreviousBatch',                     time_stamps_energy{ 1 },           ...
                                    'PreviousWorkflow',                  time_stamps_edges{ k_index, vertex_set_index },   ...
                                    'StartWorkflow',                    'network',                          ...
                                    'FinalWorkflow',                    'one',                              ...
                                    'Visual',                           'productive',                       ...
                                    'SpecialOutput',                    'none',                             ...
                                    'NewBatch',                         'no',                               ...
                                    'Presumptive',                       true,                              };

        time_stamps_network{ k_index, vertex_set_index } = vectorize_V200( name_value_pair_inputs{ 1, : });        
        
    end % FOR vertex set
end % FOR k range

clear( 'ground_truth_image' )

%% Compare network statistics across input images and parameters                                    
path_to_network = cell( number_of_contrast_indices, number_of_noise_indices, number_of_k_indices );

clear network_statistics_output

for k_index = k_index_range

    for vertex_set_index = vertex_set_index_range
        
        for contrast_index = contrast_index_range

            for noise_index = noise_index_range

                network_name = [ 'network_', time_stamps_network{ k_index, vertex_set_index }, '_', input_tif_names{ contrast_index, noise_index }( 1 : end - 4 )];
                
            	path_to_network{ contrast_index, noise_index, k_index, vertex_set_index } = [ vector_directory, network_name, '.mat' ];              
                
                load( path_to_network{ contrast_index, noise_index, k_index, vertex_set_index }, 'network_statistics' )
                
                if k_index * vertex_set_index * contrast_index * noise_index == 1
                   
                    % declare struct output to hold network statistics
                    network_statistics_output = network_statistics ;
                    
                else
                
                    network_statistics_output( contrast_index, noise_index, k_index, vertex_set_index ) = network_statistics ;
                
                end % IF first time through loop
            end % FOR noise
        end % FOR contrast
    end % FOR vertex set
end % FOR k range

% contrast_to_noise_range = true_contrast * contrast_range ./ (( 2 * true_background + contrast_range * true_contrast ) .^ 0.5 .* noise_range' );

[ contrast_to_noise_range, indices_sorted_by_CNR ] = sort( contrast_to_noise_range( : ));

number_of_input_images = number_of_contrast_indices * number_of_noise_indices ;

error_in_statistics = zeros( 4, number_of_input_images, number_of_k_indices, number_of_vertex_sets );

error_in_statistics( 1, :, :, : ) = get_network_statistic_error( original_network_statistics, network_statistics_output, 'length'           );
error_in_statistics( 2, :, :, : ) = get_network_statistic_error( original_network_statistics, network_statistics_output, 'area'             );
error_in_statistics( 3, :, :, : ) = get_network_statistic_error( original_network_statistics, network_statistics_output, 'volume'           );
error_in_statistics( 4, :, :, : ) = get_network_statistic_error( original_network_statistics, network_statistics_output, 'num_bifurcations' );

energy_colors = { 'r', 'g', 'b', 'y', 'm', 'c' };

set_style     = { '-', '--', ':', '-.' };

energy_labels = cell( number_of_vertex_sets, number_of_k_indices );

for energy_index = k_index_range
    
    for vertex_set_index = vertex_set_index_range

        if number_of_vertex_sets == 1

            energy_labels{ vertex_set_index, energy_index } = [ 'Gauss to Ideal: ', num2str( k_range( energy_index ), 2 )];

        else

            energy_labels{ vertex_set_index, energy_index } = [ 'Gauss to Ideal: ', num2str( k_range( energy_index ), 2 ), ' set ', num2str( vertex_set_index )];
                            
        end % IF single vertex set
    end % FOR vertex set
end % FOR classifier

if is_calculating_intensity_thesholded_statistics

    CNR_index = 0 ;

%     error_in_area_intensity   = zeros( number_of_contrast_indices * number_of_noise_indices, 1 );
    error_in_volume_intensity = zeros( number_of_contrast_indices * number_of_noise_indices, 1 );

    for noise_index = noise_index_range

        for contrast_index = contrast_index_range

            CNR_index = CNR_index + 1 ;

            input_image = double( tif2mat( input_tif_paths{ contrast_index, noise_index }));        

            binary_input_image = input_image > best_intensity_thresholds( contrast_index, noise_index );

%             area = area_from_binary_image( binary_input_image, microns_per_voxel );

            % summing up the ones in the binary image and then scaling by the voxel volume in microns^3
            volume = sum( binary_input_image( : )) * prod( microns_per_voxel );

%               error_in_area_intensity( indices_sorted_by_CNR( CNR_index )) = ( area   - original_area   ) / original_area   ;
            error_in_volume_intensity( indices_sorted_by_CNR( CNR_index )) = ( volume - original_volume ) / original_volume ;        

        end % FOR contrast
    end % FOR noise
end

clear( 'binary_input_image' )

is_calculating_intensity_thesholded_statistics = false ;

max_accuracies_energy = reshape( max_accuracies_energy, [ number_of_input_images, number_of_k_indices, number_of_vertex_sets ]);

figure

for statistic_index = 1 : 5

    subplot( 5, 1, statistic_index )

    hold on
    
    switch statistic_index
        
        case 1, labels_temp = [{ 'intensity' }; energy_labels( : )]; plot( contrast_to_noise_range, 100 * max_accuracies_intensity( indices_sorted_by_CNR ), 'k-o' ), xticks({}), xticklabels({}), ylabel(   'maximum accuracy [%]' )
        case 2, labels_temp =                   energy_labels( : ) ;                                                                                                  xticks({}), xticklabels({}), ylabel(       'length error [%]' )
        case 3, labels_temp =                   energy_labels( : ) ;                                                                                                  xticks({}), xticklabels({}), ylabel(         'area error [%]' )
        case 4, labels_temp = [{ 'intensity' }; energy_labels( : )]; plot( contrast_to_noise_range, 100 * error_in_volume_intensity,                         'k-o' ), xticks({}), xticklabels({}), ylabel(       'volume error [%]' )
        case 5, labels_temp =                   energy_labels( : ) ;                                                                                                                               ylabel( 'bifurcations error [%]' ), xlabel( 'Contrast / Noise' )

    end % SWITCH statistic
    
    for energy_index = k_index_range

        for vertex_set_index = vertex_set_index_range

            if statistic_index == 1 % max accuracy statistic

                plot( contrast_to_noise_range, 100 * max_accuracies_energy( indices_sorted_by_CNR, energy_index, vertex_set_index ), ...
                      [ energy_colors{ energy_index }, set_style{ vertex_set_index }, 'o' ]                                          )                
                
            else
            
                plot( contrast_to_noise_range, 100 * error_in_statistics( statistic_index - 1, indices_sorted_by_CNR, energy_index, vertex_set_index ), ...
                      [ energy_colors{ energy_index }, set_style{ vertex_set_index }, 'o' ]                                                             )
              
            end
                                  
        end % FOR vertex set
    end % FOR classifier

%     legend( labels_temp )

end % FOR statistic

clear( 'input_image' );

% Get a list of all variables
allvars = whos;

% Identify the variables that ARE NOT graphics handles. This uses a regular
% expression on the class of each variable to check if it's a graphics object
tosave = cellfun(@isempty, regexp({allvars.class}, '^matlab\.(ui|graphics)\.'));

time_stamp = char( datetime('now', 'TimeZone', 'local', 'Format', 'yyMMdd-HHmmss' )); 

% Pass these variable names to save
save([ batch_directory, 'noise_study_save_' time_stamp, '.mat' ], allvars(tosave).name)

%% Functions -------------------------------------------------------------------------------------- 

function [ out ] = num2stripped( x, y )

out = erase( erase( num2str( x, y ), '.' ), '+' );

end

function [ batch_directory, visual_data_directory, visual_vector_directory, vector_directory ] = get_directories( output_directory, time_stamp )

            batch_directory = [ output_directory, 'batch_', time_stamp, '\' ];

      visual_data_directory = [  batch_directory, 'visual_data\'    ];
    visual_vector_directory = [  batch_directory, 'visual_vectors\' ];
           vector_directory = [  batch_directory, 'vectors\'        ];

end

function [ sensitivities, specificities, thresholds, accuracies ] = get_sens_and_spec( ground_truth, classifier, number_of_thresholds )

% !!!! replace '3.5' with network_statistics.volume_density * prod( microns_per_voxel )
volume_frac_of_interest = 3.5 ;

threshold_percentiles = volume_frac_of_interest * logspace( log(( 0.25 / volume_frac_of_interest ) ^ 2 ) / log( 10 ), log(( 100 / volume_frac_of_interest ) ^ 2 ) / log( 10 ), number_of_thresholds ) .^ 0.5 ; 
% threshold_percentiles = 100 * logspace( - 3 / log( 2 ) * log( 10 ), 0 / log( 2 ) * log( 10 ), number_of_thresholds ) .^ 2 ; 
% threshold_percentiles = logspace( - 1 , 2 , number_of_thresholds );
% threshold_percentiles = linspace( 0, 100, number_of_thresholds );

% threshold_percentiles = fliplr( 100 - logspace( -1, 2, number_of_thresholds ));

thresholds = reshape( prctile( classifier( classifier( : ) < max( classifier( : ))), threshold_percentiles ), 1, number_of_thresholds );
% thresholds = reshape( prctile( classifier( : ), threshold_percentiles ), 1, number_of_thresholds );

thresholds(  1  ) = - Inf ;                                                       % Ideal best object (unattainable)
thresholds( end ) =   max( classifier( classifier( : ) < max( classifier( : )))); % worst non-background object

predictions = classifier( : ) <= thresholds ;

 true_positives =   predictions &   ground_truth( : );
% false_positives =   predictions & ~ ground_truth ;
 true_negatives = ~ predictions & ~ ground_truth( : );
% false_negatives = ~ predictions &   ground_truth ;

clear( 'predictions' )
% condition_positives = true_positives + false_negatives ;
% condition_negatives = true_negatives + false_positives ;

condition_positives    = sum( single( ground_truth( : )));

total_number_of_voxels = single( numel( ground_truth ));

clear( 'ground_truth' )

condition_negatives = total_number_of_voxels - condition_positives ;

threshold_indices = 1 : number_of_thresholds ;

accuracies    = zeros( 1, number_of_thresholds );
sensitivities = zeros( 1, number_of_thresholds );
specificities = zeros( 1, number_of_thresholds );

for threshold_index = threshold_indices

    sensitivities( threshold_index ) = sum( single( true_positives( :, threshold_index ))) ./ condition_positives ;
    specificities( threshold_index ) = sum( single( true_negatives( :, threshold_index ))) ./ condition_negatives ;

    accuracies( threshold_index ) = sum( single(   true_positives( :, threshold_index )                             ...
                                                 + true_negatives( :, threshold_index ))) / total_number_of_voxels  ;

end % FOR threshold_indices

end % FUNCTION get_sens_and_spec

function [ areas_under_ROC_curves_energy,                                                           ...
           areas_under_ROC_curves_intensity ]                                                       ...
                                       = get_ROC( sensitivities_intensity, specificities_intensity, ...
                                                  sensitivities_energies,  specificities_energies,  ...
                                                  contrast_index, noise_index,                      ...
                                                  max_accuracies_energy, max_accuracies_indices,    ...
                                                  max_accuracy_intensity, max_accuracy_index, k_range  )

number_of_energy_classifiers = size( specificities_energies, 1 );
number_of_vertex_sets        = size( specificities_energies, 2 );
                                              
figure

axes

hold on

xlabel( '1 - specificity' )
ylabel(     'sensitivity' )

intensity_label = { 'intensity' };

plot( 1 - specificities_intensity, sensitivities_intensity, '-ko', 'MarkerIndices', max_accuracy_index )

if ~isempty( max_accuracy_intensity )

    %display max accuracy at corresponding coordinate on plot
    text( 1 - specificities_intensity( max_accuracy_index ),  ...
              sensitivities_intensity( max_accuracy_index ),  ...
      num2str( max_accuracy_intensity, 4 )                    )

end

title([ 'Contrast #', num2str( contrast_index ), ', Noise #', num2str( noise_index )])

energy_colors = { 'r', 'g', 'b', 'y', 'm', 'c' };

set_style     = { '-', '--', ':', '-.' };

energy_labels = cell( number_of_vertex_sets, number_of_energy_classifiers );

for energy_index = 1 : number_of_energy_classifiers
    
    for vertex_set_index = 1 : number_of_vertex_sets

        if number_of_vertex_sets == 1

            energy_labels{ vertex_set_index, energy_index } = [ 'Gauss to Ideal: ', num2str( k_range( energy_index ), 2 )];

        else

            energy_labels{ vertex_set_index, energy_index } = [ 'Gauss to Ideal: ', num2str( k_range( energy_index ), 2 ), ' set ', num2str( vertex_set_index )];
                            
        end % IF single vertex set

        % display max accuracy text at corresponding coordinate on plot shown with a marker        
        plot( 1 - squeeze( specificities_energies( energy_index, vertex_set_index, : )), ...
                  squeeze( sensitivities_energies( energy_index, vertex_set_index, : )), ...
              [ energy_colors{ energy_index }, set_style{ vertex_set_index }, 'o' ],     ...
              'MarkerIndices', max_accuracies_indices( energy_index, vertex_set_index )  )
            
        text( 1 - specificities_energies( energy_index, vertex_set_index, max_accuracies_indices( energy_index, vertex_set_index )),  ...
                  sensitivities_energies( energy_index, vertex_set_index, max_accuracies_indices( energy_index, vertex_set_index )),  ...
          num2str( max_accuracies_energy( 1, 1, energy_index, vertex_set_index ), 4 )                                                 )

    end % FOR vertex set
end % FOR classifier

xlim([ 0, 0.05 ]) % only consider above 95 % specificity

legend([ intensity_label; energy_labels( : )], 'Location', 'southeast' )

% calculate the area under the receiver operating curve with the trapezoid method
areas_under_ROC_curves_intensity                                                                                                             ...
                     = squeeze( sum(    (     sensitivities_intensity(      2 : end )     + sensitivities_intensity(      1 : end - 1 )) / 2 ...
                                     .* ( 1 - specificities_intensity(      2 : end ) - 1 + specificities_intensity(      1 : end - 1 ))     ));

areas_under_ROC_curves_energy                                                                                                                ...
                     = squeeze( sum(    (     sensitivities_energies( :, :, 2 : end )     + sensitivities_energies( :, :, 1 : end - 1 )) / 2 ...
                                     .* ( 1 - specificities_energies( :, :, 2 : end ) - 1 + specificities_energies( :, :, 1 : end - 1 )),    ...
                                     3                                                                                                       ));

end % FUNCTION get_ROC

function [ energy_thresholds, energy_threshold_indices ] = apply_specificity_constraints( specificities, thresholds, specificity_constraints, method )

number_of_classifiers = size( specificities, 1 );

number_of_vertex_sets = numel( specificity_constraints );

energy_thresholds        = zeros( number_of_classifiers, number_of_vertex_sets );

energy_threshold_indices = zeros( number_of_classifiers, number_of_vertex_sets );

for classifier_index = 1 : number_of_classifiers

    for vertex_set_index = 1 : number_of_vertex_sets
        
        switch method
            
            case 'vertices'
                
                energy_threshold_indices( classifier_index, vertex_set_index ) = find( specificities( classifier_index,                   : ) >= specificity_constraints( vertex_set_index ), 1, 'last' );
        
            case 'edges'
                
                energy_threshold_indices( classifier_index, vertex_set_index ) = find( specificities( classifier_index, vertex_set_index, : ) >= specificity_constraints( vertex_set_index ), 1, 'last' );

        end
        
        energy_thresholds( classifier_index, vertex_set_index ) = thresholds( classifier_index, 1, energy_threshold_indices( classifier_index, vertex_set_index ) );
                
    end % FOR vertex set
end % FOR classifier
end % FUNCTION apply_specificity_constraints

function error = get_network_statistic_error( original_network_statistics, network_statistics_output, field )

% extracting the entries from each field, comparing to the original, catching these in a vector,
% and then reshaping to match the dimensions of the input, except with the first and second
% dimensions combined into one.
error = reshape( eval([ '([ network_statistics_output.', field, ' ] - original_network_statistics.', field, ' ) ./ original_network_statistics.', field, ' ;' ]), ...
                   size( network_statistics_output, 1 )  ...
                 * size( network_statistics_output, 2 ), ...
                   size( network_statistics_output, 3 ), ...
                   size( network_statistics_output, 4 )  );

end % get_network_statistic_error

function area = area_from_binary_image( binary_input_image, microns_per_voxel )

figure_handle = figure ;

fv = isosurface( double( binary_input_image ), 0.5 );

% run 3D rendering and patch area analysis
p = patch( fv, 'visible', 'off' );
%         isonormals(x,y,z,v,p)
%         set(p,'FaceColor','red','EdgeColor','none');
%         daspect([1 1 1])
%         view(3); 
%         camlight 
%         lighting gouraud
verts = get(p, 'Vertices');

% scaling verts from units of voxels to microns
verts = verts .* microns_per_voxel ;

faces = get(p, 'Faces');
a = verts(faces(:, 2), :) - verts(faces(:, 1), :);
b = verts(faces(:, 3), :) - verts(faces(:, 1), :);
c = cross(a, b, 2);
area = 1/2 * sum(sqrt(sum(c.^2, 2)));
%         fprintf('\nThe surface area is %f\n\n', area)

close( figure_handle )

end % FUNCTION area_from_binary_image

%% Scratch

% get_ROC 
% expecting 1 x n row vectors for the ground_truth and intensity inputs, and expecting m x n for the
% energy.  Each row of the energy is treated as a new classification method to have its own ROC
% curve and summary area.  Expecting lower energy to correspond to a higher probability of true in
% the ground_truth vector.
%
% SAM 4/9/19

% function [ areas_under_ROC_curves_energy, areas_under_ROC_curves_intensity, energy_thresholds ] = get_ROC( ground_truth, energy, intensity, specificity_constraints )
% 
% number_of_thresholds = 31 ;
% 
% energy = [ - intensity ; ...
%               energy     ];
% 
% number_of_classifiers = size( energy, 1 );       
%        
% threshold_percentiles = logspace( -2, 2, number_of_thresholds );
% 
% thresholds = reshape( prctile( energy, threshold_percentiles, 2 ), number_of_classifiers, 1, number_of_thresholds );
% 
% % thresholds( :, 1 ) = thresholds( :, 1 ) - 1 ;
% thresholds( :, 1 ) = - Inf ;
% 
% predictions = energy <= thresholds ;
% 
%  true_positives =   predictions &   ground_truth ;
% % false_positives =   predictions & ~ ground_truth ;
%  true_negatives = ~ predictions & ~ ground_truth ;
% % false_negatives = ~ predictions &   ground_truth ;
% 
% % condition_positives = true_positives + false_negatives ;
% % condition_negatives = true_negatives + false_positives ;
% 
% condition_positives = sum(   ground_truth, 2 );
% % condition_negatives = sum( ~ ground_truth, 2 );
% 
% condition_negatives = numel( ground_truth ) - condition_positives ;
% 
% sensitivities = squeeze( sum( true_positives, 2 )) ./ condition_positives ;
% specificities = squeeze( sum( true_negatives, 2 )) ./ condition_negatives ;
% 
% % excluding the intensity classifier.  Computing an energy threshold for each specificity constraint
% energy_thresholds = zeros( number_of_classifiers - 1, number_of_vertex_sets );
% 
% number_of_vertex_sets = numel( specificity_constraints );
% 
% for classifier_index = 2 : number_of_classifiers
% 
%     for vertex_set_index = vertex_set_index_range
%         
%         last_good_specificity_index = find( specificities( classifier_index, : ) >= specificity_constraints( vertex_set_index ), 1, 'last' );
%         
%         energy_thresholds( classifier_index - 1, vertex_set_index ) = thresholds( classifier_index, last_good_specificity_index );
%         
%     end % FOR vertex set
% end % FOR classifier
% 
% figure
% 
% axes
% 
% hold on
% 
% xlabel( '1 - specificity' )
% ylabel(     'sensitivity' )
% 
% labels = cell( number_of_classifiers, 1 );
% 
% for classifier_index = 1 : number_of_classifiers
%     
%     labels{ classifier_index } = [ 'energy ', num2str( classifier_index - 1 )];
%     
%     plot( 1 - specificities( classifier_index, : ), sensitivities( classifier_index, : ))
%     
%     % % label different threshold values
%     % text(.2,.1,'example minimum ratings for positive diagnosis shown')
% 
% %     number_of_labels = 4;
% % 
% %     for k=0:number_of_labels-1
% % 
% %         threshold_index = 1 + floor(( number_of_thresholds - 1 ) * k / ( number_of_labels - 1 ));
% % 
% %         %display actual threshold at corresponding coordinate on plot
% %         text( 1 - specificities( threshold_index ),       ...
% %                   sensitivities( threshold_index ) + .02, ...
% %             num2str( thresholds( threshold_index ), 2 )   )
% % 
% %     end
% 
% end % FOR classifier
% 
% labels{ 1 } = 'intensity' ;
% 
% legend( labels )
% 
% % calculate the area under the receiver operating curve with the trapezoid method
% areas_under_ROC_curves = sum(    (     sensitivities( :, 2 : end )     + sensitivities( :, 1 : end - 1 )) / 2 ...
%                               .* ( 1 - specificities( :, 2 : end ) - 1 + specificities( :, 1 : end - 1 )),    ...
%                               2                                                                               );
%                           
% areas_under_ROC_curves_energy    = areas_under_ROC_curves( 2 : end );
% areas_under_ROC_curves_intensity = areas_under_ROC_curves(    1    );
% 
% end % FUNCTION get_ROC