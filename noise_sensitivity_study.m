%% Noise Sensitivity Study
% The purpose of this study is to investigate the robustness of the vectorization algorithm under
% various amounts of signal to noise.  Noise will be added onto the original on a pixel to pixel
% basis, and each pixel intensity will be drawn from a Gaussian distribution.  The width of the
% Gaussian compared to the contrast between foreground and background will be varied on a
% logarithmic scale.  Vectorization will be attempted with many different values for two inputs, the
% symmetry_sensitivity parameter and the temperature of the edge trajectory random walks, each one
% swept on a logarithmic scale.  The best vectorization classifier will be recorded.  Classifier
% strength will be judged by the area under its ROC curve with respect to the input ground truth
% image.
%
% SAM 11/2/18

%% buid ground truth from raw data using manual vectorization
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
% symmetry_ratio_factor                =   0.1
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

%% blur, add noise, and automatically vectorize ground truth image
time_stamp = '190415-174607' ;

output_directory = 'C:\Users\sam7343\Box Sync\AA\vectorization\vectorization_V200\' ;

[ batch_directory, visual_data_directory, visual_vector_directory ] = get_directories( output_directory, time_stamp );

strands_spheres_image = tif2mat([ visual_vector_directory, 'network_', time_stamp, '_strands_spheres_PMT01_Red_Raw_Data_16bit.tif' ]);

size_of_image_upsampled = size( strands_spheres_image );

% downsampling back to the original dimensions
strands_bifurcations_image = tif2mat([ visual_vector_directory, 'network_', time_stamp, '_strands_bifurcations_PMT01_Red_Raw_Data_16bit.tif' ]);

size_of_image = size( strands_bifurcations_image );

% resolution_factors = [ 1.5, 1.5, 5 ];

% resolution_factors = [ 2, 2, 10 ];

resolution_factors = round( 60 * ( size_of_image_upsampled - 1 ) ./ ( size_of_image - 1 )) / 60 ;

ground_truth_image_upsampled = logical( strands_spheres_image );

signal = 1000 ;

background = 10 ;

ground_truth_image_double = signal * double( ground_truth_image_upsampled ) + background ;

mat2tif( ground_truth_image_double, [ visual_data_directory, 'original_ground_truth_binary_', num2str( signal ), '_and_', num2str( background ), '.tif' ])

% microns_per_voxel                    = [ 1.07, 1.07, 5 ]
% microns_per_sigma_microscope_PSF     = [ 0.31303, 0.31303, 1.2251 ]
% PSF_fudge_factor = 3 ;
blurred_ground_truth_image_upsampled = gaussian_blur( ground_truth_image_double, 3 * [ 0.31303, 0.31303, 1.2251 ] ./ [ 1.07, 1.07, 5 ] .* resolution_factors );

mat2tif( blurred_ground_truth_image_upsampled, [ visual_data_directory, 'blurred_ground_truth_image_upsampled.tif' ])

% downsampling back to the original dimensions
[ mesh_Y, mesh_X, mesh_Z ] = ndgrid( 1 + ( 0 : resolution_factors( 1 ) : resolution_factors( 1 ) * ( size_of_image( 1 ) - 1 )), ... 
                                     1 + ( 0 : resolution_factors( 2 ) : resolution_factors( 2 ) * ( size_of_image( 2 ) - 1 )), ...
                                     1 + ( 0 : resolution_factors( 3 ) : resolution_factors( 3 ) * ( size_of_image( 3 ) - 1 ))  ); 

blurred_ground_truth_image = interp3( blurred_ground_truth_image_upsampled, mesh_X, mesh_Y, mesh_Z );

mat2tif( blurred_ground_truth_image, [ visual_data_directory, 'blurred_ground_truth_image.tif' ])

ground_truth_image = logical( round( interp3( double( ground_truth_image_upsampled ), mesh_X, mesh_Y, mesh_Z, 'nearest' )));

mat2tif( ground_truth_image, [ visual_data_directory, 'ground_truth_image.tif' ])

% noise_range = 10 .^ linspace( -1, 3, 5 );
% k_range     = 10 .^ linspace( -1, 1, 3 );

% noise_range = [ 0, 10 .^ linspace(  0, 3, 4 )]';
% k_range     = [ 0, 10 .^ linspace( -1, 1, 3 )]';

noise_range = [ 0; 5 ];
k_range     = [ 0; 1 ];

number_of_noise_indices = numel( noise_range );
number_of_k_indices     = numel(     k_range );

noise_index_range  = 1 : number_of_noise_indices ;
k_index_range      = 1 : numel(     k_range );

time_stamps = cell( numel( k_range ), 1 );

noisy_tif_names                                                                                     ...
    = cellfun( @( x ) [ 'noisy_', erase( erase( num2str( x, '%0.1e' ),'.' ), '+' ), '_image.tif' ], ...
               mat2cell( noise_range, ones( size( noise_range ))), 'UniformOutput', false           );

noisy_tif_paths = cellfun( @( x ) [ visual_data_directory, x ], noisy_tif_names, 'UniformOutput', false );

% size_of_input_image_set = [ number_of_noise_indices, size_of_image ];
% 
% noisy_ground_truth_images = zeros( size_of_input_image_set );

noisy_ground_truth_images ...
                          = reshape(      blurred_ground_truth_image,                                    [ 1, size_of_image ]) ...
                          + reshape( abs( blurred_ground_truth_image ) .^ 0.5 .* randn( size_of_image ), [ 1, size_of_image ]) ...
                         .* noise_range                                                                                        ;

% noisy_ground_truth_images ...
%                      = reshape(                blurred_ground_truth_image,    [ 1, size_of_image ]) ...
%                      + reshape( poissrnd( abs( blurred_ground_truth_image )), [ 1, size_of_image ]) ...
%                      .* noise_range                                                                         ;

for noise_index = noise_index_range
                                                  
    mat2tif( squeeze( noisy_ground_truth_images( noise_index, :, :, : )), noisy_tif_paths{ noise_index })
           
end % FOR NOISE

for k_index = k_index_range

    if isempty( time_stamps{ 1 })

        NewBatch_value = 'yes' ;

        PreviousBatch_value    = 'none' ;
        PreviousWorkflow_value = 'none' ;

    else

        NewBatch_value = 'no'  ;

        PreviousBatch_value    = 'recent' ;
        PreviousWorkflow_value = 'recent' ;

    end
    
    name_value_pair_inputs = {  'OutputDirectory',           output_directory,                  ...
                                'PreviousBatch',            PreviousBatch_value,                ...
                                'PreviousWorkflow',      PreviousWorkflow_value,                ...
                                'VertexCuration',                   'none',                     ...
                                'EdgeCuration',                     'none',                     ...
                                'StartWorkflow',                    'energy',                   ...
                                'FinalWorkflow',                    'vertices',                 ...
                                'Visual',                           'productive',               ...
                                'NewBatch',                     NewBatch_value,                 ...
                                'Presumptive',                              true,               ...
                                'symmetry_ratio_factor',                    k_range( k_index ), ...
                                'microns_per_voxel',                        [ 1.07, 1.07, 5 ],  ...
                                'radius_of_smallest_vessel_in_microns',     1.5,                ... % lower this.
                                'radius_of_largest_vessel_in_microns',      8,                  ... % raise this.
                                'approximating_PSF',                        true,               ...
                                'scales_per_octave',                        1,                  ... % double this.
                                'max_voxels_per_node_energy',               500000              };    

    if isempty( time_stamps{ 1 })

        time_stamps{ k_index } = vectorize_V200(  noisy_tif_paths, name_value_pair_inputs{ 1, : }); 
        
    else
        
        time_stamps{ k_index } = vectorize_V200(                   name_value_pair_inputs{ 1, : });            
        
    end
end % FOR k RANGE

[ batch_directory, visual_data_directory, visual_vector_directory ] = get_directories( output_directory, time_stamps{ 1 });

size_of_output_image_set = [ number_of_noise_indices, number_of_k_indices, size_of_image ];

vertex_energy_images = zeros( size_of_output_image_set, 'uint16' );

areas_under_ROC_curves_energy    = zeros( number_of_noise_indices, number_of_k_indices );
areas_under_ROC_curves_intensity = zeros( number_of_noise_indices,          1          );

%% compare classification strength of different energy functions and of the observed intensity
for noise_index = noise_index_range

    for k_index = k_index_range

%         for f_index = f_index_range
            
%             vertex_energy_images( noise_index, :, :, : )                                            ...
%                 = tif2mat([ visual_vector_directory, 'vertices_', time_stamps{ k_index, f_index }, '_', noisy_tif_names{ noise_index }]);
            
            vertex_energy_images( noise_index, k_index, :, :, : )                                                            ...
                = tif2mat([ visual_vector_directory, 'vertices_', time_stamps{ k_index }, '_', noisy_tif_names{ noise_index }]);            

%         end % FOR f range
    end % FOR k RANGE
    
    % Compute ROC curves and the areas undereath them
    [ areas_under_ROC_curves_energy(    noise_index, : ),                          ...
      areas_under_ROC_curves_intensity( noise_index    )  ]                        ...
        = get_ROC(                      ground_truth_image( : )',                  ...
                   - double( squeeze( vertex_energy_images( noise_index, :, : ))), ...
                                 noisy_ground_truth_images( noise_index, : )       );
    
end % FOR NOISE

function [ batch_directory, visual_data_directory, visual_vector_directory ] = get_directories( output_directory, time_stamp )

            batch_directory = [ output_directory, 'batch_', time_stamp, '\' ];

      visual_data_directory = [ batch_directory, 'visual_data\'    ];
    visual_vector_directory = [ batch_directory, 'visual_vectors\' ];

end

%% get_ROC 
% expecting 1 x n row vectors for the ground_truth and intensity inputs, and expecting m x n for the
% energy.  Each row of the energy is treated as a new classification method to have its own ROC
% curve and summary area.  Expecting lower energy to correspond to a higher probability of true in
% the ground_truth vector.
%
% SAM 4/9/19

function [ areas_under_ROC_curves_energy, areas_under_ROC_curves_intensity ] = get_ROC( ground_truth, energy, intensity )

number_of_thresholds = 31 ;

energy = [ - intensity;  ...
              energy    ];

number_of_classifiers = size( energy, 1 );       
       
threshold_percentiles = logspace( -2, 2, number_of_thresholds );

thresholds = reshape( prctile( energy, threshold_percentiles, 2 ), number_of_classifiers, 1, number_of_thresholds );

% thresholds( :, 1 ) = thresholds( :, 1 ) - 1 ;
thresholds( :, 1 ) = - Inf ;

predictions = energy <= thresholds ;

 true_positives =   predictions &   ground_truth ;
% false_positives =   predictions & ~ ground_truth ;
 true_negatives = ~ predictions & ~ ground_truth ;
% false_negatives = ~ predictions &   ground_truth ;

% condition_positives = true_positives + false_negatives ;
% condition_negatives = true_negatives + false_positives ;

condition_positives = sum(   ground_truth, 2 );
% condition_negatives = sum( ~ ground_truth, 2 );

condition_negatives = numel( ground_truth ) - condition_positives ;

sensitivities = squeeze( sum( true_positives, 2 )) ./ condition_positives ;
specificities = squeeze( sum( true_negatives, 2 )) ./ condition_negatives ;

figure

axes

hold on

xlabel( '1 - specificity' )
ylabel(     'sensitivity' )

labels = cell( number_of_classifiers, 1 );

for classifier_index = 1 : number_of_classifiers
    
    labels{ classifier_index } = [ 'energy ', num2str( classifier_index - 1 )];
    
    plot( 1 - specificities( classifier_index, : ), sensitivities( classifier_index, : ))
    
    % % label different threshold values
    % text(.2,.1,'example minimum ratings for positive diagnosis shown')

%     number_of_labels = 4;
% 
%     for k=0:number_of_labels-1
% 
%         threshold_index = 1 + floor(( number_of_thresholds - 1 ) * k / ( number_of_labels - 1 ));
% 
%         %display actual threshold at corresponding coordinate on plot
%         text( 1 - specificities( threshold_index ),       ...
%                   sensitivities( threshold_index ) + .02, ...
%             num2str( thresholds( threshold_index ), 2 )   )
% 
%     end

end

labels{ 1 } = 'intensity' ;

legend( labels )

% calculate the area under the receiver operating curve with the trapezoid method
areas_under_ROC_curves = sum(    (     sensitivities( :, 2 : end )     + sensitivities( :, 1 : end - 1 )) / 2 ...
                              .* ( 1 - specificities( :, 2 : end ) - 1 + specificities( :, 1 : end - 1 )),    ...
                              2                                                                               );
                          
areas_under_ROC_curves_energy    = areas_under_ROC_curves( 2 : end );
areas_under_ROC_curves_intensity = areas_under_ROC_curves(    1    );

end
