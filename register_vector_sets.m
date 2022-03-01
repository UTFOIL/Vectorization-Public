function [ vector_before2after, vector_after2before, transformation, goodness_of_registration ] = register_vector_sets( vector_set_before, vector_set_after )
%REGISTER_VECTOR_SETS takes a before and after set of vectors and attempts to find the best
%registration transformation to take the set of vectors from the before to the after. The
%transformation acts on the before vectors to put them in the same coordinate space as the after
%vectors. The mappings (ideally invertible, but not in general) between the two vector sets are
%extracted.

% inputs are arranged in cell arrays where each entery is strand object. Strand objects are 1D
% discrete trajectories represented by 7 column (3 space, 1 size, 3 direction) matrices where each
% consecutive row is a consecutive, discrete location along the discrete, 1D object. Units for pos
% are microns. units for size are microns of the radius. directions are unit vectors in real,
% euclidean 3-space.

%% center of rotation
% find the center of surface area of the before vector set in order to substract this off before
% rotation (added back in after) so that the rotation is less likely to be coupled to the
% translation. Subtract it off here once for all. Add it back in after each rotation in the
% optimization loop.

% microns_per_voxel = [ 1, 1, 1 ] ; % !!! put strand subscripts in units of microns before this function. Also put the size subscript in microns of the radius

strand_subscripts_after = cell2mat( vector_set_after );

strand_space_subscripts_after = strand_subscripts_after( :, 1 : 3 );

% strand_subscripts_before = cell2mat( vector_set_before );
% 
% strand_space_subscripts_before = strand_subscripts_before( :, 1 : 3 );

center_of_area = calculate_center_of_area( vector_set_before );

%% subsampling
% resample the strands from the before and after to be sampled N times: once at each endpoints and N
% - 2 points in between (evenly spaced).

samples_per_strand = 3 ;

% use_case = 'registration' ;

vector_set_before_subsampled = subsample_vectors( vector_set_before, samples_per_strand );
vector_set_after_subsampled  = subsample_vectors( vector_set_after , samples_per_strand );

% subtract off center of area from the before set
vector_set_before_subsampled = vector_set_before_subsampled - center_of_area ;

%% optimization loop

is_first_iteration = true ;

is_optimizing = true ;

number_of_parameters = 9 ;

parameter_mean_vector_observed  = zeros( number_of_parameters, 1 ); 
parameter_sigma_vector_observed = zeros( number_of_parameters, 1 ); 
parameter_abs_change            = zeros( number_of_parameters, 1 ); 
transformation                  = zeros( number_of_parameters, 1 ); 

% initial standard deviations for the variations to be explored in each registration parameter.
% Initial gueses for the soft maxes are assumed to be overestimates for each variation divided by
% six (so that 99% of guesses are inside the overestimate)


% overestimate the translation parameter by looking at the side lengths of a box containing the
% extreme positions in x, y, and z. Because the translation in x can only go as far as the target
% image is long in x (rotation done first, so translation is in target domain)
max_size_of_image_after  = max( vector_set_after_subsampled(  :, 1 : 3 )) - min( vector_set_after_subsampled(  :, 1 : 3 ));
max_size_of_image_before = max( vector_set_before_subsampled( :, 1 : 3 )) - min( vector_set_before_subsampled( :, 1 : 3 ));

% max_scaling = max( max( max_size_of_image_before ) / min( max_size_of_image_after ), max( max_size_of_image_after ) / min( max_size_of_image_before ));
max_scaling = log( 2 ); % two fold scaling is soft max

max_rotation = pi / 2 ;

soft_maxes = [ max_size_of_image_after / 3, abs( max_scaling ) / 3 * [ 1, 1, 1 ], max_rotation / 3 * [ 1, 1, 1 ]]; % translation in y, x, z, scaling in y, x, z, rotation about y, x, z

number_of_perturbations_in_set = 1000 ;

while is_optimizing

    %% perturbation
    % add perturbations to each parameter gaussian distributed in rotation angle, translation, and log
    % of scaling factor. perturbations are zero-mean gaussian distributed with standard deviation given by the
    % original soft bounds, or the larger of (A) the change in the parameter from the last iteration and
    % (B) the variablea objective function-weighted variance of the parameter in the last set of
    % perturbations. The centralvalues for each functions are either the original values or the
    % objective function-weighted mean of the parameter in the last set of perturbations.

    parameter_sigma_vector = soft_maxes * is_first_iteration + parameter_sigma_vector_observed * ~ is_first_iteration ;
    
    parameter_mean_vector  = parameter_mean_vector_observed ; % original guess is all zeros, the original orientation

    is_first_iteration = false ;

    perturbation_vector_set = parameter_sigma_vector .* randn( number_of_parameters, number_of_perturbations_in_set );
    
    perturbation_index = 0 ;
    
    registration_scores = zeros( 1, number_of_perturbations_in_set );
    
    for perturbation_vector = perturbation_vector_set

        perturbation_index = perturbation_index + 1 ;
        
        %% transformation
        % transform the before vectors with the (perturbed) current best guess.
        vector_set_before_transformed_mat = transform_vector_set( vector_set_before_mat - center_of_mass, transformation + perturbation_vector, center_of_mass );

        %% evaluation
        % loop through the vectors in the before looking for the best match in the after and vice
        % versa
        registration_scores( perturbation_index ) = evaluate_registration( vector_set_before_transformed_mat, vector_set_after_mat, parameter_sigma_vector );

    end % FOR perturbation
    
    %% updating perturbation statistics

    parameter_mean_vector_observed_previous = parameter_mean_vector_observed ;

    parameter_mean_perturbation_observed = sum( registration_scores .* perturbation_vector_set, 2 ) / sum( registration_scores );
    
    parameter_mean_vector_observed = parameter_mean_vector_observed_previous + parameter_mean_perturbation_observed ;

	parameter_variance_vector_observed = sum( registration_scores .* ( perturbation_vector_set - parameter_mean_perturbation_observed ) .^ 2, 2 ) / sum( registration_scores );;

    parameter_sigma_vector_observed = parameter_variance_vector_observed .^ 0.5 ;
    
    parameter_abs_change = abs( parameter_mean_vector_observed_previous - parameter_mean_vector_observed );
       
    parameter_sigma_vector_observed = max( parameter_sigma_vector_observed, parameter_abs_change );

% %     argmax of the registration score for best perturbation
%     [ registration_score, best_perturbation_idx ] = max( registration_scores );
%     
%     transformation = transformation + perturbation_vector_set( :, best_perturbation_idx );
    
    % registration_score-weighted mean parameter vector observed most recently for best perturbation
    transformation = transformation + parameter_mean_vector_observed ;

end % WHILE optimizing

end % FUNCTION main

