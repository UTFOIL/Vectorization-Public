% pre-processing images to improve quality
% 
% SAM 9/9/21 in collaboration with Annie
% 
%%%%%% input
% original_image_address = 'E:\Annie\pre-processing\og.tif' ; % 8 bit
original_image_address = 'E:\Annie\pre-processing\BL_00001.tif' ; % 16 bit

original_image = double( tif2mat( original_image_address )); % load original image

%%%%%%%%%%%% bakcground subtraction the original image


%%%%%%%%%%!!!!!!!!!!!!!!!! convert to z score up here and divide as well as subtract

sigma_Gaussian_1 = [ 40, 40, 6 ] ; % microns
sigma_Gaussian_2 = [ 100, 100, 50 ];

voxel_lengths = [ 1.37, 1.37, 3 ] ; % for Shaun and Annie's images from year 21

% function [ background_subtracted_image ] = subtract_background( original_image, sigma_Gaussian, voxel_lengths )
% sigma_Gaussian and voxel_lengths are both the same unit of real length (e.g. microns)

pixels_per_sigma_gaussian_1 = sigma_Gaussian_1 ./ voxel_lengths ;
pixels_per_sigma_gaussian_2 = sigma_Gaussian_2 ./ voxel_lengths ;

% blur image, pad 3*sigma distance in all 3 dimensions with pad value equal to the mean image value
 blurred_image_1 = gaussian_blur( original_image, pixels_per_sigma_gaussian_1 ); 

 % background subtraction
pre_processesed_image = original_image - blurred_image_1 ;

mat2tif( pre_processesed_image, 'E:\Annie\pre-processing\backgr_subtr_1.tif' )

 blurred_image_2 = gaussian_blur( pre_processesed_image, pixels_per_sigma_gaussian_2 ); 

 % background subtraction
pre_processesed_image = pre_processesed_image - blurred_image_2 ;

mat2tif( pre_processesed_image, 'E:\Annie\pre-processing\backgr_subtr_2.tif' )


% end % FUNCTION

%%%%%%%%%%%%%% log transfomr the background subtracted image

% subtract minimum
pre_processesed_image = pre_processesed_image - min( pre_processesed_image( : ));

% % normalize by std of background noise 
%     assuming 95 percent of the voxels in the image are bacgkround voxels
% and assuming left half of background distribution is Gaussian shaped. 
% estimating sigma of background noise using the quartile and the mean from the left of the dist.
% CDF @ z-score of -0.67 is 0.25
% CDF @ z-score of  0    is 0.50000
% 0.95 * 0.5  = 0.475 ; z_score =  0    ; z2
% 0.95 * 0.25 = 0.2375; z_score = -0.67 ; z1
% z_score difference = 0.67
intensity_at_median   = quantile( pre_processesed_image( : ), 0.475  );
intensity_at_quartile = quantile( pre_processesed_image( : ), 0.2375 );

sigma_estimate_background_noise = ( intensity_at_median - intensity_at_quartile ) / 0.67 ;

% put background median at 0
pre_processesed_image = ( pre_processesed_image - intensity_at_median ) / sigma_estimate_background_noise ;

% clip any signal below z score of -2
pre_processesed_image( pre_processesed_image < -2 ) = -2 ;

% add 6 for stability (puts  the -6 z_score at 0 intensity), then log-transform
pre_processesed_image = log( pre_processesed_image + 4 );


%%%%%%%% !!!!!!!!!!!!!!!! %%%%%%%%%%%%%%%%% directed diffusion on the log-transformed image



%%%%%%% full-scale contrast stretch for uint8 output
pre_processesed_image = uint8(     256                                                                  ...
                               * (      pre_processesed_image       - min( pre_processesed_image( : ))) ...
                              ./ ( max( pre_processesed_image( : )) - min( pre_processesed_image( : ))));

%%%%%% output
mat2tif( pre_processesed_image, 'E:\Annie\pre-processing\backgr_subtr_log_transf.tif' )

mat2tif( blurred_image_1, 'E:\Annie\pre-processing\blurred_image_1.tif' )
mat2tif( blurred_image_2, 'E:\Annie\pre-processing\blurred_image_2.tif' )
