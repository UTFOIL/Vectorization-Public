function [ chunk_lattice_dimensions, number_of_chunks_in_lattice ] ...
         = get_chunking_lattice_V190( strel_size_in_pixels, max_voxels_per_node, size_of_image )
%% octave downsampling parameters
% SAM 12/17/17
%
% V02 in which the inputs are changed to be ammenable to the vectorize_V02 script 3/5/18 SAM
%
% V10: downsampling removed and octave structure removed SAM 3/28/18
%
% NOTE: this function was renamed to "get_chunking_lattice_V10" to be more descriptive.  Later it should
% be further renamed and broken up into 2 more functions, a new called something like "get_strel",
% and the existing one "get_directory_structure_V10" SAM 3/28/18
%
% NOTE: this chunking could be improved by taking into account the fact that larger scales will have
% larger overlaps input to the GET_STARTS_AND_COUNTS function and should therefore have smaller
% chuels
%
% V190: the first input was simplified. it is no longer a matrix but a vector SAM 8/26/18
%
% Note: a minimum size requirement should be computed and enforced from the largest scale and a 50 %
% efficiency requirement.  SAM 8/26/18
                                


%% chunking parameters
% note: a chuel is a chunking element for binning up the 3-D images across the worker nodes

target_voxel_per_chunk = max_voxels_per_node;

target_chunk_characteristic_length = double( target_voxel_per_chunk ) ^ ( 1 / 3 );

unit_volume_voxel_aspect_ratio =       strel_size_in_pixels               ...
                               / prod( strel_size_in_pixels ) ^ ( 1 / 3 ) ;

target_chunk_dimensions = target_chunk_characteristic_length * unit_volume_voxel_aspect_ratio;

chunk_lattice_dimensions = uint16( max( size_of_image ./ target_chunk_dimensions, 1 ));

number_of_chunks_in_lattice = uint32( prod( chunk_lattice_dimensions, 2 ));

end % function