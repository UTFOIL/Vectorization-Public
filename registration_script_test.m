%% register vector sets demonstration script

load('E:\Michael\batch_190912-113508\vectors\network_190917-232930_tie2gfp16 9juyly2018 870nm region a-082_Cycle00001_Ch3_00000.mat')

strand_numels = cellfun( @( x ) size( x, 1 ), strand_subscripts );

vector_set_before_mat = [cell2mat( strand_subscripts ), cell2mat(  )];

% convert space subscripts to position in microns
vector_set_before_mat( :, 1 : 3 ) = vector_set_before_mat( :, 1 : 3 ) .* microns_per_voxel ;

% convert size subscript to radius in microns
vector_set_before_mat( :, 4 ) = exp( interp1( log( lumen_radius_in_microns_range ), vector_set_before_mat( :, 4 )));

% transform the before with some transformation
test_transformation = [ 10, 5, 7, 1.1, 0.9, 0.95, pi/24, pi/16, pi/32 ]; % [ translation y,x,z, scaling y,x,z, rotation about y,x,z ]

vector_set_before = mat2cell( vector_set_before_mat, strand_numels, 7 );

center_of_mass = calculate_center_of_area( vector_set_before );

vector_set_after_mat = transform_vector_set( vector_set_before_mat - center_of_mass, test_transformation, center_of_mass );

vector_set_after  = mat2cell( vector_set_after_mat,  strand_numels, 7 );

[ vector_before2after, vector_after2before, transformation, goodness_of_registration ] = register_vector_sets( vector_set_before, vector_set_after );