function vector_set_after = transform_vector_set( vector_set_before, transformation, center_of_rotation )
% applies the transform to the vectors

% vectors have already been centerd on their center of rotation

% vectors are in a N x 7 matrix where N is the number of vectors, the 7 columns are pos y,x,z, size,
% direction y,x,z

% positions are microns, and size is radius in microns

% the transformation acts on the before vectors to create the afters. The transormation is a
% rotation (about the center of rotation already subtracted off the vectors) followed by a
% scaling, then returned to the center of rotation, then translated. These moves are represented by
% the transformation vector input [ translation y,x,z, scaling y,x,z, and rotation about y,x,z ]

% rotation


% scaling


% re-center


% translation


vector_set_after = [];



end % FUNCTION