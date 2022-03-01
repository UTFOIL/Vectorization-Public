function [ registration_score, vector_before2after, vector_after2before ] = evaluate_registration( vector_set_after_mat, vector_set_before_mat, parameter_sigma_vector )

vector_set_after_mat  = vector_set_after_mat'  ;    
vector_set_before_mat = vector_set_before_mat' ;

% !!!! parallelize one of these, test which one is faster, inside has less overhead but gets called
% more.
for vector_A = vector_set_after_mat

    for vector_B = vector_set_before_mat

        % decompose the rotation into rotations around x, y, and z in that order, so that we can use
        % the independent sigma's for each one: First perform the inverse (translation and) scaling
        % and rotation on both vectors to put them in the original coordinate space.
        
        % Now find angle between them in the yz plane to solve for the residual angle for the
        % rotation about the x axis
        
        % Now apply the original (not the corrected) rotation about x to both vectors and look for
        % angle between them in the xz plane for the residual angle for the roation about y.
        
        % Now apply the original rotation about y to both of those rotated vectors and look for
        % angle between them in the xy plane to solve for residual angle for rotation about z
        
        registration_matrix( vector_B, vector_A ) = [];
        
        
        
    end % FOR vector before
end % FOR vector after

% pick max's for vector in both directions of the evaluation
[ best_match_score_before2after, vector_before2after ] = max( registration_matrix, [ ], 2 );
[ best_match_score_after2before, vector_after2before ] = max( registration_matrix, [ ], 1 );

% sum up all the max's both ways individually. Product these values (want forward AND backwards to
% do well. This is the invertability regularization constraint
registration_score = sum( best_match_score_before2after ) * sum( best_match_score_after2before );

end % FUNCTION