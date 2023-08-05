%% main
function [ strand_subscripts_B, best_registration_score ] = register_strands( strand_subscripts_A, microns_per_voxel_A, lumen_radius_in_microns_range_A, ...
                                                                              strand_subscripts_B, microns_per_voxel_B, lumen_radius_in_microns_range_B, varargin )
%% register_strands                                                                                 
% registers the B vector set onto the A vector set and returns the B vector set in the coordinate system of the A vector set
% 
% % SAM 8/22/22 
%
% Note: ( if A and B are different sizes in 3-space, A should be the bigger one) !!!!! automate a warning ????
%
%
% optional initial search tolerances as varargin = { distance_tolerance, 
%                                                    angle____tolerance, 
%                                                    radius___tolerance  }
%
% Assumes: orthonormal subscript space (cubic voxel shape)
                                                                                                                            
% deltas = zeros( 1, 8 ); % 3 position (microns), 3 rotation axis, 1 rotation magnitude (radians), 1 vessel size dilation (dB)

%% resample the B vector set to match the resolution of A (each will retain its original size LUT: (lumen_radius_in_pixels_range_X) 
microns_per_voxel = microns_per_voxel_A ;

if ~ all(    microns_per_voxel_A ...
          == microns_per_voxel_B )
    
    resolution_factors = [    microns_per_voxel_B    ...
                           ./ microns_per_voxel_A, 1 ]; % 3 SAPCE, 1 SCALE
    
%     [ ~, ~, strand_subscripts_B ] = resample_vectors( lumen_radius_in_microns_range_B ./ microns_per_voxel_B, resolution_factors, ...
%                                                                   strand_subscripts_B, ones( 1, 3 )        );
    [ ~, ~, strand_subscripts_B ] = resample_vectors( ones( 5, 3 ), resolution_factors, strand_subscripts_B, ...
                                                      ones( 1, 3 )                                           );
                                        
end

%% extract directions of vessels from both vector sets                                              

% parse out the positional subscripts (from the dilational 4th component)
strand_space_subscripts_A = cellfun( @( x ) x( :, 1 : 3 ), strand_subscripts_A, 'UniformOutput', false ); 
strand_space_subscripts_B = cellfun( @( x ) x( :, 1 : 3 ), strand_subscripts_B, 'UniformOutput', false ); 

% microns_per_voxel = [ 1, 1, 1 ];
[ strand_directions_A ] = get_vessel_directions_V3( strand_space_subscripts_A, microns_per_voxel );
[ strand_directions_B ] = get_vessel_directions_V3( strand_space_subscripts_B, microns_per_voxel );

%% convert cell format to matrix format with directional quantities
matrix_directions_A = cell2mat( strand_directions_A );
matrix_directions_B = cell2mat( strand_directions_B );
matrix_subscripts_A = cell2mat( strand_subscripts_A );
matrix_subscripts_B = cell2mat( strand_subscripts_B );

% concatenate the directional quantities onto the spatial
vectors_A = [ matrix_subscripts_A, matrix_directions_A ]; % 3 position, 1 size, 3 direction
vectors_B = [ matrix_subscripts_B, matrix_directions_B ];

% remember the key for returning back to cell array from matrix format
cell_lengths_B = cellfun( @( x ) size( x, 1 ), strand_subscripts_B, 'UniformOutput', false ); 

% distance_stdevs_A = stdev( matrix_subscripts_A );
% distance_stdevs_B = stdev( matrix_subscripts_B );
% 
% distance_tolerance_default = (   sum( distance_stdevs_A .^ 2 )        ...
%                                + sum( distance_stdevs_B .^ 2 )) ^ 0.5 ;
% angle____tolerance_default = 2 * pi ;
% radius___tolerance_default = 2      ;
% 
% if ~ isempty( varargin{ 1 }), distance_tolerance = varargin{ 1 }; else, distance_tolerance = distance_tolerance_default ; end
% if ~ isempty( varargin{ 2 }), angle____tolerance = varargin{ 2 }; else, angle____tolerance = angle____tolerance_default ; end
% if ~ isempty( varargin{ 3 }), radius___tolerance = varargin{ 3 }; else, radius___tolerance = radius___tolerance_default ; end
% 
% is_optimizing = true ;
%     
% distance_tolerance = distance_tolerance ;
% angle____tolerance = angle____tolerance ;
% radius___tolerance = radius___tolerance ;

% % reset the subscripts to minimize the size of the reference image % !!!!! ????? should be sparse
% delta_position_offset = min( matrix_subscripts_B( :, 1 : 3 )) ...
%                       - min( matrix_subscripts_A( :, 1 : 3 ));
%                     
% deltas( 1 : 3 ) = - delta_position_offset ;
                    
% original_intersection_lower = max( min( matrix_subscripts_A ), min( matrix_subscripts_B ));
% original_intersection_upper = min( max( matrix_subscripts_A ), max( matrix_subscripts_B ));

% is_vector_in_original_intersection_A = matrix_subscripts_A >= original_intersection_lower & matrix_subscripts_A <= original_intersection_upper ;
% is_vector_in_original_intersection_B = matrix_subscripts_B >= original_intersection_lower & matrix_subscripts_B <= original_intersection_upper ;

%% set the origin of A to the smallest subscript, force minimal positive subscripts for A vectors 
vectors_B = vectors_B( :, 1 : 3 ) - min( vectors_A ) + 1 ;
vectors_A = vectors_A( :, 1 : 3 ) - min( vectors_A ) + 1 ;

% center_of_mass_A  =  mean( matrix_subscripts_A( is_vector_in_original_intersection_A ));
% center_of_mass_B  =  mean( matrix_subscripts_B( is_vector_in_original_intersection_B ));
% 
% % calculate initial displacement in center of mass
% deltas( 1 : 3, 1 ) = center_of_mass_B ...
%                    - center_of_mass_A ;

size_of_image_A = max( vectors_A ); %...
%                - min( matrix_subscripts_A );

% index_image_A = zeros( size_of_image_A ); % ?????? make sparse !!!!!!

locations_A = sub2ind( size_of_image_A, vectors_A( :, 1 ), ...
                                        vectors_A( :, 2 ), ...
                                        vectors_A( :, 3 )) ;

Num_vectors_B = size( vectors_B, 1 );
Num_vectors_A = size( vectors_A, 1 );

index_image_A( locations_A ) = 1 : Num_vectors_A ;

%% set tolerances for first round of searching 
positional_tolerance = 50 ; % microns
% dilational_tolerance = log(  2 ) ...
%                      / log( 10 ); % dilation factor [dB] %% !!!!! start this out much larger to decouple the optimizations which are not intrinsically coupled !!!!!
dilational_tolerance = log( 5 ) / log( 10 );
% dilational_tolerance = 1 ; % 1dB = 10 x fold change
% rotational_tolerance = pi / 2 ; % radians, 90 deg
% rotational_tolerance = pi / 4 ; % radians, 45 deg
rotational_tolerance = pi / 8 ; % radians, 22.5 deg

% !!!!!! add stretching factors ???????

sigmas = 2 * [ positional_tolerance, ...
               dilational_tolerance, ...
               rotational_tolerance ]; % 1 positional, 1 size, 1 directional
           
tolerance_octave = 0 ;

number_of_tolerance_octaves = 5 ;

tolerances_in_microns_range = positional_tolerance * 2 .^ - ( 0 : number_of_tolerance_octaves - 1 )';

% strels = generate_strels( positional_tolerance, microns_per_voxel, number_of_tolerance_octaves ); % row vectors

strels = calculate_linear_strel_range( size_of_image_A, microns_per_voxel, tolerances_in_microns_range );

% best_deltas = [ ];

%% initialization of perturbations FOR loop 
Num_perturbations = 1e6 ;
    
perturbation_Indcs = 1 : Num_perturbations ;

center_of_rotation = mean( vectors_B( :, 1 : 3 )); % updated to remain at the center of mass of the vector B set
       
is_optimizing = true ;

%% registration optimization loop 
while is_optimizing
    %% vectorized perturbation calculation:    
    sigmas = sigmas / 2 ;
    
	tolerance_octave = tolerance_octave + 1 ;
    
    is_optimizing = tolerance_octave < number_of_tolerance_octaves ;
        
                                                  %N%ormal dist.
                                                  %!%
    vector_deltas = [          ones( 1, 3 ) .* randn( Num_perturbations, 3 ) * sigmas( 1 ) , ... 3 positional [microns]
                      10 .^ (                  randn( Num_perturbations, 1 ) * sigmas( 2 )), ... 1 dilational [fold-change]
                             [ 2*pi, 1, 0 ] .* randu( Num_perturbations, 3 )               , ... 3 rotational axis (for now it is polar angle and Z-component)
                                               randn( Num_perturbations, 1 ) * sigmas( 3 )   ];% 1 rotational angle
                                                  %!%
                                                  %U%niform dist.

    % Z-component is uniformly sampled in [-1,1], then transformed to azimuthal angle to achieve isotropic sampling of spherical surface
	vector_deltas( :,   6   ) =  asin( vector_deltas( :, 6 )) .* sign( randu( Num_perturbations, 1 ) - 0.5 ); % random sign (+/-) attached here to extend [0,1]->[-1,1]
               
    % conversion from polar and azimuthal angle to x, y, z coordinates
	vector_deltas( :, 5 : 7 ) = [ cos( vector_deltas( :, 6 )) .* cos( vector_deltas( :, 5 )), ...
                                  cos( vector_deltas( :, 6 )) .* sin( vector_deltas( :, 5 )), ...
                                  sin( vector_deltas( :, 6 ))                                 ];
                                                
%     registration_scores = zeros( 1, Num_perturbations );
    
    %% resample ctrl points % ???? wrap this in another batch FOR ??????
    Num_control_points = 1e2 ;
    
    registration_scores_at_ctrl_pt = zeros( Num_perturbations,    ...
                                            Num_control_points );
    
	control_vector_Og_Idcs = randperm( Num_vectors_B, Num_control_points );

    %% iterate through random perturbations onto the current B->A mapping 
    for perturbation_Idx = perturbation_Indcs
        %% perturb the control points
        [ vectors_C, lumen_radius_in_microns_range_C ] = perturb( vectors_B( control_vector_Og_Idcs, : ), vector_deltas( perturbation_Idx, : ), microns_per_voxel, lumen_radius_in_microns_range_B, center_of_rotation );
        
        vectors_C( :, 1 : 3 )         = round( vectors_C( :, 1 : 3 ));
        
        is_ctrl_vector_within_bounds  = all(   vectors_C( :, 1 : 3 ) <= size_of_image_A   ...
                                             & vectors_C( :, 1 : 3 ) <=    ones( 1, 3 ),  ...
                                             2                                                    );
                                          
        vectors_C                     =        vectors_C( is_ctrl_vector_within_bounds, : );
        
        linears_C                     = sub2ind(                        size_of_image_A,  ...
                                               vectors_C( :,   1   ),                     ...
                                               vectors_C( :,   2   ),                     ...
                                               vectors_C( :,   3   )                      );
        
        vector__C_idcs     =    1 : size( vectors_C, 1 );
        
        Num_vectors_C_in_bounds = numel( vector__C_idcs );
                                                        
        %% loop thru control points from B vectors
        for control_vector_idx = vector__C_idcs
            %% search for matching A vector initialization
            
%             vector__C = vectors_C( control_vector_idx, : );
% 
%             is_searching_for_match = true ;
%             
%             previous_search_loop_index = 0 ;
%             
% %             search_loop_index = 0 ;
%             search_loop_index = vector__C( 4 ) - 1 ;
%             
%             vector_A_matches = [ ];
% 
%             %% loop thru sizes, ever-increasing the search radius, until finding at least 2 matches 
%             %%% !!!!!!!! convert this to a FOR loop and dynamically redefine the search strels LUT from the spatial tolerance
%             while is_searching_for_match
%                 
%                 search_loop_index = search_loop_index + 1 ;
%                 
% %                 % ???? alternative method perhaps better? ????
% %                 control_point_matches = [ control_point_matches; ...
% %                                           find( any( locations_A ==   control_point_linear       ...
% %                                                                   + strel{ search_loop_index },  ...
% %                                                      2                                          ))];
%                 
%                 % reference the index image
%                 vector_A_matches = [ vector_A_matches ;                                         ...
%                                      nonzeros( index_image_A(   linears_C( control_vector_idx ) ...
%                                                               + cell2mat( strel( 1 + previous_search_loop_index ...
%                                                                                  :            search_loop_index ))))];
                vector_A_matches = nonzeros( index_image_A(   linears_C( control_vector_idx ) ...
                                                            +    strels(  tolerance_octave  )));
           %% IF at least two matches
            if numel( vector_A_matches ) >= 2 

                positional_epsilons           =             microns_per_voxel                     ...
                                                       .*(   vector__C(                   1 : 3 ) ...
                                                           - vectors_A( vector_A_matches, 1 : 3 )); % L2 norm (microns)

                dilational_epsilon                                                                                ...
                    = log(   lumen_radius_in_microns_range_C( vector__C(                     4   )) ...
                           / lumen_radius_in_microns_range_A( vectors_A( vector_A_matches,   4   )))...
                    / log(                                          10                             ); % fold difference [dB]

                % ||A x B|| = ||A|| ||B|| ||sin(Th)||
                % rotational_epsilon = Th = arcsin(A x B); % angle [radians] % A and B should be unit vectors
                rotational_epsilon =          arcsin( sum(( cross( vector__C(                   5 : 7 ), ... % ??? may need to extend the vectorized arithmetic from compatible sizes to identical ???
                                                                   vectors_A( vector_A_matches, 5 : 7 ))) .^ 2, ...
                                                           2                                                    ) .^ 0.5 );

                positional_epsilon = sum( positional_epsilons .^ 2 ) .^ 0.5 ; % pythagorean sum

                epsilons = [ positional_epsilon, ...
                             dilational_epsilon, ...
                             rotational_epsilon  ];

                registration_scores_matches = prod( exp( - ( epsilons ./ sigmas ) .^ 2 / 2 ), 2 );

%                 % pick top two matches
%                 top_score = zeros( 2, 1 );
%                 top_match = zeros( 2, 1 );
% 
%                 for idx = 1 : 2
% 
%                     [ top_score( idx ), top_match( idx )] = max( registration_scores_matches );
% 
%                     registration_scores_matches( top_match( idx )) = -1 ; % 0 is the minimum registration score
% 
%                 end
% 
%                 registration_score_top_two_matches = geomean( top_score );

                % registration_scores_ctrl_pt = 

%                     is_searching_for_match = false ;

%                 registration_scores_at_ctrl_pt( perturbation_Idx, control_vector_idx ) = registration_score_top_two_matches / Num_vectors_C_in_bounds ; 
                registration_scores_at_ctrl_pt( perturbation_Idx, control_vector_idx ) = max( registration_scores_matches / Num_vectors_C_in_bounds ); 
                % ???? decrement Num_vectors_C_in_bounds when this statement is not reached for a ctrl vector ?????

            end % IF match(es) found                
%             end % match searching WHILE loop
        end % control point FOR loop        
    end % batch FOR loop
    
    %% sum registration scroes across ctrl points and choose the best perturbation 
           registration_scores   = sum( registration_scores_at_ctrl_pt, 2 );
    [ best_registration_score,                                            ...
      best_registration_idx    ] = max( registration_scores               );
    
%     best_previous_deltas = [ deltas( best_previous, : ),   sum( deltas .* registration_scores ) ...
%                                                          / sum(           registration_scores )];
    
    %% perturb the entire B vector set (not just the control points) toward the best registration 
    [ vectors_B, lumen_radius_in_microns_range_B, center_of_rotation ] = perturb( vectors_B, vector_deltas( best_registration_idx, : ), microns_per_voxel, lumen_radius_in_microns_range_B );
    
    % ???????????????????? recalculate registration score with all vectors (?and vice versa?) and wrap that in a batch FOR loop ???????????????????

%     best_deltas         = [ best_deltas,                deltas( best_registration_score_idx, : )];
        
end % optimization WHILE loop

%% output formatting 
strand_subscripts_B = mat2cell( vectors_B( :, 1 : 4 ), cell_lengths_B, 4 ); % output positional and dilational subscripts 

end % main

%% supporting fxns:
    %% perturb subscripts given delta vector 
function [ vectors, lumen_radius_in_microns_range, center_of_rotation ] = perturb( vectors, vector_D, microns_per_voxel, lumen_radius_in_microns_range, varargin )
    % vector_Delta is 1 x 8 vector, 3 position, 1 radius scaling, 3 rotational axis unit vector, 1 rotational axis angle
    
    % convention to subtract Delta from B vectors: A + Delta = B <=> A = B - Delta
    
    % center of roatation determination
    if ~isempty( varargin ), center_of_rotation =           varargin{ 1 };        % center of rotation is *optional input*
    else,                    center_of_rotation = mean( vectors( :, 1 : 3 )); end % center of rotation is *center of mass*
    % !!!!!!!!!!! alternatively, just translate the center of rotation with the otehr vectors instead of recalculating
    
%     rotation_matrix = [ cos( vector_D( 8 )) + vector_D( 5 ) ^ 2 * ( 1 - cos( vector_D( 8 ))), ]; 
    % where, theta = vector_D( 8 ), u = vector_D( 5 : 7 ), formula at: https://en.wikipedia.org/wiki/Rotation_matrix#Rotation_matrix_from_axis_and_angle
    rotation_matrix =       cos( vector_D( 8 ))  *                                eye( 3 )       ... %   identity    matrix
                    +       sin( vector_D( 8 ))  * cross( vector_D( 5 : 7 ),      eye( 3 )     ) ... % cross-product matrix of u
                    + ( 1 - cos( vector_D( 8 ))) * times( vector_D( 5 : 7 ), vector_D( 5 : 7 )') ;   % outer-product matrix of u
    
    %% rotational (vessel positions around center of rotation), convert into and out of real space [miocrons] dimensions
    vectors( :, 1 : 3 ) = (( microns_per_voxel .* ( vectors( :, 1 : 3 ) - center_of_rotation )) * rotation_matrix ) ./ microns_per_voxel ... * row vector on the left => rotates the opposite direction (equivalent to switching the sign of sine(theta) in rotation matrix defn.
                        +                                                 center_of_rotation                                             ;
    
    %% positional  
    vectors( :, 1 : 3 ) =   vectors( :, 1 : 3 ) - vector_D( 1 : 3 ) ...
                                               ./ microns_per_voxel ; 
    %% rotational (vessel directions)
    vectors( :, 5 : 7 ) =                                             vectors( :, 5 : 7 )  *       cos( vector_D( 8 )) ...
                        -                   cross( vector_D( 5 : 7 ), vectors( :, 5 : 7 )) *       sin( vector_D( 8 )) ...
                        + vector_D( 5 : 7 ) * dot( vector_D( 5 : 7 ), vectors( :, 5 : 7 )) * ( 1 - cos( vector_D( 8 ))); % Rodrigues' rotation Formula (with negative angle argument ( vector_D( 8 )))
    
    %% dilational
%     vectors_B( :, 4 ) = vectors_B( :, 4 ) % keep the size indices the same, only change the LUT entries
    lumen_radius_in_microns_range = lumen_radius_in_microns_range / vector_D( 4 );
    
end

%     %% generate structuring elements for searching at ever-incresing radii from an origin 
% function [ strels ] = generate_strels( position_tolerance,  microns_per_voxel, number_of_tolerance_octaves )
% 
%     strels = cell( number_of_tolerance_octaves, 1 );
%     
%     for 1 : number_of_tolerance_octaves
%         
%        strel_radii = position_tolerance ./ microns_per_voxel ;
%        
%        
%        
%        strels{ octave_idx } = 
%        
%        position_tolerance = position_tolerance / 2 ;
%         
%         
%     end
% 
%     %
% 
% end