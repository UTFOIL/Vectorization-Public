function [ edges2vertices, original_edge_indices ] = clean_edge_pairs( edges2vertices, edge_energies, is_keeping_only_mutual_edge_pairs )
%% Choosing best unique trajectories between vertices A and B.

    number_of_original_edges = size( edges2vertices, 1 );

    original_edge_indices = 1 : number_of_original_edges ;

    % % not possible to have these cases (A-A edge, A- edge, or edge with nonnegative max energy) with
    % the new deterministic edge search SAM 5/31/19
    % 
    % % remove trajectories that terminate on their original vertex
    % indices_of_non_self_referential_edges = find( edges2vertices( :, 1 ) ~= edges2vertices( :, 2 ));
    % 
    % edges2vertices          =          edges2vertices( indices_of_non_self_referential_edges, : );
    % % edge_sizes_temp       =       edge_subscripts( indices_of_non_self_referential_edges, 4 );
    % edge_energies_temp    =         edge_energies( indices_of_non_self_referential_edges    );
    % original_edge_indices = original_edge_indices( indices_of_non_self_referential_edges    );
    % 
    % % remove trajectories that don't land on any vertex
    % indices_of_terminal_edges = find( edges2vertices( :, 2 ) > 0 );
    % 
    % edges2vertices          =          edges2vertices( indices_of_terminal_edges, : );
    % % edge_sizes_temp       =       edge_sizes_temp( indices_of_terminal_edges    );
    % edge_energies_temp    =    edge_energies_temp( indices_of_terminal_edges    );
    % original_edge_indices = original_edge_indices( indices_of_terminal_edges    );
    % 
    % % remove edges that ever attained a negative energy on the trajectory
    % max_edge_energies = cellfun( @max,  edge_energies_temp );
    % 
    % indices_of_negative_energy_edges = find( max_edge_energies < 0 );

    [ mean_edge_energies ] = get_edge_metric( edge_energies );

    % 10/4/18 Possible Improvement: weight the edge means by a factor that is low when a certain
    % variance in the size coordinate of each trajectory is large.  This variance is first adjusted so
    % that any variance explained by a linear fit is removed:
    %
    % characteristic_size_stdev = scales_per_octave ; % needs to be input to this function
    %
    % % subtract off the constant term
    % edge_sizes_temp = cellfun( @(v) v - mean(v), edge_sizes_temp, 'UniformOutput', false )
    %
    % subtract off the (orthogonal) linear fits through the origin 
    % edge_sizes_temp = cellfun( @(v) v -       ((1:length(v))-(length(v)-1)/2)'    ...
    %                                     *     ((1:length(v))-(length(v)-1)/2)     ...
    %                                     * v/( ((1:length(v))-(length(v)-1)/2)     ...
    %                                          *((1:length(v))-(length(v)-1)/2)' ), ...
    %                             edge_sizes_temp, 'UniformOutput', false           )
    % 
    % % calculate the variances of the sizes, take negative exponent, and multiply this factor by the
    % % mean energy
    % mean_edge_energies = cellfun( @mean,  edge_energies_temp ) .* exp( - cellfun( @stdev, edge_sizes_temp ) / characteristic_size_stdev, edge_sizes_temp ));

    % mean_edge_energies    =     mean_edge_energies( indices_of_negative_energy_edges    );
    % original_edge_indices = original_edge_indices(  indices_of_negative_energy_edges    );
    % edges2vertices        =   edges2vertices(       indices_of_negative_energy_edges, : );

    % pre-sorting trajectories by length from shortest to longest, so that any A->B and B->A with the
    % same energy will go to the shorter one.  This is because the longer one may be double counting
    % some stretch.  Also, if you can make the trajectory in fewer steps with the same energy, why not?
    lengths_of_edges = cellfun( @length, edge_energies );

    [ ~, indices_sorted_by_length ] = sort( lengths_of_edges );

    mean_edge_energies    =    mean_edge_energies( indices_sorted_by_length    );
    original_edge_indices = original_edge_indices( indices_sorted_by_length    );
    edges2vertices        =        edges2vertices( indices_sorted_by_length, : );

    % sorting the trajectories by activation energy in ascending order, so lowest (best) energies at the
    % top of the list.
    [ ~, indices_sorted_by_mean_energy ] = sort( mean_edge_energies );

    original_edge_indices = original_edge_indices( indices_sorted_by_mean_energy    );
    edges2vertices        =        edges2vertices( indices_sorted_by_mean_energy, : );

    % only possible to have multiple trajectories from A to B with deterministic edge search after
    % adding children vertices and new edges

    % choosing only the best trajectory from vertex A to B for all A, B in the set of all vertices.
    % UNIQUE chooses the first of the multiple instances, so top of the list/lowest energy is preferred.
    [ edges2vertices, indices_of_unique_edges ] = unique( edges2vertices, 'rows', 'stable' );

    original_edge_indices = original_edge_indices(  indices_of_unique_edges );
    % mean_edge_energies    =    mean_edge_energies( indices_of_unique_edges );

    % choosing the best of the two directions of edges in the cases where we have pairs of mutual
    % trajectories (from A to B and B to A).
    [ ~, mutual_edge_indices, reverse_mutual_edge_indices ]            ...
        = intersect([ edges2vertices( :, 1 ), edges2vertices( :, 2 )], ...
                    [ edges2vertices( :, 2 ), edges2vertices( :, 1 )], ...
                    'rows', 'stable'                                   );
    
    if is_keeping_only_mutual_edge_pairs
        
        is_mutual_edge_pair_to_be_kept                          ...
            = mutual_edge_indices < reverse_mutual_edge_indices ;          
        
        original_edge_indices = original_edge_indices( is_mutual_edge_pair_to_be_kept    );
               edges2vertices =        edges2vertices( is_mutual_edge_pair_to_be_kept, : );
        
    else % tossing the worse members of each mutual edge pair
                
        is_mutual_edge_pair_to_be_erased                        ...
            = mutual_edge_indices > reverse_mutual_edge_indices ;  

        original_edge_indices( mutual_edge_indices( is_mutual_edge_pair_to_be_erased )    ) = [ ];
        %    mean_edge_energies( mutual_edge_indices( is_mutual_edge_pair_to_be_erased )    ) = [ ];
               edges2vertices( mutual_edge_indices( is_mutual_edge_pair_to_be_erased ), : ) = [ ];

    end
           
end % FUNCTION clean_edge_pairs

