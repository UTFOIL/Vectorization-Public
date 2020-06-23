function [ bifurcation_vertices, vertices_in_strands,                                               ...
           edge_indices_in_strands, end_vertices_in_strands ]                                       ...
                                                           = get_network_V190( vertex_pairs_of_edges )
%% get_probabilities:

% V160:
% The purpose of this function is essentially to clean up the output of the get_edges and
% get_vertices functions.  The strategy is to break up the image into networks based on their
% connectivity (from an adjacency matrix formed by the edges that were found before).  We will view
% each network as a stochastic process where the transition matrix is formed by weighting the
% adjacency matrix by the boltzmann factors associated with the average energies on the edges
% connecting two vertices (taken as the best trajectory found between each vertex based on average
% energy).  The energies are normalized in the same way as the get_edges algorithm, where they are
% divided by the origin vertex energy scaled by some characteristic fraction. SAM 5/11/18
%
% V161 where the edges have already been chosen before this function in the function choose_edges
% after the get_edges function, and what remains is to make the transition matrix from all of the
% edges that are passed. SAM 5/15/18
%
% function name changed from get_probabilities to get_strands. SAM 5/22/18
%
% V170 in which transitions from vertex A to A are not allowed, and the purpose of the transition
% matrix is for each vertex to select its top two vertices.  Vertices that select each other then
% belong to the same strand.  These strands are then easily found by looking at the connected
% components of transition matrix with only the top two edges from each vertex nonzero.  SAM 5/22/18
%
% function name changed from get_strands to get_network SAM 5/31/18
%
% V180 in which the junctions are also retrieved (in addition to the strands) by looking at the
% edges in the top 4 connectivity network that meet all of the following conditions
%
% 1) not found in the top 2 connectivity network
% 2) not connecting a strand to itself (all single vertex strands are considered the same strand)
%
% The connected components of this network are junctions.  The junction sets need to be trimmed to
% have only one edge connecting strand A to B again considering all single vertex strands to be the
% same strand).  After trimming this network is analyzed for connectivity. Any single vertex strands
% that dont contact at least two other strands are removed
%
% % Once we have the junctions, loop through the strands and break them up into the smaller strands
% % just going from junction to junction !!!!! not yet implemente !!!!!
%
% SAM 5/31/18
%
% V181 in which the way the junctions are found is changed.  In this version, we only look at top 3
% connectivity for vertices that had two neighbors in the the top two connectivity network.  The
% edges in this set that were not in the top two connectivity network are the first junction pieces.
% 
%
% !!! not yet implemented section !!!
% Iterating from here to now look at the top 4 connectivity for vertices that had three neighbors in
% the previous connectivity network that included top 3 connectivity for some of the vertices.
% !!! not yet implemented section !!!
%
% SAM 6/12/18
%
% V182 in which an outputs are added, end_vertices_of_junctions and bifurcation_vertices. 
%
% Also removed bugs: 
%
% extra edges were sometimes added to strands in previous version (two neighbors max was not
% guaranteed)
%
% the requirement for the opportunity to form a junction was less strict than intended (some singly
% connected vertices were granted extra edges in the top three connectivity network)
%
% the whole section about edge trimming inside junctinos was removed because I now believe that no
% edges would have been removed in that section anyway.
%
% SAM 7/11/18
%
% Attempting to decrease the runtime in this version. Taking a more direct approach in the FOR
% loops.  Building the adjacency matrices from the strand from bottom up instead of starting with
% the network-wide one and widdling it down. SAM 8/14/18
%
% V190 Complete Overhaul: The network objects of interest are changed and the method for finding
% them has changed. Now there are two objects only: strands and bifurcation vertices.  Strands are
% end-to-end collections of edges that don't bifurcate, and bifurcation vertices are vertices where
% those strands end. These are found by treating all edges as true and then manipulating the
% adjacency matrix to find the vertices with two neighbors.  The connected components when only
% considering two-neighbor vertices gives the strand vertices and associated edges.


% create sparse adjacency matrix from which to calculate the number of neighbors
vertex_pairs_of_edges = double( vertex_pairs_of_edges );

number_of_edges = size( vertex_pairs_of_edges, 1 );

number_of_vertices = max( vertex_pairs_of_edges( : ));

edge_index_range  = 1 : number_of_edges ;

% The first vertex in the pair is given by the row, second by the column.
edge_lookup_table                            ...
    = sparse( vertex_pairs_of_edges( :, 1 ), ...
              vertex_pairs_of_edges( :, 2 ), ...
              edge_index_range',             ...
              number_of_vertices,            ...
              number_of_vertices,            ...
              number_of_edges                );

adjacency_matrix                                     ...
    = sparse( vertex_pairs_of_edges( :, 1 ),         ...
              vertex_pairs_of_edges( :, 2 ),         ...
              ones( number_of_edges, 1, 'logical' ), ...
              number_of_vertices,                    ...
              number_of_vertices,                    ...
              number_of_edges                        );
          
% % fill out the backwards directions too so the matrix is symmetric (average energy is same from A to
% % B as B to A).          
adjacency_matrix = adjacency_matrix | adjacency_matrix' ;
                             
number_of_edges_at_vertex = sum( adjacency_matrix, 2 );

% If a vertex has 1 asoociated edge, it is an endpoint. If it has two it is an interior point on a
% strand. If it has n where n > 2, then it is a bifurcation of degree (n-1).  e.g. if a vertex has
% three neighbors then it is a degree two bifurcation.
vertex_is_endpoint    = number_of_edges_at_vertex == 1 ;
vertex_is_interior    = number_of_edges_at_vertex == 2 ;
vertex_is_bifurcation = number_of_edges_at_vertex >= 3 ;
vertex_is_isolated    = number_of_edges_at_vertex == 0 ;

bifurcation_vertices = find( vertex_is_bifurcation );

bifurcation_or_endpoint_vertices = find( vertex_is_bifurcation | vertex_is_endpoint );

non_strand_interior_vertices = find( vertex_is_bifurcation | vertex_is_endpoint | vertex_is_isolated )' ;
                      
% zero out the vertices that represent bifurcations of any degree or endpoints to create a new
% adjacency matrix of only strand interior vertices.
adjacency_matrix_only_strand_interiors = adjacency_matrix ;

adjacency_matrix_only_strand_interiors( bifurcation_or_endpoint_vertices, : ) = false ;
adjacency_matrix_only_strand_interiors( :, bifurcation_or_endpoint_vertices ) = false ;

% convert the adjacency matrix to a graph to find the connected components
strand_graph = graph( adjacency_matrix_only_strand_interiors );

% The connected components of this graph are the strand objects.
interior_vertices_in_strands = conncomp( strand_graph, 'OutputForm', 'cell' );

interior_vertices_in_strands_list = cell2mat( interior_vertices_in_strands );

number_of_vertices_in_strands = cellfun( 'length', interior_vertices_in_strands );

cumulative_number_of_vertices_in_strands = cumsum( number_of_vertices_in_strands );

% !!! this removal looks more complicated than it needs to be and somewhat redundant. SAM 6/11/19

% remove strands that are single nodes consisting of a bifurcation, endpoint, or isolated vertex
strand_indices_to_remove = zeros( length( non_strand_interior_vertices ), 1 );

counter = 0 ;

for vertex_index = non_strand_interior_vertices
    
    list_index = find( interior_vertices_in_strands_list == vertex_index, 1 );
    
    counter = counter + 1 ;
    
    strand_indices_to_remove( counter ) = find( cumulative_number_of_vertices_in_strands >= list_index, 1 );
    
end % FOR vertex index

interior_vertices_in_strands( strand_indices_to_remove ) = [ ];

number_of_strands = length( interior_vertices_in_strands );

  edge_indices_in_strands = cell( number_of_strands, 1 );
  end_vertices_in_strands = cell( number_of_strands, 1 );
      vertices_in_strands = cell( number_of_strands, 1 );
      
strand_index_range = 1 : number_of_strands ;
                     
for strand_index = strand_index_range
    
    possible_vertex_pairs_for_strand_interior                                                   ...
                    =   interior_vertices_in_strands{ strand_index }'                           ...
                    + ( interior_vertices_in_strands{ strand_index } - 1 ) * number_of_vertices ;
                    
    % build a small (in one dimension) adjacency matrix just for this strand
    adjacency_matrix_for_strand = adjacency_matrix( interior_vertices_in_strands{ strand_index }, : );
    
    adjacency_matrix_for_strand_exterior = adjacency_matrix_for_strand ;
    
    % zero out the interior in the other direction to look for the strand exterior, if it exists
    adjacency_matrix_for_strand_exterior( :, interior_vertices_in_strands{ strand_index }) = 0 ;
        
    % find the end vertices for this strand (vertices with only one neighbor (not two)).  This
    % won't exist for all strands in general.   Strands that form loops will not have end vertices.
    number_of_neighbors_at_strand_exterior = sum( adjacency_matrix_for_strand_exterior, 1 );
    
    if number_of_neighbors_at_strand_exterior == 0
        
        % !!! this is redundant.  In next version just put possible_vertex_pairs_for_strand_interior
        % directly into nonzeros( edge_lookup_table( vertex_pairs_for_strand_interior )) SAM 6/11/19
                

        % build a small adjacency matrix just for this strand        
        adjacency_matrix_for_strand = full( adjacency_matrix( possible_vertex_pairs_for_strand_interior ));

        [ vertices_in_strand_row, vertices_in_strand_column ] = find( adjacency_matrix_for_strand );        
        
        vertex_pairs_for_strand_interior =   interior_vertices_in_strands{ strand_index }( vertices_in_strand_row    )                            ...
                                         + ( interior_vertices_in_strands{ strand_index }( vertices_in_strand_column ) - 1 ) * number_of_vertices ;        
                
        edge_indices_in_strands{ strand_index } = nonzeros( edge_lookup_table( vertex_pairs_for_strand_interior ));

        % remove the worst edge from this strand to break the cycle. (Edge indices are sorted by
        % energy.)
        worst_edge_index = max( edge_indices_in_strands{ strand_index });
        
        % remove the two vertices associated with this worst edge from the list of interior vertices
        interior_vertices_in_strands{ strand_index }( interior_vertices_in_strands{ strand_index } == vertex_pairs_of_edges( worst_edge_index, 1 )) = [ ];
        interior_vertices_in_strands{ strand_index }( interior_vertices_in_strands{ strand_index } == vertex_pairs_of_edges( worst_edge_index, 2 )) = [ ];
        
        possible_vertex_pairs_for_strand_interior                                                   ...
                        =   interior_vertices_in_strands{ strand_index }'                           ...
                        + ( interior_vertices_in_strands{ strand_index } - 1 ) * number_of_vertices ;
    
        % build a small adjacency matrix just for this strand
        adjacency_matrix_for_strand = adjacency_matrix( interior_vertices_in_strands{ strand_index }, : );

        adjacency_matrix_for_strand_exterior = adjacency_matrix_for_strand ;

        adjacency_matrix_for_strand_exterior( :, interior_vertices_in_strands{ strand_index }) = 0 ;

        % find the end vertices for this strand (vertices with only one neighbor (not two)).  This
        % won't exist for all strands in general.   Strands that form loops will not have end vertices.
        number_of_neighbors_at_strand_exterior = sum( adjacency_matrix_for_strand_exterior, 1 );

    end % IF cyclical strand
           
    end_vertices_in_strands{ strand_index } = find( number_of_neighbors_at_strand_exterior );
    
    % build a small adjacency matrix just for this strand
    adjacency_matrix_for_strand_interior = full( adjacency_matrix( possible_vertex_pairs_for_strand_interior ));

    [ vertices_in_strand_row, vertices_in_strand_column ] = find( adjacency_matrix_for_strand_interior );
    
    % !!! this is redundant.  In next version just put possible_vertex_pairs_for_strand_interior
    % directly into nonzeros( edge_lookup_table( vertex_pairs_for_strand_interior ))  SAM 6/11/19   
    vertex_pairs_for_strand_interior =   interior_vertices_in_strands{ strand_index }( vertices_in_strand_row    )                            ...
                                     + ( interior_vertices_in_strands{ strand_index }( vertices_in_strand_column ) - 1 ) * number_of_vertices ;
                    
    number_of_end_vertices = length( end_vertices_in_strands{ strand_index });
                       
    number_of_edges_at_end_vertices = number_of_neighbors_at_strand_exterior( end_vertices_in_strands{ strand_index });
    
    vertex_pairs_for_endpoints_local = zeros( 1, sum( number_of_edges_at_end_vertices ));
    
    for end_vertex_ordinate = 1 : number_of_end_vertices
        
        for edge_ordinate = 1 : number_of_edges_at_end_vertices( end_vertex_ordinate )
            
            vertex_pairs_for_endpoints_local( end_vertex_ordinate + edge_ordinate - 1 ) = find( adjacency_matrix_for_strand( :, end_vertices_in_strands{ strand_index }( end_vertex_ordinate )), 1 );
            
            adjacency_matrix_for_strand( vertex_pairs_for_endpoints_local( end_vertex_ordinate + edge_ordinate - 1 ), end_vertices_in_strands{ strand_index }( end_vertex_ordinate )) = 0 ;
        
        end % FOR edge ordinate
    end % FOR end vertex ordinate
        
    vertex_pairs_for_endpoints = [ end_vertices_in_strands{ strand_index } + ( interior_vertices_in_strands{ strand_index }( vertex_pairs_for_endpoints_local ) - 1 ) * number_of_vertices , ...
                                   interior_vertices_in_strands{ strand_index }( vertex_pairs_for_endpoints_local ) + ( end_vertices_in_strands{ strand_index } - 1 ) * number_of_vertices   ];  
    
    vertex_pairs_for_strand = [ vertex_pairs_for_strand_interior, vertex_pairs_for_endpoints ];
                        
    % get original edge indices for the edges in this strand
    edge_indices_in_strands{ strand_index } = nonzeros( edge_lookup_table( vertex_pairs_for_strand ));
    
    vertices_in_strands{ strand_index } = [ interior_vertices_in_strands{ strand_index }, ...
                                                 end_vertices_in_strands{ strand_index }  ];
    
end % FOR strand

%% go back and assign every edge that connects (bifurcation or endpoint) to (bifurcation or endpoint) to its own strand

% zero out the vertices that represent strand interiors to create a new adjacency matrix of only
% bifurcations and endpoint vertices.
adjacency_matrix_without_strand_interiors = adjacency_matrix ;

adjacency_matrix_without_strand_interiors( vertex_is_interior, : ) = false ;
adjacency_matrix_without_strand_interiors( :, vertex_is_interior ) = false ;

vertex_pairs_in_extra_strands = find( adjacency_matrix_without_strand_interiors )';

edge_indices_in_extra_strands_vector = nonzeros( edge_lookup_table( vertex_pairs_in_extra_strands ));

vertices_in_extra_strands_matrix = vertex_pairs_of_edges( edge_indices_in_extra_strands_vector, : );

number_of_extra_strands = length( edge_indices_in_extra_strands_vector );

edge_indices_in_extra_strands = num2cell( edge_indices_in_extra_strands_vector );

    vertices_in_extra_strands = mat2cell( vertices_in_extra_strands_matrix, ones( number_of_extra_strands, 1 ), 2 );
    
end_vertices_in_extra_strands = vertices_in_extra_strands ;

edge_indices_in_strands = [ edge_indices_in_strands; edge_indices_in_extra_strands ];
end_vertices_in_strands = [ end_vertices_in_strands; end_vertices_in_extra_strands ];
    vertices_in_strands = [     vertices_in_strands;     vertices_in_extra_strands ];

end % function