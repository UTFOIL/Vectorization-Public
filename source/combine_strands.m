%% combine strands SAM 9/21/21
% this function combines the strand objects so that they are more continuous. At each bifurcation,
% it combines the two largest strands at that bifurcation, with a possible relative size condition
% between the strands involved.

% threshold = 10 ;
% 
% %% loading energy settings file
% path_to_energy_settings  = 'E:\Annie\200923 2x2 vasculature RG med filter\200923 2x2 vasculature RG med filter\batch_200925-184104\settings\energy_200925-184104.mat' ;
% path_to_energy_data      = 'E:\Annie\200923 2x2 vasculature RG med filter\200923 2x2 vasculature RG med filter\batch_200925-184104\data\energy_200925-184104_Fused_medfilt_nobg.mat' ;
% % path_to_curated_vertices = 'E:\Annie\200923 2x2 vasculature RG med filter\200923 2x2 vasculature RG med filter\batch_200925-184104\vectors\curated_vertices_200925-184104_Fused_medfilt_nobg.mat' ;
% path_to_curated_edges    = 'E:\Annie\200923 2x2 vasculature RG med filter\200923 2x2 vasculature RG med filter\batch_200925-184104\vectors\curated_edges_201111-053420_Fused_medfilt_nobg.mat' ;
% path_to_network          = 'E:\Annie\200923 2x2 vasculature RG med filter\200923 2x2 vasculature RG med filter\batch_200925-184104\vectors\network_210328-032210_Fused_medfilt_nobg.mat' ;

function combine_strands( path_to_energy_settings, path_to_energy_data, path_to_curated_edges, path_to_network, threshold )

path_to_network_backup       = [ path_to_network(       1 : end - 4 ), '_backup.mat' ];
path_to_curated_edges_backup = [ path_to_curated_edges( 1 : end - 4 ), '_backup.mat' ];

load( path_to_energy_settings  )
load( path_to_energy_data      )
% load( path_to_curated_vertices )
load( path_to_curated_edges    ) % vertices are in here

%% loading and backing-up network file

if isfile( path_to_network_backup ) % do not overwrite on this save/ instead load the backup

    load( path_to_network_backup )
    
else
    
    load( path_to_network )

    save( path_to_network_backup, ...
                  'bifurcation_vertices'       , ...
                  'strand_subscripts'          , ...
                  'strand_energies'            , ...
                  'mean_strand_energies'       , ...
                  'vessel_directions'          , ...
                  'network_runtime_in_seconds' , ...
                  'network_statistics'         , ...
                  'strands2vertices'             );

end

%% loading and backing-up curated_vertices file

if isfile( path_to_curated_edges_backup ) % do not overwrite on this save/ instead load the backup

    load( path_to_curated_edges_backup )
    
else
    
    load( path_to_curated_edges )

    save( path_to_curated_edges_backup,              ...
                'edge_curation', ...
                'edge_curation_runtime_in_seconds', ...
                'edge_energies', ...
                'edge_scale_subscripts', ...
                'edge_space_subscripts', ...
                'edges2vertices', ...
                'mean_edge_energies', ...
                'vertex_energies',                    ...
                'vertex_scale_subscripts',            ... 
                'vertex_space_subscripts'             );

end

%% check for issues

% discontinuities

strand_max_delta_length = cellfun( @( x ) max( abs(   x( 2 : end    , 1 : 3 )                 ...
                                                    - x( 1 : end - 1, 1 : 3 )), [ ], 'all' ), ...
                                   strand_subscripts                                          );
                               
if any( strand_max_delta_length > 2)
    
    warning('discontinuity in strand')
    
    strand_subscripts{strand_max_delta_length > 2}
    strand_max_delta_length(strand_max_delta_length > 2)
    find(strand_max_delta_length > 2)
    
end

% bifurcation and strands2vertices mismatch

if ~ isempty( setdiff( bifurcation_vertices, strands2vertices( : )))
    
	warning('bifurcation not in strands2vertices LUT')
    
end

% vertices not matching with STRAND endpoints in 3 space
endpoint_locations  = cell2mat([ cellfun( @( x ) x(  1 , 1 : 3 ), strand_subscripts, 'UniformOutput', false ); ...
                                 cellfun( @( x ) x( end, 1 : 3 ), strand_subscripts, 'UniformOutput', false )  ]);

distance_from_vertex = max( abs( endpoint_locations - double( vertex_space_subscripts( strands2vertices( : ), : ))), [ ], 2 );
                             
if  any( distance_from_vertex > 2 )
    
    warning('strands2vertices vertex physically away from strand end')
    
    distance_from_vertex( distance_from_vertex > 2 )
    numel(find( distance_from_vertex > 2 ))
    
    % !!!!!!!!!! found a problem: the vertices identified in the strands2vertices LUT are not
    % colocated witht the ends of the strands:
    histogram( distance_from_vertex )
    
    % compare to random pairing of vertices identified in the strands2vertices LUT:
    figure
    distance_from_vertex = max( abs( double( vertex_space_subscripts( strands2vertices( randperm( numel( strands2vertices ))), : )) - double( vertex_space_subscripts( strands2vertices( : ), : ))), [ ], 2 );
    histogram( distance_from_vertex )
    
end

% % vertices not matching with EDGE endpoints in 3 space
% endpoint_locations  = cell2mat([ cellfun( @( x ) x(  1 , 1 : 3 ), edge_subscripts, 'UniformOutput', false ); ...
%                                  cellfun( @( x ) x( end, 1 : 3 ), edge_subscripts, 'UniformOutput', false )  ]);
% 
% distance_from_vertex = max( abs( endpoint_locations - double( vertex_space_subscripts( edges2vertices( : ), : ))), [ ], 2 );
%                              
% if  any( distance_from_vertex > 2 )
%     
%     warning('strands2vertices vertex physically away from strand end')
%     
%     distance_from_vertex( distance_from_vertex > 2 )
%     numel(find( distance_from_vertex > 2 ))
%     
%     % !!!!!!!!!! found a problem: the vertices identified in the strands2vertices LUT are not
%     % colocated witht the ends of the strands:
%     histogram( distance_from_vertex )
%     
%     % compare to random pairing of vertices identified in the strands2vertices LUT:
%     figure
%     distance_from_vertex = max( abs( double( vertex_space_subscripts( edges2vertices( randperm( numel( edges2vertices ))), : )) - double( vertex_space_subscripts( edges2vertices( : ), : ))), [ ], 2 );
%     histogram( distance_from_vertex )
%     
% end


%% combining strands

number_of_strands = size( strands2vertices, 1 );

strands = 1 : number_of_strands ;

% [ original_vertices, vertices ] = unique( strands2vertices( : ));
% 
% number_of_vertices = numel( vertices );
% 
% adjacency_matrix = logical( spalloc( number_of_vertices, number_of_vertices, number_of_strands ));
% 

vertices = unique( strands2vertices( : ))';

number_of_vertices = max( strands2vertices( : ));

number_of_combined_strands = 0 ;

counter = 0 ;

for vertex = vertices
    
    % find which strands are adjacent to this vertex (and which end of the strand it is)
    strand_ends_at_vertex = find( strands2vertices( : ) == vertex );
    
    if numel( strand_ends_at_vertex ) >= 2
        if numel( strand_ends_at_vertex ) == 2
            
            warning( 'input network had interior pointn marked as bifurcation' )
            
        end

        % keep track of how many strands are combined to assign unique vertex labels to them
        number_of_combined_strands = number_of_combined_strands + 1 ;

        strands_at_vertex =    mod( strand_ends_at_vertex - 1, number_of_strands ) + 1 ;
        
        is_flipped_strand = double( strand_ends_at_vertex   >  number_of_strands );
        
        % look for strand discontinuities from the adjacent vertex
        strand_end_centroid = [ 0, 0, 0 ];
                
        number_of_strands_at_vertex = numel( strands_at_vertex );
        
        for strand_idx = 1 : number_of_strands_at_vertex
           
            strand_end_centroid = strand_end_centroid + strand_subscripts{ strands_at_vertex( strand_idx )}(         ~ is_flipped_strand( strand_idx ) ...
                                                                                                             + end *   is_flipped_strand( strand_idx ), 1 : 3 );
            
        end
        
        strand_end_centroid = strand_end_centroid / number_of_strands_at_vertex ; % take mean
        
        distance_from_centroid = zeros( number_of_strands_at_vertex, 1 );
        
        for strand_idx = 1 : number_of_strands_at_vertex
            
            distance_from_centroid( strand_idx ) = max( abs(   strand_subscripts{ strands_at_vertex( strand_idx )}(         ~ is_flipped_strand( strand_idx ) ...
                                                                                                                    + end *   is_flipped_strand( strand_idx ), 1 : 3 ) ...
                                                             - strand_end_centroid ));
            
        end
                
        distance_from_centroid = round( distance_from_centroid / threshold ) * threshold ;
        
        [ ~, index_of_size_sort ] = sort( network_statistics.strand_ave_radii( strands_at_vertex ), 'descend' );
                
        [ ~, index_of_dist_sort ] = sort( distance_from_centroid( index_of_size_sort ));
        
        
        if max( distance_from_centroid ) > 0
            
            warning( 'discontinuity' )
            
            distance_from_centroid
            
        end
        
        % combine the two edges at this bifurcation with the largest radii
        strand_ends_to_combine = strand_ends_at_vertex( index_of_size_sort( index_of_dist_sort( 1 : 2 )));
        
        % assign the strands to be combined as adjacent to a new vertex, leave the other strands be
        strands2vertices( strand_ends_to_combine ) = number_of_vertices + number_of_combined_strands ;
        
        is_flipped_strand_to_combine = double( strand_ends_to_combine > number_of_strands );
        
        strands_at_vertex_to_combine = mod( strand_ends_to_combine - 1, number_of_strands ) + 1 ;
        
        if max( abs(   strand_subscripts{ strands_at_vertex_to_combine( 1 )}( ~ is_flipped_strand_to_combine( 1 ) + end * is_flipped_strand_to_combine( 1 ), 1 : 3 ) ...
                     - strand_subscripts{ strands_at_vertex_to_combine( 2 )}( ~ is_flipped_strand_to_combine( 2 ) + end * is_flipped_strand_to_combine( 2 ), 1 : 3 ))) ...
         > threshold
     
            warning('strand vertex mismatch')
            
            counter = counter + 1 ;
     
        end
        
                
    end % IF    
end % FOR

% remove redundant strand connections between two vertices
[ ~, IA ] = unique( strands2vertices, 'rows' );
to_remove = setdiff(1:max(IA),IA); % selects later strand (hopefully sorted by energy?)

mean_strand_energies( to_remove ) = [ ];

strands2vertices( to_remove, : ) = [ ];

strand_energies(   to_remove ) = [ ];
strand_subscripts( to_remove ) = [ ];
vessel_directions( to_remove ) = [ ];

% %% recalculate network for the newly created strand XL objects

[ bifurcation_vertices, ~,                                            ...
   strand_indices_in_strands_XL, end_vertices_of_strands_XL ]      = get_network_V190( strands2vertices );

% sort the strand output
[ vertex_indices_in_strands_XL, strand_indices_in_strands_XL, edge_backwards_in_strands_XL ]           ...
                              = sort_network_V180( strands2vertices, end_vertices_of_strands_XL, ...
                                                   strand_indices_in_strands_XL                  );

strands2vertices = [ cellfun( @( x ) x(  1  ), vertex_indices_in_strands_XL ), ...
                     cellfun( @( x ) x( end ), vertex_indices_in_strands_XL )  ];

% strand_subscripts = cellfun(      @(       space_subscripts,      scale_subscripts  )  ...
%                            double([      space_subscripts,        scale_subscripts ]), ...
%                                   strand_space_subscripts, strand_scale_subscripts,    ...
%                                           'UniformOutput', false                     );

[ strand_subscripts, strand_energies ]                                                                           ...
        = get_strand_objects( strand_subscripts, strand_energies, strand_indices_in_strands_XL, edge_backwards_in_strands_XL );

% remove strands using the new vertices

to_remove = find(   strands2vertices( :, 1 ) > number_of_vertices ...
                  | strands2vertices( :, 2 ) > number_of_vertices ); % selects later strand (hopefully sorted by energy?)

mean_strand_energies( to_remove ) = [ ];

strands2vertices( to_remove, : ) = [ ];

strand_energies(   to_remove ) = [ ];
strand_subscripts( to_remove ) = [ ];
vessel_directions( to_remove ) = [ ];
    
%% check for discontinuities

strand_max_delta_length = cellfun( @( x ) max( abs(   x( 2 : end    , 1 : 3 )                 ...
                                                    - x( 1 : end - 1, 1 : 3 )), [ ], 'all' ), ...
                                   strand_subscripts                                          );

number_of_strands = length( strand_subscripts );

number_of_vertices = max( strands2vertices( : ));

if any( strand_max_delta_length > 10 )
    
    warning('discontinuity in strand')
    
    is_problem_strand = strand_max_delta_length > 10 ;
    
    strand_subscripts{ is_problem_strand }
    strand_max_delta_length( is_problem_strand )
    problem_strands = find( is_problem_strand );
    
    disp( 'breaking strands at discontinuity' )
    
    number_of_problem_strands = 0 ;
    
    length( strand_subscripts );
    
    for problem_strand = problem_strands'
       
        number_of_problem_strands = number_of_problem_strands + 1 ;
        
        [ maxval, argmaxval ] = max( max( abs(   strand_subscripts{ problem_strand }( 2 : end    , 1 : 3 )          ...
                                               - strand_subscripts{ problem_strand }( 1 : end - 1, 1 : 3 )), [ ], 2 ), [ ], 1 );
                                           
        solution_strand = number_of_strands + number_of_problem_strands ;
                                           
        strand_subscripts{ solution_strand } = strand_subscripts{ problem_strand }(       1 : argmaxval, : );
        strand_subscripts{  problem_strand } = strand_subscripts{ problem_strand }( argmaxval + 1 : end, : );
        
        strand_energies{ solution_strand } = strand_energies{ problem_strand }(       1 : argmaxval );
        strand_energies{  problem_strand } = strand_energies{ problem_strand }( argmaxval + 1 : end );
        
        strands2vertices( solution_strand, 2 ) = strands2vertices(  problem_strand, 2 );       
        strands2vertices(  problem_strand, 2 ) = number_of_vertices + number_of_problem_strands * 2 - 1 ;
        strands2vertices( solution_strand, 1 ) = number_of_vertices + number_of_problem_strands * 2     ;
        
        vertex_energies         = [ vertex_energies         ; strand_energies{    problem_strand }( end        )];
        vertex_energies         = [ vertex_energies         ; strand_energies{   solution_strand }(  1         )];
        vertex_scale_subscripts = [ vertex_scale_subscripts ; strand_subscripts{  problem_strand }( end,   4   )];
        vertex_scale_subscripts = [ vertex_scale_subscripts ; strand_subscripts{ solution_strand }(  1 ,   4   )];
        vertex_space_subscripts = [ vertex_space_subscripts ; strand_subscripts{  problem_strand }( end, 1 : 3 )];
        vertex_space_subscripts = [ vertex_space_subscripts ; strand_subscripts{ solution_strand }(  1 , 1 : 3 )];
        
    end
end

%% continue strand extraction

strand_space_subscripts = cellfun( @( x ) x( :, 1 : 3 ), strand_subscripts, 'UniformOutput', false );
strand_scale_subscripts = cellfun( @( x ) x( :,   4   ), strand_subscripts, 'UniformOutput', false );

[ vessel_directions ] = get_vessel_directions_V3( strand_space_subscripts, microns_per_voxel );

strand_subscripts = cellfun( @( x, y ) [ x, y ], strand_space_subscripts, strand_scale_subscripts, 'UniformOutput', false );

mean_strand_energies = get_edge_metric( strand_energies );

network_statistics = calculate_network_statistics( strand_subscripts, bifurcation_vertices, lumen_radius_in_microns_range, microns_per_voxel, size_of_image );

%% check for vertices not matching with STRAND endpoints in 3 space
endpoint_locations  = cell2mat([ cellfun( @( x ) x(  1 , 1 : 3 ), strand_subscripts, 'UniformOutput', false ); ...
                                 cellfun( @( x ) x( end, 1 : 3 ), strand_subscripts, 'UniformOutput', false )  ]);

distance_from_vertex = max( abs( endpoint_locations - double( vertex_space_subscripts( strands2vertices( : ), : ))), [ ], 2 );
                             
if  any( distance_from_vertex > 2 )
    
    warning('strands2vertices vertex physically away from strand end')
    
    distance_from_vertex( distance_from_vertex > 2 )
    numel(find( distance_from_vertex > 2 ))
    
    % !!!!!!!!!! found a problem: the vertices identified in the strands2vertices LUT are not
    % colocated witht the ends of the strands:
    histogram( distance_from_vertex )
    
    % compare to random pairing of vertices identified in the strands2vertices LUT:
    figure
    distance_from_vertex = max( abs( double( vertex_space_subscripts( strands2vertices( randperm( numel( strands2vertices ))), : )) - double( vertex_space_subscripts( strands2vertices( : ), : ))), [ ], 2 );
    histogram( distance_from_vertex )
    
end

%% overwriting
delete( path_to_network )

save(   path_to_network, ...
              'bifurcation_vertices'       , ...
              'strand_subscripts'          , ...
              'strand_energies'            , ...
              'mean_strand_energies'       , ...
              'vessel_directions'          , ...
              'network_runtime_in_seconds' , ...
              'network_statistics'         , ...
              'strands2vertices'             );

delete( path_to_curated_edges )

save(   path_to_curated_edges, ...
                'edge_curation', ...
                'edge_curation_runtime_in_seconds', ...
                'edge_energies', ...
                'edge_scale_subscripts', ...
                'edge_space_subscripts', ...
                'edges2vertices', ...
                'mean_edge_energies', ...
                'vertex_energies',                    ...
                'vertex_scale_subscripts',            ... 
                'vertex_space_subscripts'             );

end