% fix_strand_vertex_mismatch_again SAM 9/18/21

% the problem found in august 2020 persists in some weird way (see for backstory: %%
% fix_strand_vertex_mismatch 200811 (August 11th, 2020) SAM):s

is_good_for_bad = true ; % keep only the "bad" strands if TRUE

% path_to_network = 'E:\2P imaging\2021_Chronic_Imaging\Dan week 1\fused\batch_210902-175057\vectors\network_210915-224156_Dan_Fused_Raw_w1.mat' ;
% path_to_network = 'E:\2P imaging\2021_Chronic_Imaging\Doug week 1\batch_210626-184625\vectors\network_210901-131109_Doug_Fused_Raw_w1.mat' ;
% path_to_network = 'E:\2P imaging\2021_Chronic_Imaging\Eddy week 1\batch_210604-203819\vectors\network_210604-203819_Eddy_Fused_Raw_w1_log_stable.mat' ;
path_to_network         = 'E:\Annie\200923 2x2 vasculature RG med filter\200923 2x2 vasculature RG med filter\batch_200925-184104\vectors\network_210328-032210_Fused_medfilt_nobg.mat' ;

path_to_network_backup = [ path_to_network( 1 : end - 4 ), '_backup.mat' ];

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

strand_delta_lengths_max = cellfun( @( x ) max( max( abs( x( 2 : end, 1 : 3 )- x( 1 : end - 1, 1 : 3 )), [ ], 2 )), strand_subscripts );

% thresholoding
is_strand_bad = strand_delta_lengths_max > 1 ;

figure, histogram( log( strand_delta_lengths_max( is_strand_bad ) - 1 ))

% is_strand_bad = log( strand_delta_lengths_max - 1 ) > -15 ;
is_strand_bad = strand_delta_lengths_max > 1 + exp( -15 );

% diagnosics
figure, histogram( network_statistics.strand_lengths(          is_strand_bad ))
figure, histogram( network_statistics.strand_tortuosity(       is_strand_bad ))
figure, histogram( network_statistics.strand_ave_radii(        is_strand_bad ))
figure, histogram( network_statistics.strand_ave_directions(   is_strand_bad ))
figure, histogram( network_statistics.strand_areas(            is_strand_bad ))
figure, histogram( network_statistics.strand_volumes(          is_strand_bad ))
figure, histogram( network_statistics.strand_z_direction(      is_strand_bad ))

if is_good_for_bad, is_strand_bad = ~ is_strand_bad ; end

% !!!!!!!!!! recorde the bad vertices associated to the bad strands before erasing them in this line
% remove those bad vertices from the bifurcation vertices if they do not appear with other strands!!
strands2vertices( is_strand_bad, : ) = [ ];

mean_strand_energies( is_strand_bad ) = [ ];
     strand_energies( is_strand_bad ) = [ ];
   strand_subscripts( is_strand_bad ) = [ ];
   vessel_directions( is_strand_bad ) = [ ];

% !!!!!!!!!!!!!!!! recalculate strand statistics for accuracy of the bulk network statistics
% !!!!!!!!! calculate_strand_statistics()
network_statistics.strand_lengths(          is_strand_bad ) = [ ];
network_statistics.strand_tortuosity(       is_strand_bad ) = [ ];
network_statistics.strand_ave_radii(        is_strand_bad ) = [ ];
network_statistics.strand_ave_directions(   is_strand_bad ) = [ ];
network_statistics.strand_areas(            is_strand_bad ) = [ ];
network_statistics.strand_volumes(          is_strand_bad ) = [ ];
network_statistics.strand_z_direction(      is_strand_bad ) = [ ];

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