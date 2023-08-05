% SAM 7/10/23
input_root = 'E:\David\martinos\round 2\batch_230629-110216\' ;
input_time_stamp = '230703-000627';
% !! functionalize this code % no more hard-coding required from here !!
load([ input_root, 'settings\energy_230629-110216.mat'                                 ])
load([ input_root, '\vectors\network_', input_time_stamp, '_Angio_4X_Mouse2_Compiled_820um.mat' ])
strand_subscripts_thresholded=strand_subscripts(network_statistics.strand_z_direction>0.5^0.5);
strand_position_radius = cellfun( @(x) [ x( :, [2,1,3]) .* microns_per_voxel([2,1,3]), lumen_radius_in_microns_range( floor(x(:,4))).^(1-(x(:,4)-floor(x(:,4)))).*lumen_radius_in_microns_range( ceil(x(:,4))).^((x(:,4)-floor(x(:,4))))], strand_subscripts_thresholded, 'UniformOutput', false );
position_radius = cell2mat(strand_position_radius);
xy_position_radius = position_radius( :, [1,2,4]);
save([ input_root, 'network_export_', input_time_stamp ], 'strand_position_radius', 'position_radius', 'xy_position_radius')
