function casx_mat2file( casx_file_path, point_coordinates, arc_connectivity, arc_diameters )
%  casx_mat2file SAM 200812 (August 12th, 2020)
%  This function writes the three casX output vasriables to a text file in the .casx format.

if isfile( casx_file_path ), delete( casx_file_path ), end

number_of_points     =   size( point_coordinates, 1 );
number_of_arcs       = length( arc_diameters        );

fileID = fopen( casx_file_path, 'a' );

%% writing vmv header
affiliation = 'put_your_affiliation_here' ;
location    = 'put_your_location_here'    ;
% affiliation = 'SMihelic at FOIL' ;
% location    = 'Austin, TX'    ;
date        = char( datetime('now', 'TimeZone', 'local', 'Format', 'MM/dd/yyyy' )); 

fprintf( fileID, ...
         [ '//"3d vascular network generation" author="A.Linninger" date="2007-2011"\n',     ...
           '//File format designed by GHartung and ALinninger 10/9/2018\n',                  ... 
           '//This file was created by ', affiliation, ' in ', location, ' on: ', date, '\n' ]);

%% data formatting
% formatSpec_point_coordinate = '\t%.9E\t%.9E\t%.9E\n' ;
formatSpec_arc_connectivity =   '%x\t%x\t\n' ;
formatSpec_arc_diameter     =         '%g\n' ;       
% formatSpec_groupId          =         '%u\n' ;       
       
%% writing point values table

fprintf( fileID, '\n//point coordinates;   nPoints=%u\n', number_of_points );

for point_idx = 1 : number_of_points
    
    pow = floor(log10(abs(point_coordinates( point_idx, : ))));
    s = sprintf('\t%.9fE%+.3d', [point_coordinates( point_idx, : )./10.^pow; pow]);
    
    fprintf( fileID, [ s, '\n' ]);

end

fprintf( fileID, '//end point coordinates\n' );

%% writing connectivity table
fprintf( fileID, '\n//arc connectivity matrix;   nArcs=%u\n', number_of_arcs );

for arc_idx = 1 : number_of_arcs
   
    fprintf( fileID, formatSpec_arc_connectivity, arc_connectivity( arc_idx, : ));
    
end

fprintf( fileID, '//end arc connectivity matrix\n' );


%% writing diameters list
fprintf( fileID, '\n//diameter: vector on arc;   nArcs=%u\n', number_of_arcs );

for arc_idx = 1 : number_of_arcs
   
    fprintf( fileID, formatSpec_arc_diameter, arc_diameters( arc_idx ));
    
end

fprintf( fileID, '//end diameter\n' );

%% writing groupId list
% % !!!! use varargin for this input?

% fprintf( fileID, '\n//groupId: vector on arc;   nArcs=%u\n', number_of_arcs );
% 
% for arc_idx = 1 : number_of_arcs
%    
%     fprintf( fileID, formatSpec_groupId, groupId( arc_idx ));
%     
% end
% 
% fprintf( fileID, '//end groupId\n' );


%% close file
fclose( fileID );

end % FUNCTION casx_file2mat


