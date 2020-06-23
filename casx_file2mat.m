function [ point_coordinates, arc_connectivity, arc_diameters ] = casx_file2mat( casX_file_path )
% casX2mat SAM 7/22/19 
%  This function reads the casX text file to extract the three output variables as matlab variables
%  (arrays of doubles)

fileID = fopen( casX_file_path, 'r' );

number_of_points = find_number_after_literal( fileID, 'nPoints=' );

formatSpec_point_coordinates = '%f %f %f \n' ;
formatSpec_connectivities    =    '%x %x \n' ;
formatSpec_diameters         =       '%f \n' ;

point_coordinates = fscanf( fileID, formatSpec_point_coordinates, [ 3, number_of_points ])';

number_of_arcs = find_number_after_literal( fileID, 'nArcs=' );

arc_connectivity  = fscanf( fileID, formatSpec_connectivities,    [ 2, number_of_arcs   ])';

                 find_number_after_literal( fileID, 'nArcs=' ); % just to move the file pointer

arc_diameters     = fscanf( fileID, formatSpec_diameters,         [ 1, number_of_arcs   ])';

fclose( fileID );

end % FUNCTION casX2mat


function number = find_number_after_literal( fileID, literal )

is_searching_for_literal = true ;

while is_searching_for_literal
    
    tline = fgetl( fileID );
    
    matches = strfind( tline, literal );
    
    is_searching_for_literal = isempty( matches ) & ischar( tline );
    
end % WHILE is_searching_for_literal

number = str2double( tline( matches + length( literal ) : end ));

end % FUNCTION find_number_after_literal

% function y = litcount(filename, literal)
% % Count the number of times a given literal appears in each line.
% 
% fid = fopen(filename);
% y = 0;
% tline = fgetl(fid);
% while ischar(tline)
%    matches = strfind(tline, literal);
%    num = length(matches);
%    if num > 0
%       y = y + num;
%       fprintf(1,'%d:%s\n',num,tline);
%    end
%    tline = fgetl(fid);
% end
% fclose(fid);
% 
% end % FUNCTION