function [ number, string ] = find_number_after_literal( fileID, literal )
% string to doubles the end of the line after the first matching literal SAM 4/22/21

is_searching_for_literal = true ;

while is_searching_for_literal
    
    tline = fgetl( fileID );
    
    matches = strfind( tline, literal );
    
    is_searching_for_literal = isempty( matches ) & ischar( tline );
    
end % WHILE is_searching_for_literal

string = tline( matches + length( literal ) : end );

number = str2double( string );

end % FUNCTION find_number_after_literal