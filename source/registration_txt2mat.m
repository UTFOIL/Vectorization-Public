function [ starts, dims, number_of_images ] = registration_txt2mat( path_to_registration_txt_file )

formatSpec_coordinates = '%f,%f,%f' ;
formatSpec_dimensions  = '%d,%d,%d' ;

fid = fopen( path_to_registration_txt_file, 'r' );

number_of_images = find_number_after_literal( fid, 'num =' );

starts = zeros( 4, 3 );
dims   = zeros( 4, 3 );

for image_idx = 1 : number_of_images

    [ ~, string ] = find_number_after_literal( fid, '(' ); % just to move the file pointer and get the rest of the line

    starts( image_idx, : ) = sscanf( string, formatSpec_coordinates,[ 3, 1 ])';

end

for image_idx = 1 : number_of_images

    [ ~, string ] = find_number_after_literal( fid, '(' ); % just to move the file pointer and get the rest of the line

    dims( image_idx, : ) = sscanf( string, formatSpec_coordinates,[ 3, 1 ])';

end

fclose(fid);

end