function vmv_mat2file( vmv_file_path, point_coordinates, strand_points )
%  vmv_mat2file SAM 3/5/21
%  This function writes the a vmv text file from MATLAB variables

if isfile( vmv_file_path ), delete( vmv_file_path ), end

number_of_points       = size( point_coordinates, 1 );
number_of_point_values = size( point_coordinates, 2 ); % originally 3-space and radius

number_of_strands = length( strand_points );

number_of_points_for_strands = cellfun( @length, strand_points );

fileID = fopen( vmv_file_path, 'a' );

%% writing vmv header
fprintf( fileID, ...
         '$PARAM_BEGIN\nNUM_VERTS	%u\nNUM_STRANDS	%u\nNUM_ATTRIB_PER_VERT	%u\n$PARAM_END\n', ...
         [ number_of_points, number_of_strands, number_of_point_values ]);

%% writing point values table
formatSpec_point_coordinates = '%u\t' ;

for value_idx = 1 : number_of_point_values
            
    if value_idx < number_of_point_values
    
        formatSpec_point_coordinates = [ formatSpec_point_coordinates, '%f\t' ];
        
    else
        
        formatSpec_point_coordinates = [ formatSpec_point_coordinates, '%f\n'  ];
    
    end 
end

fprintf( fileID, '\n$VERT_LIST_BEGIN\n' );

for point_idx = 1 : number_of_points
   
    fprintf( fileID, formatSpec_point_coordinates, [ point_idx, point_coordinates( point_idx, : )]);
    
end

fprintf( fileID, '$VERT_LIST_END\n' );

%% writing strand indices list

fprintf( fileID, '\n$STRANDS_LIST_BEGIN\n' );

for strand_idx = 1 : number_of_strands
    
    formatSpec_strand_points = '%u\t' ; 
    
    for strand_point_idx = 1 : number_of_points_for_strands( strand_idx )
        
        if strand_point_idx < number_of_points_for_strands( strand_idx )
            
            formatSpec_strand_points = [ formatSpec_strand_points, '%u\t' ];
            
        else
            
            formatSpec_strand_points = [ formatSpec_strand_points, '%u\n' ];
            
        end        
    end
    
    fprintf( fileID, formatSpec_strand_points, [ strand_idx, strand_points{ strand_idx }]);
    
end

fprintf( fileID, '$STRANDS_LIST_END' );

fclose( fileID );

end % FUNCTION casx_file2mat


