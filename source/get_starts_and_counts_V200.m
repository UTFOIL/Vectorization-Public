function  [ y_reading_starts, x_reading_starts, z_reading_starts,      ...
            y_reading_counts, x_reading_counts, z_reading_counts,      ...
            y_writing_starts, x_writing_starts, z_writing_starts,      ...
            y_writing_counts, x_writing_counts, z_writing_counts,      ...
            y_offsets,        x_offsets,        z_offsets         ]    ...
                    = get_starts_and_counts_V200( chunk_lattice_dimensions, chunk_overlap_in_pixels, size_of_image, resolution_factors )
%% SAM, December of 2017
%
%
% V140 in which the chunk overlap is a vector quantity to handle the anisotropic case SAM 5/7/18
%
% V200 in which a resolution factor is input to  determine how much downsampling will occur. SAM
% 12/10/18
chunk_overlap_in_pixels = uint16( chunk_overlap_in_pixels );

y_writing_borders = uint16( linspace( 0, double( size_of_image( 1 )), double( chunk_lattice_dimensions( 1 ) + 1 )));
x_writing_borders = uint16( linspace( 0, double( size_of_image( 2 )), double( chunk_lattice_dimensions( 2 ) + 1 )));
z_writing_borders = uint16( linspace( 0, double( size_of_image( 3 )), double( chunk_lattice_dimensions( 3 ) + 1 )));

y_reading_starts  = ( y_writing_borders( 1 : end - 1 ) - chunk_overlap_in_pixels( 1 )) + 1 ;
x_reading_starts  = ( x_writing_borders( 1 : end - 1 ) - chunk_overlap_in_pixels( 2 )) + 1 ;
z_reading_starts  = ( z_writing_borders( 1 : end - 1 ) - chunk_overlap_in_pixels( 3 )) + 1 ;

size_of_image = uint16( size_of_image );

y_reading_ends    = min( y_writing_borders( 2 : end ) + chunk_overlap_in_pixels( 1 ), size_of_image( 1 ));
x_reading_ends    = min( x_writing_borders( 2 : end ) + chunk_overlap_in_pixels( 2 ), size_of_image( 2 ));
z_reading_ends    = min( z_writing_borders( 2 : end ) + chunk_overlap_in_pixels( 3 ), size_of_image( 3 ));

y_reading_counts  = ( y_reading_ends - y_reading_starts ) + 1 ;
x_reading_counts  = ( x_reading_ends - x_reading_starts ) + 1 ;
z_reading_counts  = ( z_reading_ends - z_reading_starts ) + 1 ;

% for the last chunks in each dimension ensure that the total distance (in pixels) covered by the
% reading counts is divisible by the resolution factor in that dimension (so it lands on the last
% pixel of the boundary of the image for the far edge)
y_reading_counts( end ) = 1 + resolution_factors( 1 ) * floor(( double( y_reading_counts( end ) - 1 )) / resolution_factors( 1 ));
x_reading_counts( end ) = 1 + resolution_factors( 2 ) * floor(( double( x_reading_counts( end ) - 1 )) / resolution_factors( 2 ));
z_reading_counts( end ) = 1 + resolution_factors( 3 ) * floor(( double( z_reading_counts( end ) - 1 )) / resolution_factors( 3 ));

% ensure that the chunks that will write to the end of the image in each dimension will read the
% last pixel in that dimension
y_reading_starts( end ) = y_reading_ends( end ) - y_reading_counts( end ) + 1 ;
x_reading_starts( end ) = x_reading_ends( end ) - x_reading_counts( end ) + 1 ;
z_reading_starts( end ) = z_reading_ends( end ) - z_reading_counts( end ) + 1 ;

y_writing_starts  = y_writing_borders( 1 : end - 1 ) + 1 ;
x_writing_starts  = x_writing_borders( 1 : end - 1 ) + 1 ;
z_writing_starts  = z_writing_borders( 1 : end - 1 ) + 1 ;

y_writing_ends    = y_writing_borders( 2 : end ) ;
x_writing_ends    = x_writing_borders( 2 : end ) ;
z_writing_ends    = z_writing_borders( 2 : end ) ;

y_writing_counts  =  y_writing_ends - y_writing_starts + 1 ;
x_writing_counts  =  x_writing_ends - x_writing_starts + 1 ;
z_writing_counts  =  z_writing_ends - z_writing_starts + 1 ;

y_offsets =  y_writing_starts - y_reading_starts ;
x_offsets =  x_writing_starts - x_reading_starts ;
z_offsets =  z_writing_starts - z_reading_starts ;

% convert to double because the h5 file indexing requires the indeces to be double
y_reading_starts = double( y_reading_starts );
x_reading_starts = double( x_reading_starts );
z_reading_starts = double( z_reading_starts );

y_reading_counts = double( y_reading_counts );
x_reading_counts = double( x_reading_counts );
z_reading_counts = double( z_reading_counts );

y_writing_starts = double( y_writing_starts );
x_writing_starts = double( x_writing_starts );
z_writing_starts = double( z_writing_starts );

y_writing_counts = double( y_writing_counts );
x_writing_counts = double( x_writing_counts );
z_writing_counts = double( z_writing_counts );

y_offsets = double( y_offsets );
x_offsets = double( x_offsets );
z_offsets = double( z_offsets );

end