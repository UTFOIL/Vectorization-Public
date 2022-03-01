function blurred_image = gaussian_blur_in_chunks( original_image, pixels_per_sigma_gaussian, max_voxels_per_node )
%% Gaussian blur at many scales, chunk-wise and in parallel
% SAM 12/8/19

size_of_image = size( original_image );

% original_handle = 'original_input_to_gaussian_blur_in_chunks' ;
%  blurred_handle =  'blurred_output_of_gaussian_blur_in_chunks' ;

% original_file = [ data_directory, original_handle ];  
%  blurred_file = [ data_directory,  blurred_handle ];

original_file = 'original_input_to_gaussian_blur_in_chunks' ;
 blurred_file =  'blurred_output_of_gaussian_blur_in_chunks' ;
 
chunk_directory = [ blurred_file, '_chunks\' ];

try rmdir( chunk_directory, 's' ), catch, end

try delete( blurred_file  ), catch, end

try delete( original_file ), catch, end

h5create( original_file, '/d', size_of_image )

mat2h5( original_file, original_image )

clear( 'original_image' )

mkdir( chunk_directory );

% pixels_per_radius_range         = lumen_radius_in_microns_range    ./ microns_per_voxel ;

[ chunk_lattice_dimensions, number_of_chunks ] ...
         = get_chunking_lattice_V190( pixels_per_sigma_gaussian, max_voxels_per_node, size_of_image );     
     
chunk_index_range = 1 : number_of_chunks ;     
     
[ ~, ~, ~, ~, ~, ~,                                                                         ...
  y_writing_starts, x_writing_starts, z_writing_starts,                                     ...
  y_writing_counts, x_writing_counts, z_writing_counts   ]                                  ...
        = get_starts_and_counts_V200( chunk_lattice_dimensions, [ 0, 0, 0 ], size_of_image, [ 1, 1, 1 ]);
                                                                            
%% Main chunk PARFOR: 
% Reads a chunk from the interpolated file that is larger than the volume of interest. Does
% gaussian filtering at many scales.  Writes a smaller chunk to the blurred.
parfor chunk_index = chunk_index_range

    [ y_chunk_index, x_chunk_index, z_chunk_index ] = ind2sub( chunk_lattice_dimensions, chunk_index );
    
    writing_counts  = [ y_writing_counts( y_chunk_index ), ...
                        x_writing_counts( x_chunk_index ), ...
                        z_writing_counts( z_chunk_index )  ];
    
    chunk_overlap_vector = 6 * pixels_per_sigma_gaussian ;

    [ y_reading_starts, x_reading_starts, z_reading_starts,         ...
      y_reading_counts, x_reading_counts, z_reading_counts,         ...
      ~, ~, ~, ~, ~, ~,                                             ...
      y_offset,         x_offset,         z_offset,          ]      ...
            = get_starts_and_counts_V200( chunk_lattice_dimensions, ...
                                          chunk_overlap_vector,     ...
                                          size_of_image,            ...
                                          [ 1, 1, 1 ]               );

    reading_counts  = [ y_reading_counts( y_chunk_index ), ...
                        x_reading_counts( x_chunk_index ), ...
                        z_reading_counts( z_chunk_index )  ];

    % Read original image chunk and interpolate it to downsample at higher octaves
    original_chunk   = h52mat( original_file,                         ...
                               [ y_reading_starts( y_chunk_index ),   ...
                                 x_reading_starts( x_chunk_index ),   ...
                                 z_reading_starts( z_chunk_index ) ], ...
                               reading_counts,                        ...
                               [ 1, 1, 1 ]                            ); 

%         chunk_dft = fourier_transform_V2( original_chunk );
    blurred_chunk = gaussian_blur( original_chunk, pixels_per_sigma_gaussian );

    % extract a local range that is guaranteed to include the range to be written
    y_local_range = y_offset( y_chunk_index ) + 1                               ...
                  : y_offset( y_chunk_index ) + y_writing_counts( y_chunk_index );
              
    x_local_range = x_offset( x_chunk_index ) + 1                               ...
                  : x_offset( x_chunk_index ) + x_writing_counts( x_chunk_index );
              
    z_local_range = z_offset( z_chunk_index ) + 1                               ...
                  : z_offset( z_chunk_index ) + z_writing_counts( z_chunk_index );

    blurred_chunk = blurred_chunk( y_local_range, x_local_range, z_local_range );
    
    blurred_chunk_name = [ int2str( chunk_index ), ' of ', int2str( number_of_chunks )];
    
    blurred_chunk_file = [ chunk_directory, blurred_chunk_name ];

    h5create( blurred_chunk_file, '/d', writing_counts )
    
    % Write h5 files to individual chunk files (to be combined later)   
    mat2h5( blurred_chunk_file, blurred_chunk );
        
end % chunk PARFOR

%% make tiled master files
% combine the chunk files into a master file outside of a parfor to avoid simulataneous writing
% issues

h5create( blurred_file, '/d', size_of_image )

for chunk_index = chunk_index_range

    [ y_chunk_index, ...
      x_chunk_index, ...
      z_chunk_index ] = ind2sub( chunk_lattice_dimensions, chunk_index );

    blurred_chunk_name = [ int2str( chunk_index ), ' of ', int2str( number_of_chunks )];

    writing_starts = [ y_writing_starts( y_chunk_index ), ...
                       x_writing_starts( x_chunk_index ), ...
                       z_writing_starts( z_chunk_index )  ];

    writing_counts = [ y_writing_counts( y_chunk_index ), ...
                       x_writing_counts( x_chunk_index ), ...
                       z_writing_counts( z_chunk_index )  ];                        

    blurred_chunk_file = [ chunk_directory, blurred_chunk_name ];

    mat2h5( blurred_file, h52mat( blurred_chunk_file ), writing_starts, writing_counts );
    
end % chunk FOR

%% output
blurred_image = h52mat( blurred_file );

%% Remove h52mat files that were used for data processing

try rmdir( chunk_directory, 's' ), catch, end

try delete( blurred_file  ), catch, end

try delete( original_file ), catch, end

end % FUNCTION