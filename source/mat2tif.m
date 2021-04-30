function  mat2tif( image_matrix, tif_path )
% much of this taken from Process_Data_David_03092017
%
% SAM 9/6/2017

[ tagstruct.ImageLength, tagstruct.ImageWidth, number_of_slices ] = size( image_matrix );

input_data_type = class( image_matrix );

digit_indices = regexp( input_data_type, '\d' );        

if isempty( digit_indices ) % double input goes to int16 output
    
    bits_per_sample = 16 ;
    
    image_matrix = int16( image_matrix );
    
    tagstruct.SampleFormat = Tiff.SampleFormat.Int;
    
else % assume unsigned integer input !!!!!
    
    bits_per_sample = sscanf( input_data_type( digit_indices( 1 ) : end ), '%f' );

end

tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = bits_per_sample ;
tagstruct.SamplesPerPixel = 1; %added 11/1
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB'; %added 11/1

tRaw = Tiff( tif_path, 'w' ); %Create tif file

tRaw.setTag( tagstruct ); %instantiate parameters

% temp_raw = int16( image_matrix( :, :, 1 )); %create 16 bit signed values for image slice i of raw data
temp_raw = image_matrix( :, :, 1 ); 

tRaw.write( temp_raw ); 

tRaw.close( );

for slice_index = 2 : number_of_slices
    
%     temp_raw = int16( image_matrix( :, :, slice_index)); %create 16 bit signed values for image slice i of raw data
    temp_raw = image_matrix( :, :, slice_index );

    tRaw = Tiff( tif_path, 'a' ); %Create tif file

    tRaw.setTag( tagstruct ); %instantiate parameters

    tRaw.write( temp_raw ); 
end

tRaw.close( ); %close tif file

end % FUNCTION

