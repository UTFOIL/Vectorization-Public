function  mat2tif( image_matrix, tif_path )
% much of this taken from Process_Data_David_03092017
%
% SAM 9/6/2017

[ tagstruct.ImageLength, tagstruct.ImageWidth, number_of_slices ] = size( image_matrix );

input_data_type = class( image_matrix );

digit_indices = regexp( input_data_type, '\d' );        

if isempty( digit_indices ) % double input goes to int16 output
    
%     tagstruct.SampleFormat = Tiff.SampleFormat.Int;    
    sample_format = 2 ; % Int, same as line above commented out
    
    bits_per_sample = 16 ;
        
    image_matrix = int16( image_matrix );
    
else
        
    bits_per_sample = sscanf( input_data_type( digit_indices( 1 ) : end ), '%f' );
    
    sample_format   =  input_data_type( 1 : digit_indices( 1 ) - 1 );
  
    switch sample_format
        
        case  'int', sample_format = 2 ; %  int % !!!!!! int8 not supported in imageJ
        case 'uint', sample_format = 1 ; % uint
            
    end     
end

tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
tagstruct.BitsPerSample = bits_per_sample ;
tagstruct.SampleFormat = sample_format ;
tagstruct.SamplesPerPixel = 1; %added 11/1 !!!!!!!!!!!!!!! addd color (mutliple channel) handling here ?????????????????????
tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
tagstruct.Software = 'MATLAB'; %added 11/1

tif_h = Tiff( tif_path, 'w' ); %Create tif file

tif_h.setTag( tagstruct ); %instantiate parameters

tif_h.write( image_matrix( :, :, 1 )); 

tif_h.close( );

for slice_index = 2 : number_of_slices
    
    tif_h = Tiff( tif_path, 'a' ); %Create tif file

    tif_h.setTag( tagstruct ); %instantiate parameters

    tif_h.write( image_matrix( :, :, slice_index )); 
end

tif_h.close( ); %close tif file

end % FUNCTION

