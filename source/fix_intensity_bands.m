% SAM 7/10/23

% input_file = 'E:\David\martinos\round 2\Angio_4X_Mouse2_Compiled_820um.tif'

function output_file = fix_intensity_bands(input_file)

image = tif2mat(input_file);

image_means = squeeze(mean(squeeze(mean(image,1)),1));

image = double(image) ./ reshape(image_means,[1,1,length(image_means)]) * mean(image_means);

figure
plot(image_means)

band_limits =   [ 7:   51;
                  52:  96;
                  97: 141;
                  142:186;
                  187:231;
                  232:276;  
                  277:321;
                  322:366;
                  367:410,410 ];

% band_limits(:,end)-band_limits(:,1)

bands = image_means(band_limits);

figure
plot((bands./bands(:,1))')

% avg_decay_factor = 0.7 ;

band_means = mean(bands,2);

figure
plot(band_means)


output_file = [ input_file(1:end-4), '_renormalized.tif' ];

mat2tif(uint16(image),output_file)

end
