function [ image_statistics ] = calculate_image_statistics_from_binary(   signal_binary, ...
                                                                        original_image   )

% % SAM 2/11/22
    
    signal     = original_image(   signal_binary );
    background = original_image( ~ signal_binary );
    
    image_statistics.signal     = mean(     signal );
    image_statistics.background = mean( background );
    
    image_statistics.noise      = (   var(   signal   )        ...
                                    + var( background )) ^ 0.5 ;  
    
    image_statistics.contrast = image_statistics.signal     ...
                              - image_statistics.background ;
    
    image_statistics.CNR = image_statistics.contrast ...
                         / image_statistics.noise    ;
    
end