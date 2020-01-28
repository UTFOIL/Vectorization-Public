function [ radius_centerline_mask ] = centerlines_from_morphology( vessel_mask, scale_space_sigma, visual, root_directory )
% This function does repeated morphological close and erode to the binary mask of vasculature to
% obtain a mask of centerlines.  The function works like this:  First it erodes, then it dilates the
% mask (this is a close operation yiedling a "closed mask").  Then it takes the original mask and
% subtracts the closed one from it.  Whatever was in the original but not in the closed one, must be
% a centerline for a vessel of width equal to the diameter of the structuring element.  The
% algorithm then starts from the mask that has been eroded once and repeats the process until all of
% the vessels have been eroded away.  
%
% SM 17/10/17

erosion_structuring_element_radius = 3 * scale_space_sigma( 1 );

% [ erosion_structuring_element_subscripts, erosion_structuring_element_linear_indexing ]   = construct_structuring_element(          scale_space_sigma, size( vessel_mask ));
[ erosion_structuring_element_subscripts, erosion_structuring_element_linear_indexing ]   = construct_structuring_element(          scale_space_sigma, size( vessel_mask ));


dilation_scale_space_sigma = 2 * scale_space_sigma([ 1, 1, 1 ]);

[ dilation_structuring_element_subscripts, dilation_structuring_element_linear_indexing ] = construct_structuring_element( dilation_scale_space_sigma, size( vessel_mask ));

all_vessels_eroded_flag = false;

vessel_mask = logical( vessel_mask );

% centerline_mask            = zeros( size( vessel_mask ));
radius_centerline_mask     = zeros( size( vessel_mask ));

number_of_erosions = 0;

centerline_directory = [ 'centerlines_', char( datetime('now', 'TimeZone', 'local', 'Format', 'yyMMdd-HHmmss' )), '\' ];

mkdir( root_directory, centerline_directory )

while ~all_vessels_eroded_flag
    
    eroded_vessel_mask = erosion_filter(  vessel_mask,  erosion_structuring_element_subscripts,  erosion_structuring_element_linear_indexing );
    
    closed_vessel_mask = dilation_filter( eroded_vessel_mask, dilation_structuring_element_subscripts, dilation_structuring_element_linear_indexing );
    
    difference_mask = vessel_mask - closed_vessel_mask;
    
    nonnegative_difference_mask = difference_mask;
    
    nonnegative_difference_mask( nonnegative_difference_mask < 0 ) = 0;
    
%     centerline_mask = centerline_mask + nonnegative_difference_mask;
        
    number_of_erosions = number_of_erosions + 1;
    
    radius_centerline_mask = radius_centerline_mask + nonnegative_difference_mask * number_of_erosions * erosion_structuring_element_radius;
    
    if visual
        
%         mat2tif( vessel_mask, [ root_directory, centerline_directory, int2str( number_of_erosions ), '_erosion_', 'vessel_mask.tif' ])
%         mat2tif( eroded_vessel_mask, [ root_directory, centerline_directory, int2str( number_of_erosions ), '_erosion_', 'eroded_vessel_mask.tif' ])
%         mat2tif( closed_vessel_mask, [ root_directory, centerline_directory, int2str( number_of_erosions ), '_erosion_', 'closed_vessel_mask.tif' ])
%         mat2tif( nonnegative_difference_mask, [ root_directory, centerline_directory, int2str( number_of_erosions ), '_erosion_', 'nonnegative_difference_mask.tif' ])
        mat2tif( difference_mask, [ root_directory, centerline_directory, int2str( number_of_erosions ), '_erosion_', 'difference_mask.tif' ])
%         mat2tif( centerline_mask, [ root_directory, centerline_directory, int2str( number_of_erosions ), '_erosion_', 'centerline_mask.tif' ])
        mat2tif( radius_centerline_mask, [ root_directory, centerline_directory, int2str( number_of_erosions ), '_erosion_', 'radius_mask.tif' ])

    end
    
    vessel_mask = eroded_vessel_mask;
    
    all_vessels_eroded_flag = all( vessel_mask( : ) == 0 );
end