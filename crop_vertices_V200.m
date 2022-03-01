function [ space_subscripts, scale_subscripts, energy_values ]                               ...
            = crop_vertices_V200( space_subscripts, scale_subscripts, energy_values,         ...
                                  lumen_radius_in_microns_range, microns_per_pixel, size_of_image )
%% Choose vertices_V180
% the purpose of this function is to clean up the output of the get_vertices function. It does this
% by sorting the vertices from best to worst (highest to lowest) energy and then looping through and
% eliminating vertices whose volumes overlap with vertices of lower energy.
%
% SAM 5/29/18
%
% V181 in which the elimination is one sided: volume overlap isn't sufficient to eliminate vertices.
% We look to see if the current vertex's center (not volume) is inside the volume of any previously
% accepted vertex. SAM 5/29/18
%
% V182 is somewhere between V180 and V181 SAM 5/29/18, the inner sigma volumes need to overlap (not
% the sizes as in the case of V180.
%
% V183 in which the size of the volumes needed to overlap to exclude a vertex is input as a factor
% of sigma. SAM 5/29
%
% V184 where instead of searching for previously identified vertices that will interfere with the
% current object, the method is to look at a painted image of previously identified vertices where
% previous vertices have been painted. SAM 5/30/18
%
% V190 removing vertices that overlap the image boundary before entering the for loop. SAM 8/16/18
%
% function name changed from choose_vertices_V190 to crop_vertices_V190.  Keeping only the part of
% the choose vertices function that does the cropping. eliminating references to sigma; instead
% referring to lumen radius.  SAM 8/23/18
%
% V200 in which vertices of the largest and smallest scale coordinate are also removed.  SAM 12/4/18

radii_in_pixels =  uint16( round(    lumen_radius_in_microns_range( round( scale_subscripts )) ...
                                  ./                 microns_per_pixel                 ));

% predict which vertices will overlap the outside of the image
subscript_maxs = space_subscripts + radii_in_pixels ;
subscript_mins = space_subscripts - radii_in_pixels ;

y_is_over  = subscript_maxs( :, 1 ) > size_of_image( 1 );
x_is_over  = subscript_maxs( :, 2 ) > size_of_image( 2 );
z_is_over  = subscript_maxs( :, 3 ) > size_of_image( 3 );

y_is_under = subscript_mins( :, 1 ) <         1         ;
x_is_under = subscript_mins( :, 2 ) <         1         ;
z_is_under = subscript_mins( :, 3 ) <         1         ;

scale_is_max = round( scale_subscripts ) == length( lumen_radius_in_microns_range );
scale_is_min = round( scale_subscripts ) == 1 ;

% don't include vertices that touch the image boundary
excluded_vertices_logical =  y_is_over   |  x_is_over   | z_is_over  ...
                          |  y_is_under  |  x_is_under  | z_is_under ...
                          | scale_is_max | scale_is_min ;

space_subscripts( excluded_vertices_logical, : ) = [ ];
scale_subscripts( excluded_vertices_logical    ) = [ ];
   energy_values( excluded_vertices_logical    ) = [ ];
     
end % end cropping FUNCTION