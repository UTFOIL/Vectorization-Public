function [ vessel_directions ] = get_vessel_directions_V3( strand_space_subscripts, microns_per_voxel )
%% get_vessel_directions SAM for Chakameh 7/12/18
% The purpose of this function is to loop through the strand/junction objects putting the edges
% together into strands or junctions as encoded by the inputs edge_indices_in_strands and
% edge_backwards_in_strands.  Then for each strand doing some spatial and radial smoothing and then
% computing the spatial derivative with respect the axial direction of the strand.  The output is
% then put back in the organization of the edge_subscripts input (a cell array of edges).  Note that
% the length units of the different components of the vessel direction vectors may not be the same
% due to anisotropic sampling in y, x, and z.
%
% V2, the sigma of the gaussian kernel for smoothing is an input % SAM 7/16/18
%
% V2 also, this version was merged with V4 which includes the tissue type cutoffs input and
% tissue_types output. SAM 3/5/19
%
% V3, removed the strand assembly step and put it in the function GET_STRAND_OBJECTS.  Removed the
% strand smoothing step and instead use the SMOOTH_EDGES function.  Now expecting the strand objects
% to come in as a cell array, and will output a cell array.  SAM 4/23/19

number_of_strands = numel( strand_space_subscripts );

strand_index_range = 1 : number_of_strands ;

vessel_directions = cell( size( strand_space_subscripts ));

is_symmetric_difference_good = cellfun( @( x ) size( x, 1 ) > 2, strand_space_subscripts );

for strand_index = strand_index_range

    strand_subscripts_at_strand = strand_space_subscripts{ strand_index };
    
    % Spatial Derivative Approximation:
                            
    % convert to microns in each spatial dimension before computing real-space directions
    strand_subscripts_at_strand = strand_subscripts_at_strand .* microns_per_voxel ;    
    
    if is_symmetric_difference_good( strand_index )

        % approximate the spatial derivative with respect to the axial direction of the vessel using a
        % symmetric difference.
        vessel_directions_at_strand_cropped                                                  ...
                                         = strand_subscripts_at_strand( 1 + 2 : end    , : ) ...
                                         - strand_subscripts_at_strand( 1     : end - 2, : ) ;
                                     
        % Pad the ends of the cropped output with copies of the first and last entries
        vessel_directions_at_strand = [ vessel_directions_at_strand_cropped( 1,   : ); ...
                                        vessel_directions_at_strand_cropped          ; ...
                                        vessel_directions_at_strand_cropped( end, : )  ];                                     
                                     
    else % Else aysmetric difference, for when there are only two edge locations
       
        vessel_directions_at_strand_cropped                                     ...
                                         = strand_subscripts_at_strand( 2 , : ) ...
                                         - strand_subscripts_at_strand( 1 , : ) ;

        % Pad only one end of the cropped output with a copy of the other entry
        vessel_directions_at_strand = [ vessel_directions_at_strand_cropped( 1, : ); ...
                                        vessel_directions_at_strand_cropped          ];                                     
        
    end

    % normalize the spatial derivatives with respect to the pythagorean sum (make them unit vectors)
    vessel_directions{ strand_index } =         vessel_directions_at_strand                   ...
                                      ./ ( sum( vessel_directions_at_strand .^ 2, 2 )) .^ 0.5 ;
                                 
end % strand FOR
end % FUNCTION