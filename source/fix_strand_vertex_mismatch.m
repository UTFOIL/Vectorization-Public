function [ strands2vertices ] = fix_strand_vertex_mismatch( vertex_space_subscripts, strand_subscripts, strands2vertices )
%% fix_strand_vertex_mismatch 200811 (August 11th, 2020) SAM
% This function is a patch on the vector set used in the SLAVV methods paper. The strands2vertices
% function had a bug affecting 56 strands, in which the second entry in strands2vertices pointed to
% a vertex which was not physically collocated with the end of the strand. This function replaces
% these 58 vertices with the 58 unique, closest vertices to the ends of the strands. Downstream
% effects of this patch are that data conversion to .casx and .vmv file types will be free of the
% long straight tubular artifacts crossing the entirety of the image. After running this patch,
% however two strands remained who did not perfectly coincide with a vertex. Now the conversion
% functions to casx and vmv files are robust to the strand vertex mismatch and its artifacts. the
% vertex positions are no longer inputs! artifacts may still exist in the bifurcation vertices
% visualization tif and perhaps the bifurcations themselves.
%
% SAM 200811

diffs = zeros( 2, 3, length( strand_subscripts ));
for idx = 1 : length( strand_subscripts )
diffs( :, :, idx ) = abs( strand_subscripts{ idx }([ 1, end ], 1 : 3 ) - double( vertex_space_subscripts( strands2vertices( idx, [ 1, 2 ]), : )));
end

strands_mismatched_to_vertices = find( any( any( diffs, 2 ), 1 ))';

new_mismatch_error = zeros( size( strands2vertices, 1 ), 1 );

for idx = strands_mismatched_to_vertices
    
    [ new_mismatch_error( idx ), strands2vertices( idx, 2 )] = min( sum( abs( strand_subscripts{ idx }( end, 1 : 3 ) - double( vertex_space_subscripts )), 2 ));
    
end
