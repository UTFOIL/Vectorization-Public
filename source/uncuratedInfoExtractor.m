function featurePool = uncuratedInfoExtractor(object_string, location_of_vectors, path_to_original_data, path_to_energy_settings)
 
    load(location_of_vectors)
    load(path_to_energy_settings)

    if strcmp(object_string,'vertices')

        disp([newline 'Extracting Vertex Info']);
        tic

        featurePool.vertex_space_subscripts = vertex_space_subscripts;
        featurePool.vertex_scale_subscripts = vertex_scale_subscripts;       
        featurePool.vertex_energies = vertex_energies;       

        % end spatial info

        size_of_image = h5info(path_to_original_data);
        size_of_image = size_of_image.Datasets.Dataspace.Size;

        number_of_scales = length( lumen_radius_in_microns_range );

        scale_subscript_range = 1 : number_of_scales ;

        structuring_etlement_linear_indexing_templates = cell( number_of_scales, 1 );

        radii_in_pixels_range = lumen_radius_in_microns_range ./ microns_per_voxel ;

        for scale_subscript = scale_subscript_range

            % find all pixel locations within the ellipsoid radii from the vertex position    
            structuring_element_linear_indexing_templates{ scale_subscript }                                            ...
                = int64( construct_structuring_element( radii_in_pixels_range( scale_subscript, : ), size_of_image ));

        end % constructing relative elements FOR scale

        % init structuring elements
        space_subscripts_int64 = int64( vertex_space_subscripts );

        vertex_position_linear_indices =   space_subscripts_int64( :, 1 )                                               ...
                                       + ( space_subscripts_int64( :, 2 ) - 1 ) * size_of_image( 1 )                    ...
                                       + ( space_subscripts_int64( :, 3 ) - 1 ) * size_of_image( 1 ) * size_of_image( 2 );

        vertex_position_linear_indices_cell = num2cell( vertex_position_linear_indices );

        structuring_element_linear_indexing = structuring_element_linear_indexing_templates( round( vertex_scale_subscripts ))';

        vertex_structure_positions_linear_indexing = cellfun( @plus, structuring_element_linear_indexing,   ...
                                                                       vertex_position_linear_indices_cell, ...
                                                                                     'UniformOutput', false );   
        % end initialize_structuring_elements

        orig = double(h52mat(path_to_original_data));
        orig = reshape(orig,[],1);

        average_intensity_across_original =  mean(orig,'all');
         std_of_intensity_across_original = std(orig,0,'all');

        average_intensity_across_vertices = cellfun(@(x) mean(orig(x)), vertex_structure_positions_linear_indexing, 'UniformOutput', true);
         std_of_intensity_across_vertices = cellfun(@(x)  std(orig(x)), vertex_structure_positions_linear_indexing, 'UniformOutput', true);

        featurePool.normalized_mean_intensity_across_vertices   = average_intensity_across_vertices ./ average_intensity_across_original;
        featurePool.normalized_std_of_intensity_across_vertices =  std_of_intensity_across_vertices ./  std_of_intensity_across_original;

        fprintf('     '); toc

    end %vertex info

    if strcmp(object_string,'edges')

        disp([newline 'Extracting Edge Info']);
        tic

        mean_edge_space_subscripts = cellfun(@(x) mean(x,1),edge_space_subscripts,'UniformOutput', false);
        mean_edge_space_subscripts = cell2mat(mean_edge_space_subscripts);

        mean_edge_scale_subscripts = cellfun(@(x) mean(x,1),edge_scale_subscripts);
        std_edge_scale_subscripts = cellfun(@(x) std(x,1),edge_scale_subscripts);

        featurePool.largest_scale_subscript = cellfun(@(x) max(x(:)),edge_scale_subscripts);
        featurePool.edge_scale_residuals = cellfun(@(x) [x(1), x(end)],edge_scale_subscripts,'UniformOutput', false);
        featurePool.edge_scale_residuals = cell2mat(featurePool.edge_scale_residuals);
        featurePool.edge_scale_residuals = featurePool.largest_scale_subscript - featurePool.edge_scale_residuals;
        featurePool.edge_scale_residuals = sort(featurePool.edge_scale_residuals,2,'descend');

        featurePool.mean_edge_space_subscripts = mean_edge_space_subscripts;
        featurePool.mean_edge_scale_subscripts = mean_edge_scale_subscripts;
        featurePool.std_edge_scale_subscripts = std_edge_scale_subscripts;
        featurePool.mean_edge_energies = mean_edge_energies;

        % synthesize intenity info

        size_of_image = h5info(path_to_original_data);
        size_of_image = size_of_image.Datasets.Dataspace.Size;


        number_of_scales = length( lumen_radius_in_microns_range );

        scale_subscript_range = 1 : number_of_scales ;

        vertex_element_linear_indexing_templates = cell( number_of_scales, 1 );
          edge_element_linear_indexing_templates = cell( number_of_scales, 1 );

        vertex_radii_in_pixels_range =  1 ...
                                     .* lumen_radius_in_microns_range  ...
                                     ./ microns_per_voxel              ;

          edge_radii_in_pixels_range =  1   ...
                                     .* lumen_radius_in_microns_range ...
                                     ./ microns_per_voxel             ;

        for scale_subscript = scale_subscript_range

            % find all pixel locations within the ellipsoid radii from the element position
            vertex_element_linear_indexing_templates{ scale_subscript }                                                        ...
                = int64( construct_structuring_element( vertex_radii_in_pixels_range( scale_subscript, : ), size_of_image ));

              edge_element_linear_indexing_templates{ scale_subscript }                                                        ...
                = int64( construct_structuring_element(   edge_radii_in_pixels_range( scale_subscript, : ), size_of_image ));

        end % constructing edge and vertex element templates FOR scale


        degrees_of_edges_uint_32 = uint32( cellfun( @length, edge_scale_subscripts ));

        edge_space_subscripts_mat =         cell2mat( edge_space_subscripts );
        edge_scale_subscripts_mat = uint8(  cell2mat( edge_scale_subscripts ));        

        space_subscripts_int64 = int64( edge_space_subscripts_mat );

        position_linear_indices =   space_subscripts_int64( :, 1 )                                               ...
                                + ( space_subscripts_int64( :, 2 ) - 1 ) * size_of_image( 1 )                    ...
                                + ( space_subscripts_int64( :, 3 ) - 1 ) * size_of_image( 1 ) * size_of_image( 2 );

        sphere_position_linear_indices_cell = num2cell( position_linear_indices );

        structuring_element_linear_indexing = edge_element_linear_indexing_templates( edge_scale_subscripts_mat );

        sphere_structure_positions_linear_indexing = cellfun( @plus, structuring_element_linear_indexing, ...
                                                                     sphere_position_linear_indices_cell, ...
                                                                                  'UniformOutput', false  ); 

    %         vertex_structure_positions_linear_indexing = structuring_element_linear_indexing_templates( edge_scale_subscripts_mat );

        sphere_structure_positions_linear_indexing_edge_cell = mat2cell( sphere_structure_positions_linear_indexing,  degrees_of_edges_uint_32, 1 );

    %         end_vertices_structure_positions_linear_indexing_1 = cellfun( @( x ) x{  1  }, sphere_structure_positions_linear_indexing_edge_cell, 'UniformOutput', false );
    %         end_vertices_structure_positions_linear_indexing_2 = cellfun( @( x ) x{ end }, sphere_structure_positions_linear_indexing_edge_cell, 'UniformOutput', false );

        space_subscripts_int64 = int64( vertex_space_subscripts );

        position_linear_indices =   space_subscripts_int64( :, 1 )                                               ...
                                + ( space_subscripts_int64( :, 2 ) - 1 ) * size_of_image( 1 )                    ...
                                + ( space_subscripts_int64( :, 3 ) - 1 ) * size_of_image( 1 ) * size_of_image( 2 );

        vertex_position_linear_indices_cell = num2cell( position_linear_indices );

        structuring_element_linear_indexing = vertex_element_linear_indexing_templates( round( vertex_scale_subscripts ));

        vertex_structure_positions_linear_indexing = cellfun( @plus, structuring_element_linear_indexing,   ...
                                                                       vertex_position_linear_indices_cell, ...
                                                                                     'UniformOutput', false ); 

    %         vertex_structure_positions_linear_indexing = [ end_vertices_structure_positions_linear_indexing_1, ...
    %                                                        end_vertices_structure_positions_linear_indexing_2  ];

        edge_structure_positions_linear_indexing = cellfun( @cell2mat, sphere_structure_positions_linear_indexing_edge_cell, 'UniformOutput', false );

        edge_structure_positions_linear_indexing = cellfun( @( x ) unique( x, 'rows' ), edge_structure_positions_linear_indexing, 'UniformOutput', false );

        edge_structure_positions_subscript_xy = cellfun( @( x )      1 + mod( x - 1,  prod( size_of_image([ 1, 2 ]))), edge_structure_positions_linear_indexing, 'UniformOutput', false );

        edge_structure_positions_subscript_z  = cellfun( @( x ) floor( double( x ) ./ prod( size_of_image([ 1, 2 ]))), edge_structure_positions_linear_indexing, 'UniformOutput', false );

        edge_structure_positions_subscript_y  = cellfun( @( x )      1 + mod( x - 1,        size_of_image( 1 )),       edge_structure_positions_subscript_xy,    'UniformOutput', false );

        edge_structure_positions_subscript_x  = cellfun( @( x ) floor( double( x ) ./       size_of_image( 1 )),       edge_structure_positions_subscript_xy,    'UniformOutput', false );        

        edge_structure_positions_subscripts = cellfun( @( y, x, z ) [ y, x, z ], edge_structure_positions_subscript_y,                        ...
                                                                                 edge_structure_positions_subscript_x,                        ...
                                                                                 edge_structure_positions_subscript_z, 'UniformOutput', false );

        edge_centers = cellfun( @( x ) round( mean( x, 1 )), edge_structure_positions_subscripts, 'UniformOutput', false );

        edge_centers = cell2mat( edge_centers );  
                % end initialize_structuring_elements

        orig = double(h52mat(path_to_original_data));
        orig = reshape(orig,[],1);

        average_intensity_across_original =  mean(orig,'all');
         std_of_intensity_across_original = std(orig,0,'all');

        average_intensity_across_edges = cellfun(@(x) mean(orig(x)), edge_structure_positions_linear_indexing, 'UniformOutput', true);
         std_of_intensity_across_edges = cellfun(@(x)  std(orig(x)), edge_structure_positions_linear_indexing, 'UniformOutput', true);

        featurePool.normalized_mean_intensity_across_edges   = average_intensity_across_edges ./ average_intensity_across_original;
        featurePool.normalized_std_of_intensity_across_edges =  std_of_intensity_across_edges ./  std_of_intensity_across_original;

        fprintf('     '); toc
    end %edge info
   
end %function



