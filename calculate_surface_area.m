function calculate_surface_area

[x,y,z,v] = flow;
p = patch(isosurface(x,y,z,v,-3));
isonormals(x,y,z,v,p)
set(p,'FaceColor','red','EdgeColor','none');
daspect([1 1 1])
view(3); 
camlight 
lighting gouraud
verts = get(p, 'Vertices');
faces = get(p, 'Faces');
a = verts(faces(:, 2), :) - verts(faces(:, 1), :);
b = verts(faces(:, 3), :) - verts(faces(:, 1), :);
c = cross(a, b, 2);
area = 1/2 * sum(sqrt(sum(c.^2, 2)));
fprintf('\nThe surface area is %f\n\n', area)