function draw_sphere(K,vertices, normals, x_in, y_in, z_in)

figure
for n = 1 : size(K,1)
    patch(vertices(1,K{n}),vertices(2,K{n}),vertices(3,K{n}),'red','FaceAlpha',0.5)
    hold on
end
plot3(x_in,y_in,z_in,'k.');
hold on

quiver3(x_in,y_in,z_in,normals(1,:)',normals(2,:)',normals(3,:)','k') 
view(3)
axis equal tight

end

