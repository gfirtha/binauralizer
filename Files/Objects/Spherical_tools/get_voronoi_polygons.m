function [v, K, n_in, dS] = get_voronoi_polygons( delauney_trias,x0 )

xin = x0';

%% Get voronoi vertices
v = zeros(3,size(delauney_trias,2));
for n = 1 : length(delauney_trias)
%         i1 = delauney_trias(1,n);
%         i2 = delauney_trias(2,n);
%         i3 = delauney_trias(3,n);
%         v(1,n) = ( xin(2,i2) - xin(2,i1) ) * ( xin(3,i3) - xin(3,i1) ) ...
%             - ( xin(3,i2) - xin(3,i1) ) * ( xin(2,i3) - xin(2,i1) );
%         v(2,n) = ( xin(3,i2) - xin(3,i1) ) * ( xin(1,i3) - xin(1,i1) ) ...
%             - ( xin(1,i2) - xin(1,i1) ) * ( xin(3,i3) - xin(3,i1) );
%         v(3,n) = ( xin(1,i2) - xin(1,i1) ) * ( xin(2,i3) - xin(2,i1) ) ...
%             - ( xin(2,i2) - xin(2,i1) ) * ( xin(1,i3) - xin(1,i1) );
%         v(1:3,n) = v(1:3,n) ./ sqrt ( sum ( v(1:3,n).^2 ) ).*mean(sqrt(sum(x0.^2,2)));
      [~, I, ~] = triangle_circumcircle(x0(delauney_trias(1,n),:)',x0(delauney_trias(2,n),:)', x0(delauney_trias(3,n),:)')  ;
      v(:,n) = I/norm(I)*mean(sqrt(sum(x0.^2,2)));
end
%% Get connectivity between points
for n = 1 : size(x0,1)
    % For each x0: center of Voronoi polygon find the delauney trias,
    % containing x0. Delauney indices = Voronoi vertex indices
    [rows,cols] = find(delauney_trias' == n);
    % Voronoi vertices : v(:,rows)
    [order, nin] = sort_vertices(v(:,rows'));
 %   n_in(:,n) = x0(n,:)./norm(x0(n,:));
    n_in(:,n) = nin;
    K{n} = rows(order);
end
dS = get_spherical_polygon_area(K,v);

end

