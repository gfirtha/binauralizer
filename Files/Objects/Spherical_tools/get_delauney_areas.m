        function [S] = get_daluney_areas(delauney_trias,x0)
            p0 = x0(delauney_trias(1,:)',:);
            v1 = x0(delauney_trias(2,:)',:)-p0;
            v2 = x0(delauney_trias(3,:)',:)-p0;
            S = sqrt(sum((cross(v1,v2,2)).^2,2))/2;
        end
