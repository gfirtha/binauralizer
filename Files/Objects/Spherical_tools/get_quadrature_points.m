function [quadrature_points] = get_quadrature_points(delauney_trias, x0)
            order = 10;
            [xx, ww] = gaussquad2(order, 3);
            tria_set = ShapeSet.LinearTria;
            w = tria_set.eval(xx);

            x_quad = [];
            w_quad = [];
            w_interp = [];
            ind_interp = [];
            Ssphere = 0;
            for n = 1 : size(delauney_trias,2)
                tria = delauney_trias(:,n)';

                p0 = x0(tria(1),:);
                v1 = x0(tria(2),:) - p0;
                v2 = x0(tria(3),:) - p0;
                normal = cross(v1',v2');
                S = norm(cross(v1',v2'))/2;
                Ssphere = Ssphere + S;
                normal = normal ./ sqrt( sum( normal.^2 ) );
                W = [v1',v2',normal];
                x_quad  = [x_quad ;  (W*[xx,zeros(size(xx,1),1)]')'+p0 ];
                w_quad = [w_quad; ww*S/0.5];
                w_interp = [w_interp; w];
                ind_interp = [ind_interp; repmat(tria,[size(xx,1),1])];
            end
            quadrature_points = struct('quad_pos',x_quad,'quad_w',w_quad,'interp_w',w_interp,'interp_pos',ind_interp);
end

