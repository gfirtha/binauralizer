function [quadrature_points] = get_adaptive_quadrature_points(delauney_trias, x0)
            edges = unique(sort([delauney_trias([1,2],:)';delauney_trias([1,3],:)';delauney_trias([2,3],:)'],2),'rows');
            Le = sqrt(sum((x0(edges(:,1),:)-x0(edges(:,2),:)).^2,2));
            Lmin = min(Le);
            Lmax = max(Le);
            min_order = 6;
            max_order = 20;

            x_quad = [];
            w_quad = [];
            w_interp = [];
            ind_interp = [];
            for n = 1 : size(delauney_trias,2)
                tria = delauney_trias(:,n)';
                p0 = x0(tria(1),:);
                v1 = x0(tria(2),:) - p0;
                v2 = x0(tria(3),:) - p0;
                normal = cross(v1',v2');
                Stria = norm(cross(v1',v2'))/2;

                L = max([norm(v1),norm(v2),norm(v1-v2)]);
                s0 = (L-Lmin)/(Lmax-Lmin);
                order = round(((max_order-min_order)*s0+min_order));
                %                [xx, ww] = gaussquad2(order, 3);
                Q = quadtriangle(order,'Domain',[0,0;1 0;0 1]);
                xx = Q.Points;
                ww = Q.Weights;
                tria_set = ShapeSet.LinearTria;
                w = tria_set.eval(xx);

                normal = normal ./ sqrt( sum( normal.^2 ) );
                W = [v1',v2',normal];
                x_quad  = [x_quad ;  (W*[xx,zeros(size(xx,1),1)]')'+p0 ];
                w_quad = [w_quad; ww*Stria/0.5];
                w_interp = [w_interp; w];
                ind_interp = [ind_interp; repmat(tria,[size(xx,1),1])];
            end
            quadrature_points = struct('quad_pos',x_quad,'quad_w',w_quad,'interp_w',w_interp,'interp_pos',ind_interp);
end
