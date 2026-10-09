        function [dx] = get_voronoi_edge_dx(obj,kin)
            kin = kin./sqrt(sum(kin.^2,2));
            dx = zeros(size(kin,1),1);
            for ki = 1 : size(kin,1)
                k0 = kin(ki,:);
                if norm(k0)<1e-9
                    dx(ki) = 1e-20;
                else
                    % Choose furthest vertices
                    K = obj.voronoi_cells.connectivity{ki};
                    V = obj.voronoi_cells.vertices(:,K);
                    edges_ixs = [(1:size(V,2))', circshift((1:size(V,2))',-1)];

                    [~,ix] = max(sum((V(:,edges_ixs(:,1))'-V(:,edges_ixs(:,2))').^2,2));
                    p0 = mean(V,2);
                    p1 = V(:,edges_ixs(ix(1),1));
                    p2 = V(:,edges_ixs(ix(1),2));

                    v1 = (p1-p0)/norm(p1-p0);
                    v2 = (p2-p0)/norm(p2-p0);
                    Tr = [v1,v2,cross(v1,v2)/norm(cross(v1,v2))];
                    V_tr = Tr(:,1:2)\(V-p0);
                    c_tr = Tr(:,1:2)\(obj.x0(ki,:)'-p0);
                    c = Tr*[c_tr;0] + p0;
                    kr_tr = Tr(:,1:2)\k0';
                    kr_tr = kr_tr/norm(kr_tr);
                    m = 0;
                    for n = 1 : size(edges_ixs,1)
                        l0 = V_tr(:,edges_ixs(n,1));
                        l1 = V_tr(:,edges_ixs(n,2));
                        l = l1-l0;
                        d = [kr_tr, -l]\(l0-c_tr);
                        if d(1)>=0 && abs(d(1))<=max(norm(l1-c_tr),norm(l0-c_tr)) && d(2)>=0% && abs(d(2))<=norm(l)
                            m = m+1;
                            out(:,m) = d;
                        end
                    end
                    x_out_tr = c_tr + kr_tr*out(1);
                    x_out = Tr*[x_out_tr;0] + p0;
                    dx(ki) = norm(x_out-c);
                end
            end
        end
