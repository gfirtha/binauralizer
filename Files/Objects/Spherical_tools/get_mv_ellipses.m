function obj = get_mv_ellipses(obj)
                    K = obj.voronoi_cells.connectivity;
                    vertices = obj.voronoi_cells.vertices;
                    Tr_mx = zeros(2,3,length(K));
                    R = zeros(2,length(K));
                    for n = 1 : length(K)
                        p1 = vertices(:,K{n}(1));
                        p2 = vertices(:,K{n}(2));
                        p3 = vertices(:,K{n}(end));
                        v1 = (p2-p1)/norm(p2-p1);
                        v2 = (p3-p1)/norm(p3-p1);
                        Tmx = [v1 ,v2 ];
                        points = [obj.x0(n,:);(obj.x0(obj.voronoi_cells.adjacency{n},:))]';
                        tr_points = Tmx\(points-p1);
                        [A, c0] = mvee_fit(tr_points);
                        [~, D, V] = svd(A);
                        R(:,n) = 1./sqrt(diag(D));
                        Tr_mx(:,:,n) = pinv(Tmx*V);
                        center(:,n) = Tmx*c0 + p1;
                    end
                    obj.min_vol_ellipses = struct('Tmx',Tr_mx,'R',R,'center',center);
                end
