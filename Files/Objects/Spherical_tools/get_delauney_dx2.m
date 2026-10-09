        function [dx] = get_delauney_dx2(obj,kin)
            dx = zeros(size(kin,1),1);
            for n = 1 : size(kin,1)
                [rows,cols] = find(obj.delauney_trias'==n);
                p0 = obj.x0(n,:)';
                nout = p0/norm(p0);
                trias = obj.delauney_trias(:,rows);
                w = [];
                for ti = 1 : size(trias,2)
                    ixs = setdiff((1:3),cols(ti));
                    p1 = obj.x0( trias(ixs(1),ti), :)';
                    p2 = obj.x0( trias(ixs(2),ti), :)';
                    v1 = (p1-p0) - nout*dot(nout,p1-p0);
                    v2 = (p2-p0) - nout*dot(nout,p2-p0);
                    w(:,ti) = [v1,v2]\kin(n,:)';
                end
                ix_sol = find(w(1,:)>=-1e-10&w(2,:)>=-1e-10);
                ix_sol = ix_sol(1);

                ixs = setdiff((1:3),cols(ix_sol));
                v1 = obj.x0( trias(ixs(1),ix_sol), :)'-p0;
                v2 = obj.x0( trias(ixs(2),ix_sol), :)'-p0;
                dx(n) = norm([v1,v2]*w(:,ix_sol)/norm(w(:,ix_sol),1));
                if isnan(dx(n))
                    dx(n) = 0;
                end

            end
        end
