        function dS = get_tria_areas(varargin)
            if isempty(varargin)
                dS = [];
                for n = 1 : length(obj.x0)

                    [row,col] = find(obj.delauney_trias == n);
                    trias = obj.delauney_trias(:,col);
                    area = [];
                    for m = 1 : length(col)
                        tria = trias(:,m);
                        vertices = obj.x0(tria,:)';
                        v1 = vertices(:,2) - vertices(:,1);
                        v2 = vertices(:,3) - vertices(:,1);
                        area(m) = norm(cross(v1, v2, 1)) / 2;
                    end
                    dS(n) = sum(area)/3;
                end
            else
                ix0 = varargin{1};
                [row,col] = find(obj.delauney_trias == ix0);
                trias = obj.delauney_trias(:,col);
                area = [];
                for m = 1 : length(col)
                    tria = trias(:,m);
                    vertices = obj.x0(tria,:)';
                    v1 = vertices(:,2) - vertices(:,1);
                    v2 = vertices(:,3) - vertices(:,1);
                    area(m) = norm(cross(v1, v2, 1)) / 2;
                end
                dS = sum(area)/3;
            end
        end
