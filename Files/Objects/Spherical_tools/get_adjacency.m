
        function adjacency = get_adjacency(obj,N,trias)
                    for n = 1 : N
                        ixs = (trias(1,:)==n)|(trias(2,:)==n)|(trias(3,:)==n);
                        adjacency{n} = setdiff(unique(trias(:,ixs)),n);
                    end
        end