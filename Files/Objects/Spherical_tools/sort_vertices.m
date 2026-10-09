function [order, n_in] = sort_vertices(V)

    % Translate to origin and normalize
    cent = mean(V,2);
    V = V - cent;
    V = V./sqrt(sum(V.^2,1));
    % Rotate it to xy plane by choosing most adjacent vector pairs and
    % rotate their normal to [0,0,-1]
    [r,c] = find(abs(V'*V) == min(min(abs(V'*V))));
    n_in = cross(V(:,r(1)),V(:,c(1)))/norm(cross(V(:,r(1)),V(:,c(1))));
    V = [V(:,r(1)), cross(V(:,r(1)),n_in), n_in]\V;
    if acos(dot(n_in,cent/norm(cent)))<pi/2
        n_in = -n_in;
    end
    % Sort in terms of polar angle
    [~,order] = sort(atan2(V(2,:),V(1,:)));
end

