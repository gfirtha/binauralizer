function dS = get_spherical_polygon_area(K,vert)
R = mean(sqrt( sum( vert.^2, 1) ));
for n = 1 : length(K)
    S = 0;
    for s = 1 : length(K{n})-2
        idx = K{n}([1 s+1 s+2]);
        v0 = vert(:,idx(1));
        v1 = vert(:,idx(2));
        v2 = vert(:,idx(3));
        a = (dot(v0,v1)/norm(v0)/norm(v1));
        b = (dot(v0,v2)/norm(v0)/norm(v2));
        c = (dot(v1,v2)/norm(v1)/norm(v2));
        if (abs(a-1)<1e-10) || (abs(b-1)<1e-10) || (abs(c-1)<1e-10)
            continue
        end
        A = acos( (a-b*c)./(sin(acos(b))*sin(acos(c))) );
        B = acos( (b-c*a)./(sin(acos(c))*sin(acos(a))) );
        C = acos( (c-a*b)./(sin(acos(a))*sin(acos(b))) );
        S = S + real((A + B + C - pi)*R^2);
    end
    dS(n) = S;
end
end
