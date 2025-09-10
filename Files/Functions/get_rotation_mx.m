function [R] = get_rotation_mx(u,v)
u = u/norm(u);
v = v/norm(v);
a = cross(u,v)';
sinfi = 1;
cosfi = dot(u,v);
R = [cosfi+ a(1)^2*(1-cosfi), a(1)*a(2)*(1-cosfi)-a(3)*sinfi, a(1)*a(3)*(1-cosfi)+a(2)*sinfi;
    a(1)*a(2)*(1-cosfi)+a(3)*sinfi, cosfi+a(2)^2*(1-cosfi), a(2)*a(3)*(1-cosfi)-a(1)*sinfi;
    a(1)*a(3)*(1-cosfi)-a(2)*sinfi, a(2)*a(3)*(1-cosfi)+ a(1)*sinfi, cosfi+a(3)^3*(1-cosfi)];

end

