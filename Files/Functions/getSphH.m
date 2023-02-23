function [ h ] = getSphH( n, o, x )

h = sqrt(pi./(2*x)).*besselh(n+1/2,o,x);

end

