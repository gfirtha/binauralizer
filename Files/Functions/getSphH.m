function [ h ] = getSphH( n, o, x )
h = zeros(length(n),length(x));
for ni = 1 : length(n)
    h(ni,:) = sqrt(pi./(2*x)).*besselh(n(ni)+1/2,o,x);
end
end

