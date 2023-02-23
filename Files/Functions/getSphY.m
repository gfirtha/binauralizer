function [ y ] = getSphY( n, x )

    y = (-1)^(n+1)*sqrt(pi./(2*x)).*besselj(-n-1/2,x);
    y(x == 0) = NaN;
end

