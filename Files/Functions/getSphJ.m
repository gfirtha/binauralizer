function [ j ] = getSphJ( n, x )

    j = sqrt(pi./(2*x)).*besselj(n+1/2,x);
    j(n==0 & x == 0) = 1;    
    j(n~=0 & x == 0) = 0;

end

