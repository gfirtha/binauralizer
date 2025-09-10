function [ dj ] = getDifSphJ( n, x )
    
    dj = -getSphJ(n+1,x) + (n./x).*getSphJ(n,x);

end

