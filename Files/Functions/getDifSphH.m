function [ dh ] = getDifSphH( n,o,x )
    
    dh = -getSphH(n+1,o,x) + (n./x).*getSphH(n,o,x);

end

