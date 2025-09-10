function [ out ] = SHT_new( in, sphere, N, mode, type )
% A numerikus SHT megvalósítása

switch mode
    case 'trans'
        Y = getSpherHarmMx( zenith, azim, N, type );
        out = zeros(size(Y,2),size(in,2),size(in,3));
        for m = 1 : size(in,3)
            out(:,:,m) = (4*pi/length(zenith)*Y')*squeeze(in(:,:,m));
        end
    case 'pinv'
        Y = getSpherHarmMx( zenith, azim, N, type );
        out = zeros(size(Y,2),size(in,2),size(in,3));
        Yi = (Y'*Y)\Y';
        for m = 1 : size(in,3)
            out(:,:,m) = Yi*squeeze(in(:,:,m));
        end
    case 'num_int'
        in;

end

end

