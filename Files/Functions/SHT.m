function [ out ] = SHT( in, theta, phi, N, mode, type )
% A numerikus SHT megvalósítása 
    Y = getSpherHarmMx( theta, phi, N, type );
    if strcmp(mode,'trans')
        out = (4*pi/length(theta)*Y')*in;
    elseif strcmp(mode,'pinv')
        out = ((Y'*Y)\Y')*in;
    end

end

