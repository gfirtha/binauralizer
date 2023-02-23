    function [ out ] = ISHT( in, theta, phi, mode )
% This function performs inverse Speherical Harmonic Transform
    N = sqrt(size(in,1))-1;
    Y = getSpherHarmMx( theta, phi, N, mode );
    out = Y*in;
end