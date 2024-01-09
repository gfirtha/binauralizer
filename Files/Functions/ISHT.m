function [ out ] = ISHT( in, N, inclination, azimuth, mode )
% This function performs inverse Speherical Harmonic Transform
%    N = sqrt(size(in,1))-1;
Y = getSpherHarmMx( inclination, azimuth, N, mode );
out = zeros(size(inclination,1),size(in,2),size(in,3));
for m = 1 : size(in,3)
    out(:,:,m) = Y*squeeze(in(:,:,m));
end
end