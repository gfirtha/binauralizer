function [ Y_mx] = getSpherHarmMx( theta, phi, N, type)

Y_mx = zeros(length(theta),N^2+1);

wb = waitbar(0,'Assembling Spherical Harmonic Matrix');
for n = 0:N
    for m = -n:n   
        waitbar((n^2+n+m+1)/((N+1)^2) ,wb);
        Y_mx(:,n^2+n+m+1) = getSpherHarm( theta, phi, n, m, type );
    end
end
close(wb);

end

