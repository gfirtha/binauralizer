function [N] = get_maxN(sphere)


N0 = ceil(sqrt(size(sphere.Phi0,1))-1);
w = sphere.voronoi_cells.dS;
k(1) = 1;
n = 1;
while k(n)<10
n = n+1;
Y = getSpherHarmMx( pi/2-sphere.Phi0(:,2), sphere.Phi0(:,1), n, 'real' );
Yt = pinv(Y'*diag(w)*Y)*Y'*diag(w);
k(n) = cond(Yt*Y);
end
N = n-1;

end

