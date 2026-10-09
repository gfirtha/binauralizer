function Anm = xyplane_sht_coeff(max_degree, max_order)
% Square of associated Legendre polynomial evaluated in the xy-plane.
Anm = cell(1, max_degree);
Anm{1} = 1;
Nprev = 1;
for n = 1:(max_degree + 1)
    mmax = min([n, max_order]);
    N = 2 * mmax + 1;
    Am = zeros(1, N);
    for m = (-mmax + mod(n + mmax, 2)):2:(mmax)-2
        Am(mod(m,N)+1) = Anm{n}(mod(m+1,Nprev)+1) * (n-m-1) / (n-m);
    end
    Am(mmax+1) = Anm{n}(mmax) * (n+mmax-1) / (n+mmax);
    Anm{n+1} = Am;
    Nprev = N;
end
end
