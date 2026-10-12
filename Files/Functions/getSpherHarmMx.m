function [ Y_mx] = getSpherHarmMx( inclination, azimuth, N, type)
%GETSPHERHARMMX  Orthonormal spherical harmonics up to order N at the given
%directions: Y_mx [K x (N+1)^2], column n^2 + n + m + 1.
%   type 'complex': Y_n^m = sqrt((2n+1)/4pi (n-|m|)!/(n+|m|)!) P_n^|m|(cos th) e^{j m phi}
%   (no Condon-Shortley phase, Y_n^-m = conj(Y_n^m)); 'real': the cos / sin
%   (m < 0) pairs, sqrt(2) for m ~= 0.
%   The normalised Legendre functions (legendre(..., 'norm')) keep high
%   orders accurate (the factorials of the plain form overflow above
%   n ~ 85).

K = numel(inclination);
Y_mx = zeros(K, (N+1)^2);
x = cos(inclination(:));
azimuth = azimuth(:);
for n = 0:N
    m = (-n:n);
    nlin = n^2 + n + m + 1;
    % legendre 'norm' = sqrt((n+1/2)(n-m)!/(n+m)!) P_n^m, without the
    % Condon-Shortley phase -> /sqrt(2 pi) = sqrt(F1 F2) P
    Pn = legendre(n, x', 'norm')/sqrt(2*pi);          % [n+1 x K], m = 0..n
    Pnm = [flipud(Pn(2:end,:)); Pn].';                % [K x 2n+1], |m|

    if ( strcmp( type, 'complex' ) )
        Y_mx(:,nlin) = Pnm.*exp(1i*azimuth*m);
    elseif( strcmp( type, 'real' ) )
        s = ones(size(m)); s(m ~= 0) = sqrt(2);
        nix = m < 0; pix = ~nix;
        if any(nix)
            Y_mx(:,nlin(nix)) = s(nix).*Pnm(:,nix).*sin(azimuth*m(nix));
        end
        Y_mx(:,nlin(pix)) = s(pix).*Pnm(:,pix).*cos(azimuth*m(pix));
    else
        error( 'Unknown type.' );
    end
end
end
