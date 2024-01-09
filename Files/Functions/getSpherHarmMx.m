function [ Y_mx] = getSpherHarmMx( inclination, azimuth, N, type)

Y_mx = zeros(length(inclination),N^2+1);
%wb = waitbar(0,'Assembling Spherical Harmonic Matrix');
for n = 0:N
 %   waitbar(n/N ,wb);
    m = (-n:n);
    nlin = n^2 + n + m + 1;

    F1 =  (2*n+1)/(4*pi);
    F2 = factorial(n-abs(m))./factorial(n+abs(m));
    Pnm = legendre( n, cos( inclination(:)  ) );
    Pnm = [flipud(Pnm(2:end,:));Pnm].';

    if ( strcmp( type, 'complex' ) )
        Y_mx(:,nlin) = (-1).^m.*sqrt(F1 .* F2).*Pnm.*exp(1i*azimuth*m);
    elseif( strcmp( type, 'real' ) )
        F2(m ~= 0) = 2*F2(m ~= 0);
        nix = m<0;
        pix = m>=0;
        if ~isempty(find(nix, 1))
        Y_mx(:,nlin(nix)) = (-1).^m(nix) .* sqrt( F1.*F2(nix)).*Pnm(:,nix).*sin(azimuth*m(nix));
        end
        Y_mx(:,nlin(pix)) = (-1).^m(pix) .* sqrt( F1.*F2(pix)).*Pnm(:,pix).*cos(azimuth*m(pix));

    else
        error( 'Unknown type.' );
    end
end
%close(wb);
end

