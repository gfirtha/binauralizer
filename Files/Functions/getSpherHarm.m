function [ Ynm ] = getSpherHarm( inclination, azimuth, n, m, type )

Pnm = getLegendre( n, abs(m), cos( inclination(:)  ) );

F1 = (2*n+1)/(4*pi);
F2 = factorial(n-abs(m))/factorial(n+abs(m));

if ( strcmp( type, 'complex' ) )
    Ynm = (-1).^m *sqrt( F1 .*F2 ) * reshape(Pnm,size(inclination,1),size(inclination,2)).*exp(1i*m*azimuth);
elseif( strcmp( type, 'real' ) )
    
    if ( m ~= 0 )
       F1 = 2 * F1;
    end
    
    if ( m < 0 )
        Ynm = (-1).^m .* sqrt( F1.*F2)*reshape(Pnm,size(inclination,1),size(inclination,2)).*sin( m*azimuth );
    else
        Ynm = (-1).^m .* sqrt( F1.*F2)*reshape(Pnm,size(inclination,1),size(inclination,2)).*cos( m*azimuth);
    end
else
    error( 'Unknown type.' );
end
end

