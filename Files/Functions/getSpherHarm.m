function [ Ynm ] = getSpherHarm( Theta, Phi, n, m, type )

Pnm = getLegendre( n, abs(m), cos( Theta(:) ) );
F1 = (2*n+1)/(4*pi);
F2 = factorial(n-abs(m))/factorial(n+abs(m));

if ( strcmp( type, 'complex' ) )
    Ynm = (-1).^m *sqrt( F1 .*F2 ) * reshape(Pnm,size(Theta,1),size(Theta,2)).*exp(1i*m*Phi);
elseif( strcmp( type, 'real' ) )
    
    if ( m ~= 0 )
       F1 = 2 * F1;
    end
    
    if ( m < 0 )
        Ynm = (-1).^m .* sqrt( F1.*F2)*reshape(Pnm,size(Theta,1),size(Theta,2)).*sin( m*Phi );
    else
        Ynm = (-1).^m .* sqrt( F1.*F2)*reshape(Pnm,size(Theta,1),size(Theta,2)).*cos( m*Phi);
    end
else
    error( 'Unknown type.' );
end
end

