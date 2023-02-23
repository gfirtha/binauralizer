function [ Lnm ] = getLegendre( n, m, arg )

Lnm = legendre( n, arg );

if ( n ~= 0 )
    Lnm = squeeze( Lnm( abs( m ) + 1, :, : ) );  
end

if ( m < 0 )
    Lnm = ( -1 ).^abs( m ) .* factorial( n - abs( m ) ) ./ factorial( n + abs( m ) ) .* Lnm;
end

Lnm = reshape( Lnm, size( arg ) );

end