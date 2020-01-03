function [ amp, delay ] = get_wfs_driving_function( xs, x0, n0, fs, c )

    v = bsxfun( @minus, x0, xs );
    rho_P = sqrt(sum(v.^2,2));
    rho_G = sqrt(sum(x0.^2,2));
    k_P = bsxfun(@times, v,1./rho_P);
    k_n = sum(k_P.*n0,2);
    dl = sqrt(sum((x0-circshift(x0,1)).^2,2));
    amp = ((k_n).*double(k_n>=0))/sqrt(2*pi).*sqrt(rho_P.*rho_G./(rho_P + rho_G))./rho_P.*dl;
    delay = round( rho_P / c *fs );

end

