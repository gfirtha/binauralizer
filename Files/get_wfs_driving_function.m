function [ amp, delay, focused, AAfilt ] = get_wfs_driving_function( xs, x0, n0, fs, c, antialiasing )

v = bsxfun( @minus, x0, xs );
rho_P = sqrt(sum(v.^2,2));
rho_G = sqrt(sum(x0.^2,2));
k_P = bsxfun(@times, v,1./rho_P);
k_n = sum(k_P.*n0,2);
dl = sqrt(sum((x0-circshift(x0,1)).^2,2));

if ~all(k_n<0)
    % Non-focused case
    focused = -1;
    amp = ((k_n).*double(k_n>=0))/sqrt(2*pi).*sqrt(rho_P.*rho_G./(rho_P + rho_G))./rho_P.*dl;
    delay = round( rho_P / c *fs );
else
    focused = 1;
    center = [0, 0];
    k_P0 = (xs-center)/norm(xs-center);
    amp = ((k_P*k_P0')>0).*k_n/sqrt(2*pi).*sqrt(rho_P.*rho_G./(rho_P + rho_G))./rho_P.*dl;
    delay = round( rho_P / c *fs );
    delay = max(delay) - delay;
end

switch antialiasing
    case 'on'
        wc = pi./dl.*343./abs(1-k_n);
        N = 512;
        w = [(0 : N/2 - 1)';(-N/2:-1)' ]/N*2*pi*fs;
        [Wc,W] = meshgrid(wc,w);
        AAfilt = ifft(1./sqrt(1+(W./Wc).^10),[],1);
    case 'off'
        AAfilt = [];
end

end