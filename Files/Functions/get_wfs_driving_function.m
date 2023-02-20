function [ amp, delay, focused, AAfilt ] = get_wfs_driving_function( xs, model, x0, n0, fs, c, antialiasing )
            Nbut = 64*2;
switch model
    case 'point_source'
        k = bsxfun( @minus, x0, xs );
        rho_P = sqrt(sum(k.^2,2));
        rho_G = sqrt(sum(x0.^2,2));

    
        k_P = bsxfun(@times, k,1./rho_P);
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
            win = k_P*k_P0'.*((k_P*k_P0')>0);
            k_s = (xs-center)/norm(xs-center);
 %           win1 = k_P*k_s'.*((k_P*k_s')>0);
            win = (k_P*k_s'+1)/2;

            dref = rho_P.*rho_G./(rho_G - rho_P);
            dref = rho_P.*rho_G./(rho_G + rho_P);
            amp = -win.*k_n/sqrt(2*pi).*sqrt(dref)./rho_P.*dl;
            delay = -round( rho_P / c *fs );
            delay = delay - min(delay) + Nbut/2;
        end

    case 'plane_wave'
        k = bsxfun( @minus, 0, xs );
        k = k/norm(k);
        kx = k(1);
        ky = k(2);
        rho_G = sqrt(sum(x0.^2,2));
        %  [fi,r] = cart2pol(x0(:,1),x0(:,2));
        %  [fi0,r0] = cart2pol(k(1),k(2));
        %  d_ref = -r.*( cos(fi-fi0) - 1i*sin(fi-fi0));
        k_n = sum(k.*n0,2);
        dl = sqrt(sum((x0-circshift(x0,1)).^2,2));

        amp = ((k_n).*double(k_n>=0)).*sqrt(rho_G).*dl*sqrt(8*pi);
        delay = round( (x0*k') / c *fs );
        delay = delay - min(delay);
        focused = -1;
end

switch antialiasing
    case 'on'
        if focused == -1
            wc = pi./dl.*343./abs(sqrt(1-k_n.^2));
            N = 16;
            w = [(0 : N/2 - 1)';(-N/2:-1)' ]/N*2*pi*fs;
            [Wc,W] = meshgrid(wc,w);
            transfer = 1./sqrt(1+(W./Wc).^(2*Nbut));
            AAfilt = fftshift(ifft(transfer,[],1),1);
        else

            wc = pi./dl.*343./abs(sqrt(1-k_n.^2));
            N = 512;
            w = [(0 : N/2 - 1)';(-N/2:-1)' ]/N*2*pi*fs;
            [Wc,W] = meshgrid(wc,w);
            transfer = 1./sqrt(1+(W./Wc).^(2*Nbut));
            AAfilt = fftshift(ifft(transfer,[],1),1);
            delay = delay - Nbut*2;
        end
    case 'off'
        AAfilt = [];
end


end