function [ amp, delay, focused, AAfilt ] = get_wfs_driving_function( xs, model, x0, n0, fs, c, antialiasing )
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
            delay = delay - min(delay);
        end

    case 'plane_wave'
        k = bsxfun( @minus, 0, xs );
        k = k/norm(k);

        rho_G = sqrt(sum(x0.^2,2));
        k_n = sum(k.*n0,2);
        dl = sqrt(sum((x0-circshift(x0,1)).^2,2));

        amp = 0.1*((k_n).*double(k_n>=0)).*sqrt(rho_G).*dl*sqrt(8*pi);
        delay = round( (x0*k') / c *fs );
        delay = delay - min(delay);
        focused = -1;
end

switch antialiasing
    case 'on'
        Nbut = 4;
        Naa = 128;
        if focused == -1

            kt = sqrt(1-k_n.^2);
            [km,ixm] = max(k_n);
            kt(ixm) = 0;
            wc = pi./dl.*343./abs(kt);
            w = [(0 : Naa/2 - 1)';(-Naa/2:-1)' ]/Naa*2*pi*fs;
            [Wc,W] = meshgrid(wc,w);
            transfer = 1./sqrt(1+(W./Wc).^(2*Nbut));
            AAfilt = fftshift(ifft(transfer,[],1),1).*kaiser(Naa,4);
        else
            kt = sqrt(1-k_n.^2);
            [km,ixm] = min(kt);
            kt(ixm) = 0;

            wc = pi./dl.*343./abs(kt);
            w = [(0 : Naa/2 - 1)';(-Naa/2:-1)' ]/Naa*2*pi*fs;
            [Wc,W] = meshgrid(wc,w);
            transfer = 1./sqrt(1+(W./Wc).^(2*Nbut));
            AAfilt = fftshift(ifft(transfer,[],1),1);
        end
    case 'off'
        AAfilt = [];
end


end