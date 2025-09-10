classdef wfs_renderer < base_renderer

    properties
        setup
        c
        H_pre
        amp
        delay
        dref
        WFS_filter_bank
        antialiasing_filter_bank
    end

    methods
        function obj = wfs_renderer(virtual_source,SSD,setup)
            obj = obj@base_renderer(virtual_source,SSD);
            obj.c = 343.1;
            obj.setup = setup;
            N = length(obj.virtual_source.source_signal.time_series);
            obj.H_pre = obj.get_wfs_prefilter;

            for n = 1 : length(SSD)
                obj.output_signal{n} = signal;
                obj.WFS_filter_bank{n} =  OLS_convolver(zeros(N,1),length(obj.virtual_source.source_signal.time_series));
            end
            obj.update_renderer;

        end

        function obj = update_renderer(obj, ~)
            N = length(obj.virtual_source.source_signal.time_series);
            omega = (0 : N/2)'/N*2*pi*obj.virtual_source.source_signal.fs;

            [obj.amp,obj.delay] = obj.get_delay_and_gain;
            if obj.setup.Antialiasing
                obj.get_antialiasing_filters(obj.delay);
                D_wfs = ifft(obj.antialiasing_filter_bank.* (obj.H_pre*obj.amp.').*exp(-1i*omega*obj.delay'),length(obj.virtual_source.source_signal.time_series),1,'symmetric');
            else
                D_wfs = ifft((obj.H_pre*obj.amp.').*exp(-1i*omega*obj.delay'),length(obj.virtual_source.source_signal.time_series),1,'symmetric');
            end
            for n = 1 : length(obj.WFS_filter_bank)
                obj.WFS_filter_bank{n}.update_coefficients(D_wfs(:,n));
            end
        end

        function obj = update_settings(obj,setup)
            obj.setup = setup;
            obj.update_renderer;
        end

        function render(obj)
            for n = 1 : length(obj.WFS_filter_bank)
                obj.output_signal{n}.set_signal( obj.WFS_filter_bank{n}.convolve( obj.virtual_source.source_signal.time_series ) );
            end
            obj.add_output_to_ssd_signal;
        end

        function H = get_wfs_prefilter(obj)
            % The function calculates the common prefilter present in WFS
            % driving function, given be sqrt(1i*omega)
            % Strategies:
            %   - windowing: calculates prefilters directly from evaluation
            %   the analytical formula in the frequency domain and IFFT-s
            %   into the temporal domain applying windowing
            %   - optimal: implements optimal FIR filter design as given by
            %   Schultz in Sound Field Synthesis for Line Source Array
            %   Applications in Large-Scale Sound Reinforcement PhD thesis,
            %   in Section 2.5
            Nfilt = 128;
            method = 'windowing';
            switch method
                case 'windowing'
                    N = size(obj.virtual_source.source_signal.time_series,1);
                    omega = (0 : N/2)'/N*2*pi*obj.virtual_source.source_signal.fs;
                    H = sqrt(1i*omega);
                    %                    h = fftshift(ifft(H,N,'symmetric'));
                    %                    Nfilt = N;
                    %                    h = h(round(end/2)-Nfilt/2+1:round(end/2)+Nfilt/2).*tukeywin(Nfilt,0.1);
                case 'optimal'
                    n = (-Nfilt/2:Nfilt/2-1)';
                    h = (sqrt(8*n).*cos(pi*n)+ ...
                        erfz( exp(1i*3/4*pi).*sqrt(pi*n) )-erfz( exp(1i*pi/4).*sqrt(pi*n) )).*...
                        (4*sqrt(pi)*n.^(3/2)).^(-1);
                    h(end/2+1) = sqrt(2*pi/9);
                    h = kaiser(Nfilt,4).*real(h)*Nfilt;
                case 'nara'
                    [h2, ~] = wfs_prefilter_fir('2.5D', Nfilt, 20, 20e3, obj.fs, 343);
                    h = h2(1:end-1).'*Nfilt;
            end


        end


        function [amp,delay] = get_delay_and_gain(obj)
            center = [0, 0];
            x0 = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
            n0 = cell2mat(cellfun( @(x) x.orientation, obj.secondary_source_distribution, 'UniformOutput', false)');
            dl = min( sqrt(sum((x0-circshift(x0,-1)).^2,2)), sqrt(sum((x0-circshift(x0,1)).^2,2)) );
            rhoG = sqrt(sum((x0-center).^2,2));

            switch obj.virtual_source.source_type.Shape
                case 'plane_wave'
                    A0 = .1;
                    kP = bsxfun( @minus, [0,0], obj.virtual_source.position );
                    kh = kP/norm(kP);
                    khn = sum(kh.*n0,2);
                    win = double(khn>=0);
                    dref = sqrt(rhoG);
                    delay =  (x0*kh') / obj.c;
                    amp = A0*sqrt(8*pi/obj.c).*win.*dref.*khn.*dl;
                case 'point_source'
                    kP = bsxfun( @minus, x0, obj.virtual_source.position );
                    rhoP = sqrt(sum(kP.^2,2));
                    kh = kP./sqrt(sum(kP.^2,2));
                    khn = sum(kh.*n0,2);

                    F = -1+2*all(khn<0); % Focused flag
                    if F == -1
                        win = double(khn>=0);
                    elseif F == 1
                        ks = (obj.virtual_source.position -center)/norm(obj.virtual_source.position-center);
                        win = (kh*ks').^0.*((kh*ks')>0);
                    end
                    dref = sqrt(rhoG.*rhoP./(rhoG-F*rhoP));
                    A0 = 0.1*4*pi;
                    amp = A0*sqrt(1/(obj.c*2*pi)).*win.*dref.*khn.*dl./rhoP;
                    delay = -F*rhoP / obj.c;
            end
            delay = delay - min(delay);
            [~,stat_ix] = max(abs(khn));
            win_taper = abs(circshift( win>0, -stat_ix));
            L1 = numel( find( win_taper(1:round( end/2 ))>0) );
            L2 = numel( find( win_taper(round( end/2 )+1:end)>0) );
            win1 = tukeywin(2*L1,obj.setup.Tapering);
            win2 = tukeywin(2*L2,obj.setup.Tapering);
            win_taper(1:L1) = win1(end/2+1:end);
            win_taper(end-L2+1:end) = win2(1:end/2);
            win_taper = circshift(win_taper, stat_ix-1);
            amp = amp.*win_taper;
        end

        function get_antialiasing_filters(obj,tau)

            x0 = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
            dx1 = (circshift(x0,-1)-x0);
            x1 = x0+dx1/2;
            dx2 = (circshift(x0,1)-x0);
            x2 = x0+dx2/2;
            kh1 = (x1-obj.virtual_source.position)./sqrt(sum( (x1-obj.virtual_source.position).^2,2));
            kh2 = (x2-obj.virtual_source.position)./sqrt(sum( (x2-obj.virtual_source.position).^2,2));
            dx_dir = max( [sum(kh1.*dx1,2),sum(kh2.*dx2,2)] ,[] ,2);
            R0 = mean(sqrt(sum(x0.^2,2)));
            [~,ix] = sort(x0*obj.virtual_source.position','descend');
            Nbut = 7;
            ixs_to_mod = ix(1:2);
            G = (x0(ixs_to_mod,:)')\obj.virtual_source.position'/norm(obj.virtual_source.position)*R0;
            ix0 = find(abs(G-1)<1e-9);
            if ~isempty(ix0)
                ixs_to_mod = ixs_to_mod(ix0);
                G = G(ix0);
            end

            wc = pi.*343./abs(dx_dir);

            N = length(obj.virtual_source.source_signal.time_series);
            w = (0 : N/2)'/N*2*pi*obj.virtual_source.source_signal.fs;
            [Wc,W] = meshgrid(wc,w);

            obj.antialiasing_filter_bank = 1./sqrt( 1 + ( W./Wc ).^(2*Nbut));
            H0 =1./sqrt(1+(1i*W(:,ixs_to_mod)./Wc(:,ixs_to_mod) ));

            H1 = 1./sqrt(1+(W(:,ixs_to_mod)./Wc(:,ixs_to_mod) ).^(4*Nbut));
            G1 = H1 + (1-H1).*sqrt(G').*exp(1i*W(:,ixs_to_mod).*(tau(ixs_to_mod)'-mean(tau(ixs_to_mod))))./sqrt(1-abs(dx_dir(ixs_to_mod)).^2)';

            obj.antialiasing_filter_bank(:,ixs_to_mod) = H0.*G1;
        end

        function xref = get_reference_curve(obj)

        end

    end
end

%nti xl3