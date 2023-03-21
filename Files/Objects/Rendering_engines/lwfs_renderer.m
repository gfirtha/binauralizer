classdef lwfs_renderer < base_renderer

    properties
        fs
        c
        receiver
        focused_flag
        H_pre
        WFS_filter_bank
        omega
        antialiasing_filter_bank
    end

    methods
        function obj = lwfs_renderer(virtual_source,SSD, fs)
            obj = obj@base_renderer(virtual_source,SSD);
            obj.fs = fs;
            obj.c = 343.1;
            obj.receiver = [0, 0];
            N = length(obj.virtual_source.source_signal.time_series);
            obj.omega = 2*pi*(0:N-1)'/N*obj.fs;
            obj.H_pre = fft(obj.get_wfs_prefilter,N);

            for n = 1 : length(SSD)
                obj.output_signal{n} = signal;
                obj.WFS_filter_bank{n} =  OLS_convolver(zeros(N,1),length(obj.virtual_source.source_signal.time_series));
            end
            obj.update_renderer;
        end

        function obj = update_renderer(obj, ~)
            [amp,delay] = obj.get_delay_and_gain;
            obj.get_antialiasing_filters(length(obj.omega));
            D_wfs = ifft(obj.antialiasing_filter_bank.*(obj.H_pre*amp.').*exp(-1i*obj.omega*delay'),length(obj.omega),1,'symmetric');
            for n = 1 : length(obj.WFS_filter_bank)
                obj.WFS_filter_bank{n}.update_coefficients(D_wfs(:,n));
            end
        end

        function render(obj)
            for n = 1 : length(obj.WFS_filter_bank)
                obj.output_signal{n}.set_signal( obj.WFS_filter_bank{n}.convolve( obj.virtual_source.source_signal.time_series ) );
            end
            obj.add_output_to_ssd_signal;
        end

        function h = get_wfs_prefilter(obj)
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
                    om = (0 : N/2 - 1)'/N*2*pi*obj.fs;
                    H = sqrt(1i*om);
                    h = fftshift(ifft(H,N,'symmetric'));
                    h = h(round(end/2)-Nfilt/2+1:round(end/2)+Nfilt/2).*hann(Nfilt);
                case 'optimal'
                    n = (-Nfilt/2:Nfilt/2-1)';
                    h = (sqrt(8*n).*cos(pi*n)+ ...
                        erfz( exp(1i*3/4*pi).*sqrt(pi*n) )-erfz( exp(1i*pi/4).*sqrt(pi*n) )).*...
                        (4*sqrt(pi)*n.^(3/2)).^(-1);
                    h(end/2+1) = sqrt(2*pi/9);
                    h = kaiser(Nfilt,4).*real(h);
            end
        end


        function [amp,delay] = get_delay_and_gain(obj)
            x0 = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
            n0 = cell2mat(cellfun( @(x) x.orientation, obj.secondary_source_distribution, 'UniformOutput', false)');
            dl = sqrt(sum((x0-circshift(x0,1)).^2,2));
            rho_G = sqrt(sum(x0.^2,2));

            switch obj.virtual_source.source_type.Shape
                case 'plane_wave'
                    A0 = 0.1;
                    kP = bsxfun( @minus, obj.receiver, obj.virtual_source.position );
                    kh = kP/norm(kP);
                    khn = sum(kh.*n0,2);
                    win = double(khn>=0);
                    dref = sqrt(rho_G);
                    delay =  (x0*kh') / obj.c;
                    amp = A0*sqrt(8*pi/obj.c).*win.*dref.*khn.*dl;
                case 'point_source'
                    kP = bsxfun( @minus, x0, obj.virtual_source.position );
                    rho_P = sqrt(sum(kP.^2,2));
                    dref = sqrt(rho_P.*rho_G./(rho_P + rho_G));
                    kh = kP./sqrt(sum(kP.^2,2));
                    khn = sum(kh.*n0,2);

                    F = -1+2*all(khn<0); % Focused flag
                    if F == -1
                        win = double(khn>=0);
                    elseif F == 1
                        center = [0, 0];
                        ks = (obj.virtual_source.position -center)/norm(obj.virtual_source.position-center);
                        %win = (kh*ks'+1)/2;
                        win = kh*ks'.*((kh*ks')>0);
                    end
                    A0 = 1;%*norm(bsxfun( @minus, obj.receiver, obj.virtual_source.position ))*4*pi;
                    amp = A0*sqrt(-F/(obj.c*2*pi)).*win.*dref.*khn.*dl./rho_P;
                    delay = -F*rho_P / obj.c;
            end
            delay = delay - min(delay);
        end

        function get_antialiasing_filters(obj,N)
            Nbut = 4;
            x0 = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
            n0 = cell2mat(cellfun( @(x) x.orientation, obj.secondary_source_distribution, 'UniformOutput', false)');
            dl = sqrt(sum((x0-circshift(x0,1)).^2,2));
            switch obj.virtual_source.source_type.Shape
                case 'plane_wave'
                    kP = bsxfun( @minus, obj.receiver, obj.virtual_source.position );
                case 'point_source'
                    kP = bsxfun( @minus, x0, obj.virtual_source.position );
            end
            kh = kP./sqrt(sum(kP.^2,2));
            khn = sum(kh.*n0,2);
            kht = sqrt(1-khn.^2);
            [~,ixm] = max(khn);
            kht(ixm) = 0;
            wc = pi./dl.*343.1./abs(kht);
            w = [(0 : N/2 - 1)';(-N/2:-1)' ]/N*2*pi*obj.fs;
            [Wc,W] = meshgrid(wc,w);
            obj.antialiasing_filter_bank = 1./sqrt(1+(W./Wc).^(2*Nbut));
        end

    end
end

