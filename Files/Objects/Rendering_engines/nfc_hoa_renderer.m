classdef nfc_hoa_renderer < base_renderer

    properties
        fs
        c
        receiver
        HOA_filter_bank
        omega
        M
    end

    methods
        function obj = nfc_hoa_renderer(virtual_source,SSD, fs, M)
            obj = obj@base_renderer(virtual_source,SSD);
            obj.fs = fs;
            obj.c = 343.1;
            obj.receiver = [0, 0];
            obj.M = M;
            N = length(obj.virtual_source.source_signal.time_series);
            obj.omega = 2*pi*(0:N-1)'/N*obj.fs;

            for n = 1 : length(SSD)
                obj.output_signal{n} = signal;
                obj.HOA_filter_bank{n} =  OLS_convolver(zeros(N,1),length(obj.virtual_source.source_signal.time_series));
            end
            obj.update_renderer;
        end

        function obj = update_renderer(obj, ~)
            x0 = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
            phi_ssd = cart2pol(x0(:,1),x0(:,2));

            R0 = mean(sqrt(sum(x0.^2,2)));
            dR = mean(diff(unwrap(phi_ssd)))*R0;
            switch obj.virtual_source.source_type.Shape
                case 'plane_wave'
                    A0 = 0.1;
                    k0 = bsxfun( @minus, obj.receiver, obj.virtual_source.position);
                    phi_pw = cart2pol(k0(1),k0(2));
                    D_hoa = zeros(length(obj.omega),length(obj.secondary_source_distribution));
                    for m = -obj.M:obj.M
                        Gm = getSphH( abs(m), 2, obj.omega/obj.c*R0 )./exp(-1i*obj.omega/obj.c*R0);
                        D_hoa = D_hoa - A0*2/R0*1i^(-abs(m)).*1./(1i*obj.omega/obj.c.*Gm)*exp(1i*m*(phi_pw-phi_ssd))'*dR;
                    end
                case 'point_source'
                    k0 = bsxfun( @minus, obj.receiver, obj.virtual_source.position);
                    [phi_s,rs] = cart2pol(k0(1),k0(2));
                    D_hoa = zeros(length(obj.omega),length(obj.secondary_source_distribution));
                    for m = -obj.M:obj.M
                        Gm = getSphH( abs(m), 2, obj.omega/obj.c*rs )./getSphH( abs(m), 2, obj.omega/obj.c*R0 );
                        D_hoa = D_hoa + 1/(2*pi*R0)*Gm*exp(1i*m*((phi_s-pi)-phi_ssd))'*dR;
                    end
            end
            D_hoa(isnan(D_hoa)) = 0;
            d_hoa = ifft(D_hoa,length(obj.omega),1,'symmetric');
            for n = 1 : length(obj.HOA_filter_bank)
                obj.HOA_filter_bank{n}.update_coefficients(d_hoa(:,n));
            end
        end

        function render(obj)
            for n = 1 : length(obj.HOA_filter_bank)
                obj.output_signal{n}.set_signal( obj.HOA_filter_bank{n}.convolve( obj.virtual_source.source_signal.time_series ) );
            end
            obj.add_output_to_ssd_signal;
        end

    end
end

