classdef nfc_hoa_renderer < base_renderer

    properties
        c
        receiver
        HOA_filter_bank
        setup
    end

    methods
        function obj = nfc_hoa_renderer(virtual_source,SSD, Norder)
            obj = obj@base_renderer(virtual_source,SSD);
            obj.c = 343.1;
            obj.receiver = [0, 0];
            obj.setup.Norder = Norder;
            N = length(obj.virtual_source.source_signal.time_series);

            for n = 1 : length(SSD.loudspeakers)
                obj.output_signal{n} = signal;
                obj.HOA_filter_bank{n} =  OLS_convolver(zeros(N,1),length(obj.virtual_source.source_signal.time_series));
            end
            obj.update_renderer;
        end
        function info = get_renderer_info(obj)
            info = "The renderer implements Near-field compensated Higher Order Ambisonics rendering.";
        end

        function obj = update_renderer(obj, ~)
            phi_ssd = cart2pol(obj.ssd.x0(:,1),obj.ssd.x0(:,2));

            N = length(obj.virtual_source.source_signal.time_series);
            fs = obj.virtual_source.source_signal.fs;

            omega = 2*pi*(0:N-1)'/N*fs;

            R0 = mean(sqrt(sum(obj.ssd.x0.^2,2)));
            dR = mean(diff(unwrap(phi_ssd)))*R0;
            switch obj.virtual_source.source_type.Shape
                case 'plane_wave'
                    A0 = 0.1;
                    k0 = bsxfun( @minus, obj.receiver, obj.virtual_source.position);
                    phi_pw = cart2pol(k0(1),k0(2));
                    D_hoa = zeros(length(omega),length(obj.ssd.loudspeakers));
                    for m = -obj.setup.Norder:obj.setup.Norder
                        Gm = getSphH( abs(m), 2, omega/obj.c*R0 ).'./exp(-1i*omega/obj.c*R0);
                        D_hoa = D_hoa - A0*2/R0*1i^(-abs(m)).*1./(1i*omega/obj.c.*Gm)*exp(1i*m*(phi_pw-phi_ssd))'*dR;
                    end
                case 'point_source'
                    k0 = bsxfun( @minus, obj.receiver, obj.virtual_source.position);
                    [phi_s,rs] = cart2pol(k0(1),k0(2));
                    D_hoa = zeros(length(omega),length(obj.ssd.loudspeakers));
                    for m = -obj.setup.Norder : obj.setup.Norder
                        Gm = getSphH( abs(m), 2, omega/obj.c*rs ).'./getSphH( abs(m), 2, omega/obj.c*R0 ).';
                        D_hoa = D_hoa + 1/(2*pi*R0)*Gm*exp(1i*m*((phi_s-pi)-phi_ssd))'*dR;
                    end
            end
            D_hoa(isnan(D_hoa)) = 0;
            d_hoa = fftshift( ifft(D_hoa,length(omega),1,'symmetric') , 1);
            for n = 1 : length(obj.HOA_filter_bank)
                obj.HOA_filter_bank{n}.update_coefficients(d_hoa(:,n));
            end
        end

        function obj = update_settings(obj,setup)
            obj.setup.Norder = setup.HOAorder;
            obj.update_renderer;
        end

        function render(obj)
            for n = 1 : length(obj.HOA_filter_bank)
                obj.output_signal{n}.set_signal( obj.HOA_filter_bank{n}.convolve( obj.virtual_source.source_signal.time_series ) );
            end
            obj.add_output_to_ssd_signal;
        end

    end
end






