classdef hoa_renderer < base_renderer

    properties
        c
        method
        receiver
  %      HOA_filter_bank
        setup
        Y_ls
        virtual_source_spectrum
    end

    methods
        function obj = hoa_renderer(virtual_source,SSD, N)
            obj = obj@base_renderer(virtual_source,SSD);
            obj.method = 'SAD';
            obj.c = 343.1;
            obj.receiver = [0, 0];
            obj.setup.Norder = N;
            for n = 1 : length(SSD.loudspeakers)
                obj.output_signal{n} = signal;
            end
            obj.update_renderer;
        end
        function info = get_renderer_info(obj)
            info = "The renderer implements Higher Order Ambisonics rendering, currently using SAD decoding strategy by default.";
        end

        function obj = update_renderer(obj, ~)
            x0 = obj.ssd.x0;
            z0 = zeros(size(x0,1),1);
            if size(x0,2) >= 3, z0 = x0(:,3); end
            [azim_ssd, elev_ssd, r_ssd] = cart2sph(x0(:,1), x0(:,2), z0);
            obj.Y_ls = getSpherHarmMx(pi/2-elev_ssd,azim_ssd,obj.setup.Norder,'real');

            xs = obj.virtual_source.position;
            xs_z = 0; if numel(xs) >= 3, xs_z = xs(3); end
            [azim_s, elev_s, r_s] = cart2sph(xs(1), xs(2), xs_z);
            obj.virtual_source_spectrum = 1/obj.setup.Norder*getSpherHarmMx(pi/2-elev_s,azim_s,obj.setup.Norder,'real');
        end

        function obj = update_settings(obj,setup)
            obj.setup.Norder = setup.HOAorder;
            obj.update_renderer;
        end



        function render(obj)
            yout = (obj.Y_ls*(obj.virtual_source.source_signal.time_series*obj.virtual_source_spectrum)')';

            for n = 1 : length(obj.ssd.loudspeakers)
                obj.output_signal{n}.set_signal(yout(:,n));
            end
            obj.add_output_to_ssd_signal;
        end

    end
end

