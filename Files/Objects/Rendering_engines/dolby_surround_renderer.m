classdef dolby_surround_renderer < base_renderer
    %VBAP_RENDERER Summary of this class goes here
    %   Detailed explanation goes here

    properties
        G_vec
        mode
        surround_dec_filter
        surround_enc_filter
        bp_filter
        virtual_SSD
        setup
    end

    methods
        function obj = dolby_surround_renderer(virtual_source,SSD)
            obj = obj@base_renderer(virtual_source,SSD);
            obj.mode = 'True_DS';
            obj.virtual_source = virtual_source;
            switch obj.mode
                case 'Virtual_DS'
                    if length(SSD.loudspeakers)~=5
                        error('Channel layout must be 5.1');
                        return;
                    end
                    obj.virtual_SSD = get_default_layout(5,1,'circular');

                case 'True_DS'
                    if length(SSD.loudspeakers)~=2
                        error('Channel layout must be 2.0');
                        return;
                    end
                    obj.virtual_SSD = get_default_layout(5,2,'circular');%TODO
            end
            % Create filter for surround channel
            fs = obj.virtual_source.source_signal.fs;
            Nt = length(obj.virtual_source.source_signal.time_series)/2;
            [b,a] = butter(3,[100,7e3]/(fs/2),'bandpass');
            imp = zeros(Nt,1);
            imp(1) = 1;
            h_bp = filter(b,a,imp);
            h_enc = (ifft(1i*fft(h_bp,2*Nt),2*Nt,'symmetric'));

            obj.surround_enc_filter = OLS_convolver(h_enc, length(obj.virtual_source.source_signal.time_series));

            tau = 20e-3;
            h_dec = [zeros(round(tau*fs),1);h_bp;zeros(Nt-round(tau*fs),1)];
            obj.surround_dec_filter = OLS_convolver(h_dec, length(obj.virtual_source.source_signal.time_series));

            [b,a] = butter(3,[700,1.5e3]/(fs/2),'bandpass');
            obj.bp_filter = [b;a];
            obj.update_renderer;
        end

        function info = get_renderer_info(obj)
            info = "The renderer implements classic Dolby Stereo encoding of the virtual source signal." + ...
                "The encoded stream can be played back in surround using a suitable Dolby Surround capable (typically analog) AVR.";
        end


        function obj = update_renderer(obj, ~)
            obj.G_vec = zeros(length(obj.virtual_SSD),1);
            %L = cell2mat(cellfun( @(x) x.position', obj.virtual_SSD,'UniformOutput',false));
            %R = cell2mat(cellfun( @(x) sqrt(sum(x.position.^2,2)), obj.virtual_SSD,'UniformOutput',false));
            %L = bsxfun(@times, L, 1./R);
            L = obj.virtual_SSD';

            ls_pairs = [(1:length(obj.virtual_SSD));
                circshift((1:length(obj.virtual_SSD)),1)];

            for n = 1 : size(ls_pairs,2)
                g(:,n) = ( L(:,ls_pairs(:,n)) \ obj.virtual_source.position' )/norm(obj.virtual_source.position);
            end
            [~,ix] = max(sum( (g>=0), 1));
            obj.G_vec(ls_pairs(:,ix)) = g(:,ix)/norm(g(:,ix));
            obj.G_vec(4:5) = sum(obj.G_vec(4:5));
        end

        function obj = update_settings(obj,setup)
            obj.setup = setup;
            obj.update_renderer;
        end
        function output = get_renderer_visuals(obj)

        end

        function render(obj)

            for n = 1 : length(obj.virtual_SSD)
                virtual_signals{n} = obj.G_vec(n)*obj.virtual_source.source_signal.time_series ;
            end
            Surr = obj.surround_enc_filter.convolve( virtual_signals{4} );
            Lt = virtual_signals{3} + 1/sqrt(2)*(virtual_signals{2} + Surr);
            Rt = virtual_signals{1} + 1/sqrt(2)*(virtual_signals{2} - Surr);
            switch obj.mode
                case 'Virtual_DS'

                    % Decoding
                    Ldec = Lt;
                    Rdec = Rt;
                    Cdec = Lt + Rt;
                    Sdec = obj.surround_dec_filter.convolve( Lt - Rt );

                    Le = sqrt(sum(filter( obj.bp_filter(1,:), obj.bp_filter(2,:), Lt ).^2));
                    Re = sqrt(sum(filter( obj.bp_filter(1,:), obj.bp_filter(2,:), Rt ).^2));
                    Ce = sqrt(sum(filter( obj.bp_filter(1,:), obj.bp_filter(2,:), Lt+Rt ).^2));
                    Se = sqrt(sum(filter( obj.bp_filter(1,:), obj.bp_filter(2,:),  Lt-ifft(-1i*fft(Rt) ,'symmetric')  ).^2 ) );
                    dir_vec = max ( min( 20*log10( [ Re./Le, Ce./Se ]  ), 40 ), -40 );
                    threshold = 3;
                    Al = 1;
                    Ar = 1;
                    As = 1;
                    Ac = 1;
                    if norm(dir_vec)>threshold
                        if dir_vec(2) > 0 % Front dominance
                            A1 = 1e-2;
                            A0 = ( (Le+Re+Se+Ce)-Se * A1 ) / (Le + Re + Ce);
                            As = A1;
                            Al = A0;
                            Ar = A0;
                            Ac = A0;
                        elseif dir_vec(2) < 0 % surr. dominance
                            A1 = 1e-2;
                            A0 = ( (Le+Re+Se+Ce)-Ce * A1 ) / (Le + Re + Se);
                            Ac = A1;
                            Al = A0;
                            Ar = A0;
                            As = A0;
                        end
                    end
                    output_signal{1} = Rdec*Ar;
                    output_signal{2} = Cdec*Ac;
                    output_signal{3} = Ldec*Al;
                    output_signal{4} = Sdec*As;
                    output_signal{5} = Sdec*As;
                    %%
                    for n = 1 : length(obj.ssd.loudspeakers)
                        obj.output_signal{n}.set_signal( output_signal{n} );
                    end
                case 'True_DS'
                    % 'No decoding'
                    obj.output_signal{1}.set_signal( Lt );
                    obj.output_signal{2}.set_signal( Rt );

            end
            obj.add_output_to_ssd_signal;
        end

    end
end