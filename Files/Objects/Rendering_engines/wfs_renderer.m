classdef wfs_renderer < base_renderer

    properties
        fs
        c
        source_directivity
        focused_flag
        prefilter
        delay_line
        amp
        delay
        antialiasing
        antialiasing_filters
    end

    methods
        function obj = wfs_renderer(virtual_source,SSD, fs, directivity, antialiasing)
            obj = obj@base_renderer(virtual_source,SSD);
            obj.antialiasing = antialiasing;
            obj.fs = fs;
            obj.c = 343.1;
            obj.source_directivity = directivity;
            for n = 1 : length(SSD)
                obj.output_signal{n} = signal;
            end
            x0 = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
            n0 = cell2mat(cellfun( @(x) x.orientation, obj.secondary_source_distribution, 'UniformOutput', false)');
            [ obj.amp, obj.delay, obj.focused_flag, aliasing_filter ] = get_wfs_driving_function( obj.virtual_source.position, obj.virtual_source.source_type.Shape, x0, n0, obj.fs, obj.c, obj.antialiasing );

            h_pre = obj.get_wfs_prefilter(obj.focused_flag);
            obj.prefilter = OLS_convolver(h_pre, length(obj.virtual_source.source_signal.time_series));

            N_blocks = length(obj.virtual_source.source_signal.time_series);
            N = ceil((100/obj.c*obj.fs)/N_blocks);
            obj.delay_line = delay_line(N*N_blocks,N_blocks, N);
            switch obj.antialiasing
                case 'on'
                    for n = 1 : length(obj.secondary_source_distribution)
                        obj.antialiasing_filters{n} = OLS_convolver(aliasing_filter(:,n), length(obj.virtual_source.source_signal.time_series));
                    end
                case 'off'
            end
        end

        function obj = update_renderer(obj, ~)
            x0 = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
            n0 = cell2mat(cellfun( @(x) x.orientation, obj.secondary_source_distribution, 'UniformOutput', false)');
            [ obj.amp, obj.delay, focused, aliasing_filter ] = get_wfs_driving_function( obj.virtual_source.position, obj.virtual_source.source_type.Shape, x0, n0, obj.fs, obj.c, obj.antialiasing );
            if focused ~= obj.focused_flag
                obj.prefilter.update_coefficients(obj.get_wfs_prefilter(focused));
                obj.focused_flag = focused;
            end
            switch obj.antialiasing
                case 'on'
                    for n = 1 : length(obj.secondary_source_distribution)
                        obj.antialiasing_filters{n}.update_coefficients(aliasing_filter(:,n));
                    end
                case 'off'
            end
        end


        function h = get_wfs_prefilter(obj,focused)
            Nfilt = 128;
            method = 'windowing';
            switch method
                case 'windowing'
                    N = size(obj.virtual_source.source_signal.time_series,1);
                    w = [(0 : N/2 - 1)';(-N/2:-1)' ]/N*2*pi*obj.fs;
                    H = sqrt(-focused*1i*w/obj.c);
                    H(end/2+1) = real(H(end/2+1));
                    h = fftshift(ifft(H));
                    h = h(round(end/2)-Nfilt/2+1:round(end/2)+Nfilt/2).*hann(Nfilt);
                    h = h/sum(h);
                case 'optimal'
                    n = (-Nfilt/2:Nfilt/2-1)';
                    h = (sqrt(8*n).*cos(pi*n)+ ...
                        erfz( exp(1i*3/4*pi).*sqrt(pi*n) )-erfz( exp(1i*pi/4).*sqrt(pi*n) )).*...
                        (4*sqrt(pi)*n.^(3/2)).^(-1);
                    h(end/2+1) = sqrt(2*pi/9);
                    h = kaiser(Nfilt,4).*real(h);
                    h = h/sum(h);
            end
        end

        function render(obj)
            obj.delay_line.write( obj.prefilter.convolve( obj.virtual_source.source_signal.time_series) );
            switch obj.antialiasing
                case 'on'
                    for n = 1 : length(obj.output_signal)
                        obj.output_signal{n}.set_signal( obj.amp(n)*...
                            obj.antialiasing_filters{n}.convolve( ...
                            obj.delay_line.read( obj.delay(n),size(obj.virtual_source.source_signal.time_series,1) ) ));
                    end
                case 'off'
                    for n = 1 : length(obj.output_signal)
                        obj.output_signal{n}.set_signal( obj.amp(n)*...
                            obj.delay_line.read( obj.delay(n),size(obj.virtual_source.source_signal.time_series,1) ) );
                    end
            end
            obj.add_output_to_ssd_signal;
        end
    end
end

