classdef wfs_renderer < handle
    %WFS_JOB Summary of this class goes here
    %   Detailed explanation goes here
    
    % Input: (set of?) virtual_source_object,
    % Output: n-channel signal object
    properties
        fs
        c
        virtual_source
        secondary_source_distribution
        output_signal
        prefilter
        delay_line
        amp
        delay
    end
    
    methods
        function obj = wfs_renderer(virtual_source,SSD)
            obj.fs = 48e3;
            obj.c = 343.1;
            obj.virtual_source = virtual_source;
            obj.secondary_source_distribution = SSD;
            for n = 1 : length(SSD)
                obj.output_signal{n} = zeros(size(virtual_source.source_signal,1),1);
            end
            h_pre = obj.get_wfs_prefilter;
            obj.prefilter = OLS_convolver(h_pre, length(obj.virtual_source.source_signal));
            
            x0 = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
            n0 = cell2mat(cellfun( @(x) x.orientation, obj.secondary_source_distribution, 'UniformOutput', false)');
            [ obj.amp, obj.delay ] = get_wfs_driving_function( obj.virtual_source.position, x0, n0 );
            
            N_blocks = length(obj.virtual_source.source_signal);
            N = ceil((50/obj.c*obj.fs)/N_blocks);
            obj.delay_line = delay_line(N*N_blocks,N_blocks, N);
        end
        
        function obj = update_driving_function(obj)
            x0 = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
            n0 = cell2mat(cellfun( @(x) x.orientation, obj.secondary_source_distribution, 'UniformOutput', false)');
            [ obj.amp, obj.delay ] = get_wfs_driving_function( obj.virtual_source.position, x0, n0 );            
        end
        
        function h = get_wfs_prefilter(obj)
            N = size(obj.output_signal{1},1);
            w = [(0 : N/2 - 1)';(-N/2:-1)' ]/N*2*pi*obj.fs;
            H = sqrt(1i*w/obj.c);
            H(end/2+1) = real(H(end/2+1));
            h = fftshift(ifft(H)).*tukeywin(N,0.1);
        end
        
        function render(obj)
            obj.delay_line.write( obj.prefilter.convolve(obj.virtual_source.source_signal) );
            for n = 1 : length(obj.output_signal)
                obj.output_signal{n} = obj.amp(n)*...
                    obj.delay_line.read( obj.delay(n),size(obj.virtual_source.source_signal,1) );
            end
        end
    end
end

