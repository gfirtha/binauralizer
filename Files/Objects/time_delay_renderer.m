classdef time_delay_renderer < handle
    
    properties
        fs
        c
        virtual_source
        secondary_source_distribution
        renderer_filter
        output_signal
    end
    
    methods
        function obj = time_delay_renderer(virtual_source,SSD, fs)
            obj.fs = fs;
            obj.c = 343.1;
            obj.virtual_source = virtual_source;
            obj.secondary_source_distribution = SSD;
            
            x0 = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
            n0 = cell2mat(cellfun( @(x) x.orientation, obj.secondary_source_distribution, 'UniformOutput', false)');
            delay_filters = get_td_filters( obj.virtual_source.position, x0, n0, obj.fs, obj.c );
            
            for n = 1 : length(obj.secondary_source_distribution)
                obj.output_signal{n} = signal;
                obj.renderer_filter{n} = OLS_convolver(delay_filters(:,n), length(obj.virtual_source.source_signal.time_series));
            end
        end
        
        function obj = update_renderer(obj)
            x0 = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
            n0 = cell2mat(cellfun( @(x) x.orientation, obj.secondary_source_distribution, 'UniformOutput', false)');
            delay_filters = get_td_filters( obj.virtual_source.position, x0, n0, obj.fs, obj.c );
            for n = 1 : length(obj.secondary_source_distribution)
                obj.renderer_filter{n}.update_coefficients(delay_filters(:,n));
            end
        end
        
        function render(obj)
            
            for n = 1 : length(obj.output_signal)
                obj.output_signal{n}.set_signal(  obj.renderer_filter{n}.convolve( obj.virtual_source.source_signal.time_series)  );
            end
        end
    end
end

