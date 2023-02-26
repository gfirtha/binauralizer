classdef dbap_renderer < handle
    %DBAP_RENDERER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        virtual_source
        secondary_source_distribution
        output_signal
        G_vec
    end
    
    methods
        function obj = dbap_renderer(virtual_source,SSD)
            obj.virtual_source = virtual_source;
            obj.secondary_source_distribution = SSD;
            for n = 1 : length(SSD)
                obj.output_signal{n} = signal;
            end
        obj.update_renderer;
        end
        
        function obj = update_renderer(obj)
            obj.G_vec = zeros(length(obj.secondary_source_distribution),1);
            L = cell2mat(cellfun( @(x) x.position', obj.secondary_source_distribution,'UniformOutput',false));
            rs = 0.1;
            di = sqrt(rs^2 + sum(((bsxfun( @minus, L, obj.virtual_source.position')).^2),1));
            a = 6/20*log10(2);
            k = 1/sqrt(sum( 1./di.^(2*a) ));
            g = k./di;
            obj.G_vec = g;
            
%            obj.G_vec(ind) = g/sum(g)*norm(obj.virtual_source.position);
        end

        function render(obj)
            for n = 1 : length(obj.output_signal)
                obj.output_signal{n}.set_signal( obj.G_vec(n)*obj.virtual_source.source_signal.time_series );
            end
        end
    end
end