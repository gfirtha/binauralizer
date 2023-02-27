classdef dbap_renderer < base_renderer
    %DBAP_RENDERER Summary of this class goes here
    %   Detailed explanation goes here

    properties
        G_vec
    end

    methods
        function obj = dbap_renderer(virtual_source,SSD)
            obj = obj@base_renderer(virtual_source,SSD);
            obj.update_renderer;
        end

        function obj = update_renderer(obj, type)
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
            obj.add_output_to_ssd_signal;
        end
    end
end