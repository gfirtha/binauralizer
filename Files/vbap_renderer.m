classdef vbap_renderer < handle
    %VBAP_RENDERER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        virtual_source
        secondary_source_distribution
        output_signal
        G_vec
    end
    
    methods
        function obj = vbap_renderer(virtual_source,SSD)
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
            R = cell2mat(cellfun( @(x) sqrt(sum(x.position.^2,2)), obj.secondary_source_distribution,'UniformOutput',false));
            L = bsxfun(@times, L, 1./R);
            fi = cart2pol(L(1,:),L(2,:));
            [~,i] = sort(abs( cart2pol(obj.virtual_source.position(1),obj.virtual_source.position(2)) - fi ));
            l = L(:,i(1:2));
            ind = i(1:2);
            g = l\obj.virtual_source.position'/norm(obj.virtual_source.position);
            obj.G_vec(ind) = g/norm(obj.virtual_source.position);
            
%            obj.G_vec(ind) = g/sum(g)*norm(obj.virtual_source.position);
        end

        function render(obj)
            for n = 1 : length(obj.output_signal)
                obj.output_signal{n}.set_signal( obj.G_vec(n)*obj.virtual_source.source_signal.time_series );
            end
        end
    end
end