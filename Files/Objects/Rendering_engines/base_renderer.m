classdef base_renderer < handle
    %BASE_RENDERER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        virtual_source
        secondary_source_distribution
        output_signal
    end
    methods (Abstract)
        update_renderer()
        render()
    end
    methods
        function obj = base_renderer(virtual_source,SSD)
            obj.virtual_source = virtual_source;
            obj.secondary_source_distribution = SSD;
            for n = 1 : length(SSD)
                obj.output_signal{n} = signal;
            end
        end
        
        function obj = add_output_to_ssd_signal(obj)
            for n = 1 : length(obj.output_signal)
                obj.secondary_source_distribution{n}.source_signal.add_signals(...
                    obj.output_signal{n});
            end
        end
    end
end

