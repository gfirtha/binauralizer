classdef base_renderer < handle
    %BASE_RENDERER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        virtual_source
        ssd
        output_signal
    end
    methods (Abstract)
        update_renderer()
        render()
        get_renderer_info()
    end
    methods
        function obj = base_renderer(virtual_source,SSD)
            obj.virtual_source = virtual_source;
            obj.ssd = SSD;
            for n = 1 : length(SSD.loudspeakers)
                obj.output_signal{n} = signal;
            end
        end
        
        function obj = add_output_to_ssd_signal(obj)
            for n = 1 : length(obj.output_signal)
                obj.ssd.loudspeakers{n}.source_signal.add_signals(...
                    obj.output_signal{n});
            end
        end
    end
end

