classdef direct_renderer < base_renderer
    %VBAP_RENDERER Summary of this class goes here
    %   Detailed explanation goes here

    properties
        ix
    end

    methods
        function obj = direct_renderer(virtual_source,SSD)
            obj = obj@base_renderer(virtual_source,SSD);
            obj.ix = obj.virtual_source.source_index;
            obj.update_renderer;
        end

        function obj = update_renderer(obj,~)
            obj.virtual_source.position = obj.secondary_source_distribution{obj.ix}.position ;
            obj.virtual_source.orientation = obj.secondary_source_distribution{obj.ix}.orientation;
        end

        function render(obj)
            obj.output_signal{obj.ix}.set_signal( obj.virtual_source.source_signal.time_series );
            ix_ = setdiff((1:length(obj.output_signal)),obj.ix);
            for n = 1 : length(ix_)
                obj.output_signal{ix_(n)}.set_signal( zeros(length(obj.virtual_source.source_signal.time_series),1) );
            end
            obj.add_output_to_ssd_signal;
        end
    end
end