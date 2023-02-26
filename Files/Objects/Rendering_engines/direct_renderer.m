classdef direct_renderer < handle
    %VBAP_RENDERER Summary of this class goes here
    %   Detailed explanation goes here

    properties
        virtual_source
        secondary_source_distribution
        output_signal
        ix
    end

    methods
        function obj = direct_renderer(virtual_source,SSD)
            obj.virtual_source = virtual_source;
            obj.secondary_source_distribution = SSD;
            obj.ix = obj.virtual_source.source_index;
            for n = 1 : length(SSD)
                obj.output_signal{n} = signal;
            end
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
        end
    end
end