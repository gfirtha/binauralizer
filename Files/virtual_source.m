classdef virtual_source < handle
    %VIRTUAL_SOURCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        source_index
        source_signal
        position
        orientation
        source_type
        renderer_type
    end
    
    methods
        function obj = virtual_source(idx, position, orientation, source_type, renderer_type)
            %VIRTUAL_SOURCE Construct an instance of this class
            %   Detailed explanation goes here
            obj.source_index = idx;
            obj.position = position;
            obj.source_signal = signal;
            obj.orientation = orientation;
            obj.source_type = source_type;
            obj.renderer_type = renderer_type;
        end

        function obj = set_input(obj,input_signal)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.source_signal.set_signal(input_signal);
        end
        
    end
end

