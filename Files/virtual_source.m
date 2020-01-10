classdef virtual_source < handle
    %VIRTUAL_SOURCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        source_index
        source_signal
        position
        orientation
    end
    
    methods
        function obj = virtual_source(idx, position, orientation)
            %VIRTUAL_SOURCE Construct an instance of this class
            %   Detailed explanation goes here
            obj.source_index = idx;
            obj.position = position;
            obj.source_signal = signal;
            obj.orientation = orientation;
        end
        
        function obj = set_input(obj,input_signal)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.source_signal.set_signal(input_signal);
        end
    end
end

