classdef binaural_source < handle
    %VIRTUAL_SOURCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        source_index
        source_signal
        position
        orientation
        hrtf
    end
    
    methods
        function obj = binaural_source(idx, position, orientation, input_signal, hrtf)
            %VIRTUAL_SOURCE Construct an instance of this class
            %   Detailed explanation goes here
            obj.source_index = idx;
            obj.position = position;
            obj.orientation = orientation;
            obj.source_signal = input_signal;
            obj.hrtf = hrtf;
        end
        
        function obj = set_input(obj,input_signal)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            obj.source_signal = input_signal;
        end
    end
end

