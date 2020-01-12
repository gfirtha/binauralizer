classdef binaural_source < handle
    %VIRTUAL_SOURCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        source_index
        source_signal
        position
        orientation
        source_type
        hrtf
    end
    
    methods
        function obj = binaural_source(idx, position, orientation, hrtf, type)
            obj.source_index = idx;
            obj.position = position;
            obj.orientation = orientation;
            obj.source_signal = signal;
            obj.hrtf = hrtf;
            obj.source_type = type;
        end
        
        function obj = set_input(obj,input_signal)
            obj.source_signal.set_signal(input_signal);
        end
    end
end

