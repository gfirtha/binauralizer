classdef loudspeaker < handle
    %VIRTUAL_SOURCE Summary of this class goes here
    %   Detailed explanation goes here
    events
        position_changed
        height_changed
    end
    properties
        source_index
        source_signal
        position
        height
        orientation
        source_type
    end
    
    methods
        function obj = loudspeaker(idx, position, orientation, type, height)
            obj.source_index = idx;
            obj.position = position;
            obj.height = height;
            obj.orientation = orientation;
            obj.source_signal = signal;
            obj.source_type = type;
        end
        function obj = set_height(obj, height)
            obj.height = height;
            notify(obj, 'height_changed');
        end
        function obj = set_position(obj, pos_in)
            obj.position = pos_in;
            notify(obj, 'position_changed');
        end
        function obj = set_orientation(obj, orientation_in)
            obj.orientation = orientation_in;
            notify(obj, 'position_changed');
        end
        
        function obj = set_output(obj,varargin)
            if length(varargin) == 1
                obj.source_signal.set_signal(varargin{1});
            else
                obj.source_signal.set_signal(varargin{1},varargin{2});
            end
        end
    end
end

