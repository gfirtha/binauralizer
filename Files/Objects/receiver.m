classdef receiver < handle
    %RECEIVER Summary of this class goes here
    %   Detailed explanation goes here
    
    events
        position_changed
    end
    properties
        position
        orientation
    end
    
    methods
        function obj = receiver(position, orientation)
            obj.position = position;
            obj.orientation = orientation;
        end

        function obj = set_position(obj, pos_in)
            obj.position = pos_in;
            notify(obj, 'position_changed');
        end
        function obj = set_orientation(obj, orientation_in)
            obj.orientation = orientation_in;
            notify(obj, 'position_changed');
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

