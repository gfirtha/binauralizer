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
            % Ensure 1×3
            if numel(position) < 3
                obj.position = [position(:)', zeros(1, 3-numel(position))];
            else
                obj.position = position(1:3);
            end
            if numel(orientation) < 3
                obj.orientation = [orientation(:)', zeros(1, 3-numel(orientation))];
            else
                obj.orientation = orientation(1:3);
            end
        end

        function obj = set_position(obj, pos_in)
            if numel(pos_in) >= 3
                obj.position = pos_in(1:3);
            else
                obj.position(1:numel(pos_in)) = pos_in;
            end
            notify(obj, 'position_changed');
        end
        function obj = set_orientation(obj, orientation_in)
            if numel(orientation_in) >= 3
                obj.orientation = orientation_in(1:3);
            else
                obj.orientation(1:numel(orientation_in)) = orientation_in;
            end
            notify(obj, 'position_changed');
        end
    end
end

