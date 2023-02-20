classdef receiver < handle
    %RECEIVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        position
        orientation
    end
    
    methods
        function obj = receiver(position, orientation)
            obj.position = position;
            obj.orientation = orientation;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

