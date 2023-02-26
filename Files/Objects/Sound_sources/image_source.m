classdef image_source < handle
    %IMAGE_SOURCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        index
        position
        orientation
        depth
        mother_source
        mother_wall
    end
    
    methods
        function obj = image_source(index, position, orientation, depth, mother_source, mother_wall)
            %IMAGE_SOURCE Construct an instance of this class
            %   Detailed explanation goes here
            obj.index = index;
            obj.position = position;
            obj.orientation = orientation; 
            obj.depth = depth; 
            obj.mother_source = mother_source;
            obj.mother_wall = mother_wall;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

