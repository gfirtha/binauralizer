classdef sound_source < handle
    %SOUND_SOURCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        position
        orientation
        directivity
    end
    
    methods
        function obj = sound_source(inputArg1,inputArg2)
            %SOUND_SOURCE Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

