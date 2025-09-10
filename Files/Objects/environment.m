classdef environment %< Singleton
    %ENVIRONMENT Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        rho
        c
    end
    
    methods
        function obj = environment
            %ENVIRONMENT Construct an instance of this class
            %   Detailed explanation goes here
            obj.rho = 1.2;
            obj.c = 343.1;
        end
        
    end
end

