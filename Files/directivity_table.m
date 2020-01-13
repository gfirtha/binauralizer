classdef directivity_table < handle
    %DIRECTIVITY_TABLE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        type
        theta
        directivity_mx
    end
    
    methods
        function obj = directivity_table(type,N_fft,fs)
            %DIRECTIVITY_TABLE Construct an instance of this class
            %   Detailed explanation goes here
            obj.type = type;
            N = 90;
            obj.theta = (0:N-1)/N*pi/2;
            switch type.Shape
                case 'point_source'
                    obj.directivity_mx = ones(1,N);
                case 'circular_piston'
                    obj.directivity_mx = get_piston_dir(fs,N_fft,obj.theta,type.R(1));
                case 'two_way_speaker'
                    obj.directivity_mx = get_twoway_dir(fs,N_fft,obj.theta,type.R);
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

