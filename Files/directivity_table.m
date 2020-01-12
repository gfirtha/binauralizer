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
            N = 180;
            obj.theta = (0:N-1)/N*2*pi;
            switch type.Shape
                case 'point_source'
                    obj.directivity_mx = ones(1,N);
                case 'circular_piston'
                    w = [(0 : N_fft/2 - 1)';(-N_fft/2:-1)' ]/N_fft*2*pi*fs;
                    c = 343.1;
                    k = w / c;
                    r0 = type.R;
                    [Theta, K] = meshgrid(obj.theta,k);
                    obj.directivity_mx = zeros(1,N);
                    D = 2*besselj(1,r0*K.*sin(Theta))./(r0*K.*sin(Theta)).*exp(-1i*K*c*N/fs/2);
                    D(1,:) = 1;
                    D(:,1) = 1;
                    D(end/2+1,:) = real(D(end/2+1,:));
                    obj.directivity_mx = D;
            end
        end
        
        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

