classdef nara_renderer < base_renderer

    properties
        fs
        c
        receiver
        focused_flag
        H_pre
        amp
        delay
        WFS_filter_bank
        omega
        antialiasing_filter_bank
    end

    methods
        function obj = nara_renderer(virtual_source,SSD, fs)
            obj = obj@base_renderer(virtual_source,SSD);
            obj.fs = fs;
            obj.c = 343.1;
            obj.receiver = [0, 1];

         %  fileID = fopen('data/lwfs-driving-functions-center.txt','r');
         %  fileID = fopen('data/lwfs-driving-functions-front.txt','r');
        %   fileID = fopen('data/lwfs-driving-functions-right.txt','r');
           fileID = fopen('data/lwfs-driving-functions-left-M45-N60-beta8.6.txt','r');
            d = fscanf(fileID,'%f');
            fclose(fileID);

            fileID = fopen(['data/lwfs-predelay-left-M45-N60-beta8.6.txt'],'r');
           % fileID = fopen(['data/lwfs-predelay-center.txt'],'r');
            delay = fscanf(fileID,'%f');
            fclose(fileID);


            N0 = length(obj.secondary_source_distribution);

            Nt = length(d)/N0;
            d_wfs = reshape(d,[N0,Nt]);

            Nout = length(obj.virtual_source.source_signal.time_series);
            D_wfs = zeros(Nout,N0);
            D_wfs(1:Nt,:) = d_wfs';
            omega = 2*pi*(0:Nout-1)'/Nout*obj.fs;
         %   D_wfs = ifft(fft(D_wfs,[],1).*exp(1i*omega*delay),[],1,'symmetric');
            x0 = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
            R = sqrt(mean(sum(x0.^2,2)));
            ds = 2*pi*R/N0;
            D_wfs = D_wfs.*ds;
            for n = 1 : length(SSD)
                obj.output_signal{n} = signal;
                obj.WFS_filter_bank{n} =  OLS_convolver(zeros(Nout,1),length(obj.virtual_source.source_signal.time_series));
            end
            for n = 1 : length(obj.WFS_filter_bank)
                obj.WFS_filter_bank{n}.update_coefficients(D_wfs(:,n));
            end

        end

        function obj = update_renderer(obj, ~)
        end

        function render(obj)
            for n = 1 : length(obj.WFS_filter_bank)
                obj.output_signal{n}.set_signal( obj.WFS_filter_bank{n}.convolve( obj.virtual_source.source_signal.time_series ) );
            end
            obj.add_output_to_ssd_signal;
        end
    end

end

%nti xl3