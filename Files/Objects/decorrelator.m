classdef decorrelator
    %DECORRELATOR Summary of this class goes here
    %   Detailed explanation goes here

    properties
        mode
        SSD
        alpha
        IACC_filters
    end

    methods
        function obj = decorrelator(mode,SSD, alpha)
            obj.mode = mode;
            obj.alpha = alpha;
            obj.SSD = SSD;
            Nt = length(SSD{1}.source_signal.time_series);
            Nfilt = Nt/2;
            imp_resps = obj.get_icc_filters( Nfilt, length(obj.SSD));
            for n = 1 : length(SSD)
                obj.IACC_filters{n} = OLS_convolver(imp_resps(:,n),Nt);
            end
        end

        function render(obj)
            for n = 1 : length(obj.SSD)
                obj.SSD{n}.set_output(obj.IACC_filters{n}.convolve( obj.SSD{n}.source_signal.time_series ) );
            end

        end
        function update_icc_filters(obj)
            filts = obj.get_icc_filters(obj.IACC_filters{1}.N_filt, length(obj.SSD) );
            for n = 1 : length(obj.IACC_filters)
                obj.IACC_filters{n}.update_coefficients(filts(:,n));
            end
        end

        function filters = get_icc_filters(obj, Nf, Nls)
            sequences = (rand(Nf,Nls))*pi;
            r = (obj.alpha);
            k = asin(r)*2/pi;

            phase = zeros(size(sequences));
            if k==0
                phase(:,1) = sequences(:,1);
                phase(:,2) = sequences(:,2);
            elseif k == 1
                phase(:,1) = sequences(:,1);
                phase(:,2) = sequences(:,1);
            elseif k == -1
                phase(:,1) = sequences(:,1);
                phase(:,2) = sequences(:,1) + pi;
            elseif k>0 && k<1
                phase(:,1) =  sequences(:,1) + sequences(:,2) ;
                phase(:,2) =  k*sequences(:,1) + sequences(:,2) ;
            elseif k<0 && k>-1
                phase(:,1) = sequences(:,1) +sequences(:,2) ;
                phase(:,2) = - k*(sequences(:,1)-pi) + sequences(:,2) ;
            end
         %   filters = ones(Nf,Nls);
            filters = real ( fftshift( ifft( exp(1i*phase), Nf , 1) ) );
         %  filters =  ( fftshift( ifft(  ones(Nf,2).*[ exp(1i*(1-k)/2*pi), 1] ,[],1,'symmetric')) );
            IACC = xcorr(filters(:,1),filters(:,2),'normalized');
            [~,ix] = max(abs(IACC));
            IAAC_act = IACC(ix);

        end


    end
end

