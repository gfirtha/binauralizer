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
            
            % This call now generates Nls filters
            imp_resps = obj.get_icc_filters( Nfilt, length(obj.SSD));
            
            for n = 1 : length(SSD)
                obj.IACC_filters{n} = OLS_convolver(imp_resps(:,n),Nt);
            end
        end
        
        function render(obj)
            % This loop already works for N channels
            for n = 1 : length(obj.SSD)
                obj.SSD{n}.set_output(obj.IACC_filters{n}.convolve( obj.SSD{n}.source_signal.time_series ) );
            end
        end
        
        function update_icc_filters(obj)
            % This call now gets Nls new filters based on the new alpha
            filts = obj.get_icc_filters(obj.IACC_filters{1}.N_filt, length(obj.SSD) );
            
            % This loop already works for N channels
            for n = 1 : length(obj.IACC_filters)
                obj.IACC_filters{n}.update_coefficients(filts(:,n));
            end
        end
        
        % =================================================================
        % === THIS IS THE MODIFIED HYBRID ALGORITHM =======================
        % =================================================================
        function filters = get_icc_filters(obj, Nf, Nls)
            % Get desired correlation from slider
            alpha = obj.alpha;
            
            % --- N=2 Case: Use original asin() algorithm ---
            if Nls == 2
                sequences = (rand(Nf,Nls))*pi;
                r = alpha;
                k = asin(r)*2/pi;
                phase = zeros(size(sequences));

                if k==0
                    phase(:,1) = sequences(:,1);
                    phase(:,2) = sequences(:,2);
                elseif k == 1
                    phase(:,1) = sequences(:,1);
                    phase(:,2) = sequences(:,1); % Identical
                elseif k == -1
                    phase(:,1) = sequences(:,1);
                    phase(:,2) = sequences(:,1) + pi; % Antiphase
                elseif k>0 && k<1
                    phase(:,1) =  sequences(:,1) + sequences(:,2) ;
                    phase(:,2) =  k*sequences(:,1) + sequences(:,2) ;
                elseif k<0 && k>-1
                    phase(:,1) = sequences(:,1) +sequences(:,2) ;
                    phase(:,2) = - k*(sequences(:,1)-pi) + sequences(:,2) ;
                end
                
                filters = real ( fftshift( ifft( exp(1i*phase), Nf , 1) ) );

            % --- N > 2 Case: Use N-Channel Cholesky Method ---
            elseif Nls > 2
                % Find the mathematical lower limit for N-channel correlation
                lower_limit = -1 / (Nls - 1);
                
                % Clamp the desired alpha to the valid range
                if alpha < lower_limit
                    alpha_clamped = lower_limit;
                else
                    alpha_clamped = alpha;
                end
                
                % Create the Nls x Nls target correlation matrix R
                % R = (1-alpha)*I + alpha*J, where I is identity, J is all-ones
                R = eye(Nls) * (1 - alpha_clamped);
                R = R + alpha_clamped; % Adds alpha_clamped to every element

                % Get the 'coloring' matrix L from Cholesky decomposition
                try
                    L = chol(R, 'lower');
                catch
                    % Add jitter for stability if R is not perfectly 
                    % positive semi-definite due to floating point error
                    R = R + eye(Nls) * 1e-9;
                    L = chol(R, 'lower');
                end
                
                % 1. Create Nls independent Gaussian white noise sequences
                independent_noise = randn(Nf, Nls);
                
                % 2. "Color" the noise using the correlation matrix
                % correlated_noise = independent_noise * L'
                correlated_noise = independent_noise * L'; % Result is Nf x Nls

                % 3. Get the phase spectrum from the correlated noise
                spectrum = fft(correlated_noise, Nf, 1);
                phases = angle(spectrum); % This is Nf x Nls
                
                % 4. Create the all-pass filters from this phase
                filters = real( fftshift( ifft( exp(1i*phases), Nf , 1) ) );
                
            % --- N=1 Case: Just one filter ---
            else
                % Single channel, just create one random phase filter
                phase = (rand(Nf, 1)) * pi;
                filters = real ( fftshift( ifft( exp(1i*phase), Nf , 1) ) );
            end
          %  IACC = xcorr(filters(:,1),filters(:,2),'normalized');
           % [~,ix] = max(abs(IACC));
           % IAAC_act = IACC(ix)
        end
    end
end