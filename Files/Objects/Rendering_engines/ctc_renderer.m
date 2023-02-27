classdef ctc_renderer < base_renderer
    %CTC_RENDERER Summary of this class goes here
    %   Detailed explanation goes here
    %
    %
    %   TODO: implement other models for plant_model and
    %   virtual_source_model
    %
    %
    properties
        receiver
        fs
        plant_model
        virtual_source_model    % Virtual source to receiver propagation in frequency domain
        hrtf_database
        hrtf_2d_database
        inv_plant_mx_f    % Inverse plant matrix in frequency domain
        virtual_source_coefficients
        N_filt
        r_head = 0.1 % radius of head for rigid sphere and point source model
    end
    properties (SetAccess = protected)
        CTC_filters
    end

    methods
        function obj = ctc_renderer(varargin)
            obj = obj@base_renderer(varargin{1},varargin{2});
            obj.receiver = varargin{3};
            obj.fs = varargin{4};
            for n = 1 : length(obj.secondary_source_distribution)
                obj.output_signal{n} = signal;
            end
            obj.plant_model = varargin{5};
            obj.virtual_source_model = varargin{6};
            if (strcmp(obj.plant_model,'HRTF') || strcmp(obj.virtual_source_model,'HRTF'))
                obj.hrtf_database = varargin{7};
                obj.N_filt = size(obj.hrtf_database.Data.IR,3);

                R_measurement = mean(obj.hrtf_database.SourcePosition(:,3));
                theta_measurement = obj.hrtf_database.SourcePosition(:,1:2);
                ixs = find(theta_measurement(:,2) == 0);
                theta_measurement = theta_measurement(ixs,1)*pi/180;
                hrtf_measured = obj.hrtf_database.Data.IR(ixs,:,:);
                [theta_measurement,ix] = sort(theta_measurement);
                hrtf_measured = hrtf_measured(ix,:,:);
                %  theta_measurement = theta_measurement(1:16:end);
                %  hrtf_measured = hrtf_measured(1:16:end,:,:);
                obj.hrtf_2d_database = struct('R',R_measurement,'theta',theta_measurement, 'spectrum',fft(hrtf_measured,[],3));

            else
                obj.N_filt = varargin{8};
            end

            obj.update_plant_mx;
            obj.update_vs_model;
            driving_signal = obj.get_driving_filter;
            for n = 1 : length(obj.secondary_source_distribution)
                obj.CTC_filters{n}  = OLS_convolver(driving_signal (:,n), length(obj.virtual_source.source_signal.time_series));
            end

        end


        function obj = update_plant_mx(obj)
            switch obj.plant_model
                % Plant matrix: row index: ear  , column index: loudspeaker
                case 'HRTF'
                    xs = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
                    plant_mx_t = get_hrtfs( xs, obj.receiver.position,...
                        obj.receiver.orientation, obj.hrtf_database, obj.hrtf_2d_database );
                    plant_mx_f = fft(plant_mx_t,[],3);
                    freq = reshape((0:obj.N_filt-1)'/obj.N_filt*obj.fs, [1,1,obj.N_filt] ) ;
                case 'point_source'
                    xs = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
                    x_ear = bsxfun( @plus, obj.receiver.position', fliplr((obj.receiver.orientation'*[1,-1]*obj.r_head)'));
                    Rmx = zeros(size(x_ear,1), size(xs,1));
                    for n = 1 : size(xs,1)
                        v_sr = bsxfun(@minus, xs(n,:), x_ear);
                        Rmx(:,n) = sqrt(sum(v_sr.^2,2));
                    end
                    freq = reshape((0:obj.N_filt-1)'/obj.N_filt*obj.fs, [1,1,obj.N_filt] ) ;
                    plant_mx_f = 1/(4*pi)*bsxfun( @times, exp( -1i*2*pi*bsxfun( @times, freq, Rmx/340  ) ), 1./Rmx);
                case 'rigid_sphere'
                    xs = cell2mat(cellfun( @(x) x.position,    obj.secondary_source_distribution, 'UniformOutput', false)');
                    x_ear = bsxfun( @plus, obj.receiver.position', fliplr((obj.receiver.orientation'*[1,-1]*obj.r_head)'));

                    Rmx = zeros(size(x_ear,1), size(xs,1));
                    theta_mx = zeros(size(x_ear,1), size(xs,1));
                    for n = 1 : size(xs,1)
                        v_sr = bsxfun(@minus, xs(n,:), x_ear);
                        Rmx(:,n) = sqrt(sum(v_sr.^2,2));
                        v_ls = bsxfun(@minus, xs(n,:), obj.receiver.position);
                        v_ls = v_ls/norm(v_ls);
                        v_med = obj.receiver.orientation;
                        theta_mx(:,n) = -atan2(v_ls(1)*v_med(2)-v_ls(2)*v_med(1),v_ls(1)*v_med(1)+v_ls(2)*v_med(2));
                    end
                    c = 343.1;
                    freq = reshape((0:obj.N_filt-1)'/obj.N_filt*obj.fs, [1,1,obj.N_filt] ) ;
                    k = 2*pi*freq / c;
                    Norder = 20;

                    plant_mx_f = zeros(2,size(xs,1),length(freq));
                    sign_mx = [ones(1, size(xs,1));-ones(1, size(xs,1))];

                    kR = bsxfun(@times , k, Rmx );
                    A0 = bsxfun( @times, Rmx, 1./(bsxfun(@times, k*obj.r_head^2, exp(-1i*kR) )));
                    for n = 0 : Norder
                        Pn = getLegendre(n,0,sin(theta_mx));
                        plant_mx_f = plant_mx_f + (2*n+1)*A0.*getSphH( n, 2, kR ).*bsxfun( @times, sign_mx.^n.*Pn, 1./getDifSphH( n, 2, k*obj.r_head ) );
                    end
                    plant_mx_f(isnan(plant_mx_f)) = 0;
                    plant_mx_f(isinf(plant_mx_f)) = 0;

            end


            obj.inv_plant_mx_f = zeros(size(plant_mx_f));

            %%%Tikhonov regularization
            for n = 1 : size(plant_mx_f,3)
                X = squeeze(plant_mx_f(:,:,n));
                lambda = 1e-5;
                obj.inv_plant_mx_f(:,:,n) = pinv(X.'*X + lambda*eye(size(X)))*X.';
            end
            obj.inv_plant_mx_f = 200*obj.inv_plant_mx_f / max(max(max(obj.inv_plant_mx_f)));
%             %%no regularization
%             for n = 1 : size(plant_mx_f,3) 
%                 X = squeeze(plant_mx_f(:,:,n));
%                 obj.inv_plant_mx_f(:,:,n) = pinv(X);
%             end

        end

        function obj = update_vs_model(obj)
            switch obj.virtual_source_model
                case 'HRTF'
                    obj.virtual_source_coefficients = fft(get_hrtfs( obj.virtual_source.position, obj.receiver.position, obj.receiver.orientation, obj.hrtf_database, obj.hrtf_2d_database  ),[],2);
                case 'point_source'
                    xs = obj.virtual_source.position;
                    x_ear = bsxfun( @plus, obj.receiver.position', fliplr((obj.receiver.orientation'*[1,-1]*obj.r_head)'));
                    R = sqrt(sum( (bsxfun( @plus, x_ear, -xs)).^2,2));
                    f = (0:obj.N_filt-1)'/obj.N_filt*obj.fs ;
                    obj.virtual_source_coefficients = 1/(4*pi)* bsxfun(@times, exp( -1i*2*pi* bsxfun(@times, f', R/340 )), 1./R );
                    obj.virtual_source_coefficients = obj.virtual_source_coefficients  / max(max(abs(obj.virtual_source_coefficients )));
                case 'rigid_sphere'
                    xs = obj.virtual_source.position;
                    x_ear = bsxfun( @plus, obj.receiver.position', fliplr((obj.receiver.orientation'*[1,-1]*obj.r_head)'));

                    v_sr = bsxfun(@minus, xs, x_ear);
                    R_vec = sqrt(sum(v_sr.^2,2));
                    v_ls = bsxfun(@minus, xs, obj.receiver.position);
                    v_ls = v_ls/norm(v_ls);
                    v_med = obj.receiver.orientation;
                    theta = -atan2(v_ls(1)*v_med(2)-v_ls(2)*v_med(1),v_ls(1)*v_med(1)+v_ls(2)*v_med(2));
                    % theta = acos(v_ls*v_med');
                    c = 343.1;
                    freq = (0:obj.N_filt-1)/obj.N_filt*obj.fs ;
                    k = 2*pi*freq / c;
                    Norder = 20;
                    VS_coef = zeros(2,length(freq));
                    sign_vec = [1;-1];

                    kR = bsxfun(@times , k, R_vec );
                    A0 = bsxfun( @times, R_vec, 1./(bsxfun(@times, k*obj.r_head^2, exp(-1i*kR) )));
                    for n = 0 : Norder
                        Pn = getLegendre(n,0,sin(theta));
                        VS_coef = VS_coef + (2*n+1)*A0.*getSphH( n, 2, kR ).*bsxfun( @times, sign_vec.^n.*Pn, 1./getDifSphH( n, 2, k*obj.r_head ) );
                    end
                    VS_coef(isnan(VS_coef)) = 0;
                    VS_coef(isinf(VS_coef)) = 0;
                    obj.virtual_source_coefficients = VS_coef;
            end
%             figure('Name',obj.virtual_source_model)
%             semilogx(squeeze((0:obj.N_filt-1)/obj.N_filt*obj.fs),20*log10(squeeze(abs(obj.virtual_source_coefficients(1,:)))))
%             hold on
%             semilogx(squeeze((0:obj.N_filt-1)/obj.N_filt*obj.fs),20*log10(squeeze(abs(obj.virtual_source_coefficients(2,:)))))
%             grid on
%             ylabel({'Magnitude [dB]'});
%             xlabel({'frequency [Hz]'});
%             xlim([200,20e3])
        end

        function obj = update_renderer(obj,type)
            switch type
                case 'receiver_moved'
                    obj.update_plant_mx;
                    obj.update_vs_model;
                case 'receiver_rotated'
                    obj.update_plant_mx;
                    obj.update_vs_model;
                case 'loudspeaker_moved'
                    obj.update_plant_mx;
                case 'loudspeaker_rotated'
                    obj.update_plant_mx;
                case 'virtual_source_moved'
                    obj.update_vs_model;
                case 'virtual_source_rotated'
                    obj.update_vs_model;
            end
            driving_signal = obj.get_driving_filter;
            for n = 1 : length(obj.secondary_source_distribution)
                obj.CTC_filters{n}.update_coefficients(driving_signal (:,n));
            end

        end

        function driving_filter = get_driving_filter(obj)
            N = size(obj.inv_plant_mx_f,3);
            driving_function = zeros(N, size(obj.inv_plant_mx_f,2));
            for n = 1 : N
                driving_function(n,:) = squeeze(obj.inv_plant_mx_f(:,:,n))*obj.virtual_source_coefficients(:,n);
            end
            driving_filter = fftshift(ifft(driving_function,N,1,'symmetric'),1).*tukeywin(N,0.1);
        end

        function render(obj)
            for n = 1 : length(obj.secondary_source_distribution)
                obj.output_signal{n}.set_signal( obj.CTC_filters{n}.convolve( obj.virtual_source.source_signal.time_series ) );
            end
            obj.add_output_to_ssd_signal;
        end

    end
end

