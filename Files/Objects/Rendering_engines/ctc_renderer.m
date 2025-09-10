classdef ctc_renderer < base_renderer
    % CTC renderer (cross-talk cancellation).
    %
    % Uses a (2 x Nls) plant transfer matrix (ears x loudspeakers) in the
    % frequency domain, inverts it (Tikhonov) and drives the SSD so that the
    % ears receive the virtual source signal.
    %
    % Renderer_setup fields used (defaults in parentheses):
    %   CTCPlantModel : 'HRTF' | 'point_source' | 'rigid_sphere'   ('HRTF')
    %   CTCVSModel    : 'HRTF' | 'point_source' | 'rigid_sphere'   ('HRTF')
    %   CTCNumTaps    : integer (only used when not 'HRTF')         (1024)
    %   CTCLambda     : Tikhonov regularization (>=0)               (1e-5)
    %   CTCRhead      : head radius (m) for simple models           (0.10)

    properties
        % Scene state / geometry
        receiver                % struct with .position (1x2) and .orientation (1x2 unit)

        % Audio params
        fs
        N_filt

        % Models
        plant_model             % 'HRTF' | 'point_source' | 'rigid_sphere'
        virtual_source_model    % same options
        r_head = 0.10           % m (simple models)
        tikh_lambda = 1e-5

        % HRTF resources (when requested)
        hrtf_database           % SOFA
        hrtf_2d_database        % preprocessed 2D slice (az=theta, elev=0), with FFT spectrum

        % Frequency-domain machinery
        inv_plant_mx_f              % [Nls x 2 x Nf]
        virtual_source_coefficients % [2 x Nf]
    end

    properties (SetAccess=protected)
        CTC_filters             % cell{Nls} of OLS_convolver
    end

    methods
        function obj = ctc_renderer(virtual_source, loudspeaker_array, receiver, setup)
            obj = obj@base_renderer(virtual_source, loudspeaker_array);

            % ---- Basic audio params
            obj.fs = obj.virtual_source.source_signal.fs;

            % ---- Read renderer setup with defaults
            rs = setup;
            obj.plant_model          = getOr(rs,'CTCPlantModel','HRTF');
            obj.virtual_source_model = getOr(rs,'CTCVSModel','HRTF');
            obj.N_filt               = getOr(rs,'CTCNumTaps', 2048);
            obj.tikh_lambda          = max(0, getOr(rs,'CTCLambda', 1e-5));
            obj.r_head               = max(0, getOr(rs,'CTCRhead', 0.10));

            % ---- HRTF (if requested)
            needHRTF = any(strcmpi({obj.plant_model, obj.virtual_source_model}, 'hrtf'));
            if needHRTF
                if isfield(setup,'HRTF_database') && ~isempty(setup.HRTF_database)
                    obj.hrtf_database = setup.HRTF_database;
                    obj.prepare_hrtf_2d();
                    obj.N_filt = size(obj.hrtf_database.Data.IR,3);  % adopt IR length
                else
                    warning('CTC:NoHRTF', ...
                        ['CTC requested HRTF but Setup.HRTF is empty. Falling back to point_source ' ...
                         'for both plant and virtual source models.']);
                    obj.plant_model = 'point_source';
                    obj.virtual_source_model = 'point_source';
                end
            end

            % ---- Receiver (position/orientation in scene)
            obj.receiver = normalize_receiver(receiver);

            % ---- Allocate filters & build once
            obj.rebuild_all();

            % ---------- nested helpers ----------
            function v = getOr(S, f, def)
                if isfield(S,f) && ~isempty(S.(f)), v = S.(f); else, v = def; end
            end
            function R = normalize_receiver(Rin)
                R = Rin;
                if isrow(R.position),    R.position    = R.position(:)';    end
                if isrow(R.orientation), R.orientation = R.orientation(:)'; end
                if numel(R.orientation) >= 2
                    n = norm(R.orientation(1:2));
                    if n>0, R.orientation(1:2) = R.orientation(1:2)/n; end
                end
            end
        end

        function info = get_renderer_info(obj)
            lines = [
                "CTC renderer — overview"
                ""
                "WHAT THIS RENDERER DOES"
                "• Builds an ear-by-loudspeaker plant transfer (2×Nls) in the frequency domain."
                "• Inverts it with Tikhonov regularization, and drives the array so the"
                "  ears receive the desired (binaural) virtual source."
                ""
                "MODELS"
                "• Plant model:        " + string(obj.plant_model)
                "• VS → ear model:     " + string(obj.virtual_source_model)
                "• Taps (N):           " + string(obj.N_filt)
                "• Tikhonov λ:         " + string(obj.tikh_lambda)
                "• r_head (m):         " + string(obj.r_head)
                ""
                "NOTES"
                "• If HRTF is selected, Setup.HRTF provides the database and filters adopt its IR length."
                "• If HRTF is missing, the renderer falls back to point_source models."
            ];
            info = strjoin(lines, newline);
        end

        % ------------------------------------------------------------------
        % Public API – called by your scene wrapper
        % ------------------------------------------------------------------
        function update_renderer(obj, evt)
            % evt: 'receiver_moved'|'receiver_rotated'|'loudspeaker_moved'|'loudspeaker_rotated'
            %      'virtual_source_moved'|'virtual_source_rotated'|'' (default: rebuild both)
            if nargin < 2, evt = ''; end
            switch lower(evt)
                case {'receiver_moved','receiver_rotated'}
                    obj.update_plant_mx();
                    obj.update_vs_model();
                case {'loudspeaker_moved','loudspeaker_rotated'}
                    obj.update_plant_mx();
                case {'virtual_source_moved','virtual_source_rotated'}
                    obj.update_vs_model();
                otherwise
                    obj.update_plant_mx();
                    obj.update_vs_model();
            end
            obj.update_filters();
        end

        function update_renderer_settings(obj, rs)
            % Update from Setup.Renderer_setup without reconstructing scene
            needPlant = false; needVS = false; needLen = false;

            if isfield(rs,'CTCPlantModel') && ~strcmpi(rs.CTCPlantModel, obj.plant_model)
                obj.plant_model = rs.CTCPlantModel; needPlant = true;
            end
            if isfield(rs,'CTCVSModel') && ~strcmpi(rs.CTCVSModel, obj.virtual_source_model)
                obj.virtual_source_model = rs.CTCVSModel; needVS = true;
            end
            if isfield(rs,'CTCNumTaps') && isfinite(rs.CTCNumTaps) && rs.CTCNumTaps>0 ...
                    && ~strcmpi(obj.plant_model,'hrtf') && ~strcmpi(obj.virtual_source_model,'hrtf')
                if obj.N_filt ~= rs.CTCNumTaps, obj.N_filt = rs.CTCNumTaps; needLen = true; end
            end
            if isfield(rs,'CTCLambda') && isfinite(rs.CTCLambda) && rs.CTCLambda>=0
                obj.tikh_lambda = rs.CTCLambda; needPlant = true;
            end
            if isfield(rs,'CTCRhead') && isfinite(rs.CTCRhead) && rs.CTCRhead>=0 %#ok<*STRNU>
                obj.r_head = rs.CTCRhead; needPlant = true; needVS = true;
            end

            % If either model is HRTF and we have a DB, ensure 2D slice prepared and tap-length aligned
            needHRTF = any(strcmpi({obj.plant_model, obj.virtual_source_model}, 'hrtf'));
            if needHRTF && ~isempty(obj.hrtf_database)
                obj.prepare_hrtf_2d();
                obj.N_filt = size(obj.hrtf_database.Data.IR,3);
                needLen = true;
            elseif needHRTF && isempty(obj.hrtf_database)
                warning('CTC:update_renderer_settings','HRTF model requested but no HRTF database on renderer; falling back to point_source.');
                obj.plant_model = 'point_source';
                obj.virtual_source_model = 'point_source';
            end

            if needLen
                obj.rebuild_all();
            else
                if needPlant, obj.update_plant_mx(); end
                if needVS,    obj.update_vs_model();  end
                if needPlant || needVS, obj.update_filters(); end
            end
        end

        function render(obj)
            x = obj.virtual_source.source_signal.time_series;
            for n = 1:numel(obj.CTC_filters)
                obj.output_signal{n}.set_signal( obj.CTC_filters{n}.convolve(x) );
            end
            obj.add_output_to_ssd_signal;
        end
    end

    % ===== Implementation (private) ======================================
    methods (Access=private)
        function rebuild_all(obj)
            % Allocate filters and compute initial coefficients
            Nls = numel(obj.ssd.loudspeakers);
            obj.CTC_filters = cell(1,Nls);

            obj.update_plant_mx();
            obj.update_vs_model();

            drv = obj.get_driving_filter(); % [N x Nls]
            for n = 1:Nls
                obj.CTC_filters{n} = OLS_convolver(drv(:,n), ...
                    length(obj.virtual_source.source_signal.time_series));
            end
        end

        function update_filters(obj)
            drv = obj.get_driving_filter(); % [N x Nls]
            for n = 1:numel(obj.CTC_filters)
                obj.CTC_filters{n}.update_coefficients(drv(:,n));
            end
        end

        function prepare_hrtf_2d(obj)
            % Build a 2D azimuthal slice (elev = 0°) and its FFT
            hrf = obj.hrtf_database;
            obj.N_filt = size(hrf.Data.IR,3);

            R_measurement = mean(hrf.SourcePosition(:,3));
            azel = hrf.SourcePosition(:,1:2);
            ixs  = find(azel(:,2) == 0);             % elevation 0°
            theta = azel(ixs,1)*pi/180;              % radians
            Hmeas = hrf.Data.IR(ixs,:,:);            % [numAz, ears, taps]
            [theta,ord] = sort(theta);
            Hmeas = Hmeas(ord,:,:);

            obj.hrtf_2d_database = struct( ...
                'R', R_measurement, ...
                'theta', theta, ...
                'spectrum', fft(Hmeas, [], 3) );     % [numAz, ears, N_filt]
        end

        function update_plant_mx(obj)
            % Build 2 x Nls x Nf (ears x loudspeakers x freq)
            xs = cell2mat(cellfun(@(x) x.position, obj.ssd.loudspeakers, 'uni', 0)');
            % Two ears around receiver.position; keep original convention
            x_ear = bsxfun(@plus, obj.receiver.position', ...
                fliplr((obj.receiver.orientation' * [1,-1] * obj.r_head)')); % [2 x 2]

            Nls = size(xs,1);
            Nf  = obj.N_filt;
            c   = audioapp.util.Physics.speedOfSound();
            freq = reshape((0:Nf-1)'/Nf*obj.fs, [1,1,Nf]);

            switch lower(obj.plant_model)
                case 'hrtf'
                    plant_mx_t = get_hrtfs(xs, obj.receiver.position, ...
                        obj.receiver.orientation, obj.hrtf_database, obj.hrtf_2d_database);
                    plant_mx_f = fft(plant_mx_t, [], 3); % [ears x Nls x Nf]

                case 'point_source'
                    Rmx = zeros(2, Nls);
                    for n = 1:Nls
                        v_sr = bsxfun(@minus, xs(n,:), x_ear);   % [2 x 2]
                        Rmx(:,n) = sqrt(sum(v_sr.^2,2));
                    end
                    plant_mx_f = 1/(4*pi) * bsxfun(@times, ...
                        exp(-1i * 2*pi * bsxfun(@times, freq, Rmx / c)), 1./Rmx);

                case 'rigid_sphere'
                    Rmx      = zeros(2,Nls);
                    theta_mx = zeros(2,Nls);
                    v_med    = obj.receiver.orientation;
                    for n = 1:Nls
                        v_sr = bsxfun(@minus, xs(n,:), x_ear);
                        Rmx(:,n) = sqrt(sum(v_sr.^2,2));
                        v_ls = (xs(n,:) - obj.receiver.position);
                        vv = norm(v_ls); if vv>0, v_ls = v_ls/vv; end
                        theta_mx(:,n) = -atan2( v_ls(1)*v_med(2) - v_ls(2)*v_med(1), ...
                                                v_ls(1)*v_med(1) + v_ls(2)*v_med(2) );
                    end
                    k  = 2*pi*freq/c;            % [1 1 Nf]
                    Norder = 20;

                    plant_mx_f = zeros(2, Nls, Nf);
                    sign_mx = [ones(1,Nls); -ones(1,Nls)];
                    kR = bsxfun(@times, k, Rmx);            % [2 Nls Nf]
                    A0 = 0.1 * bsxfun(@times, Rmx, 1./(bsxfun(@times, k*obj.r_head^2, exp(-1i*kR))));
                    for n = 0:Norder
                        Pn = getLegendre(n,0,sin(theta_mx));         % [2 Nls]
                        term = (2*n+1) * A0 .* getSphH(n,2,kR) ...
                             .* bsxfun(@times, sign_mx.^n .* Pn, 1./getDifSphH(n,2,k*obj.r_head));
                        plant_mx_f = plant_mx_f + term;
                    end
                    plant_mx_f(~isfinite(plant_mx_f)) = 0;

                otherwise
                    error('Unknown CTC plant model: %s', obj.plant_model);
            end

            % ---- Pseudoinverse with Tikhonov (per frequency)
            obj.inv_plant_mx_f = zeros(Nls, 2, Nf);
            lam = obj.tikh_lambda;
            for kf = 1:Nf
                X = squeeze(plant_mx_f(:,:,kf));       % [2 x Nls]
                obj.inv_plant_mx_f(:,:,kf) = (X'*X + lam*eye(size(X,2))) \ (X');
            end

            % prevent huge gains
            % mx = max(abs(obj.inv_plant_mx_f),[],'all');
            % if mx > 0
            %     obj.inv_plant_mx_f = obj.inv_plant_mx_f / mx;
            % end
        end

        function update_vs_model(obj)
            % Build ear target response (2 x Nf)
            c   = audioapp.util.Physics.speedOfSound();
            Nf  = obj.N_filt;
            f   = (0:Nf-1)/Nf*obj.fs;
            xs  = obj.virtual_source.position;
            v_med = obj.receiver.orientation;

            switch lower(obj.virtual_source_model)
                case 'hrtf'
                    coef_t = get_hrtfs(xs, obj.receiver.position, v_med, ...
                                       obj.hrtf_database, obj.hrtf_2d_database); % [ears x 1 x taps]
                    obj.virtual_source_coefficients = squeeze(fft(coef_t,[],2)); % [2 x Nf]

                case 'point_source'
                    x_ear = bsxfun(@plus, obj.receiver.position', ...
                                   fliplr((v_med' * [1,-1] * obj.r_head)'));
                    R = sqrt(sum(bsxfun(@plus, x_ear, -xs).^2, 2)); % [2 x 1]
                    tmp = 1/(4*pi) * exp(-1i*2*pi*(f(:)' .* (R'/c))) ./ (R');
                    obj.virtual_source_coefficients = tmp.';  % [2 x Nf]
                    s = max(abs(obj.virtual_source_coefficients),[],'all');
                    if s>0, obj.virtual_source_coefficients = obj.virtual_source_coefficients / s; end

                case 'rigid_sphere'
                    x_ear = bsxfun(@plus, obj.receiver.position', ...
                                   fliplr((v_med' * [1,-1] * obj.r_head)'));
                    v_sr = bsxfun(@minus, xs, x_ear);
                    R_vec = sqrt(sum(v_sr.^2,2));             % [2 x 1]
                    v_ls  = (xs - obj.receiver.position);
                    vv = norm(v_ls); if vv>0, v_ls = v_ls/vv; end
                    theta = -atan2(v_ls(1)*v_med(2) - v_ls(2)*v_med(1), ...
                                   v_ls(1)*v_med(1) + v_ls(2)*v_med(2));
                    k = 2*pi*f/c;
                    Norder = 20;
                    VS = zeros(2, numel(f));
                    sign_vec = [1; -1];

                    kR = bsxfun(@times, k, R_vec);            % [2 x Nf]
                    A0 = bsxfun(@times, R_vec, 1./(bsxfun(@times, k*obj.r_head^2, exp(-1i*kR))));
                    for n = 0:Norder
                        Pn = getLegendre(n,0,sin(theta)); % scalar -> broadcast
                        VS = VS + (2*n+1) * A0 .* getSphH(n,2,kR) ...
                             .* bsxfun(@times, sign_vec.^n .* Pn, 1./getDifSphH(n,2,k*obj.r_head));
                    end
                    VS(~isfinite(VS)) = 0;
                    obj.virtual_source_coefficients = VS;

                otherwise
                    error('Unknown CTC VS model: %s', obj.virtual_source_model);
            end
        end

        function driving_filter = get_driving_filter(obj)
            % Per frequency: d(f) [Nls x 1] = invPlant(f) [Nls x 2] * vsCoef(f) [2 x 1]
            Nf  = obj.N_filt;
            Nls = numel(obj.ssd.loudspeakers);
            Df  = zeros(Nls, Nf);
            for kf = 1:Nf
                Df(:,kf) = squeeze(obj.inv_plant_mx_f(:,:,kf)) * ...
                           obj.virtual_source_coefficients(:,kf);
            end
            % IFFT to taps; light window to reduce wrap
            driving_filter = fftshift(ifft(Df.', Nf, 'symmetric'),1) .* tukeywin(Nf,0.1);
        end
    end
end
