classdef beamformer_renderer < base_renderer
    % BEAMFORMER_RENDERER  Dual-beam CTC renderer backed by internal WFS.
    %
    % Two independent WFS beams (L and R) target separate foci to create
    % crosstalk-cancelled listening zones.
    %
    % Modes:
    %   'beam_l'  – mono input → left beam only   (audition/positioning)
    %   'beam_r'  – mono input → right beam only   (audition/positioning)
    %   'ctc'     – N VS × HRTF → earL/earR → beamL/beamR (full CTC)

    properties
        fs
        wfs_mode     = 'das'     % 'das' | 'focus' | 'directive' | 'lsfi' | 'race' | 'multiband'
        N_taps       = 256
        receiver                 % receiver object

        % LSFI beamformer settings (Hohnerlein & Ahrens method)
        lsfi_null_width   = 9    % null width around stop angle (degrees)
        lsfi_null_depth   = -40  % null depth constraint (dB), e.g. -40
        lsfi_fir_length   = 1024 % FIR filter length for LSFI
        lsfi_angle_span   = 80   % design angles from -span..+span around broadside
        lsfi_angle_step   = 1    % angular grid resolution (degrees)
        lsfi_head_diameter = 0.175 % ear-to-ear spacing used for auto LSFI angles
        lsfi_main_angle_L  = NaN  % deg, NaN = auto from left ear position
        lsfi_main_angle_R  = NaN  % deg, NaN = auto from right ear position
        lsfi_zero_angle_L  = NaN  % deg, NaN = auto from right ear position
        lsfi_zero_angle_R  = NaN  % deg, NaN = auto from left ear position
        lsfi_filters_L    = {}   % cell{Nls} of OLS_convolver for L beam
        lsfi_filters_R    = {}   % cell{Nls} of OLS_convolver for R beam
        lsfi_fir_L        = []   % [Nfir × Nls] time-domain FIR coefficients (L beam)
        lsfi_fir_R        = []   % [Nfir × Nls] time-domain FIR coefficients (R beam)

        % RACE — Recursive Ambiophonic Crosstalk Elimination
        % (Glasgal 2007; Hohnerlein & Ahrens ICASSP 2017, low-band component)
        % Self-standing version: only the two outermost loudspeakers are
        % active; the 6 middle elements stay silent. The 2-channel symmetric
        % heuristic does not generalise to N>2 — that's by design.
        race_delta_a       = 0.5    % attenuation factor (0..1), auto-computed
        race_delta_t       = 6e-4   % delay (s), auto-computed
        race_fft_length    = 1024   % FFT length used to build the 2x2 filters
        race_idx_L         = 1      % SSD index of the leftmost loudspeaker
        race_idx_R         = 1      % SSD index of the rightmost loudspeaker
        race_filter_LL     = []     % OLS_convolver: in_L -> out_L (direct)
        race_filter_LR     = []     % OLS_convolver: in_R -> out_L (cross)
        race_filter_RL     = []     % OLS_convolver: in_L -> out_R (cross)
        race_filter_RR     = []     % OLS_convolver: in_R -> out_R (direct)
        race_H_direct      = []     % [Nf × 1] complex spectrum, 1/(1-G^2)
        race_H_cross       = []     % [Nf × 1] complex spectrum, -G/(1-G^2)
        race_calibration_freq = 1000 % Hz, frequency at which Δa is calibrated
                                     % from rigid-sphere ILD

        % MULTIBAND — RACE (low band) + LSFI (high band) crossover
        % Linear-phase complementary FIR pair: h_lp + h_hp = δ exactly,
        % so reconstruction is perfect when summed at the loudspeakers.
        mb_crossover_hz   = 1000  % cross-over frequency (Hz), user-adjustable
        mb_filter_length  = 1024  % length of LP/HP filters (samples)
        mb_h_lp           = []    % [Nfft × 1] low-pass FIR coefficients
        mb_h_hp           = []    % [Nfft × 1] high-pass FIR coefficients
        mb_conv_lp_L      = []    % OLS_convolver: LP applied to bL pre-RACE
        mb_conv_lp_R      = []    % OLS_convolver: LP applied to bR pre-RACE
        mb_conv_hp_L      = []    % OLS_convolver: HP applied to bL pre-LSFI
        mb_conv_hp_R      = []    % OLS_convolver: HP applied to bR pre-LSFI

        % Two beam positions (scene coordinates)
        beamL_pos    = [0 0 0]   % left-ear beam VS position
        beamR_pos    = [0 0 0]   % right-ear beam VS position
        steerAngleL  = 0         % degrees, for display
        steerAngleR  = 0
        steerAngle   = 0         % = steerAngleL, kept for BeamformerPainter fallback

        % CTC operating mode
        ctc_mode     = 'beam_l'  % 'beam_l' | 'beam_r' | 'ctc'

        % Internal WFS renderers (both share the same SSD)
        wfsL         = []
        wfsR         = []

        % Extracted weights for visualization (from L beam)
        weights

        % WFS setup forwarded to both internal renderers
        wfs_setup    = struct( ...
            'Tapering',      0.3, ...
            'Antialiasing',  false, ...
            'ReferenceMode', 'fixed_point', ...
            'RefPoint',      [0 0 0], ...
            'RefDistance',    1, ...
            'DeltaY',        0.5, ...
            'Rref',          1)

        % Speaker equalization (post-mix)
        SpeakerEQEnabled = false
        eqConvolvers     = {}
        eqImpulseResponses = []   % [Ntaps x Nls] time-domain speaker EQ FIRs

        % Low-frequency pre-eq shelving (passed through to wfsL/wfsR)
        LowFreqCut        = false
        ear_receiver_L    = []   % receiver proxy at left  ear — set by panel
        ear_receiver_R    = []   % receiver proxy at right ear — set by panel

        % HRTF extrapolator handle — set from outside for CTC mode
        hrtf_extrap  = []
        % OLS convolvers for HRTF per VS (built on first CTC render)
        hrtf_convL   = {}
        hrtf_convR   = {}

        % Full multi-channel input (set by scene_renderer before render)
        fullInput    = []

        % Scene virtual source handles (set by scene_renderer) — used by CTC
        % to look up HRTFs at each source's actual position.
        scene_virtual_sources = {}

        % Marks the "primary" beamformer instance for CTC rendering.
        % Only the primary renderer does the full CTC mix; the others
        % skip rendering to avoid summing into the SSD twice.
        ctc_primary = false

        % CTC plant inversion
        CTCInverseEnabled = false    % apply plant inverse on ear signals
        ear_transfer_model = 'hrtf'  % 'analytic' | 'hrtf' — plant model
        ctc_beta = 0.05              % Tikhonov regularisation parameter
        ctc_inv_filters = []         % struct with 4 OLS_convolvers: LL,LR,RL,RR
    end

    methods
        function obj = beamformer_renderer(virtual_source, SSD, setup, receiver)
            obj = obj@base_renderer(virtual_source, SSD);
            obj.fs = obj.virtual_source.source_signal.fs;
            if nargin >= 4 && ~isempty(receiver)
                obj.receiver = receiver;
                rp = receiver.position(:)';
                if numel(rp) < 3, rp = [rp, zeros(1, 3-numel(rp))]; end
                obj.wfs_setup.RefPoint = rp;
                obj.wfs_setup.ReferenceMode = 'fixed_point';
            end

            % Read setup
            RS = struct();
            if isstruct(setup) && isfield(setup, 'Renderer_setup')
                RS = setup.Renderer_setup;
            elseif isstruct(setup)
                RS = setup;
            end
            if isfield(RS, 'BeamformerMode'),    obj.wfs_mode  = lower(RS.BeamformerMode); end
            if isfield(RS, 'BeamformerNumTaps'), obj.N_taps    = RS.BeamformerNumTaps; end
            if isfield(RS, 'CTCMode'),           obj.ctc_mode  = lower(RS.CTCMode); end
            if isfield(RS, 'LSFINullWidth'),     obj.lsfi_null_width = RS.LSFINullWidth; end
            if isfield(RS, 'LSFINullDepth'),     obj.lsfi_null_depth = RS.LSFINullDepth; end
            if isfield(RS, 'LSFIFIRLength'),     obj.lsfi_fir_length = RS.LSFIFIRLength; end
            if isfield(RS, 'LSFIHeadDiameter'),  obj.lsfi_head_diameter = RS.LSFIHeadDiameter; end
            if isfield(RS, 'LSFIMainAngleL'),    obj.lsfi_main_angle_L = RS.LSFIMainAngleL; end
            if isfield(RS, 'LSFIMainAngleR'),    obj.lsfi_main_angle_R = RS.LSFIMainAngleR; end
            if isfield(RS, 'LSFIZeroAngleL'),    obj.lsfi_zero_angle_L = RS.LSFIZeroAngleL; end
            if isfield(RS, 'LSFIZeroAngleR'),    obj.lsfi_zero_angle_R = RS.LSFIZeroAngleR; end
            if isfield(RS, 'LSFINullAngleL'),    obj.lsfi_zero_angle_L = RS.LSFINullAngleL; end
            if isfield(RS, 'LSFINullAngleR'),    obj.lsfi_zero_angle_R = RS.LSFINullAngleR; end
            if isfield(RS, 'BeamLPos'),          obj.beamL_pos = RS.BeamLPos; end
            if isfield(RS, 'BeamRPos'),          obj.beamR_pos = RS.BeamRPos; end

            % Forward WFS-specific fields
            if isfield(RS, 'Tapering'),      obj.wfs_setup.Tapering      = RS.Tapering; end
            if isfield(RS, 'Antialiasing'),  obj.wfs_setup.Antialiasing  = RS.Antialiasing; end
            if isfield(RS, 'ReferenceMode'), obj.wfs_setup.ReferenceMode = RS.ReferenceMode; end
            if isfield(RS, 'RefPoint'),      obj.wfs_setup.RefPoint      = RS.RefPoint; end
            if isfield(RS, 'RefDistance'),   obj.wfs_setup.RefDistance   = RS.RefDistance; end

            % Default beam positions: symmetric ±15 cm around VS
            if all(obj.beamL_pos == 0) && all(obj.beamR_pos == 0)
                vs = obj.virtual_source.position(:)';
                if numel(vs) < 3, vs = [vs, zeros(1, 3-numel(vs))]; end
                obj.beamL_pos = vs + [-0.15, 0, 0];
                obj.beamR_pos = vs + [+0.15, 0, 0];
            end

            obj.update_renderer();
        end

        % ═══════════════════════════════════════════════
        %  Info
        % ═══════════════════════════════════════════════
        function info = get_renderer_info(obj)
            if strcmpi(obj.wfs_mode, 'race')
                Nls = numel(obj.ssd.loudspeakers);
                lines = [
                    "Beamformer CTC renderer — RACE mode"
                    "(Glasgal 2007; Hohnerlein & Ahrens ICASSP 2017, low band)"
                    ""
                    sprintf("WFS mode:        %s", upper(obj.wfs_mode))
                    sprintf("CTC mode:        %s", obj.ctc_mode)
                    sprintf("Active speakers: %d / %d  (idx %d, %d — outermost)", ...
                        2, Nls, obj.race_idx_L, obj.race_idx_R)
                    sprintf("delta_a:         %.3f  (auto, rigid sphere @ %.0f Hz)", ...
                        obj.race_delta_a, obj.race_calibration_freq)
                    sprintf("delta_t:         %.4f ms  (auto, geometric)", ...
                        obj.race_delta_t * 1000)
                    sprintf("FFT length:      %d", obj.race_fft_length)
                ];
                info = strjoin(lines, newline);
                return;
            end

            if strcmpi(obj.wfs_mode, 'multiband')
                Nls = numel(obj.ssd.loudspeakers);
                fc_opt = obj.estimate_optimal_crossover();
                lines = [
                    "Beamformer CTC renderer — MULTIBAND mode"
                    "(Hohnerlein & Ahrens ICASSP 2017, RACE+LSFI combination)"
                    ""
                    sprintf("WFS mode:         %s", upper(obj.wfs_mode))
                    sprintf("CTC mode:         %s", obj.ctc_mode)
                    sprintf("Crossover:        %.0f Hz  (default 1000 Hz; optimal ≈ %.0f Hz)", ...
                        obj.mb_crossover_hz, fc_opt)
                    sprintf("FIR length:       %d  (linear-phase complementary pair)", ...
                        obj.mb_filter_length)
                    ""
                    sprintf("Low band → RACE   (speakers %d, %d only)", ...
                        obj.race_idx_L, obj.race_idx_R)
                    sprintf("    delta_a:      %.3f", obj.race_delta_a)
                    sprintf("    delta_t:      %.4f ms", obj.race_delta_t * 1000)
                    ""
                    sprintf("High band → LSFI  (all %d speakers)", Nls)
                    sprintf("    null width:   %.0f deg", obj.lsfi_null_width)
                    sprintf("    null depth:   %.0f dB",  obj.lsfi_null_depth)
                    sprintf("    FIR length:   %d",       obj.lsfi_fir_length)
                ];
                info = strjoin(lines, newline);
                return;
            end

            if strcmpi(obj.wfs_mode, 'lsfi')
                lines = [
                    "Beamformer CTC renderer - LSFI mode"
                    "(angular constrained least-squares beamforming)"
                    ""
                    sprintf("WFS mode:        %s", upper(obj.wfs_mode))
                    sprintf("CTC mode:        %s", obj.ctc_mode)
                    sprintf("Null width:      %.0f deg", obj.lsfi_null_width)
                    sprintf("Null depth:      %.0f dB", obj.lsfi_null_depth)
                    sprintf("FIR length:      %d", obj.lsfi_fir_length)
                    sprintf("Head diameter:   %.3f m", obj.lsfi_head_diameter)
                    sprintf("L main / zero:   %s / %s deg", ...
                        obj.lsfi_angle_label(obj.lsfi_main_angle_L), ...
                        obj.lsfi_angle_label(obj.lsfi_zero_angle_L))
                    sprintf("R main / zero:   %s / %s deg", ...
                        obj.lsfi_angle_label(obj.lsfi_main_angle_R), ...
                        obj.lsfi_angle_label(obj.lsfi_zero_angle_R))
                    sprintf("Angle grid:      +/-%.0f deg, %.2f deg step", ...
                        obj.lsfi_angle_span, obj.lsfi_angle_step)
                    sprintf("Elements:        %d", numel(obj.ssd.loudspeakers))
                ];
                info = strjoin(lines, newline);
                return;
            end

            lines = [
                "Beamformer CTC renderer (dual WFS beams)"
                ""
                sprintf("WFS mode:        %s", upper(obj.wfs_mode))
                sprintf("CTC mode:        %s", obj.ctc_mode)
                sprintf("Beam L angle:    %.1f deg", obj.steerAngleL)
                sprintf("Beam R angle:    %.1f deg", obj.steerAngleR)
                sprintf("FIR taps:        %d", obj.N_taps)
                sprintf("Elements:        %d", numel(obj.ssd.loudspeakers))
                ""
                sprintf("Tapering:        %.0f%%", obj.wfs_setup.Tapering*100)
                sprintf("AntiAliasing:    %s", string(obj.wfs_setup.Antialiasing))
            ];
            info = strjoin(lines, newline);
        end

        function fc_opt = estimate_optimal_crossover(obj)
            % Estimate the optimal RACE/LSFI crossover frequency for the
            % current geometry. The LSFI lower robust bound is set by the
            % spatial Nyquist of the array: roughly f = c / (2·Δx).
            % Below that, the array's superdirectivity collapses (WNG → 0).
            % Above that, RACE adds little because the head shadow itself
            % already provides separation.
            try
                x0 = obj.ssd.x0;
                spread_x = max(x0(:,1)) - min(x0(:,1));
                spread_y = max(x0(:,2)) - min(x0(:,2));
                if spread_y > spread_x, axis_col = 2; else, axis_col = 1; end
                xs = sort(x0(:, axis_col));
                d_elem = mean(diff(xs));     % mean inter-element spacing
                c = audioapp.util.Physics.speedOfSound();
                fc_opt = c / (2 * max(d_elem, 0.005));
            catch
                fc_opt = 1000;
            end
        end

        % ═══════════════════════════════════════════════
        %  Update / settings
        % ═══════════════════════════════════════════════
        function obj = update_renderer(obj, ~)
            x0 = obj.ssd.x0;
            center = mean(x0, 1);

            % Steering angles (for display)
            dirL = obj.beamL_pos(1:2) - center(1:2);
            obj.steerAngleL = atan2d(dirL(2), dirL(1));
            dirR = obj.beamR_pos(1:2) - center(1:2);
            obj.steerAngleR = atan2d(dirR(2), dirR(1));

            % Keep legacy steerAngle for backward compat
            obj.steerAngle = obj.steerAngleL;

            % Compute ear-position proxy receivers (always, not only when LowCut is on)
            obj.update_ear_receivers();

            % Invalidate stale caches BEFORE rebuilding
            obj.hrtf_convL = {};
            obj.hrtf_convR = {};
            obj.ctc_inv_filters = [];
            obj.lsfi_filters_L = {};
            obj.lsfi_filters_R = {};
            obj.lsfi_fir_L = [];
            obj.lsfi_fir_R = [];
            obj.race_filter_LL = [];
            obj.race_filter_LR = [];
            obj.race_filter_RL = [];
            obj.race_filter_RR = [];
            obj.race_H_direct  = [];
            obj.race_H_cross   = [];
            obj.mb_h_lp = [];
            obj.mb_h_hp = [];
            obj.mb_conv_lp_L = [];
            obj.mb_conv_lp_R = [];
            obj.mb_conv_hp_L = [];
            obj.mb_conv_hp_R = [];

            % Build both beams (WFS-based, LSFI, RACE, or MULTIBAND)
            if strcmpi(obj.wfs_mode, 'lsfi')
                obj.build_beam_lsfi('L');
                obj.build_beam_lsfi('R');
            elseif strcmpi(obj.wfs_mode, 'race')
                obj.build_beam_race();
            elseif strcmpi(obj.wfs_mode, 'multiband')
                % Multiband mode = RACE (low band) + LSFI (high band)
                % combined via a complementary linear-phase FIR crossover.
                obj.build_beam_race();
                obj.build_beam_lsfi('L');
                obj.build_beam_lsfi('R');
                obj.build_multiband_crossover();
            else
                obj.build_beam_wfs('L');
                obj.build_beam_wfs('R');
                % For directive mode: apply null-steering correction
                % to the WFS driving functions after they are built.
                if strcmpi(obj.wfs_mode, 'directive')
                    obj.apply_directive_correction();
                end
            end

            % Extract weights from L beam for visualization
            if ~isempty(obj.wfsL)
                obj.weights = abs(obj.wfsL.amp) .* obj.wfsL.tapering_window;
            end
        end

        function update_ear_receivers(obj)
            % Derive left/right ear proxy receivers from the scene receiver.
            if isempty(obj.receiver), return; end
            rec_pos = obj.receiver.position(:).';
            rec_ori = obj.receiver.orientation(:).';
            if numel(rec_pos) < 3, rec_pos = [rec_pos, zeros(1,3-numel(rec_pos))]; end
            if numel(rec_ori) < 3, rec_ori = [rec_ori, zeros(1,3-numel(rec_ori))]; end
            n_hat3 = rec_ori / max(norm(rec_ori), 1e-15);
            perp3 = cross(n_hat3, [0 0 1]);
            if norm(perp3) < 1e-10, perp3 = cross(n_hat3, [1 0 0]); end
            perp3 = perp3 / norm(perp3);
            % perp3 points to the right (+X when facing +Y); left ear is -perp3
            half_head = 0.5 * max(double(obj.lsfi_head_diameter), 0);
            obj.ear_receiver_L = receiver(rec_pos - half_head*perp3, rec_ori);
            obj.ear_receiver_R = receiver(rec_pos + half_head*perp3, rec_ori);
        end

        function obj = update_settings(obj, setup)
            changed = false;
            if isfield(setup, 'BeamformerMode')
                m = lower(setup.BeamformerMode);
                if ~strcmp(m, obj.wfs_mode), obj.wfs_mode = m; changed = true; end
            end
            if isfield(setup, 'BeamformerNumTaps')
                if setup.BeamformerNumTaps ~= obj.N_taps
                    obj.N_taps = setup.BeamformerNumTaps; changed = true;
                end
            end
            if isfield(setup, 'CTCMode')
                m = lower(setup.CTCMode);
                if ~strcmp(m, obj.ctc_mode), obj.ctc_mode = m; changed = true; end
            end
            if isfield(setup, 'LSFINullWidth') && setup.LSFINullWidth ~= obj.lsfi_null_width
                obj.lsfi_null_width = setup.LSFINullWidth; changed = true;
            end
            if isfield(setup, 'LSFINullDepth') && setup.LSFINullDepth ~= obj.lsfi_null_depth
                obj.lsfi_null_depth = setup.LSFINullDepth; changed = true;
            end
            if isfield(setup, 'LSFIFIRLength') && setup.LSFIFIRLength ~= obj.lsfi_fir_length
                obj.lsfi_fir_length = setup.LSFIFIRLength; changed = true;
            end
            if isfield(setup, 'LSFIHeadDiameter') && setup.LSFIHeadDiameter ~= obj.lsfi_head_diameter
                obj.lsfi_head_diameter = setup.LSFIHeadDiameter; changed = true;
            end
            if isfield(setup, 'LSFIMainAngleL') && ~isequaln(setup.LSFIMainAngleL, obj.lsfi_main_angle_L)
                obj.lsfi_main_angle_L = setup.LSFIMainAngleL; changed = true;
            end
            if isfield(setup, 'LSFIMainAngleR') && ~isequaln(setup.LSFIMainAngleR, obj.lsfi_main_angle_R)
                obj.lsfi_main_angle_R = setup.LSFIMainAngleR; changed = true;
            end
            if isfield(setup, 'LSFIZeroAngleL') && ~isequaln(setup.LSFIZeroAngleL, obj.lsfi_zero_angle_L)
                obj.lsfi_zero_angle_L = setup.LSFIZeroAngleL; changed = true;
            end
            if isfield(setup, 'LSFIZeroAngleR') && ~isequaln(setup.LSFIZeroAngleR, obj.lsfi_zero_angle_R)
                obj.lsfi_zero_angle_R = setup.LSFIZeroAngleR; changed = true;
            end
            if isfield(setup, 'LSFINullAngleL') && ~isequaln(setup.LSFINullAngleL, obj.lsfi_zero_angle_L)
                obj.lsfi_zero_angle_L = setup.LSFINullAngleL; changed = true;
            end
            if isfield(setup, 'LSFINullAngleR') && ~isequaln(setup.LSFINullAngleR, obj.lsfi_zero_angle_R)
                obj.lsfi_zero_angle_R = setup.LSFINullAngleR; changed = true;
            end
            if isfield(setup, 'BeamLPos') && any(setup.BeamLPos ~= 0)
                obj.beamL_pos = setup.BeamLPos; changed = true;
            end
            if isfield(setup, 'BeamRPos') && any(setup.BeamRPos ~= 0)
                obj.beamR_pos = setup.BeamRPos; changed = true;
            end

            % Low-freq shelving
            if isfield(setup, 'LowFreqCut')
                if setup.LowFreqCut ~= obj.LowFreqCut
                    obj.LowFreqCut = setup.LowFreqCut; changed = true;
                end
            end

            % CTC plant inversion
            if isfield(setup, 'CTCInverseEnabled')
                obj.CTCInverseEnabled = logical(setup.CTCInverseEnabled);
            end
            if isfield(setup, 'CTCBeta')
                obj.ctc_beta = setup.CTCBeta;
                obj.ctc_inv_filters = [];  % force rebuild
            end
            if isfield(setup, 'EarTransferModel')
                obj.ear_transfer_model = setup.EarTransferModel;
                obj.ctc_inv_filters = [];
            end

            % WFS passthrough
            if isfield(setup, 'Tapering'),      obj.wfs_setup.Tapering = setup.Tapering; changed = true; end
            if isfield(setup, 'Antialiasing'),  obj.wfs_setup.Antialiasing = setup.Antialiasing; changed = true; end
            if isfield(setup, 'ReferenceMode'), obj.wfs_setup.ReferenceMode = setup.ReferenceMode; changed = true; end
            if isfield(setup, 'RefPoint'),      obj.wfs_setup.RefPoint = setup.RefPoint; changed = true; end
            if isfield(setup, 'RefDistance'),    obj.wfs_setup.RefDistance = setup.RefDistance; changed = true; end

            % Multiband crossover (only relevant in 'multiband' mode)
            if isfield(setup, 'MBCrossoverHz')
                fc = setup.MBCrossoverHz;
                if abs(fc - obj.mb_crossover_hz) > 1e-6
                    obj.mb_crossover_hz = fc; changed = true;
                end
            end

            if changed, obj.update_renderer(); end
        end

        % ═══════════════════════════════════════════════
        %  Render (the core)
        % ═══════════════════════════════════════════════
        function render(obj)
            inputSig = obj.virtual_source.source_signal.time_series;

            if strcmpi(obj.wfs_mode, 'lsfi')
                % LSFI mode: render via per-element FIR filters
                switch lower(obj.ctc_mode)
                    case 'beam_l'
                        obj.render_lsfi('L', inputSig(:,1));
                    case 'beam_r'
                        obj.render_lsfi('R', inputSig(:,1));
                    case 'ctc'
                        if ~obj.ctc_primary, return; end
                        obj.render_ctc();  % CTC uses LSFI internally
                end
            elseif strcmpi(obj.wfs_mode, 'race')
                % RACE mode: only the 'ctc' output mode is meaningful here.
                % Beam L / Beam R are not exposed for RACE in the UI; if they
                % are requested anyway (programmatic use), we silently fall
                % through to ctc rendering.
                if ~obj.ctc_primary, return; end
                obj.render_ctc();  % render_ctc dispatches to render_race
            elseif strcmpi(obj.wfs_mode, 'multiband')
                % Multiband mode: only 'ctc' is meaningful (analogous to RACE).
                if ~obj.ctc_primary, return; end
                obj.render_ctc();  % render_ctc dispatches to render_multiband
            else
                % WFS mode: render via internal wfs_renderer
                switch lower(obj.ctc_mode)
                    case 'beam_l'
                        obj.render_single_beam(obj.wfsL, inputSig);
                    case 'beam_r'
                        obj.render_single_beam(obj.wfsR, inputSig);
                    case 'ctc'
                        if ~obj.ctc_primary, return; end
                        obj.render_ctc();
                end
            end

            % Speaker EQ on final SSD signals
            if obj.SpeakerEQEnabled && ~isempty(obj.eqConvolvers)
                Nls = numel(obj.ssd.loudspeakers);
                Neq = numel(obj.eqConvolvers);
                for n = 1:min(Nls, Neq)
                    if isempty(obj.eqConvolvers{n}), continue; end
                    raw = obj.ssd.loudspeakers{n}.source_signal.time_series;
                    if isempty(raw), continue; end
                    filtered = obj.eqConvolvers{n}.convolve(raw);
                    obj.ssd.loudspeakers{n}.source_signal.set_signal(filtered, obj.fs);
                end
            end
        end

        % ═══════════════════════════════════════════════
        %  Speaker EQ
        % ═══════════════════════════════════════════════
        function loadSpeakerEQ(obj, matFile)
            try
                d = load(matFile, 'hinv', 'fs');
            catch ME
                warning('beamformer:eqLoad', 'Cannot load %s: %s', matFile, ME.message);
                return;
            end
            if ~isfield(d, 'hinv') || ~isfield(d, 'fs')
                warning('beamformer:eqLoad', 'File must contain hinv and fs variables.');
                return;
            end

            hinv = d.hinv;
            fs_eq = d.fs;
            if fs_eq ~= obj.fs
                fprintf('Resampling speaker EQ from %d to %d Hz\n', fs_eq, obj.fs);
                hinv = resample(hinv, obj.fs, fs_eq);
            end

            Nch = size(hinv, 2);
            Nls = numel(obj.ssd.loudspeakers);
            blk = size(obj.virtual_source.source_signal.time_series, 1);
            if blk == 0, blk = 1024; end
            obj.eqConvolvers = cell(1, Nls);
            obj.eqImpulseResponses = zeros(size(hinv, 1), Nls);
            for n = 1:min(Nch, Nls)
                obj.eqConvolvers{n} = OLS_convolver(hinv(:, n), blk);
                obj.eqImpulseResponses(:, n) = hinv(:, n);
            end
            fprintf('Speaker EQ loaded: %d filters, %d taps each\n', min(Nch,Nls), size(hinv,1));
        end

        % ═══════════════════════════════════════════════
        %  Overlay drawing
        % ═══════════════════════════════════════════════
        function output = get_renderer_visuals(obj)
            output.x0         = obj.ssd.x0;
            output.n0         = obj.ssd.n0;
            output.vs         = obj.virtual_source.position;
            output.steerAngle = obj.steerAngleL;
            output.steerAngleL = obj.steerAngleL;
            output.steerAngleR = obj.steerAngleR;
            output.mode       = obj.wfs_mode;
            output.center     = mean(obj.ssd.x0, 1);
            output.weights    = obj.weights;
            output.beamL_pos  = obj.beamL_pos;
            output.beamR_pos  = obj.beamR_pos;
            output.ctc_mode   = obj.ctc_mode;
        end

        function store = drawOverlay(obj, ax, layers, settings, store)
            % Provide safe defaults for any missing settings fields
            if ~isfield(settings, 'showRef'),              settings.showRef              = true;  end
            if ~isfield(settings, 'showAmp'),              settings.showAmp              = false; end
            if ~isfield(settings, 'showTaper'),            settings.showTaper            = false; end
            if ~isfield(settings, 'showZones'),            settings.showZones            = true;  end
            if ~isfield(settings, 'showBeamDirection'),    settings.showBeamDirection    = true;  end
            if ~isfield(settings, 'showElementHighlight'), settings.showElementHighlight = true;  end

            wfs_settings = struct( ...
                'showRef',   settings.showRef, ...
                'showAmp',   settings.showAmp, ...
                'showTaper', settings.showTaper, ...
                'showZones', settings.showZones);

            % Each beam gets its own sub-store to avoid key collisions in WFSPainter
            if ~isfield(store, 'wfsL_store'), store.wfsL_store = struct(); end
            if ~isfield(store, 'wfsR_store'), store.wfsR_store = struct(); end

            if ~isempty(obj.wfsL)
                store.wfsL_store = audioapp.ui.draw.WFSPainter.draw( ...
                    ax, obj.wfsL.get_renderer_visuals(), obj.wfsL, wfs_settings, store.wfsL_store);
            end
            if ~isempty(obj.wfsR)
                store.wfsR_store = audioapp.ui.draw.WFSPainter.draw( ...
                    ax, obj.wfsR.get_renderer_visuals(), obj.wfsR, wfs_settings, store.wfsR_store);
            end

            % BeamformerPainter (arrows/crosses for both beams)
            rp = obj.get_renderer_visuals();
            bf_settings = struct( ...
                'showBeamDirection',    settings.showBeamDirection, ...
                'showElementHighlight', settings.showElementHighlight);
            store = audioapp.ui.draw.BeamformerPainter.draw(ax, rp, bf_settings, store);
        end

        function settings = defaultOverlaySettings(~)
            settings = struct( ...
                'showBeamDirection',    true, ...
                'showElementHighlight', true, ...
                'showRef',    true, ...
                'showAmp',    false, ...
                'showTaper',  false, ...
                'showZones',  true);
        end

        % ═══════════════════════════════════════════════
        %  Ear transfer analysis (uses L beam)
        % ═══════════════════════════════════════════════
        function [f, H_LL, H_LR, H_RL, H_RR] = compute_ear_transfer(obj, model, hrtf_extrap)
            % Compute 2x2 ear transfer matrix for both beams.
            %   H_LL : L beam -> left ear       H_LR : L beam -> right ear
            %   H_RL : R beam -> left ear       H_RR : R beam -> right ear
            % Driving spectra read directly from wfsL/wfsR.driving_spectrum
            % (the true WFS transfer stored by wfs_renderer.update_renderer).
            % All transfers normalised to free-field propagation.
            if nargin < 2, model = 'analytic'; end
            if nargin < 3, hrtf_extrap = []; end

            Nls     = numel(obj.ssd.loudspeakers);
            rec_pos = obj.receiver.position;
            rec_ori = obj.receiver.orientation;

            % Driving spectra [Nf x Nls]
            if strcmpi(obj.wfs_mode, 'lsfi')
                % LSFI: compute spectra from stored FIR matrices
                if isempty(obj.lsfi_fir_L) || isempty(obj.lsfi_fir_R)
                    error('beamformer_renderer:noLSFI', ...
                        'LSFI filters not built yet. Run update_renderer first.');
                end
                NfftL = size(obj.lsfi_fir_L, 1);
                Nf    = NfftL/2 + 1;
                DfL   = fft(obj.lsfi_fir_L, NfftL, 1);
                DfL   = DfL(1:Nf, :);
                DfR   = fft(obj.lsfi_fir_R, NfftL, 1);
                DfR   = DfR(1:Nf, :);
            elseif strcmpi(obj.wfs_mode, 'race')
                % RACE: only the two outermost loudspeakers are active.
                % L "beam" = the cross-coupled pair driven by a unit input
                % on the L channel: out_L = H_direct, out_R = H_cross.
                % R "beam" = unit input on R channel: out_L = H_cross,
                % out_R = H_direct (by symmetry of the 2x2 RACE matrix).
                if isempty(obj.race_H_direct) || isempty(obj.race_H_cross)
                    error('beamformer_renderer:noRACE', ...
                        'RACE filters not built yet. Run update_renderer first.');
                end
                Nf  = numel(obj.race_H_direct);
                DfL = zeros(Nf, Nls);
                DfR = zeros(Nf, Nls);
                DfL(:, obj.race_idx_L) = obj.race_H_direct;   % in_L -> out_L
                DfL(:, obj.race_idx_R) = obj.race_H_cross;    % in_L -> out_R
                DfR(:, obj.race_idx_L) = obj.race_H_cross;    % in_R -> out_L
                DfR(:, obj.race_idx_R) = obj.race_H_direct;   % in_R -> out_R
            elseif strcmpi(obj.wfs_mode, 'multiband')
                % Multiband: combine RACE (low band) and LSFI (high band)
                % driving spectra via the same FIR crossover used at runtime.
                if isempty(obj.race_H_direct) || isempty(obj.lsfi_fir_L) ...
                        || isempty(obj.mb_h_lp)
                    error('beamformer_renderer:noMultiband', ...
                        'Multiband filters not built yet. Run update_renderer first.');
                end
                % LSFI side ─ already on a NfftL grid
                NfftL = size(obj.lsfi_fir_L, 1);
                Nf    = NfftL/2 + 1;
                DfL_lsfi = fft(obj.lsfi_fir_L, NfftL, 1);
                DfL_lsfi = DfL_lsfi(1:Nf, :);
                DfR_lsfi = fft(obj.lsfi_fir_R, NfftL, 1);
                DfR_lsfi = DfR_lsfi(1:Nf, :);
                % RACE side ─ may be on a different FFT length, resample to NfftL
                race_dir = obj.race_H_direct;
                race_crs = obj.race_H_cross;
                if numel(race_dir) ~= Nf
                    % Re-evaluate the closed-form RACE spectra on the LSFI grid
                    fL = (0:Nf-1)' * obj.fs / NfftL;
                    omL = 2*pi*fL;
                    cL = obj.race_delta_a * exp(-1i * omL * obj.race_delta_t);
                    race_dir = 1 ./ (1 - cL.^2);
                    race_crs = -cL ./ (1 - cL.^2);
                end
                DfL_race = zeros(Nf, Nls);
                DfR_race = zeros(Nf, Nls);
                DfL_race(:, obj.race_idx_L) = race_dir;
                DfL_race(:, obj.race_idx_R) = race_crs;
                DfR_race(:, obj.race_idx_L) = race_crs;
                DfR_race(:, obj.race_idx_R) = race_dir;
                % Crossover spectra on the same NfftL grid
                H_lp = fft(obj.mb_h_lp, NfftL, 1);  H_lp = H_lp(1:Nf);
                H_hp = fft(obj.mb_h_hp, NfftL, 1);  H_hp = H_hp(1:Nf);
                % Frequency-selective combination
                DfL = (H_lp .* DfL_race) + (H_hp .* DfL_lsfi);
                DfR = (H_lp .* DfR_race) + (H_hp .* DfR_lsfi);
            else
                % WFS: read driving spectra from wfs_renderer
                DfL = obj.wfsL.driving_spectrum;
                DfR = obj.wfsR.driving_spectrum;
                Nf  = size(DfL, 1);
            end
            f = (0:Nf-1)' * obj.fs / (2*(Nf-1));

            % Optional per-loudspeaker speaker EQ, applied in render() after
            % beam synthesis. Fold the same spectra into the ear-transfer
            % model so the analysis matches the radiated field path.
            if obj.SpeakerEQEnabled && ~isempty(obj.eqImpulseResponses)
                Neq = min(size(obj.eqImpulseResponses, 2), Nls);
                if Neq > 0
                    Nfft = 2*(Nf-1);
                    H_eq = fft(obj.eqImpulseResponses(:, 1:Neq), Nfft, 1);
                    H_eq = H_eq(1:Nf, :);
                    DfL(:, 1:Neq) = DfL(:, 1:Neq) .* H_eq;
                    DfR(:, 1:Neq) = DfR(:, 1:Neq) .* H_eq;
                end
            end

            % HRTF / Diffractor matrix [Nf x Nls]
            G_earL = zeros(Nf, Nls);
            G_earR = zeros(Nf, Nls);
            switch lower(model)
                case 'analytic'
                    d = Diffractor();
                    d.scatterer_pos    = rec_pos;
                    d.scatterer_radius = 0.0875;
                    for n = 1:Nls
                        d.source_pos  = obj.ssd.x0(n, :);
                        d.source_type = 'point';
                        H_ears = d.evaluate_ears(f', rec_ori);
                        G_earL(:, n) = H_ears(1, :).';
                        G_earR(:, n) = H_ears(2, :).';
                    end
                case 'hrtf'
                    if isempty(hrtf_extrap)
                        error('beamformer_renderer:noHRTF', 'hrtf_extrapolator required');
                    end
                    Nfft = 2*(Nf-1);
                    for n = 1:Nls
                        hrir = hrtf_extrap.extrapolate_hrtf(obj.ssd.loudspeakers{n}, obj.receiver);
                        H    = fft(hrir, Nfft, 1);
                        G_earL(:, n) = H(1:Nf, 1);
                        G_earR(:, n) = H(1:Nf, 2);
                    end
                otherwise
                    error('beamformer_renderer:badModel', 'Unknown model: %s', model);
            end

            % Normalise by free-field amplitude at receiver:  1 / (4π · |beam_pos - rec_pos|)
            % Same reference as focus_normalisation_gain used in render_single_beam.
            % Scalar per beam — no frequency dependence, no per-loudspeaker weighting.
            if strcmpi(obj.wfs_mode, 'race')
                % RACE has no beam focus point — the L/R "channels" are the
                % two outermost loudspeakers themselves. Normalise by the
                % free-field propagation from each channel speaker to the
                % nominal head centre (rec_pos), giving a comparable level
                % reference to LSFI/WFS modes.
                xL = obj.ssd.x0(obj.race_idx_L, :);
                xR = obj.ssd.x0(obj.race_idx_R, :);
                norm_L = 1 / (4*pi * max(norm(xL - rec_pos(:).'), 0.01));
                norm_R = 1 / (4*pi * max(norm(xR - rec_pos(:).'), 0.01));
            else
                % LSFI, WFS-based, and Multiband all use the LSFI/beam focal
                % points as reference (in Multiband the high band dominates).
                norm_L = 1 / (4*pi * max(norm(obj.beamL_pos(:).' - rec_pos(:).'), 0.01));
                norm_R = 1 / (4*pi * max(norm(obj.beamR_pos(:).' - rec_pos(:).'), 0.01));
            end

            % Sum over loudspeakers
            H_LL_raw = sum(DfL .* G_earL, 2);  H_LR_raw = sum(DfL .* G_earR, 2);
            H_RL_raw = sum(DfR .* G_earL, 2);  H_RR_raw = sum(DfR .* G_earR, 2);

            % Normalise
            H_LL = H_LL_raw / norm_L;   H_LR = H_LR_raw / norm_L;
            H_RL = H_RL_raw / norm_R;   H_RR = H_RR_raw / norm_R;
        end
    end

    % =======================================================
    %  Private: beam building, render helpers, normalisation
    % =======================================================
    methods (Access = private)
        function build_beam_wfs(obj, side)
            % Build or update one internal WFS renderer ('L' or 'R')
            if side == 'L'
                beam_pos = obj.beamL_pos;
            else
                beam_pos = obj.beamR_pos;
            end

            x0 = obj.ssd.x0;
            D  = size(x0, 2);
            if numel(beam_pos) < D
                beam_pos = [beam_pos(:)', zeros(1, D-numel(beam_pos))];
            else
                beam_pos = beam_pos(1:D);
            end

            center = mean(x0, 1);

            switch lower(obj.wfs_mode)
                case 'das'
                    dir = beam_pos - center;
                    if norm(dir) > 1e-12, dir = dir / norm(dir);
                    else, dir = zeros(1,D); dir(2) = 1; end
                    src_R = 0.025;
                    wfs_vs_type = struct('Shape', 'plane_wave', 'R', src_R);
                    wfs_vs_pos  = -dir;
                otherwise  % 'focus' or 'directive' (directive builds on focus)
                    src_R = 0.025;
                    wfs_vs_type = struct('Shape', 'point_source', 'R', src_R);
                    wfs_vs_pos  = beam_pos;
            end

            wfs_rs = obj.wfs_setup;
            wfs_rs.ReferenceMode = 'fixed_point';
            if ~isempty(obj.receiver)
                rp = obj.receiver.position(:)';
                if numel(rp) < D, rp = [rp, zeros(1, D-numel(rp))]; end
                wfs_rs.RefPoint = rp(1:D);
            end

            if side == 'L'
                if isempty(obj.wfsL)
                    wfs_vs = virtual_source(1, wfs_vs_pos, [0,1,0], wfs_vs_type, 'wfs');
                    wfs_vs.set_input(obj.virtual_source.source_signal.time_series, obj.fs);
                    obj.wfsL = wfs_renderer(wfs_vs, obj.ssd, wfs_rs);
                else
                    obj.wfsL.virtual_source.position = wfs_vs_pos;
                    obj.wfsL.virtual_source.source_type = wfs_vs_type;
                    obj.wfsL.setup = wfs_rs;
                end
                obj.wfsL.LowCutEnabled = obj.LowFreqCut;
                if ~isempty(obj.ear_receiver_L)
                    obj.wfsL.receiver_ref = obj.ear_receiver_L;
                end
                obj.wfsL.update_renderer();
            else
                if isempty(obj.wfsR)
                    wfs_vs = virtual_source(1, wfs_vs_pos, [0,1,0], wfs_vs_type, 'wfs');
                    wfs_vs.set_input(obj.virtual_source.source_signal.time_series, obj.fs);
                    obj.wfsR = wfs_renderer(wfs_vs, obj.ssd, wfs_rs);
                else
                    obj.wfsR.virtual_source.position = wfs_vs_pos;
                    obj.wfsR.virtual_source.source_type = wfs_vs_type;
                    obj.wfsR.setup = wfs_rs;
                end
                obj.wfsR.LowCutEnabled = obj.LowFreqCut;
                if ~isempty(obj.ear_receiver_R)
                    obj.wfsR.receiver_ref = obj.ear_receiver_R;
                end
                obj.wfsR.update_renderer();
            end
        end

        function apply_directive_correction(obj)
            % Modify the WFS driving functions of both beams so that each
            % beam produces a "directive point source": unity gain at the
            % ipsilateral ear and a null at the contralateral ear.
            %
            % The WFS focus mode already synthesises a point source at the
            % beam target position. Here we apply a per-element, per-freq
            % complex correction α_n(ω) to the WFS filter bank so that the
            % resulting sound field satisfies:
            %
            %   Σ_n  α_n · D_n(ω) · G(x_target, x_n; ω) = 1    (target ear)
            %   Σ_n  α_n · D_n(ω) · G(x_null,   x_n; ω) = 0    (null ear)
            %
            % Using a matched-filter ansatz α_n = a·g_t*(n) + b·g_u*(n)
            % where g_t(n) = D_n·G_target(n) and g_u(n) = D_n·G_null(n),
            % reduces to a 2×2 system per frequency bin.

            if isempty(obj.wfsL) || isempty(obj.wfsR), return; end
            if isempty(obj.ear_receiver_L) || isempty(obj.ear_receiver_R), return; end

            c = audioapp.util.Physics.speedOfSound();
            x0 = obj.ssd.x0;
            N  = size(x0, 1);
            D  = size(x0, 2);

            earL = obj.ear_receiver_L.position(:).';
            earR = obj.ear_receiver_R.position(:).';
            if numel(earL) < D, earL = [earL, zeros(1,D-numel(earL))]; end
            if numel(earR) < D, earR = [earR, zeros(1,D-numel(earR))]; end

            for side_idx = 1:2
                if side_idx == 1
                    wfs = obj.wfsL;
                    target_pos = earL;
                    null_pos   = earR;
                else
                    wfs = obj.wfsR;
                    target_pos = earR;
                    null_pos   = earL;
                end

                Df = wfs.driving_spectrum;   % [Nf × N] one-sided
                if isempty(Df), continue; end
                Nf   = size(Df, 1);
                Nfft = 2*(Nf-1);
                freq = (0:Nf-1)' * obj.fs / Nfft;

                % Green's functions from each SSD element to target and null ears
                dist_t = sqrt(sum((x0(:,1:D) - target_pos(1:D)).^2, 2));  % [N×1]
                dist_n = sqrt(sum((x0(:,1:D) - null_pos(1:D)).^2, 2));    % [N×1]

                % Corrected driving spectrum
                Df_corrected = Df;

                for kf = 2:Nf   % skip DC
                    k = 2*pi*freq(kf) / c;

                    % Green's function vectors [N×1]
                    G_t = exp(-1i * k * dist_t) ./ (4*pi * dist_t);
                    G_n = exp(-1i * k * dist_n) ./ (4*pi * dist_n);

                    % Current driving function row [1×N]
                    d = Df(kf, :).';   % [N×1]

                    % Combined transfer: g_t = D.*G_target, g_n = D.*G_null
                    g_t = d .* G_t;    % [N×1]
                    g_n = d .* G_n;    % [N×1]

                    % Gram matrix for the matched-filter ansatz
                    A11 = g_t' * g_t;
                    A12 = g_n' * g_t;
                    A21 = g_t' * g_n;
                    A22 = g_n' * g_n;

                    A = [A11, A12; A21, A22];
                    scale = max(abs([A11, A22]));
                    if ~isfinite(scale) || scale < 1e-30
                        continue;
                    end

                    % Frequency-dependent regularisation:
                    % At low frequencies (kd_head << 1) the two ears are
                    % acoustically similar → Gram matrix near-singular →
                    % correction blows up. We regularise proportionally to
                    % how similar the two steering vectors are.
                    % Similarity measure: |g_t'g_n|² / (||g_t||²·||g_n||²)
                    similarity = abs(A12)^2 / (abs(A11)*abs(A22) + 1e-30);
                    if ~isfinite(similarity)
                        similarity = 1;
                    end
                    % When similarity → 1 (low freq), reg → large → no correction
                    % When similarity → 0 (high freq), reg → small → strong null
                    reg = 1e-3 * scale * (0.01 + similarity);
                    A = A + reg * eye(2);

                    % Preserve the original focused-WFS target-ear transfer
                    % envelope at this frequency. This keeps common spectral
                    % shaping such as the WFS prefilter / low-frequency shelf,
                    % while adding the contralateral-ear null.
                    desired_target = sum(g_t);
                    rhs = [desired_target; 0];
                    if ~all(isfinite(rhs))
                        continue;
                    end
                    if rcond(A) < 1e-12
                        ab = pinv(A) * rhs;
                    else
                        ab = A \ rhs;
                    end
                    if ~all(isfinite(ab))
                        continue;
                    end
                    a = ab(1);
                    b = ab(2);

                    % Per-element correction
                    alpha = a * conj(g_t) + b * conj(g_n);
                    if ~all(isfinite(alpha))
                        continue;
                    end

                    Df_corrected(kf, :) = (d .* alpha).';
                end

                % Convert back to time-domain FIR and update the WFS filter bank
                for n = 1:N
                    h = fftshift(ifft(Df_corrected(:,n), Nfft, 1, 'symmetric'));
                    wfs.WFS_filter_bank{n}.update_coefficients(h);
                end

                % Also update driving_spectrum for ear transfer analysis
                wfs.driving_spectrum = Df_corrected;
            end
        end

        function build_beam_lsfi(obj, side)
            obj.build_beam_lsfi_paper(side);
        end

%{
            % Design LSFI beamformer filters for one beam using near-field
            % LCMV (Linearly Constrained Minimum Variance) beamforming.
            %
            % Based on Mabande et al. (ICASSP 2009) and Hohnerlein & Ahrens
            % (ICASSP 2017), adapted for near-field:
            %
            %   Target control point: ipsilateral ear → response = 1
            %   Null control point:   contralateral ear → response = 0
            %   (+optional null width: extra points around null ear)
            %
            %   min ||w||²  subject to  C^H·w = f
            %
            % where C = [d_target, d_null1, d_null2, ...] are near-field
            % steering vectors (Green's functions), and f = [1; 0; 0; ...].
            %
            % Tikhonov regularisation (WNG constraint proxy):
            %   w = C·(C^H·C + β·I)^{-1} · f
            %
            % The result is IFFT-ed to per-element FIR filters.

            x0  = obj.ssd.x0;
            N   = size(x0, 1);
            D   = size(x0, 2);
            c   = audioapp.util.Physics.speedOfSound();

            % Target and null ear positions
            if side == 'L'
                target_pos = obj.ear_receiver_L.position(:).';
                null_pos   = obj.ear_receiver_R.position(:).';
            else
                target_pos = obj.ear_receiver_R.position(:).';
                null_pos   = obj.ear_receiver_L.position(:).';
            end
            if numel(target_pos) < D, target_pos = [target_pos, zeros(1,D-numel(target_pos))]; end
            if numel(null_pos) < D,   null_pos   = [null_pos,   zeros(1,D-numel(null_pos))];   end

            % Extra null points around the contralateral ear (null width)
            % Distribute points on a small arc around the null ear
            hw = obj.lsfi_null_width / 2;  % degrees
            n_null_pts = max(round(hw), 1);
            null_positions = null_pos;  % first row = exact null ear
            if n_null_pts > 1
                % Generate extra points in a small disc around null_pos
                center2d = mean(x0(:,1:2), 1);
                R_null = norm(null_pos(1:2) - center2d);
                base_angle = atan2(null_pos(2)-center2d(2), null_pos(1)-center2d(1));
                offsets = linspace(-hw, hw, n_null_pts) * pi/180;
                offsets(offsets == 0) = [];  % already have the center
                for ko = 1:numel(offsets)
                    a = base_angle + offsets(ko);
                    extra = [center2d(1) + R_null*cos(a), ...
                             center2d(2) + R_null*sin(a)];
                    if D >= 3
                        extra = [extra, null_pos(3)]; %#ok<AGROW>
                    end
                    null_positions = [null_positions; extra]; %#ok<AGROW>
                end
            end
            Nnull = size(null_positions, 1);

            % Frequency grid
            Nfft = obj.lsfi_fir_length;
            Nf   = Nfft/2 + 1;
            freq = (0:Nf-1)' * obj.fs / Nfft;

            % Regularisation: β controls WNG / robustness trade-off.
            % Map lsfi_null_depth (dB) to a β value.
            % Smaller β → more aggressive null (less robust).
            % Typical: β = 10^(-null_depth_dB/20) * 1/N
            beta = 10^(abs(obj.lsfi_null_depth)/20) / N^2;
            % Clamp to a reasonable range
            beta = max(beta, 1e-6);
            beta = min(beta, 1);

            W_all = zeros(Nf, N);

            for kf = 1:Nf
                omega = 2*pi*freq(kf);
                if omega < 1e-6
                    W_all(kf,:) = 1/N;
                    continue;
                end

                % Near-field steering vectors (Green's function):
                %   d_n(x) = exp(-j·k·|x - x_n|) / (4π·|x - x_n|)
                % where k = ω/c.
                k = omega / c;

                % Target steering vector [N × 1]
                dist_target = sqrt(sum((x0(:,1:D) - target_pos(1:D)).^2, 2));
                d_target = exp(-1i * k * dist_target) ./ (4*pi * dist_target);

                % Null steering vectors [N × Nnull]
                D_null = zeros(N, Nnull);
                for kn = 1:Nnull
                    dist_n = sqrt(sum((x0(:,1:D) - null_positions(kn,1:D)).^2, 2));
                    D_null(:, kn) = exp(-1i * k * dist_n) ./ (4*pi * dist_n);
                end

                % Constraint matrix: C = [d_target, D_null]  size [N × (1+Nnull)]
                C_mat = [d_target, D_null];
                % Desired response: f = [1; 0; 0; ...]
                f_vec = [1; zeros(Nnull, 1)];

                % LCMV solution with Tikhonov regularisation:
                %   w = C·(C^H·C + β·I)^{-1}·f
                Nc = size(C_mat, 2);
                CTC = C_mat' * C_mat;
                w = C_mat * ((CTC + beta * eye(Nc)) \ f_vec);

                W_all(kf, :) = w.';
            end

            % IFFT to time-domain FIR filters
            blk = size(obj.virtual_source.source_signal.time_series, 1);
            if blk == 0, blk = 1024; end
            filters = cell(1, N);
            fir_mat = zeros(Nfft, N);
            for n = 1:N
                h = fftshift(ifft(W_all(:,n), Nfft, 1, 'symmetric'));
                h = h .* tukeywin(Nfft, 0.1);
                fir_mat(:, n) = h;
                filters{n} = OLS_convolver(h, blk);
            end

            if side == 'L'
                obj.lsfi_filters_L = filters;
                obj.lsfi_fir_L = fir_mat;
            else
                obj.lsfi_filters_R = filters;
                obj.lsfi_fir_R = fir_mat;
            end
        end

%}

        function build_beam_lsfi_paper(obj, side)
            % Angular-domain LSFI after Hohnerlein/Ahrens: fit a target
            % lobe over listening angles while constraining stop-band energy.
            obj.update_ear_receivers();

            x0 = obj.ssd.x0;
            N = size(x0, 1);
            if N == 0
                error('beamformer_renderer:emptySSD', ...
                    'LSFI requires at least one loudspeaker.');
            end

            [array_axis, broadside, array_center] = obj.lsfi_array_axes();

            [target_angle, stop_angle] = obj.lsfi_effective_angles( ...
                side, array_center, array_axis, broadside);

            step = max(double(obj.lsfi_angle_step), 0.25);
            span_needed = max(abs([target_angle, stop_angle])) + double(obj.lsfi_null_width) + 20;
            span = min(89, max(double(obj.lsfi_angle_span), span_needed));
            angles = (-span:step:span).';
            if isempty(angles)
                angles = (-80:step:80).';
            end

            desired = obj.lsfi_desired_response(angles, target_angle, stop_angle);
            stop_mask = abs(angles - stop_angle) <= double(obj.lsfi_null_width) / 2;
            if ~any(stop_mask)
                [~, idx_stop] = min(abs(angles - stop_angle));
                stop_mask(idx_stop) = true;
            end

            Nfft = max(64, round(double(obj.lsfi_fir_length)));
            if mod(Nfft, 2) ~= 0
                Nfft = Nfft + 1;
            end
            Nf = Nfft / 2 + 1;
            freq = (0:Nf-1).' * obj.fs / Nfft;
            eps_stop = max(10^(-abs(double(obj.lsfi_null_depth)) / 20), 1e-6);

            W_all = zeros(Nf, N);
            for kf = 1:Nf
                f = freq(kf);
                if f <= 0
                    continue;
                end

                G = obj.lsfi_transfer_matrix(angles, f, array_axis, broadside, array_center);
                W_all(kf, :) = obj.lsfi_constrained_ls( ...
                    G, desired, G(stop_mask, :), eps_stop).';
            end

            W_all(1, :) = real(W_all(1, :));
            W_all(end, :) = real(W_all(end, :));
            fir_mat = 0.01*obj.lsfi_spectrum_to_fir(W_all, Nfft);

            blk = size(obj.virtual_source.source_signal.time_series, 1);
            if blk == 0
                blk = 1024;
            end

            filters = cell(1, N);
            for n = 1:N
                filters{n} = OLS_convolver(fir_mat(:, n), blk);
            end

            if strcmpi(side, 'L')
                obj.lsfi_filters_L = filters;
                obj.lsfi_fir_L = fir_mat;
            else
                obj.lsfi_filters_R = filters;
                obj.lsfi_fir_R = fir_mat;
            end
        end

        function [target_angle, stop_angle] = lsfi_effective_angles(obj, side, array_center, array_axis, broadside)
            has_ears = ~isempty(obj.ear_receiver_L) && ~isempty(obj.ear_receiver_R);
            if strcmpi(side, 'L')
                manual_target = obj.lsfi_main_angle_L;
                manual_stop = obj.lsfi_zero_angle_L;
                fallback_target = -6;
                fallback_stop = 6;
                if has_ears
                    target_pos = obj.ear_receiver_L.position(:).';
                    stop_pos = obj.ear_receiver_R.position(:).';
                end
            else
                manual_target = obj.lsfi_main_angle_R;
                manual_stop = obj.lsfi_zero_angle_R;
                fallback_target = 6;
                fallback_stop = -6;
                if has_ears
                    target_pos = obj.ear_receiver_R.position(:).';
                    stop_pos = obj.ear_receiver_L.position(:).';
                end
            end

            if has_ears
                target_angle = obj.lsfi_angle_from_point(target_pos, array_center, array_axis, broadside);
                stop_angle = obj.lsfi_angle_from_point(stop_pos, array_center, array_axis, broadside);
            else
                target_angle = fallback_target;
                stop_angle = fallback_stop;
            end

            if isfinite(manual_target)
                target_angle = manual_target;
            end
            if isfinite(manual_stop)
                stop_angle = manual_stop;
            end
        end

        function label = lsfi_angle_label(~, angle_deg)
            if isfinite(angle_deg)
                label = sprintf('%.1f', angle_deg);
            else
                label = 'auto';
            end
        end

        function [array_axis, broadside, array_center] = lsfi_array_axes(obj)
            x2 = obj.ssd.x0(:, 1:2);
            array_center = mean(x2, 1);
            X = x2 - array_center;

            if size(x2, 1) >= 2 && any(abs(X(:)) > 0)
                [~, ~, V] = svd(X, 'econ');
                array_axis = V(:, 1).';
            else
                array_axis = [1, 0];
            end

            [~, imax] = max(abs(array_axis));
            if array_axis(imax) < 0
                array_axis = -array_axis;
            end
            array_axis = array_axis ./ max(norm(array_axis), eps);

            if ~isempty(obj.receiver) && ~isempty(obj.receiver.position)
                ref = obj.receiver.position(1:2) - array_center;
            else
                ref = -mean(obj.ssd.n0(:, 1:2), 1);
            end
            broadside = ref - dot(ref, array_axis) * array_axis;
            if norm(broadside) < 1e-9
                broadside = [-array_axis(2), array_axis(1)];
            end
            broadside = broadside ./ max(norm(broadside), eps);
        end

        function angle_deg = lsfi_angle_from_point(~, point, array_center, array_axis, broadside)
            v = point(1:2) - array_center;
            if norm(v) < 1e-12
                angle_deg = 0;
                return;
            end
            v = v ./ norm(v);
            angle_deg = atan2d(dot(v, array_axis), dot(v, broadside));
        end

        function desired = lsfi_desired_response(~, angles, target_angle, zero_angle)
            first_zero_width = max(abs(zero_angle - target_angle), 1);
            offsets = -4:4;
            values = [0.2360, 0.4719, 0.7079, 0.9900, 1.0000, ...
                      0.9900, 0.7079, 0.4719, 0.2360];
            shape_angles = target_angle + first_zero_width * [-1, offsets / 5, 1];
            shape_values = [0, values, 0];
            desired = interp1(shape_angles, shape_values, angles, 'linear', 0);
            desired = desired(:);
        end

        function G = lsfi_transfer_matrix(obj, angles, freq, array_axis, broadside, array_center)
            c = audioapp.util.Physics.speedOfSound();
            k = 2 * pi * freq / c;
            rel = obj.ssd.x0(:, 1:2) - array_center;
            dirs = cosd(angles) .* broadside + sind(angles) .* array_axis;
            phase = dirs * rel.';
            G = exp(1i * k * phase);
        end

        function w = lsfi_constrained_ls(~, A, b, B, eps_stop)
            N = size(A, 2);
            AHA = A' * A;
            AHb = A' * b;
            scale = real(trace(AHA)) / max(N, 1);
            ridge = max(scale * 1e-8, 1e-12);

            w = (AHA + ridge * eye(N)) \ AHb;
            if isempty(B) || norm(B * w) <= eps_stop
                return;
            end

            BHB = B' * B;
            mu_hi = 1;
            for iter = 1:40
                w_hi = (AHA + mu_hi * BHB + ridge * eye(N)) \ AHb;
                if norm(B * w_hi) <= eps_stop
                    break;
                end
                mu_hi = mu_hi * 10;
            end

            mu_lo = 0;
            for iter = 1:50
                mu = 0.5 * (mu_lo + mu_hi);
                w_mu = (AHA + mu * BHB + ridge * eye(N)) \ AHb;
                if norm(B * w_mu) > eps_stop
                    mu_lo = mu;
                else
                    mu_hi = mu;
                    w_hi = w_mu;
                end
            end
            w = w_hi;
        end

        function fir_mat = lsfi_spectrum_to_fir(~, W_pos, Nfft)
            W_pos(1, :) = real(W_pos(1, :));
            W_pos(end, :) = real(W_pos(end, :));
            W_full = [W_pos; conj(W_pos(end-1:-1:2, :))];
            fir_mat = real(ifft(W_full, Nfft, 1));
            fir_mat = fftshift(fir_mat, 1);
            fir_mat = fir_mat .* tukeywin(Nfft, 0.1);
        end

        function render_lsfi(obj, side, mono)
            % Render one LSFI beam: convolve mono input with per-speaker
            % FIR filters and add to the SSD loudspeaker signals.
            if side == 'L'
                filters = obj.lsfi_filters_L;
            else
                filters = obj.lsfi_filters_R;
            end
            if isempty(filters), return; end

            Nls = numel(obj.ssd.loudspeakers);
            for n = 1:min(Nls, numel(filters))
                if isempty(filters{n}), continue; end
                sig = filters{n}.convolve(mono);
                % Add to existing SSD signal (other beam may have written already)
                existing = obj.ssd.loudspeakers{n}.source_signal.time_series;
                if ~isempty(existing) && length(existing) == length(sig)
                    obj.ssd.loudspeakers{n}.source_signal.set_signal(existing + sig, obj.fs);
                else
                    obj.ssd.loudspeakers{n}.source_signal.set_signal(sig, obj.fs);
                end
            end
        end

        function build_beam_race(obj)
            % Build the 2x2 RACE filter matrix for the two outermost
            % loudspeakers, with auto-computed delta_a and delta_t.
            %
            % RACE (Recursive Ambiophonic Crosstalk Elimination, Glasgal 2007;
            % low-band component of Hohnerlein & Ahrens ICASSP 2017):
            %
            %   The contralateral leakage from each of two stereo channels
            %   is cancelled by a delayed, attenuated, sign-flipped copy
            %   emitted from the other channel. Because that cancellation
            %   signal also leaks contralaterally, the recursion is captured
            %   by a closed-form 2x2 transfer matrix:
            %
            %       G(ω) = -delta_a * exp(-jω·delta_t)
            %
            %       [out_L]   1   [   1   -G ] [in_L]
            %       [out_R] = ─── [  -G    1 ] [in_R]
            %                1-G²
            %
            %   In the rendering chain, [in_L; in_R] = [bL; bR], the desired
            %   ear signals (already mixed from all virtual sources via HRTF
            %   in render_ctc). The outputs go to the two outermost
            %   loudspeakers; the 6 middle elements stay silent. By design.

            % ── Identify the two outermost loudspeakers ──────────────────
            % The "outermost" axis is taken to be the SSD's principal
            % horizontal axis. For a linear array this is x; we look at the
            % spread along x and y and pick the larger one.
            x0 = obj.ssd.x0;
            spread_x = max(x0(:,1)) - min(x0(:,1));
            spread_y = max(x0(:,2)) - min(x0(:,2));
            if spread_y > spread_x
                axis_col = 2;
            else
                axis_col = 1;
            end
            [~, idxL_arr] = min(x0(:, axis_col));
            [~, idxR_arr] = max(x0(:, axis_col));
            obj.race_idx_L = idxL_arr;
            obj.race_idx_R = idxR_arr;

            % ── Auto-compute delta_t (time-of-arrival difference) ────────
            c = audioapp.util.Physics.speedOfSound();
            xL_pos = x0(obj.race_idx_L, :);
            xR_pos = x0(obj.race_idx_R, :);

            if isempty(obj.ear_receiver_L) || isempty(obj.ear_receiver_R)
                obj.update_ear_receivers();
            end
            earL_pos = obj.ear_receiver_L.position(:).';
            earR_pos = obj.ear_receiver_R.position(:).';

            % Pad to 3D
            if numel(xL_pos)   < 3, xL_pos   = [xL_pos,   zeros(1,3-numel(xL_pos))]; end
            if numel(xR_pos)   < 3, xR_pos   = [xR_pos,   zeros(1,3-numel(xR_pos))]; end
            if numel(earL_pos) < 3, earL_pos = [earL_pos, zeros(1,3-numel(earL_pos))]; end
            if numel(earR_pos) < 3, earR_pos = [earR_pos, zeros(1,3-numel(earR_pos))]; end

            % Ipsilateral and contralateral path lengths for the L channel.
            % By assumption, RACE's symmetric heuristic uses one delta_t
            % and one delta_a for both channels — we average over the two
            % sides for robustness.
            d_LL = norm(xL_pos - earL_pos);    % L speaker to L ear (ipsi)
            d_LR = norm(xL_pos - earR_pos);    % L speaker to R ear (contra)
            d_RR = norm(xR_pos - earR_pos);    % R speaker to R ear (ipsi)
            d_RL = norm(xR_pos - earL_pos);    % R speaker to L ear (contra)

            dt_L = (d_LR - d_LL) / c;          % path-length excess for L channel
            dt_R = (d_RL - d_RR) / c;          % path-length excess for R channel
            obj.race_delta_t = max(0.5 * (dt_L + dt_R), 1e-5);

            % ── Auto-compute delta_a (ipsi vs contra amplitude ratio) ────
            % Use rigid-sphere head model at race_calibration_freq (1 kHz
            % default). This is consistent with the binauralizer's physical
            % model used elsewhere (CTC analytic plant, sound-field sim
            % head scattering).
            try
                d = Diffractor();
                d.scatterer_pos    = obj.receiver.position;
                d.scatterer_radius = 0.0875;
                d.source_type      = 'point';
                rec_ori = obj.receiver.orientation;

                % Drive from the L speaker, read both ears
                d.source_pos = xL_pos;
                H_from_L = d.evaluate_ears(obj.race_calibration_freq, rec_ori);
                a_L = abs(H_from_L(2)) / max(abs(H_from_L(1)), 1e-15);

                % Drive from the R speaker, read both ears
                d.source_pos = xR_pos;
                H_from_R = d.evaluate_ears(obj.race_calibration_freq, rec_ori);
                a_R = abs(H_from_R(1)) / max(abs(H_from_R(2)), 1e-15);

                % Average and clip to (0, 0.95) to keep the recursion stable
                obj.race_delta_a = min(0.95, max(0.05, 0.5 * (a_L + a_R)));
            catch
                % Fallback: literature default
                obj.race_delta_a = 0.5;
            end

            % ── Build the 2x2 transfer matrix in the frequency domain ────
            % Acoustic forward model (input ear-target signals → loudspeaker
            % outputs):  the contralateral leakage from each loudspeaker
            % attenuates by Δa and delays by Δt.  Defining
            %     c(ω) = Δa · exp(-jω·Δt)
            % the forward acoustic matrix [from out_L,out_R to ear_L,ear_R] is
            %     [1  c]
            %     [c  1]
            % Inverting this gives the RACE compensation filters:
            %     [out_L]   1     [ 1  -c] [in_L]
            %     [out_R] = ─── · [-c   1] [in_R]
            %               1-c²
            Nfft = obj.race_fft_length;
            Nf   = Nfft/2 + 1;
            freq = (0:Nf-1)' * obj.fs / Nfft;
            omega = 2*pi * freq;

            c = obj.race_delta_a * exp(-1i * omega * obj.race_delta_t);
            denom = 1 - c.^2;

            H_direct = 1 ./ denom;        % [Nf x 1]: in_L -> out_L, in_R -> out_R
            H_cross  = -c ./ denom;       % [Nf x 1]: in_L -> out_R, in_R -> out_L

            obj.race_H_direct = H_direct;
            obj.race_H_cross  = H_cross;

            % ── Convert to time-domain FIR (real, causal) ────────────────
            % We use ifft with 'symmetric' to enforce real impulse responses,
            % then circ-shift so the bulk of the energy is centred.
            h_direct = fftshift(ifft(H_direct, Nfft, 1, 'symmetric'));
            h_cross  = fftshift(ifft(H_cross,  Nfft, 1, 'symmetric'));

            % Mild Tukey windowing reduces edge ringing
            win = tukeywin(Nfft, 0.1);
            h_direct = h_direct .* win;
            h_cross  = h_cross  .* win;

            % ── Build OLS convolvers ─────────────────────────────────────
            blk = size(obj.virtual_source.source_signal.time_series, 1);
            if blk == 0, blk = 1024; end

            obj.race_filter_LL = OLS_convolver(h_direct, blk);  % in_L -> out_L
            obj.race_filter_LR = OLS_convolver(h_cross,  blk);  % in_R -> out_L
            obj.race_filter_RL = OLS_convolver(h_cross,  blk);  % in_L -> out_R
            obj.race_filter_RR = OLS_convolver(h_direct, blk);  % in_R -> out_R
        end

        function render_race(obj, bL, bR)
            % Apply the 2x2 RACE filter matrix to the [bL; bR] desired-ear
            % signal pair, and write the two outputs to the two outermost
            % loudspeakers. All other loudspeakers are silenced.

            if isempty(obj.race_filter_LL) || isempty(obj.race_filter_LR) || ...
               isempty(obj.race_filter_RL) || isempty(obj.race_filter_RR)
                return;  % filters not built yet
            end

            % Make sure the inputs are columns of equal length
            bL = bL(:);
            bR = bR(:);
            blk_len = max(numel(bL), numel(bR));
            if numel(bL) < blk_len, bL(end+1:blk_len) = 0; end
            if numel(bR) < blk_len, bR(end+1:blk_len) = 0; end

            % 2x2 mixing
            out_L = obj.race_filter_LL.convolve(bL) + ...
                    obj.race_filter_LR.convolve(bR);
            out_R = obj.race_filter_RL.convolve(bL) + ...
                    obj.race_filter_RR.convolve(bR);

            % Write to the two outermost loudspeakers; zero the rest.
            Nls = numel(obj.ssd.loudspeakers);
            zero_block = zeros(numel(out_L), 1);
            for k = 1:Nls
                if k == obj.race_idx_L
                    obj.ssd.loudspeakers{k}.source_signal.set_signal(out_L, obj.fs);
                elseif k == obj.race_idx_R
                    obj.ssd.loudspeakers{k}.source_signal.set_signal(out_R, obj.fs);
                else
                    obj.ssd.loudspeakers{k}.source_signal.set_signal(zero_block, obj.fs);
                end
            end
        end

        function build_multiband_crossover(obj)
            % Build a complementary linear-phase FIR pair (LP, HP) such that
            %     h_lp + h_hp = δ        (perfect reconstruction)
            % at the user-set crossover frequency. The HP is derived as the
            % spectral inverse of the LP via h_hp = δ - h_lp, which gives
            % exact summation in the time domain.

            Nfir = obj.mb_filter_length;
            if mod(Nfir, 2) == 1, Nfir = Nfir + 1; end   % even, so peak is well-defined
            fc = obj.mb_crossover_hz;
            fc = max(50, min(0.45*obj.fs, fc));           % keep in valid range

            % Linear-phase FIR low-pass via Hamming-windowed sinc.
            % This avoids the Signal Processing Toolbox dependency that
            % fir1 would impose, and gives identical output for our purpose.
            n = (0:Nfir-1)' - (Nfir-1)/2;                 % centred index
            wc = 2 * fc / obj.fs;                         % normalised cutoff (0..1)
            % Sinc with a small-eps guard at n=0
            ideal = wc * sinc(wc * n);
            % Hamming window
            win = 0.54 - 0.46*cos(2*pi*(0:Nfir-1)'/(Nfir-1));
            h_lp = ideal .* win;
            % Normalise DC gain to 1
            h_lp = h_lp / sum(h_lp);

            % Complementary high-pass: h_hp(n) = δ(n - Nfir/2) - h_lp(n)
            % The peak of the linear-phase LP sits at sample Nfir/2 + 1
            % (1-based MATLAB indexing, since n=0 maps to index Nfir/2 + 1).
            h_hp = -h_lp;
            mid = Nfir/2 + 1;
            h_hp(mid) = h_hp(mid) + 1;

            obj.mb_h_lp = h_lp;
            obj.mb_h_hp = h_hp;

            % Build OLS convolvers for both bands, both channels
            blk = size(obj.virtual_source.source_signal.time_series, 1);
            if blk == 0, blk = 1024; end
            obj.mb_conv_lp_L = OLS_convolver(h_lp, blk);
            obj.mb_conv_lp_R = OLS_convolver(h_lp, blk);
            obj.mb_conv_hp_L = OLS_convolver(h_hp, blk);
            obj.mb_conv_hp_R = OLS_convolver(h_hp, blk);
        end

        function render_multiband(obj, bL, bR)
            % Multiband CTC rendering:
            %   1. Split bL/bR into low and high bands via complementary FIR pair
            %   2. Low band → RACE (2x2 matrix), output to two outermost speakers
            %   3. High band → LSFI (per-element FIR), output to all 8 speakers
            %   4. Sum: outermost speakers receive RACE_low + LSFI_high
            %
            % Because the FIR pair is complementary (h_lp + h_hp = δ exactly)
            % and both filter chains are linear-phase, the band sum at each
            % loudspeaker reconstructs faithfully.

            if isempty(obj.mb_conv_lp_L) || isempty(obj.race_filter_LL) ...
                    || isempty(obj.lsfi_filters_L)
                return;
            end

            % ── Crossover split ─────────────────────────────────────────
            bL = bL(:); bR = bR(:);
            blk_len = max(numel(bL), numel(bR));
            if numel(bL) < blk_len, bL(end+1:blk_len) = 0; end
            if numel(bR) < blk_len, bR(end+1:blk_len) = 0; end

            bL_lo = obj.mb_conv_lp_L.convolve(bL);
            bR_lo = obj.mb_conv_lp_R.convolve(bR);
            bL_hi = obj.mb_conv_hp_L.convolve(bL);
            bR_hi = obj.mb_conv_hp_R.convolve(bR);

            % ── RACE on the low band (writes to outermost speakers, zeros others) ──
            % We capture the RACE outputs locally instead of writing them to the
            % SSD directly, so that we can sum them with the LSFI high-band
            % contributions later.
            race_out_L = obj.race_filter_LL.convolve(bL_lo) + ...
                         obj.race_filter_LR.convolve(bR_lo);
            race_out_R = obj.race_filter_RL.convolve(bL_lo) + ...
                         obj.race_filter_RR.convolve(bR_lo);

            % ── LSFI on the high band (writes to all 8 speakers) ────────
            % The LSFI per-element filters live in obj.lsfi_filters_L/R as
            % cell arrays of OLS_convolver. Each element n receives a sum
            % of the L-beam and R-beam contributions.
            Nls = numel(obj.ssd.loudspeakers);
            lsfi_out = cell(1, Nls);
            for n = 1:Nls
                ynL = obj.lsfi_filters_L{n}.convolve(bL_hi);
                ynR = obj.lsfi_filters_R{n}.convolve(bR_hi);
                lsfi_out{n} = ynL + ynR;
            end

            % ── Sum the two band contributions at each loudspeaker ──────
            for k = 1:Nls
                lsfi_part = lsfi_out{k};
                if k == obj.race_idx_L
                    sig = race_out_L + lsfi_part;
                elseif k == obj.race_idx_R
                    sig = race_out_R + lsfi_part;
                else
                    sig = lsfi_part;
                end
                obj.ssd.loudspeakers{k}.source_signal.set_signal(sig, obj.fs);
            end
        end

        function render_single_beam(obj, wfs, inputSig)
            % Feed input (first channel) into one beam and render.
            % For focus (point-source) mode, normalise the WFS output by the
            % free-field amplitude  1/(4π·r)  so the level stays constant
            % regardless of how close the beam target is to the listener.
            mono = inputSig(:,1);
            wfs.virtual_source.source_signal.set_signal(mono, obj.fs);
            wfs.render();

            if any(strcmpi(obj.wfs_mode, {'focus', 'directive'}))
                gain = obj.focus_normalisation_gain(wfs);
                obj.apply_gain_to_ssd(gain);
            end
        end

        function render_ctc(obj)
            % Full CTC render (no cross-talk cancellation matrix yet):
            %   earL = sum_m   HRTF_L(vs_m) * input_m
            %   earR = sum_m   HRTF_R(vs_m) * input_m
            % Then wfsL renders earL, wfsR renders earR.
            %
            % Each virtual source in the scene contributes its own
            % signal, convolved with the HRIR at that source's actual
            % position (looked up via hrtf_extrap).

            Nvs = numel(obj.scene_virtual_sources);
            if Nvs == 0
                % Fallback: route the local mono input to both ear channels.
                sig = obj.virtual_source.source_signal.time_series;
                if strcmpi(obj.wfs_mode, 'race')
                    obj.render_race(sig, sig);
                elseif strcmpi(obj.wfs_mode, 'multiband')
                    obj.render_multiband(sig, sig);
                else
                    obj.render_single_beam(obj.wfsL, sig);
                    obj.render_single_beam(obj.wfsR, sig);
                end
                return;
            end

            % Determine block size from the first VS signal
            blk = length(obj.scene_virtual_sources{1}.source_signal.time_series);

            % (Re)build HRTF convolvers if the scene changed
            if isempty(obj.hrtf_convL) || numel(obj.hrtf_convL) ~= Nvs
                obj.build_hrtf_convolvers(blk);
            end

            % Mix ear signals (desired ear signals d = [dL; dR])
            earL = zeros(blk, 1);
            earR = zeros(blk, 1);
            for m = 1:Nvs
                vs_m = obj.scene_virtual_sources{m};
                sig_m = vs_m.source_signal.time_series;
                if isempty(sig_m), continue; end
                if size(sig_m, 2) > 1, sig_m = sig_m(:, 1); end
                if length(sig_m) ~= blk
                    if length(sig_m) < blk
                        sig_m = [sig_m; zeros(blk-length(sig_m), 1)]; %#ok<AGROW>
                    else
                        sig_m = sig_m(1:blk);
                    end
                end
                if m <= numel(obj.hrtf_convL) && ~isempty(obj.hrtf_convL{m})
                    earL = earL + obj.hrtf_convL{m}.convolve(sig_m);
                end
                if m <= numel(obj.hrtf_convR) && ~isempty(obj.hrtf_convR{m})
                    earR = earR + obj.hrtf_convR{m}.convolve(sig_m);
                end
            end

            % Apply CTC plant inversion if enabled:
            %   [bL; bR] = C^{-1}(ω) * [earL; earR]
            % where C is the 2x2 ear-transfer matrix of the two beams.
            if obj.CTCInverseEnabled
                if isempty(obj.ctc_inv_filters)
                    obj.build_ctc_inverse_filters(blk);
                end
                if ~isempty(obj.ctc_inv_filters)
                    bL = obj.ctc_inv_filters.LL.convolve(earL) + ...
                         obj.ctc_inv_filters.LR.convolve(earR);
                    bR = obj.ctc_inv_filters.RL.convolve(earL) + ...
                         obj.ctc_inv_filters.RR.convolve(earR);
                else
                    bL = earL; bR = earR;
                end
            else
                bL = earL; bR = earR;
            end

            % Feed beam drive signals into the beam renderers
            if strcmpi(obj.wfs_mode, 'lsfi')
                % LSFI: convolve through per-element FIR filters
                obj.render_lsfi('L', bL);
                obj.render_lsfi('R', bR);
            elseif strcmpi(obj.wfs_mode, 'race')
                % RACE: 2x2 IIR-style filter applied to the [bL; bR] pair,
                % output goes to the two outermost loudspeakers only.
                obj.render_race(bL, bR);
            elseif strcmpi(obj.wfs_mode, 'multiband')
                % Multiband: split [bL; bR] into low and high bands via
                % complementary FIR pair. Low band → RACE (2 outermost LS),
                % high band → LSFI (all 8 LS). The two outputs are summed
                % at each loudspeaker.
                obj.render_multiband(bL, bR);
            else
                % WFS: feed into internal wfs_renderer
                G0l = obj.focus_normalisation_gain(obj.wfsL)/4/pi;
                obj.wfsL.virtual_source.source_signal.set_signal(bL*G0l, obj.fs);
                obj.wfsL.render();

                G0r = obj.focus_normalisation_gain(obj.wfsL)/4/pi;
                obj.wfsR.virtual_source.source_signal.set_signal(bR*G0r, obj.fs);
                obj.wfsR.render();
            end
        end

        function build_ctc_inverse_filters(obj, blk)
            % Build OLS convolvers for the Tikhonov-regularised 2x2 plant
            % inverse C^{-1}(ω), using the current ear transfer model.
            %
            % Plant matrix convention (frequency-domain):
            %   [dL(ω)]   [ H_LL(ω)  H_RL(ω) ] [bL(ω)]
            %   [dR(ω)] = [ H_LR(ω)  H_RR(ω) ] [bR(ω)]
            %
            %   H_LL: L beam -> left  ear   (desired on L)
            %   H_LR: L beam -> right ear   (crosstalk L->R)
            %   H_RL: R beam -> left  ear   (crosstalk R->L)
            %   H_RR: R beam -> right ear   (desired on R)
            %
            % Regularised inverse:
            %   C_inv = (C^H * C + β² I)^{-1} * C^H     (per frequency bin)
            %
            % The resulting 4 complex transfer functions are IFFT-ed into
            % length-Nfft FIR filters, wrapped in OLS_convolvers of block
            % size `blk`.
            try
                [f, H_LL, H_LR, H_RL, H_RR] = obj.compute_ear_transfer( ...
                    obj.ear_transfer_model, obj.hrtf_extrap); %#ok<ASGLU>
            catch ME
                warning('beamformer:ctcInverseBuild', ...
                    'Failed to compute ear transfer for CTC inverse: %s', ...
                    ME.message);
                obj.ctc_inv_filters = [];
                return;
            end

            Nf   = numel(H_LL);
            Nfft = 2 * (Nf - 1);
            beta = max(obj.ctc_beta, eps);

            % Build per-frequency inverse
            Ci_LL = zeros(Nf, 1);
            Ci_LR = zeros(Nf, 1);
            Ci_RL = zeros(Nf, 1);
            Ci_RR = zeros(Nf, 1);
            I2 = eye(2);
            for k = 1:Nf
                C = [H_LL(k), H_RL(k); ...
                     H_LR(k), H_RR(k)];
                % Tikhonov regularised pseudo-inverse
                Cinv = (C' * C + (beta^2) * I2) \ C';
                Ci_LL(k) = Cinv(1, 1);
                Ci_LR(k) = Cinv(1, 2);
                Ci_RL(k) = Cinv(2, 1);
                Ci_RR(k) = Cinv(2, 2);
            end

            % IFFT each one-sided spectrum to a length-Nfft FIR
            hLL = fftshift(ifft(Ci_LL, Nfft, 1, 'symmetric'), 1);
            hLR = fftshift(ifft(Ci_LR, Nfft, 1, 'symmetric'), 1);
            hRL = fftshift(ifft(Ci_RL, Nfft, 1, 'symmetric'), 1);
            hRR = fftshift(ifft(Ci_RR, Nfft, 1, 'symmetric'), 1);

            % Wrap into OLS convolvers (one per entry)
            obj.ctc_inv_filters = struct( ...
                'LL', OLS_convolver(hLL, blk), ...
                'LR', OLS_convolver(hLR, blk), ...
                'RL', OLS_convolver(hRL, blk), ...
                'RR', OLS_convolver(hRR, blk));
            fprintf('CTC plant inverse built (Nfft=%d, β=%.4g, model=%s)\n', ...
                Nfft, beta, obj.ear_transfer_model);
        end

        function gain = focus_normalisation_gain(obj, wfs)
            % Compute the inverse free-field amplitude for a focus beam:
            %   nominal_amp = 1 / (4π · |vs_pos - rec_pos|)
            %   gain        = 1 / nominal_amp  =  4π · r
            % Applying this gain to the rendered loudspeaker signals keeps
            % the perceived level at the receiver independent of distance.
            if isempty(obj.receiver)
                gain = 1;
                return;
            end
            vs_pos  = wfs.virtual_source.position(:).';
            rec_pos = obj.receiver.position(:).';
            % Pad to equal length
            n = max(numel(vs_pos), numel(rec_pos));
            if numel(vs_pos)  < n, vs_pos  = [vs_pos,  zeros(1, n-numel(vs_pos))];  end
            if numel(rec_pos) < n, rec_pos = [rec_pos, zeros(1, n-numel(rec_pos))]; end
            r    = norm(vs_pos - rec_pos);
            r    = max(r, 0.01);          % floor at 1 cm to avoid division by zero
            gain = 4 * pi * r;            % inverse of 1/(4π·r)
            %gain = 1;
        end

        function apply_gain_to_ssd(obj, gain)
            % Multiply every loudspeaker signal in the SSD by a scalar gain.
            for n = 1:numel(obj.ssd.loudspeakers)
                sig = obj.ssd.loudspeakers{n}.source_signal.time_series;
                if ~isempty(sig)
                    obj.ssd.loudspeakers{n}.source_signal.set_signal(sig *gain/4/pi, obj.fs);
                end
            end
        end

        function build_hrtf_convolvers(obj, blk)
            % Build OLS convolvers for each scene VS: HRIR from the VS
            % position to the left and right ears of the receiver.
            Nvs = numel(obj.scene_virtual_sources);
            obj.hrtf_convL = cell(1, Nvs);
            obj.hrtf_convR = cell(1, Nvs);

            if isempty(obj.hrtf_extrap)
                warning('beamformer:noHRTF', ...
                    'No HRTF extrapolator set. CTC mode needs obj.hrtf_extrap.');
                return;
            end
            if isempty(obj.receiver)
                warning('beamformer:noReceiver', ...
                    'No receiver set on beamformer_renderer.');
                return;
            end

            for m = 1:Nvs
                vs_m = obj.scene_virtual_sources{m};
                try
                    hrir = obj.hrtf_extrap.extrapolate_hrtf(vs_m, obj.receiver);
                    % hrir orientation: ensure [N x 2]
                    if size(hrir, 2) ~= 2 && size(hrir, 1) == 2
                        hrir = hrir.';
                    end
                    obj.hrtf_convL{m} = OLS_convolver(hrir(:, 1), blk);
                    obj.hrtf_convR{m} = OLS_convolver(hrir(:, 2), blk);
                catch ME
                    fprintf('HRTF convolver build failed for VS %d: %s\n', ...
                        m, ME.message);
                end
            end
            fprintf('CTC HRTF convolvers built for %d virtual sources\n', Nvs);
        end
    end

    % ═══════════════════════════════════════════════════════
    %  Static helpers
    % ═══════════════════════════════════════════════════════
    methods (Static, Access = private)
        function y = third_octave_smooth(f, x)
            y = zeros(size(x));
            for k = 1:numel(f)
                if f(k) <= 0, y(k) = x(k); continue; end
                f_lo = f(k) * 2^(-1/6);
                f_hi = f(k) * 2^(1/6);
                mask = (f >= f_lo) & (f <= f_hi);
                if any(mask), y(k) = mean(x(mask));
                else, y(k) = x(k); end
            end
        end
    end

    % ═══════════════════════════════════════════════════════
    %  Optimal position search (uses L beam)
    % ═══════════════════════════════════════════════════════
    methods
        function result = run_optimal_search(obj, f_low, f_high, acoustic_model, hrtf_extrap, grid_Nx, grid_Ny, metric)
            if nargin < 4 || isempty(acoustic_model), acoustic_model = 'analytic'; end
            if nargin < 5, hrtf_extrap = []; end
            if nargin < 6 || isempty(grid_Nx), grid_Nx = 40; end
            if nargin < 7 || isempty(grid_Ny), grid_Ny = 40; end
            if nargin < 8 || isempty(metric),  metric  = 'icld'; end

            search_mode = obj.wfs_mode;
            x0 = obj.ssd.x0;
            rec_pos = obj.receiver.position(:)';
            rec_ori = obj.receiver.orientation(:)';
            head_radius = 0.0875;

            % Search grid — right half-plane (x >= 0) only;
            % the L beam will be placed at the mirror position (-x, y).
            x_min = 0;
            x_max = max(x0(:,1));
            y_min_grid = rec_pos(2) - 0.05;
            y_max_grid = mean(x0(:,2)) + 0.5;
            margin_x = max((x_max - x_min) * 0.3, 0.1);
            x_grid = linspace(x_min, x_max + margin_x, grid_Nx);
            y_grid = linspace(y_min_grid, y_max_grid, grid_Ny);
            [X, Y] = meshgrid(x_grid, y_grid);

            rs = struct('Renderer_setup', struct( ...
                'BeamformerMode', search_mode, ...
                'Tapering',       obj.wfs_setup.Tapering, ...
                'Antialiasing',   obj.wfs_setup.Antialiasing));

            suppression_map = nan(grid_Ny, grid_Nx);
            hw = waitbar(0, 'Scanning grid (R beam, right half-plane)...');
            total = grid_Ny * grid_Nx;
            count = 0;
            block_len = size(obj.virtual_source.source_signal.time_series, 1);
            if block_len == 0, block_len = 1024; end

            for iy = 1:grid_Ny
                for ix = 1:grid_Nx
                    count = count + 1;
                    if mod(count, 50) == 0
                        waitbar(count/total, hw, sprintf('Scanning... %d/%d', count, total));
                    end

                    % R beam position on right half-plane; L beam mirror
                    beamR_pos = [X(iy,ix), Y(iy,ix), 0];
                    beamL_pos = [-X(iy,ix), Y(iy,ix), 0];

                    try
                        % Create a temporary beamformer at these beam positions
                        vs_type = struct('Shape', 'point_source', 'R', 0.025);
                        tmp_vs = virtual_source(1, beamR_pos, [0,-1,0], vs_type, 'beamformer');
                        tmp_vs.set_input(zeros(block_len,1), obj.fs);
                        r = beamformer_renderer(tmp_vs, obj.ssd, rs, obj.receiver);
                        r.beamL_pos = beamL_pos;
                        r.beamR_pos = beamR_pos;
                        r.update_renderer();

                        [f_local, H_LL, H_LR, H_RL, H_RR] = r.compute_ear_transfer(acoustic_model, hrtf_extrap);
                        f_mask_local  = (f_local >= f_low) & (f_local <= f_high);
                        f_band_local  = f_local(f_mask_local);

                        switch lower(metric)
                            case 'icld'
                                % Crosstalk suppression (ICLD):
                                %   L beam: desired = H_LL, crosstalk = H_LR
                                %   R beam: desired = H_RR, crosstalk = H_RL
                                ICLD_L = 20*log10((abs(H_LL(f_mask_local))+1e-15) ./ (abs(H_LR(f_mask_local))+1e-15));
                                ICLD_R = 20*log10((abs(H_RR(f_mask_local))+1e-15) ./ (abs(H_RL(f_mask_local))+1e-15));
                                ICLD_combined = min(abs(ICLD_L), abs(ICLD_R));
                                ICLD_smooth = beamformer_renderer.third_octave_smooth(f_band_local, ICLD_combined);
                                suppression_map(iy, ix) = min(ICLD_smooth);

                            case 'condition'
                                % Plant condition number quality:
                                % For each freq bin compute σ_min/σ_max of
                                % C = [H_LL H_RL; H_LR H_RR].
                                % The reciprocal condition number 1/κ ∈ [0,1].
                                % Higher = better conditioned = better CTC.
                                Nbin = sum(f_mask_local);
                                rcond_vals = zeros(Nbin, 1);
                                idx_f = find(f_mask_local);
                                for kf = 1:Nbin
                                    fi = idx_f(kf);
                                    C = [H_LL(fi), H_RL(fi); ...
                                         H_LR(fi), H_RR(fi)];
                                    sv = svd(C);
                                    rcond_vals(kf) = sv(end) / (sv(1) + 1e-30);
                                end
                                rcond_smooth = beamformer_renderer.third_octave_smooth( ...
                                    f_band_local, rcond_vals);
                                % Worst-case reciprocal condition in the band
                                suppression_map(iy, ix) = min(rcond_smooth);
                        end
                    catch
                        continue;
                    end
                end
            end
            close(hw);

            % Validity mask: exclude inside head, behind-array points
            dist_to_head = sqrt((X - rec_pos(1)).^2 + (Y - rec_pos(2)).^2);
            inside_head = dist_to_head < head_radius * 1.2;

            % Also check the mirrored L position is outside the head
            dist_mirror = sqrt(((-X) - rec_pos(1)).^2 + (Y - rec_pos(2)).^2);
            inside_head = inside_head | (dist_mirror < head_radius * 1.2);

            n_hat2 = rec_ori(1:2) / max(norm(rec_ori(1:2)), 1e-15);
            perp2 = [-n_hat2(2), n_hat2(1)];
            ear_L_pos = rec_pos(1:2) + head_radius * perp2;
            ear_R_pos = rec_pos(1:2) - head_radius * perp2;

            valid_mask = ~inside_head;
            search_map = suppression_map;
            search_map(~valid_mask) = -inf;
            [opt_val, idx] = max(search_map(:));
            [iy_opt, ix_opt] = ind2sub(size(search_map), idx);
            opt_pos_R = [X(iy_opt, ix_opt), Y(iy_opt, ix_opt), 0];
            opt_pos_L = [-X(iy_opt, ix_opt), Y(iy_opt, ix_opt), 0];

            result.X = X;
            result.Y = Y;
            result.suppression_map = suppression_map;
            result.valid_mask = valid_mask;
            result.opt_pos   = opt_pos_R(1:2);  % backward compat
            result.opt_pos_R = opt_pos_R;
            result.opt_pos_L = opt_pos_L;
            result.opt_val = opt_val;
            result.f_low = f_low;
            result.f_high = f_high;
            result.mode = search_mode;
            result.acoustic_model = acoustic_model;
            result.metric = metric;
            result.x0 = x0;
            result.rec_pos = rec_pos;
            result.head_radius = head_radius;
            result.ear_L_pos = ear_L_pos;
            result.ear_R_pos = ear_R_pos;
        end
    end
end
