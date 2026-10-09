classdef time_delay_renderer < base_renderer
    % TIME_DELAY_RENDERER  Two-channel panning engine.
    %
    % Picks two most frontal SSD loudspeakers as a stereo base, then
    % applies delay and/or amplitude differences according to the chosen
    % panning law.
    %
    % Methods are classified as:
    %   DELAY-BASED:     phase/delay between speakers, equal amplitude
    %   AMPLITUDE-BASED: gain difference between speakers, no delay
    %   HYBRID:          both delay and amplitude, controlled by α
    %
    % Setup fields (Renderer_setup.*):
    %   TDNumTaps      : FIR length (default 2048)
    %   TDMethod       : panning law string (default 'firtha_linear')
    %   TDRefPoint     : [x y] sweet-spot (default [0 0])
    %   TDRefRadiusR0  : 'auto' or scalar (default 'auto')
    %   TDHybridAlpha  : 0..1 mix for hybrid mode (default 0.5)

    properties
        fs
        c
        N_taps    = 2048
        method    = 'firtha_linear'
        refPoint  = [0 0 0]
        R0_override = []
        panelCallback = []
        hybridAlpha = 0.5
    end

    properties (SetAccess = protected)
        renderer_filter
        stereo_pair = []
        fi_base     = 0
    end

    properties (Constant)
        % {id, display_name, type}
        METHODS = { ...
            'firtha_exact',   'Firtha (exact — arctan phase)',        'delay'; ...
            'firtha_linear',  'Firtha (linearized — pure delay)',     'delay'; ...
            'sengpiel',       'Sengpiel (empirical curve)',           'delay'; ...
            'linear',         'Linear delay interpolation',          'delay'; ...
            'sine_law',       'Sine law (Makita)',                   'amplitude'; ...
            'tangent_law',    'Tangent law (Blumlein)',               'amplitude'; ...
            'hybrid',         'Hybrid (delay + amplitude)',           'hybrid' ...
        };
    end

    methods
        function obj = time_delay_renderer(virtual_source, SSD, setup)
            obj = obj@base_renderer(virtual_source, SSD);

            obj.fs = obj.virtual_source.source_signal.fs;
            obj.c  = audioapp.util.Physics.speedOfSound();

            RS = getfielddef(setup, 'Renderer_setup', struct());
            obj.N_taps      = getfielddef(RS, 'TDNumTaps',      2048);
            obj.method      = getfielddef(RS, 'TDMethod',       'firtha_linear');
            obj.refPoint    = getfielddef(RS, 'TDRefPoint',     [0 0 0]);
            if numel(obj.refPoint) < 3, obj.refPoint = [obj.refPoint(:)', zeros(1, 3-numel(obj.refPoint))]; end
            obj.R0_override = getfielddef(RS, 'TDRefRadiusR0',  'auto');
            obj.hybridAlpha = getfielddef(RS, 'TDHybridAlpha',  0.5);

            obj.build_or_update_filters();
        end

        function info = get_renderer_info(obj)
            ix = find(strcmp(obj.METHODS(:,1), obj.method), 1);
            if isempty(ix), mname = obj.method; mtype = '?';
            else, mname = obj.METHODS{ix,2}; mtype = obj.METHODS{ix,3}; end

            lines = [
                "TWO-CHANNEL PANNING ENGINE"
                ""
                "Current method:  " + string(mname) + "  [" + upper(string(mtype)) + "]"
                "FIR taps:        " + string(obj.N_taps)
                "Ref point:       [" + num2str(obj.refPoint(1),'%.3f') + ", " + num2str(obj.refPoint(2),'%.3f') + "]"
                "R0:              " + string(ifempty(obj.R0_override, 'auto'))
                ""
                "═══════════════════════════════════════════"
                "DELAY-BASED METHODS"
                "═══════════════════════════════════════════"
                ""
                "1) Firtha (exact — arctan phase)"
                "   φ(ω) = ± arctan( ωR₀/c · tan(φ₀)/tan(φ_b) )"
                "   Frequency-dependent: ITD-like at low f, saturates"
                "   at high f. From local wavevector matching."
                "   Ref: Firtha G., PhD dissertation, BME, 2019"
                ""
                "2) Firtha (linearized — pure delay)"
                "   Δt = ± (R₀/4c) · tan(φ₀)/tan(φ_b)"
                "   Low-freq approximation of (1). Pure time delay."
                "   Ref: ibid."
                ""
                "3) Sengpiel (empirical curve)"
                "   Lagrange-interpolated from measured data:"
                "   b=0%→0ms, 25%→0.23ms, 50%→0.48ms,"
                "   75%→0.81ms, 100%→1.5ms  (±30° base)"
                "   Scaled to actual base angle."
                "   Ref: Sengpiel E., measured 1992 (25 subjects,"
                "        Teldec studio, speech + string quartet)"
                ""
                "4) Linear delay interpolation"
                "   Δt = ± Δt_max · φ₀/φ_b"
                "   Simplest model. Δt_max = d/(2c)."
                ""
                "═══════════════════════════════════════════"
                "AMPLITUDE-BASED METHODS"
                "═══════════════════════════════════════════"
                ""
                "5) Sine law (Makita velocity vector)"
                "   g₁/g₂ = sin(φ_b+φ₀) / sin(φ_b−φ₀)"
                "   Power-normalized. Based on velocity vector."
                "   Ref: Makita Y., EBU Review Part A, No. 73, 1962"
                ""
                "6) Tangent law (Blumlein)"
                "   g₁/g₂ = tan(φ_b+φ₀) / tan(φ_b−φ₀)"
                "   Power-normalized. Classic intensity stereo."
                "   Ref: Bennett et al., JAES 33(5), 1985"
                ""
                "═══════════════════════════════════════════"
                "HYBRID METHOD"
                "═══════════════════════════════════════════"
                ""
                "7) Hybrid (delay + amplitude)"
                "   Delay:  Δt = ±(1−α)·(d/2c)·tan(φ₀)/tan(φ_b)"
                "   Amp:    tan-law with effective angle α·φ₀"
                "   α=0: pure delay, α=1: pure amplitude."
                "   Ref: Williams M., AES Conv. 82, 1987"
                "        Theile G., JAES 39(10), 1991"
                ""
                "───────────────────────────────────────────"
                "φ₀ = VS angle from bisector, φ_b = half-base angle"
                "R₀ = ref radius, d = head diameter (0.17 m)"
            ];
            info = strjoin(lines, newline);
        end

        function update_renderer(obj, ~)
            obj.build_or_update_filters();
        end

        function obj = update_settings(obj, renderer_setup)
            needRebuild = false;
            if isfield(renderer_setup,'TDNumTaps') && isfinite(renderer_setup.TDNumTaps) ...
                    && renderer_setup.TDNumTaps > 0 && renderer_setup.TDNumTaps ~= obj.N_taps
                obj.N_taps = renderer_setup.TDNumTaps; needRebuild = true;
            end
            if isfield(renderer_setup,'TDMethod')
                m = char(renderer_setup.TDMethod);
                if ~strcmp(m, obj.method), obj.method = m; needRebuild = true; end
            end
            if isfield(renderer_setup,'TDRefPoint') && numel(renderer_setup.TDRefPoint) >= 2
                rp_new = renderer_setup.TDRefPoint(:)';
                if numel(rp_new) < 3, rp_new = [rp_new, 0]; end
                if any(rp_new ~= obj.refPoint)
                    obj.refPoint = rp_new; needRebuild = true;
                end
            end
            if isfield(renderer_setup,'TDRefRadiusR0')
                v = renderer_setup.TDRefRadiusR0;
                if (ischar(v)||isstring(v)) && strcmpi(string(v),'auto')
                    if ~isempty(obj.R0_override), obj.R0_override = []; needRebuild = true; end
                elseif isempty(v)
                    if ~isempty(obj.R0_override), obj.R0_override = []; needRebuild = true; end
                elseif isnumeric(v) && isfinite(v) && v > 0
                    if isempty(obj.R0_override) || obj.R0_override ~= v
                        obj.R0_override = v; needRebuild = true;
                    end
                end
            end
            if isfield(renderer_setup,'TDHybridAlpha')
                v = renderer_setup.TDHybridAlpha;
                if isnumeric(v) && isfinite(v) && v >= 0 && v <= 1 && v ~= obj.hybridAlpha
                    obj.hybridAlpha = v; needRebuild = true;
                end
            end
            if needRebuild, obj.build_or_update_filters(); end
        end

        function render(obj)
            x = obj.virtual_source.source_signal.time_series;
            for n = 1:numel(obj.renderer_filter)
                obj.output_signal{n}.set_signal(obj.renderer_filter{n}.convolve(x));
            end
            obj.add_output_to_ssd_signal;
        end

        % ── Public: panning curve for plotting ──
        function [angles_deg, dt_ms, g1_vec, g2_vec] = get_panning_curve(obj)
            % Returns panning law curve.  For delay methods dt_ms is the
            % delay curve, g1=g2=1.  For amplitude methods dt_ms=0,
            % g1/g2 are the gain curves.  For hybrid, both are nonzero.

            if isempty(obj.stereo_pair) || obj.fi_base <= 0
                angles_deg = []; dt_ms = []; g1_vec = []; g2_vec = [];
                return;
            end

            R0 = obj.get_R0();
            fb = obj.fi_base;
            d_head = 0.17;
            M = 200;
            fi0_vec = linspace(-fb, fb, M)';
            angles_deg = fi0_vec * 180/pi;

            dt_s   = zeros(M, 1);
            g1_vec = ones(M, 1);
            g2_vec = ones(M, 1);

            for k = 1:M
                fi0 = fi0_vec(k);
                s = sign(fi0);
                a = abs(fi0);
                [dt_s(k), g1_vec(k), g2_vec(k)] = obj.panning_law_eval(a, s, fb, R0, d_head);
            end

            dt_ms = dt_s * 1e3;
        end

        % ── Visuals ──
        function output = get_renderer_visuals(obj)
            output.x0          = obj.ssd.x0;
            output.vs          = obj.virtual_source.position;
            output.stereo_pair = obj.stereo_pair;
            output.fi_base     = obj.fi_base;
            output.refPoint    = obj.refPoint;
            output.method      = obj.method;
            output.method_type = obj.get_method_type();
            [output.delays_ms, output.gains_lin] = obj.get_speaker_params();
        end

        function store = drawOverlay(obj, ax, layers, settings, store)
            rp = obj.get_renderer_visuals();
            store = audioapp.ui.draw.TDPainter.draw(ax, rp, settings, store);
        end

        function settings = defaultOverlaySettings(~)
            settings = struct('showDelayLabels', true, 'showPairLink', true);
        end
    end

    methods
        function mtype = get_method_type(obj)
            ix = find(strcmp(obj.METHODS(:,1), obj.method), 1);
            if isempty(ix), mtype = 'delay'; else, mtype = obj.METHODS{ix,3}; end
        end
    end

    % =====================================================================
    % Private implementation
    % =====================================================================
    methods (Access = private)

        function [dt_s, g1, g2] = panning_law_eval(obj, fi0_abs, side, fb, R0, d_head)
            % Core panning law evaluation for a single angle.
            %   fi0_abs : absolute angle from bisector (rad)
            %   side    : +1 or -1
            %   Returns: dt_s (signed delay), g1, g2 (linear gains)

            dt_s = 0;
            g1 = 1;
            g2 = 1;

            switch obj.method
                % ── DELAY-BASED ──
                case 'firtha_exact'
                    f0 = 500; w0 = 2*pi*f0;
                    alpha_arg = R0/obj.c * tan(fi0_abs)/tan(fb);
                    dt_s = side * atan(w0 * alpha_arg) / w0;

                case 'firtha_linear'
                    dt_s = side * (R0/(4*obj.c)) * tan(fi0_abs)/tan(fb);

                case 'sengpiel'
                    % Lagrange interpolation of empirical data (Sengpiel 1992)
                    % b = normalized position 0..1, dt_ref for ±30° base
                    b = fi0_abs / fb;  % 0 at center, 1 at edge
                    b = min(max(b, 0), 1);
                    dt_ref = sengpiel_lagrange(b) * 1e-3;  % ms → s
                    dt_s = side * dt_ref;

                case 'linear'
                    dt_max = d_head / (2*obj.c);
                    dt_s = side * dt_max * fi0_abs / fb;

                % ── AMPLITUDE-BASED ──
                case 'sine_law'
                    fi0_signed = side * fi0_abs;
                    g1 = sin(fb + fi0_signed);
                    g2 = sin(fb - fi0_signed);
                    gnorm = sqrt(g1^2 + g2^2);
                    if gnorm > 0, g1 = g1/gnorm; g2 = g2/gnorm; end

                case 'tangent_law'
                    fi0_signed = side * fi0_abs;
                    g1 = tan(fb + fi0_signed);
                    g2 = tan(fb - fi0_signed);
                    gnorm = sqrt(g1^2 + g2^2);
                    if gnorm > 0, g1 = g1/gnorm; g2 = g2/gnorm; end

                % ── HYBRID ──
                case 'hybrid'
                    alpha = obj.hybridAlpha;
                    % Delay part (tangent law, scaled by 1-α)
                    dt_s = side * (1-alpha) * (d_head/(2*obj.c)) * tan(fi0_abs)/tan(fb);
                    % Amplitude part (tangent law, scaled by α)
                    fi0_amp = side * alpha * fi0_abs;
                    g1 = tan(fb + fi0_amp);
                    g2 = tan(fb - fi0_amp);
                    gnorm = sqrt(g1^2 + g2^2);
                    if gnorm > 0, g1 = g1/gnorm; g2 = g2/gnorm;
                    else, g1 = 1/sqrt(2); g2 = 1/sqrt(2); end
            end
        end

        function [delays_ms, gains_lin] = get_speaker_params(obj)
            % Per-speaker delay (ms) and gain for overlay display.
            Nls = obj.ssd.N;
            delays_ms = zeros(Nls, 1);
            gains_lin = ones(Nls, 1);
            if isempty(obj.stereo_pair) || numel(obj.stereo_pair) < 2, return; end

            ixs = obj.stereo_pair;
            xr = obj.refPoint(:)';
            xs = obj.virtual_source.position;
            x0 = obj.ssd.x0;
            R0 = obj.get_R0();
            fb = obj.fi_base;
            if fb <= 0, return; end

            base_vec = unitrow((x0(ixs(1),:) + x0(ixs(2),:) - 2*xr));
            vec_xs   = unitrow(xs - xr);
            fi0 = acos(max(-1, min(1, dot(base_vec, vec_xs))));
            cross_z = base_vec(1)*vec_xs(2) - base_vec(2)*vec_xs(1);
            side = sign(cross_z);

            [dt_s, g1, g2] = obj.panning_law_eval(fi0, side, fb, R0, 0.17);

            delays_ms(ixs(1)) = dt_s * 1e3;
            delays_ms(ixs(2)) = -dt_s * 1e3;
            gains_lin(ixs(1)) = g1;
            gains_lin(ixs(2)) = g2;
        end

        function R0 = get_R0(obj)
            x0 = obj.ssd.x0;
            xr = obj.refPoint(:)';
            if isempty(obj.R0_override) || (ischar(obj.R0_override) && strcmpi(string(obj.R0_override),'auto'))
                R0 = mean(sqrt(sum((x0 - xr).^2, 2)));
            else
                R0 = obj.R0_override;
            end
        end

        function build_or_update_filters(obj)
            x0 = obj.ssd.x0;
            xs = obj.virtual_source.position;
            H = obj.compute_td_filters(xs, x0);

            Nls = size(x0, 1);
            if isempty(obj.renderer_filter) || numel(obj.renderer_filter) ~= Nls
                obj.renderer_filter = cell(1, Nls);
                for n = 1:Nls
                    obj.output_signal{n} = signal;
                    obj.renderer_filter{n} = OLS_convolver(H(:,n), ...
                        length(obj.virtual_source.source_signal.time_series));
                end
            else
                for n = 1:Nls
                    obj.renderer_filter{n}.update_coefficients(H(:,n));
                end
            end

            if ~isempty(obj.panelCallback)
                try obj.panelCallback(); catch, end
            end
        end

        function td_filters = compute_td_filters(obj, xs, x0)
            xr = obj.refPoint(:)';
            N  = obj.N_taps;
            w  = (0:N-1)'/N * 2*pi*obj.fs;
            R0 = obj.get_R0();
            d_head = 0.17;

            Nls = size(x0, 1);
            td_filters = zeros(N, Nls);

            % ── Find the stereo pair that brackets the source in angle ──
            % All angles measured from ref point (xy-plane projection)
            ang_ls = atan2(x0(:,2) - xr(2), x0(:,1) - xr(1));
            ang_vs = atan2(xs(2) - xr(2), xs(1) - xr(1));

            % Sort loudspeakers by angle
            [ang_sorted, ixSort] = sort(ang_ls);

            % Find the two adjacent speakers that bracket the VS angle
            % (treating the sorted angles as a circular list)
            diffs = mod(ang_vs - ang_sorted + pi, 2*pi) - pi;  % signed angular diff
            ixs = [];
            for k = 1:Nls
                k_next = mod(k, Nls) + 1;
                a1 = ang_sorted(k);
                a2 = ang_sorted(k_next);
                % Angular span from a1 to a2 (ccw)
                span = mod(a2 - a1, 2*pi);
                % VS angle relative to a1 (ccw)
                vs_rel = mod(ang_vs - a1, 2*pi);
                if vs_rel <= span + 1e-12
                    ixs = [ixSort(k), ixSort(k_next)];
                    break;
                end
            end

            % Fallback: nearest two
            if isempty(ixs)
                [~, nearest] = sort(abs(mod(ang_ls - ang_vs + pi, 2*pi) - pi));
                ixs = nearest(1:2)';
            end

            obj.stereo_pair = ixs;

            u1 = unitrow(x0(ixs(1),:) - xr);
            u2 = unitrow(x0(ixs(2),:) - xr);
            fb = acos(max(-1, min(1, dot(u1, u2)))) / 2;
            obj.fi_base = fb;

            % Degenerate: zero base angle
            if fb < 1e-10
                td_filters(:, ixs(1)) = fftshift(ifft(ones(N,1), N, 'symmetric'));
                return;
            end

            base_vec = unitrow((x0(ixs(1),:) + x0(ixs(2),:) - 2*xr));
            vec_xs   = unitrow(xs - xr);
            fi0 = acos(max(-1, min(1, dot(base_vec, vec_xs))));
            cross_z = base_vec(1)*vec_xs(2) - base_vec(2)*vec_xs(1);
            side = sign(cross_z);

            mtype = obj.get_method_type();

            % ── Compute phase ramp (delay methods) ──
            switch obj.method
                case 'firtha_exact'
                    arg = w/obj.c * R0 * tan(fi0)/tan(fb);
                    phi = atan(arg);

                case 'firtha_linear'
                    phi = (1/4) * w/obj.c * R0 * tan(fi0)/tan(fb);

                case 'sengpiel'
                    b = min(max(fi0/fb, 0), 1);
                    dt_ref = sengpiel_lagrange(b) * 1e-3;
                    phi = w * dt_ref;

                case 'linear'
                    dt = (d_head/(2*obj.c)) * fi0/fb;
                    phi = w * dt;

                case 'hybrid'
                    alpha = obj.hybridAlpha;
                    dt = (1-alpha) * (d_head/(2*obj.c)) * tan(fi0)/tan(fb);
                    phi = w * dt;

                otherwise  % amplitude methods: no phase
                    phi = zeros(N, 1);
            end

            % ── Compute amplitude gains ──
            switch obj.method
                case 'sine_law'
                    fi0s = side * fi0;
                    g1 = sin(fb + fi0s);
                    g2 = sin(fb - fi0s);
                    gnorm = sqrt(g1^2 + g2^2);
                    if gnorm > 0, g1 = g1/gnorm; g2 = g2/gnorm; end

                case 'tangent_law'
                    fi0s = side * fi0;
                    g1 = tan(fb + fi0s);
                    g2 = tan(fb - fi0s);
                    gnorm = sqrt(g1^2 + g2^2);
                    if gnorm > 0, g1 = g1/gnorm; g2 = g2/gnorm; end

                case 'hybrid'
                    alpha = obj.hybridAlpha;
                    fi0_amp = side * alpha * fi0;
                    g1 = tan(fb + fi0_amp);
                    g2 = tan(fb - fi0_amp);
                    gnorm = sqrt(g1^2 + g2^2);
                    if gnorm > 0, g1 = g1/gnorm; g2 = g2/gnorm;
                    else, g1 = 1/sqrt(2); g2 = 1/sqrt(2); end

                otherwise  % delay methods: equal gain
                    g1 = 1;
                    g2 = 1;
            end

            % ── Build transfer functions ──
            if strcmp(mtype, 'amplitude')
                % Pure amplitude: FIR is a scaled delta
                h1 = zeros(N, 1); h1(N/2+1) = g1;
                h2 = zeros(N, 1); h2(N/2+1) = g2;
                td_filters(:, ixs(1)) = h1;
                td_filters(:, ixs(2)) = h2;
            else
                % Delay or hybrid: phase ramp + gain
                H1 = g1 * exp( 1i * side * phi);
                H2 = g2 * exp(-1i * side * phi);
                td_filters(:, ixs(1)) = real(fftshift(ifft(H1, N, 'symmetric')));
                td_filters(:, ixs(2)) = real(fftshift(ifft(H2, N, 'symmetric')));
            end
        end
    end
end

% ═════════════════════════════════════════════════════════════════════════
% Helper functions
% ═════════════════════════════════════════════════════════════════════════
function v = getfielddef(S, f, d)
    if isstruct(S) && isfield(S, f) && ~isempty(S.(f)), v = S.(f); else, v = d; end
end

function u = unitrow(v)
    v = v(:)';  % ensure row
    n = norm(v); if n > 0, u = v./n; else, u = zeros(size(v)); u(1) = 1; end
end

function s = ifempty(val, repl)
    if isempty(val), s = repl; else, s = val; end
end

function dt_ms = sengpiel_lagrange(b)
    % Sengpiel (1992) empirical localization curve for time-delay stereo.
    % Input:  b = normalized phantom source position, 0 (center) to 1 (edge)
    % Output: dt_ms = interchannel time difference in milliseconds
    %
    % Data measured with 25 subjects at ±30° stereo triangle, Teldec studio.
    % Lagrange interpolation through 5 data points.
    %
    % b:     0      0.25    0.50    0.75    1.00
    % Δt:    0 ms   0.23ms  0.48ms  0.81ms  1.50ms

    bp = [0, 0.25, 0.50, 0.75, 1.00];
    tp = [0, 0.23, 0.48, 0.81, 1.50];

    % Lagrange interpolation
    b = max(0, min(1, b));
    dt_ms = 0;
    for i = 1:5
        Li = 1;
        for j = 1:5
            if j ~= i
                Li = Li * (b - bp(j)) / (bp(i) - bp(j));
            end
        end
        dt_ms = dt_ms + tp(i) * Li;
    end
end
