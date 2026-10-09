
classdef wfs_renderer < base_renderer

    properties
        setup
        kh
        khn
        H_pre
        amp
        delay
        dref
        F
        WFS_filter_bank
        tapering_window
        antialiasing_filter_bank
        zones
        zone
        driving_spectrum   % [Nfft/2+1 × Nls] one-sided freq-domain driving function

        % Filter bank: built on every update (default), or only when
        % render() needs it (LazyFilterBank, set by sound_scene_renderer:
        % the fast binaural path only uses driving_spectrum)
        LazyFilterBank     = false
        bankDirty          = false

        % Low-frequency compensation (Firtha, Limitations of WFS Part II):
        % below the two half-window cutoffs f_c1, f_c2 (Eq. 34) each half of
        % the SSD stops contributing its full stationary-phase share; the
        % deficit is equalised by 1/sum_i(1/2 s(f/f_ci)), boost capped
        LowCutEnabled      = false
        LowCutMaxBoostDb   = 15     % [dB] cap of the low-frequency boost
        receiver_ref       = []     % receiver object for cutoff freq computation
        lower_cutoff_freq  = 0      % max(f_c1, f_c2) [Hz] (display)
        lower_cutoff_freqs = []     % [f_c1 f_c2] [Hz], Eq. 34 per half window
        higher_cutoff_freq = inf
    end

    properties (SetAccess = private)
        % Driving FIRs may be longer than the block (secondary sources far
        % from the virtual source): t = 0 stays at sample N/2+1, the tail is
        % extended so no delay wraps around. FIRLength = their length.
        FIRLength = 0
    end

    methods
        function obj = wfs_renderer(virtual_source,loudspeaker_array,setup)
            obj = obj@base_renderer(virtual_source,loudspeaker_array);
            obj.setup = setup;
            N = length(obj.virtual_source.source_signal.time_series);
            obj.H_pre = obj.get_wfs_prefilter(N);

            for n = 1 : length(obj.ssd.loudspeakers)
                obj.output_signal{n} = signal;
                obj.WFS_filter_bank{n} =  OLS_convolver(zeros(N,1),length(obj.virtual_source.source_signal.time_series));
            end
            obj.FIRLength = N;
            obj.update_renderer;

        end

        function info = get_renderer_info(obj)
            %GET_RENDERER_INFO  Human-readable description of the WFS renderer.

            lines = [
                "Wave Field Synthesis (WFS) renderer — overview"
                ""
                "WHAT THIS RENDERER DOES"
                "• Implements generalized 2.5D WFS driving functions with explicit referencing."
                "  The renderer computes per-loudspeaker delays and amplitudes from the"
                "  stationary-phase point (local wavenumber match) and applies  WFS"
                "  prefilter √(jω). Amplitude is corrected to be exact on a user-chosen"
                "  reference curve ('d_ref'), following the generalized framework in [1, 2]."
                ""
                "• Referencing schemes supported"
                "  - concentric: reference curve is a scaled copy of the SSD. For linear"
                "    arrays this reduces to an offset line; for circular arrays, a smaller"
                "    concentric circle; for general enclosing contours we intersect the local"
                "    ray (±kh) with the scaled SSD (closed shapes only) — cf. [1,2]."
                "  - fixed_point: amplitude is exact at a single point (classic revisited-WFS"
                "    style), which induces a reference curve in the plane, see [2]."
                "  - fixed_distance: amplitude is exact at a constant offset from each SSD"
                "    element (useful with 2D target fields); see [1,2]."
                ""
                "• Validity & limitations (theory reminders)"
                "  - Low-frequency limit: accuracy degrades below a geometry-dependent"
                "    cut-off tied to SPA validity / Fresnel-zone arguments — see [4]."
                "    Low frequency limit is highlighted in 'Frequency Response' mode."
                "  - Diffraction & truncation: finite-aperture terms add error outside the"
                "    referenced zone(s); analysis and examples are discussed in [3]."
                "  - Discretization: spatial aliasing is mitigated here by a per-element"
                "    anti-aliasing low-pass (local spacing criterion); asymptotic guidance"
                "    and practical filtering are summarized in [2, §Thesis II.3]."
                ""
                "ADJUSTABLE PARAMETERS (Renderer_setup)"
                "• ReferenceMode  : 'concentric' | 'fixed_point' | 'fixed_distance'"
                "• Rref           : radius/scale for concentric referencing (closed SSDs)"
                "• DeltaY         : offset for linear/open-contour 'concentric' cases"
                "• RefPoint       : [x y] for 'fixed_point' referencing"
                "• RefDistance    : scalar distance for 'fixed_distance' referencing"
                "• Tapering       : 0…1 Tukey window factor (reduces edge diffraction)"
                "• Antialiasing   : true/false; enables local low-pass per element"
                ""
                "REFERENCES"
                "[1] G. Firtha, P. Fiala, F. Schultz, S. Spors, “Improved Referencing Schemes"
                "    for 2.5D Wave Field Synthesis Driving Functions,” IEEE/ACM TASLP,"
                "    25(5):1117–1127, 2017."
                "[2] G. Firtha, “A generalized Wave Field Synthesis framework with"
                "    application for moving virtual sources,” PhD dissertation, BME, 2019."
                "[3] G. Firtha, “Limitations of Wave Field Synthesis Part I: Diffraction"
                "    Effects and Amplitude Errors,” 2025."
                "[4] G. Firtha, “Limitations of Wave Field Synthesis Part II: Low-frequency"
                "    limits,” 2025."
                ];
            info = strjoin(lines, newline);
        end

        function obj = update_renderer(obj, type)
            if nargin > 1 && strcmp(type, 'receiver_rotated')
                return;   % the driving functions do not depend on head orientation
            end
            N = length(obj.virtual_source.source_signal.time_series);

            [obj.amp,obj.delay] = obj.get_delay_and_gain;

            % Cutoff frequencies (geometry-dependent, always computed for display)
            obj.lower_cutoff_freqs = [];
            obj.lower_cutoff_freq  = 0;
            obj.higher_cutoff_freq = inf;
            if ~isempty(obj.receiver_ref)
                try
                    fc_lo = obj.get_lower_cutoff_frequency(obj.receiver_ref);
                    obj.lower_cutoff_freqs = fc_lo(:).';
                    fin = fc_lo(isfinite(fc_lo) & fc_lo > 0);
                    if ~isempty(fin), obj.lower_cutoff_freq = max(fin); end
                    obj.higher_cutoff_freq = obj.get_higher_cutoff_frequency(obj.receiver_ref);
                catch
                    obj.lower_cutoff_freqs = [];
                    obj.lower_cutoff_freq  = 0;
                    obj.higher_cutoff_freq = inf;
                end
            end

            obj.driving_spectrum = obj.compute_driving_spectrum(N);  % one-sided [N/2+1 x Nls]

            if obj.LazyFilterBank
                obj.bankDirty = true;      % built by render() if it is ever needed
            else
                obj.build_filter_bank();
            end
            obj.filter_version = obj.filter_version + 1;
        end

        function D = compute_driving_spectrum(obj, Nfft)
            % One-sided driving spectra [Nfft/2+1 x Nls] on an Nfft-point
            % grid at the source sampling rate (Nfft = block length gives
            % driving_spectrum). Part I/II Eq. 19b: D = w sqrt(8 pi i k)
            % sqrt(d) k_n P_F, with sqrt(8 pi i k) = sqrt(8 pi/c) sqrt(i omega)
            % (sqrt(8 pi/c)/(4 pi) is in amp, e^{+i omega t} convention).
            fs = obj.virtual_source.source_signal.fs;
            omega = (0 : Nfft/2)'/Nfft*2*pi*fs;
            H = obj.get_wfs_prefilter(Nfft);
            if obj.LowCutEnabled
                H = H.*obj.low_frequency_compensation(omega);
            end
            D = (H*obj.amp.').*exp(-1i*omega*obj.delay');
            if obj.setup.Antialiasing
                D = obj.get_antialiasing_filters(obj.delay, omega).*D;
            end
        end

        function G = low_frequency_compensation(obj, omega)
            % Part II: each half window (Eq. 32) contributes its share of the
            % stationary-phase level above its cutoff f_ci (Eq. 34) and grows
            % as sqrt(f) below it (Eq. 24). Relative level of the field:
            %   r(f) = mean_i s(f/f_ci),   s(x) = sqrt(x)/(1 + x^2)^(1/4)
            % (s = smooth min(1, sqrt(x))). Below both cutoffs r is the
            % Eq. 26 asymptote sqrt(k) (L1 + L2) / (stationary-phase level),
            % between them the +5 dB/dec transition of Sec. IV. G = 1/r,
            % capped at LowCutMaxBoostDb (so the DC bin stays zero).
            % Non-finite cutoffs (no stationary point) are ignored.
            G = ones(size(omega));
            fc = obj.lower_cutoff_freqs;
            fc = fc(isfinite(fc) & fc > 0);
            if isempty(fc), return; end
            x = omega(:)./(2*pi*fc(:).');
            s = sqrt(x)./(1 + x.^2).^(1/4);
            gmax = 10^(obj.LowCutMaxBoostDb/20);
            G = min(1./max(mean(s, 2), eps), gmax);
        end

        function build_filter_bank(obj)
            % Per-loudspeaker FIRs of the per-loudspeaker chain
            d_wfs = obj.get_driving_filters();
            if size(d_wfs,1) ~= obj.FIRLength       % FIR length changed: new convolvers
                N = length(obj.virtual_source.source_signal.time_series);
                for n = 1 : length(obj.WFS_filter_bank)
                    obj.WFS_filter_bank{n} = OLS_convolver(zeros(size(d_wfs,1),1), N);
                end
                obj.FIRLength = size(d_wfs,1);
            end
            for n = 1 : length(obj.WFS_filter_bank)
                obj.WFS_filter_bank{n}.update_coefficients(d_wfs(:,n));
            end
            obj.bankDirty = false;
        end

        function d = get_driving_filters(obj)
            % Driving FIRs [L x Nls] in SSD order, t = 0 at sample N/2+1
            % (N = block length; as fftshift of an N-point IFFT). L = N
            % unless a secondary source is so far from the virtual source
            % that its delay would wrap around: then the spectrum is
            % evaluated on a longer grid and only the tail is extended (same
            % latency for every source, no circular wrap). N/8 samples are
            % kept free for the prefilter / anti-aliasing filter tails.
            N  = length(obj.virtual_source.source_signal.time_series);
            fs = obj.virtual_source.source_signal.fs;
            on = abs(obj.amp) > 0;
            dmax = max([0; obj.delay(on)*fs]);
            L = N;
            while dmax + N/8 > L - N/2
                L = 2*L;
            end
            if L == N
                D = obj.driving_spectrum;
            else
                D = obj.compute_driving_spectrum(L);
            end
            d = circshift(ifft(D, L, 1, 'symmetric'), N/2, 1);
        end

        function obj = update_settings(obj,setup)
            obj.setup = setup;
            obj.update_renderer;
        end

        function render(obj)
            if obj.bankDirty, obj.build_filter_bank(); end
            for n = 1 : length(obj.WFS_filter_bank)
                obj.output_signal{n}.set_signal( obj.WFS_filter_bank{n}.convolve( obj.virtual_source.source_signal.time_series ) );
            end
            obj.add_output_to_ssd_signal;
        end

        function H = get_wfs_prefilter(obj, N)
            % One-sided spectrum [N/2+1 x 1] of the common WFS prefilter
            % sqrt(1i*omega) on an N-point grid.
            % Strategies:
            %   - windowing: samples of the analytical response (default)
            %   - optimal: optimal FIR design as given by Schultz in Sound
            %   Field Synthesis for Line Source Array Applications in
            %   Large-Scale Sound Reinforcement PhD thesis, Section 2.5
            %   - nara: FIR of wfs_prefilter_fir
            % The FIR variants are aligned with their peak tap at t = 0 and
            % scaled to |sqrt(1i*omega)| at 1 kHz.
            if nargin < 2, N = size(obj.virtual_source.source_signal.time_series,1); end
            fs = obj.virtual_source.source_signal.fs;
            omega = (0 : N/2)'/N*2*pi*fs;
            Nfilt = 128;
            method = 'windowing';
            switch method
                case 'windowing'
                    H = sqrt(1i*omega);
                    return;
                case 'optimal'
                    n = (-Nfilt/2:Nfilt/2-1)';
                    h = (sqrt(8*n).*cos(pi*n)+ ...
                        erfz( exp(1i*3/4*pi).*sqrt(pi*n) )-erfz( exp(1i*pi/4).*sqrt(pi*n) )).*...
                        (4*sqrt(pi)*n.^(3/2)).^(-1);
                    h(Nfilt/2+1) = sqrt(2*pi/9);
                    h = kaiser(Nfilt,4).*real(h);
                case 'nara'
                    [h2, ~] = wfs_prefilter_fir('2.5D', Nfilt, 20, 20e3, fs, audioapp.util.Physics.speedOfSound());
                    h = h2(1:end-1).';
            end
            h = h(:);
            hp = zeros(N,1); hp(1:numel(h)) = h;
            [~, ic] = max(abs(h));
            hp = circshift(hp, -(ic-1));
            H = fft(hp); H = H(1:N/2+1);
            [~, i1k] = min(abs(omega - 2*pi*1000));
            H = H*abs(sqrt(1i*omega(i1k)))/max(abs(H(i1k)), eps);
        end


        function [amp,delay] = get_delay_and_gain(obj)
            % 2.5D WFS driving function, Firtha, Limitations of WFS Part I,
            % Eq. (driv_fun): D = w sqrt(8 pi i k) sqrt(rhoP rhoG/(rhoP + F rhoG))
            % k_n P_F, rhoG = |x_ref - x0|, F = 1 diverging, -1 focused.
            x0 = obj.ssd.x0;
            Ns = size(x0, 1);
            D  = size(x0, 2);       % 2 or 3
            n0 = obj.ssd.n0;
            dl = obj.ssd.dA;
            c  = audioapp.util.Physics.speedOfSound();
            vs_pos = obj.vs_position(D);
            alpha = obj.setup.Tapering;

            switch obj.virtual_source.source_type.Shape
                case 'plane_wave'
                    obj.F = 1;
                    kP = -vs_pos/norm(vs_pos);
                    kh3 = repmat(kP, Ns, 1);
                    [obj.kh, obj.khn] = obj.in_plane_ray(kh3, n0);
                    win = double(obj.khn >= 0);                 % Eq. (window_fun)
                    obj.get_reference_distance(inf(Ns,1));
                    delay = sum(x0.*kh3, 2)/c;
                    amp = 1/4/pi*sqrt(8*pi/c).*win.*sqrt(obj.dref).*obj.khn.*dl;
                case 'point_source'
                    kP = bsxfun(@minus, x0, vs_pos);
                    rhoP = sqrt(sum(kP.^2,2));      % 3D distance (elevated sources)
                    [obj.kh, obj.khn] = obj.in_plane_ray(kP./rhoP, n0);
                    obj.F = 1 - 2*obj.is_focused(vs_pos);
                    if obj.F == 1
                        win = double(obj.khn >= 0);             % Eq. (window_fun)
                    else
                        % Focused: x0 radiates through x_s towards the
                        % listening area x_l, k^P.(x_l - x_s) <= 0 (Part I,
                        % after Eq. window_fun), and into the SSD's front
                        win = double(obj.khn <= 0);
                        xl = obj.listener_position(D);
                        if ~isempty(xl)
                            win = win.*double(obj.kh*(xl - vs_pos).' < 0);
                        end
                        if isfield(obj.setup, 'TaperingFocused')
                            alpha = obj.setup.TaperingFocused;  % Part I: focused sources need a smooth taper
                        end
                    end
                    obj.get_reference_distance(rhoP);

                    % rhoP + F*dref -> 0 for a focused source on the
                    % reference curve: keep it finite (5 cm floor, sign
                    % preserved). Past the focus the root is imaginary
                    % (Part I footnote: the +-pi/2 phase of the focus).
                    den = obj.F*obj.dref + rhoP;
                    small = abs(den) < 0.05;
                    den(small) = 0.05*(sign(den(small)) + (den(small) == 0));
                    dcorr = sqrt(obj.dref.*rhoP./den);
                    amp = sqrt(8*pi/c).*win.*dcorr.*obj.khn.*dl./rhoP/4/pi;
                    delay = obj.F*rhoP / c;
            end
            amp(~(win > 0)) = 0;

            obj.tapering_window = obj.arc_length_taper(win, alpha);
            amp = amp.*obj.tapering_window;
        end

        function get_reference_distance(obj, rhoP)
            % dref = |x_ref(x0) - x0| along the local ray F*kh (Part I,
            % Eq. driv_fun). rhoP: distance of x0 from the virtual source
            % (inf for plane waves).
            % Focused sources: the reference point must lie past the focus,
            % otherwise the field past the focus has a -90 deg error (Part I,
            % Eq. WFS_foc_amp_dist); a reference point before the focus is
            % mirrored across it (same distance from the focus).
            Ns = size(obj.ssd.x0, 1);
            type = lower(obj.ssd.type);
            if ~any(strcmp(type, {'linear','circular','general_enclosing'}))
                if obj.ssd.isClosed, type = 'general_enclosing'; else, type = 'linear'; end
            end
            switch lower(obj.setup.ReferenceMode)
                case 'concentric'
                    switch type
                        case 'linear'
                            % parallel reference line at DeltaY: rhoG = DeltaY/k_n
                            dref = obj.setup.DeltaY ./ wfs_renderer.safe_div(obj.F*obj.khn);
                        case 'circular'
                            dref = obj.dref_concentric_circle(rhoP);
                        case 'general_enclosing'
                            if obj.ssd.isClosed
                                dref = obj.compute_dref_concentric_general(rhoP);
                            else
                                dref = obj.setup.DeltaY ./ wfs_renderer.safe_div(obj.F*obj.khn);
                            end
                    end

                case 'fixed_point'
                    rp = obj.setup.RefPoint;
                    D = size(obj.ssd.x0, 2);
                    rp = [rp(:)', zeros(1, max(0, D-numel(rp)))];
                    dref = sqrt(sum( (obj.ssd.x0-rp(1:D)).^2,2 ));

                case 'fixed_distance'
                    dref = obj.setup.RefDistance*ones(Ns, 1);
            end
            dref = abs(real(dref(:)));
            dref(~isfinite(dref)) = 0;
            if obj.F == -1
                before = dref < rhoP;
                dref(before) = 2*rhoP(before) - dref(before);
            end
            obj.dref = dref;
        end

        function bank = get_antialiasing_filters(obj, tau, omega)
            % Per-loudspeaker anti-aliasing filters [numel(omega) x Nls]
            if nargin < 3
                N = length(obj.virtual_source.source_signal.time_series);
                omega = (0 : N/2)'/N*2*pi*obj.virtual_source.source_signal.fs;
            end
            x0 = obj.ssd.x0;
            D  = size(x0, 2);
            vs = obj.vs_position(D);
            dx1 = (circshift(x0,-1)-x0);
            x1 = x0+dx1/2;
            dx2 = (circshift(x0,1)-x0);
            x2 = x0+dx2/2;
            if strcmp(obj.virtual_source.source_type.Shape, 'plane_wave')
                kpw = -vs/norm(vs);                     % plane wave: one direction
                kh1 = repmat(kpw, size(x0,1), 1); kh2 = kh1;
            else
                kh1 = (x1-vs)./sqrt(sum( (x1-vs).^2,2));
                kh2 = (x2-vs)./sqrt(sum( (x2-vs).^2,2));
            end
            dx_dir = max( [sum(kh1.*dx1,2),sum(kh2.*dx2,2)] ,[] ,2);
            R0 = mean(sqrt(sum(x0.^2,2)));
            [~,ix] = sort(x0*vs','descend');
            Nbut = 7;
            ixs_to_mod = ix(1:2);
            G = (x0(ixs_to_mod,:)')\vs'/norm(vs)*R0;
            ix0 = find(abs(G-1)<1e-9);
            if ~isempty(ix0)
                ixs_to_mod = ixs_to_mod(ix0);
                G = G(ix0);
            end

            wc = pi.*audioapp.util.Physics.speedOfSound()./abs(dx_dir);
            [Wc,W] = meshgrid(wc,omega);

            bank = 1./sqrt( 1 + ( W./Wc ).^(2*Nbut));
            H0 =1./sqrt(1+(1i*W(:,ixs_to_mod)./Wc(:,ixs_to_mod) ));

            H1 = 1./sqrt(1+(W(:,ixs_to_mod)./Wc(:,ixs_to_mod) ).^(4*Nbut));
            G1 = H1 + (1-H1).*sqrt(G').*exp(1i*W(:,ixs_to_mod).*(tau(ixs_to_mod)'-mean(tau(ixs_to_mod))))./sqrt(1-abs(dx_dir(ixs_to_mod)).^2)';

            bank(:,ixs_to_mod) = H0.*G1;
            obj.antialiasing_filter_bank = bank;
        end

        function output = get_renderer_visuals(obj)
            output.x0        = obj.ssd.x0;
            output.n0        = obj.ssd.n0;
            output.win_taper = obj.tapering_window;
            output.amp       = obj.amp;
            output.kh        = obj.kh;
            output.ref_pos   = output.x0 + obj.F*output.kh.*abs(obj.dref);
            output.ref_pos( ~(abs(obj.amp)>0),: ) = nan;
            output.isClosed  = obj.ssd.isClosed;
            output.F         = obj.F;

            % NEW:
            if ~isempty(obj.zone)
                output.zone = obj.zone;   % array of structs with .boundaries [i1 i2], .type
            end

            % Back-compat (optional): keep old indices for any legacy code still using them
            if ~isempty(obj.zones)
                output.zones = obj.zones;
            end
        end

        function fc = get_lower_cutoff_frequency(obj, receiver)
            % Low-frequency limit of 2.5D WFS at the receiver (Firtha,
            % Limitations of WFS Part II): [f_c1 f_c2] of the two half
            % windows, Eq. 34: k_c = pi/2 d / (L_eff^2 k_n(x0*)^2), L_eff by
            % Eq. 32, d by Eq. 36 (point source) / Eq. 38 (plane wave).
            D = size(obj.ssd.x0, 2);
            vs = obj.vs_position(D);
            rec_pos = receiver.position;
            rec_pos = [rec_pos(:)', zeros(1, max(0, D-numel(rec_pos)))];
            rec_pos = rec_pos(1:D);
            if strcmpi(obj.virtual_source.source_type.Shape, 'plane_wave')
                % Eq. 38: x0* on the line through x along k^PW, upstream of x
                intersection = obj.ssd.get_ray_intersection( rec_pos, vs/norm(vs) );
            else
                switch obj.F
                    case 1
                        intersection = obj.ssd.get_ray_intersection( vs, rec_pos-vs );
                    case -1
                        intersection = obj.ssd.get_ray_intersection( vs, -(rec_pos-vs) );
                end
            end

            if ~isempty(intersection)
                % Upsample the SSD quantities to a dense parametric grid so
                % that the stationary-phase split integral is not dominated
                % by the coarse discrete sampling of a few loudspeakers.
                Nup = 500;
                x0      = obj.ssd.x0;
                dA      = obj.ssd.dA(:);
                dens_o  = obj.amp(:)./max(dA, eps);    % amp per unit arc length (amp holds dA)
                khn_o   = obj.khn(:);
                N = size(x0, 1);

                if obj.ssd.isClosed
                    % Closed loop: periodic interpolation
                    t_orig = (0:N-1)';
                    t_up   = linspace(0, N, Nup+1)'; t_up(end) = [];
                    x0_up = zeros(Nup, D);
                    for k = 1:D
                        x0_up(:,k) = interp1([t_orig; N], [x0(:,k); x0(1,k)], t_up, 'linear');
                    end
                    dens_up = interp1([t_orig; N], [dens_o; dens_o(1)], t_up, 'linear');
                    khn_up  = interp1([t_orig; N], [khn_o; khn_o(1)], t_up, 'linear');
                    dA_up   = sum(dA) / Nup * ones(Nup, 1);
                else
                    % Open contour: linear interpolation in arc length
                    t_orig = (1:N)';
                    t_up   = linspace(1, N, Nup)';
                    x0_up = zeros(Nup, D);
                    for k = 1:D
                        x0_up(:,k) = interp1(t_orig, x0(:,k), t_up, 'linear');
                    end
                    dens_up = interp1(t_orig, dens_o, t_up, 'linear');
                    khn_up  = interp1(t_orig, khn_o, t_up, 'linear');
                    dA_up   = sum(dA) / Nup * ones(Nup, 1);
                end

                integrand = abs(dens_up/4/pi ./ sqrt(sum((rec_pos - x0_up).^2, 2)));   % A^{D.G}, Eq. 23a

                % Stationary phase point on the upsampled grid
                [~, ix_stat] = min(sum((x0_up - intersection(1,1:D)).^2, 2));
                if integrand(ix_stat) < 1e-15
                    fc = inf; return;
                end

                % Eq. 32: L_eff = int A / A(x0*) ds from x0* to each end; the
                % stationary sample is shared half / half
                a = integrand .* dA_up / integrand(ix_stat);
                if obj.ssd.isClosed
                    a = circshift(a, -(ix_stat-1));              % a(1) = stationary sample
                    h = floor((Nup-1)/2);
                    I = [a(1)/2 + sum(a(2:1+h)), a(1)/2 + sum(a(2+h:end))];
                else
                    I = [sum(a(1:ix_stat-1)) + a(ix_stat)/2, a(ix_stat)/2 + sum(a(ix_stat+1:end))];
                end

                rhoG = norm(rec_pos - x0_up(ix_stat,:));
                switch lower(obj.virtual_source.source_type.Shape)
                    case 'plane_wave'
                        d_rec = rhoG;                                    % Eq. 38
                    case 'point_source'
                        rhoP = norm(vs - x0_up(ix_stat,:));
                        d_rec = rhoP*rhoG/(rhoP + obj.F*rhoG);           % Eq. 36
                end
                kc = pi/2 * d_rec ./ (I.^2 * khn_up(ix_stat).^2);       % Eq. 34
                fc = abs(kc * audioapp.util.Physics.speedOfSound() / 2 / pi);
            else
                fc = inf;
            end
        end

        function fc = get_higher_cutoff_frequency(obj, receiver)
            % Spatial aliasing limit at the receiver x. The integrand's phase
            % along the SSD changes at k |F k^P_t - k^G_t| (k^G = unit(x - x0),
            % t = tangent); sampling with dx creates a spurious stationary
            % point once this reaches 2 pi/dx, so
            %   f_al = min over the active elements of c / (dx |F k^P_t - k^G_t|)
            % (the mismatch is zero at the stationary point, at most 2).
            try
                x0 = obj.ssd.x0;
                D  = size(x0, 2);
                N  = size(x0, 1);
                rec = receiver.position;
                rec = [rec(:)', zeros(1, max(0, D-numel(rec)))]; rec = rec(1:D);
                % element tangents (central differences, one-sided at open ends)
                if obj.ssd.isClosed
                    tg = circshift(x0,-1) - circshift(x0,1);
                else
                    tg = [x0(2,:)-x0(1,:); x0(3:end,:)-x0(1:end-2,:); x0(end,:)-x0(end-1,:)];
                end
                tg = tg./max(sqrt(sum(tg.^2,2)), eps);
                kG = rec - x0; kG = kG./max(sqrt(sum(kG.^2,2)), eps);
                mis = abs(sum((obj.F*obj.kh(:,1:D) - kG).*tg, 2));
                on = abs(obj.amp(:)) > 0 & mis > 1e-9;
                if ~any(on) || N < 2, fc = inf; return; end
                dx = obj.ssd.dA(:);
                fc = min(audioapp.util.Physics.speedOfSound()./(dx(on).*mis(on)));
            catch
                fc = inf;
            end
        end

        function A0 = get_intended_amplitude(obj, X)
            % Amplitude deviation A0(x) = |P_int(x)| / |P(x)| of the intended
            % field (Part I, Eq. amp_distribution / WFS_amp_dist, focused:
            % WFS_foc_amp_dist) at the points X [M x 2]: 1 = amplitude-correct.
            % x0*(x): the SSD point whose local ray F*kh passes through x
            % (linear interpolation between neighbouring elements); NaN where
            % no illuminated element radiates towards x.
            X = X(:,1:2);
            M = size(X,1);
            x0 = obj.ssd.x0(:,1:2);
            N  = size(x0,1);
            A0 = nan(M,1);
            if N < 2 || isempty(obj.kh), return; end
            k  = obj.F*obj.kh(:,1:2);
            k  = k./max(sqrt(sum(k.^2,2)), eps);
            w  = obj.tapering_window(:);
            dr = obj.dref(:);
            on = abs(obj.amp(:)) > 0;
            isPW = strcmpi(obj.virtual_source.source_type.Shape, 'plane_wave');
            vs = obj.vs_position(2);
            if obj.ssd.isClosed, nxt = [2:N 1]; else, nxt = 2:N; end
            cur = 1:numel(nxt);
            for m = 1:M
                u = X(m,:) - x0;                              % x0 -> x
                cr = k(:,1).*u(:,2) - k(:,2).*u(:,1);         % signed angle to the ray
                fr = sum(k.*u, 2) > 0;                        % x in front of x0
                ok = on(cur) & on(nxt) & fr(cur) & fr(nxt) & cr(cur).*cr(nxt) <= 0;
                if ~any(ok), continue; end
                j = find(ok);
                [~, b] = min(abs(cr(cur(j))) + abs(cr(nxt(j))));
                i1 = cur(j(b)); i2 = nxt(j(b));
                lam = cr(i1)/(cr(i1) - cr(i2) + eps*(cr(i1) == cr(i2)));
                xs  = (1-lam)*x0(i1,:) + lam*x0(i2,:);
                ws  = (1-lam)*w(i1) + lam*w(i2);
                ds  = (1-lam)*dr(i1) + lam*dr(i2);
                rhoG = norm(X(m,:) - xs);
                if isPW
                    A0(m) = ws*sqrt(ds/max(rhoG, eps));
                else
                    rhoP = norm(xs - vs);
                    S1 = sqrt(abs(rhoP*ds/(rhoP + obj.F*ds)));
                    S2 = sqrt(abs(rhoP*rhoG/(rhoP + obj.F*rhoG)));
                    A0(m) = ws*S1*S2*norm(X(m,:) - vs)/(rhoP*rhoG);
                end
            end
        end
    end

    methods (Static)
        function [kh, khn] = in_plane_ray(kh3, n0)
            % Local ray in the plane of the SSD for a (possibly elevated)
            % source. The stationary-phase condition along the contour
            % matches the tangential component of the 3D ray kh3, so the
            % in-plane ray keeps kh3.t and takes the normal component
            % sqrt(1 - (kh3.t)^2) with the sign of kh3.n0. Delay and
            % amplitude still use the 3D distance. For an in-plane source
            % kh = kh3 exactly.
            D = size(kh3, 2);
            khn = sum(kh3.*n0, 2);
            if D < 3 || all(abs(kh3(:,3)) < 1e-12)
                kh = kh3;
                return;
            end
            nh = [n0(:,1:2), zeros(size(n0,1),1)];
            ln = sqrt(sum(nh.^2,2));
            ok = ln > 1e-6;                  % vertical normals: leave as is
            nh(ok,:) = nh(ok,:)./ln(ok);
            t = [-nh(:,2), nh(:,1), zeros(size(nh,1),1)];
            kt = sum(kh3.*t, 2);
            kn = sign(sum(kh3.*nh,2)).*sqrt(max(0, 1 - kt.^2));
            kh = kh3;
            kh(ok,:) = kt(ok).*t(ok,:) + kn(ok).*nh(ok,:);
            khn(ok) = kn(ok);
        end
    end

    methods (Access = private)
        function p = vs_position(obj, D)
            % virtual source position padded / cut to D coordinates
            p = obj.virtual_source.position;
            p = [p(:)', zeros(1, max(0, D-numel(p)))];
            p = p(1:D);
        end

        function tf = is_focused(obj, vs_pos)
            % Closed SSD: a source inside the contour is focused (also for
            % non-convex contours). Open SSD: focused if every element sees
            % it in front (k_n < 0 everywhere).
            if obj.ssd.isClosed
                X = obj.ssd.x0(:,1:2);
                tf = inpolygon(vs_pos(1), vs_pos(2), X(:,1), X(:,2));
            else
                tf = all(obj.khn < 0);
            end
        end

        function xl = listener_position(obj, D)
            % Centre of the listening area x_l for the focused-source
            % selection (Part I): the reference point if referencing is
            % 'fixed_point', the SSD centroid for closed arrays, otherwise
            % [] (open arrays: the k_n <= 0 window alone).
            xl = [];
            if strcmpi(obj.setup.ReferenceMode, 'fixed_point') && isfield(obj.setup, 'RefPoint')
                rp = obj.setup.RefPoint;
                xl = [rp(:)', zeros(1, max(0, D-numel(rp)))];
                xl = xl(1:D);
            elseif obj.ssd.isClosed
                xl = mean(obj.ssd.x0, 1);
            end
        end

        function w = arc_length_taper(obj, win, alpha)
            % Tukey window over the arc length s of every contiguous active
            % segment (Part I: w(s), tapered fraction alpha, alpha/2 of the
            % segment at each end; smooth to the second derivative). Sets
            % zone / zones from the longest segment for the overlay.
            N = numel(win);
            act = win(:) > 0;
            dA = obj.ssd.dA(:);
            w = zeros(N,1);
            obj.zone = []; obj.zones = [];
            if ~any(act), return; end
            closed = obj.ssd.isClosed;

            idx = (1:N)';
            if closed && ~all(act)
                idx = circshift(idx, -find(~act, 1));   % idx(end) is inactive: no wrapped segment
            end
            a = act(idx);
            e = diff([0; a; 0]);
            st = find(e == 1); en = find(e == -1) - 1;

            bestL = -1; best = []; bestT = [];
            for q = 1:numel(st)
                seg = idx(st(q):en(q));
                ds = dA(seg);
                s = cumsum(ds) - ds/2;
                L = sum(ds);
                t = ones(numel(seg), 1);
                at = alpha*L/2;
                if at > 0
                    de = min(s, L - s);                 % distance to the nearer end
                    in = de < at;
                    t(in) = 0.5*(1 - cos(pi*de(in)/at));
                end
                w(seg) = t.*win(seg);
                if L > bestL, bestL = L; best = seg; bestT = t; end
            end

            one = find(bestT >= 1 - 1e-9);
            if isempty(one), [~, one] = max(bestT); end
            i1 = best(one(1)); i2 = best(one(end));
            if closed && ~all(act)
                b0 = mod(best(1)-2, N) + 1;             % last inactive element before
                b3 = mod(best(end), N) + 1;             % first inactive element after
            else
                b0 = best(1); b3 = best(end);
            end
            obj.zone = struct('boundaries', {[b0 i1], [i1 i2], [i2 b3]}, ...
                'type', {'tapered', 'illuminated', 'tapered'});
            obj.zones = [b0 i1 i2 b3];
        end

        function dref = dref_concentric_circle(obj, rhoP)
            % Distance along F*kh to a circle of radius Rref concentric with
            % the SSD: diverging - the first crossing in front of x0;
            % focused - the first crossing past the focus (t > rhoP). A ray
            % missing the circle takes its closest approach to the centre
            % (real, instead of the complex root of the earlier formula).
            X = obj.ssd.x0(:,1:2);
            cc = mean(X, 1);
            d = obj.F*obj.kh(:,1:2);
            d = d./max(sqrt(sum(d.^2,2)), eps);
            q = X - cc;
            b = sum(q.*d, 2);
            disc = b.^2 - (sum(q.^2,2) - obj.setup.Rref^2);
            r = sqrt(max(disc, 0));
            t1 = -b - r; t2 = -b + r;
            dref = max(-b, 0);                          % closest approach
            if obj.F == 1, tmin = zeros(size(rhoP)); else, tmin = rhoP; end
            t = t1; t(t1 <= tmin) = t2(t1 <= tmin);
            ok = disc >= 0 & t > tmin;
            dref(ok) = t(ok);
        end

        function dref = compute_dref_concentric_general(obj, rhoP)
            % Scaled (concentric) copy of the SSD polygon; each ray
            % x0 + t F kh, t >= 0, is intersected with it: diverging - the
            % first crossing; focused - the first crossing past the focus.
            % No crossing: closest approach to the centroid.
            X  = obj.ssd.x0(:,1:2);
            Np = size(X,1);
            mu = mean(X,1);
            Rssd = mean(sqrt(sum((X-mu).^2,2)));
            s  = obj.setup.Rref / max(Rssd,1e-12);
            Xref = (X - mu) * s + mu;
            if norm(Xref(end,:) - Xref(1,:)) > 1e-12
                Xref = [Xref; Xref(1,:)];
            end
            kh = obj.local_kh_for_x0();
            Dr = obj.F * kh(:,1:2);
            if obj.F == 1, tmin = zeros(Np,1); else, tmin = rhoP(:); end
            dref = zeros(Np,1);
            for i = 1:Np
                p = X(i,:); d = Dr(i,:);
                if ~all(isfinite(d)) || norm(d) < 1e-12, continue; end
                d = d/norm(d);
                t = wfs_renderer.ray_poly_hits(p, d, Xref);
                t = sort(t(t > tmin(i) + 1e-9));
                if ~isempty(t)
                    dref(i) = t(1);
                else
                    dref(i) = max(-(p - mu)*d', 0);
                end
            end
        end

        function kh = local_kh_for_x0(obj)
            % Return N x D unit vectors kh at each SSD point
            X = obj.ssd.x0;
            D = size(X, 2);
            vs_pos = obj.vs_position(D);
            switch lower(obj.virtual_source.source_type.Shape)
                case 'plane_wave'
                    kP = zeros(1,D) - vs_pos;
                    dir = kP / max(norm(kP),1e-12);
                    kh = repmat(dir, size(X,1), 1);
                otherwise  % 'point_source' (and default)
                    kP = bsxfun(@minus, X, vs_pos);
                    n  = sqrt(sum(kP.^2,2));
                    n(n<1e-12) = 1e-12;
                    kh = kP ./ n;
            end
            kh = obj.in_plane_ray(kh, obj.ssd.n0);
        end

    end

    methods (Static, Access = private)
        function v = safe_div(v)
            % divisor kept away from 0 (sign kept, 0 -> +)
            s = sign(v); s(s == 0) = 1;
            v = s.*max(abs(v), 1e-9);
        end

        function t = ray_poly_hits(p, d, P)
            % all t >= 0 with p + t d on the polygon P (M x 2, closed)
            A = P(1:end-1,:); E = P(2:end,:) - A;
            den = d(1)*E(:,2) - d(2)*E(:,1);
            w = A - p;
            t = (w(:,1).*E(:,2) - w(:,2).*E(:,1))./den;
            u = (w(:,1)*d(2) - w(:,2)*d(1))./den;
            ok = abs(den) > 1e-14 & u >= -1e-12 & u <= 1 + 1e-12 & t >= 0;
            t = t(ok);
        end
    end

    % ─── Overlay drawing (delegates to WFSPainter) ───
    methods
        function store = drawOverlay(obj, ax, layers, settings, store)
            %DRAWOVERLAY  WFS-specific overlays: reference contour, amplitude
            %  distribution, tapering window, illuminated/tapered zones,
            %  intended-field amplitude map A0 (Part I).
            rp = obj.get_renderer_visuals();
            store = audioapp.ui.draw.WFSPainter.draw(ax, rp, obj, settings, store);
        end

        function settings = defaultOverlaySettings(~)
            settings = struct('showRef',true,'showAmp',true,'showTaper',true,'showZones',true,'showA0',false);
        end
    end

end

%nti xl3
