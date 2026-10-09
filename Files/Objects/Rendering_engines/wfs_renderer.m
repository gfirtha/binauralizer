
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

        % Low-frequency pre-eq shelving
        LowCutEnabled      = false
        receiver_ref       = []     % receiver object for cutoff freq computation
        lower_cutoff_freq  = 0
        higher_cutoff_freq = inf
    end

    methods
        function obj = wfs_renderer(virtual_source,loudspeaker_array,setup)
            obj = obj@base_renderer(virtual_source,loudspeaker_array);
            obj.setup = setup;
            N = length(obj.virtual_source.source_signal.time_series);
            obj.H_pre = obj.get_wfs_prefilter;

            for n = 1 : length(obj.ssd.loudspeakers)
                obj.output_signal{n} = signal;
                obj.WFS_filter_bank{n} =  OLS_convolver(zeros(N,1),length(obj.virtual_source.source_signal.time_series));
            end
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
            omega = (0 : N/2)'/N*2*pi*obj.virtual_source.source_signal.fs;

            [obj.amp,obj.delay] = obj.get_delay_and_gain;

            % Always compute cutoff frequencies (geometry-dependent, always useful for display)
            if ~isempty(obj.receiver_ref)
                try
                    fc_lo = obj.get_lower_cutoff_frequency(obj.receiver_ref);
                    if numel(fc_lo) >= 2
                        obj.lower_cutoff_freq = max(fc_lo);
                    elseif numel(fc_lo) == 1
                        obj.lower_cutoff_freq = fc_lo;
                    else
                        obj.lower_cutoff_freq = 0;
                    end
                    obj.higher_cutoff_freq = obj.get_higher_cutoff_frequency(obj.receiver_ref);
                catch
                    obj.lower_cutoff_freq  = 0;
                    obj.higher_cutoff_freq = inf;
                end
            else
                obj.lower_cutoff_freq  = 0;
                obj.higher_cutoff_freq = inf;
            end

            % Apply low-freq pre-eq shelving if enabled:
            % WFS prefilter H_pre = sqrt(1i*omega) is flattened below fc.
            % Use H_shelf = sqrt(1i * max(omega, omega_c))
            if obj.LowCutEnabled && obj.lower_cutoff_freq > 0
                omega_c = 2*pi*obj.lower_cutoff_freq;
                omega_eff = max(omega, omega_c);
                H_pre_used = sqrt(1i*omega_eff);
            else
                H_pre_used = obj.H_pre;
            end

            if obj.setup.Antialiasing
                obj.get_antialiasing_filters(obj.delay);
                D_wfs = obj.antialiasing_filter_bank.* (H_pre_used*obj.amp.').*exp(-1i*omega*obj.delay');
            else
                D_wfs = (H_pre_used*obj.amp.').*exp(-1i*omega*obj.delay');
            end
            obj.driving_spectrum = D_wfs;  % store one-sided [Nfft/2+1 × Nls]

            if obj.LazyFilterBank
                obj.bankDirty = true;      % built by render() if it is ever needed
            else
                obj.build_filter_bank();
            end
            obj.filter_version = obj.filter_version + 1;
        end

        function build_filter_bank(obj)
            % Per-loudspeaker FIRs of the per-loudspeaker chain
            d_wfs = obj.get_driving_filters();
            for n = 1 : length(obj.WFS_filter_bank)
                obj.WFS_filter_bank{n}.update_coefficients(d_wfs(:,n));
            end
            obj.bankDirty = false;
        end

        function d = get_driving_filters(obj)
            % Same FIRs as loaded into WFS_filter_bank, [N x Nls] in SSD order
            N = length(obj.virtual_source.source_signal.time_series);
            d = fftshift(ifft(obj.driving_spectrum, N, 1, 'symmetric'), 1);
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

        function H = get_wfs_prefilter(obj)
            % The function calculates the common prefilter present in WFS
            % driving function, given be sqrt(1i*omega)
            % Strategies:
            %   - windowing: calculates prefilters directly from evaluation
            %   the analytical formula in the frequency domain and IFFT-s
            %   into the temporal domain applying windowing
            %   - optimal: implements optimal FIR filter design as given by
            %   Schultz in Sound Field Synthesis for Line Source Array
            %   Applications in Large-Scale Sound Reinforcement PhD thesis,
            %   in Section 2.5
            Nfilt = 128;
            method = 'windowing';
            switch method
                case 'windowing'
                    N = size(obj.virtual_source.source_signal.time_series,1);
                    omega = (0 : N/2)'/N*2*pi*obj.virtual_source.source_signal.fs;
                    H = sqrt(1i*omega);
                    %                    h = fftshift(ifft(H,N,'symmetric'));
                    %                    Nfilt = N;
                    %                    h = h(round(end/2)-Nfilt/2+1:round(end/2)+Nfilt/2).*tukeywin(Nfilt,0.1);
                case 'optimal'
                    n = (-Nfilt/2:Nfilt/2-1)';
                    h = (sqrt(8*n).*cos(pi*n)+ ...
                        erfz( exp(1i*3/4*pi).*sqrt(pi*n) )-erfz( exp(1i*pi/4).*sqrt(pi*n) )).*...
                        (4*sqrt(pi)*n.^(3/2)).^(-1);
                    h(end/2+1) = sqrt(2*pi/9);
                    h = kaiser(Nfilt,4).*real(h)*Nfilt;
                case 'nara'
                    [h2, ~] = wfs_prefilter_fir('2.5D', Nfilt, 20, 20e3, obj.fs, 343);
                    h = h2(1:end-1).'*Nfilt;
            end


        end


        function [amp,delay] = get_delay_and_gain(obj)
            x0 = obj.ssd.x0;
            D  = size(x0, 2);       % 2 or 3
            center = zeros(1, D);
            n0 = obj.ssd.n0;
            dl = obj.ssd.dA;
            c =  audioapp.util.Physics.speedOfSound();

            % Ensure VS position matches dimension
            vs_pos = obj.virtual_source.position;
            if numel(vs_pos) < D
                vs_pos = [vs_pos(:)', zeros(1, D-numel(vs_pos))];
            else
                vs_pos = vs_pos(1:D);
            end

            switch obj.virtual_source.source_type.Shape
                case 'plane_wave'
                    obj.F = 1;
                    kP = bsxfun( @minus, zeros(1,D), vs_pos );
                    kh3 = repmat( kP/norm(kP) , size(x0,1),1);
                    [obj.kh, obj.khn] = obj.in_plane_ray(kh3, n0);
                    win = double(obj.khn>=0);

                    obj.get_reference_distance;
                    dcorr = sqrt(obj.dref);
                    delay =  sum( (x0.*kh3) / c , 2);
                    amp = 1/4/pi*sqrt(8*pi/c).*win.*dcorr.*obj.khn.*dl;
                    ks = 1;
                case 'point_source'
                    kP = bsxfun( @minus, x0, vs_pos );
                    rhoP = sqrt(sum(kP.^2,2));      % 3D distance (elevated sources)
                    [obj.kh, obj.khn] = obj.in_plane_ray(kP./rhoP, n0);

                    ks = (vs_pos - center)/norm(vs_pos - center);
                    obj.F = 1-2*all(obj.khn<0); % Focused flag
                    if obj.F == 1
                        win = double(obj.khn>=0);
                    elseif obj.F == -1
                        if obj.ssd.isClosed
                            win = (obj.kh*ks').^0.*((obj.kh*ks')>0);
                        else
                            win = double(obj.khn<=0);
                        end
                    end
                    obj.get_reference_distance;

                    % F*dref + rhoP -> 0 for a focused source on the
                    % reference curve (e.g. at the reference point): keep
                    % it finite (5 cm floor, sign preserved)
                    den = obj.F*obj.dref + rhoP;
                    small = abs(den) < 0.05;
                    den(small) = 0.05*(sign(den(small)) + (den(small) == 0));
                    dcorr = sqrt(obj.dref.*rhoP./den);
                    amp = sqrt(8*pi/c).*win.*dcorr.*obj.khn.*dl./rhoP/4/pi;
                    delay = obj.F*rhoP / c;
            end
            % delay = delay - min(delay);

            if obj.ssd.isClosed
                [~,stat_ix] = max(abs(obj.khn).*(n0*ks'<0));
                stat_ix = stat_ix(1);
                win_taper = abs(circshift( win>0, -stat_ix));
                L1 = numel( find( win_taper(1:round( end/2 ))>0) );
                L2 = numel( find( win_taper(round( end/2 )+1:end)>0) );
                win1 = tukeywin(2*L1,obj.setup.Tapering);
                win2 = tukeywin(2*L2,obj.setup.Tapering);
                win_taper(1:L1) = win1(end/2+1:end);
                win_taper(end-L2+1:end) = win2(1:end/2);

                obj.zones = [find(win_taper == 0,1,'last'),...
                    find(win_taper(ceil(end/2)+1:end) == 1,1,'first') + ceil(length(win_taper)/2),...
                    find(win_taper(1:ceil(end/2)) == 1,1,'last') ,...
                    find(win_taper == 0,1,'first')];
                obj.zones = 1 + mod(obj.zones + stat_ix-1, size(obj.ssd.x0,1) );

                obj.zone(1).boundaries = [find(win_taper == 0,1,'last'),  find(win_taper(ceil(end/2)+1:end) == 1,1,'first') + ceil(length(win_taper)/2)];
                obj.zone(1).type = 'tapered';
                obj.zone(2).boundaries = [find(win_taper(ceil(end/2)+1:end) == 1,1,'first') + ceil(length(win_taper)/2), find(win_taper(1:ceil(end/2)) == 1,1,'last') ];
                obj.zone(2).type = 'illuminated';
                obj.zone(3).boundaries = [find(win_taper(1:ceil(end/2)) == 1,1,'last'), find(win_taper == 0,1,'first')];
                obj.zone(3).type = 'tapered';
                obj.zone(1).boundaries = 1 + mod(obj.zone(1).boundaries + stat_ix-1, size(obj.ssd.x0,1) );
                obj.zone(2).boundaries = 1 + mod(obj.zone(2).boundaries + stat_ix-1, size(obj.ssd.x0,1) );
                obj.zone(3).boundaries = 1 + mod(obj.zone(3).boundaries + stat_ix-1, size(obj.ssd.x0,1) );

                win_taper = circshift(win_taper, stat_ix-1);
                obj.tapering_window = win_taper.*win;


            else
                obj.tapering_window = tukeywin(size(obj.ssd.x0,1),obj.setup.Tapering);
                obj.zone(1).boundaries = [1, find(obj.tapering_window == 1,1,'first')];
                obj.zone(1).type = 'tapered';
                obj.zone(2).boundaries = [find(obj.tapering_window == 1,1,'first'), find(obj.tapering_window == 1,1,'last')];
                obj.zone(2).type = 'illuminated';
                obj.zone(3).boundaries = [find(obj.tapering_window == 1,1,'last'), length(obj.tapering_window)];
                obj.zone(3).type = 'tapered';
                obj.zones = [1, find(obj.tapering_window == 1,1,'first'), find(obj.tapering_window == 1,1,'last'), length(obj.tapering_window)];

            end
            amp = amp.*obj.tapering_window;


        end

        function get_reference_distance(obj)
            switch lower(obj.setup.ReferenceMode)
                case 'concentric'
                    switch lower(obj.ssd.type)
                        case 'linear'
                            % same as your current code
                            obj.dref = obj.F*obj.setup.DeltaY ./ obj.khn;

                        case 'circular'
                            % same as your current code
                            Rssd = mean( sqrt( sum( (obj.ssd.x0-obj.ssd.center).^2,2 ) ) );
                            obj.dref = obj.F*min( (Rssd.*(obj.khn + sqrt( obj.khn.^2 + (obj.setup.Rref/Rssd)^2 - 1) )),...
                                (Rssd.*(obj.khn - sqrt( obj.khn.^2 + (obj.setup.Rref/Rssd)^2 - 1) )) );

                        case 'general_enclosing'
                            % NEW: ray -> scaled copy of the SSD
                            if obj.ssd.isClosed
                                obj.dref = obj.compute_dref_concentric_general();
                            else
                                % Open general geometry: use the linear-style offset
                                obj.dref = obj.F*obj.setup.DeltaY ./ max(obj.khn, 1e-12);
                            end
                    end

                case 'fixed_point'
                    rp = obj.setup.RefPoint;
                    if numel(rp) < size(obj.ssd.x0, 2)
                        rp = [rp(:)', zeros(1, size(obj.ssd.x0,2)-numel(rp))];
                    end
                    obj.dref = sqrt(sum( (obj.ssd.x0-rp).^2,2 ));

                case 'fixed_distance'
                    obj.dref = obj.setup.RefDistance;
            end
        end

        function get_antialiasing_filters(obj,tau)

            x0 = obj.ssd.x0;
            dx1 = (circshift(x0,-1)-x0);
            x1 = x0+dx1/2;
            dx2 = (circshift(x0,1)-x0);
            x2 = x0+dx2/2;
            kh1 = (x1-obj.virtual_source.position)./sqrt(sum( (x1-obj.virtual_source.position).^2,2));
            kh2 = (x2-obj.virtual_source.position)./sqrt(sum( (x2-obj.virtual_source.position).^2,2));
            dx_dir = max( [sum(kh1.*dx1,2),sum(kh2.*dx2,2)] ,[] ,2);
            R0 = mean(sqrt(sum(x0.^2,2)));
            [~,ix] = sort(x0*obj.virtual_source.position','descend');
            Nbut = 7;
            ixs_to_mod = ix(1:2);
            G = (x0(ixs_to_mod,:)')\obj.virtual_source.position'/norm(obj.virtual_source.position)*R0;
            ix0 = find(abs(G-1)<1e-9);
            if ~isempty(ix0)
                ixs_to_mod = ixs_to_mod(ix0);
                G = G(ix0);
            end

            wc = pi.*audioapp.util.Physics.speedOfSound()./abs(dx_dir);

            N = length(obj.virtual_source.source_signal.time_series);
            w = (0 : N/2)'/N*2*pi*obj.virtual_source.source_signal.fs;
            [Wc,W] = meshgrid(wc,w);

            obj.antialiasing_filter_bank = 1./sqrt( 1 + ( W./Wc ).^(2*Nbut));
            H0 =1./sqrt(1+(1i*W(:,ixs_to_mod)./Wc(:,ixs_to_mod) ));

            H1 = 1./sqrt(1+(W(:,ixs_to_mod)./Wc(:,ixs_to_mod) ).^(4*Nbut));
            G1 = H1 + (1-H1).*sqrt(G').*exp(1i*W(:,ixs_to_mod).*(tau(ixs_to_mod)'-mean(tau(ixs_to_mod))))./sqrt(1-abs(dx_dir(ixs_to_mod)).^2)';

            obj.antialiasing_filter_bank(:,ixs_to_mod) = H0.*G1;
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
            switch obj.F
                case 1
                    intersection = obj.ssd.get_ray_intersection( obj.virtual_source.position, receiver.position-obj.virtual_source.position );
                case -1
                    intersection = obj.ssd.get_ray_intersection( obj.virtual_source.position,-(receiver.position-obj.virtual_source.position)  );
            end

            if ~isempty(intersection)
                % Upsample the SSD quantities to a dense parametric grid so
                % that the stationary-phase split integral is not dominated
                % by the coarse discrete sampling of a few loudspeakers.
                Nup = 500;
                x0      = obj.ssd.x0;
                dA      = obj.ssd.dA(:);
                amp_o   = obj.amp(:);
                khn_o   = obj.khn(:);
                N = size(x0, 1);

                if obj.ssd.isClosed
                    % Closed loop: periodic interpolation via circshift
                    t_orig = (0:N-1)';
                    t_up   = linspace(0, N, Nup+1)'; t_up(end) = [];
                    x0_up   = [interp1([t_orig; N], [x0(:,1); x0(1,1)], t_up, 'linear'), ...
                               interp1([t_orig; N], [x0(:,2); x0(1,2)], t_up, 'linear')];
                    if size(x0,2) >= 3
                        x0_up = [x0_up, interp1([t_orig; N], [x0(:,3); x0(1,3)], t_up, 'linear')];
                    end
                    amp_up  = interp1([t_orig; N], [amp_o; amp_o(1)], t_up, 'linear');
                    khn_up  = interp1([t_orig; N], [khn_o; khn_o(1)], t_up, 'linear');
                    % dA rescales to keep the total arc-length constant
                    dA_up   = sum(dA) / Nup * ones(Nup, 1);
                else
                    % Open contour: linear interpolation in arc length
                    t_orig = (1:N)';
                    t_up   = linspace(1, N, Nup)';
                    x0_up  = [interp1(t_orig, x0(:,1), t_up, 'linear'), ...
                              interp1(t_orig, x0(:,2), t_up, 'linear')];
                    if size(x0,2) >= 3
                        x0_up = [x0_up, interp1(t_orig, x0(:,3), t_up, 'linear')];
                    end
                    amp_up = interp1(t_orig, amp_o, t_up, 'linear');
                    khn_up = interp1(t_orig, khn_o, t_up, 'linear');
                    dA_up  = sum(dA) / Nup * ones(Nup, 1);
                end

                rec_pos = receiver.position;
                if size(x0_up,2) < numel(rec_pos)
                    rec_pos = rec_pos(1:size(x0_up,2));
                elseif size(x0_up,2) > numel(rec_pos)
                    rec_pos = [rec_pos(:)', zeros(1, size(x0_up,2)-numel(rec_pos))];
                end

                integrand = abs(amp_up/4/pi ./ sqrt(sum((rec_pos - x0_up).^2, 2)));

                % Stationary phase point on the upsampled grid
                [~, ix_stat] = min(sum((x0_up - intersection(1,1:size(x0_up,2))).^2, 2));
                if integrand(ix_stat) < 1e-15
                    fc = inf; return;
                end

                if obj.ssd.isClosed
                    integrand = circshift(integrand .* dA_up, -(ix_stat-1)) / integrand(ix_stat);
                    I = [sum(integrand(1:round(end/2))), sum(integrand(round(end/2)+1:end))];
                else
                    I = [sum(integrand(1:ix_stat)    .* dA_up(1:ix_stat)), ...
                         sum(integrand(ix_stat+1:end).* dA_up(ix_stat+1:end))] / integrand(ix_stat);
                end

                rhoP = norm(obj.virtual_source.position(1:size(x0_up,2)) - x0_up(ix_stat,:));
                rhoG = norm(rec_pos - x0_up(ix_stat,:));
                switch lower(obj.virtual_source.source_type.Shape)
                    case 'plane_wave'
                        d_rec = rhoG;
                    case 'point_source'
                        d_rec = rhoP*rhoG/(rhoP + obj.F*rhoG);
                end
                kc = pi/2 * d_rec ./ (I.^2 * khn_up(ix_stat).^2);
                fc = abs(kc * audioapp.util.Physics.speedOfSound() / 2 / pi);
            else
                fc = inf;
            end
        end

        function fc = get_higher_cutoff_frequency(obj, receiver)
            % Upper aliasing cutoff at the receiver position.
            % Follows the same derivation as in get_antialiasing_filters:
            %   wc = pi * c / |dx_dir|
            % where dx_dir is the projection of the inter-element spacing
            % onto the radiation direction (kh). Among all SSD elements we
            % pick the one closest to the listener-VS ray (stationary phase
            % point), which dominates the aliasing at that listener position.
            try
                x0 = obj.ssd.x0;
                dx1 = (circshift(x0,-1)-x0);
                x1 = x0 + dx1/2;
                dx2 = (circshift(x0,1)-x0);
                x2 = x0 + dx2/2;
                kh1 = (x1-obj.virtual_source.position) ./ sqrt(sum((x1-obj.virtual_source.position).^2,2));
                kh2 = (x2-obj.virtual_source.position) ./ sqrt(sum((x2-obj.virtual_source.position).^2,2));
                dx_dir = max([sum(kh1.*dx1,2), sum(kh2.*dx2,2)], [], 2);

                % Stationary phase point: SSD element closest to the
                % listener-VS ray
                switch obj.F
                    case 1
                        intersection = obj.ssd.get_ray_intersection( ...
                            obj.virtual_source.position, ...
                            receiver.position - obj.virtual_source.position);
                    case -1
                        intersection = obj.ssd.get_ray_intersection( ...
                            obj.virtual_source.position, ...
                            -(receiver.position - obj.virtual_source.position));
                end

                if ~isempty(intersection)
                    [~, ix_stat] = min(sum((obj.ssd.x0 - intersection(1,:)).^2, 2));
                else
                    % Fallback: element with the smallest dx_dir (densest)
                    [~, ix_stat] = max(abs(dx_dir));
                end

                wc = pi * audioapp.util.Physics.speedOfSound() / abs(dx_dir(ix_stat));
                fc = wc / (2*pi);
            catch
                fc = inf;
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
        function dref = compute_dref_concentric_general(obj)
            % Build a concentric (scaled) copy of the SSD polygon and intersect
            % each ray x0 + (F*kh)*t, t>=0, with that polygon. Return t as dref.

            X  = obj.ssd.x0;                 % N x 2
            Np = size(X,1);
            mu = mean(X,1);                  % should be ~[0 0], but compute anyway
            r  = sqrt(sum((X-mu).^2,2));
            Rssd = mean(r);                  % "average radius"
            s  = obj.setup.Rref / max(Rssd,1e-12);   % scale factor

            % Concentric reference contour (close loop)
            Xref = (X - mu) * s + mu;
            if norm(Xref(end,:) - Xref(1,:)) > 1e-12
                Xref = [Xref; Xref(1,:)];
            end

            % Build ray directions: F * kh at each SSD point
            kh = obj.local_kh_for_x0();      % N x 2 (unit)
            D  = obj.F * kh;                 % N x 2

            dref = nan(Np,1);
            tol = 1e-12;

            % Intersect each ray with the polygon
            for i = 1:Np
                p = X(i,:); d = D(i,:);
                if ~isfinite(d(1)) || ~isfinite(d(2)) || norm(d) < tol
                    dref(i) = 0; continue;
                end
                [hit,tmin, tnear] = obj.ssd.rayPolyIntersect(p,d,Xref);
                if hit
                    dref(i) = tmin;          % non-negative distance along F*kh
                else
                    % Fallback: try tiny numerical nudge along ray, otherwise 0
                    dref(i) = tnear;
                end
            end

            % Guarantee nonnegative (we march along F*kh)
            dref(~isfinite(dref) | dref<0) = 0;
        end

        function kh = local_kh_for_x0(obj)
            % Return N x D unit vectors kh at each SSD point
            X = obj.ssd.x0;
            D = size(X, 2);
            vs_pos = obj.virtual_source.position;
            if numel(vs_pos) < D, vs_pos = [vs_pos(:)', zeros(1, D-numel(vs_pos))]; end
            vs_pos = vs_pos(1:D);
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



    % ─── Overlay drawing (delegates to WFSPainter) ───
    methods
        function store = drawOverlay(obj, ax, layers, settings, store)
            %DRAWOVERLAY  WFS-specific overlays: reference contour, amplitude
            %  distribution, tapering window, illuminated/tapered zones.
            rp = obj.get_renderer_visuals();
            store = audioapp.ui.draw.WFSPainter.draw(ax, rp, obj, settings, store);
        end

        function settings = defaultOverlaySettings(~)
            settings = struct('showRef',true,'showAmp',true,'showTaper',true,'showZones',true);
        end
    end

end

%nti xl3
