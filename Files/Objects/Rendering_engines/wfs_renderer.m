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

        function obj = update_renderer(obj, ~)
            N = length(obj.virtual_source.source_signal.time_series);
            omega = (0 : N/2)'/N*2*pi*obj.virtual_source.source_signal.fs;

            [obj.amp,obj.delay] = obj.get_delay_and_gain;
            if obj.setup.Antialiasing
                obj.get_antialiasing_filters(obj.delay);
                D_wfs = obj.antialiasing_filter_bank.* (obj.H_pre*obj.amp.').*exp(-1i*omega*obj.delay');length(obj.virtual_source.source_signal.time_series);
                d_wfs =  fftshift( ifft( D_wfs,length(obj.virtual_source.source_signal.time_series),1,'symmetric'),1 );
            else
                D_wfs = (obj.H_pre*obj.amp.').*exp(-1i*omega*obj.delay');length(obj.virtual_source.source_signal.time_series);
                d_wfs = fftshift(ifft( D_wfs,length(obj.virtual_source.source_signal.time_series),1,'symmetric'),1);
            end

            for n = 1 : length(obj.WFS_filter_bank)
                obj.WFS_filter_bank{n}.update_coefficients(d_wfs(:,n));
            end
        end

        function obj = update_settings(obj,setup)
            obj.setup = setup;
            obj.update_renderer;
        end

        function render(obj)
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
            center = [0, 0];
            x0 = obj.ssd.x0;
            n0 = obj.ssd.n0;
            dl = obj.ssd.dA;
            c =  audioapp.util.Physics.speedOfSound();
            switch obj.virtual_source.source_type.Shape
                case 'plane_wave'
                    kP = bsxfun( @minus, [0,0], obj.virtual_source.position );
                    obj.kh = repmat( kP/norm(kP) , size(obj.ssd.x0,1),1);
                    obj.khn = sum(obj.kh.*n0,2);
                    win = double(obj.khn>=0);

                    obj.get_reference_distance;
                    dcorr = sqrt(obj.dref);
                    delay =  sum( (x0.*obj.kh) / c , 2);
                    amp = 1/4/pi*sqrt(8*pi/c).*win.*dcorr.*obj.khn.*dl;
                    ks = 1;
                    obj.F = 1;
                case 'point_source'
                    kP = bsxfun( @minus, x0, obj.virtual_source.position );
                    rhoP = sqrt(sum(kP.^2,2));
                    obj.kh = kP./sqrt(sum(kP.^2,2));
                    obj.khn = sum(obj.kh.*n0,2);

                    ks = (obj.virtual_source.position -center)/norm(obj.virtual_source.position-center);
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

                    dcorr = sqrt(obj.dref.*rhoP./(obj.F*obj.dref+rhoP));
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
                    obj.dref = sqrt(sum( (obj.ssd.x0-obj.setup.RefPoint).^2,2 ));

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
                integrand = abs(obj.amp/4/pi./sqrt(sum( (receiver.position-obj.ssd.x0).^2,2 )  ));
                switch obj.ssd.isClosed
                    case true
                        [~,ix_stat] = min( sum((obj.ssd.x0-intersection(1,:)).^2,2) );
                        integrand = circshift(integrand.*obj.ssd.dA, -(ix_stat-1))/integrand(ix_stat);
                        I = [sum(integrand(1:round(end/2))), sum(integrand(round(end/2)+1:end))];
                    case false
                        [~,ix_stat] = min( sum((obj.ssd.x0-intersection(1,:)).^2,2) );
                        I = [sum(integrand(1:ix_stat).*obj.ssd.dA(1:ix_stat)), sum(integrand(ix_stat+1:end).*obj.ssd.dA(ix_stat+1:end))]/integrand(ix_stat);
                end
                rhoP = norm(obj.virtual_source.position-obj.ssd.x0(ix_stat,:));
                rhoG = norm(receiver.position-obj.ssd.x0(ix_stat,:));
                switch lower(obj.virtual_source.source_type.Shape)
                    case 'plane_wave'
                        d_rec = rhoG;
                    case 'point_source'
                        d_rec = rhoP*rhoG/(rhoP+obj.F*rhoG);
                end
                kc = pi/2*d_rec./(I.^2*obj.khn(ix_stat).^2);
                fc = abs( kc*audioapp.util.Physics.speedOfSound()/2/pi ) ;
            else
                fc = inf;
            end
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
            % Return N x 2 unit vectors kh at each SSD point (same logic you use in get_delay_and_gain)
            X = obj.ssd.x0;
            switch lower(obj.virtual_source.source_type.Shape)
                case 'plane_wave'
                    % direction taken from your code
                    kP = [0,0] - obj.virtual_source.position(:).';
                    dir = kP / max(norm(kP),1e-12);
                    kh = repmat(dir, size(X,1), 1);
                otherwise  % 'point_source' (and default)
                    kP = bsxfun(@minus, X, obj.virtual_source.position);
                    n  = sqrt(sum(kP.^2,2));
                    n(n<1e-12) = 1e-12;
                    kh = kP ./ n;        % N x 2
            end
        end

    end



end

%nti xl3