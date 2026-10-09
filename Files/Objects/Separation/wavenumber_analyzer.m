classdef wavenumber_analyzer < handle
    %WAVENUMBER_ANALYZER  Local wavenumber vector and phase Hessian of the
    %sound field that a channel mix creates at the listening position
    %(ray-theory description of the local wavefront, after Firtha,
    %"Limitations of WFS" Part I/II), per time window (block switching
    %as in soundfield_analyzer) and per frequency bin or critical band.
    %
    %   Reproduction model: channel i is a point source (3D Green's
    %   function G = e^{-ikr}/(4 pi r)) at x_i = R_i u_i around the listener
    %   at the origin. At the origin, per bin (closed form):
    %     P = sum X_i G_i,   grad P = sum X_i G_i (ik + 1/R_i) u_i,
    %     Hess P = sum X_i G_i { [(ik + 1/R_i)^2 + 1/R_i^2] u_i u_i' -
    %                            (ik + 1/R_i)(I - u_i u_i')/R_i }
    %   With P = A e^{ik phi} and g = grad P / P:
    %     grad A / A = Re g,   k grad phi = Im g
    %     local wavenumber   k^ = -grad phi = -Im(g)/k   (|k^| = 1 for a
    %                        simple field; |k^| < 1 or > 1: interference)
    %     Hess phi = Im(Hess P / P - g g')/k
    %     wavefront curvature across the ray (horizontal):
    %                        kappa_h = -t' (Hess phi) t,  t = k^ turned 90 deg
    %     vertical:          kappa_v = -(Hess phi)_zz
    %     transport check:   kappa_tr = -k^ . (grad A / A) / |k^|^2
    %     source position (ray traced back): x_s = -rho_h k^/|k^|,
    %                        rho_h = 1/kappa_h (< 0: converging, a focus
    %                        between the loudspeakers and the listener)
    %   One loudspeaker gives k^ = -u_i, |k^| = 1, kappa_h = kappa_v =
    %   1/R_i: the source at x_i.
    %
    %   Aggregation: 'bands' (critical bands, Zwicker) - power (|P|^2)
    %   weighted means of the per-bin k^ and Hess phi, plus the band's
    %   direction spread (1 - |mean unit k^|, 0 = one direction); 'bins' -
    %   every bin (above MinHz) kept.
    %
    %   A = analyze(x, fs): per frame t, dur, isShort; per band (or bin)
    %   and frame: kx, ky (k^), kmag, kappaH, kappaV, kappaTr, E (|P|^2),
    %   spread; band edges; block switching flags.

    properties
        R = 2               % loudspeaker distance(s) [m] (scalar or per channel)
        Layout = []         % [nCh x 3] unit vectors ([] = from channel_layout)
        Mode = 'bands'      % 'bands' (critical bands) | 'bins'
        BandEdges = [100 200 300 400 510 630 770 920 1080 1270 1480 1720 2000 2320 ...
                     2700 3150 3700 4400 5300 6400 7700 9500 12000 15500]   % Zwicker
        MinHz = 100
        Plan                % soundfield_analyzer: block switching settings
        c = 343
    end

    methods
        function obj = wavenumber_analyzer(varargin)
            obj.Plan = soundfield_analyzer();
            for k = 1:2:numel(varargin), obj.(varargin{k}) = varargin{k+1}; end
        end

        function A = analyze(obj, x, fs)
            nCh = size(x, 2);
            U = obj.Layout;
            if isempty(U)
                az = separator_base.channel_layout(nCh);
                keep = ~isnan(az);
                x = x(:, keep); az = az(keep);
                U = separator_base.dir2pos(az);
            end
            Rv = obj.R(:)'.*ones(1, size(U,1));
            pl = obj.Plan;
            [fr, gAtt, gSpat, gShort] = pl.block_plan(x, U, fs);
            NL = pl.NLong; NS = 2*pl.Hop/pl.NSub; H = pl.Hop;
            nF = size(fr, 1);
            xp = [zeros(NL, size(x,2)); x; zeros(NL, size(x,2))];
            win = @(N) 0.5*(1 - cos(2*pi*((0:N-1)' + 0.5)/N));
            be = obj.BandEdges(:)';
            lf = sqrt(be(1:end-1).*be(2:end))' < pl.XoverHz;   % bands kept at long resolution (short frames)
            % per-window-length precomputation
            S = struct();
            for N = unique(fr(:,2))'
                S.(sprintf('n%d', N)) = wavenumber_analyzer.precompute(N, fs, U, Rv, win(N), be, obj.MinHz, obj.c, obj.Mode);
            end
            if strcmpi(obj.Mode, 'bands')
                nQ = numel(be) - 1;
            else
                nQ = size(S.(sprintf('n%d', NL)).Bin, 1);    % bins of the long window
            end
            Z = nan(nQ, nF);
            A = struct('kx', Z, 'ky', Z, 'kmag', Z, 'kappaH', Z, 'kappaV', Z, 'kappaTr', Z, ...
                'E', zeros(nQ, nF), 'spread', Z);
            for f = 1:nF
                N = fr(f,2); c0 = round(fr(f,1));
                Q = obj.frame_stats(xp, NL + c0 - N/2 + (1:N), S.(sprintf('n%d', N)));
                if strcmpi(obj.Mode, 'bands') && N ~= NL
                    % short frame: low bands from the granule's long frame
                    g = floor(c0/H) + 1;
                    QL = obj.frame_stats(xp, NL + (g-1)*H + H/2 - NL/2 + (1:NL), S.(sprintf('n%d', NL)));
                    fn = fieldnames(Q);
                    for i = 1:numel(fn), Q.(fn{i})(lf) = QL.(fn{i})(lf); end
                end
                if strcmpi(obj.Mode, 'bins') && N ~= NL
                    continue;                                % bins: long frames only (one grid)
                end
                fn = fieldnames(Q);
                for i = 1:numel(fn), A.(fn{i})(:, f) = Q.(fn{i}); end
            end
            A.t = fr(:,1)'/fs; A.dur = fr(:,2)'/fs; A.isShort = (fr(:,2) == NS)';
            A.gAttack = gAtt; A.gSpatial = gSpat; A.gShort = gShort;
            A.bandEdges = be; A.fs = fs; A.U = U; A.R = Rv;
            A.az = atan2d(-A.kx, A.ky);                     % direction the wave travels to
            A.srcAz = mod(A.az, 360) - 180;                  % direction it comes from (+ = left)
            A.rhoH = 1./A.kappaH;
        end
    end

    methods (Static)
        function S = precompute(N, fs, U, Rv, w, be, minHz, c, mode)
            % Green's function factors per bin and loudspeaker for an
            % N-point window w: P, grad P, Hess P at the origin are sums of
            % X_i times these (closed form, point sources at R_i u_i)
            nb = N/2 + 1;
            f = (0:nb-1)'*fs/N;
            k = 2*pi*f/c;
            use = f >= minHz;
            k = k(use);
            G = exp(-1i*k*Rv)./(4*pi*Rv);                   % [bins x nCh]
            a = 1i*k + 1./Rv;                               % ik + 1/R
            S.w = w/sqrt(sum(w.^2));
            S.use = use; S.k = k; S.G = G; S.f = f(use);
            S.Gx = G.*a.*U(:,1)'; S.Gy = G.*a.*U(:,2)'; S.Gz = G.*a.*U(:,3)';
            % Hessian: G {[a^2 + 1/R^2] u u' - a (I - u u')/R}
            b1 = G.*(a.^2 + 1./Rv.^2); b2 = G.*a./Rv;
            uu = @(p, q) U(:,p)'.*U(:,q)';
            S.Hxx = b1.*uu(1,1) - b2.*(1 - uu(1,1));
            S.Hyy = b1.*uu(2,2) - b2.*(1 - uu(2,2));
            S.Hxy = b1.*uu(1,2) - b2.*(-uu(1,2));
            S.Hzz = b1.*uu(3,3) - b2.*(1 - uu(3,3));
            fu = f(use);
            if strcmpi(mode, 'bands')
                [~, ~, bi] = histcounts(fu, be);
                ok = bi > 0;
                S.Bin = sparse(bi(ok), find(ok), 1, numel(be) - 1, numel(fu));
            else
                S.Bin = speye(numel(fu));
            end
        end

        function Q = bin_stats(X, S)
            % Per-bin local wavefront of the spectrum X [bins (S.use) x nCh]:
            % kx, ky (k^), km = |k^|, kapH, kapV, kapTr, E = |P|^2, P
            P = sum(X.*S.G, 2);
            gx = sum(X.*S.Gx, 2)./P; gy = sum(X.*S.Gy, 2)./P;
            hxx = sum(X.*S.Hxx, 2)./P - gx.^2;
            hyy = sum(X.*S.Hyy, 2)./P - gy.^2;
            hxy = sum(X.*S.Hxy, 2)./P - gx.*gy;
            hzz = sum(X.*S.Hzz, 2)./P;                      % g_z = 0 (horizontal layout)
            k = S.k;
            Q.kx = -imag(gx)./k; Q.ky = -imag(gy)./k;       % local wavenumber vector
            Pxx = imag(hxx)./k; Pyy = imag(hyy)./k; Pxy = imag(hxy)./k;   % Hess phi
            Q.km = sqrt(Q.kx.^2 + Q.ky.^2);
            tx = -Q.ky./max(Q.km, eps); ty = Q.kx./max(Q.km, eps);   % transverse unit vector
            Q.kapH = -(tx.^2.*Pxx + 2*tx.*ty.*Pxy + ty.^2.*Pyy);
            Q.kapV = -imag(hzz)./k;
            Q.kapTr = -(Q.kx.*real(gx) + Q.ky.*real(gy))./max(Q.km.^2, eps);
            Q.E = abs(P).^2; Q.P = P;
            bad = ~isfinite(Q.kx) | ~isfinite(Q.kapH);
            if any(bad)
                fn = {'kx','ky','km','kapH','kapV','kapTr'};
                for i = 1:numel(fn), Q.(fn{i})(bad) = 0; end
                Q.E(bad) = 0;
            end
        end
    end

    methods (Access = private)
        function Q = frame_stats(~, xp, idx, S)
            X = fft(xp(idx, :) .* S.w, [], 1);
            X = X(1:numel(S.use), :); X = X(S.use, :);
            q = wavenumber_analyzer.bin_stats(X, S);
            kx = q.kx; ky = q.ky; km = q.km; kapH = q.kapH; kapV = q.kapV; kapTr = q.kapTr; E = q.E;
            % aggregation (power weighted)
            B = S.Bin; sw = B*E; sw(sw == 0) = eps;
            m = @(q) (B*(E.*q))./sw;
            Q.kx = m(kx); Q.ky = m(ky);
            Q.kmag = m(km);
            Q.kappaH = m(kapH); Q.kappaV = m(kapV); Q.kappaTr = m(kapTr);
            Q.E = B*E;
            ux = kx./max(km, eps); uy = ky./max(km, eps);
            Q.spread = 1 - sqrt(m(ux).^2 + m(uy).^2);
        end
    end
end
