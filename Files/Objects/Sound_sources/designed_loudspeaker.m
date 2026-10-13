classdef designed_loudspeaker < handle
    %DESIGNED_LOUDSPEAKER  A loudspeaker of the Loudspeaker designer as a
    %source model of the Binauralizer (ls_binaural_export).
    %
    %   One channel per driver: channel k is driver k driven alone with 1 V
    %   (the other drivers shorted through the amplifier, i.e. the box
    %   coupling and the vents it drives are included).
    %
    %   Local frame: the loudspeaker faces +x, +y is its left, z up; the
    %   origin is the centre of the model's bounding box. In the scene the
    %   model is placed at a position and turned by a yaw angle about z
    %   (speaker_placement).
    %
    %   Data:
    %     Mesh     V [Nv x 3] (local), F [Nf x 3], C [Nf x 3] face colours
    %              with the shading baked in (as scene_glyph)
    %     Drivers  name, label, pos [1 x 3] and axis [1 x 3] (local)
    %     Az, El   direction grid [deg] (azimuth from +x towards +y,
    %              elevation from the horizontal plane)
    %     f        frequencies [Hz]
    %     H        [nAz x nEl x nCh x nf] response of channel k towards the
    %              direction (on the export sphere), propagation from the
    %              driver removed: p(r) = Gref H e^{-jkr}/(4 pi r) per volt
    %              (r from the driver), i.e. H = 1 is a point source of
    %              strength Gref. Used where the SH model does not apply.
    %     SH       the field as outgoing spherical waves about the origin,
    %              exact at any distance r >= a (near field included):
    %                p(r, dir)/Gref = sum_nm B_nm h_n(kr)/h_n(k r0) Y_nm(dir)
    %              (h_n = h_n^(2), e^{jwt}; Y_nm = getSpherHarmMx 'complex')
    %              fields: B [(Nmax+1)^2 x nCh x nf] (SHT of p/Gref on the
    %              sphere r0), N [nf x 1] order per frequency (0: not
    %              stored, above fTop), Nmax, r0, a (radius of the sphere
    %              enclosing the model), fTop, residual [nf x nCh] (relative
    %              SHT fit error)
    %     Gref     [Pa m / V * 4 pi] the normalisation (on-axis level of the
    %              loudest channel, 300 Hz - 3 kHz); Sensitivity() in dB
    %
    %   G = s.response(k, x, fq): channel k at the listener position x
    %   (local frame, from the origin), frequencies fq, normalised as H;
    %   h = s.fir(k, x, fs, N): the same as an FIR filter (bulk delay
    %   FirDelay samples); p = s.field(k, X, f0): normalised pressure p/Gref
    %   at the points X [M x 3] (local), one frequency.
    %   SH below fTop (blended into the H model over the top third octave),
    %   H above it and for points closer than a.
    %   s.save(file) adds the variable 'speaker' to file (a design file of
    %   the Loudspeaker designer); designed_loudspeaker.load(file).

    properties
        Name = 'loudspeaker'
        Mesh = struct('V', zeros(0,3), 'F', zeros(0,3), 'C', zeros(0,3))
        Drivers = struct('name', {}, 'label', {}, 'pos', {}, 'axis', {})
        Az = []
        El = []
        f = []
        H = []
        SH = struct([])
        Gref = 1
        Info = struct()
        File = ''
        PreEQ = []              % its Pre-EQ (speaker_eq data: biquads / FIR), [] = none: applied by the
                                % Binauralizer on its driving signals (Speaker EQ, switchable), not in H
    end

    properties (Constant)
        FirDelay = 48           % [samples] bulk delay of fir() (keeps the low-frequency interpolation causal)
        c0 = 343                % [m/s] as in ls_binaural_export
    end

    methods
        function obj = designed_loudspeaker(varargin)
            for k = 1:2:numel(varargin), obj.(varargin{k}) = varargin{k+1}; end
        end

        function n = nChannels(obj), n = numel(obj.Drivers); end

        function s = Sensitivity(obj)
            % on-axis level of the loudest channel [dB SPL @ 1 V, 1 m]
            s = 20*log10(obj.Gref/(4*pi)/20e-6);
        end

        function tf = hasSH(obj), tf = ~isempty(obj.SH) && any(obj.SH.N > 0); end

        function Hm = directions(obj, k, D)
            % [M x nf] channel k (H model) towards the local directions D
            % [M x 3]; bilinear on the az/el grid, az periodic
            D = D./max(vecnorm(D, 2, 2), eps);
            az = mod(atan2d(D(:,2), D(:,1)), 360);
            el = asind(max(-1, min(1, D(:,3))));
            A = obj.Az(:); E = obj.El(:); nA = numel(A); nE = numel(E);
            u = (az - A(1))/(A(2) - A(1)); i0 = floor(u); wa = u - i0;
            i1 = mod(i0, nA) + 1; i2 = mod(i0 + 1, nA) + 1;
            j0 = min(max(sum(E' <= el, 2), 1), nE - 1);
            we = min(max((el - E(j0))./(E(j0+1) - E(j0)), 0), 1);
            Hk = reshape(obj.H(:, :, k, :), nA*nE, []);
            ix = @(i, j) i + (j - 1)*nA;
            Hm = (1-wa).*(1-we).*Hk(ix(i1, j0), :) + wa.*(1-we).*Hk(ix(i2, j0), :) + ...
                 (1-wa).*we.*Hk(ix(i1, j0+1), :) + wa.*we.*Hk(ix(i2, j0+1), :);
        end

        function Hd = direction(obj, k, d)
            % [nf x 1] channel k (H model) towards one local direction d
            Hd = obj.directions(k, d(:)').';
        end

        function Gs = stored(obj, k, x)
            % [nf x 1] channel k at the local position x (from the origin),
            % stored frequencies, normalised as H (propagation from the
            % driver removed): SH where it applies, else the H model
            x = x(:)';
            dk = x - obj.Drivers(k).pos;
            Gs = obj.direction(k, dk);
            r = norm(x);
            if ~obj.hasSH() || r < obj.SH.a, return; end
            % all frequencies at once: the angular part once, the radial
            % ratios h_n(kr)/h_n(k r0) by the upward recursion (radialAll),
            % one product - a few ms per point instead of a besselh call per
            % order and frequency
            S = obj.SH; fv = obj.f(:);
            w = obj.shWeight(fv); ms = find(w > 0);
            if isempty(ms), return; end
            Y = designed_loudspeaker.shRow(acos(x(3)/r), atan2(x(2), x(1)), S.Nmax);
            rk = norm(dk); kv = 2*pi*fv(ms).'/obj.c0;
            L = (S.Nmax + 1)^2; nIdx = floor(sqrt(0:L - 1))' + 1;
            HR = designed_loudspeaker.radialAll(S.Nmax, kv, r, S.r0);
            Bk = reshape(S.B(1:L, k, ms), L, numel(ms));
            Nm = S.N(ms); mask = (nIdx - 1) <= Nm(:).';
            p = Y(1:L)*(Bk.*HR(nIdx, :).*mask);
            Gsh = (p.*4*pi*rk.*exp(1i*kv*rk)).';
            Gs(ms) = w(ms).*Gsh + (1 - w(ms)).*Gs(ms);
        end

        function p = field(obj, k, X, f0)
            % normalised pressure p/Gref of channel k at local points X
            % [M x 3], frequency f0 (point-source scale: G e^{-jkr}/(4 pi r))
            M = size(X, 1); kk = 2*pi*f0/obj.c0; fv = obj.f(:);
            % H model: per point, from the driver, interpolated in f
            Dk = X - obj.Drivers(k).pos; rk = vecnorm(Dk, 2, 2);
            [m1, m2, t] = bracket(fv, f0);
            Hm = obj.directions(k, Dk);
            p = ((1-t)*Hm(:, m1) + t*Hm(:, m2)).*exp(-1i*kk*rk)./(4*pi*rk);
            if ~obj.hasSH(), return; end
            w = obj.shWeight(f0);
            r = vecnorm(X, 2, 2); out = r >= obj.SH.a;
            if w == 0 || ~any(out), return; end
            % SH: exact at the neighbouring stored frequencies; between
            % them the response with the propagation from the driver
            % removed is interpolated (smooth in f), then put back at f0
            ok = find(obj.SH.N > 0); [q1, q2, tq] = bracket(fv(ok), f0); q1 = ok(q1); q2 = ok(q2);
            idx = find(out); rj = rk(idx);
            g = 0;
            for q = [q1 q2; 1 - tq, tq]
                if q(2) == 0, continue; end
                kq = 2*pi*fv(q(1))/obj.c0;
                g = g + q(2)*obj.shEval(k, X(idx, :), q(1)).*(4*pi*rj).*exp(1i*kq*rj);
            end
            p(idx) = w*g.*exp(-1i*kk*rj)./(4*pi*rj) + (1 - w)*p(idx);
            p = reshape(p, M, 1);
        end

        function G = response(obj, k, x, fq)
            % channel k at the local listener position x, frequencies fq:
            % magnitude (log) and unwrapped phase interpolated in log f;
            % below the first frequency the magnitude falls 12 dB/octave,
            % above the last it is held
            G = obj.interpF(obj.stored(k, x), fq);
        end

        function G = interpF(obj, Hd, fq)
            % stored-frequency responses Hd [nf x M] at the frequencies fq
            % (as response(): log magnitude pchip, unwrapped phase linear in
            % log f, 12 dB/octave below the first, held above the last)
            fs0 = obj.f(:); fq = fq(:);
            mag = log(abs(Hd) + 1e-12); ph = unwrap(angle(Hd), [], 1);
            lf = min(max(log(max(fq, eps)), log(fs0(1))), log(fs0(end)));
            M = interp1(log(fs0), mag, lf, 'pchip');
            P = interp1(log(fs0), ph, lf, 'linear');
            lo = fq < fs0(1);
            if any(lo), M(lo, :) = mag(1, :) + 2*log(max(fq(lo), eps)/fs0(1)); end
            G = exp(M + 1i*P);
            G(fq <= 0, :) = 0;
        end

        function h = fir(obj, k, x, fs, N)
            % N-tap FIR of channel k at the local listener position x (fs
            % [Hz]); delayed by FirDelay samples, the last quarter faded out
            fk = (0:N/2)'*fs/N;
            X = obj.response(k, x, fk).*exp(-1i*2*pi*fk/fs*obj.FirDelay);
            X(end) = real(X(end));
            h = real(ifft([X; conj(X(end-1:-1:2))]));
            nF = round(N/4);
            h(end-nF+1:end) = h(end-nF+1:end).*(0.5 + 0.5*cos(pi*(1:nF)'/nF));
        end

        function s = toStruct(obj)
            s = struct('Name', obj.Name, 'Mesh', obj.Mesh, 'Drivers', obj.Drivers, 'Az', obj.Az, ...
                'El', obj.El, 'f', obj.f, 'H', obj.H, 'SH', obj.SH, 'Gref', obj.Gref, 'Info', obj.Info, 'PreEQ', obj.PreEQ);
        end

        function save(obj, file)
            % adds / replaces the variable 'speaker' in file
            speaker = obj.toStruct(); %#ok<NASGU>
            if isfile(file), save(file, 'speaker', '-append'); else, save(file, 'speaker'); end
            obj.File = file;
        end
    end

    methods (Access = private)
        function w = shWeight(obj, fq)
            % 1 in the SH band, falling to 0 over its top third octave
            S = obj.SH; w = zeros(size(fq));
            if isempty(S) || ~(S.fTop > 0), return; end
            w = min(max(log(S.fTop./fq)/log(2^(1/3)), 0), 1);
            if numel(fq) == numel(obj.f)               % stored frequencies: only where B exists
                w(S.N(:) == 0) = 0;
            end
        end

        function p = shEval(obj, k, X, m)
            % SH series of channel k at stored frequency m, points X [M x 3]
            % (all at r >= a); chunks keep Y [chunk x (N+1)^2] small
            S = obj.SH; N = S.N(m); L = (N+1)^2; kk = 2*pi*obj.f(m)/obj.c0;
            [~, nIdx] = obj.radial(N, kk, 1, 1);
            h0 = getSphH((0:N)', 2, kk*S.r0);
            c = S.B(1:L, k, m)./h0(nIdx); c(~isfinite(c)) = 0;
            r = vecnorm(X, 2, 2); p = zeros(size(X, 1), 1);
            for s = 1:4000:size(X, 1)
                j = s:min(s + 3999, size(X, 1));
                Y = getSpherHarmMx(acos(X(j,3)./r(j)), atan2(X(j,2), X(j,1)), N, 'complex');
                Hn = getSphH((0:N)', 2, kk*r(j)');     % [N+1 x chunk]
                p(j) = sum(Y.*(c.*Hn(nIdx, :)).', 2);
            end
        end
    end

    methods (Static)
        function obj = fromStruct(s)
            obj = designed_loudspeaker();
            for fn = fieldnames(s)'
                if isprop(obj, fn{1}), obj.(fn{1}) = s.(fn{1}); end
            end
        end

        function obj = load(file)
            if ~designed_loudspeaker.isSpeakerFile(file)
                error('designed_loudspeaker:noSpeaker', ...
                    '%s has no exported loudspeaker (Loudspeaker designer: Send to Binauralizer).', file);
            end
            L = load(file, 'speaker');
            obj = designed_loudspeaker.fromStruct(L.speaker);
            obj.File = file;
        end

        function tf = isSpeakerFile(file)
            tf = false;
            try
                w = whos('-file', file);
                tf = any(strcmp({w.name}, 'speaker'));
            catch
            end
        end

        function S = fitSH(P, D, f, a, r0, Nmax)
            % SH model from the normalised pressure P [K x nCh x nf] on the
            % sphere r0 at the local directions D [K x 3] (least squares,
            % getSpherHarmMx). Order per frequency N = ceil(k a) + 4; where
            % that exceeds Nmax the frequency is left to the H model. r0
            % should be just outside a: the inward extrapolation (a <= r <
            % r0) amplifies order n by up to (r0/a)^(n+1).
            f = f(:); nf = numel(f); nCh = size(P, 2);
            Y = getSpherHarmMx(acos(max(-1, min(1, D(:,3)))), atan2(D(:,2), D(:,1)), Nmax, 'complex');
            N = ceil(2*pi*f/designed_loudspeaker.c0*a) + 4; N(N > Nmax) = 0;
            B = zeros((Nmax+1)^2, nCh, nf); res = nan(nf, nCh);
            for m = find(N > 0)'
                L = (N(m)+1)^2;
                B(1:L, :, m) = Y(:, 1:L)\P(:, :, m);
                res(m, :) = vecnorm(Y(:, 1:L)*B(1:L, :, m) - P(:, :, m))./vecnorm(P(:, :, m));
            end
            fTop = max([0; f(N > 0)]);
            S = struct('B', B, 'N', N, 'Nmax', Nmax, 'r0', r0, 'a', a, 'fTop', fTop, 'residual', res);
        end

        function Y = shRow(th, ph, N)
            % [1 x (N+1)^2] complex spherical harmonics at one direction, the
            % convention and column order of getSpherHarmMx ('complex':
            % (-1)^m sqrt((2n+1)/(4 pi) (n-|m|)!/(n+|m|)!) P_n^|m| e^{i m ph},
            % P with the Condon-Shortley phase) by the normalised Legendre
            % recursion - no factorials, ~50x faster
            x = cos(th); sn = sin(th);
            P = zeros(N + 1, N + 1);                       % P(n+1, m+1): normalised, CS phase
            P(1, 1) = 1/sqrt(4*pi);
            for m = 1:N, P(m+1, m+1) = -sqrt((2*m + 1)/(2*m))*sn*P(m, m); end
            for m = 0:N - 1, P(m+2, m+1) = sqrt(2*m + 3)*x*P(m+1, m+1); end
            for m = 0:N
                for n = m + 2:N
                    a = sqrt((4*n^2 - 1)/(n^2 - m^2)); b = sqrt(((n - 1)^2 - m^2)/(4*(n - 1)^2 - 1));
                    P(n+1, m+1) = a*(x*P(n, m+1) - b*P(n-1, m+1));
                end
            end
            Y = zeros(1, (N + 1)^2);
            for n = 0:N
                m = -n:n;
                Y(n^2 + n + m + 1) = (-1).^m.*P(n+1, abs(m) + 1).*exp(1i*ph*m);
            end
        end

        function HR = radialAll(N, kv, r, r0)
            % [N+1 x nk] h_n(k r)/h_n(k r0) (spherical Hankel, 2nd kind),
            % n = 0..N, all wavenumbers kv at once: h_0 = i e^{-ix}/x,
            % h_1 = e^{-ix}(i/x^2 - 1/x), h_{n+1} = (2n+1)/x h_n - h_{n-1}
            % (stable upward for the Hankel functions)
            HR = designed_loudspeaker.sphH2(N, kv(:).'*r)./designed_loudspeaker.sphH2(N, kv(:).'*r0);
            HR(~isfinite(HR)) = 0;
        end

        function H = sphH2(N, x)
            % [N+1 x numel(x)] spherical Hankel functions h_n^(2)(x), n = 0..N
            x = x(:).'; H = zeros(N + 1, numel(x)); e = exp(-1i*x);
            H(1, :) = 1i*e./x;
            if N >= 1, H(2, :) = e.*(1i./x.^2 - 1./x); end
            for n = 1:N - 1, H(n + 2, :) = (2*n + 1)./x.*H(n + 1, :) - H(n, :); end
        end

        function [hr, nIdx] = radial(N, k, r, r0)
            % h_n(kr)/h_n(k r0), n = 0..N, and the order of each SH column
            hr = getSphH((0:N)', 2, k*r)./getSphH((0:N)', 2, k*r0);
            hr(~isfinite(hr)) = 0;
            nIdx = floor(sqrt(0:(N+1)^2 - 1))' + 1;
        end
    end
end

function [m1, m2, t] = bracket(fv, f0)
% neighbouring indices of f0 in fv and the log-f weight of m2 (held outside)
if f0 <= fv(1), m1 = 1; m2 = 1; t = 0; return; end
if f0 >= fv(end), m1 = numel(fv); m2 = m1; t = 0; return; end
m1 = find(fv <= f0, 1, 'last'); m2 = m1 + 1;
t = log(f0/fv(m1))/log(fv(m2)/fv(m1));
end
