classdef stft_separator < separator_base
    %STFT_SEPARATOR  Real-time primary/ambient + directional decomposition
    %of stereo (PCA-based, Goodwin-Jot / He-Tan-Gan), output as
    %fixed-direction objects.
    %
    %   Per STFT bin, from the recursively smoothed 2x2 covariance
    %   (time constant Tau) with eigenvalues l1 >= l2 and principal unit
    %   eigenvector u (the panning direction of the dominant source):
    %     primary    P = sqrt(g) (u' x),  g = (l1 - l2)/l1  (ambience is
    %                assumed uncorrelated with equal power l2 per channel)
    %     ambience   A = x - u (1 - sqrt(1-g)) (u' x)   (power complementary)
    %     direction  stereo tangent law from |u1|/|u2| (source base +-Base),
    %                stretched to the output width +-Width
    %   A hard-panned source has u = [1 0] (primary, not ambience), a
    %   diffuse field l1 = l2 (all ambience).
    %   P is panned (constant power) between the two neighbouring of K
    %   front "slots" spread over +-Width; the ambience pair goes to
    %   +-AmbAz. Objects therefore sit still (no renderer updates per
    %   block) while the direction per bin stays continuous.
    %   Mode 'objects': instead of fixed slots, up to K moving objects
    %   that follow the sources: the primary energy's direction histogram
    %   (smoothed in time and azimuth) is peak-picked every frame, the peaks
    %   are tracked (nearest track within Gate, smoothed azimuth, activity
    %   with attack/release), and each bin's primary part goes to the
    %   active tracks near its direction (Gaussian, power complementary).
    %   Pos follows the tracks; TrackLog records them (azimuth per frame,
    %   NaN = inactive).
    %   Informed demixing (Demix): where a bin lies between two active
    %   sources (second-best mask weight >= DemixRatio x best) both are
    %   assumed present; with their known panning vectors a(theta) the 2x2
    %   mixing is inverted (regularised least squares), which separates
    %   two overlapping panned sources that the (orthogonal) PCA cannot.
    %   5.1/7.1 input: L/R are decomposed, C is added to the centre slot,
    %   the surround channels join the ambience objects, LFE is dropped.
    %
    %   Analysis: sqrt-Hann, N = 2048, hop N/2 (perfect reconstruction);
    %   latency N - hop samples.

    properties
        K = 5               % front slots
        Width = 30          % output half-width of the front stage [deg]
        Base = 30           % stereo source base of the content [deg]
        AmbAz = 110         % ambience object azimuth [deg]
        AmbGain = 1
        Tau = 0.05          % statistics smoothing [s]
        N = 2048            % FFT size
        Diagnostics = false % keep per-frame direction histograms
        Mode = 'slots'      % 'slots' (fixed directions) | 'objects' (tracked)
        TrackSigma = 8      % [deg] width of an object's direction mask
        TrackTau = 0.25     % [s] direction-histogram smoothing
        PeakFloorDb = -20   % peaks below max - this are ignored
        MinSep = 12         % [deg] minimum separation of two sources
        Gate = 20           % [deg] a peak continues a track within this
        HistBand = [150 8000] % [Hz] bins used to find directions (LF/HF directions are unreliable)
        Demix = false       % informed 2x2 demixing of overlapping sources (experimental:
                            % phantom tracks between two real sources make it worse)
        DemixRatio = 0.1    % second / first mask weight that triggers it
        DemixReg = 0.05     % Tikhonov regularisation (relative)
    end

    properties (SetAccess = private)
        DiagAz = []         % [nAzBins x frames] primary energy per direction
        DiagPA = []         % [2 x frames] primary / ambient energy
        AzEdges = -90:2:90
        TrackLog = []       % [K x frames] tracked azimuths (NaN = inactive)
    end

    properties (Access = private)
        hop; win; frame; acc; pending; outq
        Sll; Srr; Slr
        slotAz; centre; layoutIdx
        trAz; trAct; hist; kern
    end

    methods
        function obj = stft_separator(varargin)
            for k = 1:2:numel(varargin), obj.(varargin{k}) = varargin{k+1}; end
        end

        function init(obj, fs, nIn)
            obj.fs = fs; obj.nIn = nIn;
            obj.K = max(2, round(obj.K));
            obj.slotAz = linspace(obj.Width, -obj.Width, obj.K);      % left -> right
            [~, obj.centre] = min(abs(obj.slotAz));
            if strcmpi(obj.Mode, 'objects')
                obj.Names = [arrayfun(@(k) sprintf('Obj %d', k), 1:obj.K, 'UniformOutput', false), ...
                    {'Amb L', 'Amb R'}];
            else
                obj.Names = [arrayfun(@(a) sprintf('Front %+.0f', a), obj.slotAz, 'UniformOutput', false), ...
                    {'Amb L', 'Amb R'}];
            end
            g = -6:6;                                   % azimuth smoothing kernel (2 deg bins)
            obj.kern = exp(-(2*g).^2/(2*3^2))'; obj.kern = obj.kern/sum(obj.kern);
            obj.Pos = separator_base.dir2pos([obj.slotAz(:); obj.AmbAz; -obj.AmbAz]);
            obj.hop = obj.N/2;
            obj.win = sqrt(0.5*(1 - cos(2*pi*(0:obj.N-1)'/obj.N)));   % periodic sqrt-Hann
            obj.Latency = obj.N - obj.hop;
            % input channel roles
            [az, names] = separator_base.channel_layout(nIn);
            obj.layoutIdx = struct('L', find(strcmp(names,'L')), 'R', find(strcmp(names,'R')), ...
                'C', find(strcmp(names,'C')), 'SL', find(ismember(names, {'Ls','Lb'})), ...
                'SR', find(ismember(names, {'Rs','Rb'})));
            if nIn == 1, obj.layoutIdx.L = 1; obj.layoutIdx.R = 1; end
            if nIn > 2 && isempty(obj.layoutIdx.L)               % unknown layout: first two
                obj.layoutIdx.L = 1; obj.layoutIdx.R = 2;
            end
            obj.reset();
            %#ok<*NASGU>
        end

        function reset(obj)
            nb = obj.N/2 + 1;
            obj.frame = zeros(obj.N, obj.nIn);
            obj.acc = zeros(obj.N, obj.nOut());
            obj.pending = zeros(0, obj.nIn);
            obj.outq = zeros(0, obj.nOut());
            obj.Sll = 1e-9*ones(nb,1); obj.Srr = 1e-9*ones(nb,1); obj.Slr = zeros(nb,1);
            obj.DiagAz = []; obj.DiagPA = [];
            obj.trAz = obj.slotAz; obj.trAct = zeros(1, obj.K);
            obj.hist = zeros(numel(obj.AzEdges)-1, 1);
            obj.TrackLog = [];
            obj.Pos(1:obj.K,:) = separator_base.dir2pos(obj.slotAz(:));
        end

        function Y = process(obj, X)
            n = size(X,1);
            obj.pending = [obj.pending; X];
            while size(obj.pending,1) >= obj.hop
                obj.frame = [obj.frame(obj.hop+1:end,:); obj.pending(1:obj.hop,:)];
                obj.pending = obj.pending(obj.hop+1:end,:);
                obj.acc = obj.acc + obj.synth_frame();
                obj.outq = [obj.outq; obj.acc(1:obj.hop,:)];
                obj.acc = [obj.acc(obj.hop+1:end,:); zeros(obj.hop, obj.nOut())];
            end
            m = min(n, size(obj.outq,1));
            Y = [zeros(n-m, obj.nOut()); obj.outq(1:m,:)];   % start-up: zeros first
            obj.outq = obj.outq(m+1:end,:);
        end

        function Ys = separate_offline(obj, Xs)
            % Whole signals at once, aligned (no latency). Xs{1} drives all
            % decisions; every Xs{k} gets the same operations, so for
            % Xs = {mix, part1, part2} with mix = part1 + part2 the outputs
            % add up too (oracle evaluation of what goes where).
            obj.reset();
            h = obj.hop; N = obj.N; n = size(Xs{1},1); nO = obj.nOut();
            pad = @(x) [zeros(h, size(x,2)); x; zeros(N, size(x,2))];
            Xp = cellfun(pad, Xs, 'UniformOutput', false);
            Yp = cellfun(@(x) zeros(size(x,1), nO), Xp, 'UniformOutput', false);
            nb = N/2 + 1;
            for k = 0:floor((size(Xp{1},1) - N)/h)
                idx = k*h + (1:N);
                obj.frame = Xp{1}(idx,:);
                Fsh = cellfun(@(x) subsref(fft(x(idx,:).*obj.win, [], 1), struct('type','()','subs',{{1:nb, ':'}})), ...
                    Xp(2:end), 'UniformOutput', false);
                [Yf, Ysh] = obj.synth_frame(Fsh);
                Yp{1}(idx,:) = Yp{1}(idx,:) + Yf;
                for j = 1:numel(Ysh), Yp{j+1}(idx,:) = Yp{j+1}(idx,:) + Ysh{j}; end
            end
            Ys = cellfun(@(y) y(h + (1:n), :), Yp, 'UniformOutput', false);
        end
    end

    methods (Access = private)
        function [Yf, Ysh] = synth_frame(obj, Fsh)
            % One STFT frame. Decisions (statistics, directions, masks,
            % demixing) come from the input frame; apply_ops applies the
            % same linear operations to it and to optional shadow spectra
            % Fsh (cell) - e.g. a dialogue / M+E split for evaluation.
            I = obj.layoutIdx;
            F = fft(obj.frame .* obj.win, [], 1);
            nb = obj.N/2 + 1;
            F = F(1:nb,:);
            L = F(:, I.L); R = F(:, I.R);
            a = exp(-obj.hop/(obj.fs*obj.Tau));
            obj.Sll = a*obj.Sll + (1-a)*abs(L).^2;
            obj.Srr = a*obj.Srr + (1-a)*abs(R).^2;
            obj.Slr = a*obj.Slr + (1-a)*L.*conj(R);
            % eigen-decomposition of [Sll Slr; Slr' Srr]
            d = sqrt((obj.Sll - obj.Srr).^2 + 4*abs(obj.Slr).^2);
            l1 = (obj.Sll + obj.Srr + d)/2;
            l2 = max(0, (obj.Sll + obj.Srr - d)/2);
            leftDom = obj.Sll >= obj.Srr;
            u1 = l1 - obj.Srr; u2 = conj(obj.Slr);           % (C - l2 I) column, stable when L dominates
            u1(~leftDom) = obj.Slr(~leftDom); u2(~leftDom) = l1(~leftDom) - obj.Sll(~leftDom);
            nu = sqrt(abs(u1).^2 + abs(u2).^2) + eps;
            ops.u1 = u1./nu; ops.u2 = u2./nu;
            ops.g = (l1 - l2)./(l1 + eps);
            P = sqrt(ops.g) .* (conj(ops.u1).*L + conj(ops.u2).*R);
            gL = abs(ops.u1); gR = abs(ops.u2);
            % direction: tangent law on the content's base, stretched to Width
            th = atand((gL - gR)./(gL + gR + eps) * tand(obj.Base));
            th = th * obj.Width/obj.Base;
            ops.G = [];
            if strcmpi(obj.Mode, 'objects')
                [ops.W, ops.G] = obj.track_masks(th, abs(P).^2);
            else
                % constant-power panning between neighbouring slots
                u = (obj.Width - th)/(2*obj.Width)*(obj.K - 1);
                u = min(max(u, 0), obj.K - 1);
                i0 = min(floor(u), obj.K - 2); f = u - i0;
                ops.W = zeros(nb, obj.K);
                ops.W(sub2ind([nb obj.K], (1:nb)', i0+1)) = cos(f*pi/2);
                ops.W(sub2ind([nb obj.K], (1:nb)', i0+2)) = sin(f*pi/2);
            end
            S = obj.apply_ops(F, ops);
            AL = S(:, end-1); AR = S(:, end);
            if obj.Diagnostics
                e = abs(P).^2;
                [~, ~, bin] = histcounts(th, obj.AzEdges);
                ok = bin > 0;
                obj.DiagAz(:, end+1) = accumarray(bin(ok), e(ok), [numel(obj.AzEdges)-1, 1]);
                obj.DiagPA(:, end+1) = [sum(e); sum(abs(AL).^2 + abs(AR).^2)];
            end
            ifs = @(S) real(ifft([S; conj(S(end-1:-1:2,:))], [], 1)) .* obj.win;
            Yf = ifs(S);
            Ysh = {};
            if nargin > 1
                Ysh = cellfun(@(Fs) ifs(obj.apply_ops(Fs, ops)), Fsh, 'UniformOutput', false);
            end
        end

        function S = apply_ops(obj, F, ops)
            % The frame's linear operations on a spectrum F [bins x nIn]
            I = obj.layoutIdx;
            L = F(:, I.L); R = F(:, I.R);
            % power-complementary split of the principal component:
            % primary sqrt(g), its ambient share sqrt(1-g) (no level dips)
            Pu = conj(ops.u1).*L + conj(ops.u2).*R;
            P = sqrt(ops.g) .* Pu;
            r = (1 - sqrt(1 - ops.g)) .* Pu;
            AL = obj.AmbGain*(L - ops.u1.*r); AR = obj.AmbGain*(R - ops.u2.*r);
            S = [ops.W.*P, AL, AR];
            if strcmpi(obj.Mode, 'objects') && obj.Demix
                S(:, 1:obj.K) = obj.demix(S(:, 1:obj.K), ops.G, L, R, ops.g);
            end
            if ~isempty(I.C), S(:, obj.centre) = S(:, obj.centre) + F(:, I.C); end
            if ~isempty(I.SL), S(:, end-1) = S(:, end-1) + sum(F(:, I.SL), 2); end
            if ~isempty(I.SR), S(:, end)   = S(:, end)   + sum(F(:, I.SR), 2); end
        end

        function So = demix(obj, So, G, L, R, g)
            % Informed demixing of the two most likely sources per bin
            [Gs, ix] = sort(G, 2, 'descend');
            on = Gs(:,2) >= obj.DemixRatio*Gs(:,1) & Gs(:,2) > 0;
            if ~any(on), return; end
            b = find(on); i = ix(b,1); j = ix(b,2);
            % panning vectors of the two tracks (output azimuth -> content
            % azimuth -> tangent law -> unit-power gains)
            av = @(az) obj.pan_vec(az*obj.Base/obj.Width);
            [ai1, ai2] = av(obj.trAz(i)'); [aj1, aj2] = av(obj.trAz(j)');
            % regularised LS: s = (A'A + eps I)^-1 A' x,  A = [ai aj] (real)
            m11 = ai1.^2 + ai2.^2; m22 = aj1.^2 + aj2.^2; m12 = ai1.*aj1 + ai2.*aj2;
            ep = obj.DemixReg*(m11 + m22)/2;
            m11 = m11 + ep; m22 = m22 + ep;
            det = m11.*m22 - m12.^2;
            xl = L(b); xr = R(b);
            y1 = ai1.*xl + ai2.*xr; y2 = aj1.*xl + aj2.*xr;
            s1 = ( m22.*y1 - m12.*y2)./det;
            s2 = (-m12.*y1 + m11.*y2)./det;
            sg = sqrt(g(b));                        % keep the ambience out
            So(b,:) = 0;
            So(sub2ind(size(So), b, i)) = sg.*s1;
            So(sub2ind(size(So), b, j)) = sg.*s2;
        end

        function [a1, a2] = pan_vec(obj, az)
            r = tand(max(min(az, obj.Base), -obj.Base))/tand(obj.Base);
            a1 = (1 + r)/2; a2 = (1 - r)/2;
            n = sqrt(a1.^2 + a2.^2);
            a1 = a1./n; a2 = a2./n;
        end

        function [W, G] = track_masks(obj, th, e)
            % Track the sources' directions; amplitude masks [bins x K]
            hop = obj.hop; fs = obj.fs;
            [~, ~, bin] = histcounts(th, obj.AzEdges);
            f = (0:numel(th)-1)'*fs/obj.N;
            ok = bin > 0 & f >= obj.HistBand(1) & f <= obj.HistBand(2);
            inst = accumarray(bin(ok), e(ok), size(obj.hist));
            a = exp(-hop/(fs*obj.TrackTau));
            obj.hist = a*obj.hist + (1-a)*inst;
            h = conv(obj.hist, obj.kern, 'same');
            azc = (obj.AzEdges(1:end-1) + obj.AzEdges(2:end))'/2;
            % peaks: local maxima above the floor, strongest first, MinSep apart
            pk = find(h(2:end-1) >= h(1:end-2) & h(2:end-1) > h(3:end)) + 1;
            pk = pk(h(pk) > max(h)*10^(obj.PeakFloorDb/10) & h(pk) > 0);
            [~, o] = sort(h(pk), 'descend'); pk = pk(o);
            peaks = [];
            for i = pk(:)'
                if all(abs(azc(i) - peaks) >= obj.MinSep), peaks(end+1) = azc(i); end %#ok<AGROW>
                if numel(peaks) == obj.K, break; end
            end
            % association: nearest pairs first, within the gate
            aOn  = 1 - exp(-hop/(fs*0.15));             % activity attack
            aOff = exp(-hop/(fs*0.5));                  % release (position held)
            matched = false(1, obj.K); used = false(size(peaks));
            if ~isempty(peaks)
                D = abs(obj.trAz(:) - peaks(:)');
                D(obj.trAct(:) < 0.05, :) = inf;        % dead tracks do not continue
                while true
                    [dmin, ix] = min(D(:));
                    if isempty(dmin) || dmin > obj.Gate, break; end
                    [j, q] = ind2sub(size(D), ix);
                    obj.trAz(j) = obj.trAz(j) + 0.3*(peaks(q) - obj.trAz(j));
                    matched(j) = true; used(q) = true;
                    D(j,:) = inf; D(:,q) = inf;
                end
            end
            obj.trAct(matched) = obj.trAct(matched) + aOn*(1 - obj.trAct(matched));
            obj.trAct(~matched) = aOff*obj.trAct(~matched);
            for q = find(~used)                         % new source: the weakest unmatched track
                [amin, j] = min(obj.trAct + 2*matched);
                if amin < 0.2
                    obj.trAz(j) = peaks(q); obj.trAct(j) = 0.3; matched(j) = true;
                end
            end
            % two tracks on the same source: keep the more active one
            for j = 1:obj.K
                for q = j+1:obj.K
                    if obj.trAct(j) > 0.2 && obj.trAct(q) > 0.2 && abs(obj.trAz(j) - obj.trAz(q)) < obj.MinSep/2
                        if obj.trAct(j) >= obj.trAct(q), obj.trAct(q) = 0; else, obj.trAct(j) = 0; end
                    end
                end
            end
            active = obj.trAct > 0.2;
            G = exp(-(th - obj.trAz).^2/(2*obj.TrackSigma^2)) .* active;
            W = sqrt(G ./ max(sum(G, 2), eps));
            G = G ./ max(sum(G, 2), eps);
            obj.Pos(1:obj.K,:) = separator_base.dir2pos(obj.trAz(:));
            if obj.Diagnostics
                lg = obj.trAz(:); lg(~active) = NaN;
                obj.TrackLog(:, end+1) = lg;
            end
        end
    end
end
