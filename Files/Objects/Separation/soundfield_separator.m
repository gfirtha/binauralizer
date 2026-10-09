classdef soundfield_separator < separator_base
    %SOUNDFIELD_SEPARATOR  Real-time separation of a channel mix into
    %moving objects from the virtual sound field at the listener: sector
    %analysis (several simultaneous directions per band, see
    %soundfield_analyzer) + multi-object tracking (soundfield_tracker),
    %plus an ambience bed.
    %
    %   Per STFT frame (sqrt-Hann, N, hop N/2) and band:
    %   - every channel is a plane wave from its loudspeaker direction;
    %     sector pressures P_s / velocities V_s (order SectorOrder max-rE
    %     beams, sum_s w_s = 1), coherent neighbour sectors merged into
    %     groups (MergeCoh); per group direction, energy E_G and
    %     diffuseness psi_G (statistics smoothed over Tau);
    %   - detection: groups above DetBand(1) Hz with psi < PsiMax ->
    %     direction histogram -> peaks -> tracks (Kalman, novelty births,
    %     see soundfield_tracker); confirmed tracks hold the K object slots
    %     (a slot keeps its track until the track dies);
    %   - synthesis: the direct part of each group is shared among the
    %     objects within AssignDeg (Gaussian weights AssignSigma x the
    %     object's fade gain^2, power complementary, the rest to the bed),
    %     as its pressure P_G = sum_{s in G} P_s scaled so that the group's
    %     share of the band's channel power is kept:
    %     |Y|^2 = sum_i |X_i|^2 E_G/sum E (1 - psi_G) w; the diffuse parts
    %     and unassigned direct parts stay in the bed: the input channels
    %     scaled by sqrt(their power fraction).
    %   - objects: a slot takes a confirmed track that is novel or MinAge
    %     old, fades in (FadeIn), follows it (PosTau), and after the track
    %     died fades out (FadeOut) at its last position before it can be
    %     reused - no jumps of a sounding object.
    %   Objects move with their tracks (Pos per block). Bed objects sit at
    %   the input layout's directions (stereo: +-AmbAz); LFE is dropped.
    %   Latency N - hop. No block switching in the synthesis (the
    %   analyzer's short windows are analysis only); N = 1024 by default
    %   (10.7 ms hop at 48 kHz).

    properties
        K = 6               % object slots
        N = 1024            % FFT size (hop N/2)
        SectorOrder = 2
        NSectors = []       % [] = 2 SectorOrder + 2
        MergeCoh = 0.8
        Tau = 0.02          % [s] statistics smoothing (diffuseness, coherence)
        BandsPerOct = 3     % analysis bands (above 150 Hz; one band below)
        DetBand = [1500 Inf] % [Hz] bands that detect / track directions
        PsiMax = 0.5        % more diffuse groups do not detect
        AssignDeg = 40      % a group joins objects within this ...
        AssignSigma = 15    % ... with Gaussian weights (deg), power complementary
        MinAge = 0.1        % [s] a non-novel track needs this age to get an object
        FadeIn = 0.02       % [s] object fade-in when its track starts
        FadeOut = 0.1       % [s] fade-out after its track died (slot free afterwards)
        PosTau = 0.015      % [s] smoothing of the object position (moving track)
        PosTauStill = 0.5   % [s] ... of a standing track
        AmbGain = 1
        AmbAz = 110         % bed direction of a stereo input
        Width = 180         % (panel: direction map range)
        Base = 30           % (unused; common separator settings)
        Tracker             % soundfield_tracker (settings can be changed before init)
        Diagnostics = false
    end

    properties (SetAccess = private)
        DiagAz = []         % [nAzBins x frames] direct power per direction
        DiagPA = []         % [2 x frames] object / bed power
        AzEdges = -180:4:180
        TrackLog = []       % [K x frames] slot azimuths (NaN = empty)
    end

    properties (Access = private)
        hop; win; frame; acc; pending; outq
        chIdx               % input channels used (LFE out)
        W; WU; Ks; M; nB; bandOf; bandLo
        Sa                  % smoothed band cross statistics
        Slots               % object_slots
        t                   % [s]
        nBed
    end

    methods
        function obj = soundfield_separator(varargin)
            % objects closer than 20 deg merge; standing tracks barely accelerate
            % standing = within 20 deg for 0.1 s (a dialogue track wanders ~10 deg
            % with the mix, a fly-by leaves 20 deg in ~20 ms)
            obj.Tracker = soundfield_tracker('DupDeg', 20, 'AccNoiseStable', 100, ...
                'StableDeg', 20, 'StableTime', 0.1);
            for k = 1:2:numel(varargin), obj.(varargin{k}) = varargin{k+1}; end
        end

        function init(obj, fs, nIn)
            obj.fs = fs; obj.nIn = nIn;
            obj.K = max(1, round(obj.K));
            [az, names] = separator_base.channel_layout(nIn);
            obj.chIdx = find(~isnan(az));
            az = az(obj.chIdx); names = names(obj.chIdx);
            if nIn == 1, az = 0; end
            U = separator_base.dir2pos(az);
            [obj.W, obj.WU] = soundfield_analyzer.sector_weights(U, obj.SectorOrder, obj.NSectors);
            obj.Ks = size(obj.W, 2);
            % bands: one below 150 Hz, BandsPerOct per octave up to fs/2
            nb = obj.N/2 + 1; f = (0:nb-1)'*fs/obj.N;
            be = [0, 150*2.^(0:1/obj.BandsPerOct:log2(fs/2/150)), fs];
            [~, ~, obj.bandOf] = histcounts(f, be);
            used = unique(obj.bandOf);
            obj.bandOf = arrayfun(@(b) find(used == b), obj.bandOf);   % no empty bands
            obj.nB = numel(used);
            obj.M = sparse(obj.bandOf, (1:nb)', 1, obj.nB, nb);
            obj.bandLo = accumarray(obj.bandOf(:), f, [obj.nB 1], @min);   % lowest bin frequency per band
            % outputs: K objects + bed
            bedAz = az(:);
            if numel(az) == 2, bedAz = [obj.AmbAz; -obj.AmbAz]; end
            obj.nBed = numel(bedAz);
            obj.Names = [arrayfun(@(k) sprintf('Obj %d', k), 1:obj.K, 'UniformOutput', false), ...
                strcat('Amb', {' '}, names(:)')];
            obj.Pos = separator_base.dir2pos([linspace(0, 360, obj.K + 1)'; bedAz]);
            obj.Pos = obj.Pos([1:obj.K, obj.K + 2:end], :);
            obj.hop = obj.N/2;
            obj.win = sqrt(0.5*(1 - cos(2*pi*(0:obj.N-1)'/obj.N)));   % periodic sqrt-Hann
            obj.Latency = obj.N - obj.hop;
            obj.Tracker.Band = obj.DetBand;
            obj.Tracker.PsiMax = obj.PsiMax;
            obj.reset();
        end

        function reset(obj)
            obj.frame = zeros(obj.N, obj.nIn);
            obj.acc = zeros(obj.N, obj.nOut());
            obj.pending = zeros(0, obj.nIn);
            obj.outq = zeros(0, obj.nOut());
            obj.Sa = [];
            obj.Slots = object_slots('K', obj.K, 'MinAge', obj.MinAge, 'FadeIn', obj.FadeIn, ...
                'FadeOut', obj.FadeOut, 'PosTau', obj.PosTau, 'PosTauStill', obj.PosTauStill, ...
                'AssignSigma', obj.AssignSigma, 'AssignDeg', obj.AssignDeg, ...
                'DupDeg', obj.Tracker.DupDeg, 'StableTime', obj.Tracker.StableTime);
            obj.t = 0;
            obj.Tracker.reset();
            obj.DiagAz = []; obj.DiagPA = []; obj.TrackLog = [];
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

        function Y = separate_offline(obj, Xs)
            % Whole signal, aligned (no latency): Xs{1} only (the decisions
            % are not linear operations that could be shared with parts)
            obj.reset();
            x = Xs{1}; n = size(x,1); L = obj.Latency;
            Y = obj.process([x; zeros(L + obj.hop, size(x,2))]);
            Y = {Y(L + (1:n), :)};
        end
    end

    methods (Access = private)
        function Yf = synth_frame(obj)
            nb = obj.N/2 + 1; K = obj.K; Ks = obj.Ks; dt = obj.hop/obj.fs;
            obj.t = obj.t + dt;
            F = fft(obj.frame .* obj.win, [], 1);
            X = F(1:nb, obj.chIdx);
            Sb = soundfield_analyzer.sector_cross(X, obj.W, obj.WU, obj.M);
            a = exp(-obj.hop/(obj.fs*obj.Tau));
            if isempty(obj.Sa), obj.Sa = Sb; else, obj.Sa = a*obj.Sa + (1-a)*Sb; end
            Pc = obj.M*sum(abs(X).^2, 2);                   % band channel power
            P = X*obj.W;                                    % sector pressures
            % ---- groups per band ----
            G = cell(obj.nB, 1);
            dAz = []; dPw = [];
            for b = 1:obj.nB
                [grp, I, E, Ism, Esm] = soundfield_analyzer.group_sectors(Sb(b,:), obj.Sa(b,:), Ks, obj.MergeCoh);
                psi = min(max(1 - sqrt(sum(Ism.^2, 2))./max(Esm, eps), 0), 1);
                az = atan2d(-I(:,1), I(:,2));               % instantaneous: no lag for fast movers
                azS = atan2d(-Ism(:,1), Ism(:,2));          % smoothed: no flicker for standing ones
                K2 = Ks*Ks;
                pg = zeros(numel(grp), 1);                  % |P_G|^2 band sums
                for gi = 1:numel(grp)
                    m = grp{gi}; lin = reshape(m(:) + (m(:)' - 1)*Ks, [], 1);
                    pg(gi) = sum(Sb(b, 3*K2 + lin));
                end
                G{b} = struct('grp', {grp}, 'az', az, 'azS', azS, 'E', E, 'psi', psi, 'pg', pg);
                if obj.bandLo(b) >= obj.DetBand(1) && obj.bandLo(b) <= obj.DetBand(2)
                    ok = psi < obj.PsiMax & E > 0;
                    dAz = [dAz; az(ok)]; dPw = [dPw; E(ok).*(1 - psi(ok))]; %#ok<AGROW>
                end
            end
            % ---- detection + tracking -> object slots ----
            det = obj.Tracker.detect_frame(dAz, dPw, dt);
            tr = obj.Tracker.step(det, obj.t, dt);
            obj.Slots.update(tr, obj.t, dt);
            active = obj.Slots.active;
            % ---- synthesis: soft, power-complementary assignment ----
            S = zeros(nb, K + obj.nBed);
            bedFrac = zeros(obj.nB, 1);
            for b = 1:obj.nB
                g = G{b}; Et = sum(g.E);
                if Et <= 0, bedFrac(b) = 1; continue; end
                bins = obj.bandOf == b;
                amb = sum(g.E .* g.psi);
                for gi = 1:numel(g.grp)
                    direct = g.E(gi)*(1 - g.psi(gi));
                    % instantaneous direction for a moving object, smoothed
                    % for a standing one (see object_slots.weights)
                    [w, sw] = obj.Slots.weights(g.azS(gi), g.az(gi));
                    amb = amb + (1 - sw)*direct;            % the rest of the direct part: bed
                    if sw > 0 && g.pg(gi) > 0
                        pG = sum(P(bins, g.grp{gi}), 2);
                        for k = find(w > 0)
                            S(bins, k) = S(bins, k) + sqrt(Pc(b)*direct*w(k)/Et/g.pg(gi))*pG;
                        end
                    end
                end
                bedFrac(b) = amb/Et;
            end
            gb = obj.AmbGain*sqrt(min(bedFrac(obj.bandOf), 1));
            Xb = X;
            if size(Xb, 2) ~= obj.nBed, Xb = Xb(:, 1:obj.nBed); end
            S(:, K+1:end) = gb .* Xb;
            % ---- positions, diagnostics ----
            obj.Pos(1:K,:) = separator_base.dir2pos(obj.Slots.az(:));
            if obj.Diagnostics
                allAz = []; allPw = [];
                for b = 1:obj.nB
                    allAz = [allAz; G{b}.az]; allPw = [allPw; G{b}.E.*(1 - G{b}.psi)]; %#ok<AGROW>
                end
                [~, ~, bi] = histcounts(allAz, obj.AzEdges);
                ok = bi > 0;
                obj.DiagAz(:, end+1) = accumarray(bi(ok), allPw(ok), [numel(obj.AzEdges)-1, 1]);
                obj.DiagPA(:, end+1) = [sum(sum(abs(S(:,1:K)).^2)); sum(sum(abs(S(:,K+1:end)).^2))];
                lg = obj.Slots.az(:); lg(~active) = NaN;
                obj.TrackLog(:, end+1) = lg;
            end
            ifs = @(S) real(ifft([S; conj(S(end-1:-1:2,:))], [], 1)) .* obj.win;
            Yf = ifs(S);
        end
    end
end
