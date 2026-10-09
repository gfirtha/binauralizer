classdef cluster_separator < separator_base
    %CLUSTER_SEPARATOR  Sources by clustering time-frequency bins on their
    %local wavenumber direction, then a local wavenumber vector per source.
    %
    %   1. Per STFT bin (N, hop N/2): local wavenumber vector k^ of the field
    %      the channels create at the listener (point-source model, see
    %      wavenumber_analyzer); stereo is Lt/Rt-decoded first
    %      (ltrt_decoder, steering DecTau). Feature: the direction it comes
    %      from, theta; weight |P|^2.
    %   2. Segmentation: online EM of a mixture of von Mises components on
    %      theta plus one uniform 'diffuse' component (the bed), sufficient
    %      statistics with forgetting (TauEM), so a source is what keeps
    %      explaining energy from one direction, not a single frame's peak.
    %      - birth: the diffuse-explained energy's direction histogram
    %        (smoothed TauBirth) has a peak holding BirthFrac of the
    %        energy, away from the components - or an onset there (OnsetDb
    %        above its slow level, at least OnsetFrac of the frame);
    %      - death: weight below MinWeight for DeathTime;
    %      - merge: two means closer than MergeDeg.
    %   3. Per source: its local wavenumber vector from its own bins
    %      (posterior x |P|^2 weighted k^): direction and |k^|. A Kalman
    %      filter per source ([az; angular velocity], two noise levels as
    %      in soundfield_tracker) filters it and PREDICTS the source's
    %      direction for the next frame's E-step (fast fly-bys stay in
    %      their component).
    %   4. Synthesis: soft masks = posteriors; a source = an object (fixed
    %      output slot while it lives, fade in/out): the omni pressure at
    %      the origin masked and scaled to the masked channel power per
    %      band; the diffuse posterior keeps the channels (bed).

    properties
        K = 6               % max sources (object outputs)
        N = 1024
        R = 2; RefDist = 3; c = 343
        MinHz = 60
        DecTau = 0.2
        TauEM = 0.5         % [s] forgetting of the component statistics
        TauBirth = 0.2      % [s] smoothing of the diffuse-residual histogram
        BirthFrac = 0.05    % residual peak share of the energy for a birth
        OnsetDb = 10        % ... or an onset this far above the slow residual level
        OnsetFrac = 0.1     % ... carrying this share of the frame energy
        MinWeight = 0.02
        DeathTime = 0.5     % [s]
        MergeDeg = 15
        SigmaInit = 15      % [deg] width of a new component
        SigmaRange = [4 15] % [deg] a source is narrow; wider spread is diffuse
        DiffuseMin = 0.2    % minimum prior of the diffuse component
        DiffuseInit = 0.5   % prior weight of the diffuse component
        HistDeg = 5
        AccNoise = 2e4; AccNoiseStable = 300   % [deg/s^2] Kalman (moving / standing)
        MeasDeg = 3
        StableDeg = 20; StableTime = 0.1
        FadeIn = 0.02; FadeOut = 0.1
        PosTau = 0.015; PosTauStill = 0.5   % [s] object position smoothing (moving / standing)
        AmbGain = 1; AmbAz = 110; Width = 180; Base = 30
        Diagnostics = false
    end

    properties (SetAccess = private)
        DiagAz = []; DiagPA = []
        AzEdges = -180:4:180
        TrackLog = []       % [K x frames] source directions (NaN = no source)
        CountLog = []       % [1 x frames] live sources
        KmagLog = []        % [K x frames] |k^| of each source
    end

    properties (Access = private)
        hop; win; frame; acc; pending; outq
        stereo; Dec; chIdx; nBed; S0; nbU; useIdx; bandOfFull; Bmap; nB
        comp                % struct array of live components
        diffS               % diffuse statistic
        resF; resS          % residual histograms (fast / slow)
        slotOf; slotG; slotAz; nextId; t
    end

    methods
        function obj = cluster_separator(varargin)
            for k = 1:2:numel(varargin), obj.(varargin{k}) = varargin{k+1}; end
        end

        function init(obj, fs, nIn)
            obj.fs = fs; obj.nIn = nIn;
            obj.stereo = nIn == 2;
            if obj.stereo
                obj.Dec = ltrt_decoder('TauSteer', obj.DecTau); obj.Dec.setup_bands(obj.N, fs);
                az = [30 0 -30 110 -110]; names = {'L','C','R','Ls','Rs'};
            else
                [az, names] = separator_base.channel_layout(nIn);
                obj.chIdx = find(~isnan(az));
                az = az(obj.chIdx); names = names(obj.chIdx);
                if nIn == 1, az = 0; end
            end
            U = separator_base.dir2pos(az);
            obj.nBed = numel(az);
            obj.hop = obj.N/2;
            obj.win = sqrt(0.5*(1 - cos(2*pi*(0:obj.N-1)'/obj.N)));
            be = [0 100 200 300 400 510 630 770 920 1080 1270 1480 1720 2000 2320 ...
                  2700 3150 3700 4400 5300 6400 7700 9500 12000 15500 inf];
            obj.S0 = wavenumber_analyzer.precompute(obj.N, fs, U, obj.R*ones(1, numel(az)), obj.win, be, obj.MinHz, obj.c, 'bands');
            obj.useIdx = find(obj.S0.use); obj.nbU = numel(obj.useIdx);
            obj.Bmap = obj.S0.Bin; obj.nB = size(obj.Bmap, 1);
            [bi, col] = find(obj.Bmap);
            obj.bandOfFull = zeros(obj.N/2 + 1, 1); obj.bandOfFull(obj.useIdx(col)) = bi;
            obj.Names = [arrayfun(@(k) sprintf('Src %d', k), 1:obj.K, 'UniformOutput', false), ...
                strcat('Amb', {' '}, names(:)')];
            obj.Pos = [separator_base.dir2pos(linspace(0, 300, obj.K)')*obj.R/obj.RefDist; separator_base.dir2pos(az(:))];
            obj.Latency = obj.N - obj.hop;
            obj.reset();
        end

        function reset(obj)
            obj.frame = zeros(obj.N, obj.nIn);
            obj.acc = zeros(obj.N, obj.nOut());
            obj.pending = zeros(0, obj.nIn);
            obj.outq = zeros(0, obj.nOut());
            if obj.stereo, obj.Dec.reset_steering(); end
            obj.comp = struct('id', {}, 'C', {}, 'S', {}, 'W', {}, 'mu', {}, 'kap', {}, ...
                'x', {}, 'P', {}, 'low', {}, 'still', {}, 'anchor', {}, 'kmag', {});
            obj.diffS = 0;
            nA = 360/obj.HistDeg; obj.resF = zeros(nA, 1); obj.resS = [];
            obj.slotOf = zeros(1, obj.K); obj.slotG = zeros(1, obj.K);
            obj.slotAz = linspace(0, 300, obj.K); obj.nextId = 1; obj.t = 0;
            obj.DiagAz = []; obj.DiagPA = []; obj.TrackLog = []; obj.CountLog = []; obj.KmagLog = [];
        end

        function Y = process(obj, X)
            n = size(X,1);
            obj.pending = [obj.pending; X];
            p = 0;
            while size(obj.pending,1) - p >= obj.hop
                obj.frame = [obj.frame(obj.hop+1:end,:); obj.pending(p + (1:obj.hop),:)];
                p = p + obj.hop;
                obj.acc = obj.acc + obj.synth_frame();
                obj.outq = [obj.outq; obj.acc(1:obj.hop,:)];
                obj.acc = [obj.acc(obj.hop+1:end,:); zeros(obj.hop, obj.nOut())];
            end
            obj.pending = obj.pending(p+1:end,:);
            m = min(n, size(obj.outq,1));
            Y = [zeros(n-m, obj.nOut()); obj.outq(1:m,:)];
            obj.outq = obj.outq(m+1:end,:);
        end

        function Y = separate_offline(obj, Xs)
            obj.reset();
            x = [Xs{1}; zeros(obj.Latency + obj.hop, size(Xs{1},2))];
            n = size(Xs{1},1); Y = zeros(size(x,1), obj.nOut());
            for b = 1:4096:size(x,1)
                r = b:min(size(x,1), b+4095); Y(r,:) = obj.process(x(r,:));
            end
            Y = {Y(obj.Latency + (1:n), :)};
        end
    end

    methods (Access = private)
        function Yf = synth_frame(obj)
            nb = obj.N/2 + 1; dt = obj.hop/obj.fs; K = obj.K;
            obj.t = obj.t + dt;
            F = fft(obj.frame .* obj.win, [], 1); F = F(1:nb, :);
            if obj.stereo
                Xall = obj.Dec.decode_spectrum(F(:,1), F(:,2), dt);
            else
                Xall = F(:, obj.chIdx);
            end
            X = Xall(obj.useIdx, :);
            q = wavenumber_analyzer.bin_stats(X, obj.S0);
            th = atan2d(q.kx, -q.ky);                       % direction each bin comes from
            w = q.E;
            Et = sum(w);
            a = exp(-dt/obj.TauEM);

            % ---- predict the components (Kalman) ----
            for j = 1:numel(obj.comp)
                obj.comp(j) = obj.kalman_predict(obj.comp(j), dt);
            end

            % ---- E-step: posteriors of the components and the diffuse part ----
            nC = numel(obj.comp);
            tot = obj.diffS + sum([obj.comp.W]) + eps;
            L = zeros(obj.nbU, nC + 1);
            L(:, end) = max(obj.diffS/tot, obj.DiffuseMin)/(2*pi);
            for j = 1:nC
                cj = obj.comp(j);
                % von Mises: exp(kap (cos - 1)) / (2 pi I0(kap) e^-kap)
                L(:, j) = (cj.W/tot) * exp(cj.kap*(cosd(th - cj.mu) - 1)) / (2*pi*besseli(0, cj.kap, 1));
            end
            Gm = L./max(sum(L, 2), realmin);                % posteriors [bins x (nC+1)]
            Gm(w == 0, :) = 0; Gm(w == 0, end) = 1;

            % ---- M-step with forgetting; the source's own wavenumber vector ----
            for j = 1:nC
                gw = Gm(:, j).*w;
                cj = obj.comp(j);
                cj.C = a*cj.C + (1-a)*sum(gw.*cosd(th));
                cj.S = a*cj.S + (1-a)*sum(gw.*sind(th));
                cj.W = a*cj.W + (1-a)*sum(gw);
                Rb = min(hypot(cj.C, cj.S)/max(cj.W, eps), 0.999);
                kap = Rb*(2 - Rb^2)/(1 - Rb^2);
                sig = min(max(sqrt(1/max(kap, eps))*180/pi, obj.SigmaRange(1)), obj.SigmaRange(2));
                cj.kap = 1/(sig*pi/180)^2;
                % the source's local wavenumber vector from its own bins
                % its core bins only (within 2 sigma of the prediction)
                gc = gw.*(abs(mod(th - cj.mu + 180, 360) - 180) < 2*sqrt(1/cj.kap)*180/pi);
                sw = sum(gc);
                if sw > 1e-3*Et && sw > 0
                    kx = sum(gc.*q.kx)/sw; ky = sum(gc.*q.ky)/sw;
                    cj.kmag = hypot(kx, ky);
                    % measurement noise: the component's width over its effective bin count
                    neff = sw^2/max(sum(gc.^2), eps);
                    cj = obj.kalman_update(cj, atan2d(kx, -ky), max(obj.MeasDeg, sqrt(1/cj.kap)*180/pi/sqrt(neff)));
                end
                cj.mu = cj.x(1);
                obj.comp(j) = cj;
            end
            obj.diffS = a*obj.diffS + (1-a)*sum(Gm(:, end).*w);

            % ---- birth from the diffuse residual ----
            nA = 360/obj.HistDeg; edges = -180:obj.HistDeg:180;
            [~, ~, bi] = histcounts(th, edges); ok = bi > 0;
            inst = accumarray(bi(ok), Gm(ok, end).*w(ok), [nA 1]);
            inst = conv([inst(end-1:end); inst; inst(1:2)], [1 2 3 2 1]'/9, 'valid');
            ab = exp(-dt/obj.TauBirth);
            obj.resF = ab*obj.resF + (1-ab)*inst;
            if isempty(obj.resS), obj.resS = obj.resF; end
            obj.resS = exp(-dt/1)*obj.resS + (1 - exp(-dt/1))*inst;
            EtS = sum(obj.resF) + sum([obj.comp.W]) + eps;   % smoothed total energy rate
            [pk, ip] = max(obj.resF);
            az0 = edges(ip) + obj.HistDeg/2;
            far = all(abs(mod([obj.comp.mu] - az0 + 180, 360) - 180) > 2*obj.SigmaInit);
            onset = inst(ip) > 10^(obj.OnsetDb/10)*obj.resS(ip) && inst(ip) > obj.OnsetFrac*Et;
            if numel(obj.comp) < obj.K && far && (pk > obj.BirthFrac*EtS || onset)
                kap = 1/(obj.SigmaInit*pi/180)^2;
                obj.comp(end+1) = struct('id', obj.nextId, 'C', pk*cosd(az0), 'S', pk*sind(az0), 'W', pk, ...
                    'mu', az0, 'kap', kap, 'x', [az0; 0], 'P', diag([obj.MeasDeg^2, 300^2]), ...
                    'low', 0, 'still', 0, 'anchor', az0, 'kmag', 1);
                obj.nextId = obj.nextId + 1;
                obj.resF(:) = 0;                             % this residual now has a component
            end

            % ---- death and merge ----
            Wall = sum([obj.comp.W]) + obj.diffS + eps;
            keep = true(1, numel(obj.comp));
            for j = 1:numel(obj.comp)
                if obj.comp(j).W/Wall < obj.MinWeight, obj.comp(j).low = obj.comp(j).low + dt; else, obj.comp(j).low = 0; end
                if obj.comp(j).low > obj.DeathTime, keep(j) = false; end
            end
            for i = 1:numel(obj.comp)
                for j = i+1:numel(obj.comp)
                    if keep(i) && keep(j) && abs(mod(obj.comp(i).mu - obj.comp(j).mu + 180, 360) - 180) < obj.MergeDeg
                        % the younger one joins the older one
                        obj.comp(i).C = obj.comp(i).C + obj.comp(j).C; obj.comp(i).S = obj.comp(i).S + obj.comp(j).S;
                        obj.comp(i).W = obj.comp(i).W + obj.comp(j).W; keep(j) = false;
                    end
                end
            end
            obj.comp = obj.comp(keep);

            % ---- object slots: a component keeps its slot while it lives ----
            ids = [obj.comp.id];
            target = zeros(1, K);
            for k = 1:K
                j = find(ids == obj.slotOf(k), 1);
                if ~isempty(j)
                    target(k) = 1;
                    % standing source: its object stays put; moving: follows at once
                    if obj.comp(j).still >= obj.StableTime, tp = obj.PosTauStill; else, tp = obj.PosTau; end
                    da = mod(obj.comp(j).mu - obj.slotAz(k) + 180, 360) - 180;
                    obj.slotAz(k) = mod(obj.slotAz(k) + da*(1 - exp(-dt/tp)) + 180, 360) - 180;
                elseif obj.slotOf(k) && obj.slotG(k) < 0.01
                    obj.slotOf(k) = 0;
                end
            end
            for j = 1:numel(obj.comp)
                if ~any(obj.slotOf == ids(j))
                    k = find(obj.slotOf == 0, 1);
                    if isempty(k), continue; end
                    obj.slotOf(k) = ids(j); obj.slotAz(k) = obj.comp(j).mu; target(k) = 1;
                end
            end
            up = target > obj.slotG;
            tau = obj.FadeOut*ones(1, K); tau(up) = obj.FadeIn;
            obj.slotG = obj.slotG + (target - obj.slotG).*(1 - exp(-dt./tau));

            % ---- synthesis: masked omni per source, diffuse -> bed ----
            P0 = sum(Xall, 2);
            Pc = sum(abs(X).^2, 2);                          % channel power per used bin
            S = zeros(nb, K + obj.nBed);
            gDiff = Gm(:, end);
            for k = 1:K
                j = find(ids == obj.slotOf(k), 1);
                if isempty(j) || obj.slotG(k) < 0.01, continue; end
                gk = Gm(:, j)*obj.slotG(k)^2;
                % band-wise power match: masked omni -> masked channel power
                num = obj.Bmap*(gk.*Pc); den = obj.Bmap*(gk.*abs(P0(obj.useIdx)).^2);
                gb = sqrt(num./max(den, eps)); gb(den < 0.1*num) = 0;   % cancelling omni: leave to the bed
                gBin = gb(obj.bandOfFull(obj.useIdx));
                S(obj.useIdx, k) = sqrt(gk).*gBin.*P0(obj.useIdx);
                gDiff = gDiff + Gm(:, j).*(1 - obj.slotG(k)^2) + Gm(:, j).*obj.slotG(k)^2.*(gBin == 0);
            end
            bf = ones(nb, 1); bf(obj.useIdx) = min(gDiff, 1);
            S(:, K+1:end) = obj.AmbGain*sqrt(bf).*Xall;
            obj.Pos(1:K,:) = separator_base.dir2pos(obj.slotAz(:))*obj.R/obj.RefDist;

            if obj.Diagnostics
                [~, ~, b2] = histcounts(th, obj.AzEdges); ok = b2 > 0;
                obj.DiagAz(:, end+1) = accumarray(b2(ok), w(ok), [numel(obj.AzEdges)-1, 1]);
                obj.DiagPA(:, end+1) = [sum(sum(abs(S(:,1:K)).^2)); sum(sum(abs(S(:,K+1:end)).^2))];
                lg = obj.slotAz(:); lg(obj.slotG < 0.01) = NaN; obj.TrackLog(:, end+1) = lg;
                km = nan(K, 1);
                for k = 1:K, j = find(ids == obj.slotOf(k), 1); if ~isempty(j), km(k) = obj.comp(j).kmag; end, end
                obj.KmagLog(:, end+1) = km;
                obj.CountLog(end+1) = numel(obj.comp);
            end
            Yf = real(ifft([S; conj(S(end-1:-1:2,:))], [], 1)) .* obj.win;
        end

        function cj = kalman_predict(obj, cj, dt)
            if abs(mod(cj.x(1) - cj.anchor + 180, 360) - 180) <= obj.StableDeg
                cj.still = cj.still + dt;
            else
                cj.anchor = cj.x(1); cj.still = 0;
            end
            if cj.still >= obj.StableTime, q = obj.AccNoiseStable^2; else, q = obj.AccNoise^2; end
            Fm = [1 dt; 0 1];
            cj.x = Fm*cj.x; cj.x(1) = mod(cj.x(1) + 180, 360) - 180;
            cj.P = Fm*cj.P*Fm' + q*[dt^3/3 dt^2/2; dt^2/2 dt];
            cj.mu = cj.x(1);
        end

        function cj = kalman_update(~, cj, z, sdev)
            v = mod(z - cj.x(1) + 180, 360) - 180;
            Sv = cj.P(1,1) + sdev^2;
            Kg = cj.P(:,1)/Sv;
            cj.x = cj.x + Kg*v; cj.x(1) = mod(cj.x(1) + 180, 360) - 180;
            cj.P = cj.P - Kg*cj.P(1,:);
        end
    end
end
