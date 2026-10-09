classdef phantom_object_separator < separator_base
    %PHANTOM_OBJECT_SEPARATOR  Phantom objects + fast trajectories + music bed.
    %
    %   Speaker-feed covariance finds the current signal subspace. Rank-1
    %   content is matched directly to a phantom-gain dictionary; higher-rank
    %   content is rotated with a small ICA step, then each component vector
    %   is matched to the same dictionary. A high-rate tracker turns the
    %   frame candidates into moving object slots. The slot signals are
    %   projected out of the loudspeaker feed and the residual is emitted as
    %   one music bed, panned by its current dominant phantom direction.

    properties
        K = 5
        N = 512
        Hop = 128
        MusicBed = true
        MaxComponents = 4
        AmbGain = 1
        Tau = 0.03
        EnergyTau = 0.4
        OnsetDb = 5
        SharpnessDb = 5
        MinEnergyDb = -50
        CoherentThreshold = 0.25
        IcaThreshold = 0.35
        AssignSigma = 18
        AssignDeg = 55
        Width = 180
        Base = 30
        AmbAz = 110
        MaxFrameDet = 2
        MinCandidateSep = 35
        Tracker
        Diagnostics = false
    end

    properties (SetAccess = private)
        DiagAz = []
        DiagPA = []
        AzEdges = -180:4:180
        TrackLog = []
        CountLog = []
    end

    properties (Access = private)
        hop
        win
        frame
        acc
        pending
        outq
        chIdx
        az
        templates
        Slots
        t = 0
        refEnergy = 0
        slowEnergy = 0
        lastBedAz = 0
    end

    methods
        function obj = phantom_object_separator(varargin)
            obj.Tracker = soundfield_tracker('BinDeg', 2, 'SmoothDeg', 5, ...
                'TauHist', 0.012, 'TauSlow', 0.25, 'FloorDb', -24, ...
                'AbsFloorDb', -55, 'MinSepDeg', 16, 'MaxDet', obj.K, ...
                'AccNoise', 5e4, 'AccNoiseStable', 500, 'StableDeg', 10, ...
                'StableTime', 0.08, 'MeasDeg', 4, 'GateDeg', 22, ...
                'MaxRate', 7000, 'ConfirmTime', 0.012, 'MissTime', 0.18, ...
                'LevTau', 0.02, 'DupDeg', 25, 'DupTime', 0.08);
            for k = 1:2:numel(varargin), obj.(varargin{k}) = varargin{k+1}; end
        end

        function init(obj, fs, nIn)
            obj.fs = fs;
            obj.nIn = nIn;
            obj.K = max(1, round(obj.K));
            obj.N = max(128, round(obj.N));
            obj.hop = max(32, min(obj.N, round(obj.Hop)));
            obj.win = sqrt(0.5*(1 - cos(2*pi*(0:obj.N-1)'/obj.N)));
            [az0, ~] = separator_base.channel_layout(nIn);
            obj.chIdx = find(~isnan(az0));
            obj.az = az0(obj.chIdx);
            if nIn == 1
                obj.chIdx = 1;
                obj.az = 0;
            end
            obj.templates = obj.build_templates(obj.az);
            obj.Slots = object_slots('K', obj.K, 'MinAge', 0.03, ...
                'FadeIn', 0.01, 'FadeOut', 0.12, 'PosTau', 0.01, ...
                'PosTauStill', 0.18, 'AssignSigma', obj.AssignSigma, ...
                'AssignDeg', obj.AssignDeg, 'DupDeg', obj.Tracker.DupDeg, ...
                'StableTime', obj.Tracker.StableTime);
            bedNames = {};
            bedPos = zeros(0, 3);
            if obj.MusicBed
                bedNames = {'Music bed'};
                bedPos = separator_base.dir2pos(0);
            end
            obj.Names = [arrayfun(@(k) sprintf('Obj %d', k), 1:obj.K, 'UniformOutput', false), bedNames];
            obj.Pos = [separator_base.dir2pos(linspace(0, 360, obj.K + 1)'); bedPos];
            obj.Pos = obj.Pos([1:obj.K, obj.K + 2:end], :);
            obj.Latency = obj.N - obj.hop;
            obj.reset();
        end

        function reset(obj)
            obj.frame = zeros(obj.N, obj.nIn);
            obj.acc = zeros(obj.N, obj.nOut());
            obj.pending = zeros(0, obj.nIn);
            obj.outq = zeros(0, obj.nOut());
            obj.Slots.reset();
            obj.Tracker.MaxDet = obj.K;
            obj.Tracker.reset();
            obj.t = 0;
            obj.refEnergy = 0;
            obj.slowEnergy = 0;
            obj.lastBedAz = 0;
            obj.DiagAz = [];
            obj.DiagPA = [];
            obj.TrackLog = [];
            obj.CountLog = [];
        end

        function Y = process(obj, X)
            n = size(X, 1);
            obj.pending = [obj.pending; X];
            p = 0;
            while size(obj.pending, 1) - p >= obj.hop
                obj.frame = [obj.frame(obj.hop+1:end,:); obj.pending(p + (1:obj.hop), :)];
                p = p + obj.hop;
                obj.acc = obj.acc + obj.synth_frame();
                obj.outq = [obj.outq; obj.acc(1:obj.hop, :)];
                obj.acc = [obj.acc(obj.hop+1:end, :); zeros(obj.hop, obj.nOut())];
            end
            obj.pending = obj.pending(p+1:end, :);
            m = min(n, size(obj.outq, 1));
            Y = [zeros(n-m, obj.nOut()); obj.outq(1:m, :)];
            obj.outq = obj.outq(m+1:end, :);
        end

        function Y = separate_offline(obj, Xs)
            obj.reset();
            x0 = Xs{1};
            x = [x0; zeros(obj.Latency + obj.hop, size(x0, 2))];
            Y = zeros(size(x, 1), obj.nOut());
            for b = 1:4096:size(x, 1)
                r = b:min(size(x, 1), b+4095);
                Y(r, :) = obj.process(x(r, :));
            end
            Y = {Y(obj.Latency + (1:size(x0, 1)), :)};
        end
    end

    methods (Access = private)
        function Yf = synth_frame(obj)
            X = obj.frame(:, obj.chIdx);
            X = X - mean(X, 1);
            Xw = X .* obj.win;
            dt = obj.hop / obj.fs;
            obj.t = obj.t + dt;

            energy = mean(sum(Xw.^2, 2));
            obj.refEnergy = max(energy, obj.refEnergy*exp(-dt/2));
            [cand, ~, pattern] = obj.phantom_candidates(Xw);
            if obj.slowEnergy <= 0
                obj.slowEnergy = max(energy, eps);
            end
            onsetDb = 10*log10(max(energy, eps)/max(obj.slowEnergy, eps));
            aE = exp(-dt/max(obj.EnergyTau, eps));
            obj.slowEnergy = aE*obj.slowEnergy + (1-aE)*energy;
            sharpnessDb = obj.pattern_sharpness(pattern);
            activeFrame = energy > obj.refEnergy*10^(obj.MinEnergyDb/10) && ...
                (onsetDb >= obj.OnsetDb || sharpnessDb >= obj.SharpnessDb);
            det = zeros(0, 3);
            if activeFrame && ~isempty(cand)
                score = cand(:, 2).*cand(:, 3);
                [~, ord] = sort(score, 'descend');
                for kk = 1:numel(ord)
                    k = ord(kk);
                    if cand(k, 2) < obj.CoherentThreshold, continue; end
                    if ~isempty(det) && any(abs(obj.wrap180(det(:,1) - cand(k,1))) < obj.MinCandidateSep)
                        continue;
                    end
                    lev = 10*log10(max(cand(k, 3)*energy, eps)/max(obj.refEnergy, eps));
                    nov = 20*log10(cand(k, 2)/max(obj.CoherentThreshold, eps));
                    det(end+1, :) = [cand(k, 1), lev, nov]; %#ok<AGROW>
                    if size(det, 1) >= obj.MaxFrameDet
                        break;
                    end
                end
            end
            tr = obj.Tracker.step(det, obj.t, dt);
            obj.Slots.update(tr, obj.t, dt);

            S = zeros(obj.N, obj.nOut());
            residual = Xw;
            for k = find(obj.Slots.active)
                gk = obj.template_at(obj.Slots.az(k));
                if isempty(gk), continue; end
                sk = residual * gk(:) / (gk*gk.' + eps);
                sk = obj.Slots.g(k) * sk;
                S(:, k) = sk;
                residual = residual - sk * gk;
            end

            if obj.MusicBed
                [bedAz, ~, gBed] = obj.residual_bed(residual);
                if isfinite(bedAz)
                    obj.lastBedAz = obj.wrap180(obj.lastBedAz + ...
                        obj.wrap180(bedAz - obj.lastBedAz)*(1 - exp(-dt/max(obj.Tau, eps))));
                end
                S(:, end) = obj.AmbGain * residual * gBed(:) / (gBed*gBed.' + eps);
                obj.Pos(end, :) = separator_base.dir2pos(obj.lastBedAz);
            end

            obj.Pos(1:obj.K, :) = separator_base.dir2pos(obj.Slots.az(:));
            Yf = S .* obj.win;

            if obj.Diagnostics
                [~, ~, bi] = histcounts(obj.az_grid(), obj.AzEdges);
                ok = bi > 0;
                obj.DiagAz(:, end+1) = accumarray(bi(ok).', pattern(ok).', [numel(obj.AzEdges)-1, 1]);
                obj.DiagPA(:, end+1) = [sum(sum(S(:,1:obj.K).^2)); sum(S(:,end).^2)];
                lg = obj.Slots.az(:);
                lg(~obj.Slots.active) = NaN;
                obj.TrackLog(:, end+1) = lg;
                obj.CountLog(end+1) = nnz(obj.Slots.active);
            end
        end

        function [cand, gains, pattern] = phantom_candidates(obj, X)
            cand = zeros(0, 3);
            gains = zeros(0, size(X, 2));
            pattern = zeros(size(obj.az_grid()));
            if size(X, 1) < 8 || size(X, 2) < 1
                return;
            end

            R = (X' * X) / size(X, 1);
            R = (R + R') / 2;
            [V, D] = eig(R);
            lambda = max(real(diag(D)), 0);
            [lambda, idx] = sort(lambda, 'descend');
            V = V(:, idx);
            nSrc = obj.estimate_count(lambda);
            if nSrc <= 0
                return;
            end

            chPower = max(real(diag(R)), 0);
            C = abs(R) ./ sqrt(max(chPower * chPower.', eps));
            C(1:size(C,1)+1:end) = 0;
            if nSrc == 1
                rows = V(:, 1).';
            elseif max(C(:)) < obj.IcaThreshold
                [~, chOrder] = sort(chPower, 'descend');
                rows = eye(size(X, 2));
                rows = rows(chOrder(1:min(nSrc, size(rows, 1))), :);
            else
                rows = obj.ica_rows(X, V(:, 1:nSrc), lambda(1:nSrc));
                if isempty(rows), rows = V(:, 1:nSrc).'; end
            end

            weights = sum(rows.^2, 2);
            weights = weights ./ (sum(weights) + eps);
            usedAz = [];
            for c = 1:size(rows, 1)
                [azCand, conf, g, pat] = obj.match_vector(rows(c, :));
                pattern = max(pattern, weights(c) * pat);
                if conf <= 0 || any(abs(obj.wrap180(usedAz - azCand)) < 8)
                    continue;
                end
                usedAz(end+1) = azCand; %#ok<AGROW>
                cand(end+1, :) = [azCand, conf, weights(c)]; %#ok<AGROW>
                gains(end+1, :) = g; %#ok<AGROW>
            end
        end

        function db = pattern_sharpness(~, pattern)
            p = pattern(:);
            p = p(isfinite(p) & p > 0);
            if numel(p) < 4
                db = 0;
                return;
            end
            db = 20*log10(max(p)/max(prctile(p, 80), eps));
        end

        function [az, conf, g, pat] = match_vector(obj, v)
            v = v(:).';
            v = v ./ (norm(v) + eps);
            azs = obj.az_grid();
            pat = zeros(size(azs));
            bestG = [];
            for k = 1:numel(azs)
                G = obj.templates{k};
                if isempty(G), continue; end
                s = abs(G * v.');
                [pat(k), ii] = max(s);
                if isempty(bestG) || pat(k) > conf_or_zero(pat, k-1)
                    bestG = G(ii, :);
                end
            end
            [conf, im] = max(pat);
            az = azs(im);
            g = bestG;
            if isempty(g), g = obj.template_at(az); end

            function c = conf_or_zero(p, kmax)
                if kmax < 1, c = 0; else, c = max(p(1:kmax)); end
            end
        end

        function [az, conf, g] = residual_bed(obj, X)
            [cand, gains] = obj.phantom_candidates(X);
            if isempty(cand)
                az = obj.lastBedAz;
                conf = 0;
                g = obj.template_at(az);
                return;
            end
            [conf, im] = max(cand(:, 2).*cand(:, 3));
            az = cand(im, 1);
            g = gains(im, :);
        end

        function nSrc = estimate_count(obj, lambda)
            if isempty(lambda) || lambda(1) <= 0
                nSrc = 0;
                return;
            end
            maxSrc = min([obj.MaxComponents, numel(lambda)]);
            tail = lambda(max(2, numel(lambda)-1):end);
            noise = max(median(tail), lambda(1)*1e-4);
            nSrc = sum(lambda(1:maxSrc) > 4*noise);
            nSrc = min(max(nSrc, 1), maxSrc);
        end

        function rows = ica_rows(~, X, E, lambda)
            rows = [];
            q = size(E, 2);
            Z = X * E * diag(1 ./ sqrt(lambda(:).' + eps));
            Z = Z - mean(Z, 1);
            W = zeros(q, q);
            for c = 1:q
                w = zeros(q, 1);
                w(c) = 1;
                for it = 1:80
                    y = Z * w;
                    gy = tanh(y);
                    gp = 1 - gy.^2;
                    wNew = (Z' * gy) / size(Z, 1) - mean(gp) * w;
                    if c > 1
                        wNew = wNew - W(1:c-1,:).' * (W(1:c-1,:) * wNew);
                    end
                    nw = norm(wNew);
                    if nw <= eps, break; end
                    wNew = wNew ./ nw;
                    if abs(abs(wNew' * w) - 1) < 1e-5
                        w = wNew;
                        break;
                    end
                    w = wNew;
                end
                W(c,:) = w.';
            end
            for c = 1:q
                y = Z * W(c,:).';
                den = y' * y;
                if den <= eps, continue; end
                rows(end+1, :) = (y' * X) / den; %#ok<AGROW>
            end
        end

        function g = template_at(obj, az)
            [~, im] = min(abs(obj.wrap180(obj.az_grid() - az)));
            G = obj.templates{im};
            if isempty(G)
                g = zeros(1, numel(obj.az));
                g(1) = 1;
            else
                g = G(1, :);
            end
        end

        function templates = build_templates(obj, az)
            grid = obj.az_grid();
            templates = cell(size(grid));
            pairs = obj.phantom_pairs(az);
            nCh = numel(az);
            for pidx = 1:size(pairs, 1)
                i = pairs(pidx, 1);
                j = pairs(pidx, 2);
                d = obj.wrap180(az(j) - az(i));
                if abs(d) < 1e-9, continue; end
                for k = 1:numel(grid)
                    dt = obj.wrap180(grid(k) - az(i));
                    if d > 0
                        inside = dt >= -1e-9 && dt <= d + 1e-9;
                    else
                        inside = dt <= 1e-9 && dt >= d - 1e-9;
                    end
                    if ~inside, continue; end
                    u = min(1, max(0, dt/d));
                    g = zeros(1, nCh);
                    g(i) = cos(u*pi/2);
                    g(j) = sin(u*pi/2);
                    g = g ./ (norm(g) + eps);
                    if isempty(templates{k}), templates{k} = g; else, templates{k}(end+1,:) = g; end
                end
            end
            for k = 1:numel(grid)
                for i = 1:nCh
                    if abs(obj.wrap180(grid(k) - az(i))) <= 1
                        g = zeros(1, nCh);
                        g(i) = 1;
                        if isempty(templates{k}), templates{k} = g; else, templates{k}(end+1,:) = g; end
                    end
                end
            end
        end

        function pairs = phantom_pairs(~, az)
            nCh = numel(az);
            pairs = zeros(0, 2);
            for i = 1:nCh
                for j = i+1:nCh
                    d = phantom_object_separator.wrap180(az(j) - az(i));
                    if abs(d) <= 150
                        pairs(end+1,:) = [i, j]; %#ok<AGROW>
                    end
                end
            end
        end
    end

    methods (Static)
        function a = wrap180(a)
            a = mod(a + 180, 360) - 180;
        end

        function az = az_grid()
            az = -180:2:180;
        end
    end
end
