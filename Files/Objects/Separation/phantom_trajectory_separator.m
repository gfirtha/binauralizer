classdef phantom_trajectory_separator < separator_base
    %PHANTOM_TRAJECTORY_SEPARATOR  Fast phantom-pan trajectory extractor.
    %
    %   Short-window tracker for transient objects whose apparent position is
    %   encoded by coherent loudspeaker gains. The object signal is projected
    %   from the current gain vector; the reconstruction is subtracted and the
    %   residual is emitted as a music bed at its current dominant phantom
    %   direction.

    properties
        K = 3
        N = 256
        Hop = 64
        MusicBed = true
        AmbGain = 1
        MinEnergyDb = -45
        CoherentThreshold = 0.18
        AssignSigma = 18
        AssignDeg = 55
        Width = 180
        Base = 30
        AmbAz = 110
        Tau = 0.02
        Tracker
        Diagnostics = false
    end

    properties (SetAccess = private)
        DiagAz = []
        DiagPA = []
        AzEdges = -180:4:180
        TrackLog = []
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
        lastBedAz = 0
    end

    methods
        function obj = phantom_trajectory_separator(varargin)
            obj.Tracker = soundfield_tracker('BinDeg', 2, 'SmoothDeg', 4, ...
                'TauHist', 0.004, 'TauSlow', 0.12, 'FloorDb', -20, ...
                'AbsFloorDb', -50, 'MinSepDeg', 18, 'MaxDet', obj.K, ...
                'AccNoise', 8e4, 'AccNoiseStable', 1e4, 'StableDeg', 10, ...
                'StableTime', 0.04, 'MeasDeg', 4, 'GateDeg', 25, ...
                'MaxRate', 9000, 'ConfirmTime', 0.004, 'MissTime', 0.08, ...
                'LevTau', 0.01, 'DupDeg', 14, 'DupTime', 0.04);
            for k = 1:2:numel(varargin), obj.(varargin{k}) = varargin{k+1}; end
        end

        function init(obj, fs, nIn)
            obj.fs = fs;
            obj.nIn = nIn;
            obj.K = max(1, round(obj.K));
            obj.N = max(64, round(obj.N));
            obj.hop = max(16, min(obj.N, round(obj.Hop)));
            obj.win = sqrt(0.5*(1 - cos(2*pi*(0:obj.N-1)'/obj.N)));
            [az0, ~] = separator_base.channel_layout(nIn);
            obj.chIdx = find(~isnan(az0));
            obj.az = az0(obj.chIdx);
            if nIn == 1
                obj.az = 0;
            end
            obj.templates = obj.build_templates(obj.az);
            obj.Slots = object_slots('K', obj.K, 'MinAge', 0, 'FadeIn', 0.004, ...
                'FadeOut', 0.08, 'PosTau', 0.004, 'PosTauStill', 0.08, ...
                'AssignSigma', obj.AssignSigma, 'AssignDeg', obj.AssignDeg, ...
                'DupDeg', obj.Tracker.DupDeg, 'StableTime', obj.Tracker.StableTime);
            bedNames = {};
            bedPos = zeros(0, 3);
            if obj.MusicBed
                bedNames = {'Music bed'};
                bedPos = separator_base.dir2pos(0);
            end
            obj.Names = [arrayfun(@(k) sprintf('Traj %d', k), 1:obj.K, 'UniformOutput', false), bedNames];
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
            obj.lastBedAz = 0;
            obj.DiagAz = [];
            obj.DiagPA = [];
            obj.TrackLog = [];
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

            [azPeak, conf, gPeak, pattern] = obj.phantom_peak(Xw);
            energy = mean(sum(Xw.^2, 2));
            obj.refEnergy = max(energy, obj.refEnergy*exp(-dt/2));
            active = energy > obj.refEnergy*10^(obj.MinEnergyDb/10) && conf > obj.CoherentThreshold;
            if active
                det = [azPeak, 10*log10(energy/max(obj.refEnergy, eps)), 20*log10(conf + eps)];
            else
                det = zeros(0, 3);
            end
            tr = obj.Tracker.step(det, obj.t, dt);
            obj.Slots.update(tr, obj.t, dt);

            nOut = obj.nOut();
            S = zeros(obj.N, nOut);
            residual = Xw;
            for k = find(obj.Slots.active)
                azk = obj.Slots.az(k);
                gk = obj.template_at(azk);
                if isempty(gk)
                    gk = gPeak;
                end
                sk = Xw * gk(:) / (gk*gk.' + eps);
                sk = obj.Slots.g(k) * sk;
                S(:, k) = sk;
                residual = residual - sk * gk;
            end
            if obj.MusicBed
                [bedAz, ~, gBed] = obj.phantom_peak(residual);
                if isfinite(bedAz)
                    obj.lastBedAz = obj.wrap180(obj.lastBedAz + obj.wrap180(bedAz - obj.lastBedAz)*(1 - exp(-dt/max(obj.Tau, eps))));
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
            end
        end

        function [azPeak, conf, gPeak, pattern] = phantom_peak(obj, X)
            nA = numel(obj.templates);
            pattern = zeros(1, nA);
            R = (X' * X) / max(size(X, 1), 1);
            R = (R + R') / 2;
            [V, D] = eig(R);
            lambda = max(real(diag(D)), 0);
            [lambda, idx] = sort(lambda, 'descend');
            V = V(:, idx);
            if isempty(lambda) || lambda(1) <= 0
                azPeak = 0;
                conf = 0;
                gPeak = obj.template_at(0);
                return;
            end
            v = V(:, 1).';
            v = v ./ (norm(v) + eps);
            for k = 1:nA
                G = obj.templates{k};
                if isempty(G), continue; end
                pattern(k) = max(abs(G * v.'));
            end
            [conf, im] = max(pattern);
            azs = obj.az_grid();
            azPeak = azs(im);
            gPeak = obj.best_template(im, v);
        end

        function g = template_at(obj, az)
            [~, im] = min(abs(obj.wrap180(obj.az_grid() - az)));
            g = obj.best_template(im, []);
        end

        function g = best_template(obj, im, v)
            G = obj.templates{im};
            if isempty(G)
                g = zeros(1, numel(obj.az));
                g(1) = 1;
                return;
            end
            if isempty(v)
                g = G(1, :);
                return;
            end
            [~, ii] = max(abs(G * v(:)));
            g = G(ii, :);
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
                    d = phantom_trajectory_separator.wrap180(az(j) - az(i));
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
