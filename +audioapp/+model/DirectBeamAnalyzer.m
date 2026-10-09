classdef DirectBeamAnalyzer < handle
    %DIRECTBEAMANALYZER  Lightweight horizontal beam energy view.
    % Encodes known channel/source directions into a 2-D circular-harmonic
    % beamspace and scans it around the listener for display and later
    % tracking.

    properties
        Order double = 2
        ScanAnglesDeg double = -180:2:180
        SmoothPattern double = []
        SmoothMusicPattern double = []
        SmoothPhantomPattern double = []
        SmoothEnergy double = 0
        MusicMaxSources double = 3
        LastSourceNames cell = {}
        LastSourceAnglesDeg double = []
        LastSelectedChannels double = []
    end

    methods
        function obj = DirectBeamAnalyzer(order)
            if nargin >= 1 && ~isempty(order)
                obj.Order = order;
            end
        end

        function reset(obj)
            obj.SmoothPattern = [];
            obj.SmoothMusicPattern = [];
            obj.SmoothPhantomPattern = [];
            obj.SmoothEnergy = 0;
            obj.LastSourceNames = {};
            obj.LastSourceAnglesDeg = [];
            obj.LastSelectedChannels = [];
        end

        function data = update(obj, X, sourceNames)
            data = struct('active', false, 'anglesDeg', obj.ScanAnglesDeg, ...
                'pattern', zeros(size(obj.ScanAnglesDeg)), ...
                'musicPattern', zeros(size(obj.ScanAnglesDeg)), ...
                'musicSourceCount', 0, ...
                'phantomPattern', zeros(size(obj.ScanAnglesDeg)), ...
                'phantomSourceCount', 0, ...
                'sourceAnglesDeg', [], 'sourceNames', {{}}, 'selectedChannels', []);
            if nargin < 3
                sourceNames = {};
            end
            if isempty(X) || size(X, 2) < 2
                obj.decay();
                data.pattern = obj.normalizedPattern();
                data.musicPattern = obj.normalizedMusicPattern();
                data.phantomPattern = obj.normalizedPhantomPattern();
                return;
            end

            [ch, az, names] = obj.selectDirectionalChannels(size(X,2), sourceNames);
            if numel(ch) < 2
                obj.decay();
                data.pattern = obj.normalizedPattern();
                data.musicPattern = obj.normalizedMusicPattern();
                data.phantomPattern = obj.normalizedPhantomPattern();
                return;
            end

            order = min(obj.Order, floor((numel(ch)-1)/2));
            order = max(order, 1);
            X = double(X(:, ch));
            X = X - mean(X, 1);

            Bsrc = audioapp.model.DirectBeamAnalyzer.circularBasis(az, order);
            Bscan = audioapp.model.DirectBeamAnalyzer.circularBasis(obj.ScanAnglesDeg, order);

            hoa = X * Bsrc;
            Y = hoa * Bscan.';
            p = mean(abs(Y).^2, 1);
            if ~all(isfinite(p)) || max(p) <= 0
                obj.decay();
                data.pattern = obj.normalizedPattern();
                data.musicPattern = obj.normalizedMusicPattern();
                data.phantomPattern = obj.normalizedPhantomPattern();
                return;
            end

            [pMusic, nSrc] = obj.musicPattern(hoa, Bscan);
            [pPhantom, nPhantom] = obj.phantomSubspacePattern(X, az);

            obj.smooth(p);
            obj.smoothMusic(pMusic);
            obj.smoothPhantom(pPhantom);
            data.active = max(obj.SmoothPattern) > 1e-12 || max(obj.SmoothMusicPattern) > 1e-12 || ...
                max(obj.SmoothPhantomPattern) > 1e-12;
            data.pattern = obj.normalizedPattern();
            data.musicPattern = obj.normalizedMusicPattern();
            data.musicSourceCount = nSrc;
            data.phantomPattern = obj.normalizedPhantomPattern();
            data.phantomSourceCount = nPhantom;
            data.sourceAnglesDeg = az;
            data.sourceNames = names;
            data.selectedChannels = ch;
            obj.LastSourceAnglesDeg = az;
            obj.LastSourceNames = names;
            obj.LastSelectedChannels = ch;
        end
    end

    methods (Access = private)
        function smooth(obj, p)
            p = p(:).';
            if isempty(obj.SmoothPattern) || numel(obj.SmoothPattern) ~= numel(p)
                obj.SmoothPattern = p;
                return;
            end
            rising = p > obj.SmoothPattern;
            a = 0.08 * ones(size(p));   % slow release
            a(rising) = 0.45;           % faster attack
            obj.SmoothPattern = (1-a).*obj.SmoothPattern + a.*p;
        end

        function decay(obj)
            if ~isempty(obj.SmoothPattern)
                obj.SmoothPattern = 0.90 * obj.SmoothPattern;
            end
            if ~isempty(obj.SmoothMusicPattern)
                obj.SmoothMusicPattern = 0.85 * obj.SmoothMusicPattern;
            end
            if ~isempty(obj.SmoothPhantomPattern)
                obj.SmoothPhantomPattern = 0.82 * obj.SmoothPhantomPattern;
            end
        end

        function p = normalizedPattern(obj)
            if isempty(obj.SmoothPattern)
                p = zeros(size(obj.ScanAnglesDeg));
                return;
            end
            p = obj.SmoothPattern;
            p = p ./ (max(p) + eps);
            p(p < 0.02) = 0;
        end

        function smoothMusic(obj, p)
            p = p(:).';
            if isempty(p) || ~all(isfinite(p))
                obj.decay();
                return;
            end
            if isempty(obj.SmoothMusicPattern) || numel(obj.SmoothMusicPattern) ~= numel(p)
                obj.SmoothMusicPattern = p;
                return;
            end
            rising = p > obj.SmoothMusicPattern;
            a = 0.12 * ones(size(p));
            a(rising) = 0.55;
            obj.SmoothMusicPattern = (1-a).*obj.SmoothMusicPattern + a.*p;
        end

        function p = normalizedMusicPattern(obj)
            if isempty(obj.SmoothMusicPattern)
                p = zeros(size(obj.ScanAnglesDeg));
                return;
            end
            p = obj.SmoothMusicPattern;
            p = p ./ (max(p) + eps);
            p(p < 0.08) = 0;
        end

        function smoothPhantom(obj, p)
            p = p(:).';
            if isempty(p) || ~all(isfinite(p))
                obj.decay();
                return;
            end
            if isempty(obj.SmoothPhantomPattern) || numel(obj.SmoothPhantomPattern) ~= numel(p)
                obj.SmoothPhantomPattern = p;
                return;
            end
            rising = p > obj.SmoothPhantomPattern;
            a = 0.10 * ones(size(p));
            a(rising) = 0.50;
            obj.SmoothPhantomPattern = (1-a).*obj.SmoothPhantomPattern + a.*p;
        end

        function p = normalizedPhantomPattern(obj)
            if isempty(obj.SmoothPhantomPattern)
                p = zeros(size(obj.ScanAnglesDeg));
                return;
            end
            p = obj.SmoothPhantomPattern;
            p = p ./ (max(p) + eps);
            p(p < 0.10) = 0;
        end

        function [p, nSrc] = phantomSubspacePattern(obj, X, az)
            p = zeros(size(obj.ScanAnglesDeg));
            nSrc = 0;
            if size(X, 1) < 8 || size(X, 2) < 2
                return;
            end

            R = (X' * X) / size(X, 1);
            R = (R + R') / 2;
            [V, D] = eig(R);
            lambda = max(real(diag(D)), 0);
            [lambda, idx] = sort(lambda, 'descend');
            V = V(:, idx);
            nSrc = obj.estimateSignalCount(lambda);
            if nSrc <= 0
                return;
            end
            nSrc = min(nSrc, size(V, 2));
            templates = obj.phantomGainTemplates(az);
            if nSrc == 1
                p = obj.phantomVectorPattern(V(:, 1).', templates);
                return;
            end

            chPower = max(real(diag(R)), 0);
            C = abs(R) ./ sqrt(max(chPower * chPower.', eps));
            C(1:size(C,1)+1:end) = 0;
            if max(C(:)) < 0.35
                [~, chIdx] = sort(chPower, 'descend');
                chIdx = chIdx(1:min(nSrc, numel(chIdx)));
                mixRows = eye(size(X, 2));
                mixRows = mixRows(chIdx, :);
            else
                mixRows = obj.icaMixingRows(X, V(:, 1:nSrc), lambda(1:nSrc));
            end
            if isempty(mixRows)
                mixRows = V(:, 1:nSrc).';
            end
            weights = sum(mixRows.^2, 2);
            weights = weights ./ (sum(weights) + eps);
            for c = 1:size(mixRows, 1)
                p = p + weights(c) * obj.phantomVectorPattern(mixRows(c, :), templates);
            end
        end

        function p = phantomVectorPattern(obj, v, templates)
            p = zeros(size(obj.ScanAnglesDeg));
            v = v(:).';
            v = v ./ (norm(v) + eps);
            for k = 1:numel(obj.ScanAnglesDeg)
                G = templates{k};
                if isempty(G)
                    continue;
                end
                p(k) = max(abs(G * v.'));
            end
            p = p.^6;
        end

        function mixRows = icaMixingRows(~, X, E, lambda)
            mixRows = [];
            q = size(E, 2);
            if q < 1
                return;
            end
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
                    if nw <= eps
                        break;
                    end
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
                if den <= eps
                    continue;
                end
                mixRows(end+1, :) = (y' * X) / den; %#ok<AGROW>
            end
        end

        function templates = phantomGainTemplates(obj, az)
            az = az(:).';
            scan = obj.ScanAnglesDeg(:).';
            nCh = numel(az);
            templates = cell(size(scan));
            if nCh < 2
                return;
            end

            pairs = obj.phantomSpeakerPairs(az);
            for pidx = 1:size(pairs, 1)
                i = pairs(pidx, 1);
                j = pairs(pidx, 2);
                    d = audioapp.model.DirectBeamAnalyzer.wrap180(az(j) - az(i));
                    if abs(d) < 1e-9
                        continue;
                    end
                    for k = 1:numel(scan)
                        dt = audioapp.model.DirectBeamAnalyzer.wrap180(scan(k) - az(i));
                        if d > 0
                            inside = dt >= -1e-9 && dt <= d + 1e-9;
                        else
                            inside = dt <= 1e-9 && dt >= d - 1e-9;
                        end
                        if ~inside
                            continue;
                        end
                        u = min(1, max(0, dt / d));
                        g = zeros(1, nCh);
                        g(i) = cos(u*pi/2);
                        g(j) = sin(u*pi/2);
                        g = g ./ (norm(g) + eps);
                        if isempty(templates{k})
                            templates{k} = g;
                        else
                            templates{k}(end+1, :) = g;
                        end
                    end
            end

            for k = 1:numel(scan)
                for i = 1:nCh
                    if abs(audioapp.model.DirectBeamAnalyzer.wrap180(scan(k) - az(i))) <= 1
                        g = zeros(1, nCh);
                        g(i) = 1;
                        if isempty(templates{k})
                            templates{k} = g;
                        else
                            templates{k}(end+1, :) = g;
                        end
                    end
                end
            end
        end

        function pairs = phantomSpeakerPairs(~, az)
            [azSorted, ord] = sort(az(:).');
            nCh = numel(azSorted);
            pairs = zeros(0, 2);
            for k = 1:nCh-1
                pairs(end+1, :) = [ord(k), ord(k+1)]; %#ok<AGROW>
            end
            pairs(end+1, :) = [ord(end), ord(1)];

            front = find(abs(az) <= 90);
            for a = 1:numel(front)
                for b = a+1:numel(front)
                    i = front(a);
                    j = front(b);
                    if sign(az(i)) ~= sign(az(j))
                        pairs(end+1, :) = [i, j]; %#ok<AGROW>
                    end
                end
            end
            pairs = unique(sort(pairs, 2), 'rows', 'stable');
        end

        function [p, nSrc] = musicPattern(obj, hoa, Bscan)
            m = size(hoa, 2);
            p = zeros(1, size(Bscan, 1));
            nSrc = 0;
            if size(hoa, 1) < 8 || m < 3
                return;
            end

            R = (hoa' * hoa) / size(hoa, 1);
            R = (R + R') / 2;
            [V, D] = eig(R);
            lambda = max(real(diag(D)), 0);
            [lambda, idx] = sort(lambda, 'descend');
            V = V(:, idx);
            nSrc = obj.estimateSignalCount(lambda);
            if nSrc <= 0 || nSrc >= m
                return;
            end

            Un = V(:, nSrc+1:end);
            A = Bscan;
            A = A ./ max(vecnorm(A, 2, 2), eps);
            for k = 1:size(A, 1)
                a = A(k,:).';
                den = real(a' * (Un * Un') * a);
                p(k) = 1 / max(den, eps);
            end
        end

        function nSrc = estimateSignalCount(obj, lambda)
            if isempty(lambda) || lambda(1) <= 0
                nSrc = 0;
                return;
            end
            maxSrc = min([obj.MusicMaxSources, numel(lambda)-2, obj.Order+1]);
            if maxSrc < 1
                nSrc = 0;
                return;
            end
            tail = lambda(max(2, numel(lambda)-1):end);
            noise = max(median(tail), lambda(1)*1e-4);
            nSrc = sum(lambda(1:maxSrc) > 4*noise);
            nSrc = min(max(nSrc, 1), maxSrc);
        end

        function [ch, az, names] = selectDirectionalChannels(~, nCh, sourceNames)
            if isstring(sourceNames), sourceNames = cellstr(sourceNames); end
            if ischar(sourceNames), sourceNames = {sourceNames}; end
            sourceNames = sourceNames(:).';
            if numel(sourceNames) < nCh
                sourceNames(end+1:nCh) = {''};
            end

            ch = [];
            az = [];
            names = {};
            for k = 1:nCh
                [ok, a, nm] = audioapp.model.DirectBeamAnalyzer.nameToAzimuth(sourceNames{k});
                if ok
                    ch(end+1) = k; %#ok<AGROW>
                    az(end+1) = a; %#ok<AGROW>
                    names{end+1} = nm; %#ok<AGROW>
                end
            end

            if isempty(ch)
                [layoutAz, layoutNames] = separator_base.channel_layout(nCh);
                keep = isfinite(layoutAz);
                ch = find(keep);
                az = layoutAz(keep);
                names = layoutNames(keep);
            end
        end
    end

    methods (Static)
        function B = circularBasis(azDeg, order)
            az = deg2rad(azDeg(:));
            B = ones(numel(az), 1);
            for n = 1:order
                B = [B, cos(n*az), sin(n*az)]; %#ok<AGROW>
            end
            colNorm = sqrt(sum(B.^2, 1));
            colNorm(colNorm == 0) = 1;
            B = B ./ colNorm;
        end

        function a = wrap180(a)
            a = mod(a + 180, 360) - 180;
        end

        function [ok, az, nm] = nameToAzimuth(name)
            ok = false;
            az = NaN;
            nm = '';
            if isempty(name), return; end
            name = strtrim(char(name));
            if startsWith(name, 'Direct ', 'IgnoreCase', true)
                name = strtrim(name(8:end));
            elseif startsWith(name, 'Dir ', 'IgnoreCase', true)
                return;
            elseif startsWith(name, 'Amb ', 'IgnoreCase', true)
                return;
            end

            switch upper(name)
                case 'L'
                    az = 30; ok = true; nm = 'L';
                case 'R'
                    az = -30; ok = true; nm = 'R';
                case 'C'
                    az = 0; ok = true; nm = 'C';
                case {'LS', 'SL'}
                    az = 110; ok = true; nm = 'Ls';
                case {'RS', 'SR'}
                    az = -110; ok = true; nm = 'Rs';
                case {'LB', 'BL'}
                    az = 150; ok = true; nm = 'Lb';
                case {'RB', 'BR'}
                    az = -150; ok = true; nm = 'Rb';
                otherwise
                    ok = false;
            end
        end
    end
end
