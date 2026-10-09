classdef soundfield_tracker < handle
    %SOUNDFIELD_TRACKER  Multi-object direction tracker on the output of
    %soundfield_analyzer (preferably Mode 'sectors': several simultaneous
    %directions per band).
    %
    %   Detection, per analysis frame (variable length: block switching):
    %   circular histogram of the direct power E (1 - psi) of the bands in
    %   Band over azimuth (BinDeg bins, Gaussian smoothing SmoothDeg),
    %   recursively averaged in time (TauHist, frame-length aware); peaks
    %   above the strongest - FloorDb and above the excerpt's AbsFloorDb,
    %   at least MinSepDeg apart (parabolic refinement) are detections
    %   (azimuth, level).
    %   Tracking: one Kalman filter per object, state [azimuth; angular
    %   velocity] (deg, deg/s), constant-velocity model with white
    %   angular-acceleration noise AccNoise (a fly-by can sweep 180 deg in
    %   0.1 s) - except for tracks that have stood still for StableTime:
    %   they get AccNoiseStable (a two-model, IMM-like switch), so a
    %   stationary source keeps its track and a fly-by starting from its
    %   direction opens a new one instead of dragging it away (identity
    %   swap); measurement noise grows for weak detections. Global nearest
    %   neighbour association within a Mahalanobis gate (GateSigma) and
    %   GateDeg + MaxRate x frame step (short frames: small steps, so a
    %   track cannot jump between two separate sources and fake a sweep);
    %   unassociated detections start tentative tracks, confirmed
    %   after ConfirmTime of (nearly) continuous hits, deleted after
    %   MissTime without hits (position coasts on the velocity, decaying).
    %
    %   Gating and births: sector groups more diffuse than PsiMax give no
    %   detections (ambience); every detection carries its novelty, the
    %   level above a slow (TauSlow) background histogram in the same
    %   direction. A novel detection (> OnsetDb) may not continue a
    %   standing track - a fly-by entering the Ls direction opens its own
    %   object instead of disappearing in the Ls ambience track. Tracks
    %   running within DupDeg of an older one for DupTime are dropped.
    %
    %   Streaming: reset(), then per frame det = detect_frame(az, pw, dt)
    %   and tr = step(det, t, dt) (live tracks). Offline:
    %   T = track(A) -> struct array, one per confirmed track: id, t
    %   (frame times [s]), az, vel, lev (smoothed level, dB re excerpt max
    %   of the histogram), hit (detection associated in that frame), born
    %   (start time), novel (born from a novel detection).

    properties
        Band = [500 Inf]    % [Hz] bands used for detection
        BinDeg = 2
        SmoothDeg = 6       % Gaussian sigma of the azimuth smoothing
        TauHist = 0.01      % [s] time smoothing of the histogram
        FloorDb = -30       % peaks within this of the strongest one
        AbsFloorDb = -60    % ... and of the excerpt's strongest frame
        MinSepDeg = 25
        MaxDet = 6
        AccNoise = 2e4      % [deg/s^2] std of the angular acceleration
        AccNoiseStable = 300 % ... of a track that has stood still for StableTime
        StableDeg = 6       % a track staying within this of an anchor counts as standing
        StableTime = 0.2    % [s]
        MeasDeg = 3         % measurement std of a strong detection [deg]
        MeasDegPerDb = 0.25 % ... + this per dB below the strongest
        GateSigma = 3
        GateDeg = 8         % gate: GateDeg + MaxRate x frame step (no jumps between
        MaxRate = 1500      % [deg/s]  separate sources - a sweep needs detections along its way)
        ConfirmTime = 0.015 % [s]
        MissTime = 0.15     % [s]
        VelDecay = 0.05     % [s] velocity decay while coasting
        LevTau = 0.02       % [s] level smoothing
        PsiMax = 0.5        % sector groups more diffuse than this give no detections
        TauSlow = 0.3       % [s] background histogram (novelty reference)
        OnsetDb = 6         % detection this far above the background in its direction is
                            % 'new': it may not continue a standing track, it opens its own
        NovelProtect = 1    % [s] the novelty rule applies to standing tracks older than this
                            % (a source's own onset stays novel for ~TauSlow)
        DupDeg = 10         % two tracks this close ...
        DupTime = 0.1       % [s] ... for this long: the younger one is dropped (a novel
                            % track only after 3 DupTime)
        RefDecay = 5        % [s] streaming: decay of the running level reference
        KeepLog = false     % keep per-frame logs (offline track() sets it)
    end

    properties (Access = private)
        tr; done; nextId = 1
        hF = []; hS = []; gref = 0      % streaming detector state
    end

    methods
        function obj = soundfield_tracker(varargin)
            for k = 1:2:numel(varargin), obj.(varargin{k}) = varargin{k+1}; end
        end

        function reset(obj)
            % Clear all tracks (and the streaming detector state)
            obj.tr = struct('id', {}, 'x', {}, 'P', {}, 'lev', {}, 'age', {}, 'hitTime', {}, ...
                'lastHit', {}, 'conf', {}, 'still', {}, 'anchor', {}, 'novel', {}, 'dup', {}, 'log', {});
            obj.done = obj.tr; obj.nextId = 1;
            obj.hF = []; obj.hS = []; obj.gref = 0;
        end

        function tr = step(obj, det, t, dt)
            % One frame: det = [az, level dB, novelty dB] (rows), frame
            % time t and step dt [s]. Returns the live tracks (fields id,
            % x = [az; vel], lev, conf, novel, lastHit, ...).
            tr = obj.tr;
            % ---- predict ----
            Fm = [1 dt; 0 1]; Q0 = [dt^3/3 dt^2/2; dt^2/2 dt];
            for k = 1:numel(tr)
                if abs(wrap(tr(k).x(1) - tr(k).anchor)) <= obj.StableDeg
                    tr(k).still = tr(k).still + dt;
                else
                    tr(k).anchor = tr(k).x(1); tr(k).still = 0;
                end
                if tr(k).still >= obj.StableTime, q = obj.AccNoiseStable^2; else, q = obj.AccNoise^2; end
                coast = t - tr(k).lastHit > 1.5*dt;
                if coast, tr(k).x(2) = tr(k).x(2)*exp(-dt/obj.VelDecay); end
                tr(k).x = Fm*tr(k).x; tr(k).x(1) = wrap(tr(k).x(1));
                tr(k).P = Fm*tr(k).P*Fm' + q*Q0;
            end
            % ---- associate (global nearest neighbour, gated) ----
            nd = size(det, 1);
            used = false(nd, 1); hit = false(1, numel(tr));
            if nd > 0 && ~isempty(tr)
                C = inf(numel(tr), nd);
                for k = 1:numel(tr)
                    for j = 1:nd
                        r = obj.meas_var(det(j,2));
                        v = wrap(det(j,1) - tr(k).x(1)); S = tr(k).P(1,1) + r;
                        if abs(v) <= obj.GateDeg + obj.MaxRate*dt && v^2/S <= obj.GateSigma^2
                            C(k,j) = v^2/S + log(S);
                        end
                        if det(j,3) > obj.OnsetDb && tr(k).still >= obj.StableTime && t - tr(k).age > obj.NovelProtect
                            C(k,j) = inf;                   % novel: not for an established standing track
                        end
                    end
                end
                while true
                    [cmin, ix] = min(C(:));
                    if ~isfinite(cmin), break; end
                    [k, j] = ind2sub(size(C), ix);
                    r = obj.meas_var(det(j,2));
                    v = wrap(det(j,1) - tr(k).x(1)); S = tr(k).P(1,1) + r;
                    Kg = tr(k).P(:,1)/S;
                    tr(k).x = tr(k).x + Kg*v; tr(k).x(1) = wrap(tr(k).x(1));
                    tr(k).P = tr(k).P - Kg*tr(k).P(1,:);
                    al = 1 - exp(-dt/obj.LevTau);
                    tr(k).lev = tr(k).lev + al*(det(j,2) - tr(k).lev);
                    tr(k).hitTime = tr(k).hitTime + dt; tr(k).lastHit = t;
                    hit(k) = true; used(j) = true;
                    C(k,:) = inf; C(:,j) = inf;
                end
            end
            % ---- manage ----
            for k = 1:numel(tr)
                if ~hit(k) && ~tr(k).conf, tr(k).hitTime = 0; end   % tentative: needs continuity
                if ~tr(k).conf && tr(k).hitTime >= obj.ConfirmTime, tr(k).conf = true; end
                if obj.KeepLog, tr(k).log(end+1,:) = [t, tr(k).x(1), tr(k).x(2), tr(k).lev, hit(k)]; end
            end
            if ~isempty(tr)
                azs = arrayfun(@(s) s.x(1), tr); ids = [tr.id];
                for k = 1:numel(tr)
                    near = ids < tr(k).id & abs(wrap(azs - tr(k).x(1))) <= obj.DupDeg;
                    young = tr(k).novel && t - tr(k).age < 3*obj.DupTime;   % a new object keeps its own track first
                    if any(near) && ~young, tr(k).dup = tr(k).dup + dt; else, tr(k).dup = 0; end
                end
                dead = arrayfun(@(s) t - s.lastHit > obj.MissTime || (~s.conf && t - s.lastHit > 2*obj.ConfirmTime) ...
                    || s.dup > obj.DupTime, tr);
                if obj.KeepLog, obj.done = [obj.done, tr(dead & [tr.conf])]; end
                tr = tr(~dead);
            end
            for j = find(~used)'
                r = obj.meas_var(det(j,2));
                lg = zeros(0, 5); if obj.KeepLog, lg = [t, det(j,1), 0, det(j,2), 1]; end
                tr(end+1) = struct('id', obj.nextId, 'x', [det(j,1); 0], 'P', diag([r, 500^2]), ...
                    'lev', det(j,2), 'age', t, 'hitTime', dt, 'lastHit', t, 'conf', false, 'still', 0, ...
                    'anchor', det(j,1), 'novel', det(j,3) > obj.OnsetDb, 'dup', 0, 'log', lg); %#ok<AGROW>
                obj.nextId = obj.nextId + 1;
            end
            obj.tr = tr;

            function a = wrap(a), a = mod(a + 180, 360) - 180; end
        end

        function det = detect_frame(obj, az, pw, dt)
            % Streaming detection: az / pw = directions and direct powers
            % of this frame's (band, sector group) components (already
            % band-limited and diffuseness-gated). Levels are relative to a
            % running maximum (slow decay), novelty to the TauSlow
            % background.
            edges = -180:obj.BinDeg:180; nA = numel(edges) - 1;
            if isempty(obj.hF), obj.hF = zeros(nA, 1); obj.hS = []; end
            inst = zeros(nA, 1);
            if ~isempty(az)
                [~, ~, bi] = histcounts(az(:), edges);
                ok = bi > 0; pw = pw(:);
                inst = obj.cconv_same(accumarray(bi(ok), pw(ok), [nA 1]), obj.kernel());
            end
            a = exp(-dt/obj.TauHist);
            obj.hF = a*obj.hF + (1-a)*inst;
            if isempty(obj.hS), obj.hS = obj.hF; end
            as = exp(-dt/obj.TauSlow);
            obj.hS = as*obj.hS + (1-as)*inst;
            obj.gref = max(max(obj.hF), obj.gref*exp(-dt/obj.RefDecay));
            det = obj.peaks(obj.hF, obj.hS, obj.gref);
        end

        function [T, D] = track(obj, A)
            % Offline: detections from a soundfield_analyzer result, all
            % frames tracked; confirmed tracks with their logs
            [D, tf] = obj.detect(A);
            nF = numel(tf);
            dts = [median(diff(tf)), diff(tf)];
            kl = obj.KeepLog; obj.KeepLog = true;
            obj.reset();
            for f = 1:nF, obj.step(D{f}, tf(f), dts(f)); end
            done = [obj.done, obj.tr([obj.tr.conf])];
            obj.KeepLog = kl;
            T = struct('id', {}, 't', {}, 'az', {}, 'vel', {}, 'lev', {}, 'hit', {}, 'born', {}, 'novel', {});
            for k = 1:numel(done)
                L = done(k).log;
                % trim the coasting tail after the last hit
                last = find(L(:,5), 1, 'last'); L = L(1:last,:);
                T(end+1) = struct('id', done(k).id, 't', L(:,1)', 'az', L(:,2)', 'vel', L(:,3)', ...
                    'lev', L(:,4)', 'hit', logical(L(:,5)'), 'born', L(1,1), 'novel', done(k).novel); %#ok<AGROW>
            end
            [~, o] = sort([T.id]); T = T(o);
        end

        function [D, tf] = detect(obj, A)
            % Detections per analysis frame: D{f} = [az, level dB re max, novelty dB]
            fc = sqrt(A.bandEdges(1:end-1).*A.bandEdges(2:end));
            bands = fc >= obj.Band(1) & fc <= obj.Band(2);
            edges = -180:obj.BinDeg:180; azc = edges(1:end-1) + obj.BinDeg/2; nA = numel(azc);
            ker = obj.kernel();
            nF = numel(A.t); tf = A.t;
            dur = [median(diff(tf)), diff(tf)];
            H = zeros(nA, nF); HS = zeros(nA, nF); h = zeros(nA, 1); hs = [];
            for f = 1:nF
                ps = A.psi(bands,:,f);
                d = A.E(bands,:,f) .* (1 - ps); d(isnan(d) | ps > obj.PsiMax) = 0;
                az = A.az(bands,:,f);
                ok = d > 0 & ~isnan(az);
                [~, ~, bi] = histcounts(az(ok), edges);
                dd = d(ok); inst = accumarray(bi(bi > 0), dd(bi > 0), [nA 1]);
                inst = obj.cconv_same(inst, ker);
                a = exp(-dur(f)/obj.TauHist);
                h = a*h + (1-a)*inst;
                if isempty(hs), hs = h; end
                as = exp(-dur(f)/obj.TauSlow);
                hs = as*hs + (1-as)*inst;
                H(:,f) = h; HS(:,f) = hs;
            end
            gmax = max(H(:));
            D = cell(1, nF);
            for f = 1:nF, D{f} = obj.peaks(H(:,f), HS(:,f), gmax); end
        end
    end

    methods (Static, Access = private)
        function y = cconv_same(x, k)
            n = numel(x); m = (numel(k) - 1)/2;
            xp = [x(end-m+1:end); x; x(1:m)];
            y = conv(xp, k(:), 'valid');
            y = y(1:n);
        end
    end

    methods (Access = private)
        function det = peaks(obj, h, hs, gref)
            % Peaks of a direction histogram h -> [az, level dB re gref, novelty dB re hs]
            edges = -180:obj.BinDeg:180; azc = edges(1:end-1) + obj.BinDeg/2;
            det = zeros(0, 3);
            m = max(h);
            if m <= 0 || m <= gref*10^(obj.AbsFloorDb/10), return; end
            hp = [h(end); h; h(1)];
            pk = find(hp(2:end-1) > hp(1:end-2) & hp(2:end-1) >= hp(3:end));
            pk = pk(h(pk) > m*10^(obj.FloorDb/10));
            [~, o] = sort(h(pk), 'descend'); pk = pk(o);
            for i = pk(:)'
                y0 = hp(i); y1 = hp(i+1); y2 = hp(i+2);         % parabolic refinement
                dlt = 0.5*(y0 - y2)/(y0 - 2*y1 + y2); if ~isfinite(dlt) || abs(dlt) > 1, dlt = 0; end
                a0 = mod(azc(i) + dlt*obj.BinDeg + 180, 360) - 180;
                if isempty(det) || all(abs(mod(det(:,1) - a0 + 180, 360) - 180) >= obj.MinSepDeg)
                    det(end+1,:) = [a0, 10*log10(h(i)/gref), 10*log10(h(i)/max(hs(i), eps))]; %#ok<AGROW>
                end
                if size(det,1) == obj.MaxDet, break; end
            end
        end

        function ker = kernel(obj)
            g = -ceil(3*obj.SmoothDeg/obj.BinDeg):ceil(3*obj.SmoothDeg/obj.BinDeg);
            ker = exp(-(g*obj.BinDeg).^2/(2*obj.SmoothDeg^2)); ker = ker/sum(ker);
        end

        function r = meas_var(obj, levDb)
            % weaker detections are less certain (level re the excerpt max;
            % relative to the frame max would need the frame - use excerpt)
            s = obj.MeasDeg + obj.MeasDegPerDb*max(0, -levDb - 10);
            r = s^2;
        end
    end
end
