classdef wavenumber_separator < separator_base
    %WAVENUMBER_SEPARATOR  Real-time separation into moving objects from
    %the local wavefront at the listening position: local wavenumber
    %vector and phase Hessian of the field the mix creates there (ray
    %theory, see wavenumber_analyzer), per critical band, tracked by a
    %Kalman filter that predicts the change of the local wavenumber
    %direction (soundfield_tracker), objects in object_slots; the
    %wavefront curvature gives each object a distance.
    %
    %   Input: stereo is first decoded Lt/Rt -> L C R Ls Rs (ltrt_decoder,
    %   in the same STFT); 5.0/5.1/7.1 channels are used as they are (LFE
    %   dropped). Channels = point sources at R m (closed-form P, grad P,
    %   Hess P at the origin).
    %   Per hop (N/2) and band (Zwicker critical bands):
    %   - direction: power-weighted mean k^ (instantaneous, and smoothed
    %     over TauEst for standing objects); direct share fd = directional
    %     consistency^2 = |power-weighted mean of the unit k^ of the bins|^2
    %     (1 when all bins of the band point the same way - one source or
    %     a coherent phantom, whose |k^| < 1 is geometry, not diffuseness;
    %     < 1 for a spread / incoherent mixture);
    %   - distance: phase-Hessian curvature across the ray, weighted by a
    %     simple-field quality per bin, q = exp(-((|k^|-1)/KmagSigma)^2) x
    %     exp(-(R ((kappaH+kappaV)/2 - kappaTr)/TrErrSigma)^2) (the
    %     transport equation holds for one wavefront only); rho = 1/kappaH
    %     where the band's mean q > MinQ and rho is in RhoRange;
    %   - block switching of the analysis: an energy attack (high-passed,
    %     sub-blocks of N/8) or a turn of the channel energy vector within
    %     the hop -> four short (ShortN) analysis steps for the tracker
    %     (low bands from the long frame), held HoldHops hops;
    %   - detection (bands above DetBand(1), power x fd) -> tracker
    %     (Kalman, [az; angular velocity]) -> object_slots;
    %   - synthesis: band content (the omni pressure at the origin, P0 =
    %     sum_i X_i, scaled to the band's channel power x fd x weight) to
    %     the objects by object_slots.weights; the rest stays in the bed
    %     (the analysed channels, at their loudspeaker directions);
    %   - object distance: log-domain smoothing of the assigned bands'
    %     valid rho (DistTau moving, DistTauStill standing; R when none);
    %     Pos = direction x rho/RefDist (the scene maps |Pos| = 1 to its
    %     object room radius).

    properties
        K = 6
        N = 1024
        R = 2               % loudspeaker distance of the model [m]
        RefDist = 3         % distance mapped to |Pos| = 1 [m] (scene object radius)
        c = 343
        BandEdges = [0 100 200 300 400 510 630 770 920 1080 1270 1480 1720 2000 2320 ...
                     2700 3150 3700 4400 5300 6400 7700 9500 12000 15500 inf]
        MinHz = 20
        DetBand = [1500 Inf]
        XoverHz = 1500      % short analysis frames: bands below from the long frame
        TauEst = 0.02       % [s] smoothing of band statistics (standing objects, distance)
        KmagSigma = 0.15
        TrErrSigma = 0.3
        MinQ = 0.5
        RhoRange = [0.2 20] % [m]
        UseDistance = false % objects get the curvature distance (else all at R); off by
                            % default: on mixes it is ~R anyway and only adds radial jitter
        DecTau = 0.2        % [s] Lt/Rt decoder steering (a fast steering moves the decoded
                            % field and the objects with it)
        DistTau = 0.1       % [s] distance smoothing, moving object
        DistTauStill = 0.5  % [s] ... standing object
        ShortN = 256
        AttackRatio = 8
        AttackHP = 2000
        SpatialDeg = 15
        HoldHops = 2
        AmbGain = 1
        AmbAz = 110         % (common separator settings)
        Width = 180
        Base = 30
        MinAge = 0.12; FadeIn = 0.04; FadeOut = 0.15
        PosTau = 0.06; PosTauStill = 0.7
        AssignSigma = 12; AssignDeg = 32
        DetMinConsistency = 0.35   % minimum directional consistency^2 for tracking
        DetFrameFloorDb = -24      % ignore tracking candidates this far below frame max
        Tracker
        Diagnostics = false
    end

    properties (SetAccess = private)
        DiagAz = []; DiagPA = []
        AzEdges = -180:4:180
        TrackLog = []       % [K x frames] object azimuths (NaN = inactive)
        RhoLog = []         % [K x frames] object distances [m]
        ShortLog = []       % [1 x frames] short analysis used
        BandDiag = []       % last frame per band: [srcAz, fd, sum of object weights, |P0|^2/Pc]
        SlotDiag = []       % [3 x K x frames]: track id, fade gain, moving
    end

    properties (Access = private)
        hop; win; winS; frame; acc; pending; outq
        stereo; Dec; DecS; chIdx; Uc; Rv; nBed
        SL; SS; nB; bandOfFull; bandLo; lfS
        Sa                  % smoothed band statistics
        Slots; slotRho; t
        hpB; hpA; hpZ; eHist; hold; Uin
    end

    methods
        function obj = wavenumber_separator(varargin)
            obj.Tracker = soundfield_tracker('DupDeg', 20, 'AccNoiseStable', 100, ...
                'StableDeg', 20, 'StableTime', 0.12, ...
                'TauHist', 0.03, 'SmoothDeg', 10, 'FloorDb', -18, ...
                'AbsFloorDb', -45, 'MinSepDeg', 35, 'MaxDet', 4, ...
                'ConfirmTime', 0.04, 'MissTime', 0.25, 'GateDeg', 12);
            for k = 1:2:numel(varargin), obj.(varargin{k}) = varargin{k+1}; end
        end

        function init(obj, fs, nIn)
            obj.fs = fs; obj.nIn = nIn;
            obj.stereo = nIn == 2;
            if obj.stereo
                obj.Dec = ltrt_decoder('TauSteer', obj.DecTau); obj.Dec.setup_bands(obj.N, fs);
                obj.DecS = ltrt_decoder('TauSteer', obj.DecTau); obj.DecS.setup_bands(obj.ShortN, fs);
                az = [30 0 -30 110 -110]; names = {'L','C','R','Ls','Rs'};
                obj.chIdx = 1:2;
                obj.Uin = separator_base.dir2pos([30 -30]);
            else
                [az, names] = separator_base.channel_layout(nIn);
                obj.chIdx = find(~isnan(az));
                az = az(obj.chIdx); names = names(obj.chIdx);
                if nIn == 1, az = 0; end
                obj.Uin = separator_base.dir2pos(az);
            end
            obj.Uc = separator_base.dir2pos(az);
            obj.Rv = obj.R*ones(1, numel(az));
            obj.nBed = numel(az);
            obj.hop = obj.N/2;
            obj.win = sqrt(0.5*(1 - cos(2*pi*(0:obj.N-1)'/obj.N)));
            obj.winS = 0.5*(1 - cos(2*pi*((0:obj.ShortN-1)' + 0.5)/obj.ShortN));
            be = obj.BandEdges;
            obj.SL = wavenumber_analyzer.precompute(obj.N, fs, obj.Uc, obj.Rv, obj.win, be, obj.MinHz, obj.c, 'bands');
            obj.SS = wavenumber_analyzer.precompute(obj.ShortN, fs, obj.Uc, obj.Rv, obj.winS, be, obj.MinHz, obj.c, 'bands');
            obj.nB = size(obj.SL.Bin, 1);
            nb = obj.N/2 + 1;
            [bi, col] = find(obj.SL.Bin);
            idxUse = find(obj.SL.use);
            obj.bandOfFull = zeros(nb, 1); obj.bandOfFull(idxUse(col)) = bi;
            obj.bandLo = accumarray(bi, obj.SL.f(col), [obj.nB 1], @min, inf);
            fc = sqrt(max(be(1:end-1), 1).*min(be(2:end), obj.fs/2))';
            obj.lfS = fc(1:obj.nB) < obj.XoverHz;
            obj.Names = [arrayfun(@(k) sprintf('Obj %d', k), 1:obj.K, 'UniformOutput', false), ...
                strcat('Amb', {' '}, names(:)')];
            obj.Pos = [separator_base.dir2pos(linspace(0, 300, obj.K)'); separator_base.dir2pos(az(:))];
            obj.Latency = obj.N - obj.hop;
            [obj.hpB, obj.hpA] = butter(2, obj.AttackHP/(fs/2), 'high');
            obj.Tracker.MaxDet = min(obj.Tracker.MaxDet, obj.K);
            obj.Slots = object_slots('K', obj.K, 'MinAge', obj.MinAge, 'FadeIn', obj.FadeIn, ...
                'FadeOut', obj.FadeOut, 'PosTau', obj.PosTau, 'PosTauStill', obj.PosTauStill, ...
                'AssignSigma', obj.AssignSigma, 'AssignDeg', obj.AssignDeg, ...
                'DupDeg', obj.Tracker.DupDeg, 'StableTime', obj.Tracker.StableTime);
            obj.reset();
        end

        function setDecoderTau(obj, tau)
            % steering time constant of the Lt/Rt decoders (stereo input)
            if obj.stereo, obj.Dec.TauSteer = tau; obj.DecS.TauSteer = tau; end
        end

        function reset(obj)
            obj.frame = zeros(obj.N, obj.nIn);
            obj.acc = zeros(obj.N, obj.nOut());
            obj.pending = zeros(0, obj.nIn);
            obj.outq = zeros(0, obj.nOut());
            obj.Sa = [];
            obj.Slots.reset(); obj.Tracker.reset();
            if obj.stereo, obj.Dec.reset_steering(); obj.DecS.reset_steering(); end
            obj.slotRho = obj.R*ones(1, obj.K);
            obj.t = 0;
            obj.hpZ = zeros(2, numel(obj.chIdx));
            obj.eHist = zeros(1, 8); obj.hold = 0;
            obj.DiagAz = []; obj.DiagPA = []; obj.TrackLog = []; obj.RhoLog = []; obj.ShortLog = [];
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
            % whole signal, aligned (no latency); Xs{1} only
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
            nb = obj.N/2 + 1; K = obj.K; dt = obj.hop/obj.fs; hop = obj.hop;
            obj.t = obj.t + dt;
            F = fft(obj.frame .* obj.win, [], 1); F = F(1:nb, :);
            if obj.stereo
                Xall = obj.Dec.decode_spectrum(F(:,1), F(:,2), dt);
            else
                Xall = F(:, obj.chIdx);
            end
            use = obj.SL.use;
            q = wavenumber_analyzer.bin_stats(Xall(use,:), obj.SL);
            Bl = obj.band_agg(q, obj.SL);
            a = exp(-dt/obj.TauEst);
            if isempty(obj.Sa), obj.Sa = Bl; else
                fn = fieldnames(Bl);
                for i = 1:numel(fn), obj.Sa.(fn{i}) = a*obj.Sa.(fn{i}) + (1-a)*Bl.(fn{i}); end
            end
            azI = atan2d(Bl.sx, -Bl.sy);                    % source directions (+ = left)
            azS = atan2d(obj.Sa.sx, -obj.Sa.sy);
            fd = min(hypot(Bl.ux, Bl.uy)./max(Bl.E, eps), 1).^2;   % directional consistency^2
            qual = obj.Sa.sq./max(obj.Sa.E, eps);
            rho = obj.Sa.sq./obj.Sa.sH;                      % 1 / quality-weighted kappaH
            rhoOk = obj.UseDistance & qual > obj.MinQ & rho >= obj.RhoRange(1) & rho <= obj.RhoRange(2);
            det = obj.bandLo >= obj.DetBand(1) & obj.bandLo <= obj.DetBand(2);

            % ---- transient detection on the newest hop -> analysis steps ----
            [xh, obj.hpZ] = filter(obj.hpB, obj.hpA, obj.frame(end-hop+1:end, obj.chIdx), obj.hpZ);
            ns = hop/4;
            es = sum(reshape(sum(xh.^2, 2), ns, 4), 1);
            att = false;
            for j = 1:4
                if es(j) > obj.AttackRatio*mean(obj.eHist) && es(j) > 1e-10, att = true; end
                obj.eHist = [obj.eHist(2:end), es(j)];
            end
            if size(obj.Uin, 1) > 2                          % spatial transient (multichannel)
                e1 = sum(xh(1:hop/2,:).^2, 1); e2 = sum(xh(hop/2+1:end,:).^2, 1);
                r1 = e1*obj.Uin/max(sum(e1), eps); r2 = e2*obj.Uin/max(sum(e2), eps);
                if min(norm(r1), norm(r2)) > 0.2 && ...
                        acosd(min(1, dot(r1, r2)/max(norm(r1)*norm(r2), eps))) > obj.SpatialDeg
                    att = true;
                end
            end
            if att, obj.hold = obj.HoldHops + 1; end
            short = obj.hold > 0; obj.hold = max(0, obj.hold - 1);
            % long-frame band powers on the short frames' scale (window energy, bins per band)
            fL = sum(obj.win.^2)*obj.N/obj.ShortN;
            if short
                % four short windows ending every hop/4 within the newest hop
                L = obj.ShortN; nbS = L/2 + 1;
                for j = 1:4
                    e = obj.N - (4-j)*ns;
                    Fs = fft(obj.frame(e-L+1:e, :) .* obj.SS.w, [], 1); Fs = Fs(1:nbS, :);
                    if obj.stereo
                        Xs = obj.DecS.decode_spectrum(Fs(:,1), Fs(:,2), ns/obj.fs);
                    else
                        Xs = Fs(:, obj.chIdx);
                    end
                    qs = wavenumber_analyzer.bin_stats(Xs(obj.SS.use,:), obj.SS);
                    Bs = obj.band_agg(qs, obj.SS);
                    lf = obj.lfS | Bs.E == 0;
                    fn = {'E','sx','sy','ux','uy'};
                    for i = 1:numel(fn), Bs.(fn{i})(lf) = Bl.(fn{i})(lf)/fL; end
                    aI = atan2d(Bs.sx, -Bs.sy);
                    pw = Bs.E.*min(hypot(Bs.ux, Bs.uy)./max(Bs.E, eps), 1).^2;
                    obj.track_step(aI(det), pw(det), obj.t - dt + j*ns/obj.fs, ns/obj.fs);
                end
            else
                obj.track_step(azI(det), Bl.E(det).*fd(det)/fL, obj.t, dt);
            end
            sl = obj.Slots;

            % ---- object distances ----
            for k = find(sl.active)
                wk = zeros(obj.nB, 1);
                for b = find(rhoOk)'
                    w = sl.weights(azS(b), azI(b)); wk(b) = w(k)*Bl.E(b);
                end
                if sl.renewed(k), obj.slotRho(k) = obj.R; end
                if sum(wk) > 0
                    m = exp(sum(wk.*log(rho))/sum(wk));
                    if sl.moving(k), tau = obj.DistTau; else, tau = obj.DistTauStill; end
                    obj.slotRho(k) = exp(log(obj.slotRho(k)) + (log(m) - log(obj.slotRho(k)))*(1 - exp(-dt/tau)));
                end
            end

            % ---- synthesis ----
            P0 = sum(Xall, 2);                               % omni pressure (equal R: common G)
            Pc = obj.SL.Bin*sum(abs(Xall(use,:)).^2, 2);
            P0b = obj.SL.Bin*abs(P0(use)).^2;
            S = zeros(nb, K + obj.nBed);
            bedFrac = ones(obj.nB, 1);
            for b = 1:obj.nB
                if Pc(b) <= 0 || P0b(b) < 0.1*Pc(b), continue; end   % cancelling omni: all to the bed
                [w, sw] = sl.weights(azS(b), azI(b));
                if sw <= 0, continue; end
                bins = obj.bandOfFull == b;
                for k = find(w > 0)
                    S(bins, k) = S(bins, k) + sqrt(Pc(b)*fd(b)*w(k)/P0b(b))*P0(bins);
                end
                bedFrac(b) = 1 - fd(b)*sw;
            end
            if obj.Diagnostics
                swAll = zeros(obj.nB, 1);
                for b = 1:obj.nB, [~, swAll(b)] = sl.weights(azS(b), azI(b)); end
                obj.BandDiag = [azI, fd, swAll, P0b./max(Pc, eps)];
                obj.SlotDiag(:, :, end+1) = [sl.id; sl.g; sl.moving];
            end
            bf = ones(nb, 1); okb = obj.bandOfFull > 0; bf(okb) = bedFrac(obj.bandOfFull(okb));
            S(:, K+1:end) = obj.AmbGain*sqrt(max(bf, 0)).*Xall;

            % ---- positions (direction x distance), diagnostics ----
            sc = min(max(obj.slotRho(:)/obj.RefDist, 0.1), 2);
            obj.Pos(1:K,:) = separator_base.dir2pos(sl.az(:)).*sc;
            if obj.Diagnostics
                [~, ~, bi] = histcounts(azI, obj.AzEdges); ok = bi > 0;
                obj.DiagAz(:, end+1) = accumarray(bi(ok), Bl.E(ok).*fd(ok), [numel(obj.AzEdges)-1, 1]);
                obj.DiagPA(:, end+1) = [sum(sum(abs(S(:,1:K)).^2)); sum(sum(abs(S(:,K+1:end)).^2))];
                lg = sl.az(:); lg(~sl.active) = NaN; obj.TrackLog(:, end+1) = lg;
                rl = obj.slotRho(:); rl(~sl.active) = NaN; obj.RhoLog(:, end+1) = rl;
                obj.ShortLog(end+1) = short;
            end
            Yf = real(ifft([S; conj(S(end-1:-1:2,:))], [], 1)) .* obj.win;
        end

        function track_step(obj, az, pw, t, dt)
            det = obj.Tracker.detect_frame(az, pw, dt);
            tr = obj.Tracker.step(det, t, dt);
            obj.Slots.update(tr, t, dt);
        end

        function Bq = band_agg(obj, q, S)
            % band sums: power, power-weighted k^, simple-field quality and curvature
            B = S.Bin;
            qb = exp(-((q.km - 1)/obj.KmagSigma).^2) .* ...
                 exp(-(obj.R*((q.kapH + q.kapV)/2 - q.kapTr)/obj.TrErrSigma).^2);
            Bq.E = B*q.E;
            Bq.sx = B*(q.E.*q.kx); Bq.sy = B*(q.E.*q.ky);
            Bq.ux = B*(q.E.*q.kx./max(q.km, eps)); Bq.uy = B*(q.E.*q.ky./max(q.km, eps));   % unit directions
            Bq.sq = B*(q.E.*qb); Bq.sH = B*(q.E.*qb.*q.kapH);
        end
    end
end
