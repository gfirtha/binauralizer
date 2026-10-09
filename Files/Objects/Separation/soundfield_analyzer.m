classdef soundfield_analyzer < handle
    %SOUNDFIELD_ANALYZER  Direction / diffuseness of the sound field that a
    %channel mix would create at the listener (virtual DirAC analysis),
    %with transient-adaptive time resolution (block switching as in
    %MP3 / AAC encoders).
    %
    %   Every channel is taken as a plane wave from its loudspeaker
    %   direction u_i (unit vector, [right front up]). Per STFT bin:
    %     pressure   P = sum_i X_i,     velocity  V = sum_i X_i u_i   (rho c = 1)
    %     intensity  I = Re{P* V},      energy    E = (|P|^2 + |V|^2)/2
    %   summed into bands, then per band:
    %     direction  of I (azimuth 0 = front, + = left; elevation)
    %     diffuseness  psi = 1 - |<I>| / <E>   (<> = recursive average, Tau)
    %   A single loudspeaker / plane wave gives psi = 0 and its direction; a
    %   coherent phantom between two loudspeakers the Gerzon velocity vector
    %   direction; incoherent channels the energy vector (sum |X_i|^2 u_i).
    %
    %   Block switching: the signal runs in granules of Hop samples. A
    %   granule is analysed with one long window (NLong, centred on the
    %   granule) unless an attack is detected in it (high-passed energy of
    %   a sub-block > AttackRatio x the mean of the previous NSub
    %   sub-blocks) or a spatial transient (the high-passed energy vector
    %   sum_i e_i u_i / sum_i e_i turns by more than SpatialDeg between the
    %   two halves of the granule - a fly-by is a fast crescendo with a
    %   fast direction change, not an onset an MP3-style detector sees);
    %   then with NSub short windows (2 x Hop/NSub long, centred on the
    %   sub-blocks), but only for the bands above XoverHz (the low bands
    %   keep the granule's long window: a 256-sample window cannot resolve
    %   them anyway). The short mode is held for HoldGran further granules.
    %
    %   Sectors (Mode 'sectors', after Politis & Pulkki's sector-based
    %   parametric analysis): one intensity vector cannot show two
    %   simultaneous directions (front shot + rear fly-by average into a
    %   front-ish vector of reduced length). K sector beams of order N,
    %   w_s(theta) = (1 + 2 sum_m a_m cos(m (theta - theta_s)))/K with 2D
    %   max-rE weights a_m = cos(m pi/(2N+2)), K = 2N+2 regular directions
    %   (a circular design: sum_s w_s = 1, so the sectors partition P and V):
    %     P_s = sum_i w_s(theta_i) X_i,   V_s = sum_i w_s(theta_i) X_i u_i
    %   and the analysis above per sector. A single plane wave gives its
    %   true direction in every sector (weight w_s^2); a phantom keeps the
    %   intensity (Gerzon) direction instead of splitting into loudspeaker
    %   peaks, as a sharp high-order beamformer would - but each sector
    %   weights the two loudspeakers differently, so a phantom is smeared
    %   over its pair. Hence neighbouring sectors whose (smoothed, band)
    %   pressures are coherent (> MergeCoh: they see the same source) are
    %   merged and analysed together (sum of their P and V); uncorrelated
    %   sources stay in separate sectors. Horizontal only;
    %   3D layouts would use a spherical t-design.
    %
    %   A = analyze(x, fs) -> struct: per analysis frame t (centre [s]),
    %   dur, isShort; per band, sector and frame az, el, psi, E; column grid
    %   tc / Hist (direct power E (1 - psi) per azimuth bin AzEdges, all
    %   sectors, on a uniform Hop/NSub grid), SecE / SecAz (direct power and
    %   direction per sector), psiCol, attack flags per granule.

    properties
        Layout = []         % [nCh x 3] unit vectors ([] = from channel_layout)
        NLong = 2048        % long window
        Hop = 1024          % granule
        NSub = 8            % short windows per granule (short window 2 Hop/NSub)
        Window = 'hann'
        AttackRatio = 8     % sub-block energy jump that triggers short blocks (~9 dB)
        AttackHP = 2000     % [Hz] high-pass of the attack detector
        SpatialDeg = 15     % energy-vector turn within a granule that triggers short blocks
        SpatialMinR = 0.2   % ... only if the energy vector is this directional
        SpatialFloorDb = -40 % ... and the granule within this of the loudest one
        HoldGran = 2        % granules kept short after an attack
        XoverHz = 1500      % in short granules, bands below stay at long resolution
        Tau = 0.02          % [s] averaging of the diffuseness statistics
        BandEdges = []      % [Hz] ([] = 20 log-spaced bands, 180 Hz .. 16 kHz)
        AzEdges = -180:5:180
        Adaptive = true     % false: long windows only (for comparison)
        Mode = 'single'     % 'single' (one intensity vector) | 'sectors'
        SectorOrder = 2     % beam order of the sector patterns
        NSectors = []       % [] = 2 SectorOrder + 2 (regular, horizontal)
        MergeCoh = 0.8      % neighbouring sectors with this pressure coherence are merged
        HistBand = [0 Inf]  % [Hz] bands that enter Hist / SecE / SecAz
    end

    methods
        function obj = soundfield_analyzer(varargin)
            for k = 1:2:numel(varargin), obj.(varargin{k}) = varargin{k+1}; end
        end

        function A = analyze(obj, x, fs)
            nCh = size(x, 2);
            U = obj.Layout;
            if isempty(U)
                az = separator_base.channel_layout(nCh);
                keep = ~isnan(az);                          % LFE out
                x = x(:, keep); az = az(keep);
                U = separator_base.dir2pos(az);
            end
            be = obj.BandEdges;
            if isempty(be), be = logspace(log10(180), log10(16000), 21); end
            nB = numel(be) - 1;
            H = obj.Hop; ns = H/obj.NSub; NS = 2*ns; NL = obj.NLong;
            n = size(x, 1); nG = floor(n/H);

            [fr, gAtt, gSpat, gShort] = obj.block_plan(x, U, fs);
            nF = size(fr, 1);
            xp = [zeros(NL, size(x,2)); x; zeros(NL, size(x,2))];
            win = @(N) 0.5*(1 - cos(2*pi*((0:N-1)' + 0.5)/N));
            wL = win(NL); wS = win(NS);
            wL = wL/sqrt(sum(wL.^2)); wS = wS/sqrt(sum(wS.^2));   % same energy scale
            mapL = band_map(NL); mapS = band_map(NS);

            % ---- sectors: beams w_s(theta) of order N (max-rE), sum_s w_s = 1 ----
            if strcmpi(obj.Mode, 'sectors')
                [W, WU, secAz] = soundfield_analyzer.sector_weights(U, obj.SectorOrder, obj.NSectors);
                K = size(W, 2);
            else
                K = 1; secAz = NaN; W = ones(size(U,1), 1);
                WU = cat(3, W.*U(:,1), W.*U(:,2), W.*U(:,3));
            end

            A.az = nan(nB, K, nF); A.el = nan(nB, K, nF); A.psi = nan(nB, K, nF); A.E = zeros(nB, K, nF);
            A.nGroups = zeros(nB, nF);
            Sa = [];                                        % smoothed band cross statistics
            lf = sqrt(be(1:end-1).*be(2:end))' < obj.XoverHz;   % bands kept at long resolution
            aL = exp(-(NL/2)/(fs*obj.Tau)); aS = exp(-(NS/2)/(fs*obj.Tau));
            gL = 0;                                         % granule of the cached long frame
            for f = 1:nF
                N = fr(f,2); c = round(fr(f,1));
                if N == NL
                    Sb = band_stats(c, NL, wL, mapL);
                    if isempty(Sa), Sa = Sb; else, Sa = aL*Sa + (1-aL)*Sb; end
                else
                    % short frame: high bands from it, low bands from the
                    % granule's long frame (multi-resolution, no LF speckle)
                    g = floor(c/H) + 1;
                    if g ~= gL
                        SbL = band_stats((g-1)*H + H/2, NL, wL, mapL);
                        if isempty(Sa), Sa = SbL; end
                        Sa(lf,:) = aL*Sa(lf,:) + (1-aL)*SbL(lf,:);
                        gL = g;
                    end
                    Sb = band_stats(c, NS, wS, mapS);
                    Sb(lf,:) = SbL(lf,:);
                    Sa(~lf,:) = aS*Sa(~lf,:) + (1-aS)*Sb(~lf,:);
                end
                for b = 1:nB
                    [I, E, Ism, Esm, slot] = group_stats(Sb(b,:), Sa(b,:));
                    A.nGroups(b,f) = numel(slot);
                    A.psi(b,slot,f) = 1 - sqrt(sum(Ism.^2, 2))./max(Esm, eps);
                    A.az(b,slot,f) = atan2d(-I(:,1), I(:,2));
                    A.el(b,slot,f) = atan2d(I(:,3), sqrt(I(:,1).^2 + I(:,2).^2));
                    A.E(b,slot,f) = E;
                end
            end
            A.psi = min(max(A.psi, 0), 1);
            A.t = fr(:,1)'/fs; A.dur = fr(:,2)'/fs; A.isShort = (fr(:,2) == NS)';
            A.bandEdges = be; A.gAttack = gAtt; A.gSpatial = gSpat; A.gShort = gShort; A.fs = fs;
            A.AzEdges = obj.AzEdges; A.sectorAz = secAz; A.W = W;

            % ---- uniform column grid (sub-block resolution) for maps ----
            nC = nG*obj.NSub;
            A.tc = ((1:nC) - 0.5)*ns/fs;
            nA = numel(obj.AzEdges) - 1;
            A.Hist = zeros(nA, nC);                         % direct power per azimuth (all sectors)
            A.SecE = zeros(K, nC);                          % direct power per sector
            A.SecAz = nan(K, nC);                           % its power-weighted direction
            A.psiCol = zeros(1, nC);
            [~, ~, ab] = histcounts(A.az, obj.AzEdges);
            wd = A.E .* (1 - A.psi);                        % direct power
            fc = sqrt(be(1:end-1).*be(2:end));
            wd(fc < obj.HistBand(1) | fc > obj.HistBand(2), :, :) = 0;
            col = 0;
            for f = 1:nF
                if A.isShort(f), cols = col + 1; else, cols = col + (1:obj.NSub); end
                col = cols(end);
                a = ab(:,:,f); d = wd(:,:,f); d(isnan(d)) = 0; ok = a > 0;
                h = accumarray(a(ok), d(ok), [nA 1]);
                A.Hist(:, cols) = repmat(h, 1, numel(cols));
                se = sum(d, 1)';
                z = sum(d .* exp(1i*A.az(:,:,f)*pi/180), 1, 'omitnan').';
                A.SecE(:, cols) = repmat(se, 1, numel(cols));
                A.SecAz(:, cols) = repmat(angle(z)*180/pi, 1, numel(cols));
                Ef = A.E(:,:,f);
                pf = A.psi(:,:,f); pf(isnan(pf)) = 0;
                A.psiCol(cols) = sum(Ef(:).*pf(:))/max(sum(Ef(:)), eps);
            end

            function Sb = band_stats(c, N, w, M)
                % band cross statistics of the sector signals of the frame
                % of length N centred at c
                X = fft(xp(NL + c - N/2 + (1:N), :) .* w, [], 1);
                Sb = soundfield_analyzer.sector_cross(X(1:N/2+1, :), W, WU, M)/N;
            end

            function [I, E, Ism, Esm, slot] = group_stats(sb, sa)
                [~, I, E, Ism, Esm, slot] = soundfield_analyzer.group_sectors(sb, sa, K, obj.MergeCoh);
            end

            function M = band_map(N)
                fb = (0:N/2)*fs/N;
                M = zeros(nB, N/2+1);
                for b = 1:nB
                    M(b, fb >= be(b) & fb < be(b+1)) = 1;
                end
                % bands narrower than a short-window bin: nearest bin, shared
                for b = find(~any(M, 2))'
                    [~, k] = min(abs(fb - sqrt(be(b)*be(b+1)))); M(b,k) = 1;
                end
                M = M ./ max(sum(M, 1), 1);
            end
        end

        function [fr, gAtt, gSpat, gShort] = block_plan(obj, x, U, fs)
            % Block switching plan: frames [centre sample, window length]
            % (long NLong per granule, or NSub short windows where an energy
            % attack or a spatial transient is detected, held HoldGran
            % granules); flags per granule. U = loudspeaker unit vectors.
            H = obj.Hop; ns = H/obj.NSub; NS = 2*ns; NL = obj.NLong;
            n = size(x, 1); nG = floor(n/H);
            % ---- attack detection (high-passed energy in sub-blocks) ----
            [bh, ah] = butter(2, obj.AttackHP/(fs/2), 'high');
            xh = filter(bh, ah, x).^2;
            eh = sum(xh, 2);
            es = sum(reshape(eh(1:nG*H), ns, []), 1);       % sub-block energies
            fl = 1e-8*max(es);
            att = false(1, numel(es));
            for j = obj.NSub+1:numel(es)
                att(j) = es(j) > obj.AttackRatio*mean(es(j-obj.NSub:j-1)) && es(j) > fl;
            end
            gAtt = any(reshape(att, obj.NSub, []), 1);
            % spatial transients: energy vector of the granule halves
            eHalf = squeeze(sum(reshape(xh(1:nG*H,:), H/2, 2*nG, []), 1));   % [2 nG x nCh]
            rv = eHalf*U ./ max(sum(eHalf, 2), eps);
            r1 = rv(1:2:end,:); r2 = rv(2:2:end,:);
            nr1 = sqrt(sum(r1.^2, 2)); nr2 = sqrt(sum(r2.^2, 2));
            turn = acosd(min(1, sum(r1.*r2, 2)./max(nr1.*nr2, eps)))';
            eg = sum(reshape(eh(1:nG*H), H, []), 1);
            gSpat = turn > obj.SpatialDeg & min(nr1, nr2)' > obj.SpatialMinR & ...
                eg > max(eg)*10^(obj.SpatialFloorDb/10);
            gShort = gAtt | gSpat;
            for g = find(gShort), gShort(g:min(nG, g+obj.HoldGran)) = true; end
            if ~obj.Adaptive, gShort(:) = false; end

            % ---- frames: [centre sample, window length] ----
            fr = zeros(0, 2);
            for g = 1:nG
                if gShort(g)
                    c = (g-1)*H + ns*((1:obj.NSub) - 0.5);
                    fr = [fr; c(:), NS*ones(obj.NSub,1)]; %#ok<AGROW>
                else
                    fr(end+1,:) = [(g-1)*H + H/2, NL]; %#ok<AGROW>
                end
            end
        end
    end

    methods (Static)
        function [W, WU, secAz] = sector_weights(U, order, K)
            % Sector beams w_s(theta_i) of order 'order' (2D max-rE) for
            % loudspeaker unit vectors U [nCh x 3]; K regular sectors
            % (default 2 order + 2); sum_s w_s = 1. WU = W .* U (velocity).
            if nargin < 3 || isempty(K), K = 2*order + 2; end
            secAz = (0:K-1)*360/K;
            chAz = atan2d(-U(:,1), U(:,2));
            am = cos((1:order)*pi/(2*order + 2));
            W = ones(numel(chAz), K)/K;
            for m = 1:order, W = W + 2/K*am(m)*cosd(m*(chAz - secAz)); end
            WU = cat(3, W.*U(:,1), W.*U(:,2), W.*U(:,3));
        end

        function Sb = sector_cross(X, W, WU, M)
            % Band sums (rows of M) of the sector cross statistics of the
            % spectrum X [bins x nCh]: [nB x 6 K^2] =
            % [Re{P_s* V_t} (x, y, z), Re{P_s* P_t}, Im{P_s* P_t}, Re{V_s . V_t}]
            nb = size(X, 1); K = size(W, 2);
            P = X*W;                                        % sector pressures [nb x K]
            V = reshape(X*reshape(WU, size(WU,1), []), nb, K, 3);   % sector velocities
            cP = reshape(conj(P), nb, K, 1);
            PV = real(cP .* reshape(V, nb, 1, K, 3));       % [nb x K x K x 3]
            PP = cP .* reshape(P, nb, 1, K);
            VV = real(sum(conj(reshape(V, nb, K, 1, 3)) .* reshape(V, nb, 1, K, 3), 4));
            Sb = M*[reshape(PV, nb, []), reshape(real(PP), nb, []), ...
                    reshape(imag(PP), nb, []), reshape(VV, nb, [])];
        end

        function [grp, I, E, Ism, Esm, slot] = group_sectors(sb, sa, K, mergeCoh)
            % Merge ring-neighbour sectors whose (smoothed) pressures are
            % coherent (> mergeCoh: they see the same source), then
            % intensity / energy per group from the instantaneous sb and
            % smoothed sa statistics of one band (rows of sector_cross).
            % slot = strongest sector of each group.
            K2 = K*K;
            pp = reshape(sa(3*K2 + (1:K2)) + 1i*sa(4*K2 + (1:K2)), K, K);
            grp = num2cell(1:K);
            while numel(grp) > 1
                ng = numel(grp); best = 0; bi = 0;
                for gi = 1:ng
                    gj = mod(gi, ng) + 1;                   % ring neighbour
                    if ng == 2 && gi == 2, break; end
                    a = grp{gi}; b2 = grp{gj};
                    coh = abs(sum(pp(a,b2), 'all'))/sqrt(max(real(sum(pp(a,a), 'all'))*real(sum(pp(b2,b2), 'all')), eps));
                    if coh > best, best = coh; bi = gi; end
                end
                if best < mergeCoh, break; end
                gj = mod(bi, ng) + 1;
                grp{bi} = [grp{bi}, grp{gj}]; grp(gj) = [];
            end
            ng = numel(grp);
            I = zeros(ng, 3); Ism = zeros(ng, 3); E = zeros(ng, 1); Esm = zeros(ng, 1); slot = zeros(1, ng);
            for gi = 1:ng
                a = grp{gi};
                lin = reshape(a(:) + (a(:)' - 1)*K, [], 1);
                for d = 1:3
                    I(gi,d) = sum(sb((d-1)*K2 + lin)); Ism(gi,d) = sum(sa((d-1)*K2 + lin));
                end
                E(gi) = 0.5*(sum(sb(3*K2 + lin)) + sum(sb(5*K2 + lin)));
                Esm(gi) = 0.5*(sum(sa(3*K2 + lin)) + sum(sa(5*K2 + lin)));
                [~, k] = max(real(diag(pp(a,a)))); slot(gi) = a(k);   % strongest sector holds the group
            end
        end
    end
end
