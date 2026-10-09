classdef avendano_separator < separator_base
    %AVENDANO_SEPARATOR  Real-time stereo upmix after Avendano & Jot
    %(coherence-based ambience extraction + similarity panning index),
    %output as fixed-direction objects.
    %
    %   Per STFT bin, from the recursively smoothed auto / cross spectra
    %   (time constant Tau)  Pll = <|L|^2>, Prr = <|R|^2>, Plr = <L R*>:
    %     coherence   phi = |Plr| / sqrt(Pll Prr)
    %     ambience    Gam = (Mu1+Mu0)/2 + (Mu1-Mu0)/2 tanh(Sigma pi (Phi0 - phi))
    %                 (nonlinear mask: low coherence -> ambience)
    %     similarity  psi = 2|L R*| / (|L|^2 + |R|^2)   (instantaneous, as in
    %                 the patents: a pure level cue, unbiased by ambience),
    %                 partial psiL = |L R*|/|L|^2, psiR = |L R*|/|R|^2
    %     panning     Psi = (1 - psi) sign(psiL - psiR)   (-1 left .. +1 right;
    %                 the papers' text says "positive = left", their
    %                 formulas give positive = right - the formulas are used)
    %   Ambience  A = Gam X (both channels), direct D = sqrt(1-Gam^2) X
    %   (power complementary). The panning index is inverted to the
    %   unit-power amplitude-panning vector a = [cos t, sin t] with
    %   sin 2t = 1 - |Psi| (exact for a single amplitude-panned source),
    %   the direct part is downmixed onto it (m = a' D) and panned
    %   (constant power) between the two neighbouring of K fixed front
    %   slots over +-Width (direction: stereo tangent law of a on the
    %   source base +-Base, stretched to Width). What the downmix misses
    %   (D - a m: inter-channel phase, model mismatch) joins the ambience
    %   pair, so the objects keep the input power.
    %
    %   Deviation from the papers (HardPanGuard): a hard-panned source has
    %   phi ~ 0 (the other channel is empty) and would be taken as
    %   ambience; the mask is therefore scaled by 1 - b^BalExp with the
    %   channel balance b = |Pll - Prr|/(Pll + Prr).
    %   5.1/7.1 input: L/R are decomposed, C is added to the centre slot,
    %   the surround channels join the ambience objects, LFE is dropped.
    %   The experimental avendano_ring_separator enables a four-pair 5.1
    %   ring L/R, R/Rs, Rs/Ls, Ls/L while keeping C out of the output.
    %
    %   Refs: C. Avendano, J.-M. Jot, ICASSP 2002 (ambience), AES 22nd
    %   Conf. 2002, JAES 52(7/8) 2004 (upmix); C. Avendano, WASPAA 2003
    %   (panning index). See Documentation/Separation/AvendanoJot.
    %
    %   Analysis: sqrt-Hann, N = 2048, hop N/2 (perfect reconstruction);
    %   latency N - hop samples.

    properties
        K = 5               % front slots; in pair-bank mode, circular direct slots
        Width = 30          % output half-width of the front stage [deg]
        Base = 30           % stereo source base of the content [deg]
        AmbAz = 110         % ambience object azimuth [deg]
        AmbGain = 1
        Tau = 0.05          % statistics smoothing [s]
        N = 2048            % FFT size
        Phi0 = 0.7          % coherence threshold of the ambience mask (patent: 0.15
                            % at hop 256 / lambda 0.9; our coarser hop gives a
                            % higher coherence floor for diffuse sound)
        Sigma = 8           % slope of the ambience mask (patent: 8)
        Mu0 = 0.05          % mask floor (coherent bins; patent: 0.05)
        Mu1 = 1             % mask ceiling (incoherent bins)
        HardPanGuard = true % do not take hard-panned sources as ambience
        BalExp = 4          % exponent of the balance guard
        EnablePairBank = false % experimental multichannel ring mode
        PairBankMode = 'ring4' % 'ring4' excludes C; 'five_channel' includes C
        Diagnostics = false % keep per-frame direction histograms
    end

    properties (SetAccess = private)
        DiagAz = []         % [nAzBins x frames] direct energy per direction
        DiagPA = []         % [2 x frames] direct / ambient energy
        AzEdges = -90:2:90
        TrackLog = []       % (no tracking; kept for the panel's interface)
    end

    properties (Access = private)
        hop; win; frame; acc; pending; outq
        Pll; Prr; Plr
        Paa; Pbb; Pab
        slotAz; centre; layoutIdx
        pairBank = false
        pairCh = zeros(0, 2)
        pairAz = zeros(0, 2)
        pairGain = zeros(0, 2)
        pairDirAz = []
        pairMidAz = []
        pairNames = {}
        directCh = []
        directNames = {}
        directToChannels = false
        ambCh = []
        ambNames = {}
    end

    methods
        function obj = avendano_separator(varargin)
            for k = 1:2:numel(varargin), obj.(varargin{k}) = varargin{k+1}; end
        end

        function init(obj, fs, nIn)
            obj.fs = fs; obj.nIn = nIn;
            obj.K = max(2, round(obj.K));
            obj.hop = obj.N/2;
            obj.win = sqrt(0.5*(1 - cos(2*pi*(0:obj.N-1)'/obj.N)));   % periodic sqrt-Hann
            obj.Latency = obj.N - obj.hop;
            [az, names] = separator_base.channel_layout(nIn);
            obj.layoutIdx = struct('L', find(strcmp(names,'L')), 'R', find(strcmp(names,'R')), ...
                'C', find(strcmp(names,'C')), 'SL', find(ismember(names, {'Ls','Lb'})), ...
                'SR', find(ismember(names, {'Rs','Rb'})));
            if nIn == 1, obj.layoutIdx.L = 1; obj.layoutIdx.R = 1; end
            if nIn > 2 && isempty(obj.layoutIdx.L)               % unknown layout: first two
                obj.layoutIdx.L = 1; obj.layoutIdx.R = 2;
            end
            obj.setup_outputs(az, names);
            obj.reset();
        end

        function reset(obj)
            nb = obj.N/2 + 1;
            obj.frame = zeros(obj.N, obj.nIn);
            obj.acc = zeros(obj.N, obj.nOut());
            obj.pending = zeros(0, obj.nIn);
            obj.outq = zeros(0, obj.nOut());
            obj.Pll = 1e-9*ones(nb,1); obj.Prr = 1e-9*ones(nb,1); obj.Plr = zeros(nb,1);
            obj.Paa = 1e-9*ones(nb, size(obj.pairCh, 1));
            obj.Pbb = 1e-9*ones(nb, size(obj.pairCh, 1));
            obj.Pab = zeros(nb, size(obj.pairCh, 1));
            obj.DiagAz = []; obj.DiagPA = [];
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
            % decisions; every Xs{k} gets the same operations (oracle
            % evaluation, see stft_separator.separate_offline).
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
            % One STFT frame: decisions from the input frame, applied by
            % apply_ops to it and to optional shadow spectra Fsh (cell)
            I = obj.layoutIdx;
            F = fft(obj.frame .* obj.win, [], 1);
            nb = obj.N/2 + 1;
            F = F(1:nb,:);
            if obj.pairBank
                [S, th, e] = obj.synth_pair_bank(F);
                if obj.Diagnostics
                    [~, ~, bin] = histcounts(th, obj.AzEdges);
                    ok = bin > 0;
                    obj.DiagAz(:, end+1) = accumarray(bin(ok), e(ok), [numel(obj.AzEdges)-1, 1]);
                    obj.DiagPA(:, end+1) = [sum(e); sum(sum(abs(S(:, obj.K+1:end)).^2))];
                end
                ifs = @(S0) real(ifft([S0; conj(S0(end-1:-1:2,:))], [], 1)) .* obj.win;
                Yf = ifs(S);
                Ysh = {};
                if nargin > 1
                    Ysh = cellfun(@(Fs) ifs(obj.apply_pair_bank(Fs)), Fsh, 'UniformOutput', false);
                end
                return;
            end
            L = F(:, I.L); R = F(:, I.R);
            a = exp(-obj.hop/(obj.fs*obj.Tau));
            obj.Pll = a*obj.Pll + (1-a)*abs(L).^2;
            obj.Prr = a*obj.Prr + (1-a)*abs(R).^2;
            obj.Plr = a*obj.Plr + (1-a)*L.*conj(R);
            c = abs(obj.Plr);
            % ambience mask from the inter-channel coherence
            phi = min(1, c./sqrt(obj.Pll.*obj.Prr));
            gam = (obj.Mu1 + obj.Mu0)/2 + (obj.Mu1 - obj.Mu0)/2*tanh(obj.Sigma*pi*(obj.Phi0 - phi));
            if obj.HardPanGuard
                b = abs(obj.Pll - obj.Prr)./(obj.Pll + obj.Prr);
                gam = gam .* (1 - b.^obj.BalExp);
            end
            ops.gam = min(max(gam, 0), 1);
            % panning index (Avendano 2003), instantaneous: -1 = left .. +1 = right
            el = abs(L).^2 + 1e-12; er = abs(R).^2 + 1e-12; ci = abs(L).*abs(R);
            psi = 2*ci./(el + er);
            Psi = (1 - psi) .* sign(ci./el - ci./er);
            % -> unit-power panning vector [cos t, sin t], sin 2t = 1 - |Psi|
            t = asin(1 - abs(Psi))/2;                      % 0 .. pi/4 (hard .. centre)
            right = Psi > 0;
            t(right) = pi/2 - t(right);
            ops.a1 = cos(t); ops.a2 = sin(t);
            % direction: tangent law on the content's base, stretched to Width
            th = atand((ops.a1 - ops.a2)./(ops.a1 + ops.a2) * tand(obj.Base));
            th = th * obj.Width/obj.Base;
            % constant-power panning between neighbouring slots
            u = (obj.Width - th)/(2*obj.Width)*(obj.K - 1);
            u = min(max(u, 0), obj.K - 1);
            i0 = min(floor(u), obj.K - 2); f = u - i0;
            ops.W = zeros(nb, obj.K);
            ops.W(sub2ind([nb obj.K], (1:nb)', i0+1)) = cos(f*pi/2);
            ops.W(sub2ind([nb obj.K], (1:nb)', i0+2)) = sin(f*pi/2);
            S = obj.apply_ops(F, ops);
            if obj.Diagnostics
                e = sum(abs(S(:, 1:obj.K)).^2, 2);
                [~, ~, bin] = histcounts(th, obj.AzEdges);
                ok = bin > 0;
                obj.DiagAz(:, end+1) = accumarray(bin(ok), e(ok), [numel(obj.AzEdges)-1, 1]);
                obj.DiagPA(:, end+1) = [sum(e); sum(sum(abs(S(:, end-1:end)).^2))];
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
            gd = sqrt(1 - ops.gam.^2);                     % direct share
            DL = gd.*L; DR = gd.*R;
            m = ops.a1.*DL + ops.a2.*DR;                   % downmix onto the panning vector
            AL = ops.gam.*L + (DL - ops.a1.*m);            % ambience + downmix residual
            AR = ops.gam.*R + (DR - ops.a2.*m);
            S = [ops.W.*m, obj.AmbGain*AL, obj.AmbGain*AR];
            if ~isempty(I.C), S(:, obj.centre) = S(:, obj.centre) + F(:, I.C); end
            if ~isempty(I.SL), S(:, end-1) = S(:, end-1) + sum(F(:, I.SL), 2); end
            if ~isempty(I.SR), S(:, end)   = S(:, end)   + sum(F(:, I.SR), 2); end
        end

        function setup_outputs(obj, az, names)
            I = obj.layoutIdx;
            obj.pairBank = obj.EnablePairBank && ...
                numel(I.L) == 1 && numel(I.R) == 1 && numel(I.SL) == 1 && numel(I.SR) == 1;
            if obj.pairBank
                obj.directToChannels = strcmpi(obj.PairBankMode, 'five_channel') && numel(I.C) == 1;
                if obj.directToChannels
                    obj.pairCh = [I.L I.C; I.C I.R; I.R I.SR; I.SR I.SL; I.SL I.L];
                    obj.pairAz = [az(I.L) az(I.C); az(I.C) az(I.R); az(I.R) az(I.SR); ...
                        az(I.SR) az(I.SL); az(I.SL) az(I.L)];
                    obj.pairNames = {'L-C', 'C-R', 'R-Rs', 'Rs-Ls', 'Ls-L'};
                    obj.directCh = [I.L I.R I.C I.SL I.SR];
                    obj.directNames = names(obj.directCh);
                else
                    obj.pairCh = [I.L I.R; I.R I.SR; I.SR I.SL; I.SL I.L];
                    obj.pairAz = [az(I.L) az(I.R); az(I.R) az(I.SR); az(I.SR) az(I.SL); az(I.SL) az(I.L)];
                    obj.pairNames = {'L-R', 'R-Rs', 'Rs-Ls', 'Ls-L'};
                    obj.directCh = [];
                    obj.directNames = {};
                end
                if obj.directToChannels, obj.K = numel(obj.directCh); else, obj.K = size(obj.pairCh, 1); end
                obj.pairMidAz = obj.wrap_deg(obj.pairAz(:,1) + 0.5*obj.wrap_deg(obj.pairAz(:,2) - obj.pairAz(:,1)));
                obj.pairDirAz = obj.pairMidAz(:);
                obj.slotAz = obj.pairMidAz(:).';
                obj.ambCh = unique(obj.pairCh(:).', 'stable');
                obj.ambNames = cellfun(@(s) sprintf('Amb %s', s), obj.pairNames, 'UniformOutput', false);
                deg = zeros(1, max(obj.ambCh));
                for k = obj.ambCh, deg(k) = sum(obj.pairCh(:) == k); end
                obj.pairGain = zeros(size(obj.pairCh));
                for p = 1:size(obj.pairCh, 1)
                    obj.pairGain(p,1) = 1/sqrt(deg(obj.pairCh(p,1)));
                    obj.pairGain(p,2) = 1/sqrt(deg(obj.pairCh(p,2)));
                end
                if obj.directToChannels
                    obj.Names = [cellfun(@(s) sprintf('Direct %s', s), obj.directNames, 'UniformOutput', false), obj.ambNames];
                    directAz = az(obj.directCh);
                    obj.Pos = separator_base.dir2pos([directAz(:); obj.pairMidAz(:)]);
                else
                    obj.Names = [cellfun(@(s) sprintf('Dir %s', s), obj.pairNames, 'UniformOutput', false), obj.ambNames];
                    obj.Pos = separator_base.dir2pos([obj.pairDirAz(:); obj.pairMidAz(:)]);
                end
                obj.AzEdges = -180:2:180;
            else
                obj.slotAz = linspace(obj.Width, -obj.Width, obj.K);      % left -> right
                [~, obj.centre] = min(abs(obj.slotAz));
                obj.Names = [arrayfun(@(a) sprintf('Front %+.0f', a), obj.slotAz, 'UniformOutput', false), ...
                    {'Amb L', 'Amb R'}];
                obj.Pos = separator_base.dir2pos([obj.slotAz(:); obj.AmbAz; -obj.AmbAz]);
                obj.AzEdges = -90:2:90;
            end
        end

        function [S, thAll, eAll] = synth_pair_bank(obj, F)
            a = exp(-obj.hop/(obj.fs*obj.Tau));
            for p = 1:size(obj.pairCh, 1)
                A = obj.pairGain(p,1)*F(:, obj.pairCh(p,1));
                B = obj.pairGain(p,2)*F(:, obj.pairCh(p,2));
                obj.Paa(:,p) = a*obj.Paa(:,p) + (1-a)*abs(A).^2;
                obj.Pbb(:,p) = a*obj.Pbb(:,p) + (1-a)*abs(B).^2;
                obj.Pab(:,p) = a*obj.Pab(:,p) + (1-a)*A.*conj(B);
            end
            [S, thAll, eAll] = obj.apply_pair_bank(F);
        end

        function [S, thAll, eAll] = apply_pair_bank(obj, F)
            nb = obj.N/2 + 1;
            nPair = size(obj.pairCh, 1);
            if obj.directToChannels, nDirect = numel(obj.directCh); else, nDirect = nPair; end
            S = zeros(nb, nDirect + nPair);
            thAll = zeros(nb*size(obj.pairCh, 1), 1);
            eAll = zeros(size(thAll));
            row = 0;
            for p = 1:nPair
                A = obj.pairGain(p,1)*F(:, obj.pairCh(p,1));
                B = obj.pairGain(p,2)*F(:, obj.pairCh(p,2));
                ops = obj.pair_ops(A, B, obj.Paa(:,p), obj.Pbb(:,p), obj.Pab(:,p), obj.pairAz(p,:));
                gd = sqrt(1 - ops.gam.^2);
                DA = gd.*A; DB = gd.*B;
                m = ops.a1.*DA + ops.a2.*DB;
                AA = ops.gam.*A + (DA - ops.a1.*m);
                AB = ops.gam.*B + (DB - ops.a2.*m);
                if obj.directToChannels
                    ia = find(obj.directCh == obj.pairCh(p,1), 1);
                    ib = find(obj.directCh == obj.pairCh(p,2), 1);
                    S(:, ia) = S(:, ia) + ops.a1.*m;
                    S(:, ib) = S(:, ib) + ops.a2.*m;
                else
                    S(:, p) = m;
                end
                S(:, nDirect + p) = obj.AmbGain*(AA + AB)/sqrt(2);
                idx = row + (1:nb);
                thAll(idx) = ops.th;
                eAll(idx) = abs(m).^2;
                if sum(eAll(idx)) > 0
                    obj.pairDirAz(p) = atan2d(sum(eAll(idx).*sind(ops.th)), sum(eAll(idx).*cosd(ops.th)));
                else
                    obj.pairDirAz(p) = obj.pairMidAz(p);
                end
                row = row + nb;
            end
            if ~obj.directToChannels
                obj.Pos(1:nPair,:) = separator_base.dir2pos(obj.pairDirAz(:));
            end
        end

        function ops = pair_ops(obj, A, B, Paa, Pbb, Pab, azPair)
            c = abs(Pab);
            phi = min(1, c./sqrt(Paa.*Pbb));
            gam = (obj.Mu1 + obj.Mu0)/2 + (obj.Mu1 - obj.Mu0)/2*tanh(obj.Sigma*pi*(obj.Phi0 - phi));
            if obj.HardPanGuard
                bal = abs(Paa - Pbb)./(Paa + Pbb);
                gam = gam .* (1 - bal.^obj.BalExp);
            end
            ops.gam = min(max(gam, 0), 1);

            ea = abs(A).^2 + 1e-12; eb = abs(B).^2 + 1e-12; ci = abs(A).*abs(B);
            psi = 2*ci./(ea + eb);
            Psi = (1 - psi) .* sign(ci./ea - ci./eb);
            t = asin(1 - abs(Psi))/2;
            towardB = Psi > 0;
            t(towardB) = pi/2 - t(towardB);
            ops.a1 = cos(t); ops.a2 = sin(t);
            frac = t/(pi/2);
            ops.th = obj.wrap_deg(azPair(1) + frac*obj.wrap_deg(azPair(2) - azPair(1)));
        end
    end

    methods (Static, Access = private)
        function a = wrap_deg(a)
            a = mod(a + 180, 360) - 180;
        end
    end
end
