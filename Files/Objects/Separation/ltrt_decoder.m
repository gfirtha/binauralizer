classdef ltrt_decoder < separator_base
    %LTRT_DECODER  Frequency-domain active Lt/Rt matrix decoder:
    %stereo -> L C R Ls Rs, steered per critical band (a 4-2-4 matrix
    %decoder in the spirit of the classic film surround decoders, without
    %their trademarks). Only for 2-channel input; other widths pass through.
    %
    %   Per bin, mid / side:  M = (L + R)/sqrt2,  D = (L - R)/sqrt2.
    %   Per band, from smoothed powers (TauSteer):
    %     lr = (|L|^2 - |R|^2)/(|L|^2 + |R|^2)   left / right dominance
    %     cs = (|M|^2 - |D|^2)/(|M|^2 + |D|^2)   centre / surround dominance
    %   Centre  C   = gC M,  gC = max(cs, 0)^Steep
    %   Surround Sur = gS D,  gS = max(-cs, 0)^Steep  (Lt/Rt surround is the
    %   anti-phase part; its +-90 deg encoding phase is a common rotation)
    %   Fronts  L' = L - (C + Sur)/sqrt2,  R' = R - (C - Sur)/sqrt2
    %           (so the passive re-encode L' + (C+Sur)/sqrt2 gives L back)
    %   Surround pair IN PHASE, steered by lr (power complementary):
    %           Ls = Sur sqrt((1 + lr)/2),  Rs = Sur sqrt((1 - lr)/2)
    %   - never an anti-phase Ls/Rs pair: at the listener (equidistant
    %   loudspeakers) that would be a pressure node where the local
    %   wavenumber is undefined; and no surround delay (it would fake
    %   wavefront curvature in the wavenumber analysis).
    %
    %   As a separator: N-point sqrt-Hann STFT, outputs L C R Ls Rs at
    %   30 0 -30 110 -110 deg. decode_spectrum() is used inside other
    %   STFT separators (same frames, no extra latency).

    properties
        N = 1024
        TauSteer = 0.02     % [s] steering statistics smoothing
        Steep = 1           % steering exponent (> 1: gentler, < 1: more aggressive)
        BandEdges = [0 100 200 300 400 510 630 770 920 1080 1270 1480 1720 2000 2320 ...
                     2700 3150 3700 4400 5300 6400 7700 9500 12000 15500 inf]
    end

    properties (Access = private)
        hop; win; frame; acc; pending; outq
        bandOf; M; nB
        pL; pR; pM; pD                  % smoothed band powers
        pass = []
    end

    methods
        function obj = ltrt_decoder(varargin)
            for k = 1:2:numel(varargin), obj.(varargin{k}) = varargin{k+1}; end
        end

        function init(obj, fs, nIn)
            obj.fs = fs; obj.nIn = nIn;
            if nIn ~= 2
                obj.pass = passthrough_separator(); obj.pass.init(fs, nIn);
                obj.Names = obj.pass.Names; obj.Pos = obj.pass.Pos; obj.Latency = 0;
                return;
            end
            obj.pass = [];
            obj.Names = {'L','C','R','Ls','Rs'};
            obj.Pos = separator_base.dir2pos([30 0 -30 110 -110]);
            obj.hop = obj.N/2;
            obj.win = sqrt(0.5*(1 - cos(2*pi*(0:obj.N-1)'/obj.N)));
            obj.Latency = obj.N - obj.hop;
            obj.setup_bands(obj.N, fs);
            obj.reset();
        end

        function setup_bands(obj, N, fs)
            % band map for an N-point spectrum (also when used inside
            % another separator's STFT)
            nb = N/2 + 1; f = (0:nb-1)'*fs/N;
            [~, ~, b] = histcounts(f, obj.BandEdges);
            used = unique(b);
            obj.bandOf = arrayfun(@(x) find(used == x), b);
            obj.nB = numel(used);
            obj.M = sparse(obj.bandOf, (1:nb)', 1, obj.nB, nb);
            obj.fs = fs;
            obj.reset_steering();
        end

        function reset_steering(obj)
            obj.pL = []; obj.pR = []; obj.pM = []; obj.pD = [];
        end

        function reset(obj)
            if ~isempty(obj.pass), obj.pass.reset(); return; end
            obj.frame = zeros(obj.N, 2);
            obj.acc = zeros(obj.N, 5);
            obj.pending = zeros(0, 2);
            obj.outq = zeros(0, 5);
            obj.reset_steering();
        end

        function Y = decode_spectrum(obj, L, R, dt)
            % L, R: [bins x 1] spectra of one frame (bins matching
            % setup_bands); dt = frame step [s]. Returns [bins x 5]:
            % L' C R' Ls Rs.
            Mid = (L + R)/sqrt(2); Sid = (L - R)/sqrt(2);
            B = obj.M;
            iL = B*abs(L).^2; iR = B*abs(R).^2; iM = B*abs(Mid).^2; iD = B*abs(Sid).^2;
            a = exp(-dt/obj.TauSteer);
            if isempty(obj.pL)
                obj.pL = iL; obj.pR = iR; obj.pM = iM; obj.pD = iD;
            else
                obj.pL = a*obj.pL + (1-a)*iL; obj.pR = a*obj.pR + (1-a)*iR;
                obj.pM = a*obj.pM + (1-a)*iM; obj.pD = a*obj.pD + (1-a)*iD;
            end
            lr = (obj.pL - obj.pR)./max(obj.pL + obj.pR, eps);
            cs = (obj.pM - obj.pD)./max(obj.pM + obj.pD, eps);
            gC = max(cs, 0).^obj.Steep; gS = max(-cs, 0).^obj.Steep;
            gC = gC(obj.bandOf); gS = gS(obj.bandOf); lr = lr(obj.bandOf);
            C = gC.*Mid; Sur = gS.*Sid;
            Lf = L - (C + Sur)/sqrt(2);
            Rf = R - (C - Sur)/sqrt(2);
            Y = [Lf, C, Rf, Sur.*sqrt((1 + lr)/2), Sur.*sqrt((1 - lr)/2)];
        end

        function Y = process(obj, X)
            if ~isempty(obj.pass), Y = obj.pass.process(X); return; end
            n = size(X,1);
            obj.pending = [obj.pending; X];
            nb = obj.N/2 + 1;
            while size(obj.pending,1) >= obj.hop
                obj.frame = [obj.frame(obj.hop+1:end,:); obj.pending(1:obj.hop,:)];
                obj.pending = obj.pending(obj.hop+1:end,:);
                F = fft(obj.frame .* obj.win, [], 1);
                S = obj.decode_spectrum(F(1:nb,1), F(1:nb,2), obj.hop/obj.fs);
                yf = real(ifft([S; conj(S(end-1:-1:2,:))], [], 1)) .* obj.win;
                obj.acc = obj.acc + yf;
                obj.outq = [obj.outq; obj.acc(1:obj.hop,:)];
                obj.acc = [obj.acc(obj.hop+1:end,:); zeros(obj.hop, 5)];
            end
            m = min(n, size(obj.outq,1));
            Y = [zeros(n-m, 5); obj.outq(1:m,:)];
            obj.outq = obj.outq(m+1:end,:);
        end
    end
end
