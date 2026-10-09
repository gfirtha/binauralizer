classdef matrix_separator < separator_base
    %MATRIX_SEPARATOR  Broadband matrix upmix of stereo (Lt/Rt) to L C R Ls Rs.
    %
    %   Passive Dolby-Surround-style matrix with correlation steering:
    %     rho  = smoothed normalised correlation of L and R (broadband)
    %     C    = max(rho,0) * (L+R)/2,   L' = L - C,   R' = R - C
    %     S    = (L-R)/2 * (1 - max(rho,0)), delayed (SurroundDelay),
    %            Ls = S, Rs = -S (decorrelated pair)
    %   A baseline for the STFT upmix: zero spectral latency, no
    %   frequency resolution. Other input widths are passed through.

    properties
        Tau = 0.2                   % correlation smoothing [s]
        SurroundDelay = 0.012       % [s]
    end

    properties (Access = private)
        sLL = 1e-9; sRR = 1e-9; sLR = 0
        dline = []
        pass = []
    end

    methods
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
            obj.Latency = 0;
            obj.reset();
        end

        function reset(obj)
            obj.sLL = 1e-9; obj.sRR = 1e-9; obj.sLR = 0;
            obj.dline = zeros(round(obj.SurroundDelay*obj.fs), 1);
        end

        function Y = process(obj, X)
            if ~isempty(obj.pass), Y = obj.pass.process(X); return; end
            L = X(:,1); R = X(:,2); n = size(X,1);
            a = exp(-n/(obj.fs*obj.Tau));
            obj.sLL = a*obj.sLL + (1-a)*mean(L.^2);
            obj.sRR = a*obj.sRR + (1-a)*mean(R.^2);
            obj.sLR = a*obj.sLR + (1-a)*mean(L.*R);
            rho = max(0, obj.sLR / sqrt(obj.sLL*obj.sRR + eps));
            C = rho*(L + R)/2;
            S = (1 - rho)*(L - R)/2;
            buf = [obj.dline; S];
            Sd = buf(1:n);
            obj.dline = buf(n+1:end);
            Y = [L - C, C, R - C, Sd, -Sd];
        end
    end
end
