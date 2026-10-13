classdef speaker_eq < handle
    % SPEAKER_EQ  Per-loudspeaker equalization of the SSD driving signals
    % (the loudspeakers' own Pre-EQ, e.g. from the Loudspeaker designer).
    %   Applied by sound_scene_renderer after the driving functions, on every
    %   loudspeaker's signal: it reaches the sound-card outputs, the
    %   binaural rendering of the loudspeakers, and (folded into the system
    %   BRIRs) the fast binaural path.
    %
    %   E = speaker_eq(data, fs, Nls, blk)
    %     data   struct (a speaker-EQ file, speaker_eq.read):
    %            fs    the design rate [Hz]
    %            sos   biquad cascade [nb x 6] rows [b0 b1 b2 1 a1 a2] (IIR,
    %                  exact), or a cell of them per loudspeaker
    %            peq   (optional) the PEQ design (fc, gain, Q per entry) and hp
    %                  ([order fc], Butterworth): redesigned at fs when it differs
    %            hinv  FIR taps [N x nCh] (used when there is no sos; the
    %                  legacy hinv.mat format)
    %            name, source (text)
    %            One filter (1 column / one cascade) is used for every
    %            loudspeaker; nCh >= Nls: column n for loudspeaker n (SSD
    %            order); fewer: the rest unfiltered.
    %     fs, Nls, blk  the scene's rate, loudspeaker count, block size
    %   Y = E.process(X)        X [n x Nls] one block, state kept per block
    %   G = E.response(f)       [nf x Nls] complex
    %   h = E.impulse(N)        [N x Nls] impulse responses (N = [] : E.irLength)

    properties
        fs
        Nls
        kind = 'sos'            % 'sos' | 'fir'
        sos = {}                % {Nls} [nb x 6] (empty: unfiltered)
        fir = []                % [N x Nls]
        name = ''
        data = []               % the loaded struct (Setup.SpeakerEQ)
        irLength = 0            % samples where the impulse responses have decayed (-80 dB)
    end

    properties (Access = private)
        zi = {}                 % {Nls} [2 x nb] biquad states
        conv = {}               % {Nls} OLS_convolver (FIR)
    end

    methods
        function obj = speaker_eq(data, fs, Nls, blk)
            obj.fs = fs; obj.Nls = Nls; obj.data = data;
            if isfield(data, 'name'), obj.name = data.name; end
            fsD = fs; if isfield(data, 'fs') && ~isempty(data.fs), fsD = data.fs; end
            S = {};
            if isfield(data, 'peq') && fsD ~= fs && ~isempty(data.peq)
                S = {speaker_eq.design(data.peq, getf(data, 'hp', []), fs)};          % exact at the new rate
            elseif isfield(data, 'sos') && ~isempty(data.sos)
                S = data.sos; if ~iscell(S), S = {S}; end
                if fsD ~= fs
                    error('speaker_eq:fs', 'The biquads were designed at %g Hz, the scene runs at %g Hz (no PEQ parameters to redesign).', fsD, fs);
                end
            end
            if ~isempty(S)
                obj.kind = 'sos';
                obj.sos = cell(1, Nls);
                for n = 1:Nls
                    if numel(S) == 1, obj.sos{n} = S{1}; elseif n <= numel(S), obj.sos{n} = S{n}; end
                end
                obj.zi = cellfun(@(s) zeros(2, size(s, 1)), obj.sos, 'UniformOutput', false);
            elseif isfield(data, 'hinv') && ~isempty(data.hinv)
                h = data.hinv;
                if fsD ~= fs
                    if exist('resample', 'file') == 2, h = resample(h, fs, fsD);
                    else, error('speaker_eq:fs', 'The FIR was designed at %g Hz, the scene runs at %g Hz.', fsD, fs); end
                end
                obj.kind = 'fir';
                if size(h, 2) == 1, h = repmat(h, 1, Nls); end
                obj.fir = zeros(size(h, 1), Nls);
                m = min(size(h, 2), Nls); obj.fir(:, 1:m) = h(:, 1:m);
                if m < Nls, obj.fir(1, m+1:end) = 1; end                         % the rest unfiltered
                obj.conv = cell(1, Nls);
                for n = 1:Nls, obj.conv{n} = OLS_convolver(obj.fir(:, n), blk); end
            else
                error('speaker_eq:empty', 'The speaker EQ has neither biquads (sos) nor FIR taps (hinv).');
            end
            h = obj.impulse(2^16); e = max(abs(h), [], 2);
            obj.irLength = max(find(e > 1e-4*max(e), 1, 'last'), 1);
        end

        function Y = process(obj, X)
            Y = X;
            for n = 1:min(obj.Nls, size(X, 2))
                switch obj.kind
                    case 'sos'
                        s = obj.sos{n}; if isempty(s), continue; end
                        y = X(:, n); z = obj.zi{n};
                        for b = 1:size(s, 1), [y, z(:, b)] = filter(s(b, 1:3), s(b, 4:6), y, z(:, b)); end
                        obj.zi{n} = z; Y(:, n) = y;
                    case 'fir'
                        Y(:, n) = obj.conv{n}.convolve(X(:, n));
                end
            end
        end

        function G = response(obj, f)
            f = f(:); G = ones(numel(f), obj.Nls); z = exp(-1i*2*pi*f/obj.fs);
            for n = 1:obj.Nls
                switch obj.kind
                    case 'sos'
                        s = obj.sos{n};
                        for b = 1:size(s, 1)
                            G(:, n) = G(:, n).*(s(b, 1) + s(b, 2)*z + s(b, 3)*z.^2)./(s(b, 4) + s(b, 5)*z + s(b, 6)*z.^2);
                        end
                    case 'fir'
                        G(:, n) = polyval(flipud(obj.fir(:, n)), z);
                end
            end
        end

        function h = impulse(obj, N)
            if nargin < 2 || isempty(N), N = obj.irLength; end
            d = [1; zeros(N - 1, 1)]; h = zeros(N, obj.Nls);
            for n = 1:obj.Nls
                switch obj.kind
                    case 'sos'
                        y = d; s = obj.sos{n};
                        for b = 1:size(s, 1), y = filter(s(b, 1:3), s(b, 4:6), y); end
                        h(:, n) = y;
                    case 'fir'
                        m = min(N, size(obj.fir, 1)); h(1:m, n) = obj.fir(1:m, n);
                end
            end
        end

        function reset(obj)
            obj.zi = cellfun(@(z) 0*z, obj.zi, 'UniformOutput', false);
        end
    end

    methods (Static)
        function data = read(file)
            % a speaker-EQ .mat (Loudspeaker designer > Pre-EQ > Export, or a
            % legacy hinv.mat with hinv + fs)
            L = load(file);
            if isfield(L, 'speakerEQ'), data = L.speakerEQ;
            elseif isfield(L, 'hinv'), data = struct('fs', L.fs, 'hinv', L.hinv, 'sos', [], 'name', file);
            else, error('speaker_eq:file', '%s holds no speakerEQ struct and no hinv taps.', file);
            end
            if ~isfield(data, 'name') || isempty(data.name), [~, data.name] = fileparts(file); end
            data.file = file;
        end

        function sos = design(peq, hp, fs)
            % biquads of a PEQ design (RBJ peaking, as eq_design) + a
            % Butterworth high-pass ([order fc]) at fs
            sos = zeros(0, 6);
            if ~isempty(hp)
                n = hp(1); K = tan(pi*min(hp(2), 0.499*fs)/fs);
                for k = 1:floor(n/2)
                    Q = 1/(2*sin(pi*(2*k - 1)/(2*n))); D = 1 + K/Q + K^2;
                    sos(end+1, :) = [[1 -2 1]/D, 1, 2*(K^2 - 1)/D, (1 - K/Q + K^2)/D]; %#ok<AGROW>
                end
                if mod(n, 2), sos(end+1, :) = [[1 -1 0]/(K + 1), 1, (K - 1)/(K + 1), 0]; end
            end
            for i = 1:numel(peq)
                A = 10^(peq(i).gain/40); w0 = 2*pi*min(peq(i).fc, 0.49*fs)/fs; al = sin(w0)/(2*peq(i).Q);
                b = [1 + al*A, -2*cos(w0), 1 - al*A]; a = [1 + al/A, -2*cos(w0), 1 - al/A];
                sos(end+1, :) = [b/a(1), 1, a(2:3)/a(1)]; %#ok<AGROW>
            end
        end
    end
end

function v = getf(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
