classdef headphone_db
    % HEADPHONE_DB  Headphone frequency responses from the AutoEq project
    % (github.com/jaakkopasanen/AutoEq, MIT; measurements by oratory1990,
    % Rtings, crinacle, Innerfidelity and ~20 more sources: over-ear,
    % in-ear, earbuds, Bluetooth models included) and the headphone
    % equalisation of the Binauralizer's playback (playback_calibration).
    %
    %   T = headphone_db.index()          the measurement list (cached in
    %       Data/Headphones/autoeq_index.mat; refresh: index(true)): table
    %       source, form ('over-ear' | 'in-ear' | 'earbud'), rig (e.g.
    %       'Bruel & Kjaer 5128', 'HMS II.3', '' = the source's own), name,
    %       path (in the repository)
    %   T = headphone_db.search(T, 'airpods pro', form)   every word must match
    %   [row, score, model] = headphone_db.match(T, deviceName)   the
    %       measurement of the headphone an audio output device is (e.g.
    %       'Fejhallgató (Technics EAH-AZ100)', 'Headphones (WH-1000XM4
    %       Stereo)'): the generic words of the device name are dropped
    %       (Headphones, Stereo, Hands-Free, ...; Hungarian too), every
    %       remaining model word with a digit must be in the measurement's
    %       name; best: most of the words matched, then the shortest name
    %       (the plain model, not a variant), then the source (oratory1990,
    %       crinacle, Rtings, ...). row: one table row or empty; score 0..1
    %   M = headphone_db.fetch(T(i, :))   f [Hz], raw [dB] (frequency,raw CSV;
    %       the response at the rig's eardrum / ear-simulator microphone),
    %       name, source, rig, form - cached in Data/Headphones/cache
    %   [f, D] = headphone_db.target(rig, form)   AutoEq's diffuse-field target
    %       of the rig (Diffuse field 5128 / GRAS KEMAR / ISO 11904-1)
    %   E = headphone_db.eq(M, opt)       the correction [dB] on f:
    %       opt.mode 'eardrum' (default: -raw - a flat response at the
    %       eardrum, which is what the Binauralizer's HRTFs assume: measured
    %       with a head-and-torso simulator's ear simulator) | 'diffuse'
    %       (diffuse-field target - raw); smooth 1/6 [octave], maxBoost 10,
    %       | 'hrtf' (binaural playback with an HRTF set: the rig's
    %       diffuse-field target - the set's diffuse-field average
    %       opt.hrtfDF (hrtf_df_average) - raw: the playback gives in a
    %       diffuse field what the rig's ear would get; an eardrum set ->
    %       'eardrum', a diffuse-field equalised set -> 'diffuse', a
    %       blocked-ear-canal set -> its canal resonance restored);
    %       maxCut 15 [dB], fLow 20, fHigh 10000 [Hz] (above: faded to 0 -
    %       placement-sensitive, unreliable for a correction); normalised to
    %       0 dB at 1 kHz (the level calibration stays valid)
    %   Network: two GitHub API calls for the index, one raw download per
    %   headphone (webread / websave).
    %   Sensitivity / impedance are not in AutoEq (manufacturer data,
    %   Reference Audio Analyzer, Stereophile): the level calibration takes
    %   them as input, or a reference tone + SPL meter (always needed for
    %   Bluetooth headphones: their level is set digitally).

    properties (Constant)
        Repo = 'jaakkopasanen/AutoEq'
        Branch = 'master'
    end

    methods (Static)
        function d = dataDir()
            d = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))), 'Data', 'Headphones');
            if ~isfolder(d), mkdir(d); end
        end

        function T = index(refresh)
            if nargin < 1, refresh = false; end
            f = fullfile(headphone_db.dataDir(), 'autoeq_index.mat');
            if ~refresh && isfile(f), L = load(f, 'T'); T = L.T; return; end
            o = weboptions('Timeout', 60, 'ContentType', 'json');
            api = ['https://api.github.com/repos/' headphone_db.Repo];
            root = webread([api '/contents/?ref=' headphone_db.Branch], o);
            if iscell(root), root = [root{:}]; end
            sha = root(strcmp({root.name}, 'measurements')).sha;
            tr = webread(sprintf('%s/git/trees/%s?recursive=1', api, sha), o);
            if iscell(tr.tree), P = cellfun(@(x) x.path, tr.tree, 'UniformOutput', false);   % mixed entries (blob / tree)
            else, P = {tr.tree.path}; end
            P = P(endsWith(P, '.csv') & contains(P, '/data/'));
            n = numel(P); src = cell(n, 1); form = src; rig = src; name = src;
            for i = 1:n
                s = strsplit(P{i}, '/');                       % source/data/form[/rig]/name.csv
                src{i} = s{1}; form{i} = s{3}; name{i} = erase(s{end}, '.csv');
                if numel(s) >= 5, rig{i} = s{4}; else, rig{i} = ''; end
            end
            T = table(string(src), string(form), string(rig), string(name), string(P(:)), ...
                'VariableNames', {'source', 'form', 'rig', 'name', 'path'});
            T = sortrows(T, {'name', 'source'});
            save(f, 'T');
        end

        function R = search(T, q, form)
            % rows whose name (or source) contains every word of q
            if nargin < 3, form = ''; end
            w = strsplit(lower(strtrim(q)));
            ok = true(height(T), 1); hay = lower(T.name + " " + T.source);
            for k = 1:numel(w), if ~isempty(w{k}), ok = ok & contains(hay, w{k}); end, end
            if ~isempty(form), ok = ok & T.form == string(form); end
            R = T(ok, :);
        end

        function [row, score, model] = match(T, dev)
            row = T([], :); score = 0;
            d = char(dev);
            t = regexp(d, '\(([^()]*)\)\s*$', 'tokens', 'once');      % the part in the last brackets
            if ~isempty(t), d = t{1}; end
            d = regexprep(d, '^\s*\d+\s*-\s*', '');                   % '2- ' prefixes
            model = strtrim(d);
            stop = {'headphones', 'headphone', 'headset', 'earphones', 'earbuds', 'stereo', 'hands', 'free', ...
                'handsfree', 'hands-free', 'ag', 'audio', 'hf', 'bluetooth', 'bt', 'le', 'a2dp', 'sound', 'device', ...
                'speakers', 'speaker', 'output', 'realtek', 'r', 'high', 'definition', 'intel', 'nvidia', 'amd', ...
                'display', 'digital', 'usb', 'generic', 'default', 'primary', 'driver', 'monitor', 'hd', 'for', ...
                'fejhallgató', 'fejhallgato', 'fejbeszélő', 'fejbeszelo', 'hangszórók', 'hangszóró', 'hangszorok', ...
                'elsődleges', 'hangillesztő', 'mikrofon', 'fülhallgató', 'fulhallgato'};
            w = regexp(lower(model), '[a-z0-9áéíóöőúüű]+', 'match');
            w = w(~ismember(w, stop));
            if isempty(w), return; end
            tok = [w, {[w{:}]}];                                         % + the joined form ('eahaz100')
            hasDigit = ~cellfun(@isempty, regexp(w, '\d', 'once'));
            if ~any(hasDigit) && max(cellfun(@numel, w)) < 5, return; end   % nothing model-like
            nm = regexprep(lower(T.name), '[^a-z0-9]', '');
            hit = false(height(T), numel(tok));
            for k = 1:numel(tok), hit(:, k) = contains(nm, regexprep(tok{k}, '[^a-z0-9]', '')); end
            need = hit(:, [hasDigit, false]);
            ok = all(need, 2) & (any(hit(:, 1:end-1), 2) | hit(:, end));
            if ~any(ok), return; end
            len = cellfun(@numel, w);
            sc = (hit(:, 1:end-1)*len(:))/sum(len);
            sc(hit(:, end)) = max(sc(hit(:, end)), 1);                   % the joined model word matches
            sc(~ok) = 0;
            best = max(sc);
            if best < 0.6, return; end
            c = find(sc == best);
            srcRank = @(s) find(strcmpi(s, {'oratory1990', 'crinacle', 'Rtings', 'Innerfidelity', 'Headphone.com Legacy'}), 1);
            r = zeros(numel(c), 1);
            for i = 1:numel(c)
                q = srcRank(char(T.source(c(i)))); if isempty(q), q = 9; end
                r(i) = strlength(T.name(c(i)))*10 + q;
            end
            [~, i] = min(r); row = T(c(i), :); score = best;
        end

        function M = fetch(row)
            p = char(row.path);
            c = fullfile(headphone_db.dataDir(), 'cache', strrep(p, '/', filesep));
            if ~isfile(c)
                if ~isfolder(fileparts(c)), mkdir(fileparts(c)); end
                websave(c, headphone_db.rawUrl(['measurements/' p]), weboptions('Timeout', 60));
            end
            [f, raw] = headphone_db.readCsv(c);
            M = struct('f', f, 'raw', raw, 'name', char(row.name), 'source', char(row.source), ...
                'rig', char(row.rig), 'form', char(row.form), 'file', c);
        end

        function [f, D, name] = target(rig, form)
            % the rig's diffuse-field target (AutoEq targets/)
            rig = lower(char(rig));
            if contains(rig, '5128'), name = 'Diffuse field 5128';
            elseif contains(rig, 'kemar') || contains(rig, '711') || contains(rig, 'gras'), name = 'Diffuse field GRAS KEMAR';
            else, name = 'Diffuse field ISO 11904-1';
            end
            if nargin > 1 && strcmp(form, 'in-ear') && contains(rig, '5128'), name = 'Diffuse field 5128'; end
            c = fullfile(headphone_db.dataDir(), 'cache', 'targets', [name '.csv']);
            if ~isfile(c)
                if ~isfolder(fileparts(c)), mkdir(fileparts(c)); end
                websave(c, headphone_db.rawUrl(['targets/' name '.csv']), weboptions('Timeout', 60));
            end
            [f, D] = headphone_db.readCsv(c);
        end

        function E = eq(M, opt)
            if nargin < 2, opt = struct(); end
            g = @(n, d) headphone_db.getopt(opt, n, d);
            mode = g('mode', 'eardrum');
            f = M.f(:); r = M.raw(:);
            switch mode
                case 'diffuse'
                    [ft, D] = headphone_db.target(M.rig, M.form);
                    C = interp1(log(ft), D, log(f), 'linear', 'extrap') - r;
                case 'hrtf'
                    [ft, D] = headphone_db.target(M.rig, M.form); H = opt.hrtfDF;
                    C = interp1(log(ft), D, log(f), 'linear', 'extrap') - r ...
                        - interp1(log(max(H.f(2:end), 1)), H.dB(2:end), log(f), 'linear', 'extrap');
                otherwise
                    C = -r;
            end
            C = headphone_db.smooth(f, C, g('smooth', 1/6));
            C = C - interp1(log(f), C, log(1000));                  % 0 dB at 1 kHz
            C = min(max(C, -g('maxCut', 15)), g('maxBoost', 10));
            fH = g('fHigh', 10000); fL = g('fLow', 20);
            w = ones(size(f));
            hi = f > fH; w(hi) = max(0, 1 - log2(f(hi)/fH)/0.5);     % half-octave fade above fHigh
            lo = f < fL; w(lo) = 0;
            C = C.*w;
            E = struct('f', f, 'dB', C, 'mode', mode, 'name', M.name, 'source', M.source, 'rig', M.rig, 'opt', opt);
        end
    end

    methods (Static, Access = private)
        function u = rawUrl(p)
            u = sprintf('https://raw.githubusercontent.com/%s/%s/%s', headphone_db.Repo, headphone_db.Branch, ...
                strjoin(cellfun(@(s) urlencodeSeg(s), strsplit(p, '/'), 'UniformOutput', false), '/'));
            function e = urlencodeSeg(s)
                e = char(java.net.URLEncoder.encode(s, 'UTF-8')); e = strrep(e, '+', '%20');
            end
        end

        function [f, v] = readCsv(file)
            A = readmatrix(file, 'NumHeaderLines', 1, 'Delimiter', ',');
            A = A(all(isfinite(A(:, 1:2)), 2), :);
            f = A(:, 1); v = A(:, 2);
        end

        function Ls = smooth(f, L, bw)
            if ~(bw > 0), Ls = L; return; end
            Ls = L; lf = log2(f);
            for k = 1:numel(f)
                in = abs(lf - lf(k)) <= bw/2;
                Ls(k) = mean(L(in));
            end
        end

        function v = getopt(s, n, d)
            if isfield(s, n) && ~isempty(s.(n)), v = s.(n); else, v = d; end
        end
    end
end
