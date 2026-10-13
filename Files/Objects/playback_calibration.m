classdef playback_calibration < handle
    % PLAYBACK_CALIBRATION  Calibrated headphone output of the real-time
    % Binauralizer: the ear signals in pascal are played at the simulated
    % SPL, or - when that is too loud or a quieter level is chosen - quieter
    % with the timbre of the simulated level (ISO 226 loudness compensation).
    %
    %   Calibration (the same fields as nl_playback, Loudspeaker designer):
    %     mode 'headphone': sens [dB SPL] per sensUnit 'mW' | 'V', Z [Ohm],
    %       Vfs [V rms of a full-scale sine at the output, with the
    %       computer's volume as used]
    %     mode 'reference': Lref [dB SPL measured at the ear / coupler] for the
    %       1 kHz test tone at refDBFS [dBFS] (testTone)
    %   -> SPLfs: the SPL of a full-scale sine; digital = p * kPa with
    %      kPa = (1/sqrt(2)) / (20e-6 * 10^(SPLfs/20)).
    %   Level: a running Leq of the ear pressure (time constant tau [s],
    %   unweighted, both ears' power mean) -> Lsim. Attenuation att = max(
    %   offset, Lsim - maxL, peak headroom); the compensation (when att > 0.5
    %   dB and compensate) is the ISO 226 contour difference between the
    %   played and the simulated loudness level (0 dB at 1 kHz), as a
    %   minimum-phase FIR (nTaps) re-designed when the level changes by more
    %   than 1 dB, cross-faded over one block. A pure-tone contour applied to
    %   a broadband signal at its Leq is an approximation; the headphone's own
    %   free / diffuse-field equalisation is not included.
    %   process(Pa) is the per-block stage: [nBlock x 2] pascal -> [nBlock x 2]
    %   digital, clipped at -headroom dBFS (a sample-peak limiter with
    %   release, never hard clipping).
    %   Headphone EQ (hpEq, headphone_db.eq: the chosen headphone's
    %   correction, 0 dB at 1 kHz) joins the compensation in the same
    %   minimum-phase FIR (setHeadphoneEq; 1024 taps with an EQ).
    %   eqDirect(x, fs, diffuse): the headphone EQ alone (no level
    %   calibration) for any other stereo output - direct playback with the
    %   binauralization off (music: the diffuse-field EQ hpEqDirect, the ear's
    %   own resonances kept), or uncalibrated binaural playback (the eardrum
    %   EQ hpEq); own minimum-phase FIR, cross-faded on a change, the same
    %   sample-peak limiter.
    %   Persisted in prefdir (playback_calibration.mat): load() / save().

    properties
        mode = 'headphone'
        sens = 100              % dB SPL per sensUnit
        sensUnit = 'mW'
        Z = 32                  % Ohm
        Vfs = 1.0               % V rms of a full-scale sine
        Lref = 80               % reference mode: dB SPL of the test tone
        refDBFS = -20           % reference mode: the test tone's level
        maxL = 95               % dB SPL Leq cap (hearing safety)
        offset = 0              % dB deliberate attenuation
        compensate = true       % ISO 226 when attenuated
        headroom = 1            % dB below full scale
        tau = 1.0               % s, Leq time constant
        nTaps = 256             % compensation FIR
        enabled = true          % false: the ear pressure is scaled 1:1 to digital (legacy, uncalibrated)
        hpEq = []               % headphone equalisation (headphone_db.eq: f, dB, name, ...) or []
        hpEqDirect = []         % its diffuse-field version for direct (non-binaural) playback, [] = hpEq
        hrtfName = ''           % the HRTF set chosen with it (display / persistence)
    end

    properties (SetAccess = private)
        fs = 48000
        Lsim = -Inf             % running Leq of the simulated ear pressure [dB SPL]
        Lplay = -Inf            % what is played [dB SPL]
        att = 0                 % dB
        peak = 0                % dBFS of the last block (after the stage)
        reason = ''
        limGain = 1             % peak limiter state
    end

    properties (Access = private)
        pw = 0                  % running mean square [Pa^2]
        h = 1                   % current compensation FIR
        hNew = []               % FIR to cross-fade to
        zi = []                 % FIR state [nTaps-1 x 2]
        ziNew = []
        Ldesign = NaN           % (Lsim, att) of the current FIR
        attDesign = 0
        compDesign = false      % the current FIR has the loudness compensation
        eqVer = 0               % bumped by setHeadphoneEq
        eqDesign = -1
        hD = 1                  % eqDirect: current FIR, the one to fade to, states, design key
        hDNew = []
        ziD = []
        ziDNew = []
        dKey = ''
    end

    methods
        function obj = playback_calibration(varargin)
            for k = 1:2:numel(varargin), obj.(varargin{k}) = varargin{k+1}; end
        end

        function s = SPLfs(obj)
            % SPL of a full-scale sine [dB]
            switch obj.mode
                case 'headphone'
                    sV = obj.sens; if strcmpi(obj.sensUnit, 'mW'), sV = obj.sens + 10*log10(1000/obj.Z); end
                    s = sV + 20*log10(obj.Vfs);
                otherwise
                    s = obj.Lref - obj.refDBFS;
            end
        end

        function k = kPa(obj)
            % digital units per pascal
            k = (1/sqrt(2))/(20e-6*10^(obj.SPLfs()/20));
        end

        function reset(obj, fs)
            if nargin > 1, obj.fs = fs; end
            obj.pw = 0; obj.h = 1; obj.hNew = []; obj.zi = []; obj.ziNew = []; obj.Ldesign = NaN; obj.attDesign = 0;
            obj.Lsim = -Inf; obj.Lplay = -Inf; obj.att = 0; obj.limGain = 1; obj.peak = -Inf;
        end

        function y = process(obj, p, fs)
            % one block: ear pressure [Pa] (n x 2) -> digital (n x 2)
            if nargin > 2 && fs ~= obj.fs, obj.reset(fs); end
            if ~obj.enabled, y = p; obj.peak = 20*log10(max(abs(y(:))) + eps); return; end
            n = size(p, 1);
            % running Leq (one-pole on the block mean square)
            a = exp(-n/(obj.fs*obj.tau));
            ms = mean(p(:).^2);
            if obj.pw == 0, obj.pw = ms; else, obj.pw = a*obj.pw + (1 - a)*ms; end
            obj.Lsim = 10*log10(obj.pw/(20e-6)^2 + eps);
            % attenuation: chosen offset, safety cap
            att_ = max(obj.offset, 0); why = {};
            if att_ > 0, why{end+1} = sprintf('%.0f dB quieter (chosen)', att_); end
            if obj.Lsim - att_ > obj.maxL
                att_ = obj.Lsim - obj.maxL; why{end+1} = sprintf('capped at %.0f dB SPL (hearing safety)', obj.maxL);
            end
            obj.att = att_; obj.reason = strjoin(why, '; ');
            g = obj.kPa()*10^(-att_/20);
            x = p*g;
            % loudness compensation, re-designed when the levels moved
            compOn = obj.compensate && att_ > 0.5 && isfinite(obj.Lsim);
            eqOn = ~isempty(obj.hpEq);
            if (compOn || eqOn) && isempty(obj.hNew)
                stale = compOn ~= obj.compDesign || obj.eqVer ~= obj.eqDesign || ...
                    (compOn && ~(abs(obj.Lsim - obj.Ldesign) < 1 && abs(att_ - obj.attDesign) < 1));
                if stale
                    obj.hNew = obj.design(compOn, obj.Lsim, att_); obj.ziNew = obj.zi;
                    if numel(obj.hNew) ~= numel(obj.h), obj.ziNew = []; end
                    obj.Ldesign = obj.Lsim; obj.attDesign = att_; obj.compDesign = compOn; obj.eqDesign = obj.eqVer;
                end
            elseif ~(compOn || eqOn) && ~isequal(obj.h, 1) && isempty(obj.hNew)
                obj.hNew = 1; obj.ziNew = []; obj.Ldesign = NaN; obj.compDesign = false; obj.eqDesign = -1;
            end
            [y, obj.zi] = obj.fir(obj.h, x, obj.zi);
            if ~isempty(obj.hNew)                       % cross-fade to the new FIR over this block
                [y2, obj.ziNew] = obj.fir(obj.hNew, x, obj.ziNew);
                w = (0.5 - 0.5*cos(pi*(1:n)'/n));
                y = (1 - w).*y + w.*y2;
                obj.h = obj.hNew; obj.zi = obj.ziNew; obj.hNew = []; obj.ziNew = [];
            end
            [y, gl] = obj.limit(y);
            if min(gl) < 0.999, obj.reason = strtrim([obj.reason, iff(isempty(obj.reason), '', '; '), ...
                    sprintf('output limited (%.1f dB)', -20*log10(min(gl)))]); end
            obj.Lplay = obj.Lsim - att_ + 20*log10(mean(gl));
            obj.peak = 20*log10(max(abs(y(:))) + eps);
        end

        function setHeadphoneEq(obj, E, Ed)
            % headphone_db.eq result, or [] (none); Ed: its diffuse-field
            % version for direct playback ([] = E); re-designed on the next block
            if nargin < 3, Ed = []; end
            obj.hpEq = E; obj.hpEqDirect = Ed; obj.eqVer = obj.eqVer + 1;
        end

        function y = eqDirect(obj, x, fs, diffuse)
            % one block (n x channels, digital): the headphone EQ only - diffuse
            % (true): the diffuse-field version (direct playback of music),
            % else the eardrum one (uncalibrated binaural)
            if nargin > 2 && fs ~= obj.fs, obj.fs = fs; obj.hD = 1; obj.ziD = []; obj.dKey = ''; end
            E = obj.hpEq; tag = 'ear';
            if diffuse && ~isempty(obj.hpEqDirect), E = obj.hpEqDirect; tag = 'dif'; end
            key = sprintf('%d|%s|%g', obj.eqVer, tag, obj.fs); if isempty(E), key = 'none'; end
            if ~strcmp(key, obj.dKey) && isempty(obj.hDNew)
                if isempty(E)
                    obj.hDNew = 1;
                else
                    C = obj.eqCurve(E, obj.fs);
                    C = C - max(0, max(C));            % no boost above 0 dB: music near full scale is not limited
                    obj.hDNew = obj.mphase(C, 1024, obj.fs);
                end
                obj.ziDNew = []; obj.dKey = key;
            end
            if ~isempty(obj.ziD) && size(obj.ziD, 2) ~= size(x, 2), obj.ziD = []; end   % the channel count changed
            [y, obj.ziD] = obj.fir(obj.hD, x, obj.ziD);
            if ~isempty(obj.hDNew)
                [y2, obj.ziDNew] = obj.fir(obj.hDNew, x, obj.ziDNew);
                n = size(x, 1); w = (0.5 - 0.5*cos(pi*(1:n)'/n));
                y = (1 - w).*y + w.*y2;
                obj.hD = obj.hDNew; obj.ziD = obj.ziDNew; obj.hDNew = []; obj.ziDNew = [];
            end
            y = obj.limit(y);
            obj.peak = 20*log10(max(abs(y(:))) + eps);
        end

        function C = totalCurve(obj, f, compOn, Lsim, att)
            % the output filter [dB] at f: loudness compensation + headphone EQ
            C = zeros(size(f));
            if compOn, C = obj.curve(f, Lsim, att); end
            if ~isempty(obj.hpEq)
                E = obj.hpEq; C = C + interp1(log(E.f(:)), E.dB(:), log(min(max(f, E.f(1)), E.f(end))), 'linear');
            end
        end

        function [x, fsT] = testTone(obj, dur)
            % the 1 kHz reference tone at refDBFS (both ears), dur [s]
            if nargin < 2, dur = 3; end
            fsT = obj.fs; t = (0:round(dur*fsT) - 1)'/fsT;
            x = sqrt(2)*10^(obj.refDBFS/20)*sin(2*pi*1000*t).*min(1, min(t, dur - t)/0.02);
            x = [x x];
        end

        function C = curve(obj, f, Lsim, att)
            % the compensation [dB] at f for the simulated level Lsim and attenuation att
            N = min(max(Lsim, 20), 90); N2 = max(N - att, 20);
            C = nl_iso226(N2, f) - nl_iso226(N, f) + att;
            C = C - nl_iso226(N2, 1000) + nl_iso226(N, 1000) - att;     % 0 dB at 1 kHz
        end

        function s = toStruct(obj)
            fn = {'mode', 'sens', 'sensUnit', 'Z', 'Vfs', 'Lref', 'refDBFS', 'maxL', 'offset', 'compensate', 'headroom', ...
                'tau', 'nTaps', 'enabled', 'hpEq', 'hpEqDirect', 'hrtfName'};
            s = struct(); for k = 1:numel(fn), s.(fn{k}) = obj.(fn{k}); end
        end

        function fromStruct(obj, s)
            for fn = fieldnames(s)', if isprop(obj, fn{1}), obj.(fn{1}) = s.(fn{1}); end, end
            obj.eqVer = obj.eqVer + 1;
        end

        function save(obj, file)
            if nargin < 2, file = playback_calibration.defaultFile(); end
            cal = obj.toStruct(); %#ok<NASGU>
            save(file, 'cal');
        end
    end

    methods (Access = private)
        function [y, gl] = limit(obj, y)
            % sample-peak limiter (fast attack, 50 ms release): never clip
            n = size(y, 1); hr = 10^(-obj.headroom/20);
            pk = max(abs(y), [], 2);
            rel = exp(-1/(0.05*obj.fs)); gl = zeros(n, 1); gc = obj.limGain;
            for i = 1:n
                need = min(1, hr/max(pk(i), eps));
                if need < gc, gc = need; else, gc = rel*gc + (1 - rel)*min(1, need); end
                gl(i) = gc;
            end
            obj.limGain = gc;
            y = y.*gl;
        end

        function C = eqCurve(~, E, fs)
            % a headphone EQ [dB] on the FIR design grid (1024 taps x 4)
            N = 4096; f = max((0:N/2)'*fs/N, 20);
            C = interp1(log(E.f(:)), E.dB(:), log(min(max(f, E.f(1)), E.f(end))), 'linear');
        end

        function hm = mphase(~, C, nT, ~)
            % minimum-phase FIR (nT taps) of the magnitude C [dB] on the
            % grid 0..fs/2 (numel(C) = N/2 + 1)
            N = 2*(numel(C) - 1); mag = 10.^(C/20);
            c = real(ifft(log([mag; flipud(mag(2:end-1))])));
            fold = [c(1); 2*c(2:N/2); c(N/2+1); zeros(N/2 - 1, 1)];
            hm = real(ifft(exp(fft(fold))));
            hm = hm(1:nT).*(0.5 + 0.5*cos(pi*(0:nT-1)'/nT));
        end

        function hm = design(obj, compOn, Lsim, att)
            % minimum-phase FIR of the ISO 226 compensation and the headphone
            % EQ (longer with the EQ: its low-frequency detail)
            nT = obj.nTaps; if ~isempty(obj.hpEq), nT = max(nT, 1024); end
            N = nT*4; fk = (0:N/2)'*obj.fs/N;
            C = obj.totalCurve(max(fk, 20), compOn, Lsim, att);
            if compOn && isempty(obj.hpEq), C(fk > 12500) = C(find(fk <= 12500, 1, 'last')); end
            mag = 10.^(C/20);
            c = real(ifft(log([mag; flipud(mag(2:end-1))])));
            fold = [c(1); 2*c(2:N/2); c(N/2+1); zeros(N/2 - 1, 1)];
            hm = real(ifft(exp(fft(fold))));
            hm = hm(1:nT).*(0.5 + 0.5*cos(pi*(0:nT-1)'/nT));
        end
    end

    methods (Static)
        function [y, zf] = fir(h, x, zi)
            if isequal(h, 1), y = x; zf = zi; return; end
            if isempty(zi), zi = zeros(numel(h) - 1, size(x, 2)); end
            [y, zf] = filter(h, 1, x, zi);
        end

        function f = defaultFile()
            f = fullfile(prefdir, 'binauralizer_playback_calibration.mat');
        end

        function obj = load(file)
            % the saved calibration, or the defaults
            if nargin < 1, file = playback_calibration.defaultFile(); end
            obj = playback_calibration();
            if isfile(file)
                L = load(file, 'cal'); if isfield(L, 'cal'), obj.fromStruct(L.cal); end
            end
        end
    end
end

function out = iff(c, a, b)
if c, out = a; else, out = b; end
end
