classdef HeadphoneLevelPanel < handle
    % HeadphoneLevelPanel  "Headphones & level" (Output setup): the ear (HRTF
    % set), the headphone and its equalisation, and the calibrated playback
    % level of the binaural output.
    %   Ear: a local SOFA set (Data/HRTFs) or one downloaded from the
    %   sofacoustics.org catalogue (hrtf_db) - dummy heads and individual
    %   listeners of ~60 databases.
    %   Headphone: search the AutoEq measurements (headphone_db: ~6800
    %   frequency responses, over-ear / in-ear / earbuds incl. Bluetooth
    %   models); EQ 'eardrum' (flat at the eardrum - matches HRTFs measured
    %   with an ear simulator, such as the default B&K HATS set) or
    %   'diffuse' (to the rig's diffuse-field target - for diffuse-field
    %   equalised HRTF sets); 0 dB at 1 kHz, so the level stays calibrated.
    %   Recognised output devices (autoDetect: a Bluetooth / USB headphone
    %   whose name is an AutoEq model, headphone_db.match) select their
    %   measurement automatically; a popup tells it and lets the EQ off
    %   ("No EQ"); "Detect from output" repeats it. "Auto": the active audio
    %   endpoints (Windows) are watched every 3 s - a headphone connected
    %   while running is selected as the output (SwitchTo, OutputTab) and
    %   recognised; when it disappears the output goes back to Default and
    %   the headphone EQ off.
    %   Level: with "Calibrated playback" the ear signals are pascal
    %   (AppContext.levelScales: the SSD amplifier 0 dBFS = V rms, designed
    %   loudspeakers by their model, plain ones by a sensitivity) and are
    %   played at the simulated SPL (playback_calibration): headphone data
    %   (sensitivity, impedance, the output's full-scale voltage) or a
    %   reference tone measured with an SPL meter (the only way for
    %   Bluetooth headphones: their level is set digitally). Above the
    %   hearing-safety cap, or with a chosen attenuation, quieter with the
    %   ISO 226 loudness compensation. Live meter: simulated / played SPL,
    %   attenuation, output peak.

    properties
        C
        Root
        PopLocal
        PopDb
        PopDbFile
        EdtSearch
        PopForm
        LstHp
        PopEq
        AxHp
        TxtHp
        ChkCal
        PopMode
        Edt = struct()
        PopUnit
        ChkComp
        TxtMeter
        Tmr
        Hits = []               % search result rows (headphone_db table)
        Meas = []               % the chosen headphone's measurement
        LastRow = ''            % the headphone (AutoEq path) last recognised: one popup per headphone
        ChkAuto
        TmrAuto
        Seen = []               % the active audio endpoints at the last watch (cell)
        SwitchTo = []           % fn(deviceName) -> true when the output was switched to it (OutputTab)
        SetBinaural = []        % fn(on) -> the binauralization state after the request (OutputTab: channel check)
        ChkBin
        DFcache = struct('key', '', 'D', [])   % the HRTF set's diffuse-field average (hrtf_df_average)
    end

    methods
        function obj = HeadphoneLevelPanel(parent, ctx)
            obj.C = ctx;
            u = {'Units', 'normalized'};
            bg = get(ancestor(parent, 'figure'), 'Color');
            obj.Root = uipanel('Parent', parent, u{:}, 'Position', [0 0 1 1], 'BorderType', 'none', ...
                'DeleteFcn', @(~,~) obj.stopTimers());
            lab = @(P, s, x, y, w) uicontrol(P, 'Style', 'text', 'String', s, u{:}, 'Position', [x y w 0.06], ...
                'HorizontalAlignment', 'left', 'BackgroundColor', bg);
            % ---- ear ----
            E = uipanel('Parent', obj.Root, u{:}, 'Position', [0.01 0.80 0.98 0.195], 'Title', 'Binauralization (ear: HRTF set)');
            obj.ChkBin = uicontrol(E, 'Style', 'checkbox', 'String', 'Binauralization', 'Value', obj.C.Setup.Binauralization, ...
                u{:}, 'Position', [0.02 0.62 0.24 0.30], 'BackgroundColor', bg, 'FontWeight', 'bold', 'Callback', @(~,~) obj.onBin(), ...
                'TooltipString', 'On: the scene through the HRTF set (headphones); off: the loudspeaker signals directly (multichannel / stereo)');
            L = hrtf_db.local();
            obj.PopLocal = uicontrol(E, 'Style', 'popupmenu', 'String', [{'(local sets)'}, {L.name}], u{:}, ...
                'Position', [0.27 0.62 0.40 0.28], 'UserData', {L.file});
            uicontrol(E, 'Style', 'pushbutton', 'String', '...', u{:}, 'Position', [0.68 0.62 0.07 0.30], ...
                'TooltipString', 'A SOFA file from anywhere', 'Callback', @(~,~) obj.onBrowse());
            uicontrol(E, 'Style', 'pushbutton', 'String', 'Use', u{:}, 'Position', [0.76 0.62 0.22 0.30], ...
                'Callback', @(~,~) obj.onUseLocal());
            obj.PopDb = uicontrol(E, 'Style', 'popupmenu', 'String', {'(online: sofacoustics.org - click to list)'}, u{:}, ...
                'Position', [0.02 0.18 0.36 0.28], 'Callback', @(~,~) obj.onDb());
            obj.PopDbFile = uicontrol(E, 'Style', 'popupmenu', 'String', {'-'}, u{:}, 'Position', [0.40 0.18 0.34 0.28]);
            uicontrol(E, 'Style', 'pushbutton', 'String', 'Download + use', u{:}, 'Position', [0.76 0.18 0.22 0.30], ...
                'Callback', @(~,~) obj.onUseOnline());
            set(obj.PopDb, 'ButtonDownFcn', @(~,~) obj.listDbs());
            % ---- headphone ----
            H = uipanel('Parent', obj.Root, u{:}, 'Position', [0.01 0.47 0.98 0.325], 'Title', 'Headphone (AutoEq measurements)');
            obj.EdtSearch = uicontrol(H, 'Style', 'edit', 'String', '', u{:}, 'Position', [0.02 0.83 0.40 0.13], ...
                'TooltipString', 'Words of the model name, e.g. "airpods pro" or "hd 600"; Enter', 'Callback', @(~,~) obj.onSearch());
            obj.PopForm = uicontrol(H, 'Style', 'popupmenu', 'String', {'all', 'over-ear', 'in-ear', 'earbud'}, u{:}, ...
                'Position', [0.43 0.83 0.17 0.13], 'Callback', @(~,~) obj.onSearch());
            uicontrol(H, 'Style', 'pushbutton', 'String', 'Search', u{:}, 'Position', [0.61 0.83 0.14 0.14], ...
                'Callback', @(~,~) obj.onSearch());
            obj.ChkAuto = uicontrol(H, 'Style', 'checkbox', 'String', 'Auto', 'Value', 1, u{:}, 'Position', [0.02 0.01 0.08 0.14], ...
                'BackgroundColor', bg, 'Callback', @(~,~) obj.onAuto(), 'TooltipString', ...
                'Watch the audio devices: a headphone connected while running becomes the output and is recognised');
            uicontrol(H, 'Style', 'pushbutton', 'String', 'Detect from output', u{:}, 'Position', [0.10 0.01 0.17 0.14], ...
                'TooltipString', 'Recognise the headphone of the current output device (its name) and select its measurement', ...
                'Callback', @(~,~) obj.autoDetect(obj.C.Engine.getCurrentDeviceName(), true));
            obj.PopEq = uicontrol(H, 'Style', 'popupmenu', 'String', {'EQ: off', 'EQ: auto', 'EQ: eardrum (flat)', 'EQ: diffuse field'}, ...
                'Value', 2, u{:}, 'Position', [0.76 0.83 0.22 0.13], 'Callback', @(~,~) obj.applyEq(), 'TooltipString', ...
                ['auto: binaural playback - from the HRTF set (its diffuse-field average: right for eardrum, blocked-canal ' ...
                'and diffuse-field equalised sets alike); binauralization off - diffuse field (music)']);
            obj.LstHp = uicontrol(H, 'Style', 'listbox', 'String', {}, u{:}, 'Position', [0.02 0.18 0.44 0.62], ...
                'Callback', @(~,~) obj.onPickHp());
            obj.AxHp = axes('Parent', H, u{:}, 'Position', [0.55 0.24 0.43 0.52]);
            obj.TxtHp = uicontrol(H, 'Style', 'text', 'String', 'No headphone EQ.', u{:}, 'Position', [0.28 0.01 0.70 0.14], ...
                'HorizontalAlignment', 'left', 'BackgroundColor', bg, 'ForegroundColor', [0.1 0.3 0.6]);
            % ---- level ----
            V = uipanel('Parent', obj.Root, u{:}, 'Position', [0.01 0.13 0.98 0.335], 'Title', 'Playback level');
            P = obj.C.Playback;
            obj.ChkCal = uicontrol(V, 'Style', 'checkbox', 'String', 'Calibrated playback (the simulated SPL)', 'Value', P.enabled, ...
                u{:}, 'Position', [0.02 0.87 0.60 0.11], 'BackgroundColor', bg, 'Callback', @(~,~) obj.onCal());
            obj.PopMode = uicontrol(V, 'Style', 'popupmenu', 'String', {'from headphone data', 'reference tone + SPL meter'}, ...
                'Value', 1 + strcmp(P.mode, 'reference'), u{:}, 'Position', [0.63 0.87 0.35 0.11], 'Callback', @(~,~) obj.onCal());
            row = @(y) y; ew = 0.14;
            mk = @(nm, s, x, y, w, tip) uicontrol(V, 'Style', 'edit', 'String', s, u{:}, 'Position', [x y w 0.10], ...
                'TooltipString', tip, 'Callback', @(~,~) obj.onCal());
            lab(V, 'Sensitivity', 0.02, row(0.70), 0.16);
            obj.Edt.sens = mk('sens', num2str(P.sens), 0.18, 0.72, ew, 'dB SPL per mW or per V (datasheet)');
            obj.PopUnit = uicontrol(V, 'Style', 'popupmenu', 'String', {'dB/mW', 'dB/V'}, 'Value', 1 + strcmp(P.sensUnit, 'V'), ...
                u{:}, 'Position', [0.33 0.72 0.13 0.10], 'Callback', @(~,~) obj.onCal());
            lab(V, 'Z [Ohm]', 0.49, 0.70, 0.12); obj.Edt.Z = mk('Z', num2str(P.Z), 0.61, 0.72, 0.1, 'Impedance');
            lab(V, 'Out FS [V]', 0.73, 0.70, 0.13); obj.Edt.Vfs = mk('Vfs', num2str(P.Vfs), 0.86, 0.72, 0.12, ...
                'V rms of a full-scale sine at the headphone output, with the computer volume as used');
            lab(V, 'Ref. tone SPL', 0.02, 0.55, 0.16); obj.Edt.Lref = mk('Lref', num2str(P.Lref), 0.18, 0.57, ew, ...
                'dB SPL measured at the ear (coupler / SPL meter) for the 1 kHz test tone');
            lab(V, 'at dBFS', 0.33, 0.55, 0.10); obj.Edt.refDBFS = mk('refDBFS', num2str(P.refDBFS), 0.43, 0.57, 0.08, 'Test tone level');
            uicontrol(V, 'Style', 'pushbutton', 'String', 'Play 1 kHz tone', u{:}, 'Position', [0.53 0.565 0.20 0.11], ...
                'Callback', @(~,~) obj.onTone());
            lab(V, 'Max SPL', 0.75, 0.55, 0.11); obj.Edt.maxL = mk('maxL', num2str(P.maxL), 0.86, 0.57, 0.12, 'Hearing-safety cap, dB SPL Leq');
            lab(V, 'Quieter by [dB]', 0.02, 0.40, 0.17); obj.Edt.offset = mk('offset', num2str(P.offset), 0.18, 0.42, ew, ...
                'Deliberate attenuation (the timbre stays: ISO 226)');
            obj.ChkComp = uicontrol(V, 'Style', 'checkbox', 'String', 'ISO 226 compensation', 'Value', P.compensate, u{:}, ...
                'Position', [0.34 0.42 0.30 0.10], 'BackgroundColor', bg, 'Callback', @(~,~) obj.onCal());
            lab(V, 'SSD amp 0 dBFS [V]', 0.02, 0.25, 0.22); obj.Edt.amp = mk('amp', num2str(obj.C.AmpVfs), 0.24, 0.27, 0.10, ...
                'V rms of a full-scale sine at the loudspeaker amplifier output (the SSD signals drive the loudspeakers through it)');
            lab(V, 'Plain speakers [dB @1V,1m]', 0.37, 0.25, 0.30); obj.Edt.plain = mk('plain', num2str(obj.C.PlainSens), 0.67, 0.27, 0.10, ...
                'Sensitivity of the SSD loudspeakers without a designed model');
            uicontrol(V, 'Style', 'pushbutton', 'String', 'Save', u{:}, 'Position', [0.80 0.265 0.18 0.11], ...
                'Callback', @(~,~) obj.onSave(), 'TooltipString', 'Keep these settings (prefdir) for the next start');
            obj.TxtMeter = uicontrol(V, 'Style', 'text', 'String', '', u{:}, 'Position', [0.02 0.01 0.96 0.22], ...
                'HorizontalAlignment', 'left', 'BackgroundColor', bg, 'FontName', 'Consolas', 'FontSize', 9);
            uicontrol(obj.Root, 'Style', 'text', u{:}, 'Position', [0.01 0.005 0.98 0.12], 'BackgroundColor', bg, ...
                'HorizontalAlignment', 'left', 'ForegroundColor', [0.35 0.35 0.35], 'String', ...
                ['Frequency responses: AutoEq (github.com/jaakkopasanen/AutoEq) - oratory1990, Rtings, crinacle, Innerfidelity and ' ...
                'others. HRTFs: sofacoustics.org. Sensitivity / impedance are not in AutoEq: datasheet, Reference Audio Analyzer, ' ...
                'or the reference tone. Bluetooth: use the reference tone (digital volume).']);
            if ~isempty(P.hpEq), obj.showEq(P.hpEq); end
            obj.startMeter();
            obj.onAuto();
        end

        % ---- automatic headphone detection ----
        function onAuto(obj)
            if get(obj.ChkAuto, 'Value') == 1
                if isempty(obj.TmrAuto) || ~isvalid(obj.TmrAuto)
                    obj.Seen = [];
                    obj.TmrAuto = timer('ExecutionMode', 'fixedSpacing', 'Period', 3, 'BusyMode', 'drop', ...
                        'TimerFcn', @(~,~) obj.watch(), 'ObjectVisibility', 'off', 'Name', 'HeadphoneWatch');
                    start(obj.TmrAuto);
                end
            else
                try stop(obj.TmrAuto); delete(obj.TmrAuto); catch, end
                obj.TmrAuto = [];
            end
        end

        function watch(obj)
            % new / gone audio endpoints since the last look
            if ~isvalid(obj.TxtHp), obj.stopTimers(); return; end
            names = audioapp.ui.controls.HeadphoneLevelPanel.endpoints();
            if isempty(names) && isempty(obj.Seen), return; end
            if isempty(obj.Seen), obj.Seen = names; return; end            % the first look: the state at start
            new = setdiff(names, obj.Seen, 'stable'); gone = setdiff(obj.Seen, names, 'stable'); obj.Seen = names;
            if isempty(new) && isempty(gone), return; end
            try
                T = headphone_db.index();
            catch
                return;
            end
            cur = obj.C.Engine.getCurrentDeviceName();
            if ~isempty(gone) && (any(strcmp(gone, cur)) || (~isempty(obj.LastRow) && ...
                    ~any(cellfun(@(n) isMatch(T, n, obj.LastRow), names))))
                % the recognised headphone is gone: default output, no EQ
                if ~isempty(obj.SwitchTo) && any(strcmp(gone, cur)), obj.SwitchTo('Default'); end
                obj.LastRow = ''; set(obj.PopEq, 'Value', 1); obj.applyEq();
                set(obj.TxtHp, 'String', sprintf('Headphone disconnected (%s): default output, no headphone EQ.', strjoin(gone, ', ')));
            end
            new = new(cellfun(@isempty, regexpi(new, 'hands-?free', 'once')));   % the telephone (HFP) endpoint: never
            % a headphone has several endpoints (stereo output, headset mic):
            % the output one is switched to (-> autoDetect via the device
            % change); only when none of them is an output, the EQ alone
            hit = {};
            for k = 1:numel(new)
                row = headphone_db.match(T, new{k});
                if ~isempty(row) && ~strcmp(obj.LastRow, char(row.path)), hit{end+1} = new{k}; end %#ok<AGROW>
            end
            if isempty(hit), return; end
            for k = 1:numel(hit)
                if ~isempty(obj.SwitchTo) && obj.SwitchTo(hit{k}), return; end
            end
            obj.autoDetect(hit{1});
        end

        function stopTimers(obj)
            obj.stopMeter();
            try stop(obj.TmrAuto); delete(obj.TmrAuto); catch, end
        end

        % ---- ear ----
        function onUseLocal(obj)
            v = get(obj.PopLocal, 'Value'); F = get(obj.PopLocal, 'UserData');
            if v < 2, return; end
            obj.useHrtf(F{v - 1});
        end

        function listDbs(obj)
            if numel(get(obj.PopDb, 'String')) > 1, return; end
            set(obj.TxtHp, 'String', 'Listing sofacoustics.org...'); drawnow;
            try
                D = hrtf_db.databases(); set(obj.PopDb, 'String', [{'(online database)'}, D(:)'], 'Value', 1);
                set(obj.TxtHp, 'String', sprintf('%d HRTF databases online.', numel(D)));
            catch ME
                set(obj.TxtHp, 'String', ['sofacoustics.org: ' ME.message]);
            end
        end

        function onDb(obj)
            obj.listDbs(); s = get(obj.PopDb, 'String'); v = get(obj.PopDb, 'Value');
            if v < 2, return; end
            try
                F = hrtf_db.files(s{v});
                if isempty(F), F = {'(no SOFA files at the top level)'}; end
                set(obj.PopDbFile, 'String', F, 'Value', 1);
            catch ME
                set(obj.PopDbFile, 'String', {['error: ' ME.message]}, 'Value', 1);
            end
        end

        function onUseOnline(obj)
            s = get(obj.PopDb, 'String'); v = get(obj.PopDb, 'Value'); F = get(obj.PopDbFile, 'String');
            if v < 2 || ~endsWith(lower(F{get(obj.PopDbFile, 'Value')}), '.sofa'), return; end
            set(obj.TxtHp, 'String', 'Downloading the HRTF set...'); drawnow;
            try
                file = hrtf_db.download(s{v}, F{get(obj.PopDbFile, 'Value')});
            catch ME
                set(obj.TxtHp, 'String', ['Download failed: ' ME.message]); return;
            end
            L = hrtf_db.local(); set(obj.PopLocal, 'String', [{'(local sets)'}, {L.name}], 'UserData', {L.file});
            obj.useHrtf(file);
        end

        function useHrtf(obj, file)
            set(obj.TxtHp, 'String', 'Loading the HRTF set and rebuilding the scene...'); drawnow;
            try
                obj.C.loadHRTF(file);
                [~, n, e] = fileparts(file); obj.C.Playback.hrtfName = [n e];
                set(obj.TxtHp, 'String', sprintf('Ear: %s%s (fs %g Hz).', n, e, obj.C.Setup.SampleRate));
                if ~isempty(obj.Meas), obj.applyEq(); end            % the binaural EQ follows the set
            catch ME
                set(obj.TxtHp, 'String', ['HRTF: ' ME.message]);
            end
        end

        function onBrowse(obj)
            [fn, p] = uigetfile({'*.sofa', 'HRTF set (*.sofa)'}, 'HRTF set', fullfile(pwd, 'Data', 'HRTFs'));
            if isequal(fn, 0), return; end
            obj.useHrtf(fullfile(p, fn));
        end

        function onBin(obj)
            % binauralization on / off (the channel check: OutputTab); the
            % headphone EQ shown switches with it
            on = get(obj.ChkBin, 'Value') == 1;
            if ~isempty(obj.SetBinaural)
                on = obj.SetBinaural(on);
            else
                obj.C.Setup.Binauralization = on; obj.C.rebuildSoundScene();
            end
            set(obj.ChkBin, 'Value', on);
            if ~isempty(obj.Meas), obj.applyEq(); end
        end

        function D = hrtfDF(obj)
            % the current HRTF set's diffuse-field average (cached), [] when
            % there is none (BRIR set, no set)
            D = [];
            H = []; try, H = obj.C.Setup.HRTF; catch, end
            if isempty(H) || ~isstruct(H) || ~isfield(H, 'Data'), return; end
            key = sprintf('%s|%g|%s', mat2str(size(H.Data.IR)), H.Data.SamplingRate, obj.C.Playback.hrtfName);
            if strcmp(obj.DFcache.key, key), D = obj.DFcache.D; return; end
            try, D = hrtf_df_average(H); catch, D = []; end
            obj.DFcache = struct('key', key, 'D', D);
        end

        % ---- headphone ----
        function onSearch(obj)
            q = strtrim(get(obj.EdtSearch, 'String'));
            forms = {'', 'over-ear', 'in-ear', 'earbud'}; fm = forms{get(obj.PopForm, 'Value')};
            set(obj.TxtHp, 'String', 'Searching...'); drawnow;
            try
                T = headphone_db.index();
            catch ME
                set(obj.TxtHp, 'String', ['AutoEq index: ' ME.message]); return;
            end
            R = headphone_db.search(T, q, fm); obj.Hits = R;
            rg = R.rig; rg(rg == "") = "-";
            set(obj.LstHp, 'String', cellstr(R.name + "   [" + R.source + ", " + R.form + ", " + rg + "]"), 'Value', min(1, height(R)));
            set(obj.TxtHp, 'String', sprintf('%d of %d measurements.', height(R), height(T)));
        end

        function onPickHp(obj)
            i = get(obj.LstHp, 'Value'); if isempty(obj.Hits) || isempty(i) || i > height(obj.Hits), return; end
            set(obj.TxtHp, 'String', 'Downloading the measurement...'); drawnow;
            try
                obj.Meas = headphone_db.fetch(obj.Hits(i, :));
            catch ME
                set(obj.TxtHp, 'String', ['Measurement: ' ME.message]); return;
            end
            obj.applyEq();
        end

        function found = autoDetect(obj, dev, force)
            % the headphone of the output device dev, if its name is a model
            % of the AutoEq list: selected (search box, list, measurement),
            % its EQ applied; a popup tells it and offers "No EQ". force:
            % also for the device recognised last time (the button)
            if nargin < 3, force = false; end
            found = false;
            try
                T = headphone_db.index();
            catch
                return;                                         % no index (offline, first start)
            end
            row = headphone_db.match(T, dev);
            if isempty(row)
                if force, set(obj.TxtHp, 'String', sprintf('"%s": not a headphone of the AutoEq list (search it above).', dev)); end
                return;
            end
            found = true;
            if ~force && strcmp(obj.LastRow, char(row.path)), return; end
            obj.LastRow = char(row.path);
            set(obj.EdtSearch, 'String', char(row.name)); set(obj.PopForm, 'Value', 1);
            obj.onSearch();
            i = find(obj.Hits.path == row.path, 1); if isempty(i), i = 1; end
            set(obj.LstHp, 'Value', i);
            if get(obj.PopEq, 'Value') == 1, set(obj.PopEq, 'Value', 2); end   % recognised: EQ on (eardrum) by default
            obj.onPickHp();
            ok = ~isempty(obj.Meas) && strcmp(obj.Meas.name, char(row.name));
            msg = sprintf('Recognised headphone on the output "%s":\n\n    %s   (%s, measured by %s)\n\n', ...
                dev, char(row.name), char(row.form), char(row.source));
            if ok
                msg = [msg 'Its equalisation is on: ' strrep(get(obj.TxtHp, 'String'), [char(row.name) ' '], '')];
            else
                msg = [msg 'Its measurement could not be loaded (offline?): no EQ yet.'];
            end
            if ~(usejava('desktop') && strcmp(get(0, 'DefaultFigureVisible'), 'on')), return; end   % scripts: no dialog
            if ok
                a = questdlg([msg newline newline 'You can change this any time under Output setup > Headphones & level.'], ...
                    'Headphone recognised', 'Use the EQ', 'No EQ', 'Use the EQ');
                if strcmp(a, 'No EQ'), set(obj.PopEq, 'Value', 1); obj.applyEq(); end
            else
                msgbox(msg, 'Headphone recognised', 'modal');
            end
        end

        function applyEq(obj)
            P = obj.C.Playback; m = get(obj.PopEq, 'Value');
            if m == 1 || isempty(obj.Meas)
                if m == 1, P.setHeadphoneEq([]); cla(obj.AxHp); set(obj.TxtHp, 'String', 'No headphone EQ.'); end
                return;
            end
            % binaural playback: auto = from the HRTF set (its diffuse-field
            % average), else the chosen target; direct playback (binauralization
            % off): diffuse field - music keeps the ear's own resonances
            D = []; how = '';
            switch m
                case 2
                    D = obj.hrtfDF();
                    if isempty(D), o = struct('mode', 'eardrum'); how = 'eardrum (no HRTF average)';
                    else, o = struct('mode', 'hrtf', 'hrtfDF', D); how = sprintf('from the HRTF set: %s', D.text);
                    end
                case 3, o = struct('mode', 'eardrum'); how = 'eardrum (flat)';
                otherwise, o = struct('mode', 'diffuse'); how = 'diffuse field';
            end
            try
                E = headphone_db.eq(obj.Meas, o);
            catch ME
                set(obj.TxtHp, 'String', ['EQ: ' ME.message]); return;
            end
            Ed = [];
            try, Ed = headphone_db.eq(obj.Meas, struct('mode', 'diffuse')); catch, end
            P.setHeadphoneEq(E, Ed);
            bin = obj.C.Setup.Binauralization;
            if bin || isempty(Ed), obj.showEq(E, obj.Meas); else, obj.showEq(Ed, obj.Meas); end
            set(obj.TxtHp, 'String', sprintf('%s [%s]. Binaural: %s. Direct: diffuse field. Active now: %s.', ...
                E.name, E.source, how, iff(bin, 'binaural', 'direct')));
        end

        function showEq(obj, E, M)
            ax = obj.AxHp; cla(ax); hold(ax, 'on');
            if nargin > 2 && ~isempty(M)
                r = M.raw - interp1(log(M.f), M.raw, log(1000));
                semilogx(ax, M.f, r, 'Color', [0.6 0.6 0.6]);
            end
            semilogx(ax, E.f, E.dB, 'b', 'LineWidth', 1.3);
            hold(ax, 'off'); set(ax, 'XScale', 'log', 'FontSize', 7); grid(ax, 'on'); xlim(ax, [20 20000]); ylim(ax, [-20 20]);
            title(ax, 'measured (grey), EQ (blue) [dB re 1 kHz]', 'FontSize', 7);
            set(obj.TxtHp, 'String', sprintf('%s [%s, %s]: EQ %s.', E.name, E.source, iff(isempty(E.rig), '-', E.rig), E.mode));
        end

        % ---- level ----
        function onCal(obj)
            P = obj.C.Playback;
            num = @(e, d) iff(isnan(str2double(get(e, 'String'))), d, str2double(get(e, 'String')));
            P.enabled = get(obj.ChkCal, 'Value') == 1;
            P.mode = iff(get(obj.PopMode, 'Value') == 1, 'headphone', 'reference');
            P.sens = num(obj.Edt.sens, P.sens); P.sensUnit = iff(get(obj.PopUnit, 'Value') == 1, 'mW', 'V');
            P.Z = num(obj.Edt.Z, P.Z); P.Vfs = num(obj.Edt.Vfs, P.Vfs);
            P.Lref = num(obj.Edt.Lref, P.Lref); P.refDBFS = num(obj.Edt.refDBFS, P.refDBFS);
            P.maxL = num(obj.Edt.maxL, P.maxL); P.offset = max(num(obj.Edt.offset, P.offset), 0);
            P.compensate = get(obj.ChkComp, 'Value') == 1;
            obj.C.AmpVfs = num(obj.Edt.amp, obj.C.AmpVfs); obj.C.PlainSens = num(obj.Edt.plain, obj.C.PlainSens);
            P.reset(obj.C.Setup.SampleRate);
            obj.C.applyLevelScales();
        end

        function onTone(obj)
            obj.onCal();
            obj.C.ensureEngine();
            if ~obj.C.Engine.playCalibrationTone()
                set(obj.TxtMeter, 'String', 'Stop the playback first (the tone uses the same output).'); return;
            end
            set(obj.TxtMeter, 'String', sprintf(['1 kHz at %g dBFS played: measure it at the ear (SPL meter / coupler) and ' ...
                'enter the reading as "Ref. tone SPL", mode "reference tone".'], obj.C.Playback.refDBFS));
        end

        function onSave(obj)
            obj.onCal(); obj.C.savePlayback();
            set(obj.TxtMeter, 'String', ['Saved: ' playback_calibration.defaultFile()]);
        end

        function startMeter(obj)
            obj.Tmr = timer('ExecutionMode', 'fixedSpacing', 'Period', 0.3, 'BusyMode', 'drop', ...
                'TimerFcn', @(~,~) obj.updateMeter(), 'ObjectVisibility', 'off', 'Name', 'HeadphoneLevelMeter');
            start(obj.Tmr);
        end

        function stopMeter(obj)
            % (the meter timer only)
            try stop(obj.Tmr); delete(obj.Tmr); catch, end
        end

        function updateMeter(obj)
            if ~isvalid(obj.TxtMeter), obj.stopMeter(); return; end
            if get(obj.ChkBin, 'Value') ~= obj.C.Setup.Binauralization   % switched elsewhere (renderer choice)
                set(obj.ChkBin, 'Value', obj.C.Setup.Binauralization);
                if ~isempty(obj.Meas), obj.applyEq(); end
            end
            P = obj.C.Playback;
            if ~P.enabled
                s = sprintf('Uncalibrated (the ear signals go to the output 1:1).   0 dBFS would be %.1f dB SPL with these data.', P.SPLfs());
            else
                s = sprintf('Simulated at the ears: %5.1f dB SPL   played: %5.1f dB SPL   attenuation %4.1f dB   peak %5.1f dBFS\n0 dBFS = %.1f dB SPL.  %s', ...
                    P.Lsim, P.Lplay, P.att, P.peak, P.SPLfs(), P.reason);
            end
            set(obj.TxtMeter, 'String', s);
        end
    end

    methods (Static)
        function names = endpoints()
            % the active audio endpoints of Windows (friendly names, as the
            % output devices are listed), {} elsewhere / on failure
            names = {};
            if ~ispc, return; end
            try
                NET.addAssembly('System.Management');
                q = System.Management.ManagementObjectSearcher(['SELECT Name, Status FROM Win32_PnPEntity ' ...
                    'WHERE ClassGuid = ''{c166523c-fe0c-4a94-a586-f1a80cfbbf3e}''']);          % AudioEndpoint
                e = q.Get().GetEnumerator();
                while e.MoveNext()
                    o = e.Current;
                    if strcmpi(char(o.GetPropertyValue('Status')), 'OK')              % connected
                        names{end+1} = char(o.GetPropertyValue('Name')); %#ok<AGROW>
                    end
                end
            catch
                names = {};
            end
        end
    end
end

function tf = isMatch(T, name, path)
% the endpoint name is the headphone of the AutoEq row path
r = headphone_db.match(T, name); tf = ~isempty(r) && strcmp(char(r.path), path);
end

function out = iff(c, a, b)
if c, out = a; else, out = b; end
end
