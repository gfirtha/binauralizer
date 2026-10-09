classdef SeparationPanel < handle
    % SeparationPanel  Top-level "Source separation" tab.
    %
    %   Uses the media loaded in the Binauralizer. Pick a method and its
    %   parameters, choose an excerpt (start, length) and "Analyze excerpt":
    %     - direction-time map of the primary energy with the tracked
    %       source paths,
    %     - object levels over time,
    %     - scene (top view): where the separator put each source; click a
    %       source to hear its separated signal alone (also: the original
    %       excerpt, all objects summed, stop).
    %   "Send to Binauralizer" plays the media through the separator (its
    %   objects become the virtual sources); "Bypass" returns to the
    %   media's own channels.

    properties
        C
        Root
        TxtMedia
        PopMethod
        Edt = struct()
        ChkDemix
        EdtStart
        EdtLen
        SldStart
        SldLen
        TxtStatus
        AxMap
        AxLev
        AxScene
        Lsn = {}
        Ex = []             % last analysed excerpt: X, Y, fs, Names, t0, colours, ...
        Player = []
        LevLines = gobjects(0)
        SceneMarks = gobjects(0)
    end

    properties (Constant)
        % legacy methods (Files/Objects/Separation/legacy) are not offered;
        % titles and descriptions: separator_info
        Methods = {'stft_objects', 'stft', 'cluster', 'wavenumber', 'soundfield', 'avendano', 'ltrt', 'passthrough'}
        MaxLen = 60             % longest excerpt [s] on the length slider
        Params = {'K', 'Objects / front slots', 1; 'Width', 'Front half-width [deg]', 1; 'Base', 'Content stereo base [deg]', 1; ...
                  'AmbAz', 'Ambience azimuth [deg]', 1; 'AmbGain', 'Ambience gain', 1; 'Tau', 'Smoothing [ms]', 1000}
    end

    methods
        function obj = SeparationPanel(parent, C)
            obj.C = C;
            obj.build(parent);
            obj.Lsn{end+1} = addlistener(C, 'SetupChanged', @(~,~) obj.refresh());
            obj.refresh();
        end

        function delete(obj)
            cellfun(@delete, obj.Lsn);
            obj.stopPreview();
        end

        function build(obj, parent)
            u = {'Units', 'normalized'};
            bg = get(ancestor(parent, 'figure'), 'Color');
            obj.Root = uipanel('Parent', parent, u{:}, 'Position', [0 0 1 1], 'BorderType', 'none');
            L = uipanel('Parent', obj.Root, u{:}, 'Position', [0.01 0.02 0.28 0.96], 'Title', 'Separation / upmix');
            txt = @(s, y, h) uicontrol(L, 'Style', 'text', 'String', s, u{:}, 'Position', [0.04 y 0.92 h], ...
                'HorizontalAlignment', 'left', 'BackgroundColor', bg);
            txt('Media (loaded in the Binauralizer):', 0.94, 0.035);
            obj.TxtMedia = txt('-', 0.88, 0.06);
            set(obj.TxtMedia, 'FontAngle', 'italic');
            txt('Method:', 0.84, 0.035);
            names = cellfun(@(m) getfield(separator_info(m), 'Title'), obj.Methods, 'UniformOutput', false);
            obj.PopMethod = uicontrol(L, 'Style', 'popupmenu', 'String', names, u{:}, ...
                'Position', [0.04 0.795 0.70 0.045], 'Callback', @(~,~) obj.onMethod());
            uicontrol(L, 'Style', 'pushbutton', 'String', 'Info', u{:}, 'Position', [0.76 0.797 0.20 0.043], ...
                'TooltipString', 'What the selected separator does', 'Callback', @(~,~) obj.onInfo());
            y = 0.745;
            for k = 1:size(obj.Params, 1)
                uicontrol(L, 'Style', 'text', 'String', obj.Params{k,2}, u{:}, 'Position', [0.04 y 0.60 0.035], ...
                    'HorizontalAlignment', 'left', 'BackgroundColor', bg);
                obj.Edt.(obj.Params{k,1}) = uicontrol(L, 'Style', 'edit', u{:}, 'Position', [0.66 y 0.30 0.04], ...
                    'Callback', @(s,~) obj.onParam(s, k));
                y = y - 0.045;
            end
            obj.ChkDemix = uicontrol(L, 'Style', 'checkbox', 'String', 'Informed demixing (experimental)', u{:}, ...
                'Position', [0.04 y 0.92 0.04], 'BackgroundColor', bg, 'Callback', @(s,~) obj.onDemix(s));
            % excerpt: start / length, edit box + slider each
            y = y - 0.05;
            uicontrol(L, 'Style', 'text', 'String', 'Excerpt start [s]:', u{:}, 'Position', [0.04 y 0.56 0.035], ...
                'HorizontalAlignment', 'left', 'BackgroundColor', bg);
            obj.EdtStart = uicontrol(L, 'Style', 'edit', 'String', '0', u{:}, 'Position', [0.74 y 0.22 0.04], ...
                'TooltipString', 'Start time in the media [s] (empty: current playback time)', ...
                'Callback', @(~,~) obj.onExcerptEdit());
            obj.SldStart = uicontrol(L, 'Style', 'slider', u{:}, 'Position', [0.04 y-0.036 0.92 0.03], ...
                'Min', 0, 'Max', 1, 'Value', 0, 'Callback', @(~,~) obj.onExcerptSlider());
            y = y - 0.075;
            uicontrol(L, 'Style', 'text', 'String', 'Excerpt length [s]:', u{:}, 'Position', [0.04 y 0.56 0.035], ...
                'HorizontalAlignment', 'left', 'BackgroundColor', bg);
            obj.EdtLen = uicontrol(L, 'Style', 'edit', 'String', '5', u{:}, 'Position', [0.74 y 0.22 0.04], ...
                'Callback', @(~,~) obj.onExcerptEdit());
            obj.SldLen = uicontrol(L, 'Style', 'slider', u{:}, 'Position', [0.04 y-0.036 0.92 0.03], ...
                'Min', 0.5, 'Max', obj.MaxLen, 'Value', 5, 'Callback', @(~,~) obj.onExcerptSlider());
            uicontrol(L, 'Style', 'pushbutton', 'String', 'Analyze excerpt', u{:}, ...
                'Position', [0.04 0.235 0.92 0.055], 'Callback', @(~,~) obj.onAnalyze());
            uicontrol(L, 'Style', 'pushbutton', 'String', 'Send to Binauralizer', u{:}, ...
                'Position', [0.04 0.17 0.92 0.055], 'FontWeight', 'bold', 'Callback', @(~,~) obj.onSend());
            uicontrol(L, 'Style', 'pushbutton', 'String', 'Bypass (plain channels)', u{:}, ...
                'Position', [0.04 0.105 0.92 0.055], 'Callback', @(~,~) obj.onBypass());
            obj.TxtStatus = txt('', 0.005, 0.095);
            set(obj.TxtStatus, 'ForegroundColor', [0.1 0.3 0.6]);

            obj.AxMap = axes('Parent', obj.Root, u{:}, 'Position', [0.35 0.60 0.53 0.35]);
            title(obj.AxMap, 'Primary energy: direction over time');
            xlabel(obj.AxMap, 't [s]'); ylabel(obj.AxMap, 'azimuth [deg]  (+ = left)');
            obj.AxLev = axes('Parent', obj.Root, u{:}, 'Position', [0.35 0.08 0.36 0.38]);
            title(obj.AxLev, 'Object levels'); xlabel(obj.AxLev, 't [s]'); ylabel(obj.AxLev, 'dBFS');
            obj.AxScene = axes('Parent', obj.Root, u{:}, 'Position', [0.75 0.08 0.23 0.38]);
            title(obj.AxScene, 'Sources (click to listen)'); axis(obj.AxScene, 'off');
            uicontrol(obj.Root, 'Style', 'pushbutton', 'String', 'Play original', u{:}, ...
                'Position', [0.75 0.495 0.075 0.04], 'Callback', @(~,~) obj.playOriginal());
            uicontrol(obj.Root, 'Style', 'pushbutton', 'String', 'Play all objects', u{:}, ...
                'Position', [0.828 0.495 0.085 0.04], 'Callback', @(~,~) obj.playAll());
            uicontrol(obj.Root, 'Style', 'pushbutton', 'String', 'Stop', u{:}, ...
                'Position', [0.916 0.495 0.06 0.04], 'Callback', @(~,~) obj.stopPreview());
        end

        function refresh(obj)
            if isempty(obj.TxtMedia) || ~isgraphics(obj.TxtMedia), return; end
            E = obj.C.Engine; f = '';
            try, f = E.MediaFilename; catch, end
            [~, nm, ex] = fileparts(f);
            state = '';
            try, if E.isSeparated(), state = '   [playing separated]'; end, catch, end
            set(obj.TxtMedia, 'String', [nm ex state]);
            cfg = obj.C.Setup.Separation;
            v = find(strcmpi(obj.Methods, cfg.Method), 1);
            if isempty(v)                                   % legacy / unknown method in the setup
                v = 1; obj.C.Setup.Separation.Method = obj.Methods{1}; cfg = obj.C.Setup.Separation;
            end
            set(obj.PopMethod, 'Value', v);
            I = separator_info(cfg.Method); used = I.Params;   % only the parameters the method reads
            for k = 1:size(obj.Params, 1)
                f = obj.Params{k,1};
                set(obj.Edt.(f), 'String', num2str(cfg.(f) * obj.Params{k,3}), ...
                    'Enable', char(matlab.lang.OnOffSwitchState(any(strcmp(used, f)))));
            end
            set(obj.ChkDemix, 'Value', isfield(cfg, 'Demix') && cfg.Demix, ...
                'Enable', char(matlab.lang.OnOffSwitchState(strcmpi(cfg.Method, 'stft_objects'))));
            obj.syncExcerpt();
        end
    end

    methods (Access = private)
        % ---------------- settings ----------------
        function onMethod(obj)
            obj.C.Setup.Separation.Method = obj.Methods{get(obj.PopMethod, 'Value')};
            obj.refresh();
        end

        function onParam(obj, src, k)
            v = str2double(get(src, 'String'));
            f = obj.Params{k,1};
            if isnan(v)
                set(src, 'String', num2str(obj.C.Setup.Separation.(f) * obj.Params{k,3}));
                return;
            end
            obj.C.Setup.Separation.(f) = v / obj.Params{k,3};
        end

        function onDemix(obj, src)
            obj.C.Setup.Separation.Demix = logical(get(src, 'Value'));
        end

        function onInfo(obj)
            m = obj.Methods{get(obj.PopMethod, 'Value')};
            try
                I = separator_info(m); ttl = I.Title; info = I.Text;
            catch ME
                ttl = m; info = sprintf('Unable to retrieve separator info.\n\nError:\n%s', ME.message);
            end
            h = figure('Name', 'Separator info', 'NumberTitle', 'off', 'MenuBar', 'none', ...
                'ToolBar', 'none', 'Units', 'normalized', 'Position', [0.22 0.15 0.56 0.70], ...
                'Color', [1 1 1], 'WindowStyle', 'modal');
            uicontrol(h, 'Style', 'text', 'String', ttl, ...
                'Units', 'normalized', 'Position', [0.04 0.93 0.92 0.05], ...
                'HorizontalAlignment', 'left', 'BackgroundColor', [1 1 1], 'FontWeight', 'bold');
            uicontrol(h, 'Style', 'edit', 'Max', 2, 'Min', 0, 'Enable', 'inactive', ...
                'Units', 'normalized', 'Position', [0.04 0.10 0.92 0.82], ...
                'HorizontalAlignment', 'left', 'String', splitlines(info), 'FontName', 'Consolas');
            uicontrol(h, 'Style', 'pushbutton', 'String', 'Close', ...
                'Units', 'normalized', 'Position', [0.80 0.025 0.16 0.055], ...
                'Callback', @(~,~) close(h));
        end

        % ---------------- excerpt (start / length) ----------------
        function onExcerptEdit(obj)
            obj.syncExcerpt();
        end

        function onExcerptSlider(obj)
            set(obj.EdtStart, 'String', sprintf('%.1f', get(obj.SldStart, 'Value')));
            set(obj.EdtLen, 'String', sprintf('%.1f', get(obj.SldLen, 'Value')));
            obj.syncExcerpt();
        end

        function syncExcerpt(obj)
            % clamp start / length to the media and move the sliders there;
            % an empty start box (= current playback time) is left empty
            if isempty(obj.SldStart) || ~isgraphics(obj.SldStart), return; end
            d = 0;
            try, d = obj.C.Engine.getDuration(); catch, end
            if ~(d > 0), d = obj.MaxLen; end
            lenMax = max(0.5, min(obj.MaxLen, d));
            len = str2double(get(obj.EdtLen, 'String'));
            if isnan(len), len = 5; end
            len = min(max(len, 0.5), lenMax);
            set(obj.SldLen, 'Max', lenMax, 'Value', len, ...
                'SliderStep', min(1, [0.5 5]/max(lenMax - 0.5, eps)));
            set(obj.EdtLen, 'String', sprintf('%.1f', len));
            sMax = max(d - len, 0);
            t0 = str2double(get(obj.EdtStart, 'String'));
            if ~isnan(t0)
                t0 = min(max(t0, 0), sMax);
                set(obj.EdtStart, 'String', sprintf('%.1f', t0));
            else
                t0 = min(get(obj.SldStart, 'Value'), sMax);
            end
            set(obj.SldStart, 'Max', max(sMax, eps), 'Value', t0, ...
                'SliderStep', min(1, [1 10]/max(sMax, 1)), 'Enable', char(matlab.lang.OnOffSwitchState(sMax > 0)));
        end

        function status(obj, msg)
            set(obj.TxtStatus, 'String', msg); drawnow;
        end

        % ---------------- binauralizer ----------------
        function onSend(obj)
            cfg = obj.C.Setup.Separation;
            obj.stopPreview();
            obj.status('Opening the media through the separator...');
            try
                obj.C.Engine.openSeparation(cfg);
            catch ME
                obj.status(['Failed: ' ME.message]); return;
            end
            obj.C.rebuildSoundScene();                       % N_in = number of objects
            rd = obj.C.Engine.Objects;
            obj.C.Scene.ObjectRoomRadius = obj.C.Setup.Loudspeaker_setup.R + 1;
            obj.C.Scene.apply_object_positions(rd.positions(), ones(rd.nObj, 1));
            obj.C.Scene.flush_updates();
            obj.status(sprintf('Playing through "%s": %d objects (%s). Press Play in the Binauralizer.', ...
                cfg.Method, rd.nObj, strjoin(rd.Names, ', ')));
            obj.refresh();
        end

        function onBypass(obj)
            try
                obj.C.Engine.bypassSeparation();
                obj.C.rebuildSoundScene();
                obj.status('Back to the media''s own channels.');
            catch ME
                obj.status(['Failed: ' ME.message]);
            end
            obj.refresh();
        end

        % ---------------- analysis ----------------
        function onAnalyze(obj)
            E = obj.C.Engine;
            if isempty(E.MediaFilename), obj.status('Load a media file in the Binauralizer first.'); return; end
            cfg = obj.C.Setup.Separation;
            fs = obj.C.Setup.SampleRate;
            obj.stopPreview();
            obj.status('Analyzing (the first run on a film extracts its audio)...');
            sep = make_separator(cfg);
            if isprop(sep, 'Diagnostics'), sep.Diagnostics = true; end
            if isa(sep, 'stft_separator') && isfield(cfg, 'Demix'), sep.Demix = cfg.Demix; end
            try
                rd = separation_reader(E.MediaFilename, sep, fs);   % converts + inits sep
            catch ME
                obj.status(['Failed: ' ME.message]); return;
            end
            len = max(0.5, str2double(get(obj.EdtLen, 'String')));
            t0 = str2double(get(obj.EdtStart, 'String'));
            if isnan(t0), t0 = E.getCurrentTime(); end
            t0 = max(0, min(t0, rd.duration() - len));
            set(obj.EdtStart, 'String', sprintf('%.1f', t0)); set(obj.EdtLen, 'String', num2str(len));
            obj.syncExcerpt();
            a = round(t0*fs) + 1; b = min(rd.total, a + round(len*fs) - 1);
            X = audioread(rd.wav, [a b]);
            tic;
            if ismethod(sep, 'separate_offline')
                Y = sep.separate_offline({X}); Y = Y{1};            % aligned with X
            else
                Y = sep.process(X);                                  % zero-latency methods
            end
            tc = toc;
            obj.Ex = struct('X', X, 'Y', Y, 'fs', fs, 't0', t0, 'Names', {sep.Names}, ...
                'colors', lines(sep.nOut()));
            obj.drawLevels();
            obj.drawMap(sep);
            obj.drawScene(sep);
            obj.status(sprintf('Analyzed %.1f s from %.1f s in %.1f s. Click a source in the scene to hear it alone.', ...
                len, t0, tc));
        end

        function drawLevels(obj)
            Ex = obj.Ex; fr = round(0.05*Ex.fs); nf = floor(size(Ex.Y,1)/fr);
            lev = 10*log10(squeeze(mean(reshape(Ex.Y(1:nf*fr,:).^2, fr, nf, []), 1)) + 1e-12);
            if isvector(lev), lev = lev(:); end
            tl = Ex.t0 + ((1:nf) - 0.5)*fr/Ex.fs;
            cla(obj.AxLev); hold(obj.AxLev, 'on');
            obj.LevLines = gobjects(1, size(lev,2));
            for k = 1:size(lev,2)
                obj.LevLines(k) = plot(obj.AxLev, tl, lev(:,k), 'Color', Ex.colors(k,:), 'LineWidth', 1, ...
                    'ButtonDownFcn', @(~,~) obj.playObject(k));
            end
            hold(obj.AxLev, 'off'); grid(obj.AxLev, 'on');
            legend(obj.AxLev, Ex.Names, 'Location', 'eastoutside', 'FontSize', 7);
            ylim(obj.AxLev, [max(min(lev(:)), max(lev(:)) - 60), max(lev(:)) + 3]);
            title(obj.AxLev, 'Object levels'); xlabel(obj.AxLev, 't [s]'); ylabel(obj.AxLev, 'dBFS');
        end

        function drawMap(obj, sep)
            Ex = obj.Ex; cla(obj.AxMap);
            if isprop(sep, 'DiagAz') && ~isempty(sep.DiagAz)
                M = 10*log10(sep.DiagAz + 1e-12);
                az = (sep.AzEdges(1:end-1) + sep.AzEdges(2:end))/2;
                tm = Ex.t0 + (0:size(M,2)-1)*sep.N/2/Ex.fs;
                imagesc(obj.AxMap, tm, az, M); axis(obj.AxMap, 'xy');
                clim(obj.AxMap, max(M(:)) + [-50 0]); colormap(obj.AxMap, parula);
                ylim(obj.AxMap, [min(sep.AzEdges), max(sep.AzEdges)]);
                pa = 10*log10(sum(sep.DiagPA(1,:))/max(sum(sep.DiagPA(2,:)), eps));
                title(obj.AxMap, sprintf('Primary energy: direction over time   (primary / ambience %.1f dB)', pa));
                if ~isempty(sep.TrackLog)
                    hold(obj.AxMap, 'on');
                    for k = 1:size(sep.TrackLog,1)
                        plot(obj.AxMap, tm, sep.TrackLog(k,:), 'Color', Ex.colors(k,:), 'LineWidth', 2);
                    end
                    hold(obj.AxMap, 'off');
                end
            else
                text(obj.AxMap, 0.5, 0.5, 'Direction map: STFT methods only', 'Units', 'normalized', ...
                    'HorizontalAlignment', 'center');
                title(obj.AxMap, 'Primary energy: direction over time');
            end
            xlabel(obj.AxMap, 't [s]'); ylabel(obj.AxMap, 'azimuth [deg]  (+ = left)');
        end

        function drawScene(obj, sep)
            % Top view: listener at the centre; each object in its direction
            % (energy-weighted mean over the excerpt, with its path),
            % objects on their own ring so close directions stay clickable;
            % marker size ~ level
            Ex = obj.Ex; ax = obj.AxScene; cla(ax); hold(ax, 'on');
            t = linspace(0, 2*pi, 100);
            plot(ax, cos(t), sin(t), ':', 'Color', [0.6 0.6 0.6], 'HitTest', 'off');
            plot(ax, 0, 0, 'o', 'MarkerSize', 12, 'MarkerFaceColor', [1 0.81 0.71], 'MarkerEdgeColor', [0.3 0.16 0.05], 'HitTest', 'off');
            plot(ax, [0 0], [0.12 0.3], 'k-', 'LineWidth', 1.5, 'HitTest', 'off');           % facing: front = up
            text(ax, 0, 1.18, 'front', 'HorizontalAlignment', 'center', 'FontSize', 7, 'Color', [0.4 0.4 0.4]);
            nO = size(Ex.Y, 2);
            e = sum(Ex.Y.^2, 1); edb = 10*log10(e/max(e) + 1e-12);
            pos = sep.Pos(:, 1:2);
            rad = ones(1, nO);
            if isprop(sep, 'K'), rad(1:sep.K) = 1 - 0.17*(0:sep.K-1); end
            if isprop(sep, 'TrackLog') && ~isempty(sep.TrackLog)
                TL = sep.TrackLog; fr = sep.N/2; nf = size(TL,2);
                for k = 1:size(TL,1)
                    a = TL(k,:); ok = ~isnan(a);
                    if ~any(ok), continue; end
                    yk = Ex.Y(:,k); w = zeros(1, nf);
                    for f = 1:nf
                        r = (f-1)*fr + (1:fr); r = r(r <= numel(yk)); w(f) = sum(yk(r).^2);
                    end
                    am = sum(a(ok).*w(ok))/max(sum(w(ok)), eps);
                    pos(k,:) = [-sind(am), cosd(am)];
                    plot(ax, -rad(k)*sind(a), rad(k)*cosd(a), '-', 'Color', [Ex.colors(k,:) 0.35], 'LineWidth', 3, 'HitTest', 'off');
                end
            end
            obj.SceneMarks = gobjects(1, nO);
            for k = 1:nO
                p = pos(k,:)/max(norm(pos(k,:)), eps)*rad(k);
                if ~any(pos(k,:)), p = [0 -0.4]; end                % unpositioned (LFE)
                ms = 6 + 14*max(0, 1 + edb(k)/40);
                obj.SceneMarks(k) = plot(ax, p(1), p(2), 'o', 'MarkerSize', ms, ...
                    'MarkerFaceColor', Ex.colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1, ...
                    'ButtonDownFcn', @(~,~) obj.playObject(k));
                if p(1) < -0.05, side = -1; ha = 'right'; else, side = 1; ha = 'left'; end   % label beside
                text(ax, p(1) + side*0.16, p(2), Ex.Names{k}, 'FontSize', 7, 'HorizontalAlignment', ha, ...
                    'ButtonDownFcn', @(~,~) obj.playObject(k));
            end
            hold(ax, 'off'); axis(ax, 'equal'); axis(ax, [-1.5 1.5 -1.5 1.5]); axis(ax, 'off');
            title(ax, 'Sources (click to listen)');
        end

        % ---------------- listening ----------------
        function playObject(obj, k)
            if isempty(obj.Ex), return; end
            obj.highlight(k);
            obj.preview(obj.Ex.Y(:,k), sprintf('Playing %s alone (%.1f s).', obj.Ex.Names{k}, size(obj.Ex.Y,1)/obj.Ex.fs));
        end

        function playOriginal(obj)
            if isempty(obj.Ex), obj.status('Analyze an excerpt first.'); return; end
            obj.highlight(0);
            obj.preview(obj.Ex.X, 'Playing the original excerpt.');
        end

        function playAll(obj)
            if isempty(obj.Ex), obj.status('Analyze an excerpt first.'); return; end
            obj.highlight(0);
            obj.preview(sum(obj.Ex.Y, 2), 'Playing all objects summed (mono, no spatialisation).');
        end

        function preview(obj, y, msg)
            obj.stopPreview();
            try, if obj.C.Engine.Running, obj.C.Engine.stop(); end, catch, end
            g = 0.9/max([abs(obj.Ex.X(:)); abs(obj.Ex.Y(:)); eps]);   % one gain for all: levels comparable
            obj.Player = audioplayer(g*y, obj.Ex.fs);
            play(obj.Player);
            obj.status(msg);
        end

        function stopPreview(obj)
            if ~isempty(obj.Player) && isvalid(obj.Player), stop(obj.Player); end
        end

        function highlight(obj, k)
            for j = 1:numel(obj.LevLines)
                if isgraphics(obj.LevLines(j)), set(obj.LevLines(j), 'LineWidth', 1 + 2*(j == k)); end
                if j <= numel(obj.SceneMarks) && isgraphics(obj.SceneMarks(j))
                    set(obj.SceneMarks(j), 'LineWidth', 1 + 2*(j == k));
                end
            end
        end
    end
end
