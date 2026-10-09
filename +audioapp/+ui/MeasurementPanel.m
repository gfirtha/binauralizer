classdef MeasurementPanel < handle
    % MeasurementPanel  Top-level "Measurement-based reproduction" tab
    % (module 5).
    %
    %   Spherical-array room impulse responses (24-ch Firthenmike rigid
    %   sphere) are analysed per time window - the content is M. Bubrik's
    %   thesis work (Files/Objects/MeasurementReproduction):
    %     - conventional max-rE SH beamforming map
    %     - SH-MUSIC (signal subspace L, band)
    %     - Herglotz echo detection and localisation (Massé et al.)
    %   Select a window by dragging on the RIR (or type start / length) and
    %   press Analyze. Planned (disabled): other microphone layouts, and
    %   the reproduction in the Binauralizer (parametric: direct sound +
    %   detected reflections as virtual sources; or SH / Ambisonic
    %   binaural decoding).

    properties
        C
        Root
        TxtFile
        PopArray
        PopMethod
        Edt = struct()
        EdtStart
        EdtLen
        TxtStatus
        AxRIR
        AxMap
        AxAux
        TxtResult
        WinPatch
        Data = []           % struct: h, fs, file, arr
        Result = []
    end

    properties (Constant)
        Methods = {'beamforming', 'music', 'herglotz'}
        MethodNames = {'Conventional SH beamforming (max-rE)', 'SH-MUSIC', 'Herglotz echo detection (Massé)'}
        % name, label, default, methods using it
        Params = {'L', 'MUSIC signal subspace L', 6, {'music'}; ...
                  'fLo', 'MUSIC band low [Hz]', 800, {'music'}; ...
                  'fHi', 'MUSIC band high [Hz]', 1600, {'music'}; ...
                  'kappa', 'Model penalty kappa', 0.75, {'herglotz'}; ...
                  'Pmax', 'Max. echoes', 10, {'herglotz'}; ...
                  'fcBF', 'Beamforming low-pass [Hz]', 3000, {'beamforming'}; ...
                  'fcHG', 'Herglotz low-pass [Hz]', 2100, {'herglotz'}}
        Examples = {'Measured: front (rir_front_trunc)', 'rir_front_trunc.mat'; ...
                    'Simulated: Massé IV-A (rir_mass4_sim)', 'rir_mass4_sim.mat'}
    end

    methods
        function obj = MeasurementPanel(parent, C)
            obj.C = C;
            obj.build(parent);
        end

        function build(obj, parent)
            u = {'Units', 'normalized'};
            bg = get(ancestor(parent, 'figure'), 'Color');
            obj.Root = uipanel('Parent', parent, u{:}, 'Position', [0 0 1 1], 'BorderType', 'none');
            L = uipanel('Parent', obj.Root, u{:}, 'Position', [0.01 0.02 0.27 0.96], ...
                'Title', 'Measurement-based reproduction');
            txt = @(s, y, h) uicontrol(L, 'Style', 'text', 'String', s, u{:}, 'Position', [0.04 y 0.92 h], ...
                'HorizontalAlignment', 'left', 'BackgroundColor', bg);

            % --- data ---
            txt('Room impulse responses:', 0.955, 0.03);
            uicontrol(L, 'Style', 'pushbutton', 'String', 'Load RIR (.mat, h [T x Nch])...', u{:}, ...
                'Position', [0.04 0.905 0.60 0.045], 'Callback', @(~,~) obj.onLoad());
            uicontrol(L, 'Style', 'popupmenu', 'String', [{'Examples...'}, obj.Examples(:,1)'], u{:}, ...
                'Position', [0.66 0.905 0.30 0.045], 'Callback', @(s,~) obj.onExample(s));
            obj.TxtFile = txt('-', 0.865, 0.035); set(obj.TxtFile, 'FontAngle', 'italic');
            txt('Microphone array:', 0.825, 0.03);
            obj.PopArray = uicontrol(L, 'Style', 'popupmenu', ...
                'String', {'Firthenmike 24 ch, rigid sphere r = 10 cm'}, u{:}, ...
                'Position', [0.04 0.785 0.92 0.04], ...
                'TooltipString', 'Other layouts: planned (mbr_mic_array.from_file)');

            % --- analysis ---
            txt('Analysis:', 0.74, 0.03);
            obj.PopMethod = uicontrol(L, 'Style', 'popupmenu', 'String', obj.MethodNames, u{:}, ...
                'Position', [0.04 0.70 0.92 0.04], 'Value', 3, 'Callback', @(~,~) obj.refreshParams());
            y = 0.655;
            for k = 1:size(obj.Params, 1)
                uicontrol(L, 'Style', 'text', 'String', obj.Params{k,2}, u{:}, 'Position', [0.04 y 0.60 0.03], ...
                    'HorizontalAlignment', 'left', 'BackgroundColor', bg);
                obj.Edt.(obj.Params{k,1}) = uicontrol(L, 'Style', 'edit', 'String', num2str(obj.Params{k,3}), ...
                    u{:}, 'Position', [0.66 y 0.30 0.035]);
                y = y - 0.042;
            end
            y = y - 0.01;
            uicontrol(L, 'Style', 'text', 'String', 'Window start / length [ms]:', u{:}, ...
                'Position', [0.04 y 0.56 0.03], 'HorizontalAlignment', 'left', 'BackgroundColor', bg);
            obj.EdtStart = uicontrol(L, 'Style', 'edit', 'String', '0', u{:}, 'Position', [0.62 y 0.16 0.035]);
            obj.EdtLen = uicontrol(L, 'Style', 'edit', 'String', '5', u{:}, 'Position', [0.80 y 0.16 0.035]);
            uicontrol(L, 'Style', 'pushbutton', 'String', 'Analyze window', u{:}, ...
                'Position', [0.04 y-0.06 0.92 0.05], 'FontWeight', 'bold', 'Callback', @(~,~) obj.onAnalyze());
            txt('Tip: drag on the RIR plot to select a window.', y-0.10, 0.03);

            % --- reproduction (planned) ---
            txt('Reproduction in the Binauralizer (planned):', 0.20, 0.03);
            uicontrol(L, 'Style', 'popupmenu', 'String', ...
                {'Parametric: direct + reflections as virtual sources', 'SH / Ambisonic binaural decoding'}, ...
                u{:}, 'Position', [0.04 0.155 0.92 0.04], 'Enable', 'off');
            uicontrol(L, 'Style', 'pushbutton', 'String', 'Send to Binauralizer', u{:}, ...
                'Position', [0.04 0.095 0.92 0.05], 'Enable', 'off');
            obj.TxtStatus = txt('Load a RIR set (or pick an example).', 0.005, 0.08);
            set(obj.TxtStatus, 'ForegroundColor', [0.1 0.3 0.6]);

            % --- plots ---
            obj.AxRIR = axes('Parent', obj.Root, u{:}, 'Position', [0.33 0.72 0.64 0.24], ...
                'ButtonDownFcn', @(~,~) obj.onDragStart());
            title(obj.AxRIR, 'RIR (channel 1) and RMS envelope - drag to select a window');
            xlabel(obj.AxRIR, 't [ms]');
            obj.AxMap = axes('Parent', obj.Root, u{:}, 'Position', [0.33 0.08 0.40 0.54]);
            title(obj.AxMap, 'Direction map'); xlabel(obj.AxMap, 'azimuth [deg]'); ylabel(obj.AxMap, 'zenith [deg]');
            obj.AxAux = axes('Parent', obj.Root, u{:}, 'Position', [0.78 0.38 0.20 0.24]);
            obj.TxtResult = uicontrol(obj.Root, 'Style', 'text', 'String', '', u{:}, ...
                'Position', [0.77 0.04 0.22 0.28], 'HorizontalAlignment', 'left', 'FontName', 'Consolas', ...
                'BackgroundColor', bg);
            obj.refreshParams();
        end
    end

    methods (Access = private)
        % ---------------- data ----------------
        function onLoad(obj)
            [f, p] = uigetfile({'*.mat', 'RIR set (*.mat with h [T x Nch])'}, 'Load room impulse responses');
            if isequal(f, 0), return; end
            obj.loadFile(fullfile(p, f));
        end

        function onExample(obj, src)
            v = get(src, 'Value') - 1;
            if v < 1, return; end
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));   % binauralizer root
            obj.loadFile(fullfile(root, 'Data', 'Measurements', 'firthenmike', obj.Examples{v,2}));
            set(src, 'Value', 1);
        end

        function loadFile(obj, file)
            try
                S = load(file);
                if isfield(S, 'h'), h = S.h; else
                    fn = fieldnames(S); h = S.(fn{find(cellfun(@(n) isnumeric(S.(n)) && ismatrix(S.(n)), fn), 1)});
                end
                fs = 48000;
                if isfield(S, 'fs'), fs = S.fs; elseif isfield(S, 'sim_truth') && isfield(S.sim_truth, 'fs'), fs = S.sim_truth.fs; end
                arr = mbr_mic_array.firthenmike(0.1);
                if size(h, 2) ~= arr.nCh(), h = h.'; end
                if size(h, 2) ~= arr.nCh()
                    error('The RIR has %d channels, the array %d.', size(h,2), arr.nCh());
                end
            catch ME
                obj.status(['Load failed: ' ME.message]); return;
            end
            obj.Data = struct('h', double(h), 'fs', fs, 'file', file, 'arr', arr);
            [~, nm, ex] = fileparts(file);
            set(obj.TxtFile, 'String', sprintf('%s%s  (%d x %d, %.1f ms @ %d Hz)', nm, ex, size(h,1), size(h,2), ...
                size(h,1)/fs*1000, fs));
            obj.drawRIR();
            obj.status('Select a window (drag on the RIR) and press Analyze.');
        end

        function drawRIR(obj)
            d = obj.Data; t = (0:size(d.h,1)-1)'/d.fs*1000;
            ax = obj.AxRIR; cla(ax); hold(ax, 'on');
            plot(ax, t, d.h(:,1), 'Color', [0.3 0.3 0.3], 'HitTest', 'off');
            plot(ax, t, sqrt(mean(d.h.^2, 2)), 'Color', [0.1 0.5 0.8], 'LineWidth', 1.2, 'HitTest', 'off');
            yl = max(abs(d.h(:)))*[-1.1 1.1];
            obj.WinPatch = patch(ax, [0 0 0 0], [yl(1) yl(1) yl(2) yl(2)], [0.9 0.2 0.2], ...
                'FaceAlpha', 0.2, 'EdgeColor', [0.9 0.2 0.2], 'HitTest', 'off', 'Visible', 'off');
            hold(ax, 'off'); grid(ax, 'on'); xlim(ax, [0 t(end)]); ylim(ax, yl);
            set(ax, 'ButtonDownFcn', @(~,~) obj.onDragStart());
            title(ax, 'RIR (channel 1) and RMS envelope - drag to select a window'); xlabel(ax, 't [ms]');
            obj.showWindow();
        end

        % ---------------- window selection ----------------
        function onDragStart(obj)
            if isempty(obj.Data), return; end
            fig = ancestor(obj.AxRIR, 'figure');
            cp = get(obj.AxRIR, 'CurrentPoint'); t0 = cp(1,1);
            old = {get(fig, 'WindowButtonMotionFcn'), get(fig, 'WindowButtonUpFcn')};
            set(fig, 'WindowButtonMotionFcn', @(~,~) obj.onDragMove(t0), ...
                'WindowButtonUpFcn', @(~,~) obj.onDragEnd(t0, old));
        end

        function onDragMove(obj, t0)
            cp = get(obj.AxRIR, 'CurrentPoint');
            obj.setWindow(t0, cp(1,1));
        end

        function onDragEnd(obj, t0, old)
            fig = ancestor(obj.AxRIR, 'figure');
            set(fig, 'WindowButtonMotionFcn', old{1}, 'WindowButtonUpFcn', old{2});
            cp = get(obj.AxRIR, 'CurrentPoint');
            obj.setWindow(t0, cp(1,1));
        end

        function setWindow(obj, ta, tb)
            xl = xlim(obj.AxRIR);
            ta = min(max(ta, xl(1)), xl(2)); tb = min(max(tb, xl(1)), xl(2));
            set(obj.EdtStart, 'String', sprintf('%.2f', min(ta, tb)));
            set(obj.EdtLen, 'String', sprintf('%.2f', max(abs(tb - ta), 0.1)));
            obj.showWindow();
        end

        function win = getWindow(obj)
            t0 = str2double(get(obj.EdtStart, 'String')); len = str2double(get(obj.EdtLen, 'String'));
            if isnan(t0), t0 = 0; end
            if isnan(len) || len <= 0, len = 5; end
            win = [t0, t0 + len];
        end

        function showWindow(obj)
            if isempty(obj.WinPatch) || ~isgraphics(obj.WinPatch), return; end
            w = obj.getWindow();
            set(obj.WinPatch, 'XData', [w(1) w(2) w(2) w(1)], 'Visible', 'on');
        end

        % ---------------- analysis ----------------
        function refreshParams(obj)
            m = obj.Methods{get(obj.PopMethod, 'Value')};
            for k = 1:size(obj.Params, 1)
                set(obj.Edt.(obj.Params{k,1}), 'Enable', ...
                    char(matlab.lang.OnOffSwitchState(any(strcmp(obj.Params{k,4}, m)))));
            end
        end

        function v = param(obj, name)
            v = str2double(get(obj.Edt.(name), 'String'));
        end

        function onAnalyze(obj)
            if isempty(obj.Data), obj.status('Load a RIR set first.'); return; end
            d = obj.Data; win = obj.getWindow(); obj.showWindow();
            m = obj.Methods{get(obj.PopMethod, 'Value')};
            obj.status(sprintf('Analysing %.2f - %.2f ms (%s)...', win, m));
            try
                tic;
                switch m
                    case 'beamforming'
                        R = mbr_beamforming_map(d.h, d.arr, d.fs, win, struct('fcLP', obj.param('fcBF')));
                    case 'music'
                        R = mbr_music_map(d.h, d.arr, d.fs, win, struct('L', obj.param('L'), ...
                            'fBand', [obj.param('fLo') obj.param('fHi')]));
                    case 'herglotz'
                        R = mbr_herglotz_echoes(d.h, d.arr, d.fs, win, struct('kappa', obj.param('kappa'), ...
                            'Pmax', obj.param('Pmax'), 'fcLP', obj.param('fcHG')));
                end
                tc = toc;
            catch ME
                obj.status(['Failed: ' ME.message]); return;
            end
            obj.Result = R;
            obj.drawResult(R);
            obj.status(sprintf('Done in %.1f s: %.2f - %.2f ms.', tc, R.window));
        end

        function drawResult(obj, R)
            ax = obj.AxMap; cla(ax); cla(obj.AxAux);
            switch R.method
                case 'beamforming'
                    imagesc(ax, R.azim*180/pi, R.zen*180/pi, 10*log10(R.map + eps)); clim(ax, [-20 0]);
                    hold(ax, 'on'); plot(ax, R.peak(1), R.peak(2), 'wp', 'MarkerSize', 12, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
                    plot(obj.AxAux, R.polarAz*180/pi, 10*log10(R.polar + eps), 'k');
                    ylim(obj.AxAux, [-20 0]); title(obj.AxAux, 'Horizontal plane [dB]');
                    res = sprintf('Peak:\n az  %.0f deg\n zen %.0f deg', R.peak);
                case 'music'
                    imagesc(ax, R.phi, R.theta, R.map); clim(ax, [-15 0]);
                    hold(ax, 'on'); plot(ax, R.peak(1), R.peak(2), 'kv', 'MarkerSize', 10, 'MarkerFaceColor', [0.9 0.2 0.2]);
                    e = 10*log10(R.eig/max(R.eig(1), eps));
                    stem(obj.AxAux, 1:numel(e), e, 'filled'); hold(obj.AxAux, 'on');
                    plot(obj.AxAux, [R.L R.L] + 0.5, [min(e) 0], 'r--'); hold(obj.AxAux, 'off');
                    title(obj.AxAux, sprintf('Eigenvalues [dB], L = %d', R.L));
                    res = sprintf('DOA:\n az    %.1f deg\n polar %.1f deg\neig1/eig2 %.1f dB\nband %d-%d Hz (%d bins)', ...
                        R.peak, R.eigRatio, R.fBand, R.bins);
                case 'herglotz'
                    imagesc(ax, R.azim*180/pi, R.zen*180/pi, R.map); clim(ax, [-20 0]);
                    hold(ax, 'on');
                    plot(ax, [R.echoes.az], [R.echoes.zen], 'wp', 'MarkerSize', 12, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'k');
                    plot(obj.AxAux, R.polarAz*180/pi, R.polar, 'k', R.polarAz*180/pi, R.polarModel, 'r--');
                    title(obj.AxAux, 'Horizontal: m (k), model (r)');
                    lines = cell(1, R.P);
                    for k = 1:R.P
                        e = R.echoes(k);
                        lines{k} = sprintf('#%-2d az %4.0f zen %4.0f %5.1f dB', k, e.az, e.zen, e.relDb);
                    end
                    res = sprintf('P = %d echoes\nband %.0f-%.0f Hz\n%s', R.P, R.band, strjoin(lines, newline));
            end
            hold(ax, 'off'); axis(ax, 'xy'); set(ax, 'YDir', 'reverse'); colormap(ax, parula); colorbar(ax);
            xlabel(ax, 'azimuth [deg]'); ylabel(ax, 'zenith [deg]');
            title(ax, sprintf('%s, %.2f - %.2f ms', R.method, R.window));
            xlabel(obj.AxAux, ''); grid(obj.AxAux, 'on');
            set(obj.TxtResult, 'String', res);
        end

        function status(obj, msg)
            set(obj.TxtStatus, 'String', msg); drawnow;
        end
    end
end
