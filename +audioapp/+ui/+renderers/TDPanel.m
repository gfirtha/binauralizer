classdef TDPanel < audioapp.ui.renderers.RendererPanel
    % TDPanel  Time-Delay Stereo renderer panel with live panning curve.

    properties (Constant)
        Mode = 'time_delay';
    end

    properties
        PopMethod
        EdtTaps
        EdtRefX
        EdtRefY
        PopR0Mode
        EdtR0
        SldAlpha        % hybrid alpha slider
        LblAlphaVal     % alpha value label
        LblAlphaRow     % "Mix (α):" label
        AxPanning       % axes for panning curve
        LineCurve       % line handle
        LineMarker      % marker for current VS position
        LblDelay        % current delay readout
    end

    methods
        function obj = TDPanel(C)
            obj@audioapp.ui.renderers.RendererPanel(C);
        end

        function build(obj, parent)
            obj.Root = uipanel('Parent', parent, 'Units', 'normalized', ...
                'Position', [0 0 1 1], 'BorderType', 'none');

            fs = 9; x1 = 0.04; x2 = 0.40; w1 = 0.35; w2 = 0.22; h = 0.050; dy = 0.058;
            y = 0.94;

            % ── Method selector ──
            uicontrol(obj.Root, 'Style', 'text', 'String', 'Panning law:', ...
                'Units', 'normalized', 'Position', [x1 y w1 h], ...
                'HorizontalAlignment', 'left', 'FontSize', fs);
            methods = time_delay_renderer.METHODS;
            obj.PopMethod = uicontrol(obj.Root, 'Style', 'popupmenu', ...
                'Units', 'normalized', 'Position', [x2 y 0.56 h], ...
                'String', methods(:,2), 'FontSize', fs, ...
                'Callback', @(s,~) obj.onMethod(s));
            y = y - dy;

            % ── FIR taps ──
            uicontrol(obj.Root, 'Style', 'text', 'String', 'FIR taps:', ...
                'Units', 'normalized', 'Position', [x1 y w1 h], ...
                'HorizontalAlignment', 'left', 'FontSize', fs);
            obj.EdtTaps = uicontrol(obj.Root, 'Style', 'edit', ...
                'Units', 'normalized', 'Position', [x2 y w2 h], ...
                'FontSize', fs, 'Callback', @(s,~) obj.onTaps(s));
            y = y - dy;

            % ── Ref point ──
            uicontrol(obj.Root, 'Style', 'text', 'String', 'Ref point [x, y]:', ...
                'Units', 'normalized', 'Position', [x1 y w1 h], ...
                'HorizontalAlignment', 'left', 'FontSize', fs);
            obj.EdtRefX = uicontrol(obj.Root, 'Style', 'edit', ...
                'Units', 'normalized', 'Position', [x2 y 0.12 h], ...
                'FontSize', fs, 'Callback', @(s,~) obj.onRefX(s));
            obj.EdtRefY = uicontrol(obj.Root, 'Style', 'edit', ...
                'Units', 'normalized', 'Position', [x2+0.14 y 0.12 h], ...
                'FontSize', fs, 'Callback', @(s,~) obj.onRefY(s));
            y = y - dy;

            % ── R0 ──
            uicontrol(obj.Root, 'Style', 'text', 'String', 'R0 (ref radius):', ...
                'Units', 'normalized', 'Position', [x1 y w1 h], ...
                'HorizontalAlignment', 'left', 'FontSize', fs);
            obj.PopR0Mode = uicontrol(obj.Root, 'Style', 'popupmenu', ...
                'Units', 'normalized', 'Position', [x2 y 0.16 h], ...
                'String', {'auto', 'manual'}, 'FontSize', fs, ...
                'Callback', @(s,~) obj.onR0Mode(s));
            obj.EdtR0 = uicontrol(obj.Root, 'Style', 'edit', ...
                'Units', 'normalized', 'Position', [x2+0.18 y 0.12 h], ...
                'FontSize', fs, 'Callback', @(s,~) obj.onR0Value(s));
            y = y - dy;

            % ── Hybrid alpha slider (visible only for 'hybrid' method) ──
            obj.LblAlphaRow = uicontrol(obj.Root, 'Style', 'text', ...
                'String', 'Mix (alpha):', ...
                'Units', 'normalized', 'Position', [x1 y w1 h], ...
                'HorizontalAlignment', 'left', 'FontSize', fs, 'Visible', 'off');
            obj.SldAlpha = uicontrol(obj.Root, 'Style', 'slider', ...
                'Units', 'normalized', 'Position', [x2 y 0.38 h], ...
                'Min', 0, 'Max', 1, 'Value', 0.5, 'Visible', 'off', ...
                'Callback', @(s,~) obj.onAlpha(s));
            obj.LblAlphaVal = uicontrol(obj.Root, 'Style', 'text', ...
                'Units', 'normalized', 'Position', [x2+0.40 y 0.20 h], ...
                'String', '0.50 (delay+amp)', 'FontSize', fs-1, ...
                'HorizontalAlignment', 'left', 'Visible', 'off');
            y = y - dy;

            % ── Current delay readout ──
            obj.LblDelay = uicontrol(obj.Root, 'Style', 'text', ...
                'Units', 'normalized', 'Position', [x1 y 0.92 h], ...
                'String', '', 'HorizontalAlignment', 'left', ...
                'FontSize', fs, 'FontAngle', 'italic');
            y = y - 0.01;

            % ── Panning curve axes (x = delay, y = angle) ──
            obj.AxPanning = axes('Parent', obj.Root, ...
                'Units', 'normalized', 'Position', [0.15 0.03 0.78 y - 0.03], ...
                'FontSize', 8, 'Box', 'on', 'NextPlot', 'add');
            xlabel(obj.AxPanning, '\Deltat (ms)', 'FontSize', 8);
            ylabel(obj.AxPanning, 'VS angle from bisector (°)', 'FontSize', 8);
            title(obj.AxPanning, 'Panning law', 'FontSize', 9);
            grid(obj.AxPanning, 'on');

            % Line: panning curve (x = delay, y = angle)
            obj.LineCurve = plot(obj.AxPanning, NaN, NaN, '-', ...
                'LineWidth', 2, 'Color', [0.1 0.4 0.8]);
            % Marker: current VS position
            obj.LineMarker = plot(obj.AxPanning, NaN, NaN, 'o', ...
                'MarkerSize', 10, 'MarkerFaceColor', [0.9 0.2 0.2], ...
                'MarkerEdgeColor', [0.6 0 0], 'LineWidth', 1.5);

            obj.syncFromSetup();
        end

        function syncFromSetup(obj)
            RS = obj.safeRS(obj.C.Setup);

            methods = time_delay_renderer.METHODS;
            ix = find(strcmp(methods(:,1), RS.TDMethod), 1);
            if isempty(ix), ix = 2; end
            set(obj.PopMethod, 'Value', ix);

            set(obj.EdtTaps, 'String', num2str(RS.TDNumTaps));
            set(obj.EdtRefX, 'String', num2str(RS.TDRefPoint(1)));
            set(obj.EdtRefY, 'String', num2str(RS.TDRefPoint(2)));

            if isempty(RS.TDRefRadiusR0) || (ischar(RS.TDRefRadiusR0) && strcmpi(RS.TDRefRadiusR0, 'auto'))
                set(obj.PopR0Mode, 'Value', 1);
                set(obj.EdtR0, 'String', '', 'Enable', 'off');
            else
                set(obj.PopR0Mode, 'Value', 2);
                set(obj.EdtR0, 'String', num2str(RS.TDRefRadiusR0), 'Enable', 'on');
            end

            % Alpha slider (hybrid only)
            isHybrid = strcmp(RS.TDMethod, 'hybrid');
            obj.setAlphaVisibility(isHybrid);
            alpha = RS.TDHybridAlpha;
            set(obj.SldAlpha, 'Value', alpha);
            obj.updateAlphaLabel(alpha, obj.LblAlphaVal);

            % Register panel callback on the renderer so the plot
            % updates whenever VS moves or settings change
            r = obj.getTDRenderer();
            if ~isempty(r)
                r.panelCallback = @() obj.updatePanningPlot();
            end

            obj.updatePanningPlot();
        end
    end

    methods (Access = private)
        function onMethod(obj, src)
            methods = time_delay_renderer.METHODS;
            m = methods{get(src,'Value'), 1};
            obj.C.Setup.Renderer_setup.TDMethod = m;
            obj.setAlphaVisibility(strcmp(m, 'hybrid'));
            obj.pushAndPlot();
        end
        function onTaps(obj, src)
            v = str2double(get(src, 'String'));
            if isfinite(v) && v > 8 && mod(v,1) == 0
                obj.C.Setup.Renderer_setup.TDNumTaps = v;
                obj.pushAndPlot();
            else, obj.syncFromSetup(); end
        end
        function onRefX(obj, src)
            v = str2double(get(src, 'String'));
            if isfinite(v)
                rp = obj.safeRS(obj.C.Setup).TDRefPoint; rp(1) = v;
                obj.C.Setup.Renderer_setup.TDRefPoint = rp;
                obj.pushAndPlot();
            else, obj.syncFromSetup(); end
        end
        function onRefY(obj, src)
            v = str2double(get(src, 'String'));
            if isfinite(v)
                rp = obj.safeRS(obj.C.Setup).TDRefPoint; rp(2) = v;
                obj.C.Setup.Renderer_setup.TDRefPoint = rp;
                obj.pushAndPlot();
            else, obj.syncFromSetup(); end
        end
        function onR0Mode(obj, src)
            if get(src, 'Value') == 1
                obj.C.Setup.Renderer_setup.TDRefRadiusR0 = 'auto';
                set(obj.EdtR0, 'Enable', 'off');
            else
                RS = obj.safeRS(obj.C.Setup);
                if isempty(RS.TDRefRadiusR0) || ~isnumeric(RS.TDRefRadiusR0)
                    obj.C.Setup.Renderer_setup.TDRefRadiusR0 = 1.0;
                end
                set(obj.EdtR0, 'Enable', 'on');
            end
            obj.pushAndPlot(); obj.syncFromSetup();
        end
        function onR0Value(obj, src)
            v = str2double(get(src, 'String'));
            if isfinite(v) && v > 0
                obj.C.Setup.Renderer_setup.TDRefRadiusR0 = v;
                obj.pushAndPlot();
            else, obj.syncFromSetup(); end
        end

        function onAlpha(obj, src)
            v = get(src, 'Value');
            obj.C.Setup.Renderer_setup.TDHybridAlpha = v;
            obj.updateAlphaLabel(v, obj.LblAlphaVal);
            obj.pushAndPlot();
        end

        function updateAlphaLabel(~, v, lbl)
            % Show descriptive text depending on alpha value
            if nargin < 3, return; end
            if v < 0.05
                desc = 'pure delay';
            elseif v > 0.95
                desc = 'pure amplitude';
            else
                desc = sprintf('%.0f%% delay + %.0f%% amp', (1-v)*100, v*100);
            end
            set(lbl, 'String', sprintf('%.2f (%s)', v, desc));
        end

        function setAlphaVisibility(obj, isHybrid)
            vis = 'off';
            if isHybrid, vis = 'on'; end
            set(obj.LblAlphaRow, 'Visible', vis);
            set(obj.SldAlpha,    'Visible', vis);
            set(obj.LblAlphaVal, 'Visible', vis);
        end

        function pushAndPlot(obj)
            try
                obj.C.Scene.scene_renderer.update_renderer_settings(obj.C.Setup);
            catch
            end
            % Also trigger overlay redraw
            try
                renderers = obj.C.Scene.scene_renderer.SFS_renderer;
                if ~isempty(renderers)
                    obj.C.SceneGUI.draw_renderer_properties(renderers{1});
                end
            catch
            end
            obj.updatePanningPlot();
        end

        function updatePanningPlot(obj)
            try
                r = obj.getTDRenderer();
                if isempty(r), return; end

                [ang_deg, dt_ms, g1_vec, g2_vec] = r.get_panning_curve();
                if isempty(ang_deg), return; end

                % Determine method type
                mtype = r.get_method_type();

                switch mtype
                    case 'delay'
                        % x = delay (ms), y = angle (deg)
                        xdata = dt_ms;
                        xlabel(obj.AxPanning, '\Deltat (ms)', 'FontSize', 8);
                    case 'amplitude'
                        % x = gain difference (dB), y = angle (deg)
                        g1_dB = 20*log10(abs(g1_vec) + 1e-12);
                        g2_dB = 20*log10(abs(g2_vec) + 1e-12);
                        xdata = g1_dB - g2_dB;  % ICLD
                        xlabel(obj.AxPanning, '\DeltaL (dB)', 'FontSize', 8);
                    case 'hybrid'
                        % Show delay on x-axis (primary cue)
                        xdata = dt_ms;
                        xlabel(obj.AxPanning, '\Deltat (ms)', 'FontSize', 8);
                end

                set(obj.LineCurve, 'XData', xdata, 'YData', ang_deg);

                % Current VS marker
                fi0_deg = obj.getCurrentVSAngle(r);
                if ~isnan(fi0_deg)
                    x_current = interp1(ang_deg, xdata, fi0_deg, 'linear', 0);
                    set(obj.LineMarker, 'XData', x_current, 'YData', fi0_deg);

                    % Status text
                    switch mtype
                        case 'delay'
                            set(obj.LblDelay, 'String', ...
                                sprintf('VS: %.1f°,  Dt = %.3f ms', fi0_deg, x_current));
                        case 'amplitude'
                            set(obj.LblDelay, 'String', ...
                                sprintf('VS: %.1f°,  DL = %.1f dB', fi0_deg, x_current));
                        case 'hybrid'
                            dL = interp1(ang_deg, 20*log10(abs(g1_vec)+1e-12) - 20*log10(abs(g2_vec)+1e-12), fi0_deg, 'linear', 0);
                            set(obj.LblDelay, 'String', ...
                                sprintf('VS: %.1f°,  Dt = %.3f ms,  DL = %.1f dB', fi0_deg, x_current, dL));
                    end
                end

                % Axis limits
                xl = max(abs(xdata)) * 1.2;
                if xl > 0, xlim(obj.AxPanning, [-xl xl]); end
                ylim(obj.AxPanning, [min(ang_deg) max(ang_deg)]);
            catch
            end
        end

        function r = getTDRenderer(obj)
            r = [];
            try
                renderers = obj.C.Scene.scene_renderer.SFS_renderer;
                if ~isempty(renderers) && isa(renderers{1}, 'time_delay_renderer')
                    r = renderers{1};
                end
            catch
            end
        end

        function fi0_deg = getCurrentVSAngle(~, r)
            fi0_deg = NaN;
            try
                xs = r.virtual_source.position;
                xr = r.refPoint;
                x0 = r.ssd.x0;
                ixs = r.stereo_pair;
                if numel(ixs) < 2, return; end

                base_vec = (x0(ixs(1),:) + x0(ixs(2),:) - 2*xr);
                nb = norm(base_vec); if nb > 0, base_vec = base_vec/nb; end
                vec_xs = xs - xr;
                nv = norm(vec_xs); if nv > 0, vec_xs = vec_xs/nv; end

                fi0 = acos(max(-1, min(1, dot(base_vec, vec_xs))));
                cross_z = base_vec(1)*vec_xs(2) - base_vec(2)*vec_xs(1);
                fi0_deg = sign(cross_z) * fi0 * 180/pi;
            catch
            end
        end

        function RS = safeRS(~, Setup)
            RS = Setup.Renderer_setup;
            if ~isfield(RS, 'TDNumTaps')      || isempty(RS.TDNumTaps),      RS.TDNumTaps = 2048; end
            if ~isfield(RS, 'TDMethod')       || isempty(RS.TDMethod),       RS.TDMethod = 'firtha_linear'; end
            if ~isfield(RS, 'TDRefPoint')     || numel(RS.TDRefPoint) < 2,   RS.TDRefPoint = [0 0 0]; end
            if ~isfield(RS, 'TDRefRadiusR0'),  RS.TDRefRadiusR0 = 'auto'; end
            if ~isfield(RS, 'TDHybridAlpha')  || isempty(RS.TDHybridAlpha),  RS.TDHybridAlpha = 0.5; end
        end
    end
end
