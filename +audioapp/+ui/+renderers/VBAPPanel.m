classdef VBAPPanel < audioapp.ui.renderers.RendererPanel
    % VBAPPanel  UI controls for VBAP renderer parameters.
    %   - Spread slider (0°–180°)
    %   - p-norm popup (amplitude / power / custom)
    %   - Spread source count

    properties
        SldSpread
        LblSpreadVal
        PopPNorm
        EdtPNorm
        EdtNSpread
    end

    methods
        function build(obj, parent)
            obj.Root = uipanel('Parent', parent, 'Units', 'normalized', ...
                'Position', [0 0 1 1], 'BorderType', 'none');

            y = 0.82; h = 0.10; lw = 0.28; gap = 0.03;

            % ── Spread slider ──
            uicontrol(obj.Root, 'Style', 'text', 'String', 'Spread (°):', ...
                'Units', 'normalized', 'Position', [0.05 y lw h], ...
                'HorizontalAlignment', 'left');

            obj.SldSpread = uicontrol(obj.Root, 'Style', 'slider', ...
                'Units', 'normalized', 'Position', [0.05+lw+gap y 0.45 h], ...
                'Min', 0, 'Max', 180, 'Value', 0, ...
                'Callback', @(s,~) obj.onSpreadChanged(s));

            obj.LblSpreadVal = uicontrol(obj.Root, 'Style', 'text', ...
                'Units', 'normalized', 'Position', [0.82 y 0.14 h], ...
                'String', '0°', 'HorizontalAlignment', 'left');

            y = y - h - gap;

            % ── p-norm popup ──
            uicontrol(obj.Root, 'Style', 'text', 'String', 'Normalization:', ...
                'Units', 'normalized', 'Position', [0.05 y lw h], ...
                'HorizontalAlignment', 'left');

            obj.PopPNorm = uicontrol(obj.Root, 'Style', 'popupmenu', ...
                'Units', 'normalized', 'Position', [0.05+lw+gap y 0.30 h], ...
                'String', {'Power (p=2)', 'Amplitude (p=1)', 'Custom'}, ...
                'Value', 1, ...
                'Callback', @(s,~) obj.onPNormPopup(s));

            obj.EdtPNorm = uicontrol(obj.Root, 'Style', 'edit', ...
                'Units', 'normalized', 'Position', [0.68 y 0.12 h], ...
                'String', '2.0', 'Enable', 'off', ...
                'Callback', @(s,~) obj.onPNormEdit(s));

            y = y - h - gap;

            % ── Spread sources count ──
            uicontrol(obj.Root, 'Style', 'text', 'String', 'Spread sources:', ...
                'Units', 'normalized', 'Position', [0.05 y lw h], ...
                'HorizontalAlignment', 'left');

            obj.EdtNSpread = uicontrol(obj.Root, 'Style', 'edit', ...
                'Units', 'normalized', 'Position', [0.05+lw+gap y 0.12 h], ...
                'String', '8', ...
                'Callback', @(s,~) obj.onNSpreadEdit(s));

            y = y - h - 2*gap;

            % ── Info text ──
            uicontrol(obj.Root, 'Style', 'text', ...
                'Units', 'normalized', 'Position', [0.05 y 0.90 h*2], ...
                'String', {'Spread distributes the source across multiple', ...
                           'loudspeaker pairs. p=2: constant power, p=1: constant amplitude.'}, ...
                'HorizontalAlignment', 'left', 'FontAngle', 'italic', ...
                'FontSize', 8);
        end

        function syncFromSetup(obj)
            S = obj.C.Setup.Renderer_setup;
            sp = 0; pn = 2; ns = 8;
            if isfield(S, 'VBAPspread'),  sp = S.VBAPspread;  end
            if isfield(S, 'VBAPpNorm'),   pn = S.VBAPpNorm;   end
            if isfield(S, 'VBAPnSpread'), ns = S.VBAPnSpread;  end

            set(obj.SldSpread, 'Value', sp);
            set(obj.LblSpreadVal, 'String', sprintf('%.0f°', sp));
            set(obj.EdtPNorm, 'String', sprintf('%.1f', pn));
            set(obj.EdtNSpread, 'String', sprintf('%d', ns));

            % Sync popup
            if abs(pn - 2) < 0.01
                set(obj.PopPNorm, 'Value', 1);
                set(obj.EdtPNorm, 'Enable', 'off');
            elseif abs(pn - 1) < 0.01
                set(obj.PopPNorm, 'Value', 2);
                set(obj.EdtPNorm, 'Enable', 'off');
            else
                set(obj.PopPNorm, 'Value', 3);
                set(obj.EdtPNorm, 'Enable', 'on');
            end
        end
    end

    methods (Access = private)
        function onSpreadChanged(obj, src)
            val = round(get(src, 'Value'));
            set(obj.LblSpreadVal, 'String', sprintf('%.0f°', val));
            obj.C.Setup.Renderer_setup.VBAPspread = val;
            obj.pushSettings();
        end

        function onPNormPopup(obj, src)
            switch get(src, 'Value')
                case 1  % Power
                    pn = 2;
                    set(obj.EdtPNorm, 'Enable', 'off', 'String', '2.0');
                case 2  % Amplitude
                    pn = 1;
                    set(obj.EdtPNorm, 'Enable', 'off', 'String', '1.0');
                case 3  % Custom
                    set(obj.EdtPNorm, 'Enable', 'on');
                    pn = str2double(get(obj.EdtPNorm, 'String'));
                    if isnan(pn) || pn < 1 || pn > 2
                        pn = 2;
                        set(obj.EdtPNorm, 'String', '2.0');
                    end
            end
            obj.C.Setup.Renderer_setup.VBAPpNorm = pn;
            obj.pushSettings();
        end

        function onPNormEdit(obj, src)
            pn = str2double(get(src, 'String'));
            if isnan(pn) || pn < 1 || pn > 2
                pn = 2;
                set(src, 'String', '2.0');
            end
            obj.C.Setup.Renderer_setup.VBAPpNorm = pn;
            obj.pushSettings();
        end

        function onNSpreadEdit(obj, src)
            ns = round(str2double(get(src, 'String')));
            if isnan(ns) || ns < 4
                ns = 8;
                set(src, 'String', '8');
            end
            obj.C.Setup.Renderer_setup.VBAPnSpread = ns;
            obj.pushSettings();
        end

        function pushSettings(obj)
            % Push updated setup to all renderers
            try
                obj.C.Scene.scene_renderer.update_renderer_settings(obj.C.Setup);
            catch
            end
            % Trigger GUI overlay redraw with the first SFS renderer
            try
                renderers = obj.C.Scene.scene_renderer.SFS_renderer;
                if ~isempty(renderers)
                    obj.C.SceneGUI.draw_renderer_properties(renderers{1});
                end
            catch
            end
        end
    end
end
