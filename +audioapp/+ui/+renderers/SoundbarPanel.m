classdef SoundbarPanel < audioapp.ui.renderers.RendererPanel
    % SoundbarPanel  UI of the soundbar renderer.
    %   - discrete surrounds on/off: standard 5.1 without the centre
    %     (L R Ls Rs + reserved LFE channel); toggling rebuilds the scene
    %   - WFS tapering

    properties
        ChkSurround
        SldTaper
        LblTaper
        TxtInfo
    end

    methods
        function build(obj, parent)
            obj.Root = uipanel('Parent', parent, 'Units', 'normalized', ...
                'Position', [0 0 1 1], 'BorderType', 'none');
            y = 0.82; h = 0.10; gap = 0.03;

            obj.ChkSurround = uicontrol(obj.Root, 'Style', 'checkbox', ...
                'String', 'Discrete surrounds: 5.1 without centre (L, R, Ls, Rs + LFE)', ...
                'Units', 'normalized', 'Position', [0.05 y 0.92 h], ...
                'Callback', @(s,~) obj.onSurround(s));
            y = y - h - gap;

            uicontrol(obj.Root, 'Style', 'text', 'String', 'WFS tapering:', ...
                'Units', 'normalized', 'Position', [0.05 y 0.28 h], ...
                'HorizontalAlignment', 'left');
            obj.SldTaper = uicontrol(obj.Root, 'Style', 'slider', ...
                'Units', 'normalized', 'Position', [0.36 y 0.45 h], ...
                'Min', 0, 'Max', 1, 'Value', 0.25, ...
                'Callback', @(s,~) obj.onTaper(s));
            obj.LblTaper = uicontrol(obj.Root, 'Style', 'text', 'String', '0.25', ...
                'Units', 'normalized', 'Position', [0.83 y 0.14 h], ...
                'HorizontalAlignment', 'left');
            y = y - h - 2*gap;

            obj.TxtInfo = uicontrol(obj.Root, 'Style', 'text', ...
                'Units', 'normalized', 'Position', [0.05 y-0.25 0.90 0.30], ...
                'HorizontalAlignment', 'left', 'FontAngle', 'italic', 'FontSize', 8);
            obj.syncFromSetup();
        end

        function syncFromSetup(obj)
            L = obj.C.Setup.Loudspeaker_setup;
            on = isfield(L, 'Surround51') && L.Surround51;
            set(obj.ChkSurround, 'Value', on);
            S = obj.C.Setup.Renderer_setup;
            if isfield(S, 'Tapering')
                set(obj.SldTaper, 'Value', S.Tapering);
                set(obj.LblTaper, 'String', sprintf('%.2f', S.Tapering));
            end
            info = {'Sources seen through the bar: WFS on the 16 drivers.'};
            try
                r = obj.C.Scene.scene_renderer.SFS_renderer{1};
                if isa(r, 'soundbar_renderer')
                    if on
                        info = {sprintf('WFS window: +-%.0f deg (from the listening position).', r.window), ...
                            'Outside it the source fades from the bar edge into L/R and', ...
                            'on into Ls/Rs (constant-power pairwise panning).'};
                    else
                        info = {sprintf('WFS window: +-%.0f deg. Without surrounds every', r.window), ...
                            'source is rendered by WFS on the bar.'};
                    end
                end
            catch
            end
            set(obj.TxtInfo, 'String', info);
        end
    end

    methods (Access = private)
        function onSurround(obj, src)
            obj.C.Setup.Loudspeaker_setup.Surround51 = logical(get(src, 'Value'));
            obj.C.rebuildSoundScene();          % the loudspeaker set changes
            obj.syncFromSetup();
        end

        function onTaper(obj, src)
            v = get(src, 'Value');
            set(obj.LblTaper, 'String', sprintf('%.2f', v));
            obj.C.Setup.Renderer_setup.Tapering = v;
            try
                obj.C.Scene.scene_renderer.update_renderer_settings(obj.C.Setup);
                obj.C.SceneGUI.draw_renderer_properties(obj.C.Scene.scene_renderer.SFS_renderer{1});
            catch
            end
        end
    end
end
