classdef DBAPPanel < audioapp.ui.renderers.RendererPanel
    % UI panel for Distance-Based Amplitude Panning (DBAP)

    properties
        % Controls (do NOT re-declare Root here; it's inherited)
        EdtA
        SldA
        EdtRs
        SldRs
        Note

        % Keep UI synced with app.Setup
        SetupLsnr event.listener
    end

    methods
        function obj = DBAPPanel(C)
            % Call base class constructor to store context
            obj@audioapp.ui.renderers.RendererPanel(C);
            % Listen for external changes to Setup so UI stays in sync
            obj.SetupLsnr = addlistener(obj.C, 'SetupChanged', @(~,~)obj.syncFromSetup());
        end

        function build(obj, parent)
            % Root panel (owned by base class, but we set it here)
            obj.Root = uipanel('Parent', parent, ...
                'Units','normalized','Position',[0 0 1 1], ...
                'BorderType','etchedin','Title','DBAP controls');

            % Layout helpers
            x1 = 0.07; x2 = 0.55; w1 = 0.45; w2 = 0.30; h = 0.12;
            y  = 0.76;

            % ---- Exponent a -------------------------------------------------------
            uicontrol(obj.Root,'Style','text','String','Exponent  a  (0–3):', ...
                'Units','normalized','Position',[x1 y w1 h],'HorizontalAlignment','left');

            obj.SldA = uicontrol(obj.Root,'Style','slider', ...
                'Min',0,'Max',3,'Value',1.0, ...
                'Units','normalized','Position',[x1 y-0.06 w1 0.06], ...
                'SliderStep',[0.01 0.10], ...
                'Callback',@(s,~)obj.onAChanged(get(s,'Value')));

            obj.EdtA = uicontrol(obj.Root,'Style','edit','String','1.0', ...
                'Units','normalized','Position',[x2 y w2 h], ...
                'Callback',@(s,~)obj.onAChanged(str2double(get(s,'String'))));

            % ---- Regularization r_s ----------------------------------------------
            y = y - 0.26;

            uicontrol(obj.Root,'Style','text','String','Regularization  r_s  (m):', ...
                'Units','normalized','Position',[x1 y w1 h],'HorizontalAlignment','left');

            obj.SldRs = uicontrol(obj.Root,'Style','slider', ...
                'Min',0,'Max',2,'Value',0.1, ...
                'Units','normalized','Position',[x1 y-0.06 w1 0.06], ...
                'SliderStep',[0.01 0.10], ...
                'Callback',@(s,~)obj.onRsChanged(get(s,'Value')));

            obj.EdtRs = uicontrol(obj.Root,'Style','edit','String','0.1', ...
                'Units','normalized','Position',[x2 y w2 h], ...
                'Callback',@(s,~)obj.onRsChanged(str2double(get(s,'String'))));

            % ---- Note -------------------------------------------------------------
            obj.Note = uicontrol(obj.Root,'Style','text', ...
                'String','DBAP normalizes energy: \Sigma g_i^2 = 1  (≈6 dB per distance doubling at a≈1).', ...
                'Units','normalized','Position',[x1 0.08 0.86 0.12], ...
                'HorizontalAlignment','left');

            % Initial sync from current setup
            obj.syncFromSetup();
        end

        function syncFromSetup(obj)
            % Ensure defaults exist
            rs = obj.C.Setup.Renderer_setup;
            if ~isfield(rs,'DBAPExponent'), rs.DBAPExponent = 1.0; end
            if ~isfield(rs,'DBAPRs'),       rs.DBAPRs       = 0.1; end
            obj.C.Setup.Renderer_setup = rs;  % persist

            % Reflect into UI
            if ~isempty(obj.SldA) && isvalid(obj.SldA)
                set(obj.SldA,'Value',rs.DBAPExponent);
                set(obj.EdtA,'String',num2str(rs.DBAPExponent,'%0.3g'));
            end
            if ~isempty(obj.SldRs) && isvalid(obj.SldRs)
                set(obj.SldRs,'Value',rs.DBAPRs);
                set(obj.EdtRs,'String',num2str(rs.DBAPRs,'%0.3g'));
            end
        end

        function delete(obj)
            % Let base delete Root; we just clear our listener
            if ~isempty(obj.SetupLsnr) && isvalid(obj.SetupLsnr)
                delete(obj.SetupLsnr);
            end
        end
    end

    methods (Access=private)
        function onAChanged(obj, v)
            if ~isfinite(v), v = obj.C.Setup.Renderer_setup.DBAPExponent; end
            v = max(0, min(3, v));   % clamp
            obj.C.Setup.Renderer_setup.DBAPExponent = v;
            if ~isempty(obj.SldA) && isvalid(obj.SldA), set(obj.SldA,'Value',v); end
            if ~isempty(obj.EdtA) && isvalid(obj.EdtA), set(obj.EdtA,'String',num2str(v,'%0.3g')); end
            obj.updateRenderer();
        end

        function onRsChanged(obj, v)
            if ~isfinite(v), v = obj.C.Setup.Renderer_setup.DBAPRs; end
            v = max(0, min(2, v));   % clamp
            obj.C.Setup.Renderer_setup.DBAPRs = v;
            if ~isempty(obj.SldRs) && isvalid(obj.SldRs), set(obj.SldRs,'Value',v); end
            if ~isempty(obj.EdtRs) && isvalid(obj.EdtRs), set(obj.EdtRs,'String',num2str(v,'%0.3g')); end
            obj.updateRenderer();
        end

        function updateRenderer(obj)
            % Same contract you use in WFSPanel
            try
                obj.C.Scene.scene_renderer.update_renderer_settings(obj.C.Setup.Renderer_setup);
            catch
                obj.C.rebuildSoundScene();
            end
        end
    end
end
