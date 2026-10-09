classdef RendererTab < handle
    properties
        C
        Panel
        PopRenderer
        BtnInfo
        PnlSettings
        CurrentPanel  % audioapp.ui.renderers.RendererPanel (base)
        CurrentModeId (1,:) char = ''   % internal id of the currently shown panel

        SetupLsnr event.listener
    end

    methods
        function obj = RendererTab(parentTab, ctx)
            obj.C = ctx; obj.Panel = parentTab;
            p = obj.Panel;
            obj.SetupLsnr = addlistener(obj.C, 'SetupChanged', @(~,~)obj.refreshFromSetup());

            uicontrol(p,'Style','text','String','Renderer:',...
                'Units','normalized','Position',[0.07 0.90 0.30 0.08],...
                'HorizontalAlignment','left');

            % ── Popup list comes from the registry ──
            reg = audioapp.RendererRegistry.instance();
            modes = reg.uiList();
            val = find(strcmpi(modes, reg.toUI(obj.C.Setup.Rendering_mode)), 1);
            if isempty(val), val = 1; end

            obj.PopRenderer = uicontrol(p,'Style','popupmenu','String',modes,'Value',val,...
                'Units','normalized','Position',[0.40 0.90 0.50 0.08],...
                'Callback', @(s,~)obj.onModeChanged(s));

            obj.BtnInfo = uicontrol(p,'Style','pushbutton','String','Renderer info',...
                'Units','normalized','Position',[0.40 0.88 0.50 0.06],...
                'Callback', @(~,~)obj.onRendererInfo());

            obj.PnlSettings = uipanel('Parent',p,'Units','normalized','Position',[0.02 0.05 0.96 0.8],...
                'BorderType','none');

            obj.rebuildPanel();  % first panel
        end


        function refreshFromSetup(obj)
            % Bail out early if UI not ready
            if isempty(obj.PopRenderer) || ~isvalid(obj.PopRenderer)
                return;
            end

            reg = audioapp.RendererRegistry.instance();
            modeNow = reg.resolveId(obj.C.Setup.Rendering_mode);

            % ---- 1) Sync popup with current mode ----
            list = get(obj.PopRenderer,'String');
            targetUI = reg.toUI(modeNow);
            idx = find(strcmp(list, targetUI), 1);
            if isempty(idx), idx = 1; end
            if get(obj.PopRenderer,'Value') ~= idx
                set(obj.PopRenderer,'Value', idx);
            end

            % ---- 2) Rebuild panel if mode changed ----
            if isempty(obj.CurrentPanel) || ~isvalid(obj.CurrentPanel) ...
                    || ~strcmp(obj.CurrentModeId, modeNow)
                obj.rebuildPanel();
            else
                % Forward a lightweight refresh if the panel supports it
                try
                    if ismethod(obj.CurrentPanel,'syncFromSetup')
                        obj.CurrentPanel.syncFromSetup();
                    end
                catch
                    obj.rebuildPanel();
                end
            end
        end

        function onModeChanged(obj, src)
            list = get(src,'String'); modeUI = list{ get(src,'Value') };
            obj.C.setRendererMode(modeUI);   % central compatibility + rebuild scene
            obj.rebuildPanel();              % swap UI for the new renderer
        end

        function rebuildPanel(obj)
            if ~isempty(obj.CurrentPanel) && isvalid(obj.CurrentPanel), delete(obj.CurrentPanel); end
            reg = audioapp.RendererRegistry.instance();
            modeId = reg.resolveId(obj.C.Setup.Rendering_mode);
            obj.CurrentPanel = reg.createPanel(modeId, obj.PnlSettings, obj.C);
            obj.CurrentModeId = modeId;
        end

        function onRendererInfo(obj)
            try
                info = obj.C.Scene.scene_renderer.get_renderer_info();
                if ~ischar(info) && ~isstring(info), info = evalc('disp(info)'); end
            catch ME
                info = sprintf('Unable to retrieve renderer info.\n\nError:\n%s', ME.message);
            end
            h = figure('Name','Renderer info','NumberTitle','off','MenuBar','none',...
                'ToolBar','none','Units','normalized','Position',[0.25 0.30 0.50 0.40],...
                'Color',[1 1 1], 'WindowStyle','modal');
            uicontrol(h,'Style','text','String','Current renderer information:',...
                'Units','normalized','Position',[0.05 0.90 0.90 0.06],...
                'HorizontalAlignment','left','BackgroundColor',[1 1 1],'FontWeight','bold');
            uicontrol(h,'Style','edit','Max',2,'Min',0,'Enable','inactive',...
                'Units','normalized','Position',[0.05 0.14 0.90 0.74],...
                'HorizontalAlignment','left','String',info,'FontName','Consolas');
            uicontrol(h,'Style','pushbutton','String','Close',...
                'Units','normalized','Position',[0.78 0.04 0.17 0.07],...
                'Callback',@(~,~) close(h));
        end
    end
end
