classdef MainWindow < handle
    properties
        C
        Fig
        PnlLeft
        PnlRight
        Ax
        Tabs
        Transport   audioapp.ui.TransportBar
        InputTab    audioapp.ui.tabs.InputTab
        OutputTab   audioapp.ui.tabs.OutputTab
        SSDTab      audioapp.ui.tabs.SSDTab
        RendererTab audioapp.ui.tabs.RendererTab

        Btn3D   % <-- NEW: 3D mode button handle
    end

    methods
        function obj = MainWindow(ctx)
            obj.C = ctx;

            obj.Fig = figure('Name','Binauralizer','NumberTitle','off','MenuBar','none','Toolbar','none',...
                             'Color',get(0,'defaultuicontrolbackgroundcolor'),...
                             'Units','normalized','Position',[0.12 0.08 0.76 0.84]);

            % --- NEW: top-bar 3D mode button (like the original GUI) ----
            obj.Btn3D = uicontrol(obj.Fig,'Style','pushbutton','String','3D mode', ...
                'Units','normalized','Position',[0.10 0.92 0.16 0.035], ...
                'FontSize',11,'Callback',@(s,e)obj.on3DMode());

            % Panels
            obj.PnlLeft  = uipanel('Parent',obj.Fig,'Units','normalized','Position',[0.03 0.05 0.67 0.86],'BorderType','none');
            obj.PnlRight = uipanel('Parent',obj.Fig,'Units','normalized','Position',[0.70 0.05 0.29 0.90],'BorderType','none');

            % Embed existing axes
            oldAx   = obj.C.SceneGUI.get_main_axes();
            oldFig  = ancestor(oldAx,'figure');
            phAx = axes('Parent',obj.PnlLeft,'Units','normalized','Position',[0.08 0.18 0.88 0.78]);
            set(oldAx,'Parent',obj.PnlLeft,'Units','normalized','Position',get(phAx,'Position'));
            obj.C.SceneGUI.reattachZoomUI();  % keep the +/- panel anchored to this axes

            delete(phAx);
            obj.Ax = oldAx;
            if isvalid(oldFig) && oldFig ~= obj.Fig
                try, set(oldFig,'CloseRequestFcn','delete(gcf)'); catch, end
                try, delete(oldFig); catch, end
            end

            % Transport
            obj.Transport = audioapp.ui.TransportBar(obj.PnlLeft, obj.C);

            % Tabs
            tg = uitabgroup(obj.PnlRight,'Units','normalized','Position',[0 0 1 1]);
            obj.Tabs.Input    = uitab(tg,'Title','Input setup  ');
            obj.Tabs.Output   = uitab(tg,'Title','Output setup  ');
            obj.Tabs.SSD      = uitab(tg,'Title','SSD setup  ');
            obj.Tabs.Renderer = uitab(tg,'Title','Renderer setup  ');

            obj.InputTab    = audioapp.ui.tabs.InputTab(obj.Tabs.Input,    obj.C);
            obj.OutputTab   = audioapp.ui.tabs.OutputTab(obj.Tabs.Output,  obj.C);
            obj.SSDTab      = audioapp.ui.tabs.SSDTab(obj.Tabs.SSD,        obj.C);
            obj.RendererTab = audioapp.ui.tabs.RendererTab(obj.Tabs.Renderer, obj.C);
        end

        % ---- NEW: 3D mode callback (mirrors your old onTopButton) -------
        function on3DMode(obj)
            try
                if isprop(obj.C,'SceneGUI') && ~isempty(obj.C.SceneGUI) && ...
                        isprop(obj.C.SceneGUI,'room') && ~isempty(obj.C.SceneGUI.room)
                    room_view(obj.C.Scene, obj.C.SceneGUI, obj.C.SceneGUI.room);
                else
                    room_view(obj.C.Scene, obj.C.SceneGUI);
                end
            catch ME
                warndlg(sprintf('Unable to open 3D view:\n\n%s', ME.message), ...
                        '3D mode', 'modal');
            end
        end
    end
end
