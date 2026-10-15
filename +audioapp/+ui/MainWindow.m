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
        TrajectoryTab audioapp.ui.tabs.TrajectoryTab

        ViewMenu   % scene view selector (top / rear / 3D / first person)

        % Top level modules: 1 Binauralizer, 2 Source separation,
        % 3 Measurement-based reproduction, 4 Room builder,
        % 5 Loudspeaker designer, 6 Array lab
        TopTabs
        MainTab
        SepTab
        MeasTab
        RoomTab
        SpeakerTab
        SeparationPanel audioapp.ui.SeparationPanel
        MeasurementPanel audioapp.ui.MeasurementPanel
        RoomBuilderPanel audioapp.ui.RoomBuilderPanel
    end

    properties (Constant)
        ViewModes = {'top','rear','3d','listener'}   % ViewMenu order
    end

    methods
        function obj = MainWindow(ctx)
            obj.C = ctx;

            obj.Fig = figure('Name','Binauralizer','NumberTitle','off','MenuBar','none','Toolbar','none',...
                             'Color',get(0,'defaultuicontrolbackgroundcolor'),...
                             'Units','normalized','Position',[0.12 0.08 0.76 0.84], ...
                             'WindowState','maximized');      % the app always fills the screen

            % --- Top-level tabs ---
            obj.TopTabs = uitabgroup(obj.Fig,'Units','normalized','Position',[0 0 1 1]);
            obj.MainTab = uitab(obj.TopTabs,'Title','  Binauralizer  ');
            obj.SepTab  = uitab(obj.TopTabs,'Title','  Source separation  ');
            obj.MeasTab = uitab(obj.TopTabs,'Title','  Measurement-based reproduction  ');
            obj.RoomTab = uitab(obj.TopTabs,'Title','  Room builder  ');
            obj.SpeakerTab = uitab(obj.TopTabs,'Title','  Loudspeaker designer  ');
            host = obj.MainTab;

            % --- Scene view selector (top bar) ---
            uicontrol(host,'Style','text','String','View:', ...
                'Units','normalized','Position',[0.04 0.917 0.05 0.03], ...
                'FontSize',11,'HorizontalAlignment','right', ...
                'BackgroundColor',get(obj.Fig,'Color'));
            obj.ViewMenu = uicontrol(host,'Style','popupmenu', ...
                'String',{'Top','Rear','3D','Listener (first person)'}, ...
                'Units','normalized','Position',[0.10 0.92 0.16 0.035], ...
                'FontSize',11,'Callback',@(s,e)obj.onViewChanged());

            % Panels
            obj.PnlLeft  = uipanel('Parent',host,'Units','normalized','Position',[0.03 0.05 0.67 0.86],'BorderType','none');
            obj.PnlRight = uipanel('Parent',host,'Units','normalized','Position',[0.70 0.05 0.29 0.90],'BorderType','none');

            % Embed existing axes
            oldAx   = obj.C.SceneGUI.get_main_axes();
            oldFig  = ancestor(oldAx,'figure');
            phAx = axes('Parent',obj.PnlLeft,'Units','normalized','Position',[0.08 0.18 0.88 0.78]);
            set(oldAx,'Parent',obj.PnlLeft,'Units','normalized','Position',get(phAx,'Position'));
            obj.C.SceneGUI.reattachZoomUI();  % keep the +/- panel anchored to this axes
            obj.C.SceneGUI.ViewChangedFcn = @(mode) obj.onSceneView(mode);

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
            obj.Tabs.Trajectory = uitab(tg,'Title','Trajectories  ');

            obj.InputTab    = audioapp.ui.tabs.InputTab(obj.Tabs.Input,    obj.C);
            obj.OutputTab   = audioapp.ui.tabs.OutputTab(obj.Tabs.Output,  obj.C);
            obj.SSDTab      = audioapp.ui.tabs.SSDTab(obj.Tabs.SSD,        obj.C);
            obj.RendererTab = audioapp.ui.tabs.RendererTab(obj.Tabs.Renderer, obj.C);
            obj.TrajectoryTab = audioapp.ui.tabs.TrajectoryTab(obj.Tabs.Trajectory, obj.C);

            % Source separation module (top-level tab)
            obj.SeparationPanel = audioapp.ui.SeparationPanel(obj.SepTab, obj.C);

            % Measurement-based reproduction (spherical-array RIR analysis)
            obj.MeasurementPanel = audioapp.ui.MeasurementPanel(obj.MeasTab, obj.C);

            % Room builder and Loudspeaker designer modules (placeholders)
            obj.RoomBuilderPanel = audioapp.ui.RoomBuilderPanel(obj.RoomTab, obj.C);
            % the Loudspeaker designer is its own app (USER, 2026-10-14): this tab only opens it
            % (run_loudspeaker_app does the same); its "Send to Binauralizer" registers here
            obj.buildSpeakerTab();

            % (the Array lab - arrays measured on a microphone grid - is a step of the Loudspeaker designer)
        end

        function onViewChanged(obj)
            % first person ('listener') is a view of the main axes as well:
            % drag = look around, W/A/S/D walk, arrows turn / tilt, wheel = FOV
            obj.C.SceneGUI.set_view_mode(obj.ViewModes{obj.ViewMenu.Value});
        end

        function onSceneView(obj, mode)
            % any view change (menu, trajectory editor): the menu follows;
            % a film's picture goes to the virtual TV in the 3D views
            k = find(strcmp(obj.ViewModes, mode), 1);
            if ~isempty(k), obj.ViewMenu.Value = k; end
            if ~isempty(obj.C.Video) && isvalid(obj.C.Video)
                obj.C.Video.setVirtual(any(strcmp(mode, {'3d', 'listener'})));
            end
        end
    end

    methods (Access = private)
        function buildSpeakerTab(obj)
            u = {'Units', 'normalized'}; bg = get(obj.Fig, 'Color');
            uicontrol(obj.SpeakerTab, 'Style', 'text', 'String', {'The Loudspeaker designer is a separate app.', '', ...
                ['It opens in its own window (also without the Binauralizer: run_loudspeaker_app). "Send to Binauralizer" ' ...
                'there registers the loudspeaker here (SSD setup tab) and saves it into Data/Loudspeakers.']}, u{:}, ...
                'Position', [0.3 0.5 0.4 0.15], 'BackgroundColor', bg, 'FontSize', 11);
            uicontrol(obj.SpeakerTab, 'Style', 'pushbutton', 'String', 'Open Loudspeaker designer', u{:}, 'Position', [0.4 0.42 0.2 0.06], ...
                'FontSize', 12, 'FontWeight', 'bold', 'Callback', @(~,~) lsapp.LoudspeakerApp.open());
        end
    end
end
