classdef TransportBar < handle
    properties
        C
        Panel
        BtnPlay matlab.ui.control.UIControl
        BtnStop matlab.ui.control.UIControl
        SldVol  matlab.ui.control.UIControl
        BtnHarm matlab.ui.control.UIControl
        BtnImp  matlab.ui.control.UIControl
        BtnFR   matlab.ui.control.UIControl
        BtnMon  matlab.ui.control.UIControl   % <--- NEW
    end

    methods
        function obj = TransportBar(parent, ctx)
            obj.C = ctx;
            obj.Panel = uipanel('Parent',parent,'Units','normalized','Position',[0.08 0 0.88 0.1],...
                                'BorderType','etchedin','Title','Transport');

            % Layout
            x = 0.01; w = 0.11; h = 0.60; pad = 0.012;
            w2 = 0.14;                  % base width for action buttons
            fsSmall = 9;                % slightly smaller font for more room

            % Play / Stop / Volume
            obj.BtnPlay = uicontrol(obj.Panel,'Style','pushbutton','String','Play',...
                'Units','normalized','Position',[x 0.25 w h],'Callback',@(s,e)obj.onPlay()); x = x+w+pad;

            obj.BtnStop = uicontrol(obj.Panel,'Style','pushbutton','String','Stop',...
                'Units','normalized','Position',[x 0.25 w h],'Callback',@(s,e)obj.onStop()); x = x+w+pad;

            uicontrol(obj.Panel,'Style','text','String','Master Vol','Units','normalized',...
                'Position',[x 0.80 w 0.25],'HorizontalAlignment','left','FontSize',fsSmall);

            obj.SldVol = uicontrol(obj.Panel,'Style','slider','Min',0,'Max',1,'Value',0.5,...
                'Units','normalized','Position',[x 0.25 w h*0.6],'Callback',@(s,e)obj.onVolumeChanged());  
            x = x+w+pad + 0.05;  % gap after volume

            % Field & analysis buttons (smaller text)
            obj.BtnHarm = uicontrol(obj.Panel,'Style','pushbutton','String','Harmonic field',...
                'Units','normalized','Position',[x 0.25 w2 h],'FontSize',fsSmall,'Callback',@(s,e)obj.onRenderHarmonic()); 
            x = x+w2+pad;

            obj.BtnImp = uicontrol(obj.Panel,'Style','pushbutton','String','Impulsive field',...
                'Units','normalized','Position',[x 0.25 w2 h],'FontSize',fsSmall,'Callback',@(s,e)obj.onRenderImpulsive()); 
            x = x+w2+pad;

            obj.BtnFR = uicontrol(obj.Panel,'Style','pushbutton','String','Frequency response',...
                'Units','normalized','Position',[x 0.25 w2*1.25 h],'FontSize',fsSmall,'Callback',@(s,e)obj.onGetTransfer()); 
            x = x+w2*1.25+pad;

            % NEW: Monitor button (real-time visualization)
            obj.BtnMon = uicontrol(obj.Panel,'Style','pushbutton','String','Monitor',...
                'Units','normalized','Position',[x 0.25 w2*0.5 h],'FontSize',fsSmall,...
                'Callback',@(s,e)obj.onMonitor());

            obj.updateEnabled();
        end

        function updateEnabled(obj)
            canPlay = obj.C.HasDSP && obj.C.HasAudio;
            set(obj.BtnPlay,'Enable', ternary(canPlay,'on','off'));
            set(obj.BtnStop,'Enable', ternary(canPlay,'on','off'));
            function y = ternary(cond,a,b); if cond, y=a; else, y=b; end; end
        end

        function onPlay(obj)
            obj.C.ensureEngine();
            obj.C.Engine.play();
        end

        function onStop(obj)
            if ~isempty(obj.C.Engine), obj.C.Engine.stop(); end
        end

        function onVolumeChanged(obj)
            v = get(obj.SldVol,'Value');
            obj.C.Setup.MasterVolume = v;
        end

        function onRenderHarmonic(obj)
            sim = sound_scene_simulator(obj.C.Scene, obj.C.SceneGUI, 'harmonic', 0.5e3 );
            sim.evaluate_field(0);
        end

        function onRenderImpulsive(obj)
            sim = sound_scene_simulator(obj.C.Scene, obj.C.SceneGUI, 'impulse', 3 );
            sim.evaluate_field(0);
        end

        function onGetTransfer(obj)
            sim = sound_scene_simulator(obj.C.Scene, obj.C.SceneGUI, 'impulse', 1 );
            sim.plot_transfer();
        end

        % NEW: launch monitor-mode simulator (timer-driven ~10 Hz per your class)
        function onMonitor(obj)
            sound_scene_simulator(obj.C.Scene, obj.C.SceneGUI, 'monitor');
        end
    end
end
