classdef TransportBar < handle
    properties
        C
        Panel
        BtnPlay   matlab.ui.control.UIControl
        BtnPause  matlab.ui.control.UIControl
        BtnStop   matlab.ui.control.UIControl
        SldVol    matlab.ui.control.UIControl
        PopSpeed  matlab.ui.control.UIControl
        SldSeek   matlab.ui.control.UIControl
        LblTime   matlab.ui.control.UIControl
        BtnHarm   matlab.ui.control.UIControl
        BtnImp    matlab.ui.control.UIControl
        BtnFR     matlab.ui.control.UIControl
        BtnMon    matlab.ui.control.UIControl
        SeekTimer
    end

    methods
        function obj = TransportBar(parent, ctx)
            obj.C = ctx;
            obj.Panel = uipanel('Parent',parent,'Units','normalized',...
                'Position',[0.08 0 0.88 0.14],...
                'BorderType','etchedin','Title','Transport');

            fs = 9;
            obj.C.ensureEngine();
            if isempty(obj.C.Engine.MasterEQ) || ~isvalid(obj.C.Engine.MasterEQ)
                obj.C.Engine.MasterEQ = audioapp.dsp.MasterEQ(obj.C);
            end

            % ═══════════════════════════════════════════
            %  TOP ROW: Play | Pause | Stop | Vol ____  analysis buttons (right-aligned)
            % ═══════════════════════════════════════════
            wBtn = 0.055; hBtn = 0.32; pad = 0.005;
            yTop = 0.55;
            x = 0.005;

            obj.BtnPlay = uicontrol(obj.Panel,'Style','pushbutton','String','Play',...
                'Units','normalized','Position',[x yTop wBtn hBtn],...
                'FontSize',fs,'Callback',@(s,e)obj.onPlay());
            x = x + wBtn + pad;

            obj.BtnPause = uicontrol(obj.Panel,'Style','pushbutton','String','Pause',...
                'Units','normalized','Position',[x yTop wBtn hBtn],...
                'FontSize',fs,'Callback',@(s,e)obj.onPause());
            x = x + wBtn + pad;

            obj.BtnStop = uicontrol(obj.Panel,'Style','pushbutton','String','Stop',...
                'Units','normalized','Position',[x yTop wBtn hBtn],...
                'FontSize',fs,'Callback',@(s,e)obj.onStop());
            x = x + wBtn + pad;

            % Volume (compact)
            uicontrol(obj.Panel,'Style','text','String','Vol',...
                'Units','normalized','Position',[x yTop+hBtn*0.3 0.025 hBtn*0.7],...
                'HorizontalAlignment','left','FontSize',fs-1);
            obj.SldVol = uicontrol(obj.Panel,'Style','slider','Min',0,'Max',1,'Value',0.5,...
                'Units','normalized','Position',[x+0.025 yTop+hBtn*0.05 0.055 hBtn*0.5],...
                'Callback',@(s,e)obj.onVolumeChanged());
            x = x + 0.025 + 0.055 + 2*pad;

            % Playback speed (slow motion for watching sources move)
            obj.PopSpeed = uicontrol(obj.Panel,'Style','popupmenu',...
                'String',{'1x','0.5x','0.25x','0.125x'},'Value',1,...
                'Units','normalized','Position',[x yTop+hBtn*0.05 0.05 hBtn*0.9],...
                'FontSize',fs,'TooltipString','Playback speed (slower = lower pitch, video follows)',...
                'Callback',@(s,e)obj.onSpeedChanged());

            % Analysis buttons — right-aligned
            wA = 0.12; % wide enough for full text
            xR = 0.995;  % start from right edge

            xR = xR - wA;
            obj.BtnMon = uicontrol(obj.Panel,'Style','pushbutton','String','Monitor',...
                'Units','normalized','Position',[xR yTop wA hBtn],'FontSize',fs,...
                'Callback',@(s,e)obj.onMonitor());

            xR = xR - wA - pad;
            obj.BtnFR = uicontrol(obj.Panel,'Style','pushbutton','String','Frequency response',...
                'Units','normalized','Position',[xR yTop wA hBtn],'FontSize',fs,...
                'Callback',@(s,e)obj.onGetTransfer());

            xR = xR - wA - pad;
            obj.BtnImp = uicontrol(obj.Panel,'Style','pushbutton','String','Impulsive field',...
                'Units','normalized','Position',[xR yTop wA hBtn],'FontSize',fs,...
                'Callback',@(s,e)obj.onRenderImpulsive());

            xR = xR - wA - pad;
            obj.BtnHarm = uicontrol(obj.Panel,'Style','pushbutton','String','Harmonic field',...
                'Units','normalized','Position',[xR yTop wA hBtn],'FontSize',fs,...
                'Callback',@(s,e)obj.onRenderHarmonic());

            % ═══════════════════════════════════════════
            %  BOTTOM ROW: Seek slider (full width) + time label
            % ═══════════════════════════════════════════
            yBot = 0.10;
            hBot = 0.30;
            xSeek = 0.005;
            wTime = 0.09;
            wSeek = 0.99 - wTime - pad - xSeek;  % fill entire row minus time label

            obj.SldSeek = uicontrol(obj.Panel,'Style','slider','Min',0,'Max',1,'Value',0,...
                'Units','normalized','Position',[xSeek yBot wSeek hBot],...
                'Callback',@(s,e)obj.onSeekChanged());

            obj.LblTime = uicontrol(obj.Panel,'Style','text','String','0:00 / 0:00',...
                'Units','normalized','Position',[xSeek+wSeek+pad yBot wTime hBot],...
                'HorizontalAlignment','left','FontSize',fs);

            obj.updateEnabled();

            % Timer: update seek slider + sync VLC play/pause state (~2 Hz)
            obj.SeekTimer = timer('ExecutionMode','fixedRate','Period',0.5,...
                'TimerFcn',@(~,~)obj.onTimerTick());
            start(obj.SeekTimer);

            % Keyboard: Ctrl+Enter → VLC fullscreen, Space → play/pause
            fig = ancestor(parent, 'figure');
            if ~isempty(fig)
                set(fig, 'WindowKeyPressFcn', @(s,e)obj.onKeyPress(e));
            end
        end

        function delete(obj)
            if ~isempty(obj.SeekTimer)
                try stop(obj.SeekTimer); catch, end
                try delete(obj.SeekTimer); catch, end
            end
        end

        function updateEnabled(obj)
            canPlay = obj.C.HasDSP && obj.C.HasAudio;
            onoff = @(c) ternary(c,'on','off');
            set(obj.BtnPlay,  'Enable', onoff(canPlay));
            set(obj.BtnPause, 'Enable', onoff(canPlay));
            set(obj.BtnStop,  'Enable', onoff(canPlay));
            function y = ternary(cond,a,b); if cond, y=a; else, y=b; end; end
        end

        % ── Transport callbacks ──
        function onPlay(obj)
            obj.C.ensureEngine();
            if obj.C.Engine.is_paused
                obj.C.Engine.resumePlayback();
                if ~isempty(obj.C.Video) && obj.C.Video.hasVideo()
                    obj.C.Video.vlcResume();
                end
            else
                obj.C.Engine.play();
            end
        end

        function onPause(obj)
            if ~isempty(obj.C.Engine)
                obj.C.Engine.pausePlayback();
                % VLC pause is handled inside AudioEngine play loop
            end
        end

        function onStop(obj)
            if ~isempty(obj.C.Engine), obj.C.Engine.stop(); end
        end

        function onVolumeChanged(obj)
            obj.C.Setup.MasterVolume = get(obj.SldVol, 'Value');
        end

        function onSpeedChanged(obj)
            obj.C.ensureEngine();
            str = get(obj.PopSpeed, 'String');
            obj.C.Engine.setSpeed(str2double(erase(str{get(obj.PopSpeed, 'Value')}, 'x')));
        end

        function onSeekChanged(obj)
            if isempty(obj.C.Engine), return; end
            frac = get(obj.SldSeek, 'Value');
            dur = obj.C.Engine.getDuration();
            if dur > 0
                obj.C.Engine.seekTo(frac * dur);
            end
        end

        function onKeyPress(obj, evtData)
            mods = evtData.Modifier;
            hasCtrl = any(strcmp(mods, 'control'));

            switch evtData.Key
                case 'return'
                    % Ctrl+Enter → VLC fullscreen
                    if hasCtrl && ~isempty(obj.C.Video) && obj.C.Video.hasVideo()
                        obj.C.Video.toggleFullscreen();
                    end
                case 'space'
                    % Space → toggle play/pause
                    if ~isempty(obj.C.Engine)
                        if obj.C.Engine.is_paused
                            obj.onPlay();
                        else
                            obj.onPause();
                        end
                    end
            end
        end

        function onTimerTick(obj)
            % 1) Update seek slider and time label
            try
                if isempty(obj.C.Engine), return; end
                dur = obj.C.Engine.getDuration();
                cur = obj.C.Engine.getCurrentTime();

                if dur > 0 && obj.C.Engine.seek_to < 0
                    set(obj.SldSeek, 'Value', min(1, max(0, cur/dur)));
                end
                if obj.C.Engine.Speed < 1      % slow motion: tenths of a second
                    set(obj.LblTime, 'String', sprintf('%d:%04.1f / %s  (%gx)', ...
                        floor(cur/60), mod(cur, 60), fmtTime(dur), obj.C.Engine.Speed));
                else
                    set(obj.LblTime, 'String', sprintf('%s / %s', ...
                        fmtTime(cur), fmtTime(dur)));
                end
            catch
            end

            % 2) Sync VLC external toggles (user pressed Space in VLC window)
            try
                if ~isempty(obj.C.Video) && obj.C.Video.hasVideo() ...
                        && obj.C.Video.HasVLC && ~isempty(obj.C.Video.VLCSock)
                    changed = obj.C.Video.syncFromVLC();
                    if changed
                        % VLC state changed externally — mirror to AudioEngine
                        if obj.C.Video.VLCIsPlaying && obj.C.Engine.is_paused
                            obj.C.Engine.resumePlayback();
                        elseif ~obj.C.Video.VLCIsPlaying && ~obj.C.Engine.is_paused
                            obj.C.Engine.pausePlayback();
                        end
                    end
                end
            catch
            end

            function s = fmtTime(sec)
                sec = max(0, round(sec));
                m = floor(sec / 60);
                sr = mod(sec, 60);
                s = sprintf('%d:%02d', m, sr);
            end
        end

        % ── Analysis buttons ──
        function onRenderHarmonic(obj)
            % the constructor computes and plots the field
            sound_scene_simulator(obj.C.Scene, obj.C.SceneGUI, 'harmonic', 0.5e3);
        end

        function onRenderImpulsive(obj)
            obj.C.ensureEngine();
            if isempty(obj.C.Engine.MasterEQ) || ~isvalid(obj.C.Engine.MasterEQ)
                obj.C.Engine.MasterEQ = audioapp.dsp.MasterEQ(obj.C);
            end
            sound_scene_simulator(obj.C.Scene, obj.C.SceneGUI, 'impulse', 3, obj.C.Engine.MasterEQ);
        end

        function onGetTransfer(obj)
            obj.C.ensureEngine();
            if isempty(obj.C.Engine.MasterEQ) || ~isvalid(obj.C.Engine.MasterEQ)
                obj.C.Engine.MasterEQ = audioapp.dsp.MasterEQ(obj.C);
            end
            sim = sound_scene_simulator(obj.C.Scene, obj.C.SceneGUI, 'impulse', 1, obj.C.Engine.MasterEQ);
            sim.plot_transfer();
        end

        function onMonitor(obj)
            sound_scene_simulator(obj.C.Scene, obj.C.SceneGUI, 'monitor');
        end
    end
end
