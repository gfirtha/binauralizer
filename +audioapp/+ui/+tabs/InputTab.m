classdef InputTab < handle
    properties
        C
        Panel
        BtnLoad    matlab.ui.control.UIControl
        BtnMic     matlab.ui.control.UIControl
        LblMic     matlab.ui.control.UIControl
        Table      audioapp.ui.SourceTable
    end

    methods
        function obj = InputTab(parentTab, ctx)
            obj.C = ctx;
            obj.Panel = parentTab;

            % ── Load media button (left half) ──
            obj.BtnLoad = uicontrol(obj.Panel,'Style','pushbutton','String','Load media ...',...
                'Units','normalized','Position',[0.07 0.88 0.42 0.08],...
                'Callback',@(s,e)obj.onLoadMedia());

            % ── Mic input button (right half) ──
            obj.BtnMic = uicontrol(obj.Panel,'Style','pushbutton','String','Mic input ...',...
                'Units','normalized','Position',[0.51 0.88 0.42 0.08],...
                'Callback',@(s,e)obj.onMicInput());

            % ── Label showing current mic device ──
            obj.LblMic = uicontrol(obj.Panel,'Style','text','String','',...
                'Units','normalized','Position',[0.07 0.82 0.86 0.05],...
                'HorizontalAlignment','left','FontSize',8,...
                'ForegroundColor',[0.2 0.5 0.2]);

            obj.Table = audioapp.ui.SourceTable(obj.Panel, [0.02 0.06 0.975 0.73], ctx);
        end

        function onLoadMedia(obj)
            % Clear mic mode if it was active
            obj.C.ensureEngine();
            if obj.C.Engine.UseMic
                obj.C.Engine.clearMicInput();
                set(obj.LblMic,'String','');
            end

            startDir = fullfile(pwd,'Data','SoundSamples');
            [f,p] = uigetfile( ...
                {'*.wav;*.mp3;*.aac;*.ac3;*.flac;*.m4a;*.ogg;*.wma;*.mp4;*.mkv;*.avi;*.webm;*.mov;*.m4v;*.flv;*.wmv;*.ts;*.mpg;*.mpeg', ...
                 'Audio or video files'; ...
                 '*.mp4;*.mkv;*.avi;*.webm;*.mov;*.m4v;*.flv;*.wmv;*.ts;*.mpg;*.mpeg', ...
                 'Video files'; ...
                 '*.wav;*.mp3;*.aac;*.ac3;*.flac;*.m4a;*.ogg;*.wma', ...
                 'Audio files'; ...
                 '*.*', 'All files'}, ...
                'Select media',startDir);
            if isequal(f,0), return; end
            obj.C.ensureEngine();
            kind = obj.C.Engine.loadMedia(fullfile(p,f)); %#ok<NASGU>
            obj.C.rebuildSoundScene();
            if ~isempty(obj.C.Video) && obj.C.Video.hasVideo(), obj.C.Video.showFirstFrame(); end
            obj.Table.refresh();
        end

        function onMicInput(obj)
            obj.C.ensureEngine();

            % ── Get list of available input devices ──
            devList = obj.C.Engine.getInputDevices();
            if isempty(devList) || (numel(devList)==1 && contains(devList{1},'no input'))
                errordlg('No audio input devices found.','Mic Input','modal');
                return;
            end

            % ── Let user pick a device ──
            [sel, ok] = listdlg('ListString', devList, ...
                'SelectionMode','single', ...
                'Name','Select Microphone', ...
                'PromptString','Choose audio input device:', ...
                'ListSize',[350 200]);
            if ~ok, return; end

            devName = devList{sel};

            % ── Let user pick number of channels (1 = mono mic) ──
            answer = inputdlg('Number of input channels:', ...
                'Mic Channels', 1, {'1'});
            if isempty(answer), return; end
            nCh = str2double(answer{1});
            if isnan(nCh) || nCh < 1, nCh = 1; end
            nCh = round(nCh);

            % ── Activate mic input on engine ──
            try
                obj.C.Engine.loadMicInput(devName, nCh);
            catch ME
                errordlg(sprintf('Failed to open mic:\n%s', ME.message), ...
                    'Mic Input Error','modal');
                return;
            end

            % ── Rebuild scene for the new N_in ──
            obj.C.rebuildSoundScene();
            obj.Table.refresh();

            % ── Update label ──
            set(obj.LblMic, 'String', ...
                sprintf('MIC: %s  (%d ch, %d Hz)', devName, nCh, obj.C.Setup.SampleRate));
        end
    end
end
