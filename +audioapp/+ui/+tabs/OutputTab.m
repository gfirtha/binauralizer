classdef OutputTab < handle
    % Output settings: Binauralization toggle, Load HRTF, device picker,
    % device info, block size slider. Uses AudioEngine as the single owner
    % of the audioDeviceWriter.
    properties
        C
        Panel

        TglBinaural   matlab.ui.control.UIControl
        BtnLoadHRTF   matlab.ui.control.UIControl
        PopOutputDev  matlab.ui.control.UIControl
        TxtDevName    matlab.ui.control.UIControl
        TxtMaxCh      matlab.ui.control.UIControl
        TxtFs         matlab.ui.control.UIControl
        TxtBuf        matlab.ui.control.UIControl
        SldBlock      matlab.ui.control.UIControl
        TxtBlock      matlab.ui.control.UIControl
        EQPanel   audioapp.ui.controls.MasterEQPanel
    end

    methods
        function obj = OutputTab(parentTab, ctx)
            obj.C = ctx;
            obj.Panel = parentTab;

            % Make sure the engine (and thus the output device) exists
            obj.C.ensureEngine();
            obj.C.Engine.refreshOutput();

            p = obj.Panel;

            % --- Binauralization + Load HRTF ------------------------------------------------
            obj.TglBinaural = uicontrol(p,'Style','radiobutton','String','Binauralization',...
                'Units','normalized','Position',[0.07 0.92 0.40 0.06],...
                'FontSize',10,'Value',obj.C.Setup.Binauralization,...
                'Callback',@(s,~)obj.onBinauralization(s));

            obj.BtnLoadHRTF = uicontrol(p,'Style','pushbutton','String','Load HRTF',...
                'Units','normalized','Position',[0.52 0.92 0.40 0.06],...
                'FontSize',10,'Callback',@(~,~)obj.onLoadHRTF());

            % --- Device picker ---------------------------------------------------------------
            uicontrol(p,'Style','text','String','Select audio output device:',...
                'Units','normalized','Position',[0.07 0.84 0.86 0.06],...
                'HorizontalAlignment','left','FontSize',10);

            devList = obj.C.Engine.getOutputDevices();
            curDev  = obj.C.Engine.getCurrentDeviceName();
            if ~any(strcmp(devList,curDev)), devList = [{curDev} devList]; end

            obj.PopOutputDev = uicontrol(p,'Style','popupmenu','String',devList,...
                'Units','normalized','Position',[0.07 0.76 0.86 0.06],...
                'FontSize',10,'Callback',@(s,~)obj.onSelectOutputDevice(s));
            set(obj.PopOutputDev,'Value',find(strcmp(devList,curDev),1,'first'));

            % --- Device info fields ----------------------------------------------------------
            uicontrol(p,'Style','text','String','Selected device:',...
                'Units','normalized','Position',[0.07 0.68 0.35 0.05],...
                'HorizontalAlignment','left');
            obj.TxtDevName = uicontrol(p,'Style','text','String','-',...
                'Units','normalized','Position',[0.42 0.68 0.51 0.05],...
                'HorizontalAlignment','left','FontWeight','bold');

            uicontrol(p,'Style','text','String','Max output channels:',...
                'Units','normalized','Position',[0.07 0.61 0.35 0.05],...
                'HorizontalAlignment','left');
            obj.TxtMaxCh = uicontrol(p,'Style','text','String','-',...
                'Units','normalized','Position',[0.42 0.61 0.51 0.05],...
                'HorizontalAlignment','left');

            uicontrol(p,'Style','text','String','Sample rate (Hz):',...
                'Units','normalized','Position',[0.07 0.54 0.35 0.05],...
                'HorizontalAlignment','left');
            obj.TxtFs = uicontrol(p,'Style','text','String','-',...
                'Units','normalized','Position',[0.42 0.54 0.51 0.05],...
                'HorizontalAlignment','left');

            uicontrol(p,'Style','text','String','Buffer size (frames):',...
                'Units','normalized','Position',[0.07 0.47 0.35 0.05],...
                'HorizontalAlignment','left');
            obj.TxtBuf = uicontrol(p,'Style','text','String','-',...
                'Units','normalized','Position',[0.42 0.47 0.51 0.05],...
                'HorizontalAlignment','left');

            % --- Block size slider -----------------------------------------------------------
            uicontrol(p,'Style','text','String','Convolution block size:',...
                'Units','normalized','Position',[0.07 0.37 0.50 0.06],...
                'HorizontalAlignment','left');

            obj.SldBlock = uicontrol(p,'Style','slider','Min',128,'Max',8192,...
                'Value',obj.C.Setup.Block_size,...
                'SliderStep',[128/(8192-128) , 128/(8192-128)],...
                'Units','normalized','Position',[0.07 0.31 0.65 0.05],...
                'Callback',@(s,~)obj.onBlockSize(s));

            obj.TxtBlock = uicontrol(p,'Style','text',...
                'String',num2str(obj.C.Setup.Block_size),...
                'Units','normalized','Position',[0.75 0.30 0.20 0.06],...
                'HorizontalAlignment','left','FontWeight','bold');

            % Initial info fill
            obj.updateDeviceInfo();

            % Master EQ panel (bottom)
            obj.C.ensureEngine();
            obj.EQPanel = audioapp.ui.controls.MasterEQPanel(obj.Panel, obj.C);  % keep a property if you want to refresh later

        end

        % =================================================================
        % Callbacks
        % =================================================================
        function onLoadHRTF(obj)
            startDir = fullfile(pwd,'Data','HRTFs');
            [file,path] = uigetfile(fullfile('HRTFs','*.sofa'),'Select SOFA file',startDir);
            if isequal(file,0), return; end
            sofa = SOFAload(fullfile(path,file));
            obj.C.Setup.HRTF = sofa;

            fs = sofa.Data.SamplingRate;
            obj.C.Setup.SampleRate = fs;

            % Reconfigure hardware & scene to the new Fs
            obj.C.ensureEngine();
            obj.C.Engine.refreshOutput();             % applies Fs/Blk to device
            obj.C.rebuildSoundScene();
            obj.updateDeviceInfo();
        end

        function onBinauralization(obj, src)
            wantOn = logical(get(src,'Value'));   % 1=on (binaural), 0=off (multichannel)
            if ~wantOn
                req   = obj.requiredSpeakerChannels();
                info  = obj.C.Engine.getCurrentDeviceInfo();
                avail = info.MaxOutputChannels;
                if isnan(avail) || avail < req
                    set(src,'Value',1);
                    obj.C.Setup.Binauralization = true;
                    warndlg(sprintf(['Not enough output channels on the selected device.\n\n' ...
                        'Required (speakers in scene): %d\nAvailable on device: %d\n\n' ...
                        'Select another output device or reduce the number of loudspeakers.'], ...
                        req, max(0,avail)), ...
                        'Insufficient output channels','modal');
                    return;
                end
            end
            obj.C.Setup.Binauralization = wantOn;
            obj.C.rebuildSoundScene;
        end

        function onSelectOutputDevice(obj, src)
            if ~obj.C.HasAudio
                errordlg('Audio output not available (Audio Toolbox missing).','No audio output','modal');
                return;
            end
            list = get(src,'String'); idx = get(src,'Value');
            newDev = list{idx};
            if strcmp(newDev,'(no output devices found)'), return; end
            obj.C.ensureEngine();
            obj.C.Engine.setOutputDevice(newDev);
            obj.updateDeviceInfo();
        end

        function onBlockSize(obj, src)
            val = round(get(src,'Value')/128)*128;
            set(src,'Value',val);
            obj.C.Setup.Block_size = val;
            obj.TxtBlock.String    = num2str(val);

            % Apply to device + RB via engine
            obj.C.ensureEngine();
            obj.C.Engine.initIO(obj.C.Setup.SampleRate, val);
            obj.updateDeviceInfo();
            obj.C.rebuildSoundScene;

        end

        % =================================================================
        % UI update
        % =================================================================
        function updateDeviceInfo(obj)
            info = obj.C.Engine.getCurrentDeviceInfo();
            set(obj.TxtDevName,'String', info.Name);
            set(obj.TxtMaxCh,  'String', num2str(info.MaxOutputChannels));
            set(obj.TxtFs,     'String', num2str(info.SampleRate));
            set(obj.TxtBuf,    'String', num2str(info.BufferSize));

            % Keep popup selection aligned to current device
            try
                list = get(obj.PopOutputDev,'String');
                j = find(strcmp(list, info.Name), 1);
                if ~isempty(j), set(obj.PopOutputDev,'Value', j); end
            catch
            end
        end

        % =================================================================
        % Helpers
        % =================================================================
        function n = requiredSpeakerChannels(obj)
            n = obj.C.Setup.Loudspeaker_setup.N;
            switch lower(obj.C.Setup.Loudspeaker_setup.Shape)
                case 'stereo', n = 2;
                case '5.0',    n = 5;
                case '5.1',    n = 6;
            end
        end
    end
end
