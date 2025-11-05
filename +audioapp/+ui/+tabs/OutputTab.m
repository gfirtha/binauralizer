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
        SldDecor      matlab.ui.control.UIControl % <-- New Property
        TxtDecor      matlab.ui.control.UIControl % <-- New Property
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
            
            % Define a standard height and gap for a denser layout
            h = 0.05;    % Standard height for most controls
            y = 0.98;    % Current vertical top position
            
            % --- Binauralization + Load HRTF ------------------------------------------------
            y = y - h;
            obj.TglBinaural = uicontrol(p,'Style','radiobutton','String','Binauralization',...
                'Units','normalized','Position',[0.07 y 0.40 h],...
                'FontSize',10,'Value',obj.C.Setup.Binauralization,...
                'Callback',@(s,~)obj.onBinauralization(s));
            obj.BtnLoadHRTF = uicontrol(p,'Style','pushbutton','String','Load HRTF',...
                'Units','normalized','Position',[0.52 y 0.40 h],...
                'FontSize',10,'Callback',@(~,~)obj.onLoadHRTF());
            
            % --- Device picker ---------------------------------------------------------------
            y = y - h - 0.02; % Section gap
            uicontrol(p,'Style','text','String','Select audio output device:',...
                'Units','normalized','Position',[0.07 y 0.86 h],...
                'HorizontalAlignment','left','FontSize',10);
            
            y = y - h;
            devList = obj.C.Engine.getOutputDevices();
            curDev  = obj.C.Engine.getCurrentDeviceName();
            if ~any(strcmp(devList,curDev)), devList = [{curDev} devList]; end
            obj.PopOutputDev = uicontrol(p,'Style','popupmenu','String',devList,...
                'Units','normalized','Position',[0.07 y 0.86 h],...
                'FontSize',10,'Callback',@(s,~)obj.onSelectOutputDevice(s));
            set(obj.PopOutputDev,'Value',find(strcmp(devList,curDev),1,'first'));
            
            % --- Device info fields ----------------------------------------------------------
            y = y - h - 0.02; % Section gap
            uicontrol(p,'Style','text','String','Selected device:',...
                'Units','normalized','Position',[0.07 y 0.35 h],...
                'HorizontalAlignment','left');
            obj.TxtDevName = uicontrol(p,'Style','text','String','-',...
                'Units','normalized','Position',[0.42 y 0.51 h],...
                'HorizontalAlignment','left','FontWeight','bold');
            
            y = y - h + 0.01; % Tighter gap for info
            uicontrol(p,'Style','text','String','Max output channels:',...
                'Units','normalized','Position',[0.07 y 0.35 h],...
                'HorizontalAlignment','left');
            obj.TxtMaxCh = uicontrol(p,'Style','text','String','-',...
                'Units','normalized','Position',[0.42 y 0.51 h],...
                'HorizontalAlignment','left');
            
            y = y - h + 0.01;
            uicontrol(p,'Style','text','String','Sample rate (Hz):',...
                'Units','normalized','Position',[0.07 y 0.35 h],...
                'HorizontalAlignment','left');
            obj.TxtFs = uicontrol(p,'Style','text','String','-',...
                'Units','normalized','Position',[0.42 y 0.51 h],...
                'HorizontalAlignment','left');
            
            y = y - h + 0.01;
            uicontrol(p,'Style','text','String','Buffer size (frames):',...
                'Units','normalized','Position',[0.07 y 0.35 h],...
                'HorizontalAlignment','left');
            obj.TxtBuf = uicontrol(p,'Style','text','String','-',...
                'Units','normalized','Position',[0.42 y 0.51 h],...
                'HorizontalAlignment','left');
            
            % --- Block size slider -----------------------------------------------------------
            y = y - h - 0.02; % Section gap
            uicontrol(p,'Style','text','String','Convolution block size:',...
                'Units','normalized','Position',[0.07 y 0.50 h],...
                'HorizontalAlignment','left');
            
            y = y - h;
            obj.SldBlock = uicontrol(p,'Style','slider','Min',128,'Max',8192,...
                'Value',obj.C.Setup.Block_size,...
                'SliderStep',[128/(8192-128) , 128/(8192-128)],...
                'Units','normalized','Position',[0.07 y 0.65 h],...
                'Callback',@(s,~)obj.onBlockSize(s));
            obj.TxtBlock = uicontrol(p,'Style','text',...
                'String',num2str(obj.C.Setup.Block_size),...
                'Units','normalized','Position',[0.75 y 0.20 h],...
                'HorizontalAlignment','left','FontWeight','bold');

            % --- NEW: Decorrelation slider ---------------------------------------------------
            % Check for C.Setup.Decorrelation, default to 1 if not present
            if isfield(obj.C.Setup, 'Decorrelation')
                initDecor = obj.C.Setup.Decorrelation;
            else
                initDecor = 1;
                obj.C.Setup.Decorrelation = 1; % Save default to setup
            end
            
            y = y - h - 0.02; % Section gap
            uicontrol(p,'Style','text','String','Decorrelation:',...
                'Units','normalized','Position',[0.07 y 0.50 h],...
                'HorizontalAlignment','left');
            
            y = y - h;
            obj.SldDecor = uicontrol(p,'Style','slider','Min',-1,'Max',1,...
                'Value',initDecor,...
                'SliderStep',[0.01 / (1 - (-1)), 0.1 / (1 - (-1))],... % 0.01 minor, 0.1 major
                'Units','normalized','Position',[0.07 y 0.65 h],...
                'Callback',@(s,~)obj.onDecorrelation(s));
            obj.TxtDecor = uicontrol(p,'Style','text',...
                'String',sprintf('%.2f', initDecor),...
                'Units','normalized','Position',[0.75 y 0.20 h],...
                'HorizontalAlignment','left','FontWeight','bold');

            % ---------------------------------------------------------------------------------
            
            % Initial info fill
            obj.updateDeviceInfo();
            
            % Master EQ panel (bottom)
            obj.C.ensureEngine();
            obj.EQPanel = audioapp.ui.controls.MasterEQPanel(obj.Panel, obj.C);
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
        
        % --- NEW CALLBACK ---
        function onDecorrelation(obj, src)
            val = get(src, 'Value');
            obj.C.Setup.Decorrelation = val;
            obj.C.Scene.scene_renderer.Setup.Decorrelation = val;
            obj.C.Scene.scene_renderer.decorrelator.alpha = val;
            obj.C.Scene.scene_renderer.decorrelator.update_icc_filters; 
            obj.TxtDecor.String       = sprintf('%.2f', val);
            
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