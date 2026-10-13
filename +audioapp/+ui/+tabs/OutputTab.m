classdef OutputTab < handle
    % Output settings: the "Binaural & headphones" tab (binauralization,
    % HRTF set, headphone EQ, level: HeadphoneLevelPanel) and "Device & EQ"
    % (device picker, device info, block size slider, master EQ). Uses AudioEngine as the single owner
    % of the audioDeviceWriter.
    properties
        C
        Panel
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
        PopDriver   matlab.ui.control.UIControl
        HpPanel     audioapp.ui.controls.HeadphoneLevelPanel   % "Headphones & level" sub-tab
    end

    methods
        function obj = OutputTab(parentTab, ctx)
            obj.C = ctx;
            % two sub-tabs: the device / EQ settings and the headphones,
            % ear (HRTF) and calibrated playback level
            tg = uitabgroup(parentTab, 'Units', 'normalized', 'Position', [0 0 1 1]);
            tHp = uitab(tg, 'Title', 'Binaural & headphones');      % binauralization, HRTF, headphone EQ, level
            tDev = uitab(tg, 'Title', 'Device & EQ');
            obj.Panel = tDev;
            % Make sure the engine (and thus the output device) exists
            obj.C.ensureEngine();
            obj.C.Engine.refreshOutput();
            p = obj.Panel;

            % Define a standard height and gap for a denser layout
            h = 0.05;    % Standard height for most controls
            y = 0.98;    % Current vertical top position

            % (binauralization and the HRTF set: the "Binaural & headphones" tab)

            % --- Device picker ---------------------------------------------------------------
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

            % Alapértelmezett driver, ha még nincs
            if ~isfield(obj.C.Setup, 'Driver'), obj.C.Setup.Driver = 'DirectSound'; end

            % --- Driver picker (ÚJ SZEKCIÓ) ---
            y = y - h - 0.02;
            uicontrol(p,'Style','text','String','Select Host API (Driver):',...
                'Units','normalized','Position',[0.07 y 0.86 h],...
                'HorizontalAlignment','left');
            y = y - h;
            driverOptions = {'DirectSound', 'ASIO', 'WASAPI'};
            obj.PopDriver = uicontrol(p,'Style','popupmenu','String',driverOptions,...
                'Units','normalized','Position',[0.07 y 0.86 h],...
                'Value', find(strcmp(driverOptions, obj.C.Setup.Driver)), ...
                'Callback', @(s,~)obj.onSelectDriver(s));

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
            obj.HpPanel = audioapp.ui.controls.HeadphoneLevelPanel(tHp, obj.C);
            obj.HpPanel.SwitchTo = @(n) obj.switchToDevice(n);
            obj.HpPanel.SetBinaural = @(on) obj.setBinauralization(on);
            obj.detectHeadphone();
        end

        function ok = switchToDevice(obj, name)
            % the output to the device name (the list refreshed: it may have
            % been connected just now); false when it is not an output
            ok = false;
            if ~obj.C.HasAudio, return; end
            L = obj.C.Engine.getOutputDevices(); cur = obj.C.Engine.getCurrentDeviceName();
            if ~any(strcmp(L, cur)), L = [{cur} L]; end
            set(obj.PopOutputDev, 'String', L);
            k = find(strcmp(L, name), 1);
            if isempty(k), obj.updateDeviceInfo(); return; end
            set(obj.PopOutputDev, 'Value', k); obj.onSelectOutputDevice(obj.PopOutputDev);
            ok = true;
        end

        function detectHeadphone(obj)
            % a recognised headphone (Bluetooth / USB, AutoEq model name):
            % binaural output on the system default -> the headphone's own
            % output device is selected; its measurement + EQ (HpPanel)
            if ~obj.C.HasAudio, return; end
            cur = obj.C.Engine.getCurrentDeviceName();
            if obj.HpPanel.autoDetect(cur), return; end
            if ~obj.C.Setup.Binauralization || ~any(strcmpi(cur, {'Default', '(system default)'})), return; end
            try
                T = headphone_db.index();
            catch
                return;
            end
            L = get(obj.PopOutputDev, 'String');
            for k = 1:numel(L)
                if ~isempty(headphone_db.match(T, L{k}))
                    set(obj.PopOutputDev, 'Value', k); obj.onSelectOutputDevice(obj.PopOutputDev);   % -> autoDetect
                    return;
                end
            end
        end

        % =================================================================
        % Callbacks
        % =================================================================

        function on = setBinauralization(obj, wantOn)
            % on (binaural) / off (multichannel): off needs the scene's
            % loudspeaker channels on the device; returns the state set
            on = wantOn;
            if ~wantOn
                req   = obj.requiredSpeakerChannels();
                info  = obj.C.Engine.getCurrentDeviceInfo();
                avail = info.MaxOutputChannels;
                if isnan(avail) || avail < req
                    on = true;
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

        function onSelectDriver(obj, src)
            list = get(src,'String');
            val = get(src,'Value');
            obj.C.Setup.Driver = list{val};

            % 1. Frissítjük az Engine-t az új driverrel
            obj.C.Engine.ensureOutput();

            % 2. Frissítjük a Device listát (Popup) az új driver eszközeivel
            newDevs = obj.C.Engine.getOutputDevices();
            set(obj.PopOutputDev, 'String', newDevs, 'Value', 1);

            % 3. Az első eszközt alapértelmezettnek állítjuk
            obj.onSelectOutputDevice(obj.PopOutputDev);

            % Ha ASIO-t választottunk, dobjuk fel a beállításokat (opcionális)
            if strcmp(obj.C.Setup.Driver, 'ASIO')
                asiosettings(obj.C.Engine.Out.Device);
            end
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
            if ~isempty(obj.HpPanel), obj.HpPanel.autoDetect(newDev); end
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
                case '5.1.4',  n = 10;    % incl. the (reserved) LFE channel
                case '7.1.4',  n = 12;
                case 'soundbar'
                    n = 16;
                    L = obj.C.Setup.Loudspeaker_setup;
                    if isfield(L, 'Surround51') && L.Surround51, n = 21; end   % + L R LFE Ls Rs
            end
        end
    end
end