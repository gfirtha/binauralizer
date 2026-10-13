    classdef SSDTab < handle
        properties
            C
            Panel
            PopGeom  matlab.ui.control.UIControl
            EdtSize  matlab.ui.control.UIControl
            EdtNum   matlab.ui.control.UIControl
            PopLSModel matlab.ui.control.UIControl
            EdtLSRadius matlab.ui.control.UIControl
            BtnLoadGeom matlab.ui.control.UIControl   % NEW
    
            SetupLsnr event.listener   % <--- NEW
            SpeakerLsnr event.listener
            LstElem  matlab.ui.control.UIControl      % SSD loudspeakers (multi-select)
            PopSpk   matlab.ui.control.UIControl      % designed loudspeaker models
            TxtSpk   matlab.ui.control.UIControl
            TxtEq    matlab.ui.control.UIControl
        end

        properties (Constant)
            LoadEntry = 'Load from file...'
        end
    
        methods
            function obj = SSDTab(parentTab, ctx)
                obj.C = ctx;
                obj.Panel = parentTab;
                p = obj.Panel;
                obj.SetupLsnr = addlistener(obj.C, 'SetupChanged', @(~,~)obj.refreshFromSetup());
    
                uicontrol(p,'Style','text','String','Geometry:',...
                    'Units','normalized','Position',[0.07 0.92 0.30 0.06],...
                    'HorizontalAlignment','left');
    
                % add "General (from file)"
                obj.PopGeom = uicontrol(p,'Style','popupmenu',...
                    'String',{'Linear','Circular','Stereo','5.0 (circular)','Atmos 5.1.4','Atmos 7.1.4','Soundbar (16 ch)','Genius Array','General'},...
                    'Value',obj.geomToVal(obj.C.Setup.Loudspeaker_setup.Shape),...
                    'Units','normalized','Position',[0.40 0.92 0.50 0.06],...
                    'Callback',@(s,~)obj.onGeometry(s),'UserData',1);
    
                uicontrol(p,'Style','text','String','Size (R):',...
                    'Units','normalized','Position',[0.07 0.85 0.30 0.06],...
                    'HorizontalAlignment','left');
    
                obj.EdtSize = uicontrol(p,'Style','edit',...
                    'String',num2str(obj.C.Setup.Loudspeaker_setup.R),...
                    'Units','normalized','Position',[0.40 0.85 0.50 0.06],...
                    'Callback',@(s,~)obj.onSize(s));
    
                uicontrol(p,'Style','text','String','LS number (N):',...
                    'Units','normalized','Position',[0.07 0.78 0.30 0.06],...
                    'HorizontalAlignment','left');
    
                obj.EdtNum = uicontrol(p,'Style','edit',...
                    'String',num2str(obj.C.Setup.Loudspeaker_setup.N),...
                    'Units','normalized','Position',[0.40 0.78 0.50 0.06],...
                    'Callback',@(s,~)obj.onNum(s));
    
                uicontrol(p,'Style','text','String','Loudspeaker model:',...
                    'Units','normalized','Position',[0.07 0.71 0.30 0.06],...
                    'HorizontalAlignment','left');
    
                obj.PopLSModel = uicontrol(p,'Style','popupmenu',...
                    'String',{'Point source','Piston'},...
                    'Value', 1 + strcmpi(obj.C.Setup.Loudspeaker_type.Shape,'circular_piston'),...
                    'Units','normalized','Position',[0.40 0.71 0.50 0.06],...
                    'Callback',@(s,~)obj.onLSModel(s));
    
                uicontrol(p,'Style','text','String','LS radius (m):',...
                    'Units','normalized','Position',[0.07 0.64 0.30 0.06],...
                    'HorizontalAlignment','left');
    
                obj.EdtLSRadius = uicontrol(p,'Style','edit',...
                    'String',num2str(obj.C.Setup.Loudspeaker_type.R),...
                    'Units','normalized','Position',[0.40 0.64 0.50 0.06],...
                    'Callback',@(s,~)obj.onLSRadius(s));
    
                % NEW: Load geometry button (only visible for General)
                obj.BtnLoadGeom = uicontrol(p,'Style','pushbutton','String','Load geometry ...',...
                    'Units','normalized','Position',[0.40 0.57 0.50 0.06],...
                    'Visible','off', ...
                    'Callback',@(s,~)obj.onLoadGeometry());
    
                % Designed loudspeakers (Loudspeaker designer) on SSD elements
                uicontrol(p,'Style','text','String','Designed loudspeakers on the SSD elements', ...
                    'FontWeight','bold','Units','normalized','Position',[0.07 0.49 0.86 0.05], ...
                    'HorizontalAlignment','left');
                obj.LstElem = uicontrol(p,'Style','listbox','Max',2,'Min',0, ...
                    'Units','normalized','Position',[0.07 0.03 0.48 0.45], ...
                    'TooltipString','SSD loudspeakers (Ctrl / Shift: several)');
                obj.PopSpk = uicontrol(p,'Style','popupmenu','String',{obj.LoadEntry}, ...
                    'Units','normalized','Position',[0.58 0.41 0.35 0.06]);
                uicontrol(p,'Style','pushbutton','String','Assign', ...
                    'Units','normalized','Position',[0.58 0.33 0.35 0.06], ...
                    'TooltipString',['One-driver model: one on every selected loudspeaker. ' ...
                    'N-driver model: N selected loudspeakers become its drivers (in list order).'], ...
                    'Callback',@(~,~)obj.onAssign());
                uicontrol(p,'Style','pushbutton','String','Remove', ...
                    'Units','normalized','Position',[0.58 0.26 0.35 0.06], ...
                    'TooltipString','Plain loudspeakers again (the whole model of every selected one)', ...
                    'Callback',@(~,~)obj.onRemoveSpeaker());
                uicontrol(p,'Style','pushbutton','String','Speaker EQ...', ...
                    'Units','normalized','Position',[0.58 0.19 0.22 0.06], ...
                    'TooltipString',['The loudspeakers'' own EQ on their driving signals (all renderers; sound card ' ...
                    'and binaural): a Loudspeaker designer Pre-EQ export (biquads / FIR) or a hinv.mat (FIR taps + fs); ' ...
                    'one filter for all loudspeakers or one column each (SSD order)'], ...
                    'Callback',@(~,~)obj.onSpeakerEQ());
                uicontrol(p,'Style','pushbutton','String','Off', ...
                    'Units','normalized','Position',[0.81 0.19 0.12 0.06], ...
                    'Callback',@(~,~)obj.onSpeakerEQOff());
                obj.TxtEq = uicontrol(p,'Style','text','String','Speaker EQ: off', ...
                    'Units','normalized','Position',[0.58 0.13 0.35 0.05], ...
                    'HorizontalAlignment','left');
                obj.TxtSpk = uicontrol(p,'Style','text','String','', ...
                    'Units','normalized','Position',[0.58 0.03 0.35 0.10], ...
                    'HorizontalAlignment','left','ForegroundColor',[0.1 0.3 0.6]);
                obj.refreshEq();
                obj.SpeakerLsnr = addlistener(obj.C, 'SpeakersChanged', @(~,~)obj.refreshSpeakers());
                obj.refreshSpeakers();

                % ensure proper visibility at startup
                obj.updateGeneralButtonVisibility();
            end

            function refreshSpeakers(obj)
                % element list (with the model each one is a driver of) and
                % the available models
                A = {};
                try, A = obj.C.Scene.loudspeaker_array; catch, end
                s = cell(numel(A), 1);
                for e = 1:numel(A)
                    l = A{e}; s{e} = sprintf('%d: ch %d  (%.2f, %.2f)', e, l.channelIndex, l.position(1), l.position(2));
                    P = obj.C.placementOfElement(e);
                    if ~isempty(P)
                        s{e} = sprintf('%s  - %s', s{e}, P.Model.Name);
                        if P.Model.nChannels > 1, s{e} = sprintf('%s [%d]', s{e}, find(P.Elements == e)); end
                    end
                end
                v = get(obj.LstElem, 'Value'); v = v(v <= numel(s));
                set(obj.LstElem, 'String', s, 'Value', v);
                names = obj.C.speakerNames();
                labels = cellfun(@(n) sprintf('%s (%d ch)', n, obj.C.Speakers(n).nChannels), names, 'UniformOutput', false);
                cur = get(obj.PopSpk, 'Value');
                set(obj.PopSpk, 'String', [reshape(labels, 1, []), {obj.LoadEntry}], 'UserData', names, ...
                    'Value', min(cur, numel(names) + 1));
            end

            function onAssign(obj)
                sel = get(obj.LstElem, 'Value');
                if isempty(sel), obj.status('Select SSD loudspeakers first.'); return; end
                names = get(obj.PopSpk, 'UserData'); v = get(obj.PopSpk, 'Value');
                if v > numel(names)
                    [f, pth] = uigetfile('*.mat', 'Designed loudspeaker (Loudspeaker designer file)', ...
                        [audioapp.AppContext.speakerDir() filesep]);
                    if isequal(f, 0), return; end
                    try
                        S = obj.C.loadSpeakerFile(fullfile(pth, f));
                    catch ME
                        errordlg(ME.message, 'Designed loudspeaker'); return;
                    end
                else
                    S = obj.C.Speakers(names{v});
                end
                obj.assign(S, sel);
            end

            function assign(obj, S, sel)
                % one-driver model: one on every loudspeaker in sel; an
                % n-driver model: exactly n loudspeakers, driver k = sel(k)
                n = S.nChannels;
                try
                    if n == 1
                        for e = sel(:)', obj.C.placeSpeaker(S, e); end
                        obj.status(sprintf('"%s" on %d loudspeaker(s).', S.Name, numel(sel)));
                    elseif numel(sel) == n
                        obj.C.placeSpeaker(S, sel);
                        obj.status(sprintf('"%s": drivers 1..%d = loudspeakers %s.', S.Name, n, mat2str(sel)));
                    else
                        obj.status(sprintf('"%s" has %d drivers: select %d loudspeakers (%d selected).', ...
                            S.Name, n, n, numel(sel)));
                    end
                catch ME
                    obj.status(['Failed: ' ME.message]);
                end
                obj.refreshSpeakers();
            end

            function onRemoveSpeaker(obj)
                for e = get(obj.LstElem, 'Value')
                    P = obj.C.placementOfElement(e);
                    if ~isempty(P), obj.C.removeSpeaker(P); end
                end
                obj.refreshSpeakers();
            end

            function status(obj, s)
                set(obj.TxtSpk, 'String', s);
            end

            function onSpeakerEQ(obj)
                [f, pth] = uigetfile('*.mat', 'Speaker EQ (Loudspeaker designer Pre-EQ export or hinv.mat)', ...
                    [audioapp.AppContext.speakerDir() filesep]);
                if isequal(f, 0), return; end
                try
                    obj.C.setSpeakerEQ(speaker_eq.read(fullfile(pth, f)));
                catch ME
                    errordlg(ME.message, 'Speaker EQ'); return;
                end
                obj.refreshEq();
            end

            function onSpeakerEQOff(obj)
                obj.C.setSpeakerEQ([]); obj.refreshEq();
            end

            function refreshEq(obj)
                E = [];
                try, E = obj.C.Scene.scene_renderer.SpeakerEQ; catch, end
                if isempty(E), set(obj.TxtEq, 'String', 'Speaker EQ: off'); return; end
                if strcmp(E.kind, 'sos'), k = sprintf('%d biquads', size(E.sos{1}, 1)); else, k = sprintf('FIR %d taps', size(E.fir, 1)); end
                set(obj.TxtEq, 'String', sprintf('Speaker EQ: %s (%s, %d ls)', E.name, k, E.Nls));
            end
    
            function refreshFromSetup(obj)
                S = obj.C.Setup;
    
                % geometry popup
                try
                    set(obj.PopGeom,'Value', obj.geomToVal(S.Loudspeaker_setup.Shape));
                catch
                    % if shape not in list (e.g., 'general'), clamp to 'General'
                    set(obj.PopGeom,'Value', 9);
                end
    
                % numeric edits
                set(obj.EdtSize,    'String', num2str(S.Loudspeaker_setup.R));
                set(obj.EdtNum,     'String', num2str(S.Loudspeaker_setup.N));
                set(obj.PopLSModel, 'Value', 1 + strcmpi(S.Loudspeaker_type.Shape,'circular_piston'));
                set(obj.EdtLSRadius,'String', num2str(S.Loudspeaker_type.R));
    
                % show/hide the “Load geometry …” button
                obj.updateGeneralButtonVisibility();
                obj.refreshSpeakers();
            end
    
    
            function onGeometry(obj, src)
                list  = get(src,'String'); val = get(src,'Value');
                chosenUI = list{val};
                obj.C.setGeometry(chosenUI, str2num(obj.EdtNum.String));
            end
    
    
            function onSize(obj, src)
                v = str2double(get(src,'String'));
                if ~isnan(v) && v>0
                    obj.C.Setup.Loudspeaker_setup.R = v;
                    obj.rebuildScenePreservingAxes();
                else
                    set(src,'String',num2str(obj.C.Setup.Loudspeaker_setup.R));
                end
            end
    
            function onNum(obj, src)
                v = str2double(get(src,'String'));
                if ~isnan(v) && v>2 && mod(v,1)==0
                    obj.C.Setup.Loudspeaker_setup.N = v;
                    % If already in general mode and a custom geometry is loaded,
                    % you may want to re-sample here—left to your workflow.
                    if obj.isGeneralSelected()
                        % optional: prompt to re-load/re-sample
                    end
                    obj.rebuildScenePreservingAxes();
                else
                    set(src,'String',num2str(obj.C.Setup.Loudspeaker_setup.N));
                end
            end
    
            function onLSModel(obj, src)
                if get(src,'Value') == 1
                    obj.C.Setup.Loudspeaker_type.Shape = 'point_source';
                else
                    obj.C.Setup.Loudspeaker_type.Shape = 'circular_piston';
                end
                obj.rebuildScenePreservingAxes();
            end
    
            function onLSRadius(obj, src)
                v = str2double(get(src,'String'));
                if ~isnan(v) && isfinite(v) && v>0
                    obj.C.Setup.Loudspeaker_type.R = v;
                    obj.rebuildScenePreservingAxes();
                else
                    set(src,'String',num2str(obj.C.Setup.Loudspeaker_type.R));
                end
            end
    
            % --- helpers ---
            function rebuildScenePreservingAxes(obj)
                obj.C.rebuildSoundScene();
                try
                    obj.C.SceneGUI.draw_renderer_properties(obj.C.Scene.scene_renderer.SFS_renderer{1});
                catch
                end
            end
    
            function updateGeneralButtonVisibility(obj)
                set(obj.BtnLoadGeom,'Visible', ternary(obj.isGeneralSelected(),'on','off'));
            end
    
            function tf = isGeneralSelected(obj)
                list = get(obj.PopGeom,'String'); val = get(obj.PopGeom,'Value');
                tf = strcmpi(list{val},'General') || strcmpi(obj.C.Setup.Loudspeaker_setup.Shape,'general');
            end
    
            function v = geomToVal(~,g)
                switch lower(g)
                    case 'linear',        v=1;
                    case 'circular',      v=2;
                    case 'stereo',        v=3;
                    case '5.0',           v=4;
                    case '5.1.4',         v=5;
                    case '7.1.4',         v=6;
                    case 'soundbar',      v=7;
                    case 'genius_array',  v=8;
                    case 'general',       v=9;
                    otherwise,            v=1;
                end
            end
        end
    
        
    end
    
    function y = ternary(c,a,b), if c, y=a; else, y=b; end, end
    

