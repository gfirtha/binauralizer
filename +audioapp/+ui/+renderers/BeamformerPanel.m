classdef BeamformerPanel < audioapp.ui.renderers.RendererPanel
    % BeamformerPanel  UI controls for the WFS-backed beamformer renderer.

    properties
        PopMode
        PopCTCMode
        TxtModeDesc
        BtnEarTransfer
        PopEarModel

        % WFS controls
        EdtTaper
        ChkAntiAlias
        ChkLowCut
        TxtLowCutFreq
        ChkShowWFS

        % CTC plant inversion
        ChkCTCInverse
        EdtCTCBeta

        % Symmetric beam drag
        ChkSymmetricDrag
        SymmetricDrag = false

        % Speaker EQ

        % Multiband (RACE+LSFI) crossover controls
        TxtMBLabel
        TxtMBFc
        SldMBFc
        TxtMBOptimal

        % LSFI steering controls
        TxtLSFILabel
        TxtLSFIHead
        EdtLSFIHead
        TxtLSFIAngles
        TxtLSFILineL
        TxtLSFILineR
        EdtLSFIMainL
        EdtLSFIZeroL
        EdtLSFIMainR
        EdtLSFIZeroR
        BtnLSFIAuto

        % Optimal search
        EdtFLow
        EdtFHigh
        BtnOptSearch
        PopSearchMetric

        % Scrollable viewport
        Content            % inner content uipanel (pixel units)
        ScrollSlider       % vertical slider on the right of Root
        ContentHeightPx    = 760
        ViewportHeightPx   = 500
        ResizeLsnr         % listener for Root SizeChanged

        % Beam target graphics on main axes
        hBeamL      % patch handle for L beam target
        hBeamR      % patch handle for R beam target

        % Custom drag state (replaces draggable library)
        dragSide    = ''    % 'L' | 'R' | '' (inactive)
        dragOrigin  = []    % [x0 y0] axes coords at mouse-down
        dragPosOrig = []    % beam position at mouse-down
        dragFigListeners = []  % event listener handles

        GeomLsnr    % event.listener for AppContext GeometryChanged
    end

    methods
        function obj = BeamformerPanel(C)
            obj = obj@audioapp.ui.renderers.RendererPanel(C);
        end

        function build(obj, parent)
            obj.Root = uipanel('Parent', parent, 'Units', 'normalized', ...
                'Position', [0 0 1 1], 'BorderType', 'etchedin', 'Title', 'Beamformer');

            % --- Scrollable viewport ---
            % The content panel is 40% taller than the viewport (height=1.4
            % in normalized coords of Root). A slider on the right edge
            % shifts it up/down.
            contentH = 1.65;  % normalised height (1 = viewport)
            sliderW  = 0.04;

            obj.ContentHeightPx  = contentH;  % re-using the field name for the ratio
            obj.ViewportHeightPx = 1;         % normalised

            obj.ScrollSlider = uicontrol(obj.Root, 'Style', 'slider', ...
                'Units', 'normalized', 'Position', [1-sliderW 0 sliderW 1], ...
                'Min', 0, 'Max', 1, 'Value', 1, ...
                'SliderStep', [0.05 0.2], ...
                'Callback', @(s,~)obj.onScroll(s));

            obj.Content = uipanel(obj.Root, 'Units', 'normalized', ...
                'Position', [0, 1-contentH, 1-sliderW, contentH], ...
                'BorderType', 'none');

            % All subsequent controls go into the content panel
            p = obj.Content;
            y = 0.94; h = 0.040; gap = 0.003;
            lx = 0.05; lw = 0.38; cx = 0.45; cw = 0.50;

            % --- WFS Mode selector (DAS/Focus) ---
            uicontrol(p,'Style','text','String','Beam mode:',...
                'Units','normalized','Position',[lx y lw h],'HorizontalAlignment','left');
            modes = {'DAS (Plane wave)', 'Focus (Point source)', 'Directive (Null-steering focus)', 'LSFI (Null-steering)', 'RACE (Low-freq CTC)', 'Multiband (RACE + LSFI)'};
            curMode = 'das';
            try curMode = lower(obj.C.Setup.Renderer_setup.BeamformerMode); catch, end
            if contains(curMode,'multiband'), val = 6;
            elseif contains(curMode,'race'), val = 5;
            elseif contains(curMode,'lsfi'), val = 4;
            elseif contains(curMode,'directive'), val = 3;
            elseif contains(curMode,'focus'), val = 2;
            else, val = 1; end
            obj.PopMode = uicontrol(p,'Style','popupmenu','String',modes,...
                'Value',val,'Units','normalized','Position',[cx y cw h],...
                'Callback',@(s,~)obj.onModeChange(s));

            % --- CTC Mode selector (Beam L / Beam R / CTC) ---
            y = y - h - gap;
            uicontrol(p,'Style','text','String','Output mode:',...
                'Units','normalized','Position',[lx y lw h],'HorizontalAlignment','left');
            ctcModes = {'Beam L', 'Beam R', 'CTC'};
            curCTC = 'beam_l';
            try curCTC = lower(obj.C.Setup.Renderer_setup.CTCMode); catch, end
            ctcVal = find(strcmpi(curCTC, {'beam_l','beam_r','ctc'}), 1);
            if isempty(ctcVal), ctcVal = 1; end
            obj.PopCTCMode = uicontrol(p,'Style','popupmenu','String',ctcModes,...
                'Value',ctcVal,'Units','normalized','Position',[cx y cw h],...
                'Callback',@(s,~)obj.onCTCModeChange(s));

            % --- CTC plant inversion toggle + beta ---
            y = y - h - gap;
            ctcInv = false;
            try ctcInv = obj.C.Setup.Renderer_setup.CTCInverseEnabled; catch, end
            obj.ChkCTCInverse = uicontrol(p,'Style','checkbox',...
                'String','Apply CTC plant inversion',...
                'Value',ctcInv,'Units','normalized','Position',[lx y 0.58 h],...
                'Callback',@(s,~)obj.onCTCInverseToggle(s));
            uicontrol(p,'Style','text','String','β:',...
                'Units','normalized','Position',[lx+0.60 y 0.08 h],...
                'HorizontalAlignment','right');
            betaVal = 0.05;
            try betaVal = obj.C.Setup.Renderer_setup.CTCBeta; catch, end
            obj.EdtCTCBeta = uicontrol(p,'Style','edit','String',num2str(betaVal),...
                'Units','normalized','Position',[lx+0.70 y 0.20 h],...
                'Callback',@(s,~)obj.onCTCBetaChange(s));

            % --- Symmetric drag toggle ---
            y = y - h - gap;
            obj.ChkSymmetricDrag = uicontrol(p,'Style','checkbox','String','Symmetric beam drag (mirror L↔R on x=0)',...
                'Value',0,'Units','normalized','Position',[lx y 0.90 h],...
                'Callback',@(s,~)obj.onSymmetricDragToggle(s));

            % === LSFI steering controls ===
            y = y - h - gap*3;
            obj.TxtLSFILabel = uicontrol(p,'Style','text','String','--- LSFI steering ---',...
                'Units','normalized','Position',[lx y 0.90 h],...
                'HorizontalAlignment','center','FontAngle','italic');

            y = y - h - gap;
            obj.TxtLSFIHead = uicontrol(p,'Style','text','String','Head diameter (m):',...
                'Units','normalized','Position',[lx y lw h],'HorizontalAlignment','left');
            obj.EdtLSFIHead = uicontrol(p,'Style','edit',...
                'String',obj.setupNumberString('LSFIHeadDiameter', 0.175),...
                'Units','normalized','Position',[cx y 0.25 h],...
                'Callback',@(s,~)obj.onLSFIHeadChange(s));

            y = y - h - gap;
            obj.TxtLSFIAngles = uicontrol(p,'Style','text',...
                'String','Angles are degrees relative to array broadside; use auto for receiver-derived angles.',...
                'Units','normalized','Position',[lx y 0.90 h*0.85],...
                'HorizontalAlignment','left','ForegroundColor',[0.35 0.35 0.35],...
                'FontSize',8);

            y = y - h - gap;
            obj.TxtLSFILineL = uicontrol(p,'Style','text','String','L main / zero:',...
                'Units','normalized','Position',[lx y lw h],'HorizontalAlignment','left');
            obj.EdtLSFIMainL = uicontrol(p,'Style','edit',...
                'String',obj.setupAngleString('LSFIMainAngleL'),...
                'Units','normalized','Position',[cx y 0.20 h],...
                'Callback',@(s,~)obj.onLSFIAngleChange(s,'LSFIMainAngleL'));
            obj.EdtLSFIZeroL = uicontrol(p,'Style','edit',...
                'String',obj.setupAngleString('LSFIZeroAngleL'),...
                'Units','normalized','Position',[cx+0.24 y 0.20 h],...
                'Callback',@(s,~)obj.onLSFIAngleChange(s,'LSFIZeroAngleL'));

            y = y - h - gap;
            obj.TxtLSFILineR = uicontrol(p,'Style','text','String','R main / zero:',...
                'Units','normalized','Position',[lx y lw h],'HorizontalAlignment','left');
            obj.EdtLSFIMainR = uicontrol(p,'Style','edit',...
                'String',obj.setupAngleString('LSFIMainAngleR'),...
                'Units','normalized','Position',[cx y 0.20 h],...
                'Callback',@(s,~)obj.onLSFIAngleChange(s,'LSFIMainAngleR'));
            obj.EdtLSFIZeroR = uicontrol(p,'Style','edit',...
                'String',obj.setupAngleString('LSFIZeroAngleR'),...
                'Units','normalized','Position',[cx+0.24 y 0.20 h],...
                'Callback',@(s,~)obj.onLSFIAngleChange(s,'LSFIZeroAngleR'));

            y = y - h - gap;
            obj.BtnLSFIAuto = uicontrol(p,'Style','pushbutton','String','Auto LSFI angles',...
                'Units','normalized','Position',[0.15 y 0.70 h],...
                'Callback',@(~,~)obj.onLSFIAutoAngles());

            % === WFS controls separator ===
            y = y - h - gap*3;
            uicontrol(p,'Style','text','String','--- WFS parameters ---',...
                'Units','normalized','Position',[lx y 0.90 h],...
                'HorizontalAlignment','center','FontAngle','italic');

            % --- Tapering ---
            y = y - h - gap;
            uicontrol(p,'Style','text','String','Tapering:',...
                'Units','normalized','Position',[lx y lw h],'HorizontalAlignment','left');
            taper = 0.3;
            try taper = obj.C.Setup.Renderer_setup.Tapering; catch, end
            obj.EdtTaper = uicontrol(p,'Style','edit','String',num2str(taper),...
                'Units','normalized','Position',[cx y 0.25 h],...
                'Callback',@(s,~)obj.onTaperChange(s));

            % --- AntiAliasing + Low-freq pre-eq shelving (same row) ---
            y = y - h - gap;
            aa = false;
            try aa = obj.C.Setup.Renderer_setup.Antialiasing; catch, end
            obj.ChkAntiAlias = uicontrol(p,'Style','checkbox','String','Antialiasing',...
                'Value',aa,'Units','normalized','Position',[lx y 0.42 h],...
                'Callback',@(s,~)obj.onAntiAlias(s));

            lowcut = false;
            try lowcut = obj.C.Setup.Renderer_setup.LowFreqCut; catch, end
            obj.ChkLowCut = uicontrol(p,'Style','checkbox','String','Low-freq pre-eq',...
                'Value',lowcut,'Units','normalized','Position',[lx+0.45 y 0.50 h],...
                'Callback',@(s,~)obj.onLowCut(s));

            % fc readout under the row
            y = y - h*0.9;
            obj.TxtLowCutFreq = uicontrol(p,'Style','text','String','fc: — Hz',...
                'Units','normalized','Position',[lx+0.04 y 0.90 h*0.85],...
                'HorizontalAlignment','left','ForegroundColor',[0.4 0.4 0.4],'FontSize',8);

            % --- Multiband (RACE+LSFI) crossover ---
            % Only meaningful in 'multiband' mode; visibility is toggled by
            % refreshCrossoverVisibility() based on the current beam mode.
            y = y - h - gap*2;
            obj.TxtMBLabel = uicontrol(p,'Style','text',...
                'String','--- Multiband crossover ---',...
                'Units','normalized','Position',[lx y 0.90 h],...
                'HorizontalAlignment','center','FontAngle','italic');

            y = y - h - gap;
            obj.TxtMBFc = uicontrol(p,'Style','text',...
                'String','Crossover fc: 1000 Hz',...
                'Units','normalized','Position',[lx y 0.90 h],...
                'HorizontalAlignment','left');

            y = y - h - gap;
            % Slider range 200..4000 Hz (log mapping handled in callback)
            fcVal = 1000;
            try fcVal = obj.C.Setup.Renderer_setup.MBCrossoverHz; catch, end
            obj.SldMBFc = uicontrol(p,'Style','slider',...
                'Min', log10(200), 'Max', log10(4000),...
                'Value', log10(max(200, min(4000, fcVal))),...
                'Units','normalized','Position',[lx y 0.90 h],...
                'Callback', @(s,~)obj.onMBCrossoverChange(s));

            % Sync the label text with the loaded slider value
            set(obj.TxtMBFc, 'String', sprintf('Crossover fc: %g Hz', round(fcVal/10)*10));

            y = y - h - gap;
            obj.TxtMBOptimal = uicontrol(p,'Style','text',...
                'String','(optimal: — Hz)',...
                'Units','normalized','Position',[lx y 0.90 h*0.85],...
                'HorizontalAlignment','left','ForegroundColor',[0.4 0.4 0.4],...
                'FontSize',8);

            % === Visualization ===
            y = y - h - gap*2;
            obj.ChkShowWFS = uicontrol(p,'Style','checkbox','String','Show WFS visuals',...
                'Value',1,'Units','normalized','Position',[lx y 0.90 h],...
                'Callback',@(s,~)obj.onWFSVisToggle(s));

            % === Ear Transfer ===
            y = y - h - gap*3;
            uicontrol(p,'Style','text','String','Ear transfer model:',...
                'Units','normalized','Position',[lx y lw h],'HorizontalAlignment','left');
            obj.PopEarModel = uicontrol(p,'Style','popupmenu',...
                'String',{'Analytic (rigid sphere)','HRTF (measured)'},...
                'Value',2,'Units','normalized','Position',[cx y cw h],...
                'Callback',@(s,~)obj.onEarModelChange(s));

            y = y - h - gap;
            obj.BtnEarTransfer = uicontrol(p,'Style','pushbutton',...
                'String','Ear Transfer','FontWeight','bold',...
                'Units','normalized','Position',[0.15 y 0.70 h*1.3],...
                'Callback',@(~,~)obj.onEarTransfer());

            % === Optimal Position Search ===
            y = y - h - gap*3;
            uicontrol(p,'Style','text','String','--- Optimal VS Search ---',...
                'Units','normalized','Position',[lx y 0.90 h],...
                'HorizontalAlignment','center','FontAngle','italic');

            y = y - h - gap;
            uicontrol(p,'Style','text','String','Freq band (Hz):',...
                'Units','normalized','Position',[lx y 0.30 h],'HorizontalAlignment','left');
            obj.EdtFLow = uicontrol(p,'Style','edit','String','500',...
                'Units','normalized','Position',[0.36 y 0.14 h]);
            uicontrol(p,'Style','text','String','-',...
                'Units','normalized','Position',[0.51 y 0.04 h],'HorizontalAlignment','center');
            obj.EdtFHigh = uicontrol(p,'Style','edit','String','5000',...
                'Units','normalized','Position',[0.56 y 0.14 h]);

            y = y - h - gap;
            uicontrol(p,'Style','text','String','Metric:',...
                'Units','normalized','Position',[lx y 0.15 h],'HorizontalAlignment','left');
            obj.PopSearchMetric = uicontrol(p,'Style','popupmenu',...
                'String',{'Worst-case ICLD (passive suppression)','Plant condition number (CTC quality)'},...
                'Value',1,'Units','normalized','Position',[lx+0.16 y 0.78 h]);

            y = y - h*1.3 - gap;
            obj.BtnOptSearch = uicontrol(p,'Style','pushbutton',...
                'String','Get Optimal Position','FontWeight','bold',...
                'Units','normalized','Position',[0.10 y 0.80 h*1.4],...
                'BackgroundColor',[0.85 0.93 1.0],...
                'Callback',@(~,~)obj.onOptimalSearch());

            % Mode description at bottom
            y = y - h*1.5 - gap;
            obj.TxtModeDesc = uicontrol(p,'Style','text','String','',...
                'Units','normalized','Position',[lx y 0.90 h*1.8],...
                'HorizontalAlignment','left','FontSize',8,'Max',2);

            obj.refreshAngle();
            obj.refreshModeDesc();
            obj.refreshCrossoverVisibility();
            obj.refreshLSFIControls();

            % Sync viz checkbox with overlay settings; ensure all fields present
            try
                s = obj.C.SceneGUI;
                if ~isfield(s.rendererSetting, 'beamformer_renderer')
                    s.rendererSetting.beamformer_renderer = struct();
                end
                ss = s.rendererSetting.beamformer_renderer;
                % Populate any missing fields with defaults
                if ~isfield(ss,'showRef'),              ss.showRef              = true;  end
                if ~isfield(ss,'showAmp'),              ss.showAmp              = false; end
                if ~isfield(ss,'showTaper'),            ss.showTaper            = false; end
                if ~isfield(ss,'showZones'),            ss.showZones            = true;  end
                if ~isfield(ss,'showBeamDirection'),    ss.showBeamDirection    = true;  end
                if ~isfield(ss,'showElementHighlight'), ss.showElementHighlight = true;  end
                s.rendererSetting.beamformer_renderer = ss;
                set(obj.ChkShowWFS, 'Value', ss.showRef);
            catch
            end

            % Place draggable beam target points on main axes
            obj.placeBeamTargets();

            % Listen for geometry changes (VS/receiver drag) to refresh fc display
            obj.GeomLsnr = addlistener(obj.C, 'GeometryChanged', ...
                @(~,~)obj.refreshLowCutFreq());

            % Initial cutoff display
            obj.refreshLowCutFreq();
        end

        function delete(obj)
            % Clean up beam target graphics, drag listeners, event listener
            obj.removeDragListeners();
            try if ~isempty(obj.GeomLsnr) && isvalid(obj.GeomLsnr), delete(obj.GeomLsnr); end; catch, end
            try if ~isempty(obj.hBeamL) && isgraphics(obj.hBeamL), delete(obj.hBeamL); end; catch, end
            try if ~isempty(obj.hBeamR) && isgraphics(obj.hBeamR), delete(obj.hBeamR); end; catch, end
            try
                ax = obj.C.SceneGUI.main_axes;
                delete(findobj(ax, 'Tag', 'beamL_txt'));
                delete(findobj(ax, 'Tag', 'beamR_txt'));
            catch, end
        end

        function syncFromSetup(obj)
            obj.refreshAngle();
            % Retry placing beam targets if they weren't placed during build
            if isempty(obj.hBeamL) || ~isgraphics(obj.hBeamL)
                obj.placeBeamTargets();
            end
            % Refresh low-cut state
            try
                lowcut = false;
                if isfield(obj.C.Setup.Renderer_setup, 'LowFreqCut')
                    lowcut = logical(obj.C.Setup.Renderer_setup.LowFreqCut);
                end
                obj.applyLowCutToBeams(lowcut);
            catch
            end
            obj.refreshLowCutFreq();

            % CTC plant inversion state → renderer
            try
                r = obj.getRenderer();
                if ~isempty(r)
                    if isfield(obj.C.Setup.Renderer_setup,'CTCInverseEnabled')
                        r.CTCInverseEnabled = logical(obj.C.Setup.Renderer_setup.CTCInverseEnabled);
                    end
                    if isfield(obj.C.Setup.Renderer_setup,'CTCBeta')
                        r.ctc_beta = obj.C.Setup.Renderer_setup.CTCBeta;
                    end
                    if isfield(obj.C.Setup.Renderer_setup,'EarTransferModel')
                        r.ear_transfer_model = obj.C.Setup.Renderer_setup.EarTransferModel;
                    end
                    r.ctc_inv_filters = [];  % force rebuild
                end
            catch
            end
        end

        function s = setupNumberString(obj, field, defaultValue)
            v = defaultValue;
            try
                if isfield(obj.C.Setup.Renderer_setup, field)
                    v = obj.C.Setup.Renderer_setup.(field);
                end
            catch
            end
            s = sprintf('%.4g', v);
        end

        function s = setupAngleString(obj, field)
            v = NaN;
            try
                if isfield(obj.C.Setup.Renderer_setup, field)
                    v = obj.C.Setup.Renderer_setup.(field);
                end
            catch
            end
            if isnumeric(v) && isscalar(v) && isfinite(v)
                s = sprintf('%.2f', v);
            else
                s = 'auto';
            end
        end

        function onLSFIHeadChange(obj, src)
            v = str2double(get(src, 'String'));
            if ~isnan(v) && isfinite(v) && v > 0 && v <= 1
                obj.C.Setup.Renderer_setup.LSFIHeadDiameter = v;
                set(src, 'String', obj.setupNumberString('LSFIHeadDiameter', 0.175));
                obj.updateRenderer();
            else
                set(src, 'String', obj.setupNumberString('LSFIHeadDiameter', 0.175));
            end
        end

        function onLSFIAngleChange(obj, src, field)
            txt = strtrim(get(src, 'String'));
            if isempty(txt) || strcmpi(txt, 'auto') || strcmpi(txt, 'nan')
                v = NaN;
            else
                v = str2double(txt);
                if isnan(v) || ~isfinite(v) || abs(v) > 89
                    set(src, 'String', obj.setupAngleString(field));
                    return;
                end
            end

            obj.C.Setup.Renderer_setup.(field) = v;
            set(src, 'String', obj.setupAngleString(field));
            obj.updateRenderer();
        end

        function onLSFIAutoAngles(obj)
            fields = {'LSFIMainAngleL','LSFIZeroAngleL', ...
                      'LSFIMainAngleR','LSFIZeroAngleR'};
            for k = 1:numel(fields)
                obj.C.Setup.Renderer_setup.(fields{k}) = NaN;
            end
            obj.refreshLSFIControls();
            obj.updateRenderer();
        end

        function refreshLSFIControls(obj)
            try
                if ~isempty(obj.EdtLSFIHead) && isgraphics(obj.EdtLSFIHead)
                    set(obj.EdtLSFIHead, 'String', obj.setupNumberString('LSFIHeadDiameter', 0.175));
                end
                if ~isempty(obj.EdtLSFIMainL) && isgraphics(obj.EdtLSFIMainL)
                    set(obj.EdtLSFIMainL, 'String', obj.setupAngleString('LSFIMainAngleL'));
                end
                if ~isempty(obj.EdtLSFIZeroL) && isgraphics(obj.EdtLSFIZeroL)
                    set(obj.EdtLSFIZeroL, 'String', obj.setupAngleString('LSFIZeroAngleL'));
                end
                if ~isempty(obj.EdtLSFIMainR) && isgraphics(obj.EdtLSFIMainR)
                    set(obj.EdtLSFIMainR, 'String', obj.setupAngleString('LSFIMainAngleR'));
                end
                if ~isempty(obj.EdtLSFIZeroR) && isgraphics(obj.EdtLSFIZeroR)
                    set(obj.EdtLSFIZeroR, 'String', obj.setupAngleString('LSFIZeroAngleR'));
                end
            catch
            end

            try
                mode = obj.C.Setup.Renderer_setup.BeamformerMode;
            catch
                mode = 'das';
            end
            visible = 'off';
            if strcmpi(mode, 'lsfi') || strcmpi(mode, 'multiband')
                visible = 'on';
            end

            hs = {obj.TxtLSFILabel, obj.TxtLSFIHead, obj.EdtLSFIHead, ...
                  obj.TxtLSFIAngles, obj.TxtLSFILineL, obj.EdtLSFIMainL, ...
                  obj.EdtLSFIZeroL, obj.TxtLSFILineR, obj.EdtLSFIMainR, ...
                  obj.EdtLSFIZeroR, obj.BtnLSFIAuto};
            for k = 1:numel(hs)
                if ~isempty(hs{k}) && isgraphics(hs{k})
                    set(hs{k}, 'Visible', visible);
                end
            end
        end

        function onModeChange(obj, src)
            val = get(src,'Value');
            modeNames = {'das', 'focus', 'directive', 'lsfi', 'race', 'multiband'};
            obj.C.Setup.Renderer_setup.BeamformerMode = modeNames{val};

            isRaceLike = strcmpi(modeNames{val}, 'race') || strcmpi(modeNames{val}, 'multiband');

            % In RACE/Multiband mode the only meaningful Output mode is 'ctc'
            % (the 2x2 mixing matrix is intrinsically two-channel). Force it
            % and make the popup reflect that.
            if isRaceLike
                obj.C.Setup.Renderer_setup.CTCMode = 'ctc';
                if isvalid(obj.PopCTCMode)
                    set(obj.PopCTCMode, 'String', {'CTC'}, 'Value', 1, 'Enable', 'off');
                end
                try
                    r = obj.getRenderer();
                    if ~isempty(r)
                        r.ctc_mode = 'ctc';
                        if isempty(r.hrtf_extrap)
                            try
                                r.hrtf_extrap = obj.C.Scene.scene_renderer.HRTF_extrapolators;
                            catch
                            end
                        end
                    end
                catch
                end
            else
                % Restore the full set of output options for non-RACE modes
                if isvalid(obj.PopCTCMode)
                    set(obj.PopCTCMode, 'String', {'Beam L', 'Beam R', 'CTC'}, 'Enable', 'on');
                end
            end

            obj.updateRenderer();
            obj.refreshModeDesc();
            obj.refreshCrossoverVisibility();
            obj.refreshLSFIControls();
        end

        function onCTCModeChange(obj, src)
            val = get(src, 'Value');
            ctcModes = {'beam_l', 'beam_r', 'ctc'};
            obj.C.Setup.Renderer_setup.CTCMode = ctcModes{val};
            try
                r = obj.getRenderer();
                if ~isempty(r)
                    r.ctc_mode = ctcModes{val};
                    % For CTC mode, wire the HRTF extrapolator
                    if val == 3 && isempty(r.hrtf_extrap)
                        try
                            r.hrtf_extrap = obj.C.Scene.scene_renderer.HRTF_extrapolators;
                        catch
                        end
                    end
                end
            catch
            end
        end

        function onTaperChange(obj, src)
            v = str2double(get(src,'String'));
            if ~isnan(v) && v >= 0 && v <= 1
                obj.C.Setup.Renderer_setup.Tapering = v;
                obj.updateRenderer();
            else
                set(src,'String','0.3');
            end
        end

        function onAntiAlias(obj, src)
            obj.C.Setup.Renderer_setup.Antialiasing = logical(get(src,'Value'));
            obj.updateRenderer();
        end

        function onLowCut(obj, src)
            val = logical(get(src, 'Value'));
            obj.C.Setup.Renderer_setup.LowFreqCut = val;
            obj.applyLowCutToBeams(val);
            obj.tryDrawProps();
            obj.refreshLowCutFreq();
        end

        function applyLowCutToBeams(obj, enabled)
            r = obj.getRenderer();
            if isempty(r), return; end
            r.LowFreqCut = enabled;
            % The renderer computes ear receivers in update_renderer() and
            % forwards them into wfsL/wfsR. We just trigger the update.
            r.update_renderer();
        end

        function refreshLowCutFreq(obj)
            r = obj.getRenderer();
            if isempty(r), return; end
            try
                parts = {};
                for side = {'L','R'}
                    if side{1} == 'L', w = r.wfsL; else, w = r.wfsR; end
                    if isempty(w), continue; end
                    fc_low  = w.lower_cutoff_freq;
                    fc_high = w.higher_cutoff_freq;
                    if fc_low > 0 && isfinite(fc_high)
                        parts{end+1} = sprintf('%s: fc%s%.0f  fc%s%.0f Hz', ...
                            side{1}, char(8595), fc_low, char(8593), fc_high); %#ok<AGROW>
                    end
                end
                if ~isempty(parts)
                    set(obj.TxtLowCutFreq, 'String', strjoin(parts, '  '));
                else
                    set(obj.TxtLowCutFreq, 'String', 'fc: — Hz');
                end
            catch
                set(obj.TxtLowCutFreq, 'String', 'fc: — Hz');
            end
        end

        function r = getRenderer(obj)
            r = [];
            try
                rr = obj.C.Scene.scene_renderer.SFS_renderer;
                for k = 1:numel(rr)
                    if isa(rr{k}, 'beamformer_renderer')
                        r = rr{k};
                        return;
                    end
                end
            catch
            end
        end

        function placeBeamTargets(obj)
            % Draw two beam target circles (L=blue, R=red) with a custom drag
            % engine that works reliably across repeated draw calls.
            r = obj.getRenderer();
            if isempty(r), return; end

            ax = obj.C.SceneGUI.main_axes;
            if isempty(ax) || ~isgraphics(ax), return; end

            % Clean up old patches and listeners first
            try
                if ~isempty(obj.hBeamL) && isgraphics(obj.hBeamL), delete(obj.hBeamL); end
                if ~isempty(obj.hBeamR) && isgraphics(obj.hBeamR), delete(obj.hBeamR); end
                delete(findobj(ax, 'Tag', 'beamL_txt'));
                delete(findobj(ax, 'Tag', 'beamR_txt'));
            catch, end
            obj.removeDragListeners();

            % Circle shape (2D, closed)
            th  = linspace(0, 2*pi, 25)';
            rad = 0.06;
            cx  = rad * cos(th);
            cy  = rad * sin(th);

            posL = r.beamL_pos(1:2);
            posR = r.beamR_pos(1:2);

            % L patch (blue) — simple XData/YData form, no z column
            obj.hBeamL = patch(ax, cx + posL(1), cy + posL(2), [0.2 0.4 1], ...
                'EdgeColor', [0 0 0.8], 'LineWidth', 1.5, ...
                'FaceAlpha', 0.5, 'Tag', 'beamL', ...
                'HitTest', 'on', 'PickableParts', 'all');
            obj.hBeamL.UserData = struct('side','L','Origin',posL,'Orientation',0);
            setappdata(obj.hBeamL, 'beamCircleR', rad);
            text(ax, posL(1), posL(2), 'L', 'HorizontalAlignment','center', ...
                'FontWeight','bold','FontSize',10,'Color',[0 0 0.8],'Tag','beamL_txt', ...
                'HitTest','off','PickableParts','none');

            % R patch (red)
            obj.hBeamR = patch(ax, cx + posR(1), cy + posR(2), [1 0.3 0.2], ...
                'EdgeColor', [0.8 0 0], 'LineWidth', 1.5, ...
                'FaceAlpha', 0.5, 'Tag', 'beamR', ...
                'HitTest', 'on', 'PickableParts', 'all');
            obj.hBeamR.UserData = struct('side','R','Origin',posR,'Orientation',0);
            setappdata(obj.hBeamR, 'beamCircleR', rad);
            text(ax, posR(1), posR(2), 'R', 'HorizontalAlignment','center', ...
                'FontWeight','bold','FontSize',10,'Color',[0.8 0 0],'Tag','beamR_txt', ...
                'HitTest','off','PickableParts','none');

            % ButtonDownFcn on each patch
            set(obj.hBeamL, 'ButtonDownFcn', @(~,e) obj.onBeamMouseDown('L', e));
            set(obj.hBeamR, 'ButtonDownFcn', @(~,e) obj.onBeamMouseDown('R', e));

            obj.raiseBeamTargets();
        end

        % ── Custom drag engine ────────────────────────────────────────────

        function onBeamMouseDown(obj, side, evt)
            if evt.Button ~= 1, return; end
            ax  = obj.C.SceneGUI.main_axes;
            cp  = get(ax, 'CurrentPoint');
            obj.dragSide    = side;
            obj.dragOrigin  = cp(1, 1:2);
            r = obj.getRenderer();
            if isempty(r), return; end
            if side == 'L'
                obj.dragPosOrig = r.beamL_pos(1:2);
            else
                obj.dragPosOrig = r.beamR_pos(1:2);
            end
            fig = ancestor(ax, 'figure');
            % Store previous callbacks so we can restore them on release
            obj.dragFigListeners = { ...
                get(fig, 'WindowButtonMotionFcn'), ...
                get(fig, 'WindowButtonUpFcn') };
            set(fig, 'WindowButtonMotionFcn',  @(~,~) obj.onDragMotion());
            set(fig, 'WindowButtonUpFcn',      @(~,~) obj.onDragRelease());
        end

        function onDragMotion(obj)
            if isempty(obj.dragSide), return; end
            ax = obj.C.SceneGUI.main_axes;
            if isempty(ax) || ~isgraphics(ax), return; end

            cp     = get(ax, 'CurrentPoint');
            delta  = cp(1,1:2) - obj.dragOrigin;
            newPos = obj.dragPosOrig + delta;

            if obj.dragSide == 'L', h = obj.hBeamL; else, h = obj.hBeamR; end
            if isempty(h) || ~isgraphics(h), return; end

            % Move active-side circle
            obj.moveBeamHandle(h, newPos, obj.dragSide);

            % Mirror to the other side if symmetric drag is on
            if obj.SymmetricDrag
                mirrorPos = [-newPos(1), newPos(2)];
                if obj.dragSide == 'L'
                    obj.moveBeamHandle(obj.hBeamR, mirrorPos, 'R');
                else
                    obj.moveBeamHandle(obj.hBeamL, mirrorPos, 'L');
                end
            end

            % Update renderer position(s) and overlay
            r = obj.getRenderer();
            if ~isempty(r)
                if obj.dragSide == 'L'
                    r.beamL_pos = [newPos, 0];
                    if obj.SymmetricDrag
                        r.beamR_pos = [-newPos(1), newPos(2), 0];
                    end
                else
                    r.beamR_pos = [newPos, 0];
                    if obj.SymmetricDrag
                        r.beamL_pos = [-newPos(1), newPos(2), 0];
                    end
                end
                r.update_renderer();
                obj.drawActiveBeamOverlay(obj.dragSide);
                if obj.SymmetricDrag
                    % Also redraw the mirrored beam overlay
                    if obj.dragSide == 'L'
                        obj.drawActiveBeamOverlay('R');
                    else
                        obj.drawActiveBeamOverlay('L');
                    end
                end
            end
            drawnow limitrate
        end

        function moveBeamHandle(~, h, newPos, ~)
            % Move a beam circle patch + its label to newPos in axes coords.
            if isempty(h) || ~isgraphics(h), return; end
            rad = getappdata(h, 'beamCircleR');
            if isempty(rad), rad = 0.06; end
            th = linspace(0, 2*pi, 25)';
            set(h, 'XData', rad*cos(th) + newPos(1), ...
                   'YData', rad*sin(th) + newPos(2));
            h.UserData.Origin = newPos;
            % Move label
            tagCandidate = h.Tag;
            if strcmp(tagCandidate, 'beamL'),    lbl = 'beamL_txt';
            elseif strcmp(tagCandidate,'beamR'), lbl = 'beamR_txt';
            else, lbl = '';
            end
            if ~isempty(lbl)
                ax = ancestor(h, 'axes');
                if isgraphics(ax)
                    txts = findobj(ax, 'Tag', lbl);
                    if ~isempty(txts)
                        set(txts(1), 'Position', [newPos(1), newPos(2), 0]);
                    end
                end
            end
        end

        function syncBeamGraphics(obj, r)
            % Reposition beam L/R patches to match renderer beam positions.
            if nargin < 2
                r = obj.getRenderer();
            end
            if isempty(r), return; end
            obj.moveBeamHandle(obj.hBeamL, r.beamL_pos(1:2), 'L');
            obj.moveBeamHandle(obj.hBeamR, r.beamR_pos(1:2), 'R');
        end

        function onDragRelease(obj)
            if isempty(obj.dragSide), return; end

            % Restore previous figure callbacks first
            obj.removeDragListeners();

            r = obj.getRenderer();
            if ~isempty(r)
                if obj.dragSide == 'L'
                    obj.C.Setup.Renderer_setup.BeamLPos = r.beamL_pos;
                    if obj.SymmetricDrag
                        obj.C.Setup.Renderer_setup.BeamRPos = r.beamR_pos;
                    end
                else
                    obj.C.Setup.Renderer_setup.BeamRPos = r.beamR_pos;
                    if obj.SymmetricDrag
                        obj.C.Setup.Renderer_setup.BeamLPos = r.beamL_pos;
                    end
                end
            end
            obj.dragSide = '';
            obj.tryDrawProps();
            obj.refreshLowCutFreq();
        end

        function onSymmetricDragToggle(obj, src)
            obj.SymmetricDrag = logical(get(src, 'Value'));
        end

        function onCTCInverseToggle(obj, src)
            val = logical(get(src, 'Value'));
            obj.C.Setup.Renderer_setup.CTCInverseEnabled = val;
            r = obj.getRenderer();
            if ~isempty(r)
                r.CTCInverseEnabled = val;
                r.ctc_inv_filters = [];  % force rebuild
            end
        end

        function onCTCBetaChange(obj, src)
            v = str2double(get(src, 'String'));
            if ~isnan(v) && v > 0
                obj.C.Setup.Renderer_setup.CTCBeta = v;
                r = obj.getRenderer();
                if ~isempty(r)
                    r.ctc_beta = v;
                    r.ctc_inv_filters = [];  % force rebuild
                end
            else
                % revert
                try set(src, 'String', num2str(obj.C.Setup.Renderer_setup.CTCBeta)); catch, end
            end
        end

        function onEarModelChange(obj, src)
            idx = get(src, 'Value');
            if idx == 2
                model = 'hrtf';
            else
                model = 'analytic';
            end
            obj.C.Setup.Renderer_setup.EarTransferModel = model;
            r = obj.getRenderer();
            if ~isempty(r)
                r.ear_transfer_model = model;
                r.ctc_inv_filters = [];  % force rebuild
            end
        end

        function onScroll(obj, src)
            if isempty(obj.Content) || ~isgraphics(obj.Content), return; end
            v = get(src, 'Value');  % 0..1, where 1 = scrolled to top
            contentH = obj.ContentHeightPx;  % e.g. 1.4
            % Bottom of content in Root normalised coords:
            %   v=1 → see top → content bottom = 1 - contentH
            %   v=0 → see bottom → content bottom = 0
            bottomY = v * (1 - contentH);
            pos = obj.Content.Position;
            obj.Content.Position = [pos(1), bottomY, pos(3), pos(4)];
        end

        function onRootResized(~)
            % No-op in normalised mode (everything scales automatically)
        end

        function removeDragListeners(obj)
            % Restore the saved WindowButton callbacks and clear state
            try
                ax  = obj.C.SceneGUI.main_axes;
                fig = ancestor(ax, 'figure');
                if ~isempty(obj.dragFigListeners) && iscell(obj.dragFigListeners)
                    set(fig, 'WindowButtonMotionFcn', obj.dragFigListeners{1});
                    set(fig, 'WindowButtonUpFcn',     obj.dragFigListeners{2});
                end
            catch
            end
            obj.dragFigListeners = [];
        end

        function drawActiveBeamOverlay(obj, ~)
            % Full overlay redraw into the main renderer_illustration store,
            % bypassing listener_space_axes so raiseSceneObjects doesn't fire.
            r = obj.getRenderer();
            if isempty(r), return; end
            ax = obj.C.SceneGUI.main_axes;
            if isempty(ax) || ~isgraphics(ax), return; end

            s    = obj.C.SceneGUI;
            rKey = 'beamformer_renderer';
            if isfield(s.rendererSetting, rKey)
                ss = s.rendererSetting.(rKey);
            else
                ss = r.defaultOverlaySettings();
            end
            wfs_settings = struct( ...
                'showRef',   ss.showRef,   'showAmp',   ss.showAmp, ...
                'showTaper', ss.showTaper, 'showZones', ss.showZones);

            store = s.renderer_illustration;
            if ~isfield(store, 'wfsL_store'), store.wfsL_store = struct(); end
            if ~isfield(store, 'wfsR_store'), store.wfsR_store = struct(); end

            % Update both beam WFS overlays via their own sub-stores
            if ~isempty(r.wfsL)
                store.wfsL_store = audioapp.ui.draw.WFSPainter.draw( ...
                    ax, r.wfsL.get_renderer_visuals(), r.wfsL, wfs_settings, store.wfsL_store);
            end
            if ~isempty(r.wfsR)
                store.wfsR_store = audioapp.ui.draw.WFSPainter.draw( ...
                    ax, r.wfsR.get_renderer_visuals(), r.wfsR, wfs_settings, store.wfsR_store);
            end

            % Update BeamformerPainter (both beam arrows/crosses)
            rp = r.get_renderer_visuals();
            bf_settings = struct( ...
                'showBeamDirection',    ss.showBeamDirection, ...
                'showElementHighlight', ss.showElementHighlight);
            store = audioapp.ui.draw.BeamformerPainter.draw(ax, rp, bf_settings, store);

            s.renderer_illustration = store;
            obj.raiseBeamTargets();
        end

        function onWFSVisToggle(obj, src)
            val = logical(get(src,'Value'));
            try
                s = obj.C.SceneGUI;
                ss = struct('showRef', val, 'showAmp', false, ...
                            'showTaper', false, 'showZones', val, ...
                            'showBeamDirection', val, 'showElementHighlight', val);
                s.rendererSetting.beamformer_renderer = ss;
                obj.tryDrawProps();
            catch
            end
        end

        function refreshAngle(~)
            % Kept as no-op; the UI no longer displays the beam angles.
        end

        function refreshModeDesc(obj)
            try mode = obj.C.Setup.Renderer_setup.BeamformerMode; catch, mode='das'; end
            try ctc = obj.C.Setup.Renderer_setup.CTCMode; catch, ctc='beam_l'; end
            switch lower(mode)
                case 'das',       wfsDesc = 'DAS (plane-wave beams)';
                case 'focus',     wfsDesc = 'Focus (point-source beams)';
                case 'directive', wfsDesc = 'Directive (null-steering focus)';
                case 'lsfi',      wfsDesc = 'LSFI (null-steering beamformer)';
                case 'race',      wfsDesc = 'RACE (low-freq CTC, two outermost speakers)';
                case 'multiband', wfsDesc = 'Multiband: RACE (low band, 2 LS) + LSFI (high band, all 8 LS)';
                otherwise,        wfsDesc = mode;
            end
            switch lower(ctc)
                case 'beam_l', ctcDesc = 'Beam L: mono input → left beam only.';
                case 'beam_r', ctcDesc = 'Beam R: mono input → right beam only.';
                case 'ctc',    ctcDesc = 'CTC: input × HRTF → earL/earR → dual beams.';
                otherwise,     ctcDesc = '';
            end
            desc = sprintf('%s\n%s', wfsDesc, ctcDesc);
            if ~isempty(obj.TxtModeDesc) && isvalid(obj.TxtModeDesc)
                set(obj.TxtModeDesc,'String',desc);
            end
        end

        function onMBCrossoverChange(obj, src)
            % Slider value is log10(Hz); convert and round to nearest 10 Hz
            v = get(src, 'Value');
            fc = 10^v;
            fc = round(fc / 10) * 10;
            obj.C.Setup.Renderer_setup.MBCrossoverHz = fc;
            if isvalid(obj.TxtMBFc)
                set(obj.TxtMBFc, 'String', sprintf('Crossover fc: %g Hz', fc));
            end
            obj.updateRenderer();
        end

        function refreshCrossoverVisibility(obj)
            % Show the multiband crossover controls only in 'multiband' mode.
            try mode = obj.C.Setup.Renderer_setup.BeamformerMode; catch, mode='das'; end
            on = strcmpi(mode, 'multiband');
            vis = 'off'; if on, vis = 'on'; end

            if ~isempty(obj.TxtMBLabel)   && isvalid(obj.TxtMBLabel),   set(obj.TxtMBLabel,   'Visible', vis); end
            if ~isempty(obj.TxtMBFc)      && isvalid(obj.TxtMBFc),      set(obj.TxtMBFc,      'Visible', vis); end
            if ~isempty(obj.SldMBFc)      && isvalid(obj.SldMBFc),      set(obj.SldMBFc,      'Visible', vis); end
            if ~isempty(obj.TxtMBOptimal) && isvalid(obj.TxtMBOptimal), set(obj.TxtMBOptimal, 'Visible', vis); end

            if on
                % Update the optimal fc readout from the active renderer
                try
                    r = obj.getRenderer();
                    if isa(r, 'beamformer_renderer')
                        fc_opt = r.estimate_optimal_crossover();
                        set(obj.TxtMBOptimal, 'String', ...
                            sprintf('(optimal: %.0f Hz)', fc_opt));
                    end
                catch, end
            end
        end

        function onEarTransfer(obj)
            try
                r = obj.C.Scene.scene_renderer.SFS_renderer{1};
                if ~isa(r,'beamformer_renderer')
                    warndlg('Active renderer is not a beamformer.','Ear Transfer');
                    return;
                end
                if isempty(r.receiver)
                    warndlg('No receiver available.','Ear Transfer');
                    return;
                end

                modelIdx = get(obj.PopEarModel,'Value');
                if modelIdx == 1
                    model = 'analytic'; hrtf_ext = [];
                else
                    model = 'hrtf';
                    try hrtf_ext = obj.C.Scene.scene_renderer.HRTF_extrapolators;
                    catch
                        warndlg('No HRTF extrapolator available.','Ear Transfer');
                        return;
                    end
                end

                [f, H_LL, H_LR, H_RL, H_RR] = r.compute_ear_transfer(model, hrtf_ext);

                fs_half = obj.C.Setup.SampleRate / 2;
                f_lim   = [max(f(2), 20), fs_half];

                % dB magnitudes
                dB = @(H) 20*log10(abs(H) + 1e-15);
                m_LL = dB(H_LL);  m_LR = dB(H_LR);
                m_RL = dB(H_RL);  m_RR = dB(H_RR);

                % ICLD per beam  (desired / crosstalk, positive = good suppression)
                ICLD_L = m_LL - m_LR;   % L beam: left ear desired
                ICLD_R = m_RR - m_RL;   % R beam: right ear desired

                fig = figure('Name', 'Ear Transfer Functions (2×2)', ...
                    'NumberTitle', 'off', 'Position', [150 80 900 700]);

                colL = [0.15 0.35 0.85];   % blue  – L beam
                colR = [0.85 0.20 0.15];   % red   – R beam
                colX = [0.50 0.50 0.50];   % grey  – crosstalk / zero line

                % ── Row 1: L beam transfer ──
                ax1 = subplot(3,2,1, 'Parent', fig);
                semilogx(ax1, f, m_LL, 'Color', colL, 'LineWidth', 1.8, ...
                    'DisplayName', 'L beam → left ear (desired)');
                hold(ax1,'on');
                semilogx(ax1, f, m_LR, '--', 'Color', colX, 'LineWidth', 1.4, ...
                    'DisplayName', 'L beam → right ear (crosstalk)');
                grid(ax1,'on');
                xlim(ax1, f_lim); xlabel(ax1,'Frequency (Hz)'); ylabel(ax1,'Magnitude (dB)');
                title(ax1, sprintf('L beam  [%s / %s]', upper(r.wfs_mode), model));
                legend(ax1, 'show', 'Location', 'southwest', 'FontSize', 7);

                % ── Row 1: R beam transfer ──
                ax2 = subplot(3,2,2, 'Parent', fig);
                semilogx(ax2, f, m_RR, 'Color', colR, 'LineWidth', 1.8, ...
                    'DisplayName', 'R beam → right ear (desired)');
                hold(ax2,'on');
                semilogx(ax2, f, m_RL, '--', 'Color', colX, 'LineWidth', 1.4, ...
                    'DisplayName', 'R beam → left ear (crosstalk)');
                grid(ax2,'on');
                xlim(ax2, f_lim); xlabel(ax2,'Frequency (Hz)'); ylabel(ax2,'Magnitude (dB)');
                title(ax2, sprintf('R beam  [%s / %s]', upper(r.wfs_mode), model));
                legend(ax2, 'show', 'Location', 'southwest', 'FontSize', 7);

                % ── Row 2: ICLD per beam ──
                ax3 = subplot(3,2,3, 'Parent', fig);
                semilogx(ax3, f, ICLD_L, 'Color', colL, 'LineWidth', 1.8);
                hold(ax3,'on');
                yline(ax3, 0, '--', 'Color', colX);
                grid(ax3,'on');
                xlim(ax3, f_lim); xlabel(ax3,'Frequency (Hz)'); ylabel(ax3,'ICLD (dB)');
                title(ax3, 'L beam: left/right ear ratio');

                ax4 = subplot(3,2,4, 'Parent', fig);
                semilogx(ax4, f, ICLD_R, 'Color', colR, 'LineWidth', 1.8);
                hold(ax4,'on');
                yline(ax4, 0, '--', 'Color', colX);
                grid(ax4,'on');
                xlim(ax4, f_lim); xlabel(ax4,'Frequency (Hz)'); ylabel(ax4,'ICLD (dB)');
                title(ax4, 'R beam: right/left ear ratio');

                % ── Row 3: Combined suppression (worst-case) ──
                ax5 = subplot(3,1,3, 'Parent', fig);
                ICLD_worst = min(ICLD_L, ICLD_R);
                semilogx(ax5, f, ICLD_L,     ':', 'Color', colL, 'LineWidth', 1.2, ...
                    'DisplayName', 'L beam ICLD');
                hold(ax5,'on');
                semilogx(ax5, f, ICLD_R,     ':', 'Color', colR, 'LineWidth', 1.2, ...
                    'DisplayName', 'R beam ICLD');
                semilogx(ax5, f, ICLD_worst,  '-', 'Color', [0.1 0.1 0.1], 'LineWidth', 2.0, ...
                    'DisplayName', 'Worst-case (min)');
                yline(ax5, 0, '--', 'Color', colX);
                grid(ax5,'on');
                xlim(ax5, f_lim); xlabel(ax5,'Frequency (Hz)'); ylabel(ax5,'ICLD (dB)');
                title(ax5, 'Combined crosstalk suppression');
                legend(ax5, 'show', 'Location', 'southwest', 'FontSize', 8);

            catch ME
                warndlg(sprintf('Error:\n%s', ME.message), 'Ear Transfer');
            end
        end

        function updateRenderer(obj)
            try
                obj.C.Scene.scene_renderer.update_renderer_settings(obj.C.Setup);
                obj.C.SceneGUI.draw_renderer_properties(obj.C.Scene.scene_renderer.SFS_renderer{1});
                obj.refreshAngle();
            catch
            end
        end

        function tryDrawProps(obj)
            try
                obj.C.SceneGUI.draw_renderer_properties(obj.C.Scene.scene_renderer.SFS_renderer{1});
            catch
            end
            % raiseSceneObjects() pushes VS/LS patches on top, which blocks
            % clicks on beam target patches. Re-raise the beam targets here.
            obj.raiseBeamTargets();
        end

        function raiseBeamTargets(obj)
            try
                if ~isempty(obj.hBeamL) && isgraphics(obj.hBeamL)
                    uistack(obj.hBeamL, 'top');
                end
                if ~isempty(obj.hBeamR) && isgraphics(obj.hBeamR)
                    uistack(obj.hBeamR, 'top');
                end
                % Also raise the text labels
                ax = obj.C.SceneGUI.main_axes;
                txts = findobj(ax, '-regexp', 'Tag', '^beamL_txt$|^beamR_txt$');
                for k = 1:numel(txts)
                    uistack(txts(k), 'top');
                end
            catch
            end
        end

        function onOptimalSearch(obj)
            try
                r = obj.C.Scene.scene_renderer.SFS_renderer{1};
                if ~isa(r, 'beamformer_renderer')
                    warndlg('Active renderer is not a beamformer.','Optimal Search');
                    return;
                end
                if isempty(r.receiver)
                    warndlg('No receiver available.','Optimal Search');
                    return;
                end

                f_lo = str2double(get(obj.EdtFLow,'String'));
                f_hi = str2double(get(obj.EdtFHigh,'String'));
                if isnan(f_lo) || isnan(f_hi) || f_lo >= f_hi || f_lo < 20
                    warndlg('Invalid frequency range.','Optimal Search');
                    return;
                end

                % Get acoustic model choice
                modelIdx = get(obj.PopEarModel,'Value');
                if modelIdx == 1
                    model = 'analytic'; hrtf_ext = [];
                else
                    model = 'hrtf';
                    try hrtf_ext = obj.C.Scene.scene_renderer.HRTF_extrapolators;
                    catch
                        warndlg('No HRTF extrapolator available.','Optimal Search');
                        return;
                    end
                end

                % Get search metric choice
                metricIdx = get(obj.PopSearchMetric, 'Value');
                if metricIdx == 1
                    metric = 'icld';
                else
                    metric = 'condition';
                end

                % Run the search with current mode and selected acoustic model
                res = r.run_optimal_search(f_lo, f_hi, model, hrtf_ext, [], [], metric);

                % Place R beam at optimal position (right half-plane),
                % L beam at mirror position
                r.beamR_pos = res.opt_pos_R;
                r.beamL_pos = res.opt_pos_L;
                r.update_renderer();
                obj.C.Setup.Renderer_setup.BeamRPos = res.opt_pos_R;
                obj.C.Setup.Renderer_setup.BeamLPos = res.opt_pos_L;
                obj.syncBeamGraphics(r);
                obj.tryDrawProps();
                obj.refreshLowCutFreq();

                % Show interactive map
                obj.showInteractiveMap(res, r);

            catch ME
                warndlg(sprintf('Error:\n%s', ME.message),'Optimal Search');
            end
        end

        function showInteractiveMap(obj, res, renderer)
            fig = figure('Name','Interaural Suppression Map','NumberTitle','off',...
                'Position',[100 100 800 600]);

            ax = axes('Parent',fig);
            imagesc(ax, res.X(1,:), res.Y(:,1), res.suppression_map);
            set(ax,'YDir','normal');
            hold(ax,'on');
            colormap(ax, parula);
            cb = colorbar(ax);
            if isfield(res,'metric') && strcmpi(res.metric,'condition')
                ylabel(cb, 'Worst-case 1/\kappa (reciprocal condition)');
                metricLabel = '1/κ';
            else
                ylabel(cb, 'Min 1/3-oct |ICLD| (dB)');
                metricLabel = 'ICLD';
            end
            mx = max(res.suppression_map(:));
            if isfinite(mx) && mx > 0
                caxis(ax, [0, mx*1.05]);
            end

            % Dim excluded region
            excl = zeros(size(res.valid_mask));
            excl(~res.valid_mask) = 1;
            h_excl = imagesc(ax, res.X(1,:), res.Y(:,1), excl);
            set(h_excl, 'AlphaData', excl * 0.4);

            % Array elements
            plot(ax, res.x0(:,1), res.x0(:,2), 'ws', 'MarkerSize',10,...
                'MarkerFaceColor',[0.2 0.2 0.2], 'LineWidth',1.5);

            % Head
            th = linspace(0,2*pi,100);
            plot(ax, res.rec_pos(1)+res.head_radius*cos(th),...
                     res.rec_pos(2)+res.head_radius*sin(th),...
                     'w-','LineWidth',2);

            % Ears
            plot(ax, res.ear_L_pos(1), res.ear_L_pos(2), 'bo','MarkerSize',8,'MarkerFaceColor','b');
            plot(ax, res.ear_R_pos(1), res.ear_R_pos(2), 'ro','MarkerSize',8,'MarkerFaceColor','r');

            % Optimal beam markers (R on right, L mirrored)
            if isfield(res,'opt_pos_R')
                optR = res.opt_pos_R;
                optL = res.opt_pos_L;
            else
                optR = [res.opt_pos, 0];
                optL = [-res.opt_pos(1), res.opt_pos(2), 0];
            end
            h_opt = plot(ax, optR(1), optR(2), 'p','MarkerSize',18,...
                'MarkerFaceColor',[1 0.3 0.2],'MarkerEdgeColor','k','LineWidth',1.5);
            plot(ax, optL(1), optL(2), 'p','MarkerSize',18,...
                'MarkerFaceColor',[0.2 0.4 1.0],'MarkerEdgeColor','k','LineWidth',1.5);

            % Clickable marker (current selection) — for R beam
            h_sel = plot(ax, optR(1), optR(2), 'o','MarkerSize',12,...
                'MarkerEdgeColor',[0 1 0],'LineWidth',2.5,'Visible','off');

            title(ax, sprintf('Min %s — %s mode, %s model (%d-%d Hz) — Click to place beams', ...
                metricLabel, upper(res.mode), res.acoustic_model, res.f_low, res.f_high));
            xlabel(ax, 'x [m]'); ylabel(ax, 'y [m]');
            axis(ax, 'equal','tight');

            % Store references for click callback
            mapData.res = res;
            mapData.h_sel = h_sel;
            mapData.h_opt = h_opt;
            mapData.panel = obj;
            ax.UserData = mapData;

            set(fig, 'WindowButtonDownFcn', @(src,~) obj.onMapClick(src));
        end

        function onMapClick(obj, figSrc)
            ax = findobj(figSrc, 'Type', 'axes');
            if isempty(ax), return; end
            ax = ax(1);
            cp = get(ax, 'CurrentPoint');
            x_click = cp(1,1);
            y_click = cp(1,2);

            md = ax.UserData;
            if isempty(md), return; end

            % Check if click is within map bounds
            res = md.res;
            xl = get(ax,'XLim'); yl = get(ax,'YLim');
            if x_click < xl(1) || x_click > xl(2) || y_click < yl(1) || y_click > yl(2)
                return;
            end

            % Get suppression value at clicked point (nearest grid point)
            [~, ix] = min(abs(res.X(1,:) - x_click));
            [~, iy] = min(abs(res.Y(:,1) - y_click));
            val = res.suppression_map(iy, ix);

            % Update selection marker
            set(md.h_sel, 'XData', x_click, 'YData', y_click, 'Visible','on');

            % Update title with clicked info
            if isfinite(val)
                title(ax, sprintf('Clicked: (%.3f, %.3f) — ICLD = %.1f dB | Optimal: %.1f dB', ...
                    x_click, y_click, val, res.opt_val));
            end

            % Place R beam at clicked position, L beam mirrored on x=0
            r = obj.getRenderer();
            if ~isempty(r)
                r.beamR_pos = [x_click, y_click, 0];
                r.beamL_pos = [-x_click, y_click, 0];
                r.update_renderer();
                obj.C.Setup.Renderer_setup.BeamRPos = r.beamR_pos;
                obj.C.Setup.Renderer_setup.BeamLPos = r.beamL_pos;
                obj.syncBeamGraphics(r);
                obj.tryDrawProps();
                obj.refreshLowCutFreq();
            end
        end

        function moveVSto(obj, pos2d)
            try
                vs = obj.C.Scene.virtual_sources{1};
                vs.set_position([pos2d(1:2), vs.position(3)]);

                % Redraw and update renderer
                obj.C.Scene.scene_renderer.update_SFS_renderers('source_moved');
                gui = obj.C.SceneGUI;
                if ~isempty(gui.virtual_source_points) && ~isempty(gui.virtual_source_points{1})
                    vsp = gui.virtual_source_points{1};
                    dx = pos2d(1:2) - vsp.UserData.Origin;
                    vsp.UserData.Origin = vsp.UserData.Origin + dx;
                    vsp.Vertices(:,1:2) = vsp.Vertices(:,1:2) + dx;
                    if isfield(vsp.UserData,'text') && isgraphics(vsp.UserData.text)
                        vsp.UserData.text.Position(1:2) = pos2d(1:2) + [0.1, 0.1];
                    end
                end

                % Refresh overlay
                obj.tryDrawProps();
                obj.refreshAngle();
            catch
            end
        end
    end
end
