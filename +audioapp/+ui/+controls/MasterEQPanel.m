classdef MasterEQPanel < handle
    % audioapp.ui.controls.MasterEQPanel
    % Creates a bottom panel with:
    %  - Enable toggle
    %  - Mode popup (Graphic / Low-pass)
    %  - Presets popup
    %  - If Graphic: 10 sliders (-12..+12 dB)
    %  - If Low-pass: log slider + numeric edit for cutoff

    properties
        C
        Parent
        Panel

        % top controls
        ChkEnable
        PopMode
        PopPreset

        % containers
        PnlGraphic
        PnlLowpass

        % graphic sliders
        GraphicSliders

        % lowpass widgets
        LpSlider
        LpEdit
    end

    methods
        function obj = MasterEQPanel(parent, ctx)
            obj.Parent = parent;
            obj.C = ctx;
            obj.C.ensureEngine();
            if isempty(obj.C.Engine.MasterEQ) || ~isvalid(obj.C.Engine.MasterEQ)
                obj.C.Engine.MasterEQ = audioapp.dsp.MasterEQ(obj.C);
            end

            % host panel at bottom
            eqH = 0.26;
            obj.Panel = uipanel('Parent',parent,'Units','normalized', ...
                'Position',[0.02 0.02 0.96 eqH], ...
                'BorderType','etchedin','Title','Master EQ');

            % --- Top row ---------------------------------------------------
            obj.ChkEnable = uicontrol(obj.Panel,'Style','checkbox','String','Enable', ...
                'Units','normalized','Position',[0.02 0.8 0.17 0.18], ...
                'Value', obj.C.Engine.MasterEQ.Enabled, ...
                'Callback', @(s,~) obj.onEnable(s));

            uicontrol(obj.Panel,'Style','text','String','Mode:', ...
                'Units','normalized','Position',[0.23 0.80 0.10 0.16], 'HorizontalAlignment','left');

            obj.PopMode = uicontrol(obj.Panel,'Style','popupmenu', ...
                'String',{'10-band Graphic','Low-pass'}, ...
                'Value', 1 + strcmpi(obj.C.Engine.MasterEQ.Mode,'lowpass'), ...
                'Units','normalized','Position',[0.34 0.80 0.22 0.16], ...
                'Callback', @(s,~) obj.onMode(s));

            % If Audio TB is missing, pre-select Low-pass and disable graphic choice
            if ~obj.C.HasAudio
                set(obj.PopMode,'Value',2);
            end

            uicontrol(obj.Panel,'Style','text','String','Preset:', ...
                'Units','normalized','Position',[0.59 0.80 0.12 0.16], 'HorizontalAlignment','left');

            obj.PopPreset = uicontrol(obj.Panel,'Style','popupmenu', ...
                'String',{'(none)','Flat','Bass boost','Treble cut','Vocal clarity','Lowpass 8k','Lowpass 2k'}, ...
                'Value', 1,'Units','normalized','Position',[0.75 0.80 0.23 0.16], ...
                'Callback', @(s,~) obj.onPreset(s));

            % --- Mode containers ------------------------------------------
            obj.PnlGraphic = uipanel('Parent',obj.Panel,'Units','normalized','Position',[0.01 0.01 0.98 0.74], 'BorderType','none');
            obj.PnlLowpass = uipanel('Parent',obj.Panel,'Units','normalized','Position',[0.01 0.01 0.98 0.74], 'BorderType','none');

            % Graphic sliders
            bandLabels = {'31','62','125','250','500','1k','2k','4k','8k','16k'};
            G0 = obj.C.Engine.MasterEQ.GraphicGains; if numel(G0)~=10, G0 = zeros(1,10); end
            sW = 0.07; gap = 0.013; x0 = 0.04; y0 = 0.18; h = 0.76;
            obj.GraphicSliders = gobjects(1,10);
            for k = 1:10
                x = x0 + (k-1)*(sW+gap);
                obj.GraphicSliders(k) = uicontrol(obj.PnlGraphic,'Style','slider','Min',-12,'Max',+12, ...
                    'Value', G0(k), 'Units','normalized','Position',[x y0 sW h], ...
                    'Callback', @(~,~) obj.onGraphicGains());
                uicontrol(obj.PnlGraphic,'Style','text','String',bandLabels{k}, ...
                    'Units','normalized','Position',[x 0.02 sW 0.14], 'HorizontalAlignment','center');
            end

            % Lowpass widgets
            uicontrol(obj.PnlLowpass,'Style','text','String','Cutoff (Hz):', ...
                'Units','normalized','Position',[0.06 0.66 0.22 0.20], 'HorizontalAlignment','left');

            obj.LpSlider = uicontrol(obj.PnlLowpass,'Style','slider','Min',0,'Max',1,'Value',0.9, ...
                'Units','normalized','Position',[0.28 0.70 0.50 0.16], ...
                'Callback', @(s,~) obj.onLowpassSlider(s));
            obj.LpEdit = uicontrol(obj.PnlLowpass,'Style','edit','String',num2str(obj.C.Engine.MasterEQ.LowpassFc), ...
                'Units','normalized','Position',[0.80 0.68 0.16 0.20], ...
                'Callback', @(s,~) obj.onLowpassEdit(s));

            set(obj.LpSlider,'Value', obj.fc2slider(obj.C.Engine.MasterEQ.LowpassFc));
            set(obj.LpEdit,  'String', num2str(round(obj.C.Engine.MasterEQ.LowpassFc)));

            obj.refreshModePanels();
            obj.refreshEnableState();
        end

        % ======= Callbacks ===============================================
        function onEnable(obj, src)
            on = logical(get(src,'Value'));
            obj.C.Engine.MasterEQ.setEnabled(on);
            obj.refreshEnableState();
        end

        function onMode(obj, src)
            v = get(src,'Value');
            mode = ternary(v==1,'graphic','lowpass');
            if strcmp(mode,'graphic') && ~obj.C.HasAudio
                % No Audio TB → force lowpass
                mode = 'lowpass'; set(src,'Value',2);
            end
            obj.C.Engine.MasterEQ.setMode(mode);
            obj.refreshModePanels();
            obj.refreshEnableState();
        end

        function onPreset(obj, src)
            list = get(src,'String'); v = get(src,'Value');
            if v==1, return; end
            obj.C.Engine.MasterEQ.setPreset(list{v});
            % snap UI to new settings
            if strcmpi(obj.C.Engine.MasterEQ.Mode,'graphic')
                g = obj.C.Engine.MasterEQ.GraphicGains;
                for kk=1:10, set(obj.GraphicSliders(kk),'Value',g(kk)); end
            else
                set(obj.LpSlider,'Value', obj.fc2slider(obj.C.Engine.MasterEQ.LowpassFc));
                set(obj.LpEdit,  'String', num2str(round(obj.C.Engine.MasterEQ.LowpassFc)));
            end
            set(src,'Value',1);
        end

        function onGraphicGains(obj)
            g = zeros(1,10);
            for kk=1:10, g(kk) = get(obj.GraphicSliders(kk),'Value'); end
            obj.C.Engine.MasterEQ.setGraphicGains(g);
        end

        function onLowpassSlider(obj, s)
            fc = obj.slider2fc(get(s,'Value'));
            obj.C.Engine.MasterEQ.setLowpassFc(fc);
            set(obj.LpEdit,'String',num2str(round(fc)));
        end

        function onLowpassEdit(obj, s)
            fc = str2double(get(s,'String'));
            if ~isnan(fc) && isfinite(fc)
                obj.C.Engine.MasterEQ.setLowpassFc(fc);
                set(obj.LpSlider,'Value', obj.fc2slider(obj.C.Engine.MasterEQ.LowpassFc));
                set(s,'String', num2str(round(obj.C.Engine.MasterEQ.LowpassFc)));
            else
                set(s,'String', num2str(round(obj.slider2fc(get(obj.LpSlider,'Value')))));
            end
        end

        % ======= Helpers & UI refresh ====================================
        function refreshModePanels(obj)
            if strcmpi(obj.C.Engine.MasterEQ.Mode,'graphic')
                set(obj.PnlGraphic,'Visible','on'); set(obj.PnlLowpass,'Visible','off');
                set(obj.PopMode,'Value',1);
            else
                set(obj.PnlGraphic,'Visible','off'); set(obj.PnlLowpass,'Visible','on');
                set(obj.PopMode,'Value',2);
            end
        end

        function refreshEnableState(obj)
            on = obj.C.Engine.MasterEQ.Enabled;
            kids = [obj.PnlGraphic.Children; obj.PnlLowpass.Children; obj.PopMode; obj.PopPreset];
            set(kids, 'Enable', ternary(on,'on','off'));

            % If graphic EQ is not supported, keep the popup showing Low-pass
            if ~obj.C.HasAudio
                set(obj.PopMode,'Value',2);
            end
        end

        function fc = slider2fc(~, t)
            fmin=100; fmax=20000;
            fc = 10^( log10(fmin) + t*(log10(fmax)-log10(fmin)) );
        end
        function t = fc2slider(~, fc)
            fmin=100; fmax=20000; fc = max(fmin,min(fmax,fc));
            t = (log10(fc)-log10(fmin))/(log10(fmax)-log10(fmin));
        end
    end
end

function y = ternary(cond, a, b); if cond, y=a; else, y=b; end; end
