classdef TDPanel < audioapp.ui.renderers.RendererPanel
    % Time-Delay Stereo renderer control panel
    % Exposes:
    %   TDNumTaps      (edit)
    %   TDUseLinear    (checkbox)
    %   TDRefPoint     (two edits)
    %   TDRefRadiusR0  ('auto' | numeric via popup + edit)

    properties (Constant)
        Mode = 'time_delay';
    end

    properties
        % UI
        EdtTaps
        ChkLinear
        EdtRefX
        EdtRefY
        PopR0Mode
        EdtR0

        SetupLsnr event.listener
    end

    methods
        function obj = TDPanel(C)
            obj@audioapp.ui.renderers.RendererPanel(C);
        end

        function build(obj, parent)
            obj.Root = uipanel('Parent',parent,'Units','normalized',...
                'Position',[0 0 1 1],'BorderType','etchedin','Title','Time-Delay Stereo');

            obj.SetupLsnr = addlistener(obj.C,'SetupChanged',@(~,~)obj.syncFromSetup());

            fs = 10;                      % compact font
            x1 = 0.06; x2 = 0.48; w1=0.40; w2=0.18; h=0.10; dy=0.12; y=0.82;

            % ---- FIR taps -------------------------------------------------
            uicontrol(obj.Root,'Style','text','String','FIR taps (N):',...
                'Units','normalized','Position',[x1 y w1 h],'HorizontalAlignment','left','FontSize',fs);
            obj.EdtTaps = uicontrol(obj.Root,'Style','edit','Units','normalized',...
                'Position',[x2 y w2 h],'FontSize',fs,'Callback',@(s,~)obj.onTaps(s));
            y = y - dy;

            % ---- Linear phase checkbox -----------------------------------
            obj.ChkLinear = uicontrol(obj.Root,'Style','checkbox','String','Linear phase ramp',...
                'Units','normalized','Position',[x1 y 0.60 h],'FontSize',fs,...
                'Callback',@(s,~)obj.onLinear(s));
            y = y - dy;

            % ---- Ref point ------------------------------------------------
            uicontrol(obj.Root,'Style','text','String','Ref point [x y] (m):',...
                'Units','normalized','Position',[x1 y w1 h],'HorizontalAlignment','left','FontSize',fs);
            obj.EdtRefX = uicontrol(obj.Root,'Style','edit','Units','normalized',...
                'Position',[x2 y w2 h],'FontSize',fs,'Callback',@(s,~)obj.onRefX(s));
            obj.EdtRefY = uicontrol(obj.Root,'Style','edit','Units','normalized',...
                'Position',[x2+w2+0.03 y w2 h],'FontSize',fs,'Callback',@(s,~)obj.onRefY(s));
            y = y - dy;

            % ---- R0 mode + value -----------------------------------------
            uicontrol(obj.Root,'Style','text','String','R0 (ref radius):',...
                'Units','normalized','Position',[x1 y w1 h],'HorizontalAlignment','left','FontSize',fs);

            obj.PopR0Mode = uicontrol(obj.Root,'Style','popupmenu','String',{'auto','manual'},...
                'Units','normalized','Position',[x2 y w2 h],'FontSize',fs,...
                'Callback',@(s,~)obj.onR0Mode(s));

            obj.EdtR0 = uicontrol(obj.Root,'Style','edit','Units','normalized',...
                'Position',[x2+w2+0.03 y w2 h],'FontSize',fs,'Callback',@(s,~)obj.onR0Value(s));

            obj.syncFromSetup();    % initial values
        end

        function syncFromSetup(obj)
            RS = obj.safeRS(obj.C.Setup);

            set(obj.EdtTaps,   'String', num2str(RS.TDNumTaps));
            set(obj.ChkLinear, 'Value',  logical(RS.TDUseLinear));
            set(obj.EdtRefX,   'String', num2str(RS.TDRefPoint(1)));
            set(obj.EdtRefY,   'String', num2str(RS.TDRefPoint(2)));

            % R0 mode/value
            if isempty(RS.TDRefRadiusR0) || (ischar(RS.TDRefRadiusR0) && strcmpi(RS.TDRefRadiusR0,'auto'))
                set(obj.PopR0Mode,'Value',1);   % auto
                set(obj.EdtR0,'String','','Enable','off');
            else
                set(obj.PopR0Mode,'Value',2);   % manual
                set(obj.EdtR0,'String',num2str(RS.TDRefRadiusR0),'Enable','on');
            end
        end
    end

    % ----- callbacks -------------------------------------------------------
    methods (Access=private)
        function onTaps(obj, src)
            v = str2double(get(src,'String'));
            if isfinite(v) && v>8 && mod(v,1)==0
                obj.C.Setup.Renderer_setup.TDNumTaps = v;
                obj.updateRenderer();
            else
                obj.syncFromSetup();
            end
        end

        function onLinear(obj, src)
            obj.C.Setup.Renderer_setup.TDUseLinear = logical(get(src,'Value'));
            obj.updateRenderer();
        end

        function onRefX(obj, src)
            v = str2double(get(src,'String'));
            if isfinite(v)
                rp = obj.safeRS(obj.C.Setup).TDRefPoint; rp(1) = v;
                obj.C.Setup.Renderer_setup.TDRefPoint = rp;
                obj.updateRenderer();
            else
                obj.syncFromSetup();
            end
        end

        function onRefY(obj, src)
            v = str2double(get(src,'String'));
            if isfinite(v)
                rp = obj.safeRS(obj.C.Setup).TDRefPoint; rp(2) = v;
                obj.C.Setup.Renderer_setup.TDRefPoint = rp;
                obj.updateRenderer();
            else
                obj.syncFromSetup();
            end
        end

        function onR0Mode(obj, src)
            isAuto = get(src,'Value')==1;
            if isAuto
                obj.C.Setup.Renderer_setup.TDRefRadiusR0 = 'auto';
                set(obj.EdtR0,'Enable','off');
            else
                % switch to manual; keep current value or suggest something
                RS = obj.safeRS(obj.C.Setup);
                if isempty(RS.TDRefRadiusR0) || ~isfinite(RS.TDRefRadiusR0)
                    obj.C.Setup.Renderer_setup.TDRefRadiusR0 = 1.0;
                end
                set(obj.EdtR0,'Enable','on');
            end
            obj.updateRenderer();
            obj.syncFromSetup();
        end

        function onR0Value(obj, src)
            v = str2double(get(src,'String'));
            if isfinite(v) && v>0
                obj.C.Setup.Renderer_setup.TDRefRadiusR0 = v;
                obj.updateRenderer();
            else
                obj.syncFromSetup();
            end
        end

        function updateRenderer(obj)
            % Ask current scene renderer to apply new settings (or rebuild)
            try
                obj.C.Scene.scene_renderer.update_renderer_settings(obj.C.Setup.Renderer_setup);
            catch
                obj.C.rebuildSoundScene();
            end
        end

        function RS = safeRS(~, Setup)
            % Ensure defaults so the UI is robust
            RS = Setup.Renderer_setup;
            if ~isfield(RS,'TDNumTaps')      || isempty(RS.TDNumTaps),      RS.TDNumTaps = 512; end
            if ~isfield(RS,'TDUseLinear')    || isempty(RS.TDUseLinear),    RS.TDUseLinear = true; end
            if ~isfield(RS,'TDRefPoint')     || numel(RS.TDRefPoint)<2,     RS.TDRefPoint = [0 0]; end
            if ~isfield(RS,'TDRefRadiusR0'), RS.TDRefRadiusR0 = 'auto'; end
        end
    end
end
