classdef CTCPanel < audioapp.ui.renderers.RendererPanel
    properties
        % Controls
        PopPlant
        PopVS
        EdtLambda
        SldLambda
        EdtR
        SldR
        EdtN
        Note

        SetupLsnr event.listener
    end

    methods
        function obj = CTCPanel(C)
            obj@audioapp.ui.renderers.RendererPanel(C);
            obj.SetupLsnr = addlistener(obj.C, 'SetupChanged', @(~,~)obj.syncFromSetup());
        end

        function build(obj, parent)
            obj.Root = uipanel('Parent',parent,'Units','normalized','Position',[0 0 1 1], ...
                               'BorderType','etchedin','Title','CTC controls');

            % Layout
            x1=0.06; x2=0.46; x3=0.73; w1=0.36; w2=0.24; h=0.11; y=0.82;

            % Plant model
            uicontrol(obj.Root,'Style','text','String','Plant model:', ...
                'Units','normalized','Position',[x1 y w1 h],'HorizontalAlignment','left');
            obj.PopPlant = uicontrol(obj.Root,'Style','popupmenu', ...
                'String',{'HRTF','point_source','rigid_sphere'}, ...
                'Units','normalized','Position',[x2 y w2 h], ...
                'Callback',@(s,~)obj.onPlant(getChoice(s)));

            % VS model
            y = y-0.16;
            uicontrol(obj.Root,'Style','text','String','VS→ear model:', ...
                'Units','normalized','Position',[x1 y w1 h],'HorizontalAlignment','left');
            obj.PopVS = uicontrol(obj.Root,'Style','popupmenu', ...
                'String',{'HRTF','point_source','rigid_sphere'}, ...
                'Units','normalized','Position',[x2 y w2 h], ...
                'Callback',@(s,~)obj.onVS(getChoice(s)));

            % Lambda
            y = y-0.16;
            uicontrol(obj.Root,'Style','text','String','Tikhonov λ:', ...
                'Units','normalized','Position',[x1 y w1 h],'HorizontalAlignment','left');
            obj.SldLambda = uicontrol(obj.Root,'Style','slider','Min',0,'Max',1e-2,'Value',1e-5, ...
                'Units','normalized','Position',[x1 y-0.06 w1 0.06], ...
                'SliderStep',[0.01 0.10], ...
                'Callback',@(s,~)obj.onLambda(get(s,'Value')));
            obj.EdtLambda = uicontrol(obj.Root,'Style','edit','String','1e-5', ...
                'Units','normalized','Position',[x3 y w2 h], ...
                'Callback',@(s,~)obj.onLambda(str2double(get(s,'String'))));

            % r_head
            y = y-0.16;
            uicontrol(obj.Root,'Style','text','String','Head radius r (m):', ...
                'Units','normalized','Position',[x1 y w1 h],'HorizontalAlignment','left');
            obj.SldR = uicontrol(obj.Root,'Style','slider','Min',0,'Max',0.2,'Value',0.10, ...
                'Units','normalized','Position',[x1 y-0.06 w1 0.06], ...
                'SliderStep',[0.01 0.10], ...
                'Callback',@(s,~)obj.onR(get(s,'Value')));
            obj.EdtR = uicontrol(obj.Root,'Style','edit','String','0.10', ...
                'Units','normalized','Position',[x3 y w2 h], ...
                'Callback',@(s,~)obj.onR(str2double(get(s,'String'))));

            % N taps (when not HRTF)
            y = y-0.16;
            uicontrol(obj.Root,'Style','text','String','Filter length N:', ...
                'Units','normalized','Position',[x1 y w1 h],'HorizontalAlignment','left');
            obj.EdtN = uicontrol(obj.Root,'Style','edit','String','1024', ...
                'Units','normalized','Position',[x2 y w2 h], ...
                'Callback',@(s,~)obj.onN(str2double(get(s,'String'))));

            % Note
            obj.Note = uicontrol(obj.Root,'Style','text', ...
                'String','If Plant or VS model = HRTF, filters adopt the IR length of Setup.HRTF.', ...
                'Units','normalized','Position',[x1 0.05 0.88 0.10], ...
                'HorizontalAlignment','left');

            obj.syncFromSetup();

            function v = getChoice(src)
                s = get(src,'String'); v = s{get(src,'Value')};
            end
            function v = get(h, p), v = get(h,p); end %#ok<*GETPROP>
        end

        function syncFromSetup(obj)
            rs = obj.C.Setup.Renderer_setup;
            rs = ensureDefaults(rs);

            % reflect
            set(obj.PopPlant,'Value', find(strcmpi(get(obj.PopPlant,'String'), rs.CTCPlantModel),1,'first'));
            set(obj.PopVS,   'Value', find(strcmpi(get(obj.PopVS,'String'),    rs.CTCVSModel),1,'first'));
            set(obj.SldLambda,'Value', rs.CTCLambda);
            set(obj.EdtLambda,'String', num2str(rs.CTCLambda,'%0.3g'));
            set(obj.SldR,'Value', rs.CTCRhead);
            set(obj.EdtR,'String', num2str(rs.CTCRhead,'%0.3g'));
            set(obj.EdtN,'String', num2str(rs.CTCNumTaps));

            % enable N only when neither model is HRTF
            needN = ~strcmpi(rs.CTCPlantModel,'hrtf') && ~strcmpi(rs.CTCVSModel,'hrtf');
            set(obj.EdtN,'Enable', ternary(needN,'on','off'));

            function y = ternary(c,a,b), if c, y=a; else, y=b; end; end
            function rs = ensureDefaults(rs)
                if ~isfield(rs,'CTCPlantModel'), rs.CTCPlantModel = 'HRTF'; end
                if ~isfield(rs,'CTCVSModel'),    rs.CTCVSModel    = 'HRTF'; end
                if ~isfield(rs,'CTCNumTaps'),    rs.CTCNumTaps    = 1024;   end
                if ~isfield(rs,'CTCLambda'),     rs.CTCLambda     = 1e-5;   end
                if ~isfield(rs,'CTCRhead'),      rs.CTCRhead      = 0.10;   end
                obj.C.Setup.Renderer_setup = rs;
            end
        end

        function delete(obj)
            if ~isempty(obj.SetupLsnr) && isvalid(obj.SetupLsnr), delete(obj.SetupLsnr); end
        end
    end

    methods (Access=private)
        function onPlant(obj, val)
            obj.C.Setup.Renderer_setup.CTCPlantModel = val;
            obj.push();
        end
        function onVS(obj, val)
            obj.C.Setup.Renderer_setup.CTCVSModel = val;
            obj.push();
        end
        function onLambda(obj, v)
            if ~isfinite(v) || v<0, v = obj.C.Setup.Renderer_setup.CTCLambda; end
            obj.C.Setup.Renderer_setup.CTCLambda = v;
            set(obj.SldLambda,'Value',v); set(obj.EdtLambda,'String',num2str(v,'%0.3g'));
            obj.push();
        end
        function onR(obj, v)
            if ~isfinite(v) || v<0, v = obj.C.Setup.Renderer_setup.CTCRhead; end
            obj.C.Setup.Renderer_setup.CTCRhead = v;
            set(obj.SldR,'Value',v); set(obj.EdtR,'String',num2str(v,'%0.3g'));
            obj.push();
        end
        function onN(obj, v)
            if ~isfinite(v) || v<1, v = obj.C.Setup.Renderer_setup.CTCNumTaps; end
            obj.C.Setup.Renderer_setup.CTCNumTaps = round(v);
            set(obj.EdtN,'String',num2str(obj.C.Setup.Renderer_setup.CTCNumTaps));
            obj.push();
        end

        function push(obj)
            try
                obj.C.Scene.scene_renderer.update_renderer_settings(obj.C.Setup.Renderer_setup);
            catch
                obj.C.rebuildSoundScene();
            end
        end
    end
end
