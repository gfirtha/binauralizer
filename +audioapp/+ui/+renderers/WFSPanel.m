classdef WFSPanel < audioapp.ui.renderers.RendererPanel
    properties
        ChkAA
        ChkLowCut
        TxtLowCutFreq
        EdtTaper
        PopRefMode
        PnlRef
        ChkRef
        ChkAmp
        ChkTaper
        ChkZones
        Visuals

        SetupLsnr event.listener   % <--- NEW
        GeomLsnr  event.listener   % GeometryChanged listener
    end

    methods
        function build(obj, parent)
            obj.SetupLsnr = addlistener(obj.C, 'SetupChanged', @(~,~)obj.rebuildRefArea());
            obj.GeomLsnr  = addlistener(obj.C, 'GeometryChanged', @(~,~)obj.refreshLowCutFreq());

            obj.Root = uipanel('Parent',parent,'Units','normalized','Position',[0 0.3 1 0.7],'BorderType','etchedin','Title','Controls');

            uicontrol(obj.Root,'Style','text','String','Antialiasing filter:',...
                'Units','normalized','Position',[0.07 0.88 0.50 0.08],...
                'HorizontalAlignment','left');

            obj.ChkAA = uicontrol(obj.Root,'Style','checkbox','String','Enabled',...
                'Units','normalized','Position',[0.60 0.88 0.30 0.08],...
                'Callback', @(s,~)obj.onAntialias(s));

            uicontrol(obj.Root,'Style','text','String','Tapering coefficient (0..1):',...
                'Units','normalized','Position',[0.07 0.78 0.60 0.08],...
                'HorizontalAlignment','left');

            obj.EdtTaper = uicontrol(obj.Root,'Style','edit',...
                'Units','normalized','Position',[0.70 0.78 0.20 0.08],...
                'Callback', @(s,~)obj.onTaper(s));

            % Low-freq pre-eq shelving
            obj.ChkLowCut = uicontrol(obj.Root,'Style','checkbox','String','Low-freq pre-eq shelving',...
                'Units','normalized','Position',[0.07 0.68 0.85 0.08],...
                'Callback', @(s,~)obj.onLowCut(s));
            obj.TxtLowCutFreq = uicontrol(obj.Root,'Style','text','String','fc: — Hz',...
                'Units','normalized','Position',[0.10 0.60 0.85 0.07],...
                'HorizontalAlignment','left','ForegroundColor',[0.4 0.4 0.4],'FontSize',8);

            uicontrol(obj.Root,'Style','text','String','Referencing setup:',...
                'Units','normalized','Position',[0.07 0.48 0.50 0.08],...
                'HorizontalAlignment','left');

            obj.PopRefMode = uicontrol(obj.Root,'Style','popupmenu',...
                'String',{'Concentric/parallel contour','Fixed reference point','Fixed reference distance'},...
                'Units','normalized','Position',[0.4 0.48 0.50 0.08],...
                'Callback', @(s,~)obj.onRefMode(s));

            obj.PnlRef = uipanel('Parent',obj.Root,'Units','normalized','Position',[0.04 0.18 0.92 0.28],'BorderType','none');

            obj.Visuals = uipanel('Parent',parent,'Units','normalized','Position',[0 0 1 0.3],'BorderType','etchedin','Title','Visualization');

            % uicontrol(obj.Root,'Style','text','String','Visualization:',...
            %     'Units','normalized','Position',[0.07 0.36 0.50 0.08],...
            %     'HorizontalAlignment','left','FontWeight','bold');

            obj.ChkRef = uicontrol(obj.Visuals,'Style','checkbox','String','Reference contour',...
                'Units','normalized','Position',[0.10 0.80 0.70 0.08],...
                'Callback', @(s,~)obj.onViz('ref',s));

            obj.ChkAmp = uicontrol(obj.Visuals,'Style','checkbox','String','Amplitude distribution',...
                'Units','normalized','Position',[0.10 0.60 0.70 0.08],...
                'Callback', @(s,~)obj.onViz('amp',s));

            obj.ChkTaper = uicontrol(obj.Visuals,'Style','checkbox','String','Tapering window',...
                'Units','normalized','Position',[0.10 0.40 0.70 0.08],...
                'Callback', @(s,~)obj.onViz('taper',s));

            obj.ChkZones = uicontrol(obj.Visuals,'Style','checkbox','String','Illuminated / taper zones',...
                'Units','normalized','Position',[0.10 0.2 0.70 0.08],...
                'Callback', @(s,~)obj.onViz('zones',s));

            obj.syncFromSetup();    % initial values
            obj.rebuildRefArea();   % build the subarea according to mode+geom
        end

        function syncFromSetup(obj)
            rs = obj.C.Setup.Renderer_setup;
            set(obj.ChkAA,   'Value', logical(rs.Antialiasing));
            set(obj.EdtTaper,'String', num2str(rs.Tapering));
            set(obj.PopRefMode,'Value', obj.refModeToVal(rs.ReferenceMode));

            % Low-cut sync
            lowcut = false;
            if isfield(rs, 'LowFreqCut'), lowcut = logical(rs.LowFreqCut); end
            set(obj.ChkLowCut, 'Value', lowcut);
            obj.applyLowCutToWFS(lowcut);
            obj.refreshLowCutFreq();

            % viz toggles from SceneGUI flags (if available)
            s = obj.C.SceneGUI;
            try
                set(obj.ChkRef,  'Value', s.rendererSetting.wfs_renderer.showRef);
                set(obj.ChkAmp,  'Value', s.rendererSetting.wfs_renderer.showAmp);
                set(obj.ChkTaper,'Value', s.rendererSetting.wfs_renderer.showTaper);
                set(obj.ChkZones,'Value', s.rendererSetting.wfs_renderer.showZones);
            catch
            end
        end
    end

    % ----- callbacks -----
    methods (Access=private)
        function onAntialias(obj, src)
            obj.C.Setup.Renderer_setup.Antialiasing = logical(get(src,'Value'));
            obj.updateRenderer();
        end
        function onLowCut(obj, src)
            val = logical(get(src,'Value'));
            obj.C.Setup.Renderer_setup.LowFreqCut = val;
            obj.applyLowCutToWFS(val);
            obj.updateRenderer();
            obj.refreshLowCutFreq();
        end
        function applyLowCutToWFS(obj, enabled)
            r = obj.getWfsRenderer();
            if isempty(r), return; end
            r.LowCutEnabled = enabled;
            % Always set receiver_ref so cutoff is computed and displayed
            r.receiver_ref = obj.C.Scene.receiver;
            r.update_renderer();
        end
        function refreshLowCutFreq(obj)
            r = obj.getWfsRenderer();
            if isempty(r)
                set(obj.TxtLowCutFreq, 'String', 'fc: — Hz');
                return;
            end
            try
                fc_lo = r.lower_cutoff_freq;
                fc_hi = r.higher_cutoff_freq;
                if fc_lo > 0 && isfinite(fc_hi)
                    set(obj.TxtLowCutFreq, 'String', ...
                        sprintf('fc%s%.0f  fc%s%.0f Hz', char(8595), fc_lo, char(8593), fc_hi));
                else
                    set(obj.TxtLowCutFreq, 'String', 'fc: — Hz');
                end
            catch
                set(obj.TxtLowCutFreq, 'String', 'fc: — Hz');
            end
        end
        function r = getWfsRenderer(obj)
            r = [];
            try
                rr = obj.C.Scene.scene_renderer.SFS_renderer;
                for k = 1:numel(rr)
                    if isa(rr{k}, 'wfs_renderer')
                        r = rr{k};
                        return;
                    end
                end
            catch, end
        end
        function onTaper(obj, src)
            v = str2double(get(src,'String'));
            if ~isnan(v) && isfinite(v) && v>=0 && v<=1
                obj.C.Setup.Renderer_setup.Tapering = v;
                obj.updateRenderer();
            else
                set(src,'String',num2str(obj.C.Setup.Renderer_setup.Tapering));
            end
        end
        function onRefMode(obj, src)
            modes = {'concentric','fixed_point','fixed_distance'};
            obj.C.Setup.Renderer_setup.ReferenceMode = modes{ get(src,'Value') };
            obj.rebuildRefArea();
            obj.updateRenderer();
        end
        function onViz(obj, what, src)
            val = logical(get(src,'Value'));
            s = obj.C.SceneGUI;
            try
                switch what
                    case 'ref',   s.rendererSetting.wfs_renderer.showRef   = val; set(s.renderer_illustration.wfsRefContour,'Visible',val);
                    case 'amp',   s.rendererSetting.wfs_renderer.showAmp   = val; set(s.renderer_illustration.amplitudeDistribution,'Visible',val);
                    case 'taper', s.rendererSetting.wfs_renderer.showTaper = val; set(s.renderer_illustration.windowFunction,'Visible',val);
                    case 'zones'
                        s.rendererSetting.wfs_renderer.showZones = val;
                        f = s.renderer_illustration;
                        if isfield(f,'zonePolys') && iscell(f.zonePolys)
                            cellfun(@(h) set(h,'Visible',val), f.zonePolys(cellfun(@isgraphics, f.zonePolys)));
                        end
                        if isfield(f,'zoneRays') && iscell(f.zoneRays)
                            cellfun(@(h) set(h,'Visible',val), f.zoneRays(cellfun(@isgraphics, f.zoneRays)));
                        end
                end
                obj.tryDrawProps();
            catch
            end
        end
        function rebuildRefArea(obj)
            delete(allchild(obj.PnlRef));
            geom = lower(obj.C.Setup.Loudspeaker_setup.Shape);
            mode = obj.C.Setup.Renderer_setup.ReferenceMode;

            switch mode
                case 'concentric'
                    if any(strcmp(geom,{'linear','stereo'}))
                        uicontrol(obj.PnlRef,'Style','text','String','normal distance (m):',...
                            'Units','normalized','Position',[0.05 0.55 0.40 0.30],...
                            'HorizontalAlignment','left');
                        uicontrol(obj.PnlRef,'Style','edit','String',num2str(obj.C.Setup.Renderer_setup.DeltaY),...
                            'Units','normalized','Position',[0.50 0.60 0.20 0.25],...
                            'Callback',@(s,~)obj.onDeltaY(s));
                    else
                        uicontrol(obj.PnlRef,'Style','text','String','Reference radius r_{ref} (m):',...
                            'Units','normalized','Position',[0.05 0.55 0.60 0.30],...
                            'HorizontalAlignment','left');
                        uicontrol(obj.PnlRef,'Style','edit','String',num2str(obj.C.Setup.Renderer_setup.Rref),...
                            'Units','normalized','Position',[0.68 0.60 0.20 0.25],...
                            'Callback',@(s,~)obj.onRref(s));
                    end
                case 'fixed_point'
                    uicontrol(obj.PnlRef,'Style','text','String','Reference point (x, y) [m]:',...
                        'Units','normalized','Position',[0.05 0.55 0.60 0.30],'HorizontalAlignment','left');
                    uicontrol(obj.PnlRef,'Style','edit','String',num2str(obj.C.Setup.Renderer_setup.RefPoint(1)),...
                        'Units','normalized','Position',[0.50 0.60 0.16 0.25],...
                        'Callback',@(s,~)obj.onRefPointX(s));
                    uicontrol(obj.PnlRef,'Style','edit','String',num2str(obj.C.Setup.Renderer_setup.RefPoint(2)),...
                        'Units','normalized','Position',[0.70 0.60 0.16 0.25],...
                        'Callback',@(s,~)obj.onRefPointY(s));
                case 'fixed_distance'
                    uicontrol(obj.PnlRef,'Style','text','String','Reference distance [m]:',...
                        'Units','normalized','Position',[0.05 0.55 0.60 0.30],'HorizontalAlignment','left');
                    uicontrol(obj.PnlRef,'Style','edit','String',num2str(obj.C.Setup.Renderer_setup.RefDistance),...
                        'Units','normalized','Position',[0.50 0.60 0.16 0.25],...
                        'Callback',@(s,~)obj.onRefDistance(s));
            end
        end
        function onDeltaY(obj, src), v=str2double(get(src,'String')); if ~isnan(v), obj.C.Setup.Renderer_setup.DeltaY=v; obj.updateRenderer(); else, set(src,'String',num2str(obj.C.Setup.Renderer_setup.DeltaY)); end, end
        function onRref(obj, src),   v=str2double(get(src,'String')); if ~isnan(v)&&v>0, obj.C.Setup.Renderer_setup.Rref=v; obj.updateRenderer(); else, set(src,'String',num2str(obj.C.Setup.Renderer_setup.Rref)); end, end
        function onRefPointX(obj, src), v=str2double(get(src,'String')); if ~isnan(v), rp=obj.C.Setup.Renderer_setup.RefPoint; rp(1)=v; obj.C.Setup.Renderer_setup.RefPoint=rp; obj.updateRenderer(); else, set(src,'String',num2str(obj.C.Setup.Renderer_setup.RefPoint(1))); end, end
        function onRefPointY(obj, src), v=str2double(get(src,'String')); if ~isnan(v), rp=obj.C.Setup.Renderer_setup.RefPoint; rp(2)=v; obj.C.Setup.Renderer_setup.RefPoint=rp; obj.updateRenderer(); else, set(src,'String',num2str(obj.C.Setup.Renderer_setup.RefPoint(2))); end, end
        function onRefDistance(obj, src), v=str2double(get(src,'String')); if ~isnan(v), obj.C.Setup.Renderer_setup.RefDistance=v; obj.updateRenderer(); else, set(src,'String',num2str(obj.C.Setup.Renderer_setup.RefDistance)); end, end

        function v = refModeToVal(~, m)
            switch lower(m)
                case 'concentric', v=1;
                case 'fixed_point', v=2;
                case 'fixed_distance', v=3;
                otherwise, v=1;
            end
        end
        function updateRenderer(obj)
            try
                obj.C.Scene.scene_renderer.update_renderer_settings(obj.C.Setup);
            catch
                obj.C.rebuildSoundScene();
            end
            obj.tryDrawProps();
            obj.refreshLowCutFreq();
        end
        function tryDrawProps(obj)
            try
                obj.C.SceneGUI.draw_renderer_properties(obj.C.Scene.scene_renderer.SFS_renderer{1});
            catch
            end
        end
    end
end
