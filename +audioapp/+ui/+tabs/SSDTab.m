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
        end
    
        methods
            function obj = SSDTab(parentTab, ctx)
                obj.C = ctx;
                obj.Panel = parentTab;
                p = obj.Panel;
                obj.SetupLsnr = addlistener(obj.C, 'SetupChanged', @(~,~)obj.refreshFromSetup());
    
                uicontrol(p,'Style','text','String','Geometry:',...
                    'Units','normalized','Position',[0.07 0.85 0.30 0.08],...
                    'HorizontalAlignment','left');
    
                % add "General (from file)"
                obj.PopGeom = uicontrol(p,'Style','popupmenu',...
                    'String',{'Linear','Circular','Stereo','5.0 (circular)','General'},...
                    'Value',obj.geomToVal(obj.C.Setup.Loudspeaker_setup.Shape),...
                    'Units','normalized','Position',[0.40 0.85 0.50 0.08],...
                    'Callback',@(s,~)obj.onGeometry(s),'UserData',1);
    
                uicontrol(p,'Style','text','String','Size (R):',...
                    'Units','normalized','Position',[0.07 0.72 0.30 0.08],...
                    'HorizontalAlignment','left');
    
                obj.EdtSize = uicontrol(p,'Style','edit',...
                    'String',num2str(obj.C.Setup.Loudspeaker_setup.R),...
                    'Units','normalized','Position',[0.40 0.72 0.50 0.08],...
                    'Callback',@(s,~)obj.onSize(s));
    
                uicontrol(p,'Style','text','String','LS number (N):',...
                    'Units','normalized','Position',[0.07 0.60 0.30 0.08],...
                    'HorizontalAlignment','left');
    
                obj.EdtNum = uicontrol(p,'Style','edit',...
                    'String',num2str(obj.C.Setup.Loudspeaker_setup.N),...
                    'Units','normalized','Position',[0.40 0.60 0.50 0.08],...
                    'Callback',@(s,~)obj.onNum(s));
    
                uicontrol(p,'Style','text','String','Loudspeaker model:',...
                    'Units','normalized','Position',[0.07 0.48 0.30 0.08],...
                    'HorizontalAlignment','left');
    
                obj.PopLSModel = uicontrol(p,'Style','popupmenu',...
                    'String',{'Point source','Piston'},...
                    'Value', 1 + strcmpi(obj.C.Setup.Loudspeaker_type.Shape,'circular_piston'),...
                    'Units','normalized','Position',[0.40 0.48 0.50 0.08],...
                    'Callback',@(s,~)obj.onLSModel(s));
    
                uicontrol(p,'Style','text','String','LS radius (m):',...
                    'Units','normalized','Position',[0.07 0.36 0.30 0.08],...
                    'HorizontalAlignment','left');
    
                obj.EdtLSRadius = uicontrol(p,'Style','edit',...
                    'String',num2str(obj.C.Setup.Loudspeaker_type.R),...
                    'Units','normalized','Position',[0.40 0.36 0.50 0.08],...
                    'Callback',@(s,~)obj.onLSRadius(s));
    
                % NEW: Load geometry button (only visible for General)
                obj.BtnLoadGeom = uicontrol(p,'Style','pushbutton','String','Load geometry ...',...
                    'Units','normalized','Position',[0.40 0.26 0.50 0.08],...
                    'Visible','off', ...
                    'Callback',@(s,~)obj.onLoadGeometry());
    
                % ensure proper visibility at startup
                obj.updateGeneralButtonVisibility();
            end
    
            function refreshFromSetup(obj)
                S = obj.C.Setup;
    
                % geometry popup
                try
                    set(obj.PopGeom,'Value', obj.geomToVal(S.Loudspeaker_setup.Shape));
                catch
                    % if shape not in list (e.g., 'general'), clamp to 'General'
                    set(obj.PopGeom,'Value', 5);
                end
    
                % numeric edits
                set(obj.EdtSize,    'String', num2str(S.Loudspeaker_setup.R));
                set(obj.EdtNum,     'String', num2str(S.Loudspeaker_setup.N));
                set(obj.PopLSModel, 'Value', 1 + strcmpi(S.Loudspeaker_type.Shape,'circular_piston'));
                set(obj.EdtLSRadius,'String', num2str(S.Loudspeaker_type.R));
    
                % show/hide the “Load geometry …” button
                obj.updateGeneralButtonVisibility();
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
                    case 'linear',   v=1;
                    case 'circular', v=2;
                    case 'stereo',   v=3;
                    case '5.0',      v=4;
                    case 'general',  v=5;      % NEW
                    otherwise,       v=1;
                end
            end
        end
    
        
    end
    
    function y = ternary(c,a,b), if c, y=a; else, y=b; end, end
    

