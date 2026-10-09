classdef HOAPanel < audioapp.ui.renderers.RendererPanel
    properties
        EdtOrder
    end

    methods
        function build(obj, parent)
            obj.Root = uipanel('Parent',parent,'Units','normalized','Position',[0 0 1 1],'BorderType','none');

            uicontrol(obj.Root,'Style','text','String','HOA Order:',...
                'Units','normalized','Position',[0.07 0.80 0.40 0.10],...
                'HorizontalAlignment','left');

            obj.EdtOrder = uicontrol(obj.Root,'Style','edit',...
                'Units','normalized','Position',[0.50 0.80 0.20 0.10],...
                'Callback', @(s,~)obj.onOrderChanged(s));

            obj.syncFromSetup();
        end

        function syncFromSetup(obj)
            set(obj.EdtOrder,'String',num2str(obj.C.Setup.Renderer_setup.HOAorder));
        end
    end

    methods (Access=private)
        function onOrderChanged(obj, src)
            v = str2double(get(src,'String'));
            if ~isnan(v) && isfinite(v) && v>=0 && mod(v,1)==0
                obj.C.Setup.Renderer_setup.HOAorder = v;
                obj.updateRenderer();
            else
                set(src,'String',num2str(obj.C.Setup.Renderer_setup.HOAorder));
            end
        end

        function updateRenderer(obj)
            try
                obj.C.Scene.scene_renderer.update_renderer_settings(obj.C.Setup.Renderer_setup);
            catch
                obj.C.rebuildSoundScene();
            end
            obj.tryDrawProps();
        end

        function tryDrawProps(obj)
            try
                obj.C.SceneGUI.draw_renderer_properties(obj.C.Scene.scene_renderer.SFS_renderer{1});
            catch
            end
        end
    end
end
