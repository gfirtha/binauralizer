classdef DirectPanel < audioapp.ui.renderers.RendererPanel
    methods
        function build(obj, parent)
            obj.Root = uipanel('Parent',parent,'Units','normalized','Position',[0 0 1 1],'BorderType','none');

            uicontrol(obj.Root,'Style','text',...
                'String','Direct playback settings (stub).',...
                'Units','normalized','Position',[0.07 0.84 0.86 0.10],...
                'HorizontalAlignment','left');
        end
        function syncFromSetup(~), end
    end
end
