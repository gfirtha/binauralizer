classdef ModeMap
    %MODEMAP  Legacy facade — delegates to audioapp.RendererRegistry.
    %
    % Existing code that calls ModeMap.toInternal(), .toUI(), .uiList()
    % keeps working without changes.  New code should use
    % audioapp.RendererRegistry.instance() directly.

    methods (Static)
        function list = uiList()
            list = audioapp.RendererRegistry.instance().uiList();
        end

        function id = toInternal(name)
            id = audioapp.RendererRegistry.instance().resolveId(name);
        end

        function ui = toUI(nameOrId)
            ui = audioapp.RendererRegistry.instance().toUI(nameOrId);
        end

        function idx = uiIndex(uiListCell, nameOrId)
            if isempty(uiListCell), idx = 1; return; end
            targetUI = audioapp.util.ModeMap.toUI(nameOrId);
            idx = find(strcmp(uiListCell, targetUI), 1);
            if isempty(idx), idx = 1; end
        end
    end
end
