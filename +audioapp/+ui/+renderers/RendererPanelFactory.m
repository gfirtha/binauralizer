classdef RendererPanelFactory
    % RendererPanelFactory  Creates the correct UI panel for a renderer mode.
    %
    % Delegates entirely to audioapp.RendererRegistry so that adding a new
    % renderer never requires touching this file.

    methods (Static)
        function pnl = create(modeUI, parent, C)
            pnl = audioapp.RendererRegistry.instance().createPanel(modeUI, parent, C);
        end
    end
end
