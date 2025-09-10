classdef (Abstract) RendererPanel < handle
    properties
        C
        Root matlab.ui.container.Panel   % our panel root
    end
    methods
        function obj = RendererPanel(C), obj.C = C; end
        function delete(obj)
            if ~isempty(obj.Root) && isvalid(obj.Root)
                delete(obj.Root);
            end
        end
    end
    methods (Abstract)
        build(obj, parentPanel)   % create controls under parentPanel; set obj.Root
        syncFromSetup(obj)        % refresh controls from C.Setup
    end
end
