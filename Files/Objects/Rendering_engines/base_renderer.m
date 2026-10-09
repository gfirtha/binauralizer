classdef base_renderer < handle
    %BASE_RENDERER  Abstract base class for all SFS renderers.
    
    properties
        virtual_source
        ssd
        output_signal
        filter_version = 0   % bumped whenever the driving filters change
    end
    methods (Abstract)
        update_renderer()
        render()
        get_renderer_info()
    end
    methods
        function obj = base_renderer(virtual_source,SSD)
            obj.virtual_source = virtual_source;
            obj.ssd = SSD;
            for n = 1 : length(SSD.loudspeakers)
                obj.output_signal{n} = signal;
            end
        end
        
        function obj = add_output_to_ssd_signal(obj)
            for n = 1 : length(obj.output_signal)
                obj.ssd.loudspeakers{n}.source_signal.add_signals(...
                    obj.output_signal{n});
            end
        end

        function d = get_driving_filters(obj) %#ok<MANU>
            %GET_DRIVING_FILTERS  Per-loudspeaker FIRs [Ntaps x Nls] (SSD
            %  order) of a linear time-invariant renderer, used by the fast
            %  binaural path. [] = not available (the per-loudspeaker chain
            %  is used instead). Overriding renderers must bump
            %  filter_version whenever the returned filters change.
            d = [];
        end

        function store = drawOverlay(obj, ax, layers, settings, store)
            %DRAWOVERLAY  Draw renderer-specific overlays on the main axes.
            %
            %   store = renderer.drawOverlay(ax, layers, settings, store)
            %
            %   ax       : main axes handle
            %   layers   : struct with .bg, .mid, .fg hggroup handles
            %   settings : renderer-specific overlay settings struct
            %   store    : struct of cached graphics handles (pass in, get back)
            %
            %   Default: no overlay (subclasses override as needed).
            %   The returned store is kept by the GUI and passed back next frame.
        end

        function settings = defaultOverlaySettings(obj)
            %DEFAULTOVERLAYSETTINGS  Return default overlay visibility flags.
            %   Subclasses override to provide renderer-specific defaults.
            %   Example (WFS): struct('showRef',true,'showAmp',true, ...)
            settings = struct();
        end
    end
end
