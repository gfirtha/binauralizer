classdef sound_scene_renderer < handle
    %RENDERER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        wfs_renderer
        binaural_renderer
    end
    
    methods
        function obj = sound_scene_renderer(virtual_sources,binaural_sources,receiver)
            if isempty(virtual_sources) % direct binauralization scenario
                for n = 1 : length(binaural_sources)
                    obj.binaural_renderer{n} = binaural_renderer(binaural_sources{n}, receiver);
                end
                
            else % virtual sound field synthesis scenario
                for n = 1 : length(virtual_sources)
                    obj.wfs_renderer{n} = wfs_renderer(virtual_sources{n}, binaural_sources);
                end
            end
            for n = 1 : length(binaural_sources)
                obj.binaural_renderer{n} = binaural_renderer(binaural_sources{n}, receiver);
            end
        end
        function update_wfs_renderers(obj, idx)
            obj.wfs_renderer{idx}.update_driving_function;
        end
        
        function update_binaural_renderers(obj, idx)
            obj.binaural_renderer{idx}.update_hrtf;
        end
        
        function output = render(obj, input)
            %% Only binauralization job
            if (isempty(obj.wfs_renderer))
               for n = 1 : length(obj.binaural_renderer)
                   obj.binaural_renderer{n}.binaural_source.set_input(input(:,n));
               end
                cellfun( @(x) x.render , obj.binaural_renderer );
                out_cell = cellfun( @(x) x.output_signal, obj.binaural_renderer, 'UniformOutput' , false);
                output = sum( cat( 3, out_cell{:} ) , 3);
            %% Sound field synthesis job
            else
                for m = 1 : length(obj.wfs_renderer)
                    obj.wfs_renderer{m}.virtual_source.set_input(input(:,m));
                    obj.wfs_renderer{m}.render;
                end
                for n = 1 : length(obj.binaural_renderer)
                    obj.binaural_renderer{n}.binaural_source.set_input(0);
                    for m = 1 : length(obj.wfs_renderer)
                        obj.binaural_renderer{n}.binaural_source.set_input(...
                            obj.binaural_renderer{n}.binaural_source.source_signal +...
                            obj.wfs_renderer{m}.output_signal{n});
                    end
                end
                cellfun( @(x) x.render , obj.binaural_renderer );
                out_cell = cellfun( @(x) x.output_signal, obj.binaural_renderer, 'UniformOutput' , false);
                output = sum( cat( 3, out_cell{:} ) , 3);
                
            end
        end
    end
end