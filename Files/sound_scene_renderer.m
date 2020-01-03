classdef sound_scene_renderer < handle
    %RENDERER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        wfs_renderer
        binaural_renderer
    end
    
    methods
        function obj = sound_scene_renderer(virtual_sources,binaural_sources,receiver, setup)
            if isempty(virtual_sources) % direct binauralization scenario
                for n = 1 : length(binaural_sources)
                    obj.binaural_renderer{n} = binaural_renderer(binaural_sources{n}, receiver);
                end
                
            else % virtual sound field synthesis scenario
                for n = 1 : length(virtual_sources)
                    obj.wfs_renderer{n} = wfs_renderer(virtual_sources{n}, binaural_sources, setup.Input_stream.SampleRate);
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
               output = 0;
               for n = 1 : length(obj.binaural_renderer)
                   obj.binaural_renderer{n}.binaural_source.set_input(input(:,n));
                   obj.binaural_renderer{n}.render;
                   output = output + obj.binaural_renderer{n}.output_signal;
               end
            %% Sound field synthesis job
            else
                wfs_output = 0;
                for m = 1 : length(obj.wfs_renderer)
                    obj.wfs_renderer{m}.virtual_source.set_input(input(:,m));
                    obj.wfs_renderer{m}.render;
                    wfs_output = wfs_output + cell2mat(obj.wfs_renderer{m}.output_signal);
                end
                output = 0;
                for n = 1 : length(obj.binaural_renderer)
                    obj.binaural_renderer{n}.binaural_source.set_input(wfs_output(:,n));
                    obj.binaural_renderer{n}.render;
                    output = output + obj.binaural_renderer{n}.output_signal;
                end
                
            end
        end
    end
end