classdef sound_scene < handle
    %SOUND_SCENE Summary of this class goes here
    %   Should contain description for the objects, present in the sound
    %   scene and interfaces positions with the gui
    
    properties
        
        virtual_sources % Virtual source: WFS or HOA or stereo etc source (should be a class with input signal,)
        binaural_sources % sources to be binauralized by binaural method
        receiver
        scene_renderer
        convolution_setup
        gui
        
    end
    
    methods
        function obj = sound_scene(source_signal, gui, setup)
            obj.receiver = [0,0];
            obj.gui = gui;
            switch setup.Rendering
                case 'Binaural'
                    for n = 1 : source_signal.info.NumChannels
                        pos = randn(1,2);
                        pos = pos/norm(pos);
                        obj.create_binaural_source(pos,-pos,zeros(setup.Block_size,1), obj.gui, setup.HRTF);
                    end                    
                otherwise
                    for n = 1 : source_signal.info.NumChannels
                        pos = randn(1,2);
                        pos = pos/norm(pos);
                        obj.create_virtual_source(pos*(setup.WFS_Setup.R+0.25),-pos,zeros(setup.Block_size,1), obj.gui);
                    end
                    N = setup.WFS_Setup.N;
                    R0 = setup.WFS_Setup.R;
                    fi = (0:2*pi/N:2*pi-2*pi/N);
                    for n = 1 : N
                        obj.create_binaural_source( R0*[cos(fi(n)),sin(fi(n))],-[cos(fi(n)),sin(fi(n))],...
                            zeros(setup.Block_size,1), obj.gui, setup.HRTF);
                    end
                    obj.gui.axes.XLim = (R0+1)*[-1,1];
                    obj.gui.axes.YLim = (R0+1)*[-1,1];
            end
            
            obj.scene_renderer = sound_scene_renderer(obj.virtual_sources,obj.binaural_sources,obj.receiver);
            
        end
        
        function obj = create_virtual_source(obj, position, orientation, input, gui)
            idx = length(obj.virtual_sources);
            obj.virtual_sources{idx + 1} = virtual_source(idx + 1, position, orientation, input);
            gui.virtual_source_points{idx + 1} = drawpoint(gui.axes,...
                'Position',position,'Color','red','Label',sprintf('%d',idx + 1),'LabelVisible','Off');
            addlistener(gui.virtual_source_points{idx + 1} ,'MovingROI',@allevents);
            function allevents(~,evt)
                obj.virtual_sources{str2double(evt.Source.Label)}.position = evt.CurrentPosition;
                obj.scene_renderer.update_wfs_renderers(str2double(evt.Source.Label));
            end
        end
        
        function obj = delete_virtual_source(obj, virt_source_idx)
            obj.virtual_sources{virt_source_idx} = {};
            delete(obj.gui.virtual_source_points{virt_source_idx});
        end
        
        function obj = create_binaural_source(obj, position, orientation, input, gui, hrtf)
            idx = length(obj.binaural_sources);
            obj.binaural_sources{idx + 1} = binaural_source(idx + 1, position, orientation, input, hrtf);
            gui.binaural_source_points{idx + 1} = drawpoint(gui.axes,...
                'Position',position,...
                'Label',sprintf('%d',idx + 1),'LabelVisible','Off');
            addlistener(gui.binaural_source_points{idx + 1} ,'MovingROI',@allevents);
            function allevents(~,evt)
                obj.binaural_sources{str2double(evt.Source.Label)}.position = evt.CurrentPosition;
                for n = 1 : length(obj.scene_renderer.wfs_renderer)
                    obj.scene_renderer.update_wfs_renderers(n);
                end
                obj.scene_renderer.update_binaural_renderers(str2double(evt.Source.Label));
            end
        end
        
        function obj = delete_binaural_source(obj, bin_source_idx)
            obj.binaural_sources(bin_source_idx) = [];
            delete(obj.gui.binaural_source_points{bin_source_idx});
        end
        
        function obj = delete(obj)
            N = length(obj.virtual_sources);
            for n = 1 : N
                obj.delete_virtual_source(N-n+1);
            end
            N = length(obj.binaural_sources);
            for n = 1 : N
                obj.delete_binaural_source(N-n+1);
            end
            clear obj
        end
        
        function output = binauralize_sound_scene(obj,input)
            output = obj.scene_renderer.render(input);
        end
    end
end

