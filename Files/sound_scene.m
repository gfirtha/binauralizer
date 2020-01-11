classdef sound_scene < handle
    %SOUND_SCENE Summary of this class goes here
    %   Should contain description for the objects, present in the sound
    %   scene and interfaces positions with the gui
    
    properties
        virtual_sources % Virtual source: WFS or HOA or stereo etc source (should be a class with input signal,)
        binaural_sources % sources to be binauralized by binaural method
        receiver
        environment
        scene_renderer
    end
    
    methods
        function obj = sound_scene(gui, setup)
            obj.environment = environment;
            obj.create_receiver(gui);
            switch setup.Rendering
                case 'Binaural'
                    pos = get_default_layout(setup.Input_stream.info.NumChannels, 2);
                    for n = 1 : setup.Input_stream.info.NumChannels
                        obj.create_binaural_source(pos(n,:),-pos(n,:), gui, setup.HRTF);
                        obj.binaural_sources{n}.set_input(zeros(setup.Block_size,1));
                    end
                otherwise
                    pos = get_default_layout(setup.Input_stream.info.NumChannels,  setup.renderer_setup.R + 0.5);
                    for n = 1 : setup.Input_stream.info.NumChannels
                        obj.create_virtual_source(pos(n,:),-pos(n,:), gui);
                        obj.virtual_sources{n}.set_input(zeros(setup.Block_size,1));
                    end
                    pos_ssd = get_default_layout(setup.renderer_setup.N,setup.renderer_setup.R);
                    for n = 1 : setup.renderer_setup.N
                        obj.create_binaural_source( pos_ssd(n,:),-pos_ssd(n,:), gui, setup.HRTF);
                        obj.binaural_sources{n}.set_input(zeros(setup.Block_size,1));
                    end
                    gui.axes.XLim = (setup.renderer_setup.R+1)*[-1,1];
                    gui.axes.YLim = (setup.renderer_setup.R+1)*[-1,1];
            end
            
            obj.scene_renderer = sound_scene_renderer(obj.virtual_sources,obj.binaural_sources,obj.receiver, setup);
            
        end
        
        function obj = create_receiver(obj, gui)
            pos = [0,0];
            R = 0.2;
            obj.receiver = receiver(pos,0);
            gui.receiver = gui.draw_head(pos,R);
            draggable(gui.receiver,@update_receiver_position, @update_receiver_orientation);
            function update_receiver_position(receiver)
                pos = [gui.receiver.UserData.Origin];
                obj.receiver.position = pos;
                for n = 1 : length(obj.scene_renderer.binaural_renderer)
                    obj.scene_renderer.update_binaural_renderers(n);
                end
            end
            function update_receiver_orientation(receiver)
                obj.receiver.orientation = [ gui.receiver.UserData.Orientation ];
                for n = 1 : length(obj.scene_renderer.binaural_renderer)
                    obj.scene_renderer.update_binaural_renderers(n);
                end
            end
            
        end
        
        function obj = create_virtual_source(obj, position, orientation, gui)
            idx = length(obj.virtual_sources) + 1;
            obj.virtual_sources{idx} = virtual_source(idx, position, orientation);
            gui.virtual_source_points{idx} = drawpoint(gui.axes,...
                'Position',position,'Color','red','Label',sprintf('%d',idx),'LabelVisible','Off');
            addlistener(gui.virtual_source_points{idx} ,'MovingROI',@allevents);
            function allevents(~,evt)
                obj.virtual_sources{str2double(evt.Source.Label)}.position = evt.CurrentPosition;
                obj.scene_renderer.update_wfs_renderers(str2double(evt.Source.Label));
            end
        end
        
        function obj = delete_virtual_source(obj, virt_source_idx, gui)
            obj.virtual_sources{virt_source_idx} = {};
            delete(gui.virtual_source_points{virt_source_idx});
        end
        
        function obj = create_binaural_source(obj, position, orientation, gui, hrtf)
            idx = length(obj.binaural_sources) + 1;
            obj.binaural_sources{idx} = binaural_source(idx, position, orientation, hrtf);
            gui.binaural_source_points{idx} = gui.draw_loudspeaker(position,0.05,...
                cart2pol(orientation(1),orientation(2))*180/pi,idx);
            
            draggable(gui.binaural_source_points{idx},@update_binaural_position, @update_binaural_orientation);
            function update_binaural_position(binaural_source)
                obj.binaural_sources{binaural_source.UserData.Label}.position...
                    = gui.binaural_source_points{binaural_source.UserData.Label}.UserData.Origin;
                for n = 1 : length(obj.scene_renderer.wfs_renderer)
                    obj.scene_renderer.update_wfs_renderers(n);
                end
                obj.scene_renderer.update_binaural_renderers(binaural_source.UserData.Label);
            end
            function update_binaural_orientation(binaural_source)
                sprintf('anyad')
            end
        end
        
        function obj = delete_binaural_source(obj, bin_source_idx, gui)
            obj.binaural_sources(bin_source_idx) = [];
            delete(gui.binaural_source_points{bin_source_idx});
        end
        
        function obj = delete(obj,gui)
            N = length(obj.virtual_sources);
            for n = 1 : N
                obj.delete_virtual_source(N-n+1,gui);
            end
            N = length(obj.binaural_sources);
            for n = 1 : N
                obj.delete_binaural_source(N-n+1,gui);
            end
            clear obj
        end
        
        function output = binauralize_sound_scene(obj,input)
            output = obj.scene_renderer.render(input);
        end
    end
end

