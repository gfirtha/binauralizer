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
                    pos = get_default_layout(setup.Input_stream.info.NumChannels, 1.5);
                    for n = 1 : setup.Input_stream.info.NumChannels
                        obj.create_binaural_source(pos(n,:),-pos(n,:)/norm(pos(n,:)), gui, setup.HRTF, setup.Binaural_source_type);
                        obj.binaural_sources{n}.set_input(zeros(setup.Block_size,1),setup.Input_stream.SampleRate);
                    end
                otherwise
                    %R0 = mean(setup.HRTF.SourcePosition(:,3));
                    R0 = setup.renderer_setup.R;
                    pos = get_default_layout(setup.Input_stream.info.NumChannels,  R0 + 0.5);
                    for n = 1 : setup.Input_stream.info.NumChannels
                        obj.create_virtual_source(pos(n,:),-pos(n,:)/norm(pos(n,:)), gui, setup.Virtual_source_type, setup.Rendering);
                        obj.virtual_sources{n}.set_input(zeros(setup.Block_size,1),setup.Input_stream.SampleRate);
                    end
                    pos_ssd = get_default_layout(setup.renderer_setup.N,R0);
                    for n = 1 : setup.renderer_setup.N
                        obj.create_binaural_source( pos_ssd(n,:),-pos_ssd(n,:)/norm(pos_ssd(n,:)), gui, setup.HRTF, setup.Binaural_source_type);
                        obj.binaural_sources{n}.set_input(zeros(setup.Block_size,1),setup.Input_stream.SampleRate);
                    end
                    gui.main_axes.XLim = (R0+1)*[-1,1];
                    gui.main_axes.YLim = (R0+1)*[-1,1];
            end
            obj.scene_renderer = sound_scene_renderer(obj.virtual_sources,obj.binaural_sources,obj.receiver, setup);

        end
        
        function obj = create_receiver(obj, gui)
            pos = [0,0];
            R = 0.15;
            obj.receiver = receiver(pos,[1,0]);
            gui.receiver = gui.draw_head(pos,R);
            draggable(gui.receiver,@update_receiver_position, @update_receiver_orientation);
            function update_receiver_position(receiver)
                obj.receiver.position = gui.receiver.UserData.Origin;
                for n = 1 : length(obj.scene_renderer.binaural_renderer)
                    obj.scene_renderer.update_binaural_renderers(n,'receiver_moved');
                end 
            end
            function update_receiver_orientation(receiver)
                obj.receiver.orientation = [cosd(gui.receiver.UserData.Orientation),...
                                            sind(gui.receiver.UserData.Orientation)];
                for n = 1 : length(obj.scene_renderer.binaural_renderer)
                    obj.scene_renderer.update_binaural_renderers(n, 'receiver_rotated' );
                end
            end
            
        end
        
        function obj = delete_receiver(obj, gui)
            obj.receiver = {};
            delete(gui.receiver);
        end
        
        function obj = create_virtual_source(obj, position, orientation, gui, source_type, rendering)
            idx = length(obj.virtual_sources) + 1;
            obj.virtual_sources{idx} = virtual_source(idx, position, orientation, source_type, rendering);
            gui.virtual_source_points{idx} = gui.draw_virtual_source(position,...
                cart2pol(orientation(1),orientation(2))*180/pi,idx);
            
            draggable(gui.virtual_source_points{idx},@update_virtual_position, @update_virtual_orientation);
            function update_virtual_position(virtual_source)
                obj.virtual_sources{virtual_source.UserData.Label}.position...
                    = gui.virtual_source_points{virtual_source.UserData.Label}.UserData.Origin;
                 obj.scene_renderer.update_SFS_renderers(virtual_source.UserData.Label);
            end
            function update_virtual_orientation(virtual_source)
                sprintf('not supported')
            end
        end
        
        function obj = delete_virtual_source(obj, virt_source_idx, gui)
            obj.virtual_sources{virt_source_idx} = {};
            delete(gui.virtual_source_points{virt_source_idx});
        end
        
        function obj = create_binaural_source(obj, position, orientation, gui, hrtf, type)
            idx = length(obj.binaural_sources) + 1;
            obj.binaural_sources{idx} = binaural_source(idx, position, orientation, hrtf, type);
            gui.binaural_source_points{idx} = ...
                gui.draw_loudspeaker(position,type.R(1),cart2pol(orientation(1),orientation(2))*180/pi,idx);
            
            draggable(gui.binaural_source_points{idx},@update_binaural_position, @update_binaural_orientation);
            function update_binaural_position(binaural_source)
                obj.binaural_sources{binaural_source.UserData.Label}.position...
                    = gui.binaural_source_points{binaural_source.UserData.Label}.UserData.Origin;
                for n = 1 : length(obj.scene_renderer.SFS_renderer)
                    obj.scene_renderer.update_SFS_renderers(n);
                end
                obj.scene_renderer.update_binaural_renderers(binaural_source.UserData.Label,'source_moved');
            end
            function update_binaural_orientation(binaural_source)
                obj.binaural_sources{binaural_source.UserData.Label}.orientation...
                    = [ cosd(gui.binaural_source_points{binaural_source.UserData.Label}.UserData.Orientation),...
                        sind(gui.binaural_source_points{binaural_source.UserData.Label}.UserData.Orientation)];
                for n = 1 : length(obj.scene_renderer.SFS_renderer)
                    obj.scene_renderer.update_SFS_renderers(n);
                end
                obj.scene_renderer.update_binaural_renderers(binaural_source.UserData.Label,'source_rotated');
            end
        end
        
        function obj = delete_binaural_source(obj, bin_source_idx, gui)
            obj.binaural_sources(bin_source_idx) = [];
            delete(gui.binaural_source_points{bin_source_idx});
        end

        % Reserved for multiple sound scene scenario (for comparing approaches)
        function obj = save_sound_scene(obj)
        end
        function obj = load_sound_scene(obj)
        end
            
        function duplum = duplicate_scene(obj)
            duplum = obj;
        end
        
        function obj = delete(obj,gui)
            obj.delete_receiver(gui);
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