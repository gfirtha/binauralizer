classdef sound_scene < handle
    %SOUND_SCENE Summary of this class goes here
    %   Should contain description for the objects, present in the sound
    %   scene and interfaces positions with the gui

    properties
        virtual_sources % Virtual source: WFS or HOA or stereo etc source (should be a class with input signal,)
        loudspeaker_array  % SSD for SFS
        headphone
        receiver
        environment
        scene_renderer
        downmixing_matrix
    end

    methods
        function obj = sound_scene(gui, setup)
            obj.environment = environment;
            obj.create_receiver(gui);

            R0 = setup.Loudspeaker_setup.R;
            pos = get_default_layout(setup.N_in,  R0 + 0.5);
            for n = 1 : setup.N_in
                obj.create_virtual_source(pos(n,:),-pos(n,:)/norm(pos(n,:)), gui, setup.Virtual_source_type, setup.Rendering_mode);
                obj.virtual_sources{n}.set_input(zeros(setup.Block_size,1),setup.SampleRate);
            end
            switch setup.Rendering_mode
                case 'Direct_playback'
                    for n = 1 : length(obj.virtual_sources)
                        obj.create_loudspeaker( obj.virtual_sources{n}.position,...
                            obj.virtual_sources{n}.orientation, ...
                            gui, setup.HRTF, setup.Loudspeaker_type);
                        obj.loudspeaker_array{n}.set_output(zeros(setup.Block_size,1),setup.SampleRate);
                        delete(gui.virtual_source_points{n});
                    end
                otherwise
                    pos_ssd = get_default_layout(setup.Loudspeaker_setup.N,R0);
                    for n = 1 : setup.Loudspeaker_setup.N
                        obj.create_loudspeaker( pos_ssd(n,:),-pos_ssd(n,:)/norm(pos_ssd(n,:)), gui, setup.HRTF, setup.Loudspeaker_type);
                        obj.loudspeaker_array{n}.set_output(zeros(setup.Block_size,1),setup.SampleRate);
                    end
                    gui.main_axes.XLim = (R0+1)*[-1,1];
                    gui.main_axes.YLim = (R0+1)*[-1,1];
            end
            obj.headphone = headphone(1,obj.receiver.position);
            obj.scene_renderer = sound_scene_renderer(obj.virtual_sources,obj.loudspeaker_array, obj.headphone, obj.receiver, setup);

            obj.downmixing_matrix = get_downmixing_mx(obj.loudspeaker_array,setup.N_out);
        end


        function obj = create_receiver(obj, gui)
            pos = [0,0];
            R = 0.15;
            obj.receiver = receiver(pos,[1,0]);
            gui.receiver = gui.draw_head(pos,R);
            draggable(gui.receiver,@update_receiver_position, @update_receiver_orientation);
            function update_receiver_position(receiver)
                obj.receiver.position = gui.receiver.UserData.Origin;
                obj.scene_renderer.update_binaural_renderers('receiver_moved');
            end
            function update_receiver_orientation(receiver)
                obj.receiver.orientation = [cosd(gui.receiver.UserData.Orientation),...
                    sind(gui.receiver.UserData.Orientation)];
                obj.scene_renderer.update_binaural_renderers('receiver_rotated' );
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
                obj.scene_renderer.update_SFS_renderers('virtual_source_moved');
            end
            function update_virtual_orientation(virtual_source)
                obj.scene_renderer.update_SFS_renderers('virtual_source_rotated');
            end
        end

        function obj = delete_virtual_source(obj, virt_source_idx, gui)
            obj.virtual_sources{virt_source_idx} = {};
            delete(gui.virtual_source_points{virt_source_idx});
        end

        function obj = create_loudspeaker(obj, position, orientation, gui, hrtf, type)
            idx = length(obj.loudspeaker_array) + 1;
            obj.loudspeaker_array{idx} = loudspeaker(idx, position, orientation, hrtf, type);
            gui.loudspeaker_points{idx} = ...
                gui.draw_loudspeaker(position,type.R(1),cart2pol(orientation(1),orientation(2))*180/pi,idx);

            draggable(gui.loudspeaker_points{idx},@update_loudspeaker_position, @update_loudspeaker_orientation);
            function update_loudspeaker_position(loudspeaker)
                obj.loudspeaker_array{loudspeaker.UserData.Label}.position...
                    = gui.loudspeaker_points{loudspeaker.UserData.Label}.UserData.Origin;
                obj.scene_renderer.update_SFS_renderers('loudspeaker_moved');
                obj.scene_renderer.update_binaural_renderers('loudspeaker_moved');
            end
            function update_loudspeaker_orientation(loudspeaker)
                obj.loudspeaker_array{loudspeaker.UserData.Label}.orientation...
                    = [ cosd(gui.loudspeaker_points{loudspeaker.UserData.Label}.UserData.Orientation),...
                    sind(gui.loudspeaker_points{loudspeaker.UserData.Label}.UserData.Orientation)];
                obj.scene_renderer.update_SFS_renderers('loudspeaker_rotated');
                obj.scene_renderer.update_binaural_renderers('loudspeaker_rotated');
            end
        end

        function obj = delete_loudspeaker(obj, bin_source_idx, gui)
            obj.loudspeaker_array(bin_source_idx) = [];
            delete(gui.loudspeaker_points{bin_source_idx});
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
            N = length(obj.loudspeaker_array);
            for n = 1 : N
                obj.delete_loudspeaker(N-n+1,gui);
            end
            clear obj
        end

        function output = render_sound_scene(obj,input, Binauralization, Downmixing_enabled)
            % Initizialize loudspeaker and headphone signals with zeros
            cellfun( @(x) x.source_signal.clear_signal, obj.loudspeaker_array);
            obj.headphone.source_signal.clear_signal;

            % Render
            obj.scene_renderer.render(input, Binauralization);
            if ~Binauralization
                if ~Downmixing_enabled
                    output = cell2mat(cellfun( @(x) x.source_signal.get_signal, obj.loudspeaker_array, 'UniformOutput',false));
                elseif Downmixing_enabled
                    output = cell2mat(cellfun( @(x) x.source_signal.get_signal, obj.loudspeaker_array, 'UniformOutput',false))*obj.downmixing_matrix';
                end
            elseif Binauralization
                output = obj.headphone.source_signal.get_signal;
            end
        end
    end
end