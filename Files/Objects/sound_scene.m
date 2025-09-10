
classdef sound_scene < handle
    %SOUND_SCENE Summary of this class goes here
    %   Should contain description for the objects, present in the sound
    %   scene and interfaces positions with the gui

    properties
        virtual_sources % Virtual source: WFS or HOA or stereo etc source (should be a class with input signal,)
        loudspeaker_array  % SSD for SFS
        headphone
        receiver
        env
        scene_renderer
        downmixing_matrix
    end

    methods
        function obj = sound_scene(gui, setup, varargin)
            obj.env = environment;
            R0 = setup.Loudspeaker_setup.R;
            

            if isempty(gui.room)

                if isempty(varargin)
                    rec_pos = [0,0];
                    rec_orientation = [0,1];
                    [pos_vs,n_vs] = get_default_layout(setup.N_in,  R0 + 1, 'circular');
                else
                    rec_pos = varargin{1}.rec_position;
                    rec_orientation = varargin{1}.rec_orientation;
                    pos_vs = varargin{1}.vs_position;
                    n_vs = varargin{1}.vs_orientation;
                end
                [pos_ssd,n_ssd] = get_default_layout(setup.Loudspeaker_setup.N,R0, setup.Loudspeaker_setup.Shape);
            else
                if isempty(varargin)
                    rec_pos = [1,1];
                    rec_orientation = [0,1];
                else
                    rec_pos = varargin{1}.rec_position;
                    rec_orientation = varargin{1}.rec_orientation;
                end
                R0 = setup.Loudspeaker_setup.R;
                [pos_vs,n_vs] = get_default_layout(setup.N_in,  R0 + 1,'circular');
                pos_vs = pos_vs  + [1,1];
                [pos_ssd,n_ssd] = get_default_layout(setup.Loudspeaker_setup.N,R0*0.2,'circular');
                pos_ssd = pos_ssd + [1,1];
            end

            obj.create_receiver(gui, rec_pos, rec_orientation);
            for n = 1 : setup.N_in
                obj.create_virtual_source(pos_vs(n,:),n_vs(n,:), gui, setup.Virtual_source_type, setup.Rendering_mode);
                obj.virtual_sources{n}.set_input(zeros(setup.Block_size,1),setup.SampleRate);
            end
            switch setup.Rendering_mode
                case 'direct_playback'
                    for n = 1 : length(obj.virtual_sources)
                        obj.create_loudspeaker( obj.virtual_sources{n}.position,...
                            obj.virtual_sources{n}.orientation, ...
                            gui, setup.Loudspeaker_type, setup.Loudspeaker_setup.Height);
                        obj.loudspeaker_array{n}.set_output(zeros(setup.Block_size,1),setup.SampleRate);
                        if isvalid(gui.virtual_source_points{n}.UserData.text)
                            gui.loudspeaker_points{n}.UserData.text = gui.virtual_source_points{n}.UserData.text;
                        end
                        delete(gui.virtual_source_points{n});
                    end
                otherwise
                    for n = 1 : setup.Loudspeaker_setup.N
                        obj.create_loudspeaker( pos_ssd(n,:),n_ssd(n,:), gui, setup.Loudspeaker_type, setup.Loudspeaker_setup.Height);
                        obj.loudspeaker_array{n}.set_output(zeros(setup.Block_size,1),setup.SampleRate);
                    end
                    if isempty(gui.room)
                        gui.main_axes.XLim = 1.5*(R0+0.5)*[-1,1];
                        gui.main_axes.YLim = 1.5*(R0+0.5)*[-1,1];
                    end
            end
            obj.headphone = headphone(1,obj.receiver.position);
            obj.scene_renderer = sound_scene_renderer(obj.virtual_sources,obj.loudspeaker_array, obj.headphone, obj.receiver, setup);

            obj.downmixing_matrix = get_downmixing_mx(obj.loudspeaker_array,setup.N_out);
        end

        function obj = change_renderer(obj, setup)
            obj.scene_renderer = sound_scene_renderer(obj.virtual_sources,obj.loudspeaker_array, obj.headphone, obj.receiver, setup);
        end

        function obj = create_receiver(obj, gui, pos,orientation)
            R = 0.15;
            obj.receiver = receiver(pos,orientation);
            gui.receiver = gui.draw_head(pos,orientation,R);

            draggable(gui.receiver,@update_receiver_position, @update_receiver_orientation);
            function update_receiver_position(receiver)
                obj.receiver.set_position( gui.receiver.UserData.Origin );
                obj.scene_renderer.update_binaural_renderers('receiver_moved',[]);
                obj.scene_renderer.update_SFS_renderers('receiver_moved');
            end
            function update_receiver_orientation(receiver)
                obj.receiver.set_orientation(  [cosd(gui.receiver.UserData.Orientation),...
                    sind(gui.receiver.UserData.Orientation)] );
                obj.scene_renderer.update_binaural_renderers('receiver_rotated',[]);
                obj.scene_renderer.update_SFS_renderers('receiver_rotated');
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
                obj.virtual_sources{virtual_source.UserData.Label}.set_position(...
                    gui.virtual_source_points{virtual_source.UserData.Label}.UserData.Origin);
                obj.scene_renderer.update_SFS_renderers('virtual_source_moved');
                gui.request_draw_visuals(obj.scene_renderer.SFS_renderer{virtual_source.UserData.Label});
            end
            function update_virtual_orientation(virtual_source)
                obj.scene_renderer.update_SFS_renderers('virtual_source_rotated');
            end
        end

        function obj = delete_virtual_source(obj, virt_source_idx, gui)
            obj.virtual_sources{virt_source_idx} = {};
            delete(gui.virtual_source_points{virt_source_idx});
        end

        function obj = create_loudspeaker(obj, position, orientation, gui, type, height)
            idx = length(obj.loudspeaker_array) + 1;
            obj.loudspeaker_array{idx} = loudspeaker(idx, position, orientation, type, height);
            gui.loudspeaker_points{idx} = ...
                gui.draw_loudspeaker(position,type.R(1),cart2pol(orientation(1),orientation(2))*180/pi,idx);

            draggable(gui.loudspeaker_points{idx},@update_loudspeaker_position, @update_loudspeaker_orientation);
            function update_loudspeaker_position(loudspeaker)
                obj.loudspeaker_array{loudspeaker.UserData.Label}.set_position...
                    (gui.loudspeaker_points{loudspeaker.UserData.Label}.UserData.Origin);
                obj.scene_renderer.update_SFS_renderers('loudspeaker_moved');
                obj.scene_renderer.update_binaural_renderers('loudspeaker_moved',loudspeaker.UserData.Label);
            end
            function update_loudspeaker_orientation(loudspeaker)
                obj.loudspeaker_array{loudspeaker.UserData.Label}.set_orientation...
                    ( [ cosd(gui.loudspeaker_points{loudspeaker.UserData.Label}.UserData.Orientation),...
                    sind(gui.loudspeaker_points{loudspeaker.UserData.Label}.UserData.Orientation)]);
                obj.scene_renderer.update_SFS_renderers('loudspeaker_rotated');
                obj.scene_renderer.update_binaural_renderers('loudspeaker_rotated',loudspeaker.UserData.Label);
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

        function setup = save_setup(obj)
            setup.vs_position = cell2mat(cellfun(@(x) x.position, obj.virtual_sources, 'UniformOutput',0)');
            setup.vs_orientation = cell2mat(cellfun(@(x) x.orientation, obj.virtual_sources, 'UniformOutput',0)');
            setup.rec_position = obj.receiver.position;
            setup.rec_orientation = obj.receiver.orientation;
            setup.ls_position = cell2mat(cellfun(@(x) x.position, obj.loudspeaker_array, 'UniformOutput',0)');
            setup.ls_orientation = cell2mat(cellfun(@(x) x.orientation, obj.loudspeaker_array, 'UniformOutput',0)');



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
            for n = 1 : length(gui.main_axes.Children)
                delete(gui.main_axes.Children(1))
            end
            clear obj
        end

        function output = render_sound_scene(obj,input, Binauralization, Downmixing_enabled, Decorrelation)
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