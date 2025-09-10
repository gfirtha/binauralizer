classdef room_view < handle
    %LISTENER_SPACE_AXES Summary of this class goes here
    %   Detailed explanation goes here

    properties
        scene
        main_gui
        room
        current_window
        axis
        receiver
        virtual_sources
        speakers
    end

    methods
        function obj = room_view(varargin)
            obj.scene = varargin{1};
            obj.receiver = obj.scene.receiver;
            obj.main_gui = varargin{2};

            if length(varargin)>2
                obj.room = varargin{3};
                obj.current_window = obj.room.drawRoomMesh(0);
            else
                obj.current_window = figure;
                axis equal tight
                set(gca, 'Clipping', 'Off')
                light("Style","local","Position",[2,5 2.9]);
                lighting phong
                size = 30;
                patch(size*[1,-1,-1,1],size*[1,1, -1,-1],[0,0,0,0],[1 1 1]*1)
            end
            obj.axis = gca;
            obj.axis.Clipping = "off";

            % Get the current axes handle
            hold on
            cnt = 0;
            for n = 1 : length(obj.main_gui.virtual_source_points)
                if isvalid(obj.main_gui.virtual_source_points{1})
                    cnt = cnt+1;
                    obj.virtual_sources{cnt} = draw_virtual_source(obj, obj.scene.virtual_sources{n});
                end
            end

            % Loudspeaker
            for n = 1 : length(obj.scene.loudspeaker_array)
                obj.speakers{n} = draw_speaker(obj,obj.scene.loudspeaker_array{n});
            end

            view(3)
            z0 = 1.7;
            % Set up initial camera position and view
            camPos = [obj.receiver.position'; z0]; % Initial camera position (adjust as needed)
            camTarget = camPos + [obj.receiver.orientation'; 0]; % Camera target position (center of the scene)
            camUp = [0; 1; 0]; % Up vector
            camva(obj.axis, 75); % Set the camera view angle
            campos(obj.axis, camPos);
            camtarget(obj.axis, camTarget);

            % Set perspective projection
            camproj(obj.axis, 'perspective');

            % Set orthographic projection
            %camproj(ax, 'orthographic');


            % Set up keyboard callbacks for camera movement/rotation
            set(gcf, 'KeyPressFcn', @obj.keyPressed);

            addlistener(obj.scene.receiver, 'position_changed', @obj.onReceiverPositionChanged);
            for n = 1 : length(obj.virtual_sources)
                addlistener(obj.scene.virtual_sources{n}, 'position_changed', @obj.onVirtualPositionChanged);
            end
            for n = 1 : length(obj.speakers)
                addlistener(obj.scene.loudspeaker_array{n}, 'position_changed', @obj.onLoudspeakerPositionChanged);
                addlistener(obj.scene.loudspeaker_array{n}, 'height_changed', @obj.onLoudspeakerHeightChanged);
            end
        end
        function onReceiverPositionChanged(obj, src, ~)
            z0 = 1.7;
            % Set up initial camera position and view
            camPos = [obj.receiver.position'; z0]; % Initial camera position (adjust as needed)
            camTarget = camPos + [obj.receiver.orientation'; 0]; % Camera target position (center of the scene)
            camUp = [0; 1; 0]; % Up vector
            camva(obj.axis, 75); % Set the camera view angle
            campos(obj.axis, camPos);
            camtarget(obj.axis, camTarget);
        end

        function onVirtualPositionChanged(obj, src, ~)
            ix = find(cell2mat(cellfun( @(x) (x.ID==src.source_index) ,obj.virtual_sources,'UniformOutput',0)));
            dx = src.position - obj.virtual_sources{ix}.center(1:2);
            obj.virtual_sources{ix}.patches.Vertices = obj.virtual_sources{ix}.patches.Vertices + [dx,0];
            obj.virtual_sources{ix}.center(1:2) =  obj.virtual_sources{ix}.center(1:2) +dx;
        end
        function onLoudspeakerPositionChanged(obj, src, ~)
            ix = find(cell2mat(cellfun( @(x) (x.ID==src.source_index) ,obj.speakers,'UniformOutput',0)));
            dx = src.position - obj.speakers{ix}.center(1:2);
            M = get_rotation_mx( [obj.speakers{ix}.orientation,0]', [src.orientation,0]' );
            M(3,3) = 1;
            for n = 1 : length(obj.speakers{ix}.patches)

                obj.speakers{ix}.patches{n}.Vertices = (M*(obj.speakers{ix}.patches{n}.Vertices-obj.speakers{ix}.center)')' + obj.speakers{ix}.center;
                obj.speakers{ix}.patches{n}.Vertices = obj.speakers{ix}.patches{n}.Vertices + [dx,0];
            end
            obj.speakers{ix}.orientation = src.orientation;
            obj.speakers{ix}.center(1:2) =  obj.speakers{ix}.center(1:2) +dx;

        end
        function onLoudspeakerHeightChanged(obj, src, ~)
            ix = find(cell2mat(cellfun( @(x) (x.ID==src.source_index) ,obj.speakers,'UniformOutput',0)));
            h0 = obj.speakers{ix}.center(3);
            dz = src.height - h0;
            for n = 1 : length(obj.speakers{ix}.patches)
                obj.speakers{ix}.patches{n}.Vertices = obj.speakers{ix}.patches{n}.Vertices + [0,0,dz];
            end
            obj.speakers{ix}.center(3) =  obj.speakers{ix}.center(3) + dz;

        end

        function virtual_source = draw_virtual_source(obj, vs_handle)
            sphere = spherical_design('Design', 120, 'T', 0.1);
            sphere.translate_sphere([vs_handle.position, 1.7]);
            hold on
            t = trisurf(sphere.delauney_trias',sphere.x0(:,1),sphere.x0(:,2),sphere.x0(:,3),'FaceColor','red','EdgeColor','none');
            virtual_source = struct('ID',vs_handle.source_index,'center', mean(sphere.x0,1), 'patches',t);
        end
        function speaker_out = draw_speaker(obj,speaker_handle)
            % Loudspeaker
            fileIn = fullfile('Data/Speaker_small.stl');
            [faces, vertices] = stlread(fileIn);
            [faces, vertices] = simplifyGeometry(faces, vertices);
            M = get_rotation_mx( [1,0,0]', [speaker_handle.orientation,0 ]' );
            M(3,3) = 1;
            vertices = (vertices - mean(vertices,1))/50;
            Lx = max(vertices(:,1))-min(vertices(:,1));
            vertices = vertices/Lx*2*speaker_handle.source_type.R;
            vertices =  (M*vertices')' + [speaker_handle.position,speaker_handle.height];
            hold on
            for n = 1 : length(faces)
                face = faces{n};
                if length(face) > 20
                    color = [0,0,0];
                else
                    color = [1,1,1];
                end
                p{n} = patch(vertices(face,1),vertices(face,2),vertices(face,3),color,'EdgeColor','none');
            end
            speaker_out = struct('ID',speaker_handle.source_index,'center', mean(vertices,1),'orientation',speaker_handle.orientation, 'patches',[]);
            speaker_out(1).patches = p;

        end
        function keyPressed(obj,event,keyData)
            switch keyData.Key
                case 'uparrow' % Rotate up
                    obj.rotateCamera('up');
                case 'downarrow' % Rotate down
                    obj.rotateCamera('down');
                case 'leftarrow' % Rotate left
                    obj.rotateCamera('left');
                case 'rightarrow' % Rotate right
                    obj.rotateCamera('right');
                case 'f' % Zoom in
                    obj.zoomCamera('in');
                case 'g' % Zoom out
                    obj.zoomCamera('out');
                case 'd' % Pan left
                    obj.panCamera('left');
                case 'a' % Pan right
                    obj.panCamera('right');
                case 'w' % Move forward
                    obj.moveForward();
                case 's' % Move backward
                    obj.moveBackward();
                case 'q' % Pan up
                    obj.panCamera('up');
                case 'e' % Pan down
                    obj.panCamera('down');
                case 'z' % Reset view
                    obj.resetView();
                case 'x' % Set view to axis limits
                    obj.setViewToAxisLimits();
                case 'y' % Move left
                    moveLeft();
                case 'c' % Move right
                    moveRight();
            end
            ax = gca;
            receiver_position = campos(ax); % Get current camera position
            receiver_orientation = camtarget(ax)-receiver_position; % Get current camera target
            obj.receiver.position = receiver_position(1:2);
            obj.receiver.orientation = receiver_orientation(1:2);
            obj.main_gui.update_receiver( receiver_position , receiver_orientation );
            obj.scene.scene_renderer.update_binaural_renderers('receiver_moved',[]);
            obj.scene.scene_renderer.update_SFS_renderers('receiver_moved');
            %            obj.scene.receiver.position
            %            obj.scene.receiver.orientation

        end



        function moveForward(obj)
            ax = gca;
            camPos = campos(ax); % Get current camera position
            camTarget = camtarget(ax); % Get current camera target

            % Calculate the camera direction
            camDir = camTarget - camPos;
            camDir = camDir / norm(camDir);

            % Move the camera forward along the camera direction
            camPos = camPos + camDir * 0.1; % Adjust the multiplier as needed
            camTarget = camPos + camDir;
            % Update the camera position
            campos(ax, camPos);
            camtarget(ax, camTarget);
        end

        % Function to move the camera backward
        function moveBackward(obj)
            ax = gca;

            camPos = campos(ax); % Get current camera position
            camTarget = camtarget(ax); % Get current camera target

            % Calculate the camera direction
            camDir = camTarget - camPos;
            camDir = camDir / norm(camDir);

            % Move the camera backward along the camera direction
            camPos = camPos - camDir * 0.1; % Adjust the multiplier as needed
            camTarget = camPos + camDir;
            % Update the camera position
            campos(ax, camPos);
            camtarget(ax, camTarget);
        end

        % Function to move the camera left
        function moveLeft(obj)
            ax = gca;

            camPos = campos(ax); % Get current camera position
            camTarget = camtarget(ax); % Get current camera target
            camUp = camup(ax); % Get current camera up vector

            % Calculate the camera direction
            camDir = camTarget - camPos;
            camDir = camDir / norm(camDir);

            % Calculate the vector pointing to the right relative to the camera's view
            camRight = cross(camDir, camUp);

            % Move the camera left along the camRight vector
            camPos = camPos - camRight * 0.1; % Adjust the multiplier as needed

            % Update the camera position
            campos(ax, camPos);
        end

        % Function to move the camera right
        function moveRight(obj)
            ax = gca;

            camPos = campos(ax); % Get current camera position
            camTarget = camtarget(ax); % Get current camera target
            camUp = camup(ax); % Get current camera up vector

            % Calculate the camera direction
            camDir = camTarget - camPos;
            camDir = camDir / norm(camDir);

            % Calculate the vector pointing to the right relative to the camera's view
            camRight = cross(camDir, camUp);

            % Move the camera right along the camRight vector
            camPos = camPos + camRight * 0.1; % Adjust the multiplier as needed

            % Update the camera position
            campos(ax, camPos);
        end

        % Function to rotate the camera
        function rotateCamera(obj, direction)
            ax = gca;

            switch direction
                case 'up'
                    camorbit(ax, 0, 1, 'camera');
                case 'down'
                    camorbit(ax, 0, -1, 'camera');
                case 'left'
                    camorbit(ax, 1, 0, 'camera');
                case 'right'
                    camorbit(ax, -1, 0, 'camera');
            end
        end

        % Function to zoom the camera
        function zoomCamera(obj, direction)
            ax = gca;

            switch direction
                case 'in'
                    camzoom(ax, 1.05);
                case 'out'
                    camzoom(ax, 0.95);
            end
        end

        % Function to pan the camera
        function panCamera(obj, direction)
            ax = gca;

            res = 5;
            switch direction
                case 'left'
                    campan(ax, res, 0, 'camera');
                case 'right'
                    campan(ax, -res, 0, 'camera');
                case 'up'
                    campan(ax, 0, 1, 'camera');
                case 'down'
                    campan(ax, 0, -1, 'camera');
                case 'panup'
                    campan(ax, 0, 0.5, 'camera');
                case 'pandown'
                    campan(ax, 0, -0.5, 'camera');
            end
        end

        % Function to reset the view
        function resetView(obj)
            ax = gca;

            camtarget(ax, [2 0 1]);
            campos(ax, [2 2 1]);
            camup(ax, camUp);
        end

        % Function to set the view to axis limits
        function setViewToAxisLimits(obj)
            ax = gca;
            axis(ax, 'tight');
        end

    end
end

