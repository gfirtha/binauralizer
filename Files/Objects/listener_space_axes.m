classdef listener_space_axes < handle
    %LISTENER_SPACE_AXES (optimized layering + object reuse)
    % - Stable z-order via 3 hggroups: bg (zones), mid (amp/taper), fg (lines/labels)
    % - Reuse graphics; no per-frame deletes/uistack
    % - Batched dashed zone rays in a single line (NaN-separated)
    % - Less transparency/edges; OpenGL renderer + childorder sort

    properties
        main_axes
        room
        receiver
        virtual_source_points
        loudspeaker_points
        designed_speaker_points = {}   % shaded models of designed loudspeakers (sound_scene.speakers order)
        binaural_source_points
        renderer_illustration
        beamPolar = struct()

        parentFig
        zoomPanel
        btnMinus
        btnPlus
        zoomFactorIn  = 0.80;  % <1 zooms in
        zoomFactorOut = 1.25;  % >1 zooms out
        zoomFactor

        rendererSetting
        activeRendererKey = ''
        lastDrawTic = [];     % for throttling
        minDrawDt = 1/15;     % max 10 FPS for visuals

        simulatorHandle

        % Scene view: 'top' (edit x/y), 'rear' (behind the listener, edit
        % x/height), '3d' (oblique, mouse orbit) or 'listener' (first
        % person: perspective camera in the listener's head; mouse look and
        % the keys turn / walk the listener itself). Heights live in the
        % patches' z vertices; the side/3D views add screen-facing markers
        % and stems down to the ear plane (z = 0), built by refresh_view.
        ViewMode = 'top'
        Move3DFcn            % @(kind, idx, xyz), set by sound_scene
        ReceiverMoveFcn      % @(xy, yawDeg), set by sound_scene (first-person walk / turn)
        FpPitch = 0          % first person: camera pitch [deg] (the listener stays level)
        FpFov = 75           % first person: camera view angle [deg]
        FloorZ = -1.2        % first person: floor below the ear plane [m]
        tv = []              % virtual TV: body patch, screen surface, size (show_tv)
        ViewChangedFcn = []  % @(mode) after a view change (main window: menu, video output)
        marks = struct()
        ViewLims = struct()  % per view: {xlim, ylim, zlim, camera} when it was left

        % Hooks for editors drawing into the scene (trajectory designer):
        BackgroundPressHook = []   % @() -> true if it consumed the press
        ExtraSceneHandlesFcn = []  % @() -> graphics that are scene objects, not overlays
    end

    methods
        function obj = listener_space_axes(varargin)
            if ~isempty(varargin) && isa(varargin{1},'matlab.graphics.axis.Axes')
                ax = varargin{1};
                obj.main_axes = ax;
            else
                obj.main_axes = axes('Units','normalized','Position',[0.08 0.18 0.88 0.78]);
                ax = obj.main_axes;
            end

            axis(obj.main_axes,'equal'); axis(obj.main_axes,'tight');
            obj.zoomFactor = 1;
            obj.lastDrawTic = tic;
            obj.renderer_illustration = struct();
            grid(obj.main_axes,'on');
            set(obj.main_axes,'ButtonDownFcn',[]);   % don’t rely on gca
            xlabel(obj.main_axes, 'x -> [m]');
            ylabel(obj.main_axes, 'y -> [m]');
            zlabel(obj.main_axes, 'z -> [m]');
            title('Main view')

            % Fast renderer defaults
            set(ancestor(ax,'figure'),'Renderer','opengl','GraphicsSmoothing','off');
            set(ax,'SortMethod','childorder','Clipping','on');

            if numel(varargin) > 1
                obj.room = varargin{2};
                hold(obj.main_axes,'on');
                obj.room.drawRoom(obj.main_axes);
                view(obj.main_axes, 2);
                xlim(obj.main_axes,[min(obj.room.vertices(:,2))-0.25 , max(obj.room.vertices(:,2))+0.25]);
                ylim(obj.main_axes,[min(obj.room.vertices(:,3))-0.25 , max(obj.room.vertices(:,3))+0.25]);
            else
                xlim(obj.main_axes,[-4,4]);
                ylim(obj.main_axes,[-4,4]);
            end

            % ---- Zoom UI (bottom-right of this axes) ----
            ax = obj.main_axes;
            obj.parentFig = ancestor(ax,'figure');

            % Put the floating controls into the SAME parent as the axes
            parentContainer = ax.Parent;
            bgc = ancestor(parentContainer,'figure').Color;

            obj.zoomPanel = uipanel( ...
                'Parent', parentContainer, ...
                'Units','pixels','BorderType','none', ...
                'BackgroundColor', bgc, ...
                'Tag','ZoomPanel');

            obj.btnMinus = uicontrol('Parent',obj.zoomPanel, ...
                'Style','pushbutton','String','–','Units','pixels', ...
                'Position',[0 0 28 28], 'TooltipString','Zoom out', ...
                'FontWeight','bold', ...
                'Callback',@(~,~)obj.zoomAxes(obj.zoomFactorOut), ...
                'BackgroundColor',[0.94 0.94 0.94]);

            obj.btnPlus = uicontrol('Parent',obj.zoomPanel, ...
                'Style','pushbutton','String','+','Units','pixels', ...
                'Position',[28 0 28 28], 'TooltipString','Zoom in', ...
                'FontWeight','bold', ...
                'Callback',@(~,~)obj.zoomAxes(obj.zoomFactorIn), ...
                'BackgroundColor',[0.94 0.94 0.94]);

            % Reposition on resize of BOTH the figure and the axes' parent panel
            prevFigSCF = obj.parentFig.SizeChangedFcn;
            obj.parentFig.SizeChangedFcn = @(src,evt) obj.onFigureResized(prevFigSCF, src, evt);

            try
                parentPrevSCF = parentContainer.SizeChangedFcn;
                parentContainer.SizeChangedFcn = @(src,evt) obj.onParentResized(parentPrevSCF, src, evt);
            catch
            end
            try
                prevAxSCF = obj.main_axes.SizeChangedFcn;
                obj.main_axes.SizeChangedFcn = @(src,evt) obj.onAxesResized(prevAxSCF, src, evt);
            catch
            end

            obj.positionZoomPanel();   % initial placement
            obj.install_mouse_nav();

            % Overlay settings: initialized lazily per renderer from
            % renderer.defaultOverlaySettings() on first draw call.
            obj.rendererSetting = struct();
        end

        function onFigureResized(obj, prevFcn, src, evt)
            if ~isempty(prevFcn)
                try, feval(prevFcn, src, evt); catch, end
            end
            obj.positionZoomPanel();
        end

        function onAxesResized(obj, prevFcn, src, evt)
            if ~isempty(prevFcn)
                try, feval(prevFcn, src, evt); catch, end
            end
            obj.positionZoomPanel();
        end

        function delete(obj)
            if ~isempty(obj.zoomPanel) && isvalid(obj.zoomPanel), delete(obj.zoomPanel); end
        end

        function hAx = get_main_axes(obj)
            hAx = obj.main_axes;
        end

        function receiver = draw_head(obj,pos,orientation,R)
            % 3D listener model (scene_glyph.head): from above it is the
            % former head-from-above icon, from the side a head and shoulders
            G = scene_glyph.head(R);
            receiver = patch(obj.main_axes, 'Faces', G.F, ...
                'Vertices', G.V + [pos(1) pos(2) 0], 'FaceVertexCData', G.C, ...
                'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceLighting', 'none', ...
                'Tag', 'Receiver');
            receiver.UserData = struct( 'Label', 1,...
                'Origin', [ mean(receiver.Vertices(:,1)), mean(receiver.Vertices(:,2))  ],...
                'Orientation', 0, 'Z', 0 );
            rotate(receiver, [0,0,1], asind(orientation(2)), [receiver.UserData.Origin(1),receiver.UserData.Origin(2),0]);
            receiver.UserData.Orientation = asind(orientation(2));
            listener_space_axes.place_z(receiver, pos);
            obj.receiver = receiver;
        end

        function update_receiver( obj, receiver_position, receiver_orientation )
            if isempty(obj.receiver) || ~isvalid(obj.receiver), return; end
            v = [cosd(obj.receiver.UserData.Orientation), sind(obj.receiver.UserData.Orientation)];
            w = receiver_orientation(1:2);
            dfi = atan2(w(2)*v(1)-w(1)*v(2), w(1)*v(1)+w(2)*v(2));

            rotate(obj.receiver,[0 0 1], dfi*180/pi,...
                [obj.receiver.UserData.Origin(1),obj.receiver.UserData.Origin(2),0]);
            obj.receiver.UserData.Orientation = obj.receiver.UserData.Orientation + dfi*180/pi;

            dx = receiver_position(1:2) - obj.receiver.UserData.Origin;
            obj.receiver.UserData.Origin = obj.receiver.UserData.Origin + dx;
            obj.receiver.Vertices(:,1:2) = obj.receiver.Vertices(:,1:2) + dx;
            if strcmp(obj.ViewMode, 'listener'), obj.apply_fp_camera(); end
        end

        function loudspeaker = draw_loudspeaker(obj,pos,R,orientation,idx)
            % 3D schematic loudspeaker (scene_glyph.loudspeaker): cabinet +
            % cone, the former 2D symbol when seen from above
            G = scene_glyph.loudspeaker(R);
            loudspeaker = patch(obj.main_axes, 'Faces', G.F, ...
                'Vertices', G.V + [pos(1) pos(2) 0], 'FaceVertexCData', G.C, ...
                'FaceColor', 'flat', 'EdgeColor', 'none', 'FaceLighting', 'none', ...
                'Tag', sprintf('loudspeaker_%i',idx));
            loudspeaker.UserData = struct( 'Label', idx,...
                'Origin', [ mean(loudspeaker.Vertices(:,1)), mean(loudspeaker.Vertices(:,2))  ],...
                'Orientation', orientation, 'Z', 0 );
            rotate(loudspeaker,[0 0 1], orientation,...
                [loudspeaker.UserData.Origin(1),loudspeaker.UserData.Origin(2),0]);
            listener_space_axes.place_z(loudspeaker, pos);
        end

        function h = draw_designed_speaker(obj, model, pos, yaw, idx)
            % Shaded 3D model of a designed loudspeaker (designed_loudspeaker
            % Mesh, local frame: faces +x, origin = bounding-box centre) at
            % pos, turned by yaw [deg]. The drag bookkeeping uses the vertex
            % mean (UserData.Origin); UserData.Offset is the vertex mean in
            % the model frame, so origin = Origin - Rz(yaw) Offset.
            M = model.Mesh;
            pos = [pos(:)', zeros(1, 3 - numel(pos))];
            % lowest faces first: the top view draws without depth sorting
            [~, o] = sort(M.V(M.F(:,1),3) + M.V(M.F(:,2),3) + M.V(M.F(:,3),3));
            M.F = M.F(o,:); M.C = M.C(o,:);
            h = patch(obj.main_axes, 'Faces', M.F, 'Vertices', obj.speaker_vertices(M.V, pos, yaw), ...
                'FaceVertexCData', M.C, 'FaceColor', 'flat', 'EdgeColor', 'none', ...
                'FaceLighting', 'none', 'Tag', sprintf('designed_speaker_%i', idx));
            setappdata(h, 'BaseV', M.V);
            h.UserData = struct('Label', idx, 'Origin', mean(h.Vertices(:,1:2), 1), ...
                'Orientation', yaw, 'Z', pos(3), 'Offset', mean(M.V(:,1:2), 1), 'DisplayName', model.Name);
            h.UserData.text = text(obj.main_axes, h.UserData.Origin(1) + 0.1, h.UserData.Origin(2) + 0.1, ...
                pos(3) + max(M.V(:,3)), model.Name, 'Interpreter', 'none');
            addlistener(h, 'ObjectBeingDestroyed', @(src,~) delete(src.UserData.text(isgraphics(src.UserData.text))));
            obj.designed_speaker_points{idx} = h;
        end

        function place_designed_speaker(obj, idx, pos, yaw)
            % Programmatic move / turn of a designed loudspeaker model
            h = obj.designed_speaker_points{idx};
            if ~isgraphics(h), return; end
            pos = [pos(:)', zeros(1, 3 - numel(pos))];
            h.Vertices = obj.speaker_vertices(getappdata(h, 'BaseV'), pos, yaw);
            h.UserData.Origin = mean(h.Vertices(:,1:2), 1);
            h.UserData.Orientation = yaw; h.UserData.Z = pos(3);
            if isgraphics(h.UserData.text)
                h.UserData.text.Position(1:2) = h.UserData.Origin + 0.1;
            end
        end

        function remove_designed_speaker(obj, idx)
            h = obj.designed_speaker_points{idx};
            if isgraphics(h), delete(h); end
            obj.designed_speaker_points(idx) = [];
            for k = idx:numel(obj.designed_speaker_points)       % labels follow the list
                if isgraphics(obj.designed_speaker_points{k})
                    obj.designed_speaker_points{k}.UserData.Label = k;
                end
            end
        end

        % --- Virtual TV (video playback in the 3D / first-person views) ---
        function show_tv(obj, tf)
            % The TV appears with a loaded video: first time 2.5 m in front
            % of the listener, facing it; then it keeps where it was put
            if tf && isempty(obj.tv)
                yaw = 90; xy = [0 0];
                if ~isempty(obj.receiver) && isgraphics(obj.receiver)
                    yaw = obj.receiver.UserData.Orientation; xy = obj.receiver.UserData.Origin;
                end
                obj.draw_tv(xy + 2.5*[cosd(yaw) sind(yaw)], yaw + 180, 65);
            end
            if isempty(obj.tv), return; end
            obj.tv.visible = tf;
            obj.ensure_tv();
            set([obj.tv.body obj.tv.screen], 'Visible', matlab.lang.OnOffSwitchState(tf));
        end

        function ensure_tv(obj)
            % a scene rebuild clears the axes: redraw the TV where it was
            if ~isempty(obj.tv) && ~isgraphics(obj.tv.body)
                obj.draw_tv(obj.tv.xy, obj.tv.yaw, obj.tv.diag, true);
            end
        end

        function draw_tv(obj, xy, yaw, diagInch, noRefresh)
            % TV body (draggable patch; double-click = screen size) and the
            % screen (texture-mapped surface that follows the body)
            old = obj.tv;
            if ~isempty(old), delete([old.body(isgraphics(old.body)) old.screen(isgraphics(old.screen))]); end
            d = diagInch*0.0254; w = d*16/sqrt(16^2 + 9^2); h = w*9/16;
            G = scene_glyph.tv(w, h, obj.FloorZ);
            body = patch(obj.main_axes, 'Faces', G.F, 'Vertices', obj.speaker_vertices(G.V, [xy 0], yaw), ...
                'FaceVertexCData', G.C, 'FaceColor', 'flat', 'EdgeColor', 'none', ...
                'FaceLighting', 'none', 'Tag', 'VirtualTV');
            body.UserData = struct('Label', 0, 'Origin', mean(body.Vertices(:,1:2), 1), ...
                'Orientation', yaw, 'Z', 0);
            img = zeros(9, 16, 3, 'uint8'); vis = true; locked = false;
            if ~isempty(old), img = old.img; vis = old.visible; locked = old.locked; end
            scr = surface(obj.main_axes, zeros(2), zeros(2), zeros(2), img, 'FaceColor', 'texturemap', ...
                'EdgeColor', 'none', 'CDataMapping', 'direct', 'FaceLighting', 'none', ...
                'HitTest', 'off', 'PickableParts', 'none', 'Tag', 'VirtualTVScreen');
            set([body scr], 'Visible', matlab.lang.OnOffSwitchState(vis));
            obj.tv = struct('body', body, 'screen', scr, 'w', w, 'h', h, 'diag', diagInch, ...
                'local', G.screen, 'locked', locked, 'visible', vis, 'img', img, 'xy', xy, 'yaw', yaw);
            draggable(body, @(~) obj.tv_follow(), @(~) obj.tv_follow());
            setappdata(body, 'user_dblclickfcn', @(~) obj.ask_tv_size());
            obj.tv_follow();
            obj.tv_interaction();
            if (nargin < 5 || ~noRefresh) && ~strcmp(obj.ViewMode, 'top'), obj.refresh_view(false); end
        end

        function tv_follow(obj)
            % screen corners from the body pose (row 1 = top; texturemap puts
            % image column 1 on the surface's last column, so column 1 = the
            % right edge as seen from the front)
            T = obj.tv; yaw = T.body.UserData.Orientation;
            obj.tv.xy = T.body.UserData.Origin; obj.tv.yaw = yaw;
            c = [T.body.UserData.Origin 0] + T.local*[cosd(yaw) sind(yaw) 0; -sind(yaw) cosd(yaw) 0; 0 0 1];
            r = [sind(yaw), -cosd(yaw), 0];          % the viewer's right
            L = c - r*T.w/2; R = c + r*T.w/2;
            set(T.screen, 'XData', [R(1) L(1); R(1) L(1)], 'YData', [R(2) L(2); R(2) L(2)], ...
                'ZData', T.h/2*[1 1; -1 -1]);
        end

        function ask_tv_size(obj)
            a = inputdlg('Screen diagonal [inch]:', 'Virtual TV', 1, {num2str(obj.tv.diag)});
            v = str2double(a);
            if isempty(a) || ~isfinite(v) || v < 10 || v > 200, return; end
            obj.draw_tv(obj.tv.body.UserData.Origin, obj.tv.body.UserData.Orientation, v);
        end

        function set_tv_image(obj, img)
            % img: still (idle) or current video frame
            if isempty(obj.tv), return; end
            obj.tv.img = img;
            obj.ensure_tv();
            obj.tv.screen.CData = img;
        end

        function tf = tv_on_screen(obj)
            % TV visible to the user: shown, and the scene tab is on top
            tf = ~isempty(obj.tv) && isgraphics(obj.tv.screen) && strcmp(obj.tv.screen.Visible, 'on');
            if tf
                tab = ancestor(obj.main_axes, 'uitab');
                tf = isempty(tab) || isequal(tab.Parent.SelectedTab, tab);
            end
        end

        function set_tv_locked(obj, tf)
            % during playback the TV stays where it is
            if isempty(obj.tv), return; end
            obj.tv.locked = tf;
            obj.tv_interaction();
        end

        function tv_interaction(obj)
            % top view and not playing: drag (left) / turn (right) / size
            % (double-click); otherwise clicks go to the view navigation
            if isempty(obj.tv) || ~isgraphics(obj.tv.body), return; end
            b = obj.tv.body;
            if ~isappdata(b, 'topButtonDownFcn'), setappdata(b, 'topButtonDownFcn', b.ButtonDownFcn); end
            if strcmp(obj.ViewMode, 'top') && ~obj.tv.locked
                b.ButtonDownFcn = getappdata(b, 'topButtonDownFcn');
            else
                b.ButtonDownFcn = @(~,~) obj.on_background_press();
            end
        end

        function set_object_hidden(obj, kind, k, tf)
            % Hide / show a virtual source ('vs') or loudspeaker ('ls') icon
            % (and its label, marker): drivers of a designed loudspeaker
            % are drawn by the model
            if strcmp(kind, 'vs'), L = obj.virtual_source_points; else, L = obj.loudspeaker_points; end
            if k > numel(L) || isempty(L{k}) || ~isgraphics(L{k}), return; end
            h = L{k};
            on = matlab.lang.OnOffSwitchState(~tf);
            if strcmp(kind, 'vs')
                setappdata(h, 'objVisible', ~tf);
                h.Visible = matlab.lang.OnOffSwitchState(~tf && strcmp(obj.ViewMode, 'top'));
                if isfield(obj.marks, 'vs') && k <= numel(obj.marks.vs) && isgraphics(obj.marks.vs(k).mk)
                    set([obj.marks.vs(k).mk, obj.marks.vs(k).stem], 'Visible', on);
                end
            else
                h.Visible = on;
            end
            if isfield(h.UserData, 'text') && isgraphics(h.UserData.text), h.UserData.text.Visible = on; end
        end

        function virtual_source = draw_virtual_source(obj,pos,orientation,idx,displayName)
            if nargin < 5 || isempty(displayName)
                displayName = sprintf('s_{%i}',idx);
            elseif isstring(displayName)
                displayName = char(displayName);
            end
            N = 20;
            R = [180,40,45,50];
            A = [0.03,0.08,0.12,0.16]*1.5;
            w = 0.002;
            x_source = A(1)*cosd(linspace(-R(1),R(1),N));
            y_source = A(1)*sind(linspace(-R(1),R(1),N));
            x_w1 = [A(2)*cosd(linspace(-R(2),R(2),N/2)) (A(2)-w)*cosd(linspace(R(2),-R(2),N/2))];
            y_w1 = [A(2)*sind(linspace(-R(2),R(2),N/2)) (A(2)-w)*sind(linspace(R(2),-R(2),N/2))];
            x_w2 = [A(3)*cosd(linspace(-R(3),R(3),N/2)) (A(3)-w)*cosd(linspace(R(3),-R(3),N/2))];
            y_w2 = [A(3)*sind(linspace(-R(3),R(3),N/2)) (A(3)-w)*sind(linspace(R(3),-R(3),N/2))];
            x_w3 = [A(4)*cosd(linspace(-R(4),R(4),N/2)) (A(4)-w)*cosd(linspace(R(4),-R(4),N/2))];
            y_w3 = [A(4)*sind(linspace(-R(4),R(4),N/2)) (A(4)-w)*sind(linspace(R(4),-R(4),N/2))];

            x = [x_source;x_w1;x_w2;x_w3]';
            y = [y_source;y_w1;y_w2;y_w3]';

            x = x - mean(mean(x));
            y = y - mean(mean(y));

            c = [255 0 0;
                0 0 0;
                0 0 0;
                0 0 0]/255;
            virtual_source = patch(obj.main_axes, x + pos(1), y + pos(2),[1;0;0;0],'Tag',sprintf('s_%i',idx));
            set(virtual_source,'FaceVertexCData',c,'FaceLighting','none');

            virtual_source.UserData = struct( 'Label', idx,...
                'Origin', [ mean(virtual_source.Vertices(:,1)), mean(virtual_source.Vertices(:,2))  ],...
                'Orientation', orientation, 'Z', 0, 'DisplayName', displayName );
            % unrotated icon around its origin (orient_virtual_source)
            virtual_source.UserData.Base = virtual_source.Vertices(:,1:2) - virtual_source.UserData.Origin;

            rotate(virtual_source,[0 0 1], orientation,[pos(1),pos(2),0]);
            listener_space_axes.place_z(virtual_source, pos);
            p = obj.patch_pos(virtual_source);
            virtual_source.UserData.text = text(obj.main_axes, pos(1)+0.1,pos(2)+0.1,p(3),displayName);
            addlistener(virtual_source,'ObjectBeingDestroyed', @(src,~) cleanupPatch(src));

            function cleanupPatch(src)
                if isfield(src.UserData,'text')
                    h = src.UserData.text;
                    if isgraphics(h) && isvalid(h)
                        delete(h);
                    end
                end
            end
        end

        function positionZoomPanel(obj)
            if isempty(obj.zoomPanel) || ~isvalid(obj.zoomPanel) || ...
                    isempty(obj.main_axes) || ~isvalid(obj.main_axes)
                return;
            end

            % compute plot box in PIXELS (accounts for axis equal / aspect)
            pb = obj.getPlotBoxPixels();

            margin = 1; panelW = 56; panelH = 28;
            obj.zoomPanel.Position = [pb(1) + pb(3) - panelW - margin, ...
                pb(2) + margin, ...
                panelW, panelH];
            uistack(obj.zoomPanel,'top');
        end
        function pb = getPlotBoxPixels(obj)
            % Returns [x y w h] of the actual data plot box in PIXELS,
            % i.e., the square/rect that respects aspect ratio (axis equal).
            ax = obj.main_axes;

            % Axes inner rectangle (in pixels)
            oldU = ax.Units; ax.Units = 'pixels';
            pos  = ax.Position;                % axes box (including labels)
            ti   = ax.TightInset;              % label/tick margins
            ax.Units = oldU;

            inner = [pos(1)+ti(1), pos(2)+ti(2), ...
                pos(3)-ti(1)-ti(3), pos(4)-ti(2)-ti(4)];

            % If something is off, fall back to inner
            if any(~isfinite(inner)) || inner(3) <= 0 || inner(4) <= 0
                pb = inner; return;
            end

            if ~strcmp(obj.ViewMode, 'top')
                pb = inner; return;
            end

            % Data aspect ratio from current limits
            xl = xlim(ax); yl = ylim(ax);
            if any(~isfinite([xl yl])) || diff(yl) == 0
                pb = inner; return;
            end
            dataAR = diff(xl) / diff(yl);      % width / height ratio

            % Fit a rectangle of aspect dataAR inside the inner box
            w = inner(3); h = inner(4);
            if w/h > dataAR
                % too wide → letterbox horizontally
                w2 = h * dataAR;
                x  = inner(1) + (w - w2)/2;
                y  = inner(2);
                pb = [x, y, w2, h];
            else
                % too tall → letterbox vertically
                h2 = w / dataAR;
                x  = inner(1);
                y  = inner(2) + (h - h2)/2;
                pb = [x, y, w, h2];
            end
        end

        function zoomAxes(obj, factor, c)
            % Zoom by factor (<1 = in) about point c (default: view centre).
            % Top view: x/y limits; rear view: the screen-plane axes (z and
            % the horizontal axes facing the camera); 3D view: camera zoom.
            ax = obj.main_axes;
            obj.zoomFactor = obj.zoomFactor .* factor;
            if strcmp(obj.ViewMode, 'listener')     % first person: field of view
                obj.FpFov = min(110, max(15, obj.FpFov*factor));
                camva(ax, obj.FpFov);
                return;
            end
            if strcmp(obj.ViewMode, '3d')
                camzoom(ax, 1/factor);
                return;
            end
            lims = {xlim(ax), ylim(ax), zlim(ax)};
            if nargin < 3
                c = cellfun(@mean, lims);
            end
            axs = [1 2];
            if strcmp(obj.ViewMode, 'rear')
                [az, ~] = view(ax);
                r = [cosd(az), sind(az)];          % screen-horizontal direction
                axs = [find(abs(r) > 1e-3), 3];
            end
            for a = axs
                lims{a} = c(a) + (lims{a} - c(a))*factor;
            end
            xlim(ax, lims{1}); ylim(ax, lims{2});
            if any(axs == 3), zlim(ax, lims{3}); end
        end

        function request_draw_visuals(obj, renderer)
            if toc(obj.lastDrawTic) < obj.minDrawDt
                return; % skip this frame
            end
            obj.lastDrawTic = tic;
            obj.draw_renderer_properties(renderer);
            drawnow limitrate nocallbacks
        end

        function obj = draw_renderer_properties(obj, renderer)
            % Generic dispatch: ask the renderer to draw its own overlay.
            % Each renderer gets its OWN settings sub-struct, keyed by class name.
            % On first call, settings are initialized from renderer.defaultOverlaySettings().
            if ~ismethod(renderer, 'drawOverlay'), return; end

            % Get renderer-specific settings (lazy init)
            rKey = class(renderer);  % e.g. 'wfs_renderer'
            if ~strcmp(obj.activeRendererKey, rKey)
                obj.clearRendererIllustration();
                obj.activeRendererKey = rKey;
            end

            L = obj.ensureLayers();

            if ~isfield(obj.rendererSetting, rKey)
                if ismethod(renderer, 'defaultOverlaySettings')
                    obj.rendererSetting.(rKey) = renderer.defaultOverlaySettings();
                else
                    obj.rendererSetting.(rKey) = struct();
                end
            end

            obj.renderer_illustration = renderer.drawOverlay( ...
                obj.main_axes, L, obj.rendererSetting.(rKey), obj.renderer_illustration);

            % Ensure scene objects (receiver, LS, VS) stay on top of overlays
            obj.raiseSceneObjects();
            obj.passthrough_overlays();
            if ~strcmp(obj.ViewMode, 'top')
                obj.stack_overlays();       % keep the layering under depth sorting
            end

            if isprop(obj,'simulatorHandle') && ~isempty(obj.simulatorHandle) && isvalid(obj.simulatorHandle)
                obj.simulatorHandle.refreshRendererCopies();
            end
        end

        function update_beam_polar(obj, anglesDeg, pattern, musicPattern, phantomPattern)
            % Live horizontal direction-energy overlay around the listener.
            if isempty(obj.main_axes) || ~isgraphics(obj.main_axes)
                return;
            end
            if nargin < 4
                musicPattern = [];
            end
            if nargin < 5
                phantomPattern = [];
            end
            hasDas = nargin >= 3 && ~isempty(pattern) && max(pattern(:)) > 0;
            hasMusic = ~isempty(musicPattern) && max(musicPattern(:)) > 0;
            hasPhantom = ~isempty(phantomPattern) && max(phantomPattern(:)) > 0;
            if ~hasDas && ~hasMusic && ~hasPhantom
                obj.set_beam_polar_visible(false);
                return;
            end

            anglesDeg = anglesDeg(:).';
            if hasDas, pattern = pattern(:).'; else, pattern = zeros(size(anglesDeg)); end
            if hasMusic, musicPattern = musicPattern(:).'; else, musicPattern = zeros(size(anglesDeg)); end
            if hasPhantom, phantomPattern = phantomPattern(:).'; else, phantomPattern = zeros(size(anglesDeg)); end
            n = min([numel(anglesDeg), numel(pattern), numel(musicPattern), numel(phantomPattern)]);
            anglesDeg = anglesDeg(1:n);
            pattern = max(0, pattern(1:n));
            musicPattern = max(0, musicPattern(1:n));
            phantomPattern = max(0, phantomPattern(1:n));
            if n < 3
                obj.set_beam_polar_visible(false);
                return;
            end

            L = obj.ensureLayers();
            ax = obj.main_axes;
            if ~isfield(obj.beamPolar, 'line') || ~isgraphics(obj.beamPolar.line)
                obj.beamPolar.fill = patch('Parent', L.mid, ...
                    'XData', nan, 'YData', nan, 'ZData', nan, ...
                    'FaceColor', [0.05 0.45 0.95], 'FaceAlpha', 0.12, ...
                    'EdgeColor', 'none', 'HitTest', 'off', 'PickableParts', 'none', ...
                    'Tag', 'direct_beam_polar_fill');
                obj.beamPolar.line = line('Parent', L.fg, ...
                    'XData', nan, 'YData', nan, 'ZData', nan, ...
                    'Color', [0.05 0.45 0.95], 'LineWidth', 1.7, ...
                    'HitTest', 'off', 'PickableParts', 'none', ...
                    'Tag', 'direct_beam_polar_line');
                obj.beamPolar.peak = line('Parent', L.fg, ...
                    'XData', nan, 'YData', nan, 'ZData', nan, ...
                    'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 7, ...
                    'MarkerFaceColor', [1 0.85 0.1], 'MarkerEdgeColor', [0.1 0.1 0.1], ...
                    'HitTest', 'off', 'PickableParts', 'none', ...
                    'Tag', 'direct_beam_polar_peak');
                obj.beamPolar.musicLine = line('Parent', L.fg, ...
                    'XData', nan, 'YData', nan, 'ZData', nan, ...
                    'Color', [0.9 0.05 0.85], 'LineWidth', 2.0, ...
                    'HitTest', 'off', 'PickableParts', 'none', ...
                    'Tag', 'direct_music_polar_line');
                obj.beamPolar.musicPeak = line('Parent', L.fg, ...
                    'XData', nan, 'YData', nan, 'ZData', nan, ...
                    'LineStyle', 'none', 'Marker', 'x', 'MarkerSize', 9, ...
                    'LineWidth', 1.8, 'Color', [0.9 0.05 0.85], ...
                    'HitTest', 'off', 'PickableParts', 'none', ...
                    'Tag', 'direct_music_polar_peak');
                obj.beamPolar.phantomLine = line('Parent', L.fg, ...
                    'XData', nan, 'YData', nan, 'ZData', nan, ...
                    'Color', [0.05 0.75 0.25], 'LineWidth', 2.1, ...
                    'HitTest', 'off', 'PickableParts', 'none', ...
                    'Tag', 'direct_phantom_polar_line');
                obj.beamPolar.phantomPeak = line('Parent', L.fg, ...
                    'XData', nan, 'YData', nan, 'ZData', nan, ...
                    'LineStyle', 'none', 'Marker', 'd', 'MarkerSize', 8, ...
                    'MarkerFaceColor', [0.05 0.75 0.25], 'MarkerEdgeColor', [0.02 0.25 0.08], ...
                    'HitTest', 'off', 'PickableParts', 'none', ...
                    'Tag', 'direct_phantom_polar_peak');
            end
            if ~isfield(obj.beamPolar, 'musicLine') || ~isgraphics(obj.beamPolar.musicLine)
                obj.beamPolar.musicLine = line('Parent', L.fg, ...
                    'XData', nan, 'YData', nan, 'ZData', nan, ...
                    'Color', [0.9 0.05 0.85], 'LineWidth', 2.0, ...
                    'HitTest', 'off', 'PickableParts', 'none', ...
                    'Tag', 'direct_music_polar_line');
                obj.beamPolar.musicPeak = line('Parent', L.fg, ...
                    'XData', nan, 'YData', nan, 'ZData', nan, ...
                    'LineStyle', 'none', 'Marker', 'x', 'MarkerSize', 9, ...
                    'LineWidth', 1.8, 'Color', [0.9 0.05 0.85], ...
                    'HitTest', 'off', 'PickableParts', 'none', ...
                    'Tag', 'direct_music_polar_peak');
            end
            if ~isfield(obj.beamPolar, 'phantomLine') || ~isgraphics(obj.beamPolar.phantomLine)
                obj.beamPolar.phantomLine = line('Parent', L.fg, ...
                    'XData', nan, 'YData', nan, 'ZData', nan, ...
                    'Color', [0.05 0.75 0.25], 'LineWidth', 2.1, ...
                    'HitTest', 'off', 'PickableParts', 'none', ...
                    'Tag', 'direct_phantom_polar_line');
                obj.beamPolar.phantomPeak = line('Parent', L.fg, ...
                    'XData', nan, 'YData', nan, 'ZData', nan, ...
                    'LineStyle', 'none', 'Marker', 'd', 'MarkerSize', 8, ...
                    'MarkerFaceColor', [0.05 0.75 0.25], 'MarkerEdgeColor', [0.02 0.25 0.08], ...
                    'HitTest', 'off', 'PickableParts', 'none', ...
                    'Tag', 'direct_phantom_polar_peak');
            end

            baseR = 0.45;
            scaleR = 1.45;
            if hasDas
                r = baseR + scaleR*sqrt(pattern ./ (max(pattern) + eps));
                x = -r .* sind(anglesDeg);
                y =  r .* cosd(anglesDeg);
                z = 0.035 * ones(size(x));
                x = [x, x(1)];
                y = [y, y(1)];
                z = [z, z(1)];
                set(obj.beamPolar.fill, 'XData', x, 'YData', y, 'ZData', z, 'Visible', 'on');
                set(obj.beamPolar.line, 'XData', x, 'YData', y, 'ZData', z, 'Visible', 'on');
                [~, imax] = max(pattern);
                set(obj.beamPolar.peak, 'XData', x(imax), 'YData', y(imax), 'ZData', z(imax), 'Visible', 'on');
            else
                set([obj.beamPolar.fill, obj.beamPolar.line, obj.beamPolar.peak], 'Visible', 'off');
            end

            if hasMusic
                rM = baseR + scaleR*sqrt(musicPattern ./ (max(musicPattern) + eps));
                xM = -rM .* sind(anglesDeg);
                yM =  rM .* cosd(anglesDeg);
                zM = 0.07 * ones(size(xM));
                xM = [xM, xM(1)];
                yM = [yM, yM(1)];
                zM = [zM, zM(1)];
                set(obj.beamPolar.musicLine, 'XData', xM, 'YData', yM, 'ZData', zM, 'Visible', 'on');
                [~, imaxM] = max(musicPattern);
                set(obj.beamPolar.musicPeak, 'XData', xM(imaxM), 'YData', yM(imaxM), 'ZData', zM(imaxM), 'Visible', 'on');
            else
                set([obj.beamPolar.musicLine, obj.beamPolar.musicPeak], 'Visible', 'off');
            end

            if hasPhantom
                rP = baseR + scaleR*sqrt(phantomPattern ./ (max(phantomPattern) + eps));
                xP = -rP .* sind(anglesDeg);
                yP =  rP .* cosd(anglesDeg);
                zP = 0.105 * ones(size(xP));
                xP = [xP, xP(1)];
                yP = [yP, yP(1)];
                zP = [zP, zP(1)];
                set(obj.beamPolar.phantomLine, 'XData', xP, 'YData', yP, 'ZData', zP, 'Visible', 'on');
                [~, imaxP] = max(phantomPattern);
                set(obj.beamPolar.phantomPeak, 'XData', xP(imaxP), 'YData', yP(imaxP), 'ZData', zP(imaxP), 'Visible', 'on');
            else
                set([obj.beamPolar.phantomLine, obj.beamPolar.phantomPeak], 'Visible', 'off');
            end

            obj.raiseSceneObjects();
            if ~strcmp(obj.ViewMode, 'top')
                obj.stack_overlays();
            end
            if isprop(obj,'simulatorHandle') && ~isempty(obj.simulatorHandle) && isvalid(obj.simulatorHandle)
                obj.simulatorHandle.refreshRendererCopies();
            end
            obj.apply_sort_method();
        end

        function apply_sort_method(obj)
            % top view: layers by child order (flat overlays); side / 3D /
            % first-person views: depth test (else the TV, loudspeakers and
            % labels are drawn in list order, behind ones over the screen)
            if strcmp(obj.ViewMode, 'top')
                obj.main_axes.SortMethod = 'childorder';
            else
                obj.main_axes.SortMethod = 'depth';
            end
        end

        function set_beam_polar_visible(obj, tf)
            if isempty(obj.beamPolar) || ~isstruct(obj.beamPolar), return; end
            fns = fieldnames(obj.beamPolar);
            vis = matlab.lang.OnOffSwitchState(tf);
            for k = 1:numel(fns)
                h = obj.beamPolar.(fns{k});
                if isgraphics(h)
                    h.Visible = vis;
                end
            end
        end

        function move_virtual_source(obj, k, xyz)
            % Programmatic move (object-based input), same bookkeeping as a
            % drag. xyz: [x y] or [x y z]
            h = obj.virtual_source_points{k};
            if ~isgraphics(h), return; end
            p0 = obj.patch_pos(h);
            p = p0;
            p(1:numel(xyz)) = xyz;
            if isequal(p, p0), return; end
            obj.set_patch_pos(h, p);
            obj.update_vs_mark(k, p);
        end

        function orient_virtual_source(obj, k, ang)
            % Turn the source icon to face direction ang [deg]. Rebuilt from
            % the unrotated icon: draggable re-sets the click-time vertices
            % on every mouse move, so an incremental rotate would be lost.
            h = obj.virtual_source_points{k};
            if ~isgraphics(h) || ~isfield(h.UserData, 'Base'), return; end
            Rz = [cosd(ang) -sind(ang); sind(ang) cosd(ang)];
            p = obj.patch_pos(h);
            h.Vertices = [h.UserData.Base*Rz.' + p(1:2), p(3)*ones(size(h.UserData.Base,1),1)];
            h.UserData.Orientation = ang;
        end

        function style_virtual_sources(obj, up, level)
            % Object-based input: dot colour = height (blue floor -> orange
            % ceiling). Shown only while an object really contributes: within
            % 30 dB of the loudest object and above -70 dBFS (level is
            % smoothed by the caller, ~0.5 s release). Opacity follows the
            % level within that 30 dB window. Silent objects - Atmos parks
            % them in the front-left corner - and lingering tails stay hidden.
            persistent cmap
            if isempty(cmap)
                cmap = interp1([0 0.5 1], [0.15 0.35 0.85; 0.55 0.55 0.55; 0.95 0.45 0.1], linspace(0,1,64));
            end
            dB = 20*log10(level(:) + eps);
            floorDB = max(-70, max(dB) - 30);
            L = min(max((dB - floorDB)/30, 0), 1);
            for k = 1:min(numel(up), numel(obj.virtual_source_points))
                h = obj.virtual_source_points{k};
                if ~isgraphics(h), continue; end
                c = h.FaceVertexCData;
                c(1,:) = cmap(1 + round(63*min(max(up(k),0),1)), :);
                h.FaceVertexCData = c;
                vis = dB(k) > floorDB;
                setappdata(h, 'objVisible', vis);
                h.Visible = matlab.lang.OnOffSwitchState(vis && strcmp(obj.ViewMode, 'top'));
                if isfield(h.UserData,'text') && isgraphics(h.UserData.text)
                    h.UserData.text.Visible = matlab.lang.OnOffSwitchState(vis);
                end
                h.FaceAlpha = 0.15 + 0.85*L(k);
                h.EdgeAlpha = 0.15 + 0.85*L(k);
                if isfield(obj.marks,'vs') && k <= numel(obj.marks.vs) && isgraphics(obj.marks.vs(k).mk)
                    set([obj.marks.vs(k).mk, obj.marks.vs(k).stem], 'Visible', matlab.lang.OnOffSwitchState(vis));
                    obj.marks.vs(k).mk.MarkerFaceColor = c(1,:);
                end
            end
        end

        % --- Scene view (top / rear / 3D) ---
        function set_view_mode(obj, mode)
            % Each view keeps its own camera and zoom/pan: remember the
            % limits and camera (position, target, view angle, up vector) of
            % the view we leave, restore those of the view we enter (side/3D
            % views start from the whole scene, not the top view's zoom)
            ax = obj.main_axes;
            mode = lower(mode);
            if ~strcmp(mode, obj.ViewMode) && isgraphics(ax)
                cam = {campos(ax), camtarget(ax), camva(ax), camup(ax)};
                obj.ViewLims.(obj.view_key(obj.ViewMode)) = {xlim(ax), ylim(ax), zlim(ax), cam};
            end
            obj.ViewMode = mode;
            if strcmp(mode, 'listener'), obj.install_key_nav(); end
            obj.refresh_view(true);
            if ~isempty(obj.ViewChangedFcn), obj.ViewChangedFcn(mode); end
        end

        function refresh_view(obj, restore)
            % Apply ViewMode: camera, markers, picking. Called on a mode
            % change and after the scene objects are (re)created.
            % Virtual sources: flat wave icon in the top view, a sphere
            % marker (+ stem to the ear plane) in the side/3D views.
            % Receiver and loudspeakers are 3D models in every view.
            if nargin < 2, restore = false; end
            ax = obj.main_axes;
            if isempty(ax) || ~isgraphics(ax), return; end
            obj.ensure_tv();
            obj.delete_marks();
            mode = obj.ViewMode;
            top = strcmp(mode, 'top');
            for k = 1:numel(obj.virtual_source_points)
                h = obj.virtual_source_points{k};
                if isempty(h) || ~isgraphics(h), continue; end
                vis = getappdata(h, 'objVisible');
                if isempty(vis), vis = true; end
                h.Visible = matlab.lang.OnOffSwitchState(vis && top);
            end
            % Picking: top = draggable patches (x/y); rear = loudspeaker
            % models and source markers drag in the screen plane; 3D =
            % clicks go to the axes (orbit)
            for k = 1:numel(obj.loudspeaker_points)
                h = obj.loudspeaker_points{k};
                if isempty(h) || ~isgraphics(h), continue; end
                if ~isappdata(h, 'topButtonDownFcn')
                    setappdata(h, 'topButtonDownFcn', h.ButtonDownFcn);
                end
                switch mode
                    case 'top'
                        h.ButtonDownFcn = getappdata(h, 'topButtonDownFcn');
                        h.PickableParts = 'visible';
                    case 'rear'
                        h.ButtonDownFcn = @(src,~) obj.start_side_drag('ls', src.UserData.Label);
                        h.PickableParts = 'visible';
                    otherwise
                        % 3D: double-click = test signal, otherwise orbit/pan
                        h.ButtonDownFcn = @(src,~) obj.on_press_3d(src);
                        h.PickableParts = 'visible';
                end
            end
            for k = 1:numel(obj.designed_speaker_points)
                h = obj.designed_speaker_points{k};
                if isempty(h) || ~isgraphics(h), continue; end
                if ~isappdata(h, 'topButtonDownFcn')
                    setappdata(h, 'topButtonDownFcn', h.ButtonDownFcn);
                end
                switch mode
                    case 'top',  h.ButtonDownFcn = getappdata(h, 'topButtonDownFcn');
                    case 'rear', h.ButtonDownFcn = @(src,~) obj.start_side_drag('ds', src.UserData.Label);
                    otherwise,   h.ButtonDownFcn = @(src,~) obj.on_press_3d(src);
                end
            end
            fp = strcmp(mode, 'listener');
            if ~isempty(obj.receiver) && isgraphics(obj.receiver)
                if ~top, obj.receiver.PickableParts = 'none'; else, obj.receiver.PickableParts = 'visible'; end
                obj.receiver.Visible = matlab.lang.OnOffSwitchState(~fp);   % first person: we are the head
            end
            ax.ButtonDownFcn = @(~,~) obj.on_background_press();
            obj.tv_interaction();
            obj.passthrough_overlays();
            obj.stack_overlays();
            % drop 3D camera zoom/pan before applying the new view
            set(ax, 'CameraPositionMode', 'auto', 'CameraTargetMode', 'auto', ...
                'CameraViewAngleMode', 'auto', 'CameraUpVectorMode', 'auto', 'Projection', 'orthographic');
            if top
                view(ax, 2);
                zlim(ax, 'auto');
                key = obj.view_key(mode);
                if restore && isfield(obj.ViewLims, key)
                    L = obj.ViewLims.(key);
                    xlim(ax, L{1}); ylim(ax, L{2});
                end
                ax.SortMethod = 'childorder';
                obj.positionZoomPanel();
                return;
            end

            key = obj.view_key(mode);
            if restore && isfield(obj.ViewLims, key)
                L = obj.ViewLims.(key);
                xlim(ax, L{1}); ylim(ax, L{2}); zlim(ax, L{3});
            elseif restore
                obj.fit_scene();                 % first visit: the whole scene
            else
                zlim(ax, obj.z_range());
            end
            obj.build_marks();                   % (ear plane spans the limits)
            ax.SortMethod = 'depth';
            % azimuth of the camera behind the listener (view az 0 looks
            % along +y; the head patch keeps its facing angle in degrees)
            az = 0;
            if ~isempty(obj.receiver) && isgraphics(obj.receiver)
                az = obj.receiver.UserData.Orientation - 90;
            end
            switch mode
                case 'rear'
                    view(ax, az, 0);
                case '3d'
                    view(ax, az + 37.5, 30);       % from behind-right of the listener
                case 'listener'
                    % the camera follows the listener, not a stored camera
                    zl = zlim(ax); zlim(ax, [min(zl(1), obj.FloorZ - 0.05), zl(2)]);
                    obj.apply_fp_camera();
                    obj.positionZoomPanel();
                    return;
            end
            if restore && isfield(obj.ViewLims, key) && numel(obj.ViewLims.(key)) >= 4
                % where this view's camera was left (orbit, pan, zoom)
                cam = obj.ViewLims.(key){4};
                set(ax, 'CameraPosition', cam{1}, 'CameraTarget', cam{2}, ...
                    'CameraViewAngle', cam{3}, 'CameraUpVector', cam{4});
            end
            obj.positionZoomPanel();
        end

        function stack_overlays(obj)
            % Renderer overlays are flat (z = 0) and layered by child order,
            % which only holds without depth sorting (top view). In the
            % side/3D views each overlay is lifted by a few mm per layer
            % (one hgtransform each, so the painters keep updating X/YData),
            % keeping the top-view stacking under depth sorting. Top view:
            % the overlays go back to the axes in the same order.
            ax = obj.main_axes;
            T = findobj(ax.Children, 'flat', 'Tag', 'OverlayLift');
            if ~isempty(T)
                % unwrap, lowest lift first: the last one lands on top
                dz = arrayfun(@(t) t.Matrix(3,4), T);
                [~, o] = sort(dz);
                for t = T(o)'
                    if ~isempty(t.Children), set(t.Children, 'Parent', ax); end
                end
                delete(T);
                obj.raiseSceneObjects();
            end
            if strcmp(obj.ViewMode, 'top'), return; end
            scene = obj.scene_handles();
            kids = ax.Children;                      % top-most first
            ov = kids(~ismember(kids, scene) & ...
                ~arrayfun(@(h) strcmp(h.Tag, 'EarPlane') || isa(h, 'matlab.graphics.primitive.Text'), kids));
            n = numel(ov);
            vis = ~any(strcmp(obj.ViewMode, {'rear', 'listener'}));   % edge-on in the rear / first-person view
            for k = 1:n
                t = hgtransform('Parent', ax, 'Tag', 'OverlayLift', ...
                    'HitTest', 'off', 'PickableParts', 'none', 'Visible', vis, ...
                    'Matrix', makehgtform('translate', [0 0 5e-4*(n - k + 1)]));
                ov(k).Parent = t;
            end
        end

        function h = scene_handles(obj)
            % Graphics of the scene objects (not renderer overlays)
            c = [obj.virtual_source_points(:)', obj.loudspeaker_points(:)', obj.designed_speaker_points(:)', {obj.receiver}];
            if ~isempty(obj.tv), c = [c, {obj.tv.body, obj.tv.screen}]; end
            c = c(cellfun(@(x) ~isempty(x) && isgraphics(x), c));
            h = [c{:}];
            if isempty(h), h = gobjects(0); end
            if ~isempty(obj.ExtraSceneHandlesFcn)
                try
                    e = obj.ExtraSceneHandlesFcn();
                    h = [h, reshape(e(isgraphics(e)), 1, [])];
                catch
                end
            end
            M = obj.marks;
            if isfield(M, 'vs'), h = [h, M.vs.mk, M.vs.stem]; end
            for f = {'lsStem','plane'}
                if isfield(M, f{1}) && ~isempty(M.(f{1})), h = [h, M.(f{1})]; end %#ok<AGROW>
            end
        end


        function raiseSceneObjects(obj)
            % Push receiver, loudspeaker and VS patches above overlay layers.
            % One reordering of the axes' children, and only when needed:
            % per-object uistack cost ~25 ms each (~10 s per redraw with 128
            % loudspeakers), and the overlay layers are persistent hggroups,
            % so once ordered the scene objects normally stay on top.
            try
                % top-most first: virtual sources, loudspeakers, receiver
                % (the order the former uistack sequence produced)
                scene = [obj.virtual_source_points(:)', obj.loudspeaker_points(:)', obj.designed_speaker_points(:)', {obj.receiver}];
                if ~isempty(obj.tv), scene = [scene, {obj.tv.screen, obj.tv.body}]; end
                % with the source labels
                lbl = cellfun(@(h) listener_space_axes.label_of(h), ...
                    [obj.virtual_source_points(:)', obj.designed_speaker_points(:)'], 'UniformOutput', false);
                scene = [lbl, scene];
                scene = scene(cellfun(@(h) ~isempty(h) && isgraphics(h) && ...
                    isequal(h.Parent, obj.main_axes), scene));
                if isempty(scene), return; end
                scene = [scene{:}]';
                kids = obj.main_axes.Children;
                nS = numel(scene);
                if numel(kids) >= nS && all(ismember(kids(1:nS), scene))
                    return;                       % already on top
                end
                obj.main_axes.Children = [scene; kids(~ismember(kids, scene))];
            catch
            end
        end

        function settings = get_settings(obj, renderer)
            % Get overlay settings for a specific renderer (by object or class name).
            % Without argument: returns the whole rendererSetting struct.
            if nargin < 2
                settings = obj.rendererSetting;
            else
                rKey = class(renderer);
                if isfield(obj.rendererSetting, rKey)
                    settings = obj.rendererSetting.(rKey);
                else
                    settings = struct();
                end
            end
        end

        function obj = set_settings(obj, renderer, setting)
            % Set overlay settings for a specific renderer.
            % Also accepts old 1-arg form: set_settings(settingStruct)
            if nargin == 2 && isstruct(renderer)
                % Legacy: set_settings(fullStruct)
                obj.rendererSetting = renderer;
            else
                rKey = class(renderer);
                obj.rendererSetting.(rKey) = setting;
            end
        end

        function obj = clearObject(obj)
            if ~isempty(obj.receiver) && isvalid(obj.receiver), delete(obj.receiver); end
            obj.receiver = [];

            if ~isempty(obj.virtual_source_points)
                for n = 1:numel(obj.virtual_source_points)
                    if isgraphics(obj.virtual_source_points{n}), delete(obj.virtual_source_points{n}); end
                end
            end
            obj.virtual_source_points = {};

            if ~isempty(obj.loudspeaker_points)
                for n = 1:numel(obj.loudspeaker_points)
                    if isgraphics(obj.loudspeaker_points{n}), delete(obj.loudspeaker_points{n}); end
                end
            end
            obj.loudspeaker_points = {};

            for n = 1:numel(obj.designed_speaker_points)
                if isgraphics(obj.designed_speaker_points{n}), delete(obj.designed_speaker_points{n}); end
            end
            obj.designed_speaker_points = {};

            if ~isempty(obj.binaural_source_points)
                for n = 1:numel(obj.binaural_source_points)
                    if isgraphics(obj.binaural_source_points{n}), delete(obj.binaural_source_points{n}); end
                end
            end
            obj.binaural_source_points = {};

            obj.clearRendererIllustration();
        end

        function onParentResized(obj, prevFcn, src, evt)
            if ~isempty(prevFcn)
                try, feval(prevFcn, src, evt); catch, end
            end
            obj.positionZoomPanel();
        end

        function reattachZoomUI(obj)
            % Call this after you move the axes to a new parent panel
            if isempty(obj.zoomPanel) || ~isvalid(obj.zoomPanel), return; end
            newParent = obj.main_axes.Parent;
            if obj.zoomPanel.Parent ~= newParent
                obj.zoomPanel.Parent = newParent;
            end
            % rehook parent panel resize callback
            try
                parentPrevSCF = newParent.SizeChangedFcn;
                newParent.SizeChangedFcn = @(src,evt) obj.onParentResized(parentPrevSCF, src, evt);
            catch
            end
            obj.positionZoomPanel();
            uistack(obj.zoomPanel,'top');
            obj.install_mouse_nav();        % the axes moved to a new figure
        end

        function build_marks(obj)
            ax = obj.main_axes;
            M = struct('vs', struct('mk',{},'stem',{}), 'lsStem', [], 'lsLabels', [], 'plane', []);
            grey = [0.5 0.5 0.5];
            for k = 1:numel(obj.virtual_source_points)
                h = obj.virtual_source_points{k};
                if isempty(h) || ~isgraphics(h)
                    M.vs(k).mk = gobjects(1); M.vs(k).stem = gobjects(1);
                    continue;
                end
                p = obj.patch_pos(h);
                vis = getappdata(h, 'objVisible');
                if isempty(vis), vis = true; end
                vis = matlab.lang.OnOffSwitchState(vis);
                M.vs(k).stem = line(ax, p([1 1]), p([2 2]), [0 p(3)], 'Color', grey, ...
                    'LineStyle', ':', 'HitTest', 'off', 'PickableParts', 'none', 'Visible', vis);
                M.vs(k).mk = line(ax, p(1), p(2), p(3), 'LineStyle', 'none', 'Marker', 'o', ...
                    'MarkerSize', 10, 'MarkerFaceColor', h.FaceVertexCData(1,:), ...
                    'MarkerEdgeColor', 'k', 'Visible', vis, ...
                    'ButtonDownFcn', @(~,~) obj.start_side_drag('vs', k));
            end
            ls = obj.loudspeaker_points(cellfun(@(h) ~isempty(h) && isgraphics(h) && strcmp(h.Visible, 'on'), obj.loudspeaker_points));
            if ~isempty(ls)
                P = cell2mat(cellfun(@(h) obj.patch_pos(h), ls(:), 'UniformOutput', false));
                M.lsStem = line(ax, obj.stem_data(P(:,1)), obj.stem_data(P(:,2)), ...
                    obj.stem_data(P(:,3), 0), 'Color', grey, 'LineStyle', ':', ...
                    'HitTest', 'off', 'PickableParts', 'none');
                M.lsLabels = cellfun(@(h) h.UserData.Label, ls(:));
            end
            if strcmp(obj.ViewMode, '3d')
                xl = xlim(ax); yl = ylim(ax);
                M.plane = patch(ax, xl([1 2 2 1]), yl([1 1 2 2]), -0.002*[1 1 1 1], [0.8 0.88 1], ...
                    'FaceAlpha', 0.25, 'EdgeColor', 'none', 'HitTest', 'off', ...
                    'PickableParts', 'none', 'Tag', 'EarPlane');
            elseif strcmp(obj.ViewMode, 'listener')
                % first person: a floor for depth (the ear plane is edge-on)
                xl = xlim(ax); yl = ylim(ax);
                M.plane = patch(ax, xl([1 2 2 1]), yl([1 1 2 2]), obj.FloorZ*[1 1 1 1], [0.86 0.86 0.84], ...
                    'EdgeColor', 'none', 'HitTest', 'off', 'PickableParts', 'none', 'Tag', 'EarPlane');
            end
            obj.marks = M;
        end

        function delete_marks(obj)
            if isempty(fieldnames(obj.marks)), return; end
            M = obj.marks;
            h = [M.vs.mk, M.vs.stem];
            for f = {'lsStem','plane'}
                if ~isempty(M.(f{1})), h = [h, M.(f{1})]; end %#ok<AGROW>
            end
            if ~isempty(h), delete(h(isgraphics(h))); end
            obj.marks = struct();
        end

        function update_vs_mark(obj, k, p)
            if ~isfield(obj.marks,'vs') || k > numel(obj.marks.vs) || ~isgraphics(obj.marks.vs(k).mk)
                return;
            end
            set(obj.marks.vs(k).mk, 'XData', p(1), 'YData', p(2), 'ZData', p(3));
            set(obj.marks.vs(k).stem, 'XData', p([1 1]), 'YData', p([2 2]), 'ZData', [0 p(3)]);
        end

        function on_press_3d(obj, h)
            % Loudspeaker model clicked in the 3D view
            fig = ancestor(obj.main_axes, 'figure');
            if strcmp(fig.SelectionType, 'open')
                dbl = getappdata(h, 'user_dblclickfcn');
                if ~isempty(dbl), dbl(h); end
            else
                obj.on_background_press();       % orbit / shift: pan
            end
        end

        function fit_scene(obj)
            % Limits around every scene object (+ margin), at least the
            % listener's +-1.5 m neighbourhood
            ax = obj.main_axes;
            scene = [obj.virtual_source_points(:)', obj.loudspeaker_points(:)', obj.designed_speaker_points(:)', {obj.receiver}];
            P = [0 0 0];
            for k = 1:numel(scene)
                if ~isempty(scene{k}) && isgraphics(scene{k})
                    P(end+1,:) = obj.patch_pos(scene{k}); %#ok<AGROW>
                end
            end
            m = 0.75;
            xl = [min(P(:,1)) - m, max(P(:,1)) + m];
            yl = [min(P(:,2)) - m, max(P(:,2)) + m];
            xl = mean(xl) + [-1 1]*max(diff(xl)/2, 1.5);
            yl = mean(yl) + [-1 1]*max(diff(yl)/2, 1.5);
            xlim(ax, xl); ylim(ax, yl); zlim(ax, obj.z_range());
        end

        function zl = z_range(obj)
            scene = [obj.virtual_source_points(:)', obj.loudspeaker_points(:)', obj.designed_speaker_points(:)', {obj.receiver}];
            z = 0;
            for k = 1:numel(scene)
                if ~isempty(scene{k}) && isgraphics(scene{k})
                    p = obj.patch_pos(scene{k});
                    z(end+1) = p(3); %#ok<AGROW>
                end
            end
            zl = [min(-1, min(z) - 0.3), max(2, max(z) + 0.3)];
        end

        function start_side_drag(obj, kind, idx)
            % Rear view: drag a marker in the screen plane (sideways and
            % height). The patch, marker and model follow; the renderer
            % update goes through Move3DFcn -> sound_scene.schedule_update.
            if ~strcmp(obj.ViewMode, 'rear'), return; end
            ax = obj.main_axes;
            fig = ancestor(ax, 'figure');
            if strcmp(kind, 'ls') && strcmp(fig.SelectionType, 'open')
                % double-click: loudspeaker test signal, as in the top view
                h = obj.loudspeaker_points{idx};
                dbl = getappdata(h, 'user_dblclickfcn');
                if ~isempty(dbl), dbl(h); end
                return;
            end
            cp0 = ax.CurrentPoint(1,:);
            [az, ~] = view(ax);
            v = [-sind(az), cosd(az), 0];          % viewing direction
            if strcmp(kind, 'ls')
                h = obj.loudspeaker_points{idx};
                j = find(obj.marks.lsLabels == idx, 1);
            elseif strcmp(kind, 'ds')
                h = obj.designed_speaker_points{idx};
            else
                h = obj.virtual_source_points{idx};
            end
            p0 = obj.patch_pos(h);
            prevMotion = fig.WindowButtonMotionFcn;
            prevUp = fig.WindowButtonUpFcn;
            fig.WindowButtonMotionFcn = @(~,~) moveTo();
            fig.WindowButtonUpFcn = @(~,~) stop();

            function moveTo()
                d = ax.CurrentPoint(1,:) - cp0;
                d = d - (d*v')*v;                  % stay in the screen plane
                p = p0 + d;
                obj.set_patch_pos(h, p);
                if strcmp(kind, 'ls')
                    if ~isempty(j)
                        st = obj.marks.lsStem;
                        st.XData(3*j-2:3*j-1) = p(1); st.YData(3*j-2:3*j-1) = p(2);
                        st.ZData(3*j-1) = p(3);
                    end
                elseif strcmp(kind, 'ds')
                    % the model itself is the marker
                else
                    obj.update_vs_mark(idx, p);
                end
                if ~isempty(obj.Move3DFcn)
                    obj.Move3DFcn(kind, idx, p);
                end
            end
            function stop()
                fig.WindowButtonMotionFcn = prevMotion;
                fig.WindowButtonUpFcn = prevUp;
            end
        end

        % --- Mouse navigation: wheel = zoom about the cursor (all views);
        %     drag on the background = pan (top/rear) or orbit (3D) ---
        function install_mouse_nav(obj)
            fig = ancestor(obj.main_axes, 'figure');
            if isempty(fig) || isappdata(fig, 'SceneWheelNav'), return; end
            prev = fig.WindowScrollWheelFcn;
            fig.WindowScrollWheelFcn = @(f, e) obj.on_scroll(f, e, prev);
            setappdata(fig, 'SceneWheelNav', true);
            obj.main_axes.ButtonDownFcn = @(~,~) obj.on_background_press();
        end

        function on_scroll(obj, fig, evt, prev)
            ax = obj.main_axes;
            if ~isgraphics(ax) || ~obj.pointer_over_axes(fig)
                if ~isempty(prev), feval(prev, fig, evt); end
                return;
            end
            f = 1.15^evt.VerticalScrollCount;       % wheel down = zoom out
            if strcmp(obj.ViewMode, 'top')
                obj.zoomAxes(f, [ax.CurrentPoint(1,1:2), 0]);
            elseif strcmp(obj.ViewMode, 'rear')
                obj.zoomAxes(f, ax.CurrentPoint(1,:));
            else
                obj.zoomAxes(f);
            end
            obj.positionZoomPanel();
        end

        function tf = pointer_over_axes(obj, fig)
            % only while the axes are visible (not on another top-level tab)
            tab = ancestor(obj.main_axes, 'uitab');
            if ~isempty(tab) && ~isequal(tab.Parent.SelectedTab, tab), tf = false; return; end
            u = fig.Units; fig.Units = 'pixels';
            p = fig.CurrentPoint;
            fig.Units = u;
            b = getpixelposition(obj.main_axes, true);
            tf = p(1) >= b(1) && p(1) <= b(1)+b(3) && p(2) >= b(2) && p(2) <= b(2)+b(4);
        end

        function on_background_press(obj)
            % plain drag: orbit (3D) / pan (top, rear); shift+drag: pan
            if ~isempty(obj.BackgroundPressHook) && obj.BackgroundPressHook()
                return;                          % e.g. a trajectory point was placed
            end
            fig = ancestor(obj.main_axes, 'figure');
            if strcmp(obj.ViewMode, 'listener')
                if any(strcmp(fig.SelectionType, {'normal', 'extend'})), obj.start_look(); end
                return;
            end
            switch fig.SelectionType
                case 'normal'
                    if strcmp(obj.ViewMode, '3d')
                        obj.start_orbit();
                    else
                        obj.start_pan();
                    end
                case 'extend'
                    if strcmp(obj.ViewMode, '3d')
                        obj.start_pan3d();
                    else
                        obj.start_pan();
                    end
            end
        end

        % --- First person ('listener' view) ---
        function apply_fp_camera(obj)
            % perspective camera in the listener's head, looking where the
            % listener faces (+ the camera-only pitch)
            ax = obj.main_axes;
            if isempty(obj.receiver) || ~isgraphics(obj.receiver) || ~isgraphics(ax), return; end
            p = obj.patch_pos(obj.receiver) + [0 0 0.05];          % eyes slightly above the ears
            yaw = obj.receiver.UserData.Orientation;
            d = [cosd(yaw)*cosd(obj.FpPitch), sind(yaw)*cosd(obj.FpPitch), sind(obj.FpPitch)];
            set(ax, 'Projection', 'perspective', 'CameraPosition', p, 'CameraTarget', p + d, ...
                'CameraUpVector', [0 0 1], 'CameraViewAngle', obj.FpFov);
        end

        function fp_move(obj, fwd, side, dyaw)
            % walk fwd / side [m] (in the facing direction) and turn by dyaw
            % [deg]: the listener itself moves (renderers follow)
            if isempty(obj.receiver) || ~isgraphics(obj.receiver), return; end
            yaw = obj.receiver.UserData.Orientation + dyaw;
            xy = obj.receiver.UserData.Origin + fwd*[cosd(yaw), sind(yaw)] + side*[sind(yaw), -cosd(yaw)];
            if ~isempty(obj.ReceiverMoveFcn)
                obj.ReceiverMoveFcn(xy, yaw);      % -> update_receiver -> camera
            else
                obj.update_receiver(xy, [cosd(yaw), sind(yaw)]);
            end
        end

        function start_look(obj)
            % drag: sideways turns the listener, up/down tilts the camera
            ax = obj.main_axes;
            fig = ancestor(ax, 'figure');
            pLast = get(0, 'PointerLocation');
            prevMotion = fig.WindowButtonMotionFcn;
            prevUp = fig.WindowButtonUpFcn;
            fig.WindowButtonMotionFcn = @(~,~) look();
            fig.WindowButtonUpFcn = @(~,~) stop();

            function look()
                p = get(0, 'PointerLocation');
                d = p - pLast; pLast = p;
                s = 0.25*obj.FpFov/75;             % deg per pixel, finer when zoomed in
                obj.FpPitch = min(80, max(-80, obj.FpPitch + s*d(2)));
                if d(1) ~= 0
                    obj.fp_move(0, 0, -s*d(1));    % drag right -> look right
                else
                    obj.apply_fp_camera();
                end
            end
            function stop()
                fig.WindowButtonMotionFcn = prevMotion;
                fig.WindowButtonUpFcn = prevUp;
            end
        end

        function install_key_nav(obj)
            % first person keys (chained in front of the figure's key
            % handler, e.g. the transport's space = play / pause):
            % W/S forward / back, A/D sideways, arrows turn / tilt, R level
            fig = ancestor(obj.main_axes, 'figure');
            if isempty(fig) || isappdata(fig, 'SceneKeyNav'), return; end
            prev = fig.WindowKeyPressFcn;
            fig.WindowKeyPressFcn = @(f, e) obj.on_key(f, e, prev);
            setappdata(fig, 'SceneKeyNav', true);
        end

        function on_key(obj, fig, evt, prev)
            used = false;
            if strcmp(obj.ViewMode, 'listener') && isgraphics(obj.main_axes) && ...
                    ~isa(fig.CurrentObject, 'matlab.ui.control.UIControl')
                tab = ancestor(obj.main_axes, 'uitab');
                if isempty(tab) || isequal(tab.Parent.SelectedTab, tab)
                    st = 0.1; tr = 5;              % step [m], turn [deg]
                    used = true;
                    switch evt.Key
                        case 'w',          obj.fp_move(st, 0, 0);
                        case 's',          obj.fp_move(-st, 0, 0);
                        case 'a',          obj.fp_move(0, -st, 0);
                        case 'd',          obj.fp_move(0, st, 0);
                        case 'leftarrow',  obj.fp_move(0, 0, tr);
                        case 'rightarrow', obj.fp_move(0, 0, -tr);
                        case 'uparrow',    obj.FpPitch = min(80, obj.FpPitch + tr); obj.apply_fp_camera();
                        case 'downarrow',  obj.FpPitch = max(-80, obj.FpPitch - tr); obj.apply_fp_camera();
                        case 'r',          obj.FpPitch = 0; obj.FpFov = 75; obj.apply_fp_camera();
                        otherwise,         used = false;
                    end
                end
            end
            if ~used && ~isempty(prev), feval(prev, fig, evt); end
        end

        function start_pan3d(obj)
            % 3D view: move camera and target together in the screen plane
            ax = obj.main_axes;
            fig = ancestor(ax, 'figure');
            pLast = get(0, 'PointerLocation');
            prevMotion = fig.WindowButtonMotionFcn;
            prevUp = fig.WindowButtonUpFcn;
            fig.WindowButtonMotionFcn = @(~,~) pan();
            fig.WindowButtonUpFcn = @(~,~) stop();

            function pan()
                p = get(0, 'PointerLocation');
                d = p - pLast;
                pLast = p;
                camdolly(ax, -d(1), -d(2), 0, 'movetarget', 'pixels');
            end
            function stop()
                fig.WindowButtonMotionFcn = prevMotion;
                fig.WindowButtonUpFcn = prevUp;
            end
        end

        function start_pan(obj)
            % Drag the view: the data point under the cursor stays under it
            ax = obj.main_axes;
            fig = ancestor(ax, 'figure');
            cp0 = ax.CurrentPoint(1,:);
            [az, ~] = view(ax);
            v = [-sind(az), cosd(az), 0];          % viewing direction (rear)
            top = strcmp(obj.ViewMode, 'top');
            prevMotion = fig.WindowButtonMotionFcn;
            prevUp = fig.WindowButtonUpFcn;
            fig.WindowButtonMotionFcn = @(~,~) pan();
            fig.WindowButtonUpFcn = @(~,~) stop();

            function pan()
                d = ax.CurrentPoint(1,:) - cp0;
                if top
                    d(3) = 0;
                else
                    d = d - (d*v')*v;              % screen plane only
                end
                xlim(ax, xlim(ax) - d(1)); ylim(ax, ylim(ax) - d(2));
                if d(3) ~= 0, zlim(ax, zlim(ax) - d(3)); end
            end
            function stop()
                fig.WindowButtonMotionFcn = prevMotion;
                fig.WindowButtonUpFcn = prevUp;
            end
        end

        function passthrough_overlays(obj)
            % Overlay graphics must not swallow clicks: a press on them
            % goes to the axes (pan/orbit), like a press on the background
            ax = obj.main_axes;
            kids = ax.Children;
            ov = kids(~ismember(kids, obj.scene_handles()));
            if isempty(ov), return; end
            h = findobj(ov, '-property', 'HitTest', 'HitTest', 'on');
            h = h(~arrayfun(@(x) isa(x, 'matlab.graphics.primitive.Text'), h));
            if ~isempty(h), set(h, 'HitTest', 'off'); end
        end

        function start_orbit(obj)
            % 3D view: left-drag on the axes orbits the camera
            ax = obj.main_axes;
            fig = ancestor(ax, 'figure');
            p0 = get(0, 'PointerLocation');
            [az0, el0] = view(ax);
            prevMotion = fig.WindowButtonMotionFcn;
            prevUp = fig.WindowButtonUpFcn;
            fig.WindowButtonMotionFcn = @(~,~) orbit();
            fig.WindowButtonUpFcn = @(~,~) stop();

            function orbit()
                d = get(0, 'PointerLocation') - p0;
                view(ax, az0 - 0.4*d(1), min(89, max(-89, el0 - 0.4*d(2))));
            end
            function stop()
                fig.WindowButtonMotionFcn = prevMotion;
                fig.WindowButtonUpFcn = prevUp;
            end
        end

        function L = ensureLayers(obj)
            % Create (or reuse) 3 groups with stable stacking:
            % bg (zones) < mid (amp/taper) < fg (lines/labels)
            if ~isfield(obj.renderer_illustration,'layers') || ...
               ~all(isfield(obj.renderer_illustration.layers,{'bg','mid','fg'})) || ...
               ~all(isgraphics(struct2array(obj.renderer_illustration.layers)))

                ax = obj.main_axes;
                L.bg  = hggroup('Parent',ax,'Tag','layer-bg','HitTest','off','PickableParts','none');
                L.mid = hggroup('Parent',ax,'Tag','layer-mid','HitTest','off','PickableParts','none');
                L.fg  = hggroup('Parent',ax,'Tag','layer-fg','HitTest','off','PickableParts','none');
                try
                    uistack(L.bg ,'bottom');  % run ONCE
                    uistack(L.mid,'top');
                    uistack(L.fg ,'top');
                catch
                end
                obj.renderer_illustration.layers = L;
                obj.apply_sort_method();
                set(ancestor(ax,'figure'),'Renderer','opengl','GraphicsSmoothing','off');
            else
                L = obj.renderer_illustration.layers;
            end
        end

        function clearRendererIllustration(obj)
            if isempty(obj.renderer_illustration)
                obj.activeRendererKey = '';
                return;
            end

            keepLayers = struct();
            if isfield(obj.renderer_illustration,'layers')
                keepLayers.layers = obj.renderer_illustration.layers;

                % Delete any overlay children parented to the persistent layer
                % groups. Some renderers also track their handles in nested
                % stores; the recursive pass below catches those as well.
                layerNames = fieldnames(keepLayers.layers);
                for k = 1:numel(layerNames)
                    hLayer = keepLayers.layers.(layerNames{k});
                    if isgraphics(hLayer)
                        try
                            delete(allchild(hLayer));
                        catch
                        end
                    end
                end
            end

            fns = fieldnames(obj.renderer_illustration);
            for k = 1:numel(fns)
                if strcmp(fns{k},'layers'), continue; end
                obj.deleteGraphicsRecursive(obj.renderer_illustration.(fns{k}));
            end

            obj.renderer_illustration = keepLayers;
            obj.activeRendererKey = '';
        end

        function deleteGraphicsRecursive(obj, h)
            if isempty(h), return; end

            if iscell(h)
                for k = 1:numel(h)
                    obj.deleteGraphicsRecursive(h{k});
                end
                return;
            end

            if isstruct(h)
                for n = 1:numel(h)
                    fns = fieldnames(h(n));
                    for k = 1:numel(fns)
                        obj.deleteGraphicsRecursive(h(n).(fns{k}));
                    end
                end
                return;
            end

            try
                mask = isgraphics(h);
                if any(mask(:))
                    delete(h(mask));
                end
            catch
            end
        end

    end

    methods (Static)
        function k = view_key(mode)
            k = ['v_' regexprep(mode, '\W', '')];   % '3d' -> valid field name
        end

        function t = label_of(h)
            t = [];
            if ~isempty(h) && isgraphics(h) && isfield(h.UserData, 'text')
                t = h.UserData.text;
            end
        end

        function V = speaker_vertices(V0, pos, yaw)
            % model vertices (local) turned by yaw about z and moved to pos
            Rz = [cosd(yaw) -sind(yaw) 0; sind(yaw) cosd(yaw) 0; 0 0 1];
            V = V0*Rz' + pos;
        end

        function place_z(h, pos)
            % Lift a scene patch to the object's height (UserData.Z)
            if numel(pos) >= 3 && pos(3) ~= 0
                V = h.Vertices;
                if size(V,2) < 3, V(:,3) = 0; end
                V(:,3) = V(:,3) + pos(3);
                h.Vertices = V;
                h.UserData.Z = pos(3);
            end
        end

        function p = patch_pos(h)
            % [x y z] of a scene patch: UserData.Origin and its height
            z = 0;
            if isfield(h.UserData, 'Z'), z = h.UserData.Z; end
            p = [h.UserData.Origin(1:2), z];
        end

        function set_patch_pos(h, p)
            p0 = listener_space_axes.patch_pos(h);
            V = h.Vertices;
            if size(V,2) < 3, V(:,3) = 0; end
            h.Vertices = V + (p - p0);
            h.UserData.Origin = p(1:2);
            h.UserData.Z = p(3);
            if isfield(h.UserData,'text') && isgraphics(h.UserData.text)
                h.UserData.text.Position = p + [0.1 0.1 0];
            end
        end

        function d = stem_data(c, base)
            % NaN-separated vertical stems for line(): [c c NaN] per point;
            % with base, [base c NaN] (the z data)
            c = c(:)';
            if nargin < 2
                d = [c; c; nan(size(c))];
            else
                d = [base*ones(size(c)); c; nan(size(c))];
            end
            d = d(:)';
        end
    end
end
