classdef TrajectoryEditor < handle
    %TRAJECTORYEDITOR  Draws and edits source trajectories in the scene axes.
    %
    %   Drawing (top view): startTool(kind) then click points.
    %     line: 2 clicks   arc: start, through, end   circle: centre, radius
    %     polyline / spline: clicks, finish with double-click, right-click
    %     or Enter; Esc cancels.
    %   Editing: drag the square control points of the selected trajectory.
    %   Display: path line (start dot, end arrowhead) and, for every source
    %   assigned to the path, a dot per second of travel (speed at a glance).

    properties
        C
        T                       % audioapp.model.TrajectorySet
        Tool = ''               % '' | 'line' | 'arc' | 'circle' | 'polyline' | 'spline'
        ToolPts = zeros(0,3)
        AppendTo = []           % trajectory Id the drawn primitive extends
        SelectedId = []
        Height = 0              % z of new points [m]
        OnStatus = []           % @(msg) status text for the tab
        OnSelect = []           % @(id) selection changed
        TickDt = 1              % [s] between time ticks
        TickTmax = 60           % [s] ticks shown up to
    end

    properties (Access = private)
        G = struct('id', {}, 'path', {}, 'ends', {}, 'handles', {}, 'ticks', {})
        Preview = gobjects(0)
        Lsn = {}
        prevMotion = []
        prevKey = []
    end

    properties (Constant)
        Need = struct('line', 2, 'arc', 3, 'circle', 2, 'polyline', inf, 'spline', inf)
    end

    methods
        function obj = TrajectoryEditor(C)
            obj.C = C;
            obj.T = C.Trajectories;
            obj.Lsn{end+1} = addlistener(obj.T, 'Changed', @(~,~) obj.redraw());
            % a scene rebuild deletes every graphics object of the axes
            obj.Lsn{end+1} = addlistener(C, 'SetupChanged', @(~,~) obj.reattach());
            obj.reattach();
        end

        function delete(obj)
            cellfun(@delete, obj.Lsn);
            obj.clear_graphics();
            try
                gui = obj.C.SceneGUI;
                gui.BackgroundPressHook = [];
                gui.ExtraSceneHandlesFcn = [];
            catch
            end
        end

        function reattach(obj)
            gui = obj.C.SceneGUI;
            if isempty(gui) || ~isvalid(gui), return; end
            gui.BackgroundPressHook = @() obj.onPress();
            gui.ExtraSceneHandlesFcn = @() obj.allHandles();
            obj.redraw();                        % deletes our old graphics first
        end

        % ================= drawing tools =================
        function startTool(obj, kind, appendTo)
            obj.cancel();
            if nargin < 3, appendTo = []; end
            gui = obj.C.SceneGUI;
            if ~strcmp(gui.ViewMode, 'top'), gui.set_view_mode('top'); end
            obj.Tool = lower(kind);
            obj.AppendTo = appendTo;
            obj.ToolPts = zeros(0,3);
            if ~isempty(appendTo) && ~strcmp(obj.Tool, 'circle')
                tr = obj.T.get(appendTo);
                if ~isempty(tr) && ~isempty(tr.pts)
                    obj.ToolPts = tr.pts(end,:);    % continues from the end
                end
            end
            fig = ancestor(gui.main_axes, 'figure');
            obj.prevMotion = fig.WindowButtonMotionFcn;
            obj.prevKey = fig.WindowKeyPressFcn;
            fig.WindowButtonMotionFcn = @(~,~) obj.onMotion();
            fig.WindowKeyPressFcn = @(~,e) obj.onKey(e);
            fig.Pointer = 'crosshair';
            % scene objects must not catch the clicks while drawing
            scene = [gui.virtual_source_points(:)', gui.loudspeaker_points(:)', {gui.receiver}];
            for k = 1:numel(scene)
                if ~isempty(scene{k}) && isgraphics(scene{k}), scene{k}.PickableParts = 'none'; end
            end
            obj.status(obj.instructions());
        end

        function cancel(obj)
            if isempty(obj.Tool), return; end
            obj.Tool = '';
            obj.ToolPts = zeros(0,3);
            delete(obj.Preview(isgraphics(obj.Preview)));
            obj.Preview = gobjects(0);
            gui = obj.C.SceneGUI;
            if ~isempty(gui) && isvalid(gui)
                fig = ancestor(gui.main_axes, 'figure');
                fig.WindowButtonMotionFcn = obj.prevMotion;
                fig.WindowKeyPressFcn = obj.prevKey;
                fig.Pointer = 'arrow';
                gui.refresh_view();              % scene objects pickable again
            end
            obj.status('');
        end

        function consumed = onPress(obj)
            % BackgroundPressHook of the scene axes
            consumed = false;
            if isempty(obj.Tool), return; end
            consumed = true;
            ax = obj.C.SceneGUI.main_axes;
            fig = ancestor(ax, 'figure');
            switch fig.SelectionType
                case {'open', 'alt'}            % double / right click: finish
                    obj.finish();
                    return;
            end
            cp = ax.CurrentPoint(1,1:2);
            obj.ToolPts(end+1,:) = [cp, obj.Height];
            if size(obj.ToolPts,1) >= obj.Need.(obj.Tool)
                obj.finish();
            else
                obj.drawPreview(cp);
                obj.status(obj.instructions());
            end
        end

        function finish(obj)
            P = obj.ToolPts;
            kind = obj.Tool;
            need = obj.Need.(kind);
            % drop a duplicate final point (the double-click's first click)
            if size(P,1) > 1 && norm(P(end,:) - P(end-1,:)) < 1e-3, P(end,:) = []; end
            ok = size(P,1) >= min(need, 2) && (isinf(need) || size(P,1) >= need);
            appendTo = obj.AppendTo;
            obj.cancel();
            if ~ok, obj.status('Too few points - nothing created.'); return; end
            if ~isempty(appendTo) && ~isempty(obj.T.get(appendTo))
                obj.T.append(appendTo, kind, P);
                obj.select(appendTo);
            else
                tr = obj.T.add(kind, P);
                obj.select(tr.Id);
            end
        end

        function select(obj, id)
            obj.SelectedId = id;
            obj.redraw();
            if ~isempty(obj.OnSelect), obj.OnSelect(id); end
        end

        % ================= graphics =================
        function redraw(obj)
            gui = obj.C.SceneGUI;
            if isempty(gui) || ~isvalid(gui) || ~isgraphics(gui.main_axes), return; end
            ax = gui.main_axes;
            obj.clear_graphics();
            for k = 1:numel(obj.T.Trajectories)
                tr = obj.T.Trajectories{k};
                if isempty(tr.pts), continue; end
                sel = isequal(tr.Id, obj.SelectedId);
                g.id = tr.Id;
                g.path = line(ax, tr.pts(:,1), tr.pts(:,2), tr.pts(:,3), 'Color', tr.Color, ...
                    'LineWidth', 1.5 + sel, 'Tag', 'Trajectory', ...
                    'ButtonDownFcn', @(~,~) obj.onPathPress(tr.Id));
                g.ends = [ ...
                    line(ax, tr.pts(1,1), tr.pts(1,2), tr.pts(1,3), 'Marker', 'o', 'MarkerSize', 6, ...
                        'MarkerFaceColor', tr.Color, 'MarkerEdgeColor', 'none', 'HitTest', 'off'), ...
                    obj.arrowhead(ax, tr)];
                g.handles = gobjects(0);
                if sel
                    for s = 1:numel(tr.Segments)
                        P = tr.Segments(s).P;
                        for j = 1:size(P,1)
                            g.handles(end+1) = line(ax, P(j,1), P(j,2), P(j,3), 'Marker', 's', ...
                                'MarkerSize', 8, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', tr.Color, ...
                                'LineWidth', 1.5, 'ButtonDownFcn', @(~,~) obj.onHandlePress(tr.Id, s, j));
                        end
                    end
                end
                g.ticks = gobjects(0);
                for a = obj.T.Assign([obj.T.Assign.TrajId] == tr.Id)
                    Q = obj.T.ticks(a, obj.TickDt, obj.TickTmax);
                    g.ticks(end+1) = line(ax, Q(:,1), Q(:,2), Q(:,3), 'LineStyle', 'none', ...
                        'Marker', '.', 'MarkerSize', 12, 'Color', tr.Color*0.6, 'HitTest', 'off');
                end
                obj.G(end+1) = g;
            end
            gui.raiseSceneObjects();
        end

        function h = allHandles(obj)
            h = gobjects(0);
            for g = obj.G
                h = [h, g.path, g.ends, g.handles, g.ticks]; %#ok<AGROW>
            end
            h = [h, obj.Preview];
        end
    end

    methods (Access = private)
        function clear_graphics(obj)
            h = obj.allHandles();
            delete(h(isgraphics(h) & ~ismember(h, obj.Preview)));
            obj.G = obj.G([]);
        end

        function h = arrowhead(~, ax, tr)
            % small triangle at the end, pointing along the path
            n = size(tr.pts,1);
            p = tr.pts(end,:);
            d = tr.pts(n,1:2) - tr.pts(max(1, n-5),1:2);
            if norm(d) < eps, d = [1 0]; end
            d = d/norm(d); q = [-d(2) d(1)];
            sz = 0.08;
            V = [p(1:2); p(1:2) - sz*d + 0.5*sz*q; p(1:2) - sz*d - 0.5*sz*q];
            h = patch(ax, V(:,1), V(:,2), p(3)*ones(3,1), tr.Color, 'EdgeColor', 'none', 'HitTest', 'off');
        end

        function onPathPress(obj, id)
            if ~isempty(obj.Tool), obj.onPress(); return; end
            obj.select(id);
        end

        function onHandlePress(obj, id, seg, k)
            if ~isempty(obj.Tool), obj.onPress(); return; end
            tr = obj.T.get(id);
            if isempty(tr), return; end
            ax = obj.C.SceneGUI.main_axes;
            fig = ancestor(ax, 'figure');
            gi = find([obj.G.id] == id, 1);
            prevM = fig.WindowButtonMotionFcn; prevU = fig.WindowButtonUpFcn;
            fig.WindowButtonMotionFcn = @(~,~) moveTo();
            fig.WindowButtonUpFcn = @(~,~) stop();

            function moveTo()
                cp = ax.CurrentPoint(1,1:2);
                tr.set_point(seg, k, cp);
                if ~isempty(gi) && isgraphics(obj.G(gi).path)
                    set(obj.G(gi).path, 'XData', tr.pts(:,1), 'YData', tr.pts(:,2), 'ZData', tr.pts(:,3));
                    hh = obj.G(gi).handles;
                    i0 = sum(arrayfun(@(g) size(g.P,1), tr.Segments(1:seg-1)));
                    if i0 + k <= numel(hh), set(hh(i0+k), 'XData', cp(1), 'YData', cp(2)); end
                end
            end
            function stop()
                fig.WindowButtonMotionFcn = prevM;
                fig.WindowButtonUpFcn = prevU;
                notify(obj.T, 'Changed');            % full redraw (ticks, joints), tab refresh
            end
        end

        function onMotion(obj)
            if isempty(obj.Tool), return; end
            cp = obj.C.SceneGUI.main_axes.CurrentPoint(1,1:2);
            obj.drawPreview(cp);
        end

        function onKey(obj, e)
            switch e.Key
                case {'return', 'enter'}, obj.finish();
                case 'escape', obj.cancel();
            end
        end

        function drawPreview(obj, cp)
            ax = obj.C.SceneGUI.main_axes;
            P = [obj.ToolPts; cp, obj.Height];
            q = P;
            try
                if size(P,1) >= 2
                    kind = obj.Tool;
                    if strcmp(kind, 'arc') && size(P,1) < 3, kind = 'line'; end
                    tmp = trajectory(0, '', kind, P);
                    q = tmp.pts;
                end
            catch
            end
            V = P(1:end-1,:);                    % points placed so far
            if isempty(V), V = nan(1,3); end     % (an empty line() creates nothing)
            if numel(obj.Preview) ~= 2 || ~all(isgraphics(obj.Preview))
                delete(obj.Preview(isgraphics(obj.Preview)));
                obj.Preview = [ ...
                    line(ax, q(:,1), q(:,2), q(:,3), 'LineStyle', '--', 'Color', [0.3 0.3 0.3], 'HitTest', 'off'), ...
                    line(ax, V(:,1), V(:,2), V(:,3), 'LineStyle', 'none', 'Marker', 's', ...
                        'MarkerSize', 7, 'MarkerEdgeColor', [0.3 0.3 0.3], 'HitTest', 'off')];
            else
                set(obj.Preview(1), 'XData', q(:,1), 'YData', q(:,2), 'ZData', q(:,3));
                set(obj.Preview(2), 'XData', V(:,1), 'YData', V(:,2), 'ZData', V(:,3));
            end
        end

        function msg = instructions(obj)
            n = size(obj.ToolPts,1);
            switch obj.Tool
                case 'line',   steps = {'click the start point', 'click the end point'};
                case 'arc',    steps = {'click the start point', 'click a point the arc passes through', 'click the end point'};
                case 'circle', steps = {'click the centre', 'click a point on the circle (start)'};
                otherwise
                    msg = sprintf('%s: click points (%d so far); double-click / right-click / Enter finishes, Esc cancels.', ...
                        obj.Tool, n);
                    return;
            end
            msg = sprintf('%s: %s  (Esc cancels)', obj.Tool, steps{min(n+1, numel(steps))});
        end

        function status(obj, msg)
            if ~isempty(obj.OnStatus), obj.OnStatus(msg); end
        end
    end
end
