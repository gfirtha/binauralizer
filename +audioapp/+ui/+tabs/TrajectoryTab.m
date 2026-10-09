classdef TrajectoryTab < handle
    % TrajectoryTab  Designs source trajectories and assigns sources to them.
    %
    %   Top: trajectory list, primitive + Draw new / Append / Delete, height.
    %   Middle: assignment table (editable t0, v0, a, speed profile, end mode).
    %   Bottom: preview time slider (moves the sources without audio),
    %   Save / Load. During playback the sources follow the media time.

    properties
        C
        Editor  audioapp.ui.TrajectoryEditor
        LstTraj
        PopPrim
        EdtHeight
        TxtStatus
        PopSource
        Tbl
        SldTime
        LblTime
        Lsn = {}
    end

    properties (Constant)
        Prims = {'Line', 'Arc (3 points)', 'Circle', 'Polyline', 'Spline'}
        PrimKeys = {'line', 'arc', 'circle', 'polyline', 'spline'}
        EndModes = {'loop', 'pingpong', 'stop'}
    end

    methods
        function obj = TrajectoryTab(parent, C)
            obj.C = C;
            obj.Editor = audioapp.ui.TrajectoryEditor(C);
            obj.Editor.OnStatus = @(m) obj.setStatus(m);
            obj.Editor.OnSelect = @(id) obj.refresh();
            obj.build(parent);
            obj.Lsn{end+1} = addlistener(C.Trajectories, 'Changed', @(~,~) obj.refresh());
            obj.Lsn{end+1} = addlistener(C, 'SetupChanged', @(~,~) obj.refresh());
            obj.refresh();
        end

        function delete(obj)
            cellfun(@delete, obj.Lsn);
            if ~isempty(obj.Editor) && isvalid(obj.Editor), delete(obj.Editor); end
        end

        function build(obj, p)
            u = {'Units', 'normalized'};
            % ---- trajectories ----
            uicontrol(p, 'Style', 'text', 'String', 'Trajectories', u{:}, 'Position', [0.04 0.94 0.5 0.04], ...
                'HorizontalAlignment', 'left', 'FontWeight', 'bold');
            obj.LstTraj = uicontrol(p, 'Style', 'listbox', u{:}, 'Position', [0.04 0.72 0.56 0.22], ...
                'Callback', @(s,~) obj.onListSelect(s));
            obj.PopPrim = uicontrol(p, 'Style', 'popupmenu', 'String', obj.Prims, u{:}, ...
                'Position', [0.63 0.895 0.33 0.04]);
            uicontrol(p, 'Style', 'pushbutton', 'String', 'Draw new', u{:}, 'Position', [0.63 0.845 0.33 0.045], ...
                'TooltipString', 'Draw the primitive in the top view', 'Callback', @(~,~) obj.onDraw(false));
            uicontrol(p, 'Style', 'pushbutton', 'String', 'Append to selected', u{:}, 'Position', [0.63 0.795 0.33 0.045], ...
                'TooltipString', 'The primitive continues the selected trajectory', 'Callback', @(~,~) obj.onDraw(true));
            uicontrol(p, 'Style', 'pushbutton', 'String', 'Delete', u{:}, 'Position', [0.63 0.745 0.16 0.045], ...
                'Callback', @(~,~) obj.onDelete());
            uicontrol(p, 'Style', 'text', 'String', 'Height (m):', u{:}, 'Position', [0.04 0.675 0.22 0.035], ...
                'HorizontalAlignment', 'left');
            obj.EdtHeight = uicontrol(p, 'Style', 'edit', 'String', '0', u{:}, 'Position', [0.26 0.675 0.12 0.04], ...
                'TooltipString', 'z of new points; applied to the selected trajectory', 'Callback', @(s,~) obj.onHeight(s));
            obj.TxtStatus = uicontrol(p, 'Style', 'text', 'String', '', u{:}, 'Position', [0.04 0.60 0.92 0.065], ...
                'HorizontalAlignment', 'left', 'FontAngle', 'italic', 'ForegroundColor', [0.1 0.3 0.6]);

            % ---- assignments ----
            uicontrol(p, 'Style', 'text', 'String', 'Sources on trajectories', u{:}, 'Position', [0.04 0.555 0.6 0.04], ...
                'HorizontalAlignment', 'left', 'FontWeight', 'bold');
            obj.PopSource = uicontrol(p, 'Style', 'popupmenu', 'String', {'s_1'}, u{:}, 'Position', [0.04 0.51 0.18 0.04]);
            uicontrol(p, 'Style', 'pushbutton', 'String', 'Assign to selected', u{:}, 'Position', [0.24 0.51 0.34 0.045], ...
                'Callback', @(~,~) obj.onAssign());
            uicontrol(p, 'Style', 'pushbutton', 'String', 'Remove', u{:}, 'Position', [0.60 0.51 0.18 0.045], ...
                'Callback', @(~,~) obj.onUnassign());
            obj.Tbl = uitable(p, u{:}, 'Position', [0.04 0.20 0.92 0.30], ...
                'ColumnName', {'Source', 'Trajectory', 't0 [s]', 'v0 [m/s]', 'a [m/s²]', 'Speed profile t:v', 'End'}, ...
                'ColumnFormat', {'char', 'char', 'numeric', 'numeric', 'numeric', 'char', obj.EndModes}, ...
                'ColumnEditable', [false false true true true true true], ...
                'ColumnWidth', {45, 70, 45, 55, 55, 110, 65}, 'RowName', [], ...
                'TooltipString', 'Speed profile overrides v0/a, e.g. 0:1, 5:3, 10:0 (time from t0 : speed)', ...
                'CellEditCallback', @(~,e) obj.onCellEdit(e));

            % ---- preview / files ----
            uicontrol(p, 'Style', 'text', 'String', 'Preview time:', u{:}, 'Position', [0.04 0.13 0.22 0.035], ...
                'HorizontalAlignment', 'left');
            obj.SldTime = uicontrol(p, 'Style', 'slider', 'Min', 0, 'Max', 60, 'Value', 0, u{:}, ...
                'Position', [0.26 0.135 0.54 0.035], 'Callback', @(s,~) obj.onPreview(s));
            obj.LblTime = uicontrol(p, 'Style', 'text', 'String', '0.0 s', u{:}, 'Position', [0.81 0.13 0.15 0.035], ...
                'HorizontalAlignment', 'left');
            uicontrol(p, 'Style', 'pushbutton', 'String', 'Save...', u{:}, 'Position', [0.04 0.05 0.2 0.05], ...
                'Callback', @(~,~) obj.onSave());
            uicontrol(p, 'Style', 'pushbutton', 'String', 'Load...', u{:}, 'Position', [0.26 0.05 0.2 0.05], ...
                'Callback', @(~,~) obj.onLoad());
            uicontrol(p, 'Style', 'text', 'String', 'During playback the sources follow the media time.', u{:}, ...
                'Position', [0.48 0.04 0.5 0.06], 'HorizontalAlignment', 'left', 'FontSize', 8, 'FontAngle', 'italic');
        end

        function refresh(obj)
            if isempty(obj.LstTraj) || ~isgraphics(obj.LstTraj), return; end
            T = obj.C.Trajectories;
            names = cellfun(@(t) sprintf('%s   (%.2f m)', t.Name, t.L), T.Trajectories, 'UniformOutput', false);
            ids = cellfun(@(t) t.Id, T.Trajectories);
            set(obj.LstTraj, 'String', names, 'UserData', ids);
            k = find(ids == obj.Editor.SelectedId, 1);
            if isempty(k), k = min(1, numel(ids)); end
            set(obj.LstTraj, 'Value', max(k, 1));
            if isempty(ids), set(obj.LstTraj, 'Value', 1); end
            nVs = numel(obj.C.Scene.virtual_sources);
            set(obj.PopSource, 'String', arrayfun(@(i) sprintf('s_%d', i), 1:nVs, 'UniformOutput', false), ...
                'Value', min(get(obj.PopSource, 'Value'), max(nVs, 1)));
            D = cell(numel(T.Assign), 7);
            for i = 1:numel(T.Assign)
                a = T.Assign(i); m = a.Motion; tr = T.get(a.TrajId);
                D(i,:) = {sprintf('s_%d', a.Vs), tr.Name, m.T0, m.V0, m.A, ...
                    motion_profile.format_profile(m.Profile), m.EndMode};
            end
            set(obj.Tbl, 'Data', D);
            dur = 0;
            try, dur = obj.C.Engine.getDuration(); catch, end
            if dur > 0, set(obj.SldTime, 'Max', dur); end
        end

        function setStatus(obj, msg)
            if ~isempty(obj.TxtStatus) && isgraphics(obj.TxtStatus), set(obj.TxtStatus, 'String', msg); end
        end
    end

    methods (Access = private)
        function id = selectedId(obj)
            ids = get(obj.LstTraj, 'UserData');
            id = [];
            if ~isempty(ids), id = ids(get(obj.LstTraj, 'Value')); end
        end

        function onListSelect(obj, ~)
            id = obj.selectedId();
            if ~isempty(id), obj.Editor.select(id); end
        end

        function onDraw(obj, append)
            kind = obj.PrimKeys{get(obj.PopPrim, 'Value')};
            target = [];
            if append
                target = obj.selectedId();
                if isempty(target), obj.setStatus('Select a trajectory to append to.'); return; end
            end
            obj.Editor.Height = str2double(get(obj.EdtHeight, 'String'));
            obj.Editor.startTool(kind, target);
        end

        function onDelete(obj)
            id = obj.selectedId();
            if ~isempty(id), obj.C.Trajectories.remove(id); end
        end

        function onHeight(obj, src)
            z = str2double(get(src, 'String'));
            if isnan(z), set(src, 'String', '0'); return; end
            obj.Editor.Height = z;
            tr = obj.C.Trajectories.get(obj.selectedId());
            if ~isempty(tr)
                tr.set_height(z);
                notify(obj.C.Trajectories, 'Changed');
            end
        end

        function onAssign(obj)
            id = obj.selectedId();
            if isempty(id), obj.setStatus('Draw or select a trajectory first.'); return; end
            vs = get(obj.PopSource, 'Value');
            obj.C.Trajectories.assign(vs, id, motion_profile('V0', 1, 'EndMode', 'loop'));
        end

        function onUnassign(obj)
            vs = get(obj.PopSource, 'Value');
            obj.C.Trajectories.unassign(vs);
        end

        function onCellEdit(obj, e)
            T = obj.C.Trajectories;
            i = e.Indices(1);
            if i > numel(T.Assign), return; end
            m = T.Assign(i).Motion;
            v = e.NewData;
            switch e.Indices(2)
                case 3, m.T0 = obj.num(v, m.T0);
                case 4, m.V0 = max(0, obj.num(v, m.V0));
                case 5, m.A  = obj.num(v, m.A);
                case 6, m.Profile = motion_profile.parse_profile(v);
                case 7, m.EndMode = char(v);
            end
            T.Assign(i).Motion = m;
            notify(T, 'Changed');
        end

        function onPreview(obj, src)
            t = get(src, 'Value');
            set(obj.LblTime, 'String', sprintf('%.1f s', t));
            obj.C.ensureEngine();
            obj.C.Engine.applyTrajectories(t);
        end

        function onSave(obj)
            [f, p] = uiputfile('*.json', 'Save trajectories', 'trajectories.json');
            if isequal(f, 0), return; end
            obj.C.Trajectories.save(fullfile(p, f));
        end

        function onLoad(obj)
            [f, p] = uigetfile('*.json', 'Load trajectories');
            if isequal(f, 0), return; end
            try
                obj.C.Trajectories.load(fullfile(p, f));
            catch ME
                errordlg(sprintf('Could not load trajectories:\n%s', ME.message), 'Load failed', 'modal');
            end
        end
    end

    methods (Static, Access = private)
        function x = num(v, def)
            if ischar(v) || isstring(v), v = str2double(v); end
            if isempty(v) || ~isfinite(v), x = def; else, x = v; end
        end
    end
end
