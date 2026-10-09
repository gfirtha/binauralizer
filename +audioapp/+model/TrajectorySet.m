classdef TrajectorySet < handle
    %TRAJECTORYSET  Designed source trajectories and their assignments.
    %
    %   Lives in AppContext (survives scene rebuilds). Each assignment ties
    %   a virtual source (index) to a trajectory (Id) with a motion_profile;
    %   during playback AudioEngine asks positions(t) with the media time.

    events
        Changed            % trajectories or assignments edited
    end

    properties
        Trajectories = {}  % cell of trajectory
        Assign = struct('Vs', {}, 'TrajId', {}, 'Motion', {})
        NextId = 1
    end

    methods
        % --- trajectories ---
        function tr = add(obj, type, P, name)
            if nargin < 4 || isempty(name)
                name = sprintf('%s %d', type, obj.NextId);
            end
            colors = lines(7);
            tr = trajectory(obj.NextId, name, type, P);
            tr.Color = colors(mod(obj.NextId-1, 7)+1, :);
            obj.NextId = obj.NextId + 1;
            obj.Trajectories{end+1} = tr;
            notify(obj, 'Changed');
        end

        function append(obj, id, type, P)
            tr = obj.get(id);
            if isempty(tr), return; end
            tr.add_segment(type, P);
            notify(obj, 'Changed');
        end

        function remove(obj, id)
            k = obj.index_of(id);
            if isempty(k), return; end
            obj.Trajectories(k) = [];
            obj.Assign([obj.Assign.TrajId] == id) = [];
            notify(obj, 'Changed');
        end

        function tr = get(obj, id)
            k = obj.index_of(id);
            tr = [];
            if ~isempty(k), tr = obj.Trajectories{k}; end
        end

        function k = index_of(obj, id)
            k = find(cellfun(@(t) t.Id == id, obj.Trajectories), 1);
        end

        % --- assignments ---
        function assign(obj, vs, id, motion)
            if nargin < 4, motion = motion_profile(); end
            obj.Assign([obj.Assign.Vs] == vs) = [];
            obj.Assign(end+1) = struct('Vs', vs, 'TrajId', id, 'Motion', motion);
            obj.sort_assign();
            notify(obj, 'Changed');
        end

        function unassign(obj, vs)
            obj.Assign([obj.Assign.Vs] == vs) = [];
            notify(obj, 'Changed');
        end

        function tf = hasAssignments(obj)
            tf = ~isempty(obj.Assign);
        end

        function [idx, P] = positions(obj, t, nVs)
            % Positions of the assigned sources at media time t
            idx = zeros(0,1); P = zeros(0,3);
            for a = obj.Assign
                if a.Vs > nVs, continue; end
                tr = obj.get(a.TrajId);
                if isempty(tr) || isempty(tr.pts), continue; end
                idx(end+1,1) = a.Vs; %#ok<AGROW>
                P(end+1,:) = tr.point_at(a.Motion.s_at(t, tr.L)); %#ok<AGROW>
            end
        end

        function P = ticks(obj, a, dt, tmax)
            % Positions every dt seconds from T0 (display of the speed)
            tr = obj.get(a.TrajId);
            if isempty(tr), P = zeros(0,3); return; end
            t = a.Motion.T0 + (0:dt:tmax)';
            P = tr.point_at(a.Motion.s_at(t, tr.L));
        end

        % --- persistence ---
        function save(obj, file)
            S.Trajectories = cellfun(@(t) t.to_struct(), obj.Trajectories, 'UniformOutput', false);
            S.Assign = arrayfun(@(a) struct('Vs', a.Vs, 'TrajId', a.TrajId, ...
                'Motion', a.Motion.to_struct()), obj.Assign, 'UniformOutput', false);
            fid = fopen(file, 'w');
            fwrite(fid, jsonencode(S, 'PrettyPrint', true), 'char');
            fclose(fid);
        end

        function load(obj, file)
            S = jsondecode(fileread(file));
            T = S.Trajectories; if isstruct(T), T = num2cell(T); end
            obj.Trajectories = cellfun(@(t) trajectory.from_struct(t), T(:)', 'UniformOutput', false);
            obj.NextId = 1 + max([0, cellfun(@(t) t.Id, obj.Trajectories)]);
            obj.Assign = struct('Vs', {}, 'TrajId', {}, 'Motion', {});
            A = S.Assign; if isstruct(A), A = num2cell(A); end
            for k = 1:numel(A)
                obj.Assign(end+1) = struct('Vs', A{k}.Vs, 'TrajId', A{k}.TrajId, ...
                    'Motion', motion_profile.from_struct(A{k}.Motion));
            end
            obj.sort_assign();
            notify(obj, 'Changed');
        end
    end

    methods (Access = private)
        function sort_assign(obj)
            if isempty(obj.Assign), return; end
            [~, o] = sort([obj.Assign.Vs]);
            obj.Assign = obj.Assign(o);
        end
    end
end
