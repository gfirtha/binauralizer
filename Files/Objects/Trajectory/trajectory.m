classdef trajectory < handle
    %TRAJECTORY  Source path built from geometric primitives.
    %
    %   A trajectory is a chain of segments; each segment is a primitive
    %   defined by its control points (rows [x y z]):
    %     'line'     [p1; p2]
    %     'arc'      [p1; pm; p2]   circle arc from p1 through pm to p2
    %     'circle'   [c; p]         full circle around c, starting at p (CCW)
    %     'polyline' [p1; ...; pn]
    %     'spline'   [p1; ...; pn]  Catmull-Rom spline through the points
    %   Heights: arcs/circles interpolate z linearly between their ends.
    %
    %   The path is sampled densely and parametrised by arc length s in
    %   [0, L], so speeds are real m/s whatever the primitive.

    properties
        Id      = 0
        Name    = ''
        Segments = struct('type', {}, 'P', {})
        Color   = [0.85 0.33 0.10]
    end

    properties (SetAccess = private)
        pts = zeros(0,3)        % samples along the path
        s   = zeros(0,1)        % arc length at the samples
        L   = 0                 % total length [m]
        segStart = []           % sample index where each segment starts
    end

    properties (Constant)
        Spacing = 0.01          % sampling step [m]
    end

    methods
        function obj = trajectory(id, name, type, P)
            obj.Id = id;
            obj.Name = name;
            if nargin > 2
                obj.add_segment(type, P);
            end
        end

        function add_segment(obj, type, P)
            % Append a primitive; it starts where the path currently ends
            P = obj.pad3(P);
            if ~isempty(obj.pts) && ~strcmp(type, 'circle')
                P(1,:) = obj.pts(end,:);
            end
            obj.Segments(end+1) = struct('type', type, 'P', P);
            obj.build();
        end

        function set_point(obj, seg, k, p)
            % Move control point k of segment seg (keeps segments joined)
            P = obj.Segments(seg).P;
            P(k, 1:numel(p)) = p;
            obj.Segments(seg).P = P;
            nP = size(P,1);
            if k == nP && seg < numel(obj.Segments) && ~strcmp(obj.Segments(seg+1).type, 'circle')
                obj.Segments(seg+1).P(1,:) = P(k,:);
            end
            if k == 1 && seg > 1 && ~strcmp(obj.Segments(seg).type, 'circle')
                obj.Segments(seg-1).P(end,:) = P(1,:);
            end
            obj.build();
        end

        function set_height(obj, z)
            for k = 1:numel(obj.Segments)
                obj.Segments(k).P(:,3) = z;
            end
            obj.build();
        end

        function build(obj)
            pts = zeros(0,3); obj.segStart = zeros(1, numel(obj.Segments));
            for k = 1:numel(obj.Segments)
                q = obj.sample_segment(obj.Segments(k));
                if ~isempty(pts) && norm(q(1,:) - pts(end,:)) < 1e-9
                    q = q(2:end,:);
                end
                obj.segStart(k) = size(pts,1) + 1;
                pts = [pts; q]; %#ok<AGROW>
            end
            obj.pts = pts;
            if size(pts,1) < 2
                obj.s = zeros(size(pts,1),1); obj.L = 0; return;
            end
            ds = sqrt(sum(diff(pts).^2, 2));
            obj.s = [0; cumsum(ds)];
            obj.L = obj.s(end);
        end

        function p = point_at(obj, s)
            % Points at arc lengths s (clamped to [0, L]); [numel(s) x 3]
            if isempty(obj.pts), p = zeros(numel(s), 3); return; end
            if obj.L == 0, p = repmat(obj.pts(1,:), numel(s), 1); return; end
            s = min(max(s(:), 0), obj.L);
            [su, iu] = unique(obj.s);
            p = interp1(su, obj.pts(iu,:), s, 'linear');
        end

        function S = to_struct(obj)
            S = struct('Id', obj.Id, 'Name', obj.Name, 'Color', obj.Color, ...
                'Segments', {arrayfun(@(g) struct('type', g.type, 'P', g.P), obj.Segments, 'UniformOutput', false)});
        end
    end

    methods (Static)
        function obj = from_struct(S)
            obj = trajectory(S.Id, S.Name);
            if isfield(S, 'Color'), obj.Color = S.Color(:)'; end
            seg = S.Segments;
            if isstruct(seg), seg = num2cell(seg); end
            for k = 1:numel(seg)
                obj.Segments(end+1) = struct('type', seg{k}.type, 'P', trajectory.pad3(seg{k}.P));
            end
            obj.build();
        end

        function P = pad3(P)
            if size(P,2) < 3, P(:,3) = 0; end
            P = P(:,1:3);
        end
    end

    methods (Static, Access = private)
        function q = sample_segment(g)
            P = g.P; h = trajectory.Spacing;
            switch g.type
                case 'line'
                    n = max(2, ceil(norm(P(2,:) - P(1,:))/h) + 1);
                    t = linspace(0, 1, n)';
                    q = P(1,:) + t.*(P(2,:) - P(1,:));
                case 'arc'
                    [c, r, a1, da] = trajectory.arc_params(P(1,1:2), P(2,1:2), P(3,1:2));
                    if isinf(r)                         % collinear: straight line
                        q = trajectory.sample_segment(struct('type','line','P',P([1 3],:)));
                        return;
                    end
                    n = max(8, ceil(abs(da)*r/h) + 1);
                    a = a1 + linspace(0, da, n)';
                    q = [c + r*[cos(a), sin(a)], linspace(P(1,3), P(3,3), n)'];
                case 'circle'
                    c = P(1,1:2); v = P(2,1:2) - c; r = norm(v);
                    a0 = atan2(v(2), v(1));
                    n = max(16, ceil(2*pi*r/h) + 1);
                    a = a0 + linspace(0, 2*pi, n)';
                    q = [c + r*[cos(a), sin(a)], P(2,3)*ones(n,1)];
                case 'polyline'
                    q = P(1,:);
                    for k = 2:size(P,1)
                        seg = trajectory.sample_segment(struct('type','line','P',P(k-1:k,:)));
                        q = [q; seg(2:end,:)]; %#ok<AGROW>
                    end
                case 'spline'
                    if size(P,1) < 3
                        q = trajectory.sample_segment(struct('type','polyline','P',P));
                        return;
                    end
                    Q = [2*P(1,:) - P(2,:); P; 2*P(end,:) - P(end-1,:)];   % end tangents
                    q = P(1,:);
                    for k = 2:size(Q,1)-2
                        p0 = Q(k-1,:); p1 = Q(k,:); p2 = Q(k+1,:); p3 = Q(k+2,:);
                        n = max(4, ceil(norm(p2 - p1)/h*1.2));
                        t = linspace(0, 1, n+1)'; t = t(2:end);
                        % uniform Catmull-Rom
                        q = [q; 0.5*((2*p1) + (-p0 + p2).*t + (2*p0 - 5*p1 + 4*p2 - p3).*t.^2 ...
                            + (-p0 + 3*p1 - 3*p2 + p3).*t.^3)]; %#ok<AGROW>
                    end
                otherwise
                    error('trajectory:type', 'Unknown primitive "%s".', g.type);
            end
        end

        function [c, r, a1, da] = arc_params(p1, pm, p2)
            % Circle through three points; signed sweep from p1 via pm to p2
            A = 2*[pm - p1; p2 - p1];
            b = [sum(pm.^2) - sum(p1.^2); sum(p2.^2) - sum(p1.^2)];
            if abs(det(A)) < 1e-12
                c = [0 0]; r = inf; a1 = 0; da = 0; return;
            end
            c = (A\b)';
            r = norm(p1 - c);
            a1 = atan2(p1(2)-c(2), p1(1)-c(1));
            am = atan2(pm(2)-c(2), pm(1)-c(1));
            a2 = atan2(p2(2)-c(2), p2(1)-c(1));
            ccw = @(a) mod(a - a1, 2*pi);           % angle travelled CCW from p1
            if ccw(am) <= ccw(a2)
                da = ccw(a2);                       % CCW passes pm
            else
                da = ccw(a2) - 2*pi;                % CW
            end
        end
    end
end
