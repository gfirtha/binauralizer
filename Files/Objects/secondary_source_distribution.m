classdef secondary_source_distribution < handle
    %SECONDARY_SOURCE_DISTRIBUTION Summary of this class goes here
    %   Detailed explanation goes here

    properties
        type
        loudspeakers
        dimensionality
        center
        x0
        n0
        dA
        isClosed
    end

    methods
        function obj = secondary_source_distribution(loudspeakers, center)
            %SECONDARY_SOURCE_DISTRIBUTION Construct an instance of this class
            %   Detailed explanation goes here
            obj.loudspeakers = loudspeakers;
            obj.update_geometry(center)
        end

        function obj = update_geometry(obj, center)
            obj.x0 = cell2mat(cellfun(@(x) x.position,    obj.loudspeakers, 'UniformOutput', false)');
            obj.n0 = cell2mat(cellfun(@(x) x.orientation, obj.loudspeakers, 'UniformOutput', false)');
            obj.center = center;

            if size(obj.x0,2) == 2
                obj.dimensionality = 2;
            else
                obj.dimensionality = 3;
            end

            % Degenerate guards
            if size(obj.x0,1) < 3 || size(unique(obj.x0,'rows'),1) < 3
                obj.isClosed = false;
                return;
            end

            % Classify geometry rank (point / line / plane)
            [geomType, basis, mu] = obj.classify_geometry(obj.x0);
            switch geomType
                case 'plane'  % rank 2 → project and do point-in-polygon
                    X2 = (obj.x0 - mu) * basis(:,1:2);  % local 2D coords
                    c2 = (obj.center - mu) * basis(:,1:2);
                    % Build boundary (allows mild concavity). Fallback to convex hull if needed.
                    try
                        idx = boundary(X2(:,1), X2(:,2), 0.9);
                    catch
                        idx = convhull(X2(:,1), X2(:,2));
                    end

                    % If boundary couldn’t form a loop, not closed
                    if numel(idx) < 3
                        obj.isClosed = false;
                        return;
                    end
                    pg = polyshape(X2(idx,1), X2(idx,2), 'Simplify', true);
                    obj.isClosed = isinterior(pg, c2(1), c2(2));
                    R = sqrt(sum((obj.x0-mu).^2,2));
                    tol = 1e-6;
                    if max(R) - min(R) <= tol
                        obj.type = 'circular';
                    else
                        obj.type = 'general_enclosing';
                    end
                case 'line'   % rank 1 → no area → not closed
                    obj.isClosed = false;
                    obj.type = 'linear';
                case 'point'  % rank 0 → not closed
                    obj.isClosed = false;
            end

            switch obj.isClosed
                case true
                    obj.dA = mean([sqrt(sum((obj.x0-circshift(obj.x0,-1)).^2,2)), sqrt(sum((obj.x0-circshift(obj.x0,1)).^2,2))],2);
                case false
                    obj.dA = min( sqrt(sum((obj.x0-circshift(obj.x0,-1)).^2,2)), sqrt(sum((obj.x0-circshift(obj.x0,1)).^2,2)) );

            end
        end

        function [P, t_hit] = get_ray_intersection(obj, origin, v)
            %GET_RAY_INTERSECTION Intersect SSD contour with a ray.
            %   [P, t] = obj.get_ray_intersection(origin, v)
            %   P     : K×D intersection points (K=0,1,2), D=2 or 3
            %   t_hit : K×1 ray parameters where P = origin + t*v, t>=0

            P = []; t_hit = [];
            if isempty(obj.x0) || size(obj.x0,1) < 2, return; end

            o = origin(:);
            v = v(:);
            nv = norm(v);
            if nv == 0, return; end
            v = v / nv;

            D = size(obj.x0,2);
            tolDist   = 1e-9;
            tolParam  = 1e-12;

            % Geometry classification (for circular/plane)
            [geomType, basis, mu] = obj.classify_geometry(obj.x0);

            switch lower(string(obj.type))
                case "linear"
                    % ---- Build best-fit segment from samples ----
                    [c0, d, umin, umax] = obj.line_segment_from_points(obj.x0); % d is unit
                    c0 = c0(:); d = d(:);

                    % Check parallelism WITHOUT cross (works in any D):
                    % v ⟂-component to d small => parallel/anti-parallel
                    v_perp = v - d * (d.'*v);
                    isParallel = norm(v_perp) <= 1e-12;

                    if isParallel
                        % Distance from ray origin to the infinite line
                        w0 = o - c0;
                        dist = norm(w0 - d*(d.'*w0));
                        if dist > tolDist
                            % Parallel, disjoint
                            return;
                        end
                        % Colinear: earliest hit on the segment in forward ray direction
                        u0    = d.'*(o - c0);     % ray origin's u along the line
                        alpha = d.'*v;            % = ±1 (since both unit), sign = direction
                        % If origin already on the segment, take t=0
                        if u0 >= umin - 1e-12 && u0 <= umax + 1e-12
                            P = o.'; t_hit = 0; return;
                        end
                        % Otherwise, see if moving forward reaches the segment
                        if alpha > 0
                            if u0 < umin
                                t = (umin - u0)/alpha; P = (o + t*v).'; t_hit = t; return;
                            else
                                % u0 > umax and moving increasing => segment is behind
                                return;
                            end
                        else % alpha < 0
                            if u0 > umax
                                t = (umax - u0)/alpha; P = (o + t*v).'; t_hit = t; return;
                            else
                                % u0 < umin and moving decreasing => segment is behind
                                return;
                            end
                        end
                    else
                        % General case: solve [v -d][t;u] = (c0 - o), then verify
                        M   = [v, -d];
                        rhs = (c0 - o);
                        xu  = M\rhs;  % 2×1
                        t = xu(1); u = xu(2);
                        p_ray = o + t*v;
                        p_seg = c0 + u*d;

                        if norm(p_ray - p_seg) <= tolDist && t >= -tolParam && ...
                                u >= umin - 1e-12 && u <= umax + 1e-12
                            P = p_ray.'; t_hit = max(t,0);
                        end
                        return;
                    end

                case "circular"
                    % ---- Fit circle and intersect with ray ----
                    if D == 2
                        % Work fully in 2D
                        [c2, R] = obj.fit_circle_2d(obj.x0);
                        o2 = o; v2 = v;
                        delta = o2 - c2(:);
                        a = dot(v2,v2); % = 1
                        b = 2*dot(delta,v2);
                        c = dot(delta,delta) - R^2;
                        disc = b^2 - 4*a*c;

                        if disc < -1e-12, return; end
                        if abs(disc) <= 1e-12
                            t = -b/(2*a);
                            if t >= -tolParam, P = (o2 + t*v2).'; t_hit = max(t,0); end
                            return;
                        end
                        rt = sqrt(max(0,disc));
                        t1 = (-b - rt)/(2*a);
                        t2 = (-b + rt)/(2*a);
                        if t1 >= -tolParam, P = [P; (o2 + t1*v2).']; t_hit = [t_hit; max(t1,0)]; end
                        if t2 >= -tolParam, P = [P; (o2 + t2*v2).']; t_hit = [t_hit; max(t2,0)]; end
                        % Dedup near-tangent
                        if size(P,1) == 2 && norm(P(1,:) - P(2,:)) < 1e-10
                            P = P(1,:); t_hit = t_hit(1);
                        end
                        return;

                    elseif D == 3
                        % Project to the SSD plane, fit circle there
                        if ~strcmp(geomType,'plane') || size(basis,2) < 2, return; end
                        e1 = basis(:,1); e2 = basis(:,2);
                        n  = cross(e1, e2); n = n / norm(n);

                        X2 = (obj.x0 - mu) * [e1 e2];           % N×2
                        [c2, R] = obj.fit_circle_2d(X2);         % center in local coords
                        c3 = mu(:) + e1*c2(1) + e2*c2(2);        % center in 3D

                        denom = dot(v, n);
                        if abs(denom) > 1e-12
                            % Intersect ray with plane at once point
                            t = dot((c3 - o), n) / denom;
                            if t < -tolParam, return; end
                            p = o + t*v;
                            if abs(norm(p - c3) - R) <= 1e-6
                                P = p.'; t_hit = max(t,0);
                            end
                            return;
                        else
                            % Ray parallel to plane normal. Check coplanarity.
                            if abs(dot(o - mu(:), n)) > 1e-8, return; end
                            % Coplanar: 2D intersection in plane coords
                            O2 = [(o - c3).'*e1, (o - c3).'*e2];
                            V2 = [v.'*e1, v.'*e2];

                            a = dot(V2,V2);
                            b = 2*dot(O2,V2);
                            c = dot(O2,O2) - R^2;

                            disc = b^2 - 4*a*c;
                            if disc < -1e-12, return; end
                            if abs(disc) <= 1e-12
                                t = -b/(2*a);
                                if t >= -tolParam, P = (o + t*v).'; t_hit = max(t,0); end
                                return;
                            end
                            rt = sqrt(max(0,disc));
                            t1 = (-b - rt)/(2*a);
                            t2 = (-b + rt)/(2*a);
                            if t1 >= -tolParam, P = [P; (o + t1*v).']; t_hit = [t_hit; max(t1,0)]; end
                            if t2 >= -tolParam, P = [P; (o + t2*v).']; t_hit = [t_hit; max(t2,0)]; end
                            if size(P,1) == 2 && norm(P(1,:) - P(2,:)) < 1e-10
                                P = P(1,:); t_hit = t_hit(1);
                            end
                            return;
                        end
                    else
                        % D not supported for circular (only 2D/3D make sense)
                        return;
                    end

                otherwise
                    % Not (yet) supported geometry
                    return;
            end
        end




    end
    methods (Static, Access = private)
        function [geomType, basis, mu] = classify_geometry(P)
            % Determine 'point' / 'line' / 'plane' and return ambient-space basis.
            % P is N×D (N points, D dims). basis is D×r (columns = principal dirs).

            mu = mean(P,1);
            X  = P - mu;            % N×D

            % SVD: X = U*S*V', columns of V are principal directions in ambient D-space
            [~,S,V] = svd(X, 'econ');    % U: N×r, S: r×r, V: D×r
            s = diag(S);
            if isempty(s)
                geomType = 'point';
                basis    = eye(size(P,2),1);
                return;
            end

            % Rank estimate
            tol_rank = 1e-8;
            r = sum(s > s(1)*tol_rank);

            if r <= 0
                geomType = 'point';
            elseif r == 1
                geomType = 'line';
            else
                geomType = 'plane';
            end

            % Return ambient-space basis (principal directions)
            basis = V(:, max(1,r));          % ensure at least 1 column
            if size(basis,2) < min(2,size(V,2))
                basis = V(:,1:min(2,size(V,2)));
            end
        end
        function [c2, R] = fit_circle_2d(X)
            %FIT_CIRCLE_2D  Algebraic circle fit (Taubin-like), X is N×2.
            %   Returns center c2 = [cx, cy], radius R.
            x = X(:,1); y = X(:,2);
            mx = mean(x); my = mean(y);
            x = x - mx; y = y - my;  % center to improve conditioning

            Z = [x.^2 + y.^2, x, y, ones(size(x))];
            [~,~,V] = svd(Z,0);
            p = V(:,end);

            if abs(p(1)) < 1e-14
                % Fallback: linear LS
                A = [2*x, 2*y, ones(size(x))];
                b = x.^2 + y.^2;
                sol = A\b;
                cx = sol(1); cy = sol(2);
                R  = sqrt(max(0, sol(3) + cx^2 + cy^2));
                c2 = [cx + mx, cy + my];
                return;
            end

            A = p(1); B = p(2); C = p(3); D = p(4);
            cx = -B/(2*A);
            cy = -C/(2*A);
            R  = sqrt(max(0, (B^2 + C^2)/(4*A^2) - D/A));

            % undo mean shift
            c2 = [cx + mx, cy + my];
        end

        function [c0, d, umin, umax] = line_segment_from_points(P)
            %LINE_SEGMENT_FROM_POINTS  Best-fit line segment from samples P (N×D).
            %   c0   : a point on the line (the mean)
            %   d    : unit direction (first principal component)
            %   umin, umax : scalar extents covering the sample span along d
            mu = mean(P,1);
            X  = P - mu;
            [~,~,V] = svd(X, 'econ');
            d = V(:,1); d = d/norm(d);
            u = X * d;               % scalar projection of each point on d
            umin = min(u); umax = max(u);
            c0 = mu(:);
        end

    end
    methods (Static)
        function [x0, n0] = load_eps_as_ssd(filename, N, R_target)
            % Read EPS → largest closed contour → resample to N points →
            % inward normals → (optional) scale so max pairwise distance == R_target.
            %
            % If you want a different scaling meaning:
            %  - max adjacent spacing  == R_target: after P is built, set
            %      scale = R_target / (arcLength(P)/N);
            %  - "radius" like circle: set max pairwise distance == 2*R_target.

            if nargin < 2 || N < 3 || mod(N,1) ~= 0
                error('load_eps_as_ssd: N must be an integer >= 3');
            end

            [V, closed] = secondary_source_distribution.read_eps_contour(filename);
            if isempty(V) || size(V,1) < 3
                error('No valid contour extracted from EPS.');
            end
            if ~closed
                if any(V(1,:) ~= V(end,:)), V = [V; V(1,:)]; end
            end
            if norm(V(1,:) - V(end,:)) < 1e-12
                V = V(1:end-1,:);
            end

            V = secondary_source_distribution.sanitize_polyline(V);

            P = secondary_source_distribution.sample_polygon_equidistant(V, N);

            P = P - mean(P,1);
            % ---- optional scaling -------------------------------------------------
            if nargin >= 3 && ~isempty(R_target) && isfinite(R_target) && R_target > 0
                % NEW (no toolboxes)
                K  = convhull(P(:,1), P(:,2));   % base MATLAB
                H  = P(K,:);                     % hull points (ordered, last == first)
                dmax = secondary_source_distribution.max_pairwise_distance(H);
                s = R_target / dmax;
                P = P * s ;         % scale about centroid
            end
            % ----------------------------------------------------------------------

            n = secondary_source_distribution.compute_inward_normals(P);
            x0 = P; n0 = n;
        end

        function [V, isClosed] = read_eps_contour(fn)
            % Minimal EPS path reader with multiple-subpath handling.
            % Supports: moveto/lineto/curveto/closepath and *relative* rmoveto/rlineto.
            % Curves are tessellated to 16 segments. Ignores transforms.
            txt = fileread(fn);

            tokens = regexp(txt, '([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)|([A-Za-z][A-Za-z]*)', 'match');
            stack = []; paths = {}; cur = []; isClosed = false;

            function flushCur()
                if ~isempty(cur)
                    % Keep as its own subpath; close if already closed-ish
                    if size(cur,1) >= 3
                        if norm(cur(1,:) - cur(end,:)) > 1e-9
                            % leave open; we’ll ignore open subpaths later
                        end
                        paths{end+1} = cur; %#ok<AGROW>
                    end
                    cur = [];
                end
            end

            i = 1;
            while i <= numel(tokens)
                t = tokens{i};
                val = str2double(t);
                if ~isnan(val)
                    stack(end+1) = val; %#ok<AGROW>
                else
                    op = lower(t);
                    switch op
                        case {'m','moveto'}
                            % starting a new subpath: flush previous
                            flushCur();
                            if numel(stack) >= 2
                                y = stack(end); x = stack(end-1); stack(end-1:end) = [];
                                cur = [cur; x y];
                            else
                                stack = [];
                            end

                        case {'rmoveto'}
                            flushCur();
                            if numel(stack) >= 2
                                dy=stack(end); dx=stack(end-1); stack(end-1:end)=[];
                                if isempty(cur), cur = [0 0]; end
                                cur = [cur; cur(end,:)+[dx dy]];
                            else
                                stack = [];
                            end

                        case {'l','lineto'}
                            if numel(stack) >= 2 && ~isempty(cur)
                                y = stack(end); x = stack(end-1); stack(end-1:end) = [];
                                cur = [cur; x y];
                            else
                                stack = [];
                            end

                        case {'rlineto'}
                            if numel(stack) >= 2 && ~isempty(cur)
                                dy=stack(end); dx=stack(end-1); stack(end-1:end)=[];
                                cur = [cur; cur(end,:)+[dx dy]];
                            else
                                stack = [];
                            end

                        case {'c','curveto'}
                            if numel(stack) >= 6 && ~isempty(cur)
                                y3=stack(end); x3=stack(end-1);
                                y2=stack(end-2); x2=stack(end-3);
                                y1=stack(end-4); x1=stack(end-5);
                                stack(end-5:end) = [];
                                p0 = cur(end,:);
                                seg = secondary_source_distribution.sample_bezier_cubic(p0,[x1 y1],[x2 y2],[x3 y3],16);
                                cur = [cur; seg(2:end,:)];
                            else
                                stack = [];
                            end

                        case {'h','cp','closepath'}
                            if ~isempty(cur)
                                % make sure it’s explicitly closed
                                if norm(cur(1,:) - cur(end,:)) > 1e-9
                                    cur = [cur; cur(1,:)];
                                end
                                paths{end+1} = cur; %#ok<AGROW>
                                cur = [];
                            end

                        otherwise
                            % ignore gsave/grestore/stroke/fill etc.
                    end
                end
                i = i + 1;
            end
            flushCur();  % leftover open subpath

            % Keep only closed loops
            polys = {};
            for k = 1:numel(paths)
                P = paths{k};
                if size(P,1) >= 3 && norm(P(1,:) - P(end,:)) <= 1e-9
                    polys{end+1} = secondary_source_distribution.sanitize_polyline(P); %#ok<AGROW>
                end
            end
            if isempty(polys), V = []; isClosed = false; return; end

            % Pick the largest area loop
            areas = cellfun(@(Q) polyarea(Q(:,1),Q(:,2)), polys);
            [~, j] = max(areas);
            V = polys{j};
            isClosed = true;
        end

        function V = sanitize_polyline(V)
            % Remove consecutive duplicates and tiny segments that create clusters.
            if isempty(V), return; end
            % Ensure closed
            closed = norm(V(1,:) - V(end,:)) <= 1e-12;
            if ~closed, V = [V; V(1,:)]; closed = true; end

            % Remove consecutive duplicates
            keep = true(size(V,1),1);
            keep(2:end) = any(abs(diff(V,1,1)) > 1e-12, 2);
            V = V(keep,:);

            % Remove very short segments
            if size(V,1) > 2
                seg = sqrt(sum(diff(V,1,1).^2, 2));
                L   = sum(seg);
                thr = max(1e-6, 1e-4 * (L / max(1,size(V,1))));  % relative threshold
                i = 1;
                while i <= size(V,1)-1
                    if norm(V(i+1,:) - V(i,:)) < thr
                        V(i+1,:) = [];      % collapse
                    else
                        i = i + 1;
                    end
                end
                % re-close if needed
                if norm(V(1,:) - V(end,:)) > 1e-12
                    V = [V; V(1,:)];
                end
            end
        end

        function P = sample_polygon_equidistant(V, N)
            if norm(V(1,:) - V(end,:)) > 1e-12
                V = [V; V(1,:)];
            end
            seg = sqrt(sum(diff(V,1,1).^2, 2));
            s   = [0; cumsum(seg)];
            L   = s(end);
            t   = linspace(0, L, N+1).'; t(end) = [];
            P   = zeros(N, size(V,2));
            for k = 1:N
                dk = t(k);
                i  = find(s <= dk, 1, 'last');
                if i == numel(s), i = numel(s)-1; end
                ds = dk - s(i);
                if seg(i) < 1e-15
                    P(k,:) = V(i,:);
                else
                    alpha = ds / seg(i);
                    P(k,:) = (1-alpha)*V(i,:) + alpha*V(i+1,:);
                end
            end
        end

        function n = compute_inward_normals(P)
            N = size(P,1);
            if norm(P(1,:) - P(end,:)) > 1e-12, P = [P; P(1,:)]; end
            cent = mean(P(1:end-1,:), 1);
            n = zeros(N, size(P,2));
            for i = 1:N
                i1 = i; i2 = i+1;
                e  = P(i2,:) - P(i1,:);
                if norm(e) < 1e-15, n_e = [0 0];
                else, t = e / norm(e); n_e = [t(2), -t(1)];
                end
                if i == 1, ep = P(1,:) - P(end-1,:); else, ep = P(i,:) - P(i-1,:); end
                if norm(ep) < 1e-15, n_ep = [0 0];
                else, tp = ep / norm(ep); n_ep = [tp(2), -tp(1)];
                end
                v = n_e + n_ep; if norm(v) < 1e-15, v = n_e; end
                if norm(v) > 0, v = v / norm(v); end
                if dot(v, cent - P(i,:)) < 0, v = -v; end
                n(i,:) = v;
            end
        end

        function Q = sample_bezier_cubic(p0, p1, p2, p3, K)
            t = linspace(0,1,K).';
            Q = (1-t).^3 .* p0 + 3*(1-t).^2 .* t .* p1 + 3*(1-t) .* t.^2 .* p2 + t.^3 .* p3;
        end
        function dmax = max_pairwise_distance(X)
            % Toolbox-free max pairwise distance (works for 2D/3D).
            % O(n^2) time, O(n) memory – fine for a few thousand points.
            n = size(X,1);
            dmax2 = 0;
            for i = 1:n
                d2 = sum( (X - X(i,:)).^2, 2 );
                m  = max(d2);
                if m > dmax2, dmax2 = m; end
            end
            dmax = sqrt(dmax2);
        end
    end
    % --- put these in the class body ---

    methods (Static)   % public static (remove Access=private)
     function [hit, tmin, tnear] = rayPolyIntersect(p, d, P)
        % Ray p + t d (t>=0) vs polygon P (Mx2, closed: P(end,:)=P(1,:))
        % Returns:
        %   hit   : logical, true if intersection exists
        %   tmin  : smallest non-negative intersection distance along the ray
        %   tnear : if no intersection, distance along the ray to the closest
        %           point on the polygon boundary (else equals tmin)
        %
        % Notes:
        %  - Distances are in the same units as P and p.
        %  - We internally normalize d so t is a true Euclidean distance.

        if size(p,2) ~= 2, p = p(:).'; end
        if size(d,2) ~= 2, d = d(:).'; end

        % Normalize ray direction so t is a metric distance.
        nd = norm(d);
        if nd < eps
            hit = false; tmin = NaN; tnear = NaN; return;
        end
        d = d / nd;

        M   = size(P,1);
        tol = 1e-12;

        hit   = false;
        tmin  = inf;

        % We also track the closest distance ray↔segment to provide tnear
        bestDist2 = inf;
        tnear     = inf;

        for j = 1:M-1
            a = P(j,:); 
            b = P(j+1,:);
            e = b - a;

            % --- 1) Exact intersection test (ray vs segment) ---------------
            den = secondary_source_distribution.cross2(d, e);
            if abs(den) >= tol
                ap = a - p;
                t  = secondary_source_distribution.cross2(ap, e) / den;
                u  = secondary_source_distribution.cross2(ap, d) / den;

                if t >= -tol && u >= -tol && u <= 1+tol
                    % record intersection
                    tpos = max(t,0); % force non-negative
                    if tpos < tmin
                        tmin = tpos;
                        hit  = true;
                    end
                end
            end

            % --- 2) Nearest approach (ray ↔ segment), used if no hit -------
            [tCand, dist2Cand] = secondary_source_distribution.closestRaySegmentT(p, d, a, b);
            if dist2Cand < bestDist2
                bestDist2 = dist2Cand;
                tnear     = tCand;
            end
        end

        % For consistency: if we did hit, the nearest "approach" along ray is the hit.
        if hit
            tnear = tmin;
        end
    end

    function z = cross2(a,b)
        % 2D scalar cross product
        z = a(1)*b(2) - a(2)*b(1);
    end
     function [t, dist2] = closestRaySegmentT(p, d, a, b)
        % Distance along the ray r(t) = p + t d, t>=0, to the closest point
        % on segment s(u) = a + u (b-a), u in [0,1].  Assumes d is unit.
        % Returns:
        %   t     : argmin_{t>=0, u∈[0,1]} || p + t d - (a + u e) ||, in distance units
        %   dist2 : the squared minimal distance at (t,u).

        e  = b - a;
        ee = dot(e,e);
        r  = p - a;

        % Handle degenerate (zero-length) segment
        if ee < 1e-20
            % Reduce to point a
            t    = max(0, dot(d, a - p));      % projection onto ray (>=0)
            qray = p + t*d;
            dist2 = sum((qray - a).^2);
            return;
        end

        de = dot(d, e);
        dr = dot(d, r);
        er = dot(e, r);

        denom = ee - de*de;     % since ||d||=1
        if abs(denom) > 1e-20
            % Unconstrained optimum
            t0 = (de*er - ee*dr) / denom;
            u0 = (de*dr - er)    / denom;

            % Clamp to constraints
            t = max(t0, 0);
            u = min(max(u0, 0), 1);

            % If clamping changed either variable, recompute the other at the clamp
            if t ~= t0
                % for fixed t, best u is projection of (p + t d - a) onto e
                u = min(max( dot(e, (p + t*d - a)) / ee, 0), 1);
            end
            if u ~= u0
                % for fixed u, best t is projection of (a + u e - p) onto d (>=0)
                t = max( dot(d, (a + u*e - p)), 0 );
            end
        else
            % Ray is (near) parallel to segment → choose the best among u=0/1 and t>=0
            tA = max(0, dot(d, a - p));
            tB = max(0, dot(d, b - p));

            qA = p + tA*d;  dist2A = sum((qA - a).^2);
            qB = p + tB*d;  dist2B = sum((qB - b).^2);

            if dist2A <= dist2B
                t = tA; dist2 = dist2A;
            else
                t = tB; dist2 = dist2B;
            end
            return;
        end

        % Compute minimal distance at (t,u)
        qray = p + t*d;
        qseg = a + u*e;
        dist2 = sum((qray - qseg).^2);
    end
    end





end



