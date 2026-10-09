classdef scene_glyph
    %SCENE_GLYPH  3D models of the scene objects for listener_space_axes.
    %   Seen from above they reproduce the former 2D icons (head from
    %   above, schematic loudspeaker), from the side/3D views they read as
    %   a person and a loudspeaker. Each model is one patch struct
    %   (Vertices, Faces, FaceVertexCData) facing +x, centred so that
    %   mean(Vertices(:,1:2)) = [0 0] (the drag bookkeeping uses the vertex
    %   mean as the object position), ear/acoustic centre at z = 0.
    %   Shading is baked into the face colours (no light object: a light
    %   would also shade the renderer overlays).

    methods (Static)
        function G = head(R)
            % Listener: head (hair on top/back, face in front), nose, ears,
            % neck and shoulders. R = head radius.
            hair  = [77 41 14]/255;
            skin  = [255 206 180]/255;
            shirt = [37 160 217]/255;
            G = scene_glyph.empty_glyph();
            G = scene_glyph.add(G, scene_glyph.ellipsoid([-R/7 0 -2.0*R], [0.7 1.7 0.9]*R, 20), shirt);
            G = scene_glyph.add(G, scene_glyph.ellipsoid([-0.1*R 0 -1.1*R], [0.45 0.45 0.6]*R, 12), skin);
            [V, F, N] = scene_glyph.ellipsoid([0 0 0.1*R], [1.05 0.95 1.15]*R, 20);
            Cf = scene_glyph.face_centers(V, F);
            isHair = Cf(:,3) > 0.25*R | Cf(:,1) < -0.35*R;
            col = repmat(skin, size(F,1), 1);
            col(isHair,:) = repmat(hair, nnz(isHair), 1);
            G = scene_glyph.add(G, {V, F, N}, col);
            G = scene_glyph.add(G, scene_glyph.ellipsoid([0.95*R 0 -0.1*R], [0.3 0.13 0.22]*R, 10), skin);
            G = scene_glyph.add(G, scene_glyph.ellipsoid([-0.05*R  0.95*R 0], [0.22 0.15 0.32]*R, 10), skin);
            G = scene_glyph.add(G, scene_glyph.ellipsoid([-0.05*R -0.95*R 0], [0.22 0.15 0.32]*R, 10), skin);
            G = scene_glyph.finish(G);
        end

        function G = loudspeaker(R)
            % Schematic loudspeaker: cabinet + flared horn/cone, membrane.
            % Same footprint as the 2D symbol (body -1.8R..-R, horn to 0,
            % mouth half-width 1.5R).
            body  = [0.2 0.2 0.2];
            horn  = [0.5 0.5 0.5];
            memb  = [0.12 0.12 0.12];
            G = scene_glyph.empty_glyph();
            G = scene_glyph.add(G, scene_glyph.box([-1.8*R -R -R], [-R R R]), body);
            G = scene_glyph.add(G, scene_glyph.frustum(-R, 0, R, 1.5*R, 24), horn);
            G = scene_glyph.add(G, scene_glyph.disc(-0.12*R, 1.38*R, 24), memb);
            G = scene_glyph.finish(G);
        end
    end

    methods (Static, Access = private)
        function G = empty_glyph()
            G = struct('V', zeros(0,3), 'F', zeros(0,4), 'C', zeros(0,3));
        end

        function G = add(G, part, col)
            % part = {V, F, N} (N: outward face normals)
            [V, F, N] = part{:};
            if size(col,1) == 1, col = repmat(col, size(F,1), 1); end
            L = [0.35 -0.45 0.82]; L = L/norm(L);       % light from above
            shade = 0.55 + 0.45*max(0, N*L');
            G.F = [G.F; F + size(G.V,1)];
            G.V = [G.V; V];
            G.C = [G.C; col.*shade];
        end

        function G = finish(G)
            G.V(:,1:2) = G.V(:,1:2) - mean(G.V(:,1:2), 1);
        end

        function C = face_centers(V, F)
            C = (V(F(:,1),:) + V(F(:,2),:) + V(F(:,3),:) + V(F(:,4),:))/4;
        end

        function [V, F, N] = ellipsoid(c, r, n)
            [X, Y, Z] = sphere(n);
            fv = surf2patch(c(1) + r(1)*X, c(2) + r(2)*Y, c(3) + r(3)*Z);
            V = fv.vertices; F = fv.faces;
            N = (scene_glyph.face_centers(V, F) - c)./r.^2;   % gradient
            N = N./max(vecnorm(N, 2, 2), eps);
            if nargout < 2, V = {V, F, N}; end
        end

        function out = box(lo, hi)
            [x, y, z] = ndgrid([lo(1) hi(1)], [lo(2) hi(2)], [lo(3) hi(3)]);
            V = [x(:) y(:) z(:)];
            F = [1 3 7 5; 2 6 8 4; 1 5 6 2; 3 4 8 7; 1 2 4 3; 5 7 8 6];
            N = [-1 0 0; 1 0 0; 0 -1 0; 0 1 0; 0 0 -1; 0 0 1];
            out = {V, F, N};
        end

        function out = frustum(x1, x2, r1, r2, n)
            % Open cone along +x, radius r1 at x1 to r2 at x2
            t = linspace(0, 2*pi, n+1)'; t(end) = [];
            V = [x1*ones(n,1), r1*cos(t), r1*sin(t);
                 x2*ones(n,1), r2*cos(t), r2*sin(t)];
            k = (1:n)'; k2 = mod(k, n) + 1;
            F = [k, k2, k2+n, k+n];
            tc = t + pi/n;                                % face mid-angle
            s = (r2 - r1)/(x2 - x1);
            N = [-s*ones(n,1), cos(tc), sin(tc)];
            N = N./vecnorm(N, 2, 2);
            out = {V, F, N};
        end

        function out = disc(x, r, n)
            % Disc facing +x (quads degenerate to triangles at the centre)
            t = linspace(0, 2*pi, n+1)'; t(end) = [];
            V = [x 0 0; x*ones(n,1), r*cos(t), r*sin(t)];
            k = (2:n+1)'; k2 = [3:n+1, 2]';
            F = [ones(n,1), k, k2, ones(n,1)];
            N = repmat([1 0 0], n, 1);
            out = {V, F, N};
        end
    end
end
