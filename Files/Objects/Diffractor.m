classdef Diffractor < handle
    % DIFFRACTOR  Rigid sphere acoustic scattering model.
    %
    % Computes the total (incident + scattered) pressure field around a
    % rigid sphere in 2D, for point-source or plane-wave incidence.
    % Uses the classical Legendre series expansion.
    %
    % This class is standalone: it knows nothing about the GUI, renderers,
    % or the scene. It can be used by:
    %   - sound_scene_simulator (field maps with head effect)
    %   - CTC renderers (plant matrix construction)
    %   - standalone scripts (transfer function analysis)
    %
    % Usage:
    %   d = Diffractor();
    %   d.scatterer_pos    = [0, 0];
    %   d.scatterer_radius = 0.0875;   % m
    %   d.source_pos       = [1, 0];
    %   d.source_type      = 'point';
    %
    %   % Total field at evaluation points
    %   p = d.evaluate(f, x_eval);          % [Npts x Nf]
    %
    %   % Left/right ear transfer functions
    %   H = d.evaluate_ears(f, head_orientation);  % [2 x Nf]
    %
    %   % Free-field reference (no scatterer)
    %   p0 = d.evaluate_freefield(f, x_eval);

    properties
        scatterer_pos    double = [0, 0, 0]    % sphere center [x, y, z]
        scatterer_radius (1,1) double = 0.0875 % sphere radius (m)
        source_pos       double = [1, 0, 0]    % source position [x, y, z]
        source_type      char = 'point'        % 'point' | 'plane'
        source_dir       double = [1, 0, 0]    % plane wave direction (source_type='plane')
        N_order          (1,1) double = -1     % series truncation (-1 = auto)
    end

    methods
        function p = evaluate(obj, f, x_eval)
            % EVALUATE  Total pressure (incident + scattered) at eval points.
            %
            %   f      : [1 x Nf] or scalar — frequencies (Hz)
            %   x_eval : [Npts x 2] — evaluation positions
            %   p      : [Npts x Nf] — complex pressure
            %
            % Points inside the sphere return 0.

            c = audioapp.util.Physics.speedOfSound();
            a = obj.scatterer_radius;
            f = f(:).';   % row
            Nf = numel(f);
            Npts = size(x_eval, 1);

            % Source-to-scatterer geometry
            x_src = obj.source_pos;
            x_scat = obj.scatterer_pos;

            % Evaluation points relative to scatterer center
            dx = bsxfun(@minus, x_eval, x_scat(1:size(x_eval,2)));
            R = sqrt(sum(dx.^2, 2));           % [Npts x 1]

            % Incidence angle for each eval point:
            % angle between (source → scatterer) direction and (scatterer → eval) direction
            switch lower(obj.source_type)
                case 'point'
                    v_inc = x_scat(1:size(x_eval,2)) - x_src(1:size(x_eval,2));
                case 'plane'
                    v_inc = obj.source_dir(1:size(x_eval,2));
            end
            v_inc_n = v_inc / max(norm(v_inc), 1e-15);

            % Angle theta: 0 = shadow side, pi = illuminated side
            cos_theta = sum(bsxfun(@times, dx, v_inc_n), 2) ./ max(R, 1e-15);
            cos_theta = max(-1, min(1, cos_theta));  % clamp

            % Source distance (for point source)
            R_src = norm(x_src(1:size(x_eval,2)) - x_scat(1:size(x_eval,2)));

            % Frequency loop
            k = 2*pi * f / c;                 % [1 x Nf]

            % Auto truncation order
            if obj.N_order < 0
                Nmax = ceil(max(k) * a) + 15;
            else
                Nmax = obj.N_order;
            end

            % Preallocate
            p = zeros(Npts, Nf);

            % Mask: points inside sphere → 0
            outside = R > a * 0.999;

            for fi = 1:Nf
                kk = k(fi);
                if kk < 1e-10
                    % DC: free-field only
                    p(:, fi) = obj.ff_greens(R, R_src, 0);
                    continue;
                end

                ka = kk * a;
                kR = kk * R(outside);            % [Nout x 1]

                switch lower(obj.source_type)
                    case 'point'
                        kR_src = kk * R_src;
                        p_out = obj.rigid_sphere_point(kk, ka, kR, kR_src, ...
                            cos_theta(outside), Nmax);
                    case 'plane'
                        p_out = obj.rigid_sphere_plane(kk, ka, kR, ...
                            cos_theta(outside), Nmax);
                end

                p(outside, fi) = p_out;
            end
        end

        function H = evaluate_ears(obj, f, head_orientation, ear_angle)
            % EVALUATE_EARS  Transfer function at left and right ear positions.
            %
            %   f                : [1 x Nf] — frequencies (Hz)
            %   head_orientation : [1 x 2] or [1 x 3] — unit vector pointing "nose" direction
            %   ear_angle        : scalar (rad), default pi/2
            %   H                : [2 x Nf] — H(1,:) = left ear, H(2,:) = right ear
            %
            %   Uses the standard rigid sphere series evaluated at r = a,
            %   vectorized over frequency.

            if nargin < 4, ear_angle = pi/2; end

            c = audioapp.util.Physics.speedOfSound();
            a = obj.scatterer_radius;
            f = f(:).';          % row
            Nf = numel(f);
            k = 2*pi * f / c;   % [1 x Nf]

            % ── Ear positions and incidence angles ──
            % Work in 3D internally
            if numel(head_orientation) < 3
                n_hat3 = [head_orientation(:)', 0];
            else
                n_hat3 = head_orientation(1:3);
            end
            n_hat3 = n_hat3 / max(norm(n_hat3), 1e-15);

            % Perpendicular direction (left ear) using cross product with up vector
            up = [0, 0, 1];
            perp3 = cross(n_hat3, up);
            if norm(perp3) < 1e-10
                % Degenerate case: nose pointing straight up/down
                perp3 = cross(n_hat3, [1, 0, 0]);
            end
            perp3 = perp3 / norm(perp3);

            x_scat = obj.scatterer_pos(1:3);
            if numel(obj.scatterer_pos) < 3
                x_scat = [obj.scatterer_pos(:)', 0];
            end
            x_src = obj.source_pos(1:3);
            if numel(obj.source_pos) < 3
                x_src = [obj.source_pos(:)', 0];
            end

            % Incidence direction (3D)
            switch lower(obj.source_type)
                case 'point'
                    v_inc = x_scat - x_src;
                case 'plane'
                    sd = obj.source_dir;
                    if numel(sd) < 3, sd = [sd(:)', 0]; end
                    v_inc = sd(1:3);
            end
            v_inc_n = v_inc / max(norm(v_inc), 1e-15);

            % cos(theta) for each ear
            cos_theta = zeros(2, 1);
            ear_dirs = [perp3; -perp3];   % left ear, right ear
            for ear = 1:2
                cos_theta(ear) = dot(ear_dirs(ear,:), v_inc_n);
            end

            % ── Series truncation ──
            ka_max = max(k) * a;
            if obj.N_order < 0
                Nmax = ceil(ka_max) + 15;
            else
                Nmax = obj.N_order;
            end

            % ── Precompute Legendre polynomials ──
            Pn_ears = zeros(2, Nmax+1);
            for n = 0:Nmax
                for ear = 1:2
                    Pn_ears(ear, n+1) = getLegendre(n, 0, cos_theta(ear));
                end
            end

            % ── Vectorized series ──
            ka = k * a;            % [1 x Nf]
            H = zeros(2, Nf);

            switch lower(obj.source_type)
                case 'point'
                    R_src = norm(x_src - x_scat);
                    kR_src = k * R_src;

                    for n = 0:Nmax
                        % All vectorized over frequency [1 x Nf]
                        jn_ka  = getSphJ(n, ka);
                        hn_ka  = getSphH(n, 2, ka);
                        djn_ka = Diffractor.dsphj(n, ka);
                        dhn_ka = getDifSphH(n, 2, ka);
                        hn_src = getSphH(n, 2, kR_src);

                        % Scattering coefficient
                        An = -djn_ka ./ dhn_ka;

                        % Surface pressure for this order
                        % p_n = (2n+1)*(-i)^n * hn_src * (jn_ka + An*hn_ka) * Pn / kR_src
                        surf_n = (2*n+1) .* ((-1i).^n) .* hn_src .* ...
                                 (jn_ka + An .* hn_ka) ;

                        for ear = 1:2
                            H(ear, :) = H(ear, :) + surf_n .* Pn_ears(ear, n+1);
                        end
                    end

                    % Point source prefactor: -ik (from Green's function expansion)
                    H = H .* (-1i * k);

                case 'plane'
                    for n = 0:Nmax
                        jn_ka  = getSphJ(n, ka);
                        hn_ka  = getSphH(n, 2, ka);
                        djn_ka = Diffractor.dsphj(n, ka);
                        dhn_ka = getDifSphH(n, 2, ka);

                        An = -djn_ka ./ dhn_ka;

                        surf_n = (2*n+1) .* ((1i).^n) .* (jn_ka + An .* hn_ka);

                        for ear = 1:2
                            H(ear, :) = H(ear, :) + surf_n .* Pn_ears(ear, n+1);
                        end
                    end
            end

            % Handle DC (k=0): set to free-field value
            dc_idx = (k < 1e-15);
            if any(dc_idx)
                switch lower(obj.source_type)
                    case 'point'
                        H(:, dc_idx) = 1 / (4*pi*R_src);
                    case 'plane'
                        H(:, dc_idx) = 1;
                end
            end
        end

        function p0 = evaluate_freefield(obj, f, x_eval)
            % EVALUATE_FREEFIELD  Free-field pressure (no scatterer).
            %
            %   f      : [1 x Nf] — frequencies
            %   x_eval : [Npts x D] — positions (D = 2 or 3)
            %   p0     : [Npts x Nf] — complex pressure

            c = audioapp.util.Physics.speedOfSound();
            f = f(:).';
            Nf = numel(f);
            Npts = size(x_eval, 1);
            D = size(x_eval, 2);

            switch lower(obj.source_type)
                case 'point'
                    dx = bsxfun(@minus, x_eval, obj.source_pos(1:D));
                    R = sqrt(sum(dx.^2, 2));      % [Npts x 1]
                    k = 2*pi * f / c;             % [1 x Nf]
                    % 3D Green's function (monopole)
                    p0 = (1./(4*pi)) * bsxfun(@times, 1./max(R, 1e-15), ...
                        exp(-1i * bsxfun(@times, R, k)));

                case 'plane'
                    d = obj.source_dir(1:D);
                    d = d(:).' / max(norm(d), 1e-15);
                    k = 2*pi * f / c;
                    % Plane wave: p = exp(-i k . x)
                    proj = x_eval * d';           % [Npts x 1]
                    p0 = exp(-1i * bsxfun(@times, proj, k));
            end
        end
    end

    % ===== Private: rigid sphere series =====
    methods (Access = private)

        function p = rigid_sphere_point(obj, k, ka, kR, kR_src, cos_th, Nmax)
            % Point source incident on rigid sphere.
            % Returns total field at points with |kR| >= ka.
            %
            % Using the standard exterior Helmholtz solution:
            %   p(r,theta) = sum_n (2n+1) * i^n * [ j_n(kr) - (j_n'(ka)/h_n'(ka)) h_n(kr) ]
            %                * P_n(cos theta) * h_n(kr_s) / (kr_s)
            %
            % Simplified for 2D cross-section (azimuthal symmetry).

            Npts = numel(kR);
            kR = kR(:);            % force column
            cos_th = cos_th(:);    % force column
            p = zeros(Npts, 1);

            for n = 0:Nmax
                % Legendre polynomial P_n(cos_theta)
                Pn = getLegendre(n, 0, cos_th);
                Pn = Pn(:);                              % force column

                % Spherical Hankel (2nd kind) at source distance
                hn_src = getSphH(n, 2, kR_src);          % scalar

                % Derivative of sph Hankel at sphere surface
                dhn_a = getDifSphH(n, 2, ka);            % scalar

                % Derivative of sph Bessel at sphere surface
                djn_a = obj.dsphj(n, ka);                % scalar

                % Spherical Hankel at eval points
                hn_R = getSphH(n, 2, kR);
                hn_R = hn_R(:);                          % force column

                % Spherical Bessel at eval points
                jn_R = getSphJ(n, kR);
                jn_R = jn_R(:);                          % force column

                % Scattering coefficient
                An = -djn_a / dhn_a;

                % Total field contribution for order n
                coeff = (2*n+1) * (-1i)^n * hn_src;
                p = p + coeff * (jn_R + An * hn_R) .* Pn;
            end

            % Normalization: match free-field amplitude at large r
            p = p / (kR_src);
        end

        function p = rigid_sphere_plane(~, k, ka, kR, cos_th, Nmax)
            % Plane wave incident on rigid sphere.
            %
            %   p(r,theta) = sum_n (2n+1) * i^n * [ j_n(kr) - (j_n'(ka)/h_n'(ka)) h_n(kr) ]
            %                * P_n(cos theta)

            Npts = numel(kR);
            kR = kR(:);            % force column
            cos_th = cos_th(:);    % force column
            p = zeros(Npts, 1);

            for n = 0:Nmax
                Pn = getLegendre(n, 0, cos_th);
                Pn = Pn(:);

                dhn_a = getDifSphH(n, 2, ka);
                djn_a = Diffractor.dsphj(n, ka);

                hn_R = getSphH(n, 2, kR);
                hn_R = hn_R(:);
                jn_R = getSphJ(n, kR);
                jn_R = jn_R(:);

                An = -djn_a / dhn_a;

                p = p + (2*n+1) * (1i)^n * (jn_R + An * hn_R) .* Pn;
            end
        end

        function g = ff_greens(~, R, R_src, ~)
            % Free-field Green's function (DC limit placeholder)
            g = 1./(4*pi*max(R, 1e-15));
        end
    end

    methods (Static)
        function dj = dsphj(n, z)
            % Derivative of spherical Bessel function j_n(z) w.r.t. z.
            %   dj_n(z) = j_{n-1}(z) - (n+1)/z * j_n(z)     (n >= 1)
            %   dj_0(z) = -j_1(z)
            if n == 0
                dj = -getSphJ(1, z);
            else
                dj = getSphJ(n-1, z) - (n+1)./z .* getSphJ(n, z);
            end
        end
    end
end
