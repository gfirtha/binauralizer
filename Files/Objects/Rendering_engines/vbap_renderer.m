classdef vbap_renderer < base_renderer
    % VBAP_RENDERER  Vector Base Amplitude Panning with spread and p-norm.
    %
    % Properties set from the panel / setup:
    %   spread   : source spread angle in degrees (0 = point source)
    %   pNorm    : gain normalization exponent (1 = amplitude, 2 = power)
    %   nSpread  : number of auxiliary spread sources (default 8)

    properties
        G_vec           % [N×1] gain vector (per loudspeaker)
        spread  = 0     % spread angle (degrees)
        pNorm   = 2     % normalization exponent (1..2)
        nSpread = 8     % number of auxiliary spread sources
        ls_pairs        % [2×Npairs] loudspeaker pair indices (2D)
        ls_tris         % [Ntri×3] loudspeaker triplet indices (3D)
        ls_dirs_rad     % [N×1] loudspeaker angles (radians)
        is3D    = false % true if layout has height variation
    end

    methods
        function obj = vbap_renderer(virtual_source, SSD)
            obj = obj@base_renderer(virtual_source, SSD);
            obj.precompute_layout();
            obj.update_renderer();
        end

        function info = get_renderer_info(obj)
            info = sprintf(['VBAP — Vector Base Amplitude Panning\n' ...
                'Spread: %.0f°   p-norm: %.1f   Spread sources: %d'], ...
                obj.spread, obj.pNorm, obj.nSpread);
        end

        function obj = update_settings(obj, renderer_setup)
            if isfield(renderer_setup, 'VBAPspread')
                obj.spread = renderer_setup.VBAPspread;
            end
            if isfield(renderer_setup, 'VBAPpNorm')
                obj.pNorm = renderer_setup.VBAPpNorm;
            end
            if isfield(renderer_setup, 'VBAPnSpread')
                obj.nSpread = max(4, round(renderer_setup.VBAPnSpread));
            end
            obj.update_renderer();
        end

        function obj = update_renderer(obj, ~)
            Nls = obj.ssd.N;
            obj.G_vec = zeros(Nls, 1);

            if Nls < 2, return; end

            % Virtual source direction (angle)
            p = obj.virtual_source.position;
            vs_angle = atan2(p(2), p(1)) * 180/pi;
            vs_elev = 0;                 % used by 3D layouts only
            if numel(p) >= 3
                vs_elev = atan2(p(3), norm(p(1:2))) * 180/pi;
            end

            if obj.spread <= 0
                % ── No spread: standard VBAP (single pair) ──
                obj.G_vec = obj.vbap_single(vs_angle, vs_elev);
            else
                % ── Spread: average gains over auxiliary directions ──
                spread_dirs = getSpreadSrcDirs(vs_angle, obj.spread, obj.nSpread);
                % spread_dirs: [(nSpread+1) × 2] unit vectors
                G_acc = zeros(Nls, 1);
                for k = 1:size(spread_dirs, 1)
                    dir_deg = atan2(spread_dirs(k,2), spread_dirs(k,1)) * 180/pi;
                    G_acc = G_acc + obj.vbap_single(dir_deg, vs_elev);
                end
                G_acc = G_acc / size(spread_dirs, 1);

                % Normalize with p-norm
                pn = max(abs(G_acc)) ^ (1 - 1/obj.pNorm);
                if pn > 0
                    G_acc = G_acc / (norm(G_acc, obj.pNorm) + eps) * pn;
                end
                obj.G_vec = G_acc;
            end
        end

        function render(obj)
            for n = 1:obj.ssd.N
                obj.output_signal{n}.set_signal( ...
                    obj.G_vec(n) * obj.virtual_source.source_signal.time_series);
            end
            obj.add_output_to_ssd_signal;
        end

        % ─── Visuals ───
        function output = get_renderer_visuals(obj)
            output.x0      = obj.ssd.x0;
            output.G_vec   = obj.G_vec;
            output.vs      = obj.virtual_source.position;
            output.spread  = obj.spread;
            output.pNorm   = obj.pNorm;
            output.ls_dirs = obj.ls_dirs_rad;
        end

        function store = drawOverlay(obj, ax, layers, settings, store)
            rp = obj.get_renderer_visuals();
            store = audioapp.ui.draw.VBAPPainter.draw(ax, rp, settings, store);
        end

        function settings = defaultOverlaySettings(~)
            settings = struct('showActiveLS', true, 'showPairLink', true, ...
                              'showSpreadArc', true, 'showGainRing', true);
        end
    end

    methods (Access = private)

        function precompute_layout(obj)
            % Precompute loudspeaker angles and adjacent pairs/triplets.
            Nls = obj.ssd.N;
            x0 = obj.ssd.x0;

            % Check if layout is 3D (speakers at different heights)
            if size(x0, 2) >= 3 && (max(x0(:,3)) - min(x0(:,3))) > 0.01
                obj.is3D = true;
            else
                obj.is3D = false;
            end

            if obj.is3D
                % 3D VBAP: Delaunay triangulation on unit sphere
                x0_unit = bsxfun(@rdivide, x0(:,1:3), sqrt(sum(x0(:,1:3).^2, 2)) + eps);
                try
                    obj.ls_tris = convhulln(x0_unit);  % [Ntri × 3]
                catch
                    % Fallback: use delaunay on azimuth/elevation
                    [az, el, ~] = cart2sph(x0_unit(:,1), x0_unit(:,2), x0_unit(:,3));
                    obj.ls_tris = delaunay(az, el);
                end
                obj.ls_pairs = [];
            else
                % 2D VBAP: adjacent pairs by azimuth angle
                obj.ls_dirs_rad = atan2(x0(:,2), x0(:,1));
                [~, sortIdx] = sort(obj.ls_dirs_rad);
                obj.ls_pairs = [sortIdx, circshift(sortIdx, -1)]';  % [2 × Nls]
                obj.ls_tris = [];
            end
        end

        function g = vbap_single(obj, vs_angle_deg, vs_elev_deg)
            % Compute VBAP gains for a single direction.
            % Returns [N×1] gain vector (non-negative, p-norm normalized).

            Nls = obj.ssd.N;
            g = zeros(Nls, 1);
            if Nls < 2, return; end

            x0 = obj.ssd.x0;

            if obj.is3D
                % ── 3D VBAP: triplet search ──
                % azimuth + elevation -> 3D unit vector
                if nargin < 3, vs_elev_deg = 0; end
                vs_rad = vs_angle_deg * pi / 180;
                el_rad = vs_elev_deg * pi / 180;
                vs_uv = [cos(el_rad)*cos(vs_rad), cos(el_rad)*sin(vs_rad), sin(el_rad)]';

                % Normalize loudspeaker positions to unit vectors
                R_ls = sqrt(sum(x0(:,1:3).^2, 2));
                L = bsxfun(@rdivide, x0(:,1:3), R_ls + eps)';  % [3 × Nls]

                best_g = [];
                best_tri = [];
                best_cond = inf;

                for t = 1:size(obj.ls_tris, 1)
                    idx = obj.ls_tris(t, :)';
                    Lt = L(:, idx);        % [3 × 3]
                    gt = Lt \ vs_uv;       % [3 × 1]

                    if all(gt >= -1e-10)
                        cn = cond(Lt);
                        if cn < best_cond
                            best_cond = cn;
                            best_g = max(gt, 0);
                            best_tri = idx;
                        end
                    end
                end

                if ~isempty(best_tri)
                    gn = best_g / (norm(best_g, obj.pNorm) + eps);
                    g(best_tri) = gn;
                end
            else
                % ── 2D VBAP: pair search ──
                R_ls = sqrt(sum(x0(:,1:2).^2, 2));
                L = bsxfun(@rdivide, x0(:,1:2), R_ls + eps)';  % [2 × Nls]

                vs_rad = vs_angle_deg * pi / 180;
                vs_uv = [cos(vs_rad); sin(vs_rad)];

                best_g = [];
                best_pair = [];
                best_dl = inf;

                for p = 1:size(obj.ls_pairs, 2)
                    idx = obj.ls_pairs(:, p);
                    Lp = L(:, idx);          % [2 × 2]
                    gp = Lp \ vs_uv;         % [2 × 1]

                    if all(gp >= -1e-10)
                        dl = norm(x0(idx(1),:) - x0(idx(2),:));
                        if dl < best_dl
                            best_dl = dl;
                            best_g = max(gp, 0);
                            best_pair = idx;
                        end
                    end
                end

                if ~isempty(best_pair)
                    gn = best_g / (norm(best_g, obj.pNorm) + eps);
                    g(best_pair) = gn;
                end
            end
        end
    end
end
