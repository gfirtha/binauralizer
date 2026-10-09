classdef soundbar_renderer < base_renderer
    %SOUNDBAR_RENDERER  Hybrid renderer of the 16-driver WFS soundbar with
    %optional discrete surrounds (standard 5.1 without the centre: L R Ls Rs).
    %
    %   Sources seen through the bar (|azimuth| <= window, from the reference
    %   listening position at the origin) are rendered by WFS on the bar.
    %   Outside the window the source is panned pairwise (2D VBAP, constant
    %   power) around the ring
    %       window edge (L) - L 30 - Ls 110 - Rs -110 - R -30 - window edge (R)
    %   where a window-edge node is a WFS source at the edge direction (same
    %   distance), so a source leaving the screen fades from the bar into L/R
    %   and on into the surrounds. Discrete feeds get the virtual source's
    %   distance attenuation and the extra delay (r_s - R_d)/c, matching the
    %   WFS level/arrival. Without discrete loudspeakers it is WFS on the bar.
    %
    %   The bar drivers are the loudspeakers with source_index <= Nbar (the
    %   layout creates them first); the rest are the discrete loudspeakers.

    properties
        setup
        Nbar = 16
        wfs                 % wfs_renderer on the bar only
        bar_ssd             % secondary_source_distribution of the bar
        proxy               % virtual source the WFS part renders
        barIdx              % main-SSD index of each bar_ssd loudspeaker
        discIdx             % main-SSD index of the discrete loudspeakers
        discAz              % their azimuths [deg, 0 = front, + = left]
        discR               % their distances [m]
        window = 0          % half-width of the WFS window [deg]
        g_wfs = 1           % gain of the WFS part
        g_disc              % gains of the discrete loudspeakers
        tau_disc            % extra delays of the discrete feeds [s]
        disc_bank = {}      % OLS convolvers (per-loudspeaker chain)
    end

    methods
        function obj = soundbar_renderer(vs, SSD, setup)
            obj = obj@base_renderer(vs, SSD);
            obj.setup = setup;
            if isfield(setup, 'SoundbarNbar'), obj.Nbar = setup.SoundbarNbar; end
            N = length(vs.source_signal.time_series);

            sidx = cellfun(@(l) l.source_index, SSD.loudspeakers);
            isBar = sidx <= obj.Nbar;
            obj.bar_ssd = secondary_source_distribution(SSD.loudspeakers(isBar), [0 0]);
            barS = cellfun(@(l) l.source_index, obj.bar_ssd.loudspeakers);
            obj.barIdx = arrayfun(@(s) find(sidx == s, 1), barS);
            obj.discIdx = find(~isBar);
            P = cell2mat(cellfun(@(l) l.position(1:2), SSD.loudspeakers(obj.discIdx)', 'UniformOutput', false));
            if isempty(P), P = zeros(0,2); end
            obj.discAz = atan2d(-P(:,1), P(:,2));
            obj.discR  = sqrt(sum(P.^2, 2));

            % WFS window: the bar seen from the origin, slightly inside its ends
            xb = obj.bar_ssd.x0;
            obj.window = min(25, 0.85*atan2d(max(abs(xb(:,1))), mean(xb(:,2))));

            obj.proxy = virtual_source(vs.source_index, vs.position, ...
                vs.orientation, vs.source_type, 'wfs');
            obj.proxy.source_signal.set_signal(zeros(N,1), vs.source_signal.fs);
            obj.wfs = wfs_renderer(obj.proxy, obj.bar_ssd, setup);
            obj.wfs.LazyFilterBank = true;

            for j = 1:numel(obj.discIdx)
                obj.disc_bank{j} = OLS_convolver(zeros(N,1), N);
            end
            obj.update_renderer();
        end

        function info = get_renderer_info(obj)
            if isempty(obj.discIdx)
                extra = 'No discrete loudspeakers: plain WFS on the bar.';
            else
                extra = sprintf(['Discrete L R Ls Rs: sources outside the +-%.0f deg WFS window ' ...
                    'are panned (constant power) from the window edge into L/R and the surrounds.'], obj.window);
            end
            info = sprintf(['Soundbar: %d-driver WFS bar (62.5 mm pitch).\n' ...
                'Sources seen through the bar are WFS point sources / plane waves.\n%s'], obj.Nbar, extra);
        end

        function obj = update_settings(obj, setup)
            obj.setup = setup;
            obj.wfs.setup = setup;
            obj.update_renderer('settings');
        end

        function obj = update_renderer(obj, type)
            if nargin < 2, type = 'init'; end
            if strcmp(type, 'receiver_rotated'), return; end
            if startsWith(string(type), 'loudspeaker')
                obj.bar_ssd.update_geometry();
            end
            obj.compute_panning();
            obj.wfs.update_renderer(type);
            obj.build_discrete_filters();
            obj.filter_version = obj.filter_version + 1;
        end

        function render(obj)
            x = obj.virtual_source.source_signal.time_series;
            obj.proxy.source_signal.set_signal(obj.g_wfs * x);
            if obj.g_wfs > 0
                obj.wfs.render();                   % adds to the bar loudspeakers
            end
            for j = 1:numel(obj.discIdx)
                if obj.g_disc(j) == 0, continue; end
                n = obj.discIdx(j);
                obj.output_signal{n}.set_signal(obj.disc_bank{j}.convolve(x));
                obj.ssd.loudspeakers{n}.source_signal.add_signals(obj.output_signal{n});
            end
        end

        function d = get_driving_filters(obj)
            % [N x Nls] in main SSD order (fast binaural path)
            dw = obj.wfs.get_driving_filters();
            N = size(dw, 1);
            d = zeros(N, obj.ssd.N);
            d(:, obj.barIdx) = obj.g_wfs * dw;
            for j = 1:numel(obj.discIdx)
                d(:, obj.discIdx(j)) = obj.discrete_filter(j, N);
            end
        end

        % --- overlays: the WFS part ---
        function store = drawOverlay(obj, ax, layers, settings, store)
            store = obj.wfs.drawOverlay(ax, layers, settings, store);
        end

        function settings = defaultOverlaySettings(obj)
            settings = obj.wfs.defaultOverlaySettings();
        end
    end

    methods (Access = private)
        function compute_panning(obj)
            vs = obj.virtual_source;
            p = vs.position(1:2);
            az = atan2d(-p(1), p(2));
            nD = numel(obj.discIdx);
            obj.g_disc = zeros(nD, 1);
            obj.tau_disc = zeros(nD, 1);
            obj.g_wfs = 1;
            obj.proxy.position = vs.position;
            if nD == 0 || abs(az) <= obj.window
                return;
            end
            % ring nodes, CCW angles in [0, 360): window edges + discrete
            w = obj.window;
            ang = [w; mod(obj.discAz, 360); 360 - w];
            who = [0; (1:nD)'; -1];                 % 0/-1: left/right window edge
            [ang, o] = sort(ang); who = who(o);
            a = mod(az, 360);
            i = find(ang <= a, 1, 'last');
            if isempty(i) || i == numel(ang), i = numel(ang) - 1; end
            u = @(t) [cosd(t); sind(t)];
            g = [u(ang(i)), u(ang(i+1))] \ u(a);
            g = max(g, 0);
            g = g / max(norm(g), eps);              % constant power
            rs = norm(p);
            plane = strcmp(vs.source_type.Shape, 'plane_wave');
            obj.g_wfs = 0;
            for k = 1:2
                n = who(i+k-1);
                if n <= 0
                    obj.g_wfs = g(k);
                    edge = w * (1 - 2*(n < 0));     % +w left edge, -w right edge
                    obj.proxy.position(1:2) = rs * [-sind(edge), cosd(edge)];
                else
                    if plane
                        obj.g_disc(n) = g(k) * obj.discR(n);
                    else
                        obj.g_disc(n) = g(k) * obj.discR(n) / max(rs, 0.1);
                        obj.tau_disc(n) = max(0, (rs - obj.discR(n)) / audioapp.util.Physics.speedOfSound());
                    end
                end
            end
        end

        function h = discrete_filter(obj, j, N)
            % gain x delay, same time origin as the WFS FIRs (fftshift: N/2+1)
            h = zeros(N, 1);
            if obj.g_disc(j) == 0, return; end
            fs = obj.virtual_source.source_signal.fs;
            k = min(N, N/2 + 1 + round(obj.tau_disc(j) * fs));
            h(k) = obj.g_disc(j);
        end

        function build_discrete_filters(obj)
            N = length(obj.virtual_source.source_signal.time_series);
            for j = 1:numel(obj.discIdx)
                obj.disc_bank{j}.update_coefficients(obj.discrete_filter(j, N));
            end
        end
    end
end
