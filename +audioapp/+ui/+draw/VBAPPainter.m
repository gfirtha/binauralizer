classdef VBAPPainter
    % VBAPPainter  Draws VBAP renderer overlays:
    %   1) Active loudspeaker dots (size ∝ gain)
    %   2) Lines from origin to active speakers (opacity ∝ gain)
    %   3) Spread arc showing the angular extent
    %   4) Gain ring: polar gain pattern around the SSD center
    %
    % rp fields (from vbap_renderer.get_renderer_visuals):
    %   .x0      : N×2 loudspeaker positions
    %   .G_vec   : N×1 gain vector
    %   .vs      : 1×2 virtual source position
    %   .spread  : spread angle (degrees)
    %   .pNorm   : normalization exponent
    %   .ls_dirs : N×1 loudspeaker angles (rad)
    %
    % settings fields:
    %   .showActiveLS   : gain-proportional dots on active speakers
    %   .showPairLink   : lines from center to active speakers
    %   .showSpreadArc  : arc showing spread extent
    %   .showGainRing   : polar gain pattern

    methods (Static)

        function store = draw(ax, rp, settings, store)
            if nargin < 4 || isempty(store), store = struct(); end
            hold(ax, 'on');

            activeIdx = find(abs(rp.G_vec) > 1e-12);
            gains     = abs(rp.G_vec(activeIdx));

            % Geometry: SSD center (use mean of loudspeaker positions)
            center = mean(rp.x0, 1);
            R_ssd  = mean(sqrt(sum(bsxfun(@minus, rp.x0, center).^2, 2)));

            % VS angle
            vs_vec = rp.vs - center;
            vs_angle = atan2(vs_vec(2), vs_vec(1));

            % ── 1) Active loudspeaker dots ──
            if settings.showActiveLS && ~isempty(activeIdx)
                xA = rp.x0(activeIdx, 1);
                yA = rp.x0(activeIdx, 2);
                dotSizes = 40 + 200 * (gains / max(gains(:) + eps));

                if ~isfield(store, 'activeDots') || ~isgraphics(store.activeDots)
                    store.activeDots = scatter(ax, xA, yA, dotSizes, ...
                        'filled', 'MarkerFaceColor', [0.9 0.1 0.1], ...
                        'MarkerFaceAlpha', 0.7, 'MarkerEdgeColor', [0.6 0 0], ...
                        'LineWidth', 1.5, 'HitTest', 'off', 'Tag', 'VBAPactiveDots');
                else
                    set(store.activeDots, 'XData', xA, 'YData', yA, ...
                        'SizeData', dotSizes, 'Visible', 'on');
                end
            else
                if isfield(store, 'activeDots') && isgraphics(store.activeDots)
                    set(store.activeDots, 'Visible', 'off');
                end
            end

            % ── 2) Lines from center to active speakers ──
            if settings.showPairLink && ~isempty(activeIdx)
                xLines = [];
                yLines = [];
                for k = 1:numel(activeIdx)
                    xLines = [xLines, center(1), rp.x0(activeIdx(k),1), NaN]; %#ok
                    yLines = [yLines, center(2), rp.x0(activeIdx(k),2), NaN]; %#ok
                end
                if ~isfield(store, 'pairLines') || ~isgraphics(store.pairLines)
                    store.pairLines = line('Parent', ax, ...
                        'XData', xLines, 'YData', yLines, ...
                        'Color', [0.9 0.1 0.1 0.4], 'LineWidth', 1.5, ...
                        'LineStyle', '-', 'HitTest', 'off', 'Tag', 'VBAPpairLines');
                else
                    set(store.pairLines, 'XData', xLines, 'YData', yLines, 'Visible', 'on');
                end
            else
                if isfield(store, 'pairLines') && isgraphics(store.pairLines)
                    set(store.pairLines, 'Visible', 'off');
                end
            end

            % ── 3) Spread arc ──
            if settings.showSpreadArc && rp.spread > 0
                half_spread = rp.spread / 2 * pi / 180;
                arc_angles = linspace(vs_angle - half_spread, vs_angle + half_spread, 60);
                arc_r = R_ssd * 1.08;  % slightly outside the SSD ring
                xArc = center(1) + arc_r * cos(arc_angles);
                yArc = center(2) + arc_r * sin(arc_angles);

                % Draw arc
                if ~isfield(store, 'spreadArc') || ~isgraphics(store.spreadArc)
                    store.spreadArc = line('Parent', ax, ...
                        'XData', xArc, 'YData', yArc, ...
                        'Color', [0.2 0.6 1.0 0.8], 'LineWidth', 3, ...
                        'LineStyle', '-', 'HitTest', 'off', 'Tag', 'VBAPspreadArc');
                else
                    set(store.spreadArc, 'XData', xArc, 'YData', yArc, 'Visible', 'on');
                end

                % Draw spread boundary lines (two dashed lines from center)
                xBnd = [center(1), center(1) + arc_r * cos(vs_angle - half_spread), NaN, ...
                        center(1), center(1) + arc_r * cos(vs_angle + half_spread)];
                yBnd = [center(2), center(2) + arc_r * sin(vs_angle - half_spread), NaN, ...
                        center(2), center(2) + arc_r * sin(vs_angle + half_spread)];

                if ~isfield(store, 'spreadBounds') || ~isgraphics(store.spreadBounds)
                    store.spreadBounds = line('Parent', ax, ...
                        'XData', xBnd, 'YData', yBnd, ...
                        'Color', [0.2 0.6 1.0 0.5], 'LineWidth', 1, ...
                        'LineStyle', '--', 'HitTest', 'off', 'Tag', 'VBAPspreadBounds');
                else
                    set(store.spreadBounds, 'XData', xBnd, 'YData', yBnd, 'Visible', 'on');
                end

                % Spread angle label
                lbl_angle = vs_angle;
                lx = center(1) + (arc_r + 0.12) * cos(lbl_angle);
                ly = center(2) + (arc_r + 0.12) * sin(lbl_angle);
                lbl_str = sprintf('%.0f°', rp.spread);

                if ~isfield(store, 'spreadLabel') || ~isgraphics(store.spreadLabel)
                    store.spreadLabel = text(ax, lx, ly, lbl_str, ...
                        'HorizontalAlignment', 'center', 'FontSize', 10, ...
                        'FontWeight', 'bold', 'Color', [0.1 0.4 0.8], ...
                        'HitTest', 'off', 'Tag', 'VBAPspreadLabel');
                else
                    set(store.spreadLabel, 'Position', [lx ly 0], ...
                        'String', lbl_str, 'Visible', 'on');
                end
            else
                % Hide spread visuals
                flds = {'spreadArc', 'spreadBounds', 'spreadLabel'};
                for k = 1:numel(flds)
                    if isfield(store, flds{k}) && isgraphics(store.(flds{k}))
                        set(store.(flds{k}), 'Visible', 'off');
                    end
                end
            end

            % ── 4) Gain ring: polar gain pattern ──
            if settings.showGainRing && any(abs(rp.G_vec) > 1e-12)
                Nls = numel(rp.G_vec);
                ls_ang = rp.ls_dirs;

                % Sort by angle [0, 2*pi]
                ls_ang_wrapped = mod(ls_ang, 2*pi);
                [ls_ang_sorted, sI] = sort(ls_ang_wrapped);
                G_sorted = abs(rp.G_vec(sI));

                % Close the loop by wrapping endpoints
                ang_interp = [ls_ang_sorted(end) - 2*pi; ls_ang_sorted; ls_ang_sorted(1) + 2*pi];
                G_interp   = [G_sorted(end);              G_sorted;      G_sorted(1)];

                % Fine grid
                th_fine = linspace(0, 2*pi, 360)';
                G_fine = interp1(ang_interp, G_interp, th_fine, 'linear', 0);
                G_fine = max(G_fine, 0);

                % Map gain to radius
                ring_scale = R_ssd * 0.3;
                R_ring = R_ssd + ring_scale * G_fine / (max(G_fine) + eps);

                xRing = center(1) + R_ring .* cos(th_fine);
                yRing = center(2) + R_ring .* sin(th_fine);

                if ~isfield(store, 'gainRing') || ~isgraphics(store.gainRing)
                    store.gainRing = patch(ax, ...
                        'XData', [xRing; xRing(1)], 'YData', [yRing; yRing(1)], ...
                        'FaceColor', [1.0 0.3 0.3], 'FaceAlpha', 0.15, ...
                        'EdgeColor', [0.8 0.1 0.1], 'EdgeAlpha', 0.6, ...
                        'LineWidth', 1.5, 'HitTest', 'off', 'Tag', 'VBAPgainRing');
                else
                    set(store.gainRing, ...
                        'XData', [xRing; xRing(1)], 'YData', [yRing; yRing(1)], ...
                        'Visible', 'on');
                end
            else
                if isfield(store, 'gainRing') && isgraphics(store.gainRing)
                    set(store.gainRing, 'Visible', 'off');
                end
            end

            % ── 5) Gain text labels on active speakers ──
            if settings.showActiveLS && ~isempty(activeIdx)
                if ~isfield(store, 'gainTexts') || ~iscell(store.gainTexts)
                    store.gainTexts = {};
                end
                for k = 1:numel(activeIdx)
                    lbl = sprintf('%.2f', gains(k));
                    tx = rp.x0(activeIdx(k), 1);
                    ty = rp.x0(activeIdx(k), 2) + 0.12;

                    if k > numel(store.gainTexts) || ~isgraphics(store.gainTexts{k})
                        store.gainTexts{k} = text(ax, tx, ty, lbl, ...
                            'HorizontalAlignment', 'center', 'FontSize', 8, ...
                            'FontWeight', 'bold', 'Color', [0.7 0 0], ...
                            'HitTest', 'off', 'Tag', sprintf('VBAPgain%d', k));
                    else
                        set(store.gainTexts{k}, ...
                            'Position', [tx ty 0], 'String', lbl, 'Visible', 'on');
                    end
                end
                for k = numel(activeIdx)+1 : numel(store.gainTexts)
                    if isgraphics(store.gainTexts{k})
                        set(store.gainTexts{k}, 'Visible', 'off');
                    end
                end
            else
                if isfield(store, 'gainTexts') && iscell(store.gainTexts)
                    for k = 1:numel(store.gainTexts)
                        if isgraphics(store.gainTexts{k})
                            set(store.gainTexts{k}, 'Visible', 'off');
                        end
                    end
                end
            end
        end
    end
end
