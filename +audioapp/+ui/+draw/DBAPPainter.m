classdef DBAPPainter
    % DBAPPainter  Draws DBAP renderer overlays:
    %   1) Gain-proportional dots on each loudspeaker (color intensity ∝ gain)
    %   2) Lines from VS to each loudspeaker (width ∝ gain)
    %   3) Gain labels on each active speaker (dB)
    %
    % rp fields (from dbap_renderer.get_renderer_visuals):
    %   .x0     : N×2 loudspeaker positions
    %   .G_vec  : N×1 gain vector (linear)
    %   .vs     : 1×2 virtual source position
    %
    % settings:
    %   .showGainDots   : gain-colored dots
    %   .showGainLines  : lines from VS to speakers
    %   .showGainLabels : dB labels

    methods (Static)
        function store = draw(ax, rp, settings, store)
            if nargin < 4 || isempty(store), store = struct(); end
            hold(ax, 'on');

            G = abs(rp.G_vec);
            Gmax = max(G) + eps;
            Gnorm = G / Gmax;        % 0..1
            activeIdx = find(G > 0.01 * Gmax);  % show only > 1% of max

            % ── 1) Gain-proportional dots ──
            if settings.showGainDots
                xA = rp.x0(:, 1);
                yA = rp.x0(:, 2);
                dotSizes = 30 + 200 * Gnorm;
                % Color: green → yellow → red with increasing gain
                cmap = [1-Gnorm, 0.3+0.5*Gnorm, 0.1*ones(size(Gnorm))];
                cmap = min(max(cmap, 0), 1);

                if ~isfield(store, 'gainDots') || ~isgraphics(store.gainDots)
                    store.gainDots = scatter(ax, xA, yA, dotSizes, cmap, ...
                        'filled', 'MarkerEdgeColor', 'none', ...
                        'MarkerFaceAlpha', 0.8, ...
                        'HitTest', 'off', 'Tag', 'DBAPgainDots');
                else
                    set(store.gainDots, 'XData', xA, 'YData', yA, ...
                        'SizeData', dotSizes, 'CData', cmap, 'Visible', 'on');
                end
            else
                if isfield(store, 'gainDots') && isgraphics(store.gainDots)
                    set(store.gainDots, 'Visible', 'off');
                end
            end

            % ── 2) Lines from VS to active speakers ──
            if settings.showGainLines && ~isempty(activeIdx)
                xLines = [];
                yLines = [];
                for k = 1:numel(activeIdx)
                    xLines = [xLines, rp.vs(1), rp.x0(activeIdx(k),1), NaN]; %#ok
                    yLines = [yLines, rp.vs(2), rp.x0(activeIdx(k),2), NaN]; %#ok
                end

                if ~isfield(store, 'gainLines') || ~isgraphics(store.gainLines)
                    store.gainLines = line('Parent', ax, ...
                        'XData', xLines, 'YData', yLines, ...
                        'Color', [0.4 0.7 0.3 0.4], 'LineWidth', 1.5, ...
                        'LineStyle', '-', 'HitTest', 'off', 'Tag', 'DBAPgainLines');
                else
                    set(store.gainLines, 'XData', xLines, 'YData', yLines, 'Visible', 'on');
                end
            else
                if isfield(store, 'gainLines') && isgraphics(store.gainLines)
                    set(store.gainLines, 'Visible', 'off');
                end
            end

            % ── 3) Gain labels (dB) ──
            if settings.showGainLabels && ~isempty(activeIdx)
                if ~isfield(store, 'gainTexts') || ~iscell(store.gainTexts)
                    store.gainTexts = {};
                end

                for k = 1:numel(activeIdx)
                    idx = activeIdx(k);
                    g_dB = 20*log10(G(idx) + 1e-12);
                    lbl = sprintf('%.1f dB', g_dB);

                    tx = rp.x0(idx, 1);
                    ty = rp.x0(idx, 2) + 0.15;

                    if k > numel(store.gainTexts) || ~isgraphics(store.gainTexts{k})
                        store.gainTexts{k} = text(ax, tx, ty, lbl, ...
                            'HorizontalAlignment', 'center', 'FontSize', 8, ...
                            'FontWeight', 'bold', 'Color', [0.2 0.5 0.15], ...
                            'HitTest', 'off', 'Tag', sprintf('DBAPgain%d', k));
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
