classdef TDPainter
    % TDPainter  Draws time-delay stereo renderer overlays:
    %   1) Highlight the active stereo pair
    %   2) Line connecting the pair through the ref point
    %   3) Delay labels above each active speaker (in ms)
    %
    % rp fields (from time_delay_renderer.get_renderer_visuals):
    %   .x0          : N×2 loudspeaker positions
    %   .vs          : 1×2 virtual source position
    %   .stereo_pair : [2×1] indices of the stereo pair
    %   .fi_base     : half-base angle (rad)
    %   .refPoint    : [1×2] ref point
    %   .method      : string method name
    %   .delays_ms   : [N×1] delay per speaker (ms)
    %
    % settings fields:
    %   .showDelayLabels : show delay text on active speakers
    %   .showPairLink    : show connecting line between pair

    methods (Static)

        function store = draw(ax, rp, settings, store)
            if nargin < 4 || isempty(store), store = struct(); end
            hold(ax, 'on');

            ixs = rp.stereo_pair;
            hasActivePair = ~isempty(ixs) && numel(ixs) >= 2;

            % ── 1) Active speaker dots ──
            if hasActivePair
                xA = rp.x0(ixs, 1);
                yA = rp.x0(ixs, 2);
                dotSizes = [160; 160];

                if ~isfield(store, 'activeDots') || ~isgraphics(store.activeDots)
                    store.activeDots = scatter(ax, xA, yA, dotSizes, ...
                        'filled', 'MarkerFaceColor', [0.15 0.55 0.85], ...
                        'MarkerFaceAlpha', 0.6, 'MarkerEdgeColor', [0.05 0.3 0.6], ...
                        'LineWidth', 1.5, 'HitTest', 'off', 'Tag', 'TDactiveDots');
                else
                    set(store.activeDots, 'XData', xA, 'YData', yA, ...
                        'SizeData', dotSizes, 'Visible', 'on');
                end
            else
                if isfield(store, 'activeDots') && isgraphics(store.activeDots)
                    set(store.activeDots, 'Visible', 'off');
                end
            end

            % ── 2) Line connecting pair through ref point ──
            if settings.showPairLink && hasActivePair
                xL = [rp.x0(ixs(1),1), rp.refPoint(1), rp.x0(ixs(2),1)];
                yL = [rp.x0(ixs(1),2), rp.refPoint(2), rp.x0(ixs(2),2)];

                if ~isfield(store, 'pairLine') || ~isgraphics(store.pairLine)
                    store.pairLine = line('Parent', ax, ...
                        'XData', xL, 'YData', yL, ...
                        'Color', [0.15 0.55 0.85 0.5], 'LineWidth', 2, ...
                        'LineStyle', '--', 'HitTest', 'off', 'Tag', 'TDpairLine');
                else
                    set(store.pairLine, 'XData', xL, 'YData', yL, 'Visible', 'on');
                end
            else
                if isfield(store, 'pairLine') && isgraphics(store.pairLine)
                    set(store.pairLine, 'Visible', 'off');
                end
            end

            % ── 3) Parameter labels above active speakers ──
            if settings.showDelayLabels && hasActivePair
                if ~isfield(store, 'delayTexts') || ~iscell(store.delayTexts)
                    store.delayTexts = {};
                end

                mtype = 'delay';
                if isfield(rp, 'method_type'), mtype = rp.method_type; end

                for k = 1:2
                    idx = ixs(k);
                    dt_val = rp.delays_ms(idx);
                    g_val  = rp.gains_lin(idx);
                    g_dB   = 20*log10(abs(g_val) + 1e-12);

                    switch mtype
                        case 'delay'
                            if dt_val >= 0
                                lbl = sprintf('+%.2f ms', dt_val);
                            else
                                lbl = sprintf('%.2f ms', dt_val);
                            end
                        case 'amplitude'
                            lbl = sprintf('%.1f dB', g_dB);
                        case 'hybrid'
                            if dt_val >= 0
                                lbl = sprintf('+%.2f ms | %.1f dB', dt_val, g_dB);
                            else
                                lbl = sprintf('%.2f ms | %.1f dB', dt_val, g_dB);
                            end
                        otherwise
                            lbl = sprintf('%.2f ms', dt_val);
                    end

                    tx = rp.x0(idx, 1);
                    ty = rp.x0(idx, 2) + 0.18;

                    if k > numel(store.delayTexts) || ~isgraphics(store.delayTexts{k})
                        store.delayTexts{k} = text(ax, tx, ty, lbl, ...
                            'HorizontalAlignment', 'center', 'FontSize', 10, ...
                            'FontWeight', 'bold', 'Color', [0.05 0.3 0.6], ...
                            'BackgroundColor', [1 1 1 0.7], 'Margin', 2, ...
                            'HitTest', 'off', 'Tag', sprintf('TDdelay%d', k));
                    else
                        set(store.delayTexts{k}, ...
                            'Position', [tx ty 0], 'String', lbl, 'Visible', 'on');
                    end
                end

                % Hide extra text handles
                for k = 3:numel(store.delayTexts)
                    if isgraphics(store.delayTexts{k})
                        set(store.delayTexts{k}, 'Visible', 'off');
                    end
                end
            else
                if isfield(store, 'delayTexts') && iscell(store.delayTexts)
                    for k = 1:numel(store.delayTexts)
                        if isgraphics(store.delayTexts{k})
                            set(store.delayTexts{k}, 'Visible', 'off');
                        end
                    end
                end
            end
        end
    end
end
