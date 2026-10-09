classdef BeamformerPainter
    % BeamformerPainter  Draws beamformer overlays for both DAS and Focus modes.
    %
    % DAS mode:
    %   - Arrow from array center in the steering direction
    %   - Parallel wavefront lines across the array
    %
    % Focus mode:
    %   - Focal point marker at VS position
    %   - Converging ray lines from each element to the focal point
    %   - Element dots sized by amplitude weight
    %
    % Common:
    %   - Element dots (blue, sized by weight in Focus mode)
    %   - Steering angle label

    methods (Static)

        function store = draw(ax, rp, settings, store)
            if nargin < 4 || isempty(store), store = struct(); end
            hold(ax, 'on');

            % ── 1) Element highlight dots ──
            if settings.showElementHighlight
                xE = rp.x0(:,1);
                yE = rp.x0(:,2);

                % In focus mode, size dots by weight
                if strcmpi(rp.mode, 'focus') && isfield(rp, 'weights')
                    w = abs(rp.weights);
                    dotSize = 30 + 120 * (w / max(w + eps));
                else
                    dotSize = 60 * ones(size(xE));
                end

                if ~isfield(store, 'elemDots') || ~isgraphics(store.elemDots)
                    store.elemDots = scatter(ax, xE, yE, dotSize, ...
                        'filled', 'MarkerFaceColor', [0.2 0.6 1.0], ...
                        'MarkerFaceAlpha', 0.8, 'MarkerEdgeColor', [0.1 0.3 0.7], ...
                        'LineWidth', 1.0, 'HitTest', 'off', 'Tag', 'BFelemDots');
                else
                    set(store.elemDots, 'XData', xE, 'YData', yE, ...
                        'SizeData', dotSize, 'Visible', 'on');
                end
            else
                if isfield(store, 'elemDots') && isgraphics(store.elemDots)
                    set(store.elemDots, 'Visible', 'off');
                end
            end

            % ── 2) Mode-specific visuals ──
            if settings.showBeamDirection
                switch lower(rp.mode)
                    case 'das'
                        store = audioapp.ui.draw.BeamformerPainter.drawDAS(ax, rp, store);
                        store = audioapp.ui.draw.BeamformerPainter.hideFocusVisuals(store);
                    case 'focus'
                        store = audioapp.ui.draw.BeamformerPainter.drawFocus(ax, rp, store);
                        store = audioapp.ui.draw.BeamformerPainter.hideDASVisuals(store);
                    otherwise
                        store = audioapp.ui.draw.BeamformerPainter.hideDASVisuals(store);
                        store = audioapp.ui.draw.BeamformerPainter.hideFocusVisuals(store);
                end

                % Angle label (common)
                lbl = sprintf('%.1f deg  [%s]', rp.steerAngle, upper(rp.mode));
                cx = rp.center(1);
                cy = rp.center(2);
                lblX = cx + 0.05;
                lblY = cy - 0.15;

                if ~isfield(store, 'angleLabel') || ~isgraphics(store.angleLabel)
                    store.angleLabel = text(ax, lblX, lblY, lbl, ...
                        'FontSize', 9, 'FontWeight', 'bold', ...
                        'Color', [0.1 0.3 0.7], 'HitTest', 'off', ...
                        'HorizontalAlignment', 'left', 'Tag', 'BFangleLabel');
                else
                    set(store.angleLabel, 'Position', [lblX lblY 0], ...
                        'String', lbl, 'Visible', 'on');
                end
            else
                store = audioapp.ui.draw.BeamformerPainter.hideDASVisuals(store);
                store = audioapp.ui.draw.BeamformerPainter.hideFocusVisuals(store);
                if isfield(store, 'angleLabel') && isgraphics(store.angleLabel)
                    set(store.angleLabel, 'Visible', 'off');
                end
            end
        end

        % ────────────────────────────────────────────────
        %  DAS visuals: steering arrow
        % ────────────────────────────────────────────────
        function store = drawDAS(ax, rp, store)
            cx = rp.center(1);
            cy = rp.center(2);
            dir = rp.vs(1:2) - rp.center(1:2);
            if norm(dir) > 1e-12
                dir = dir / norm(dir);
            else
                dir = [0, 1];
            end

            aperture = max(sqrt(sum((rp.x0(:,1:2) - rp.center(1:2)).^2, 2)));
            arrowLen = max(aperture * 1.5, 0.3);

            xTip = cx + dir(1) * arrowLen;
            yTip = cy + dir(2) * arrowLen;

            % Main beam line
            if ~isfield(store, 'beamLine') || ~isgraphics(store.beamLine)
                store.beamLine = line('Parent', ax, ...
                    'XData', [cx, xTip], 'YData', [cy, yTip], ...
                    'Color', [0.2 0.6 1.0], 'LineWidth', 2.5, ...
                    'HitTest', 'off', 'Tag', 'BFbeamLine');
            else
                set(store.beamLine, 'XData', [cx, xTip], 'YData', [cy, yTip], 'Visible', 'on');
            end

            % Arrowhead
            headLen = arrowLen * 0.15;
            perp = [-dir(2), dir(1)];
            p1 = [xTip, yTip];
            p2 = p1 - headLen*dir + headLen*0.4*perp;
            p3 = p1 - headLen*dir - headLen*0.4*perp;

            if ~isfield(store, 'beamHead') || ~isgraphics(store.beamHead)
                store.beamHead = fill(ax, [p1(1) p2(1) p3(1)], [p1(2) p2(2) p3(2)], ...
                    [0.2 0.6 1.0], 'EdgeColor', 'none', ...
                    'HitTest', 'off', 'Tag', 'BFbeamHead');
            else
                set(store.beamHead, ...
                    'XData', [p1(1) p2(1) p3(1)], ...
                    'YData', [p1(2) p2(2) p3(2)], 'Visible', 'on');
            end
        end

        function store = hideDASVisuals(store)
            flds = {'beamLine', 'beamHead'};
            for k = 1:numel(flds)
                if isfield(store, flds{k}) && isgraphics(store.(flds{k}))
                    set(store.(flds{k}), 'Visible', 'off');
                end
            end
        end

        % ────────────────────────────────────────────────
        %  Focus visuals: focal point + converging rays
        % ────────────────────────────────────────────────
        function store = drawFocus(ax, rp, store)
            vs = rp.vs;
            Nls = size(rp.x0, 1);

            % Focal point marker
            if ~isfield(store, 'focalMarker') || ~isgraphics(store.focalMarker)
                store.focalMarker = plot(ax, vs(1), vs(2), 'x', ...
                    'MarkerSize', 14, 'LineWidth', 2.5, ...
                    'Color', [0.9 0.2 0.2], 'HitTest', 'off', 'Tag', 'BFfocalMarker');
            else
                set(store.focalMarker, 'XData', vs(1), 'YData', vs(2), 'Visible', 'on');
            end

            % Converging rays: one NaN-separated line for all rays
            xRays = [];
            yRays = [];
            for n = 1:Nls
                xRays = [xRays, rp.x0(n,1), vs(1), NaN]; %#ok<AGROW>
                yRays = [yRays, rp.x0(n,2), vs(2), NaN]; %#ok<AGROW>
            end

            if ~isfield(store, 'focusRays') || ~isgraphics(store.focusRays)
                store.focusRays = line('Parent', ax, ...
                    'XData', xRays, 'YData', yRays, ...
                    'Color', [0.5 0.7 1.0 0.5], 'LineWidth', 1.0, ...
                    'LineStyle', ':', 'HitTest', 'off', 'Tag', 'BFfocusRays');
            else
                set(store.focusRays, 'XData', xRays, 'YData', yRays, 'Visible', 'on');
            end
        end

        function store = hideFocusVisuals(store)
            flds = {'focalMarker', 'focusRays'};
            for k = 1:numel(flds)
                if isfield(store, flds{k}) && isgraphics(store.(flds{k}))
                    set(store.(flds{k}), 'Visible', 'off');
                end
            end
        end

    end
end
