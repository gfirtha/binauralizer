classdef WFSPainter
    % WFSPainter  Draws WFS renderer overlays: reference contour,
    %   amplitude distribution, tapering window, illuminated/tapered zones,
    %   and the intended-field amplitude map A0.
    %
    % Zone visualization uses the zone structs provided by wfs_renderer
    % (rp.zone array with .boundaries and .type), NOT inferred from win_taper.

    methods (Static)

        function store = draw(ax, rp, renderer, setting, store)
            if nargin < 5 || isempty(store), store = struct(); end
            hold(ax,'on');

            

            % ── 1) Reference contour ──
            if setting.showRef
                xRef = rp.ref_pos(:,1);  yRef = rp.ref_pos(:,2);
                if rp.isClosed
                    xRef = [xRef; xRef(1)]; yRef = [yRef; yRef(1)];
                end
                store.wfsRefContour = audioapp.ui.draw.WFSPainter.updateLine(ax, store, 'wfsRefContour', xRef, yRef, ...
                    {'LineStyle','--','LineWidth',1.25,'Color',[0 0 0],'HitTest','off','Tag','WFSRefContour'});
            else
                store = audioapp.ui.draw.WFSPainter.hideField(store, 'wfsRefContour');
            end

            % ── 2) Amplitude ribbon ──
            % Stacking (uistack, ~0.1 s per call) only when an object is
            % created: updating XData/YData keeps the child order.
            if setting.showAmp
                [Xp,Yp] = audioapp.ui.draw.WFSPainter.ribbonPoly(rp.x0(:,1:2), rp.n0(:,1:2), abs(rp.amp), 0.5, rp.isClosed);
                isNew = audioapp.ui.draw.WFSPainter.isMissing(store, 'amplitudeDistribution');
                store.amplitudeDistribution = audioapp.ui.draw.WFSPainter.updateFill(ax, store, 'amplitudeDistribution', Xp, Yp, ...
                    {'FaceColor',[1 1 1]*0.95,'FaceAlpha',0.9,'EdgeColor',[100 100 100]/255,'LineWidth',1.5,'Tag','WFSamplitudeDistribution'});
                if isNew, uistack(store.amplitudeDistribution,'bottom'); end
            else
                store = audioapp.ui.draw.WFSPainter.hideField(store, 'amplitudeDistribution');
            end

            % ── 3) Taper ribbon (drawn ON TOP of amplitude) ──
            if setting.showTaper
                [Xp2,Yp2] = audioapp.ui.draw.WFSPainter.ribbonPoly(rp.x0(:,1:2), rp.n0(:,1:2), abs(rp.win_taper), 0.15, rp.isClosed);
                isNew = audioapp.ui.draw.WFSPainter.isMissing(store, 'windowFunction');
                store.windowFunction = audioapp.ui.draw.WFSPainter.updateFill(ax, store, 'windowFunction', Xp2, Yp2, ...
                    {'FaceColor',[245 154 159]/255,'FaceAlpha',0.5,'EdgeColor',[209 151 137]/255,'LineWidth',1.5,'Tag','WFSwindowFunction'});
                % Place taper above amplitude, but below other elements
                if isNew && setting.showAmp && isfield(store,'amplitudeDistribution') && all(isgraphics(store.amplitudeDistribution))
                    uistack(store.windowFunction,'up');
                end
            else
                store = audioapp.ui.draw.WFSPainter.hideField(store, 'windowFunction');
            end

            % ── 4) Zones — uses rp.zone from the renderer directly ──
            if setting.showZones
                zones = audioapp.ui.draw.WFSPainter.getZones(rp);
                store = audioapp.ui.draw.WFSPainter.drawZones(ax, store, zones, rp, renderer);
            else
                store = audioapp.ui.draw.WFSPainter.hideField(store, 'zoneRays');
                store = audioapp.ui.draw.WFSPainter.hideField(store, 'zonePolys');
            end

            % ── 5) Intended-field amplitude A0 (Part I, Eq. WFS_amp_dist) ──
            if isfield(setting, 'showA0') && setting.showA0 && ismethod(renderer, 'get_intended_amplitude')
                store = audioapp.ui.draw.WFSPainter.drawA0(ax, store, renderer);
            else
                store = audioapp.ui.draw.WFSPainter.hideField(store, 'wfsA0');
            end
        end

        function store = drawA0(ax, store, renderer)
            % Contours of 20 log10 A0 over the visible area: 0 dB (bold) is
            % where the synthesis is amplitude-correct.
            xl = xlim(ax); yl = ylim(ax); n = 70;
            [X, Y] = meshgrid(linspace(xl(1), xl(2), n), linspace(yl(1), yl(2), n));
            A0 = renderer.get_intended_amplitude([X(:) Y(:)]);
            Z = reshape(20*log10(max(A0, 1e-6)), size(X));
            Z(~isfinite(reshape(A0, size(X))) | reshape(A0, size(X)) <= 0) = NaN;
            if isfield(store, 'wfsA0'), delete(store.wfsA0(isgraphics(store.wfsA0))); end
            [~, h1] = contour(ax, X, Y, Z, [-6 -3 -1 1 3 6], 'LineColor', [0.2 0.4 0.8], ...
                'LineStyle', ':', 'ShowText', 'on', 'LabelSpacing', 400, 'HitTest', 'off', 'Tag', 'WFSA0');
            [~, h2] = contour(ax, X, Y, Z, [0 0], 'LineColor', [0.1 0.3 0.8], 'LineWidth', 1.5, ...
                'HitTest', 'off', 'Tag', 'WFSA0');
            store.wfsA0 = [h1 h2];
            xlim(ax, xl); ylim(ax, yl);
        end

        % ═══════════════════════════════════════════════════════════
        %  Zone drawing
        % ═══════════════════════════════════════════════════════════

        function zones = getZones(rp)
            % Return a struct array from the renderer's own zone data.
            if isfield(rp,'zone') && ~isempty(rp.zone)
                zones = rp.zone;
            else
                zones = [];
            end
        end

        function store = drawZones(ax, store, zones, rp, renderer)
            

            if isempty(zones)
                store = audioapp.ui.draw.WFSPainter.hideField(store, 'zoneRays');
                store = audioapp.ui.draw.WFSPainter.hideField(store, 'zonePolys');
                return;
            end

            N = size(rp.x0,1);
            vs = renderer.virtual_source.position(1:2);
            F  = rp.F;

            % ── Boundary rays ──
            % Collect unique boundary indices (each zone has 2, many are shared)
            allBnd = [];
            for k = 1:numel(zones)
                allBnd = [allBnd, zones(k).boundaries(:)']; %#ok<AGROW>
            end
            bIdx = unique(allBnd);
            % Remove invalid
            bIdx(bIdx < 1 | bIdx > N | ~isfinite(bIdx)) = [];

            if ~isfield(store,'zoneRays') || ~iscell(store.zoneRays)
                store.zoneRays = {};
            end
            nRays = numel(bIdx);
            for k = 1:nRays
                [rx, ry] = audioapp.ui.draw.WFSPainter.boundaryRay(bIdx(k), rp, vs, F);
                if k > numel(store.zoneRays) || ~all(isgraphics(store.zoneRays{k}))
                    store.zoneRays{k} = line('Parent',ax,'XData',rx,'YData',ry, ...
                        'LineStyle','--','LineWidth',0.75,'Color',[0 0 0], ...
                        'HitTest','off','Tag',sprintf('wfsZoneRay%d',k));
                else
                    set(store.zoneRays{k}, 'XData',rx, 'YData',ry, 'Visible','on');
                end
            end
            % Hide/delete extras
            for k = nRays+1 : numel(store.zoneRays)
                if all(isgraphics(store.zoneRays{k})), set(store.zoneRays{k},'Visible','off'); end
            end

            % ── Zone fill polygons ──
            if ~isfield(store,'zonePolys') || ~iscell(store.zonePolys)
                store.zonePolys = {};
            end
            nPolys = numel(zones);
            for k = 1:nPolys
                z = zones(k);
                if numel(z.boundaries) ~= 2, continue; end
                i1 = z.boundaries(1);  i2 = z.boundaries(2);
                if any(~isfinite([i1 i2])) || i1<1 || i2<1 || i1>N || i2>N
                    continue;
                end

                [Xp, Yp] = audioapp.ui.draw.WFSPainter.zonePoly(i1, i2, rp, vs, F);
                col = audioapp.ui.draw.WFSPainter.zoneColor(z.type);

                if k > numel(store.zonePolys) || ~all(isgraphics(store.zonePolys{k}))
                    store.zonePolys{k} = fill(ax, Xp, Yp, 1, ...
                        'FaceColor',col,'FaceAlpha',0.4, ...
                        'EdgeColor','none','HitTest','off', ...
                        'Tag',sprintf('wfsZone%d',k));
                    uistack(store.zonePolys{k},'bottom');   % once, at creation
                else
                    set(store.zonePolys{k}, 'XData',Xp, 'YData',Yp, ...
                        'FaceColor',col, 'Visible','on');
                end
            end
            for k = nPolys+1 : numel(store.zonePolys)
                if all(isgraphics(store.zonePolys{k})), set(store.zonePolys{k},'Visible','off'); end
            end
        end

        % ═══════════════════════════════════════════════════════════
        %  Geometry helpers
        % ═══════════════════════════════════════════════════════════

        function [rx, ry] = boundaryRay(idx, rp, vs, F)
            % Single boundary ray at SSD index idx
            len = 1e2;
            kh = rp.kh(idx,1:2);
            if F == 1
                % Ray from virtual source outward along kh
                rx = [vs(1),  vs(1) + len*kh(1)];
                ry = [vs(2),  vs(2) + len*kh(2)];
            else
                % Focused: ray from SSD point backward
                p0 = rp.x0(idx,1:2);
                rx = [p0(1),  p0(1) - len*kh(1)];
                ry = [p0(2),  p0(2) - len*kh(2)];
            end
        end

        function [Xp, Yp] = zonePoly(i1, i2, rp, vs, F)
            % Build a polygon for the zone between SSD indices i1..i2
            % bounded by rays at i1 and i2.
            N = size(rp.x0,1);
            len = 1e2;

            % Arc along SSD from i1 to i2
            if rp.isClosed
                if i2 >= i1
                    idx = i1:i2;
                else
                    idx = [i1:N, 1:i2];
                end
            else
                if i1 > i2, [i1,i2] = deal(i2,i1); end
                idx = i1:i2;
            end
            arc = rp.x0(idx, 1:2);

            % Ray endpoints at the two boundaries
            kh1 = rp.kh(i1,1:2);
            kh2 = rp.kh(i2,1:2);

            if F == 1
                p1far = vs + len*kh1;
                p2far = vs + len*kh2;
                % Polygon: arc along SSD, then ray-end at i2, then ray-end at i1
                Xp = [arc(:,1); p2far(1); p1far(1)];
                Yp = [arc(:,2); p2far(2); p1far(2)];
            else
                p1far = rp.x0(i1,1:2) - len*kh1;
                p2far = rp.x0(i2,1:2) - len*kh2;
                Xp = [arc(:,1); p2far(1); p1far(1)];
                Yp = [arc(:,2); p2far(2); p1far(2)];
            end
        end

        function c = zoneColor(ztype)
            switch lower(ztype)
                case 'illuminated', c = [255 248 174]/255;   % warm yellow
                case 'tapered',     c = [236 227 128]/255;   % darker yellow
                otherwise,          c = [0.9 0.9 0.9];      % grey
            end
        end

        % ═══════════════════════════════════════════════════════════
        %  Ribbon (amplitude / taper) polygon builder
        % ═══════════════════════════════════════════════════════════

        function [Xp, Yp] = ribbonPoly(x0, n0, vals, scale, isClosed)
            % Build a filled polygon that shows the distribution of |vals|
            % as an outward ribbon along the SSD.
            ma = max(abs(vals(:)));
            if ma < eps, normA = zeros(size(vals)); else, normA = abs(vals)/ma; end

            wx = x0(:,1) + normA.*n0(:,1)*scale;
            wy = x0(:,2) + normA.*n0(:,2)*scale;

            if isClosed
                % Find a seam where normA is near zero to avoid self-intersection
                ix0 = find(normA < 1e-10, 1, 'first');
                if isempty(ix0), ix0 = 1; end
                Xp = [x0(ix0:end,1); x0(1:ix0,1); flipud(wx(1:ix0)); flipud(wx(ix0:end))];
                Yp = [x0(ix0:end,2); x0(1:ix0,2); flipud(wy(1:ix0)); flipud(wy(ix0:end))];
            else
                Xp = [x0(:,1); flipud(wx)];
                Yp = [x0(:,2); flipud(wy)];
            end
        end

        % ═══════════════════════════════════════════════════════════
        %  Graphics utilities
        % ═══════════════════════════════════════════════════════════

        function h = updateLine(ax, store, field, x, y, props)
            if ~isfield(store,field) || ~all(isgraphics(store.(field)))
                h = line('Parent',ax,'XData',x,'YData',y,props{:});
            else
                h = store.(field);
                set(h, 'XData',x, 'YData',y, 'Visible','on');
            end
        end

        function h = updateFill(ax, store, field, x, y, props)
            if ~isfield(store,field) || ~all(isgraphics(store.(field)))
                h = fill(ax, x, y, 1, props{:});
            else
                h = store.(field);
                set(h, 'XData',x, 'YData',y, 'Visible','on');
            end
        end

        function tf = isMissing(store, field)
            tf = ~isfield(store, field) || ~all(isgraphics(store.(field)));
        end

        function store = hideField(store, field)
            if ~isfield(store, field), return; end
            h = store.(field);
            if iscell(h)
                for k = 1:numel(h)
                    if all(isgraphics(h{k})), set(h{k},'Visible','off'); end
                end
            elseif all(isgraphics(h))
                set(h,'Visible','off');
            end
        end

    end
end
