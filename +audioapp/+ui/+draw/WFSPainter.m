classdef WFSPainter
    % Painter for WFS renderer visuals. Supports variable zones (1..3).

    methods (Static)
        function store = draw(ax, rp, renderer, setting, store)
            if nargin < 5 || isempty(store), store = struct(); end
            hold(ax,'on');

            % 1) Reference contour
            if setting.WFSrenderer.showRef
                xRef = rp.ref_pos(:,1);  yRef = rp.ref_pos(:,2);
                if isfield(rp,'isClosed') && rp.isClosed
                    xRef = [xRef; xRef(1)]; yRef = [yRef; yRef(1)];
                end
                store.wfsRefContour = WFSPainter.updateLine(ax, store, 'wfsRefContour', xRef, yRef, ...
                    {'LineStyle','--','LineWidth',1.25,'Color',[0 0 0],'HitTest','off','Tag','WFSRefContour'});
            else
                store = WFSPainter.hideIfExists(store, 'wfsRefContour');
            end

            % 2) Amplitude ribbon
            if setting.WFSrenderer.showAmp
                a = abs(rp.amp); ma = max(a(:)); normA = (ma>0).*a/max(ma,eps);
                scale = 0.5;
                win_x = rp.x0(:,1) + normA.*rp.n0(:,1)*scale;
                win_y = rp.x0(:,2) + normA.*rp.n0(:,2)*scale;
                if isfield(rp,'isClosed') && rp.isClosed
                    % close without self-intersection (use seam where normA~0)
                    ix0 = find(normA==0,1,'first'); if isempty(ix0), ix0 = 1; end
                    Xpoly = [rp.x0(ix0:end,1); rp.x0(1:ix0,1); flipud(win_x(1:ix0)); flipud(win_x(ix0:end))];
                    Ypoly = [rp.x0(ix0:end,2); rp.x0(1:ix0,2); flipud(win_y(1:ix0)); flipud(win_y(ix0:end))];
                else
                    Xpoly = [rp.x0(:,1); flipud(win_x)];
                    Ypoly = [rp.x0(:,2); flipud(win_y)];
                end
                store.amplitudeDistribution = WFSPainter.updateFill(ax, store, 'amplitudeDistribution', Xpoly, Ypoly, ...
                    {'FaceColor',[1 1 1]*0.95,'FaceAlpha',0.9,'EdgeColor',[100,100,100]/255,'LineWidth',1.5,'Tag','WFSamplitudeDistribution'});
                uistack(store.amplitudeDistribution,'bottom');
            else
                store = WFSPainter.hideIfExists(store, 'amplitudeDistribution');
            end

            % 3) Taper ribbon
            if setting.WFSrenderer.showTaper
                a = abs(rp.win_taper); ma = max(a(:)); normA = (ma>0).*a/max(ma,eps);
                scale = 0.15;
                win_x2 = rp.x0(:,1) + normA.*rp.n0(:,1)*scale;
                win_y2 = rp.x0(:,2) + normA.*rp.n0(:,2)*scale;
                if isfield(rp,'isClosed') && rp.isClosed
                    ix0 = find(normA==0,1,'first'); if isempty(ix0), ix0 = 1; end
                    Xpoly2 = [rp.x0(ix0:end,1); rp.x0(1:ix0,1); flipud(win_x2(1:ix0)); flipud(win_x2(ix0:end))];
                    Ypoly2 = [rp.x0(ix0:end,2); rp.x0(1:ix0,2); flipud(win_y2(1:ix0)); flipud(win_y2(ix0:end))];
                else
                    Xpoly2 = [rp.x0(:,1); flipud(win_x2)];
                    Ypoly2 = [rp.x0(:,2); flipud(win_y2)];
                end
                store.windowFunction = WFSPainter.updateFill(ax, store, 'windowFunction', Xpoly2, Ypoly2, ...
                    {'FaceColor',[245,154,159]/255 ,'FaceAlpha',0.5,'EdgeColor',[209,151,137]/255,'LineWidth',1.5,'Tag','WFSwindowFunction'});
                uistack(store.windowFunction,'bottom');
            else
                store = WFSPainter.hideIfExists(store, 'windowFunction');
            end

            % 4) Zones (variable count)
            if setting.WFSrenderer.showZones
                segs = WFSPainter.segmentizeWinTaper(rp.win_taper(:), logical(rp.isClosed));
                % Boundary rays (at segment boundaries)
                bIdx = unique([segs.startIdx]);   % each segment start is a boundary
                store.zoneLines = WFSPainter.drawBoundaryRays(ax, store, bIdx, rp, renderer);

                % Taper polygons (for segments with type='taper')
                taperSegs = segs(strcmp({segs.type},'taper'));
                store.zonePolys = WFSPainter.drawTaperPolys(ax, store, taperSegs, rp, renderer);
            else
                store = WFSPainter.hideIfExists(store, 'zoneLines');
                store = WFSPainter.hideIfExists(store, 'zonePolys');
            end
        end

        % --------- segmenter: 1..3 zones from win_taper ----------
        function segs = segmentizeWinTaper(w, isClosed)
            eps1 = 1e-12;
            % classify each sample
            cls = zeros(size(w));                % 0:off, 1:illum, 2:taper
            cls(w < eps1) = 0;
            cls(abs(w-1) <= eps1) = 1;
            cls(w > eps1 & abs(w-1) > eps1) = 2;

            N = numel(w);
            if N==0, segs = struct('type',{},'startIdx',{},'endIdx',{}); return; end

            if isClosed
                % rotate so we start at a boundary
                k0 = find([true; diff(cls)~=0],1,'first'); if isempty(k0), k0 = 1; end
                cls = cls([k0:N, 1:k0-1]);
                perm = [k0:N, 1:k0-1];   % to map back
            else
                perm = 1:N;
            end

            % boundaries in this (possibly rotated) order
            b = [1; find(diff(cls)~=0)+1; N+1];
            segs = struct('type',{},'startIdx',{},'endIdx',{});
            for i=1:numel(b)-1
                ii = b(i); jj = b(i+1)-1;
                if ii>jj, continue; end
                typ = cls(ii);
                if typ==0, tstr='off'; elseif typ==1, tstr='illum'; else, tstr='taper'; end
                % map indices back to original order
                si = perm(ii); ei = perm(jj);
                segs(end+1) = struct('type',tstr,'startIdx',si,'endIdx',ei); %#ok<AGROW>
            end

            % merge wrap-around if needed (closed)
            if isClosed && ~isempty(segs)
                if segs(1).type == segs(end).type
                    segs(1).startIdx = segs(end).startIdx;
                    segs(end) = [];
                end
            end
        end

        % --------- boundary rays (variable count) ----------
        function zoneLines = drawBoundaryRays(ax, store, bIdx, rp, renderer)
            % Make as many lines as boundaries; reuse handles if possible
            L = numel(bIdx);
            if ~isfield(store,'zoneLines') || ~iscell(store.zoneLines)
                store.zoneLines = {};
            end
            zoneLines = store.zoneLines;

            % create/update
            for k=1:L
                i = bIdx(k);
                if rp.F==1
                    % ray from virtual source along kh(i)
                    p0 = renderer.virtual_source.position(1:2);
                    p1 = p0 + rp.kh(i,1:2)*1e2;
                else
                    % ray from boundary point against kh(i)
                    p0 = rp.x0(i,1:2);
                    p1 = p0 - rp.kh(i,1:2)*1e2;
                end
                if k>numel(zoneLines) || ~isgraphics(zoneLines{k})
                    zoneLines{k} = plot(ax, [p0(1) p1(1)], [p0(2) p1(2)], ...
                        'LineStyle','--','LineWidth',0.5,'Color',[0 0 0], ...
                        'HitTest','off','Tag',sprintf('wfsZoneRay%d',k));
                else
                    set(zoneLines{k}, 'XData',[p0(1) p1(1)], 'YData',[p0(2) p1(2)], 'Visible','on');
                end
            end
            % delete any extras
            for k=L+1:numel(zoneLines)
                if isgraphics(zoneLines{k}), delete(zoneLines{k}); end
            end
            zoneLines = zoneLines(1:L);
        end

        % --------- taper polygons (variable count) ----------
        function zonePolys = drawTaperPolys(ax, store, segs, rp, renderer)
            if ~isfield(store,'zonePolys') || ~iscell(store.zonePolys)
                store.zonePolys = {};
            end
            zonePolys = store.zonePolys;

            L = numel(segs);
            for k=1:L
                s = segs(k);
                idx = WFSPainter.arcIndices(s.startIdx, s.endIdx, size(rp.x0,1), logical(rp.isClosed));
                arc = rp.x0(idx,1:2);

                if rp.F==1
                    % close with two kh points from both ends extended outward
                    p1 = rp.x0(s.startIdx,1:2) + rp.kh(s.startIdx,1:2)*1e2;
                    p2 = rp.x0(s.endIdx,  1:2) + rp.kh(s.endIdx,  1:2)*1e2;
                    X = [arc(:,1); p2(1); p1(1)];
                    Y = [arc(:,2); p2(2); p1(2)];
                else
                    % close with the virtual source (triangle)
                    vs = renderer.virtual_source.position(1:2);
                    p1 = rp.x0(s.startIdx,1:2) - rp.kh(s.startIdx,1:2)*1e2;
                    p2 = rp.x0(s.endIdx,  1:2) - rp.kh(s.endIdx,  1:2)*1e2;
                    X = [vs(1); p2(1); p1(1)];
                    Y = [vs(2); p2(2); p1(2)];
                end

                if k>numel(zonePolys) || ~isgraphics(zonePolys{k})
                    zonePolys{k} = fill(ax, X, Y, 1, ...
                        'FaceColor',[255,248,174]/255,'FaceAlpha',0.5, ...
                        'EdgeColor','none','Tag',sprintf('taperRegion%d',k));
                else
                    set(zonePolys{k}, 'XData',X, 'YData',Y, 'Visible','on');
                end
                uistack(zonePolys{k},'bottom');
            end

            % delete any extras
            for k=L+1:numel(zonePolys)
                if isgraphics(zonePolys{k}), delete(zonePolys{k}); end
            end
            zonePolys = zonePolys(1:L);
        end

        % --------- utils ----------
        function idx = arcIndices(iStart, iEnd, N, isClosed)
            if isClosed
                if iStart <= iEnd
                    idx = iStart:iEnd;
                else
                    idx = [iStart:N, 1:iEnd];
                end
            else
                idx = iStart:iEnd;
            end
        end

        function h = updateLine(ax, store, field, x, y, props)
            if ~isfield(store,field) || ~isgraphics(store.(field))
                h = plot(ax, x, y, props{:});
            else
                h = store.(field);
                set(h, 'XData', x, 'YData', y, 'Visible','on');
            end
        end

        function h = updateFill(ax, store, field, x, y, props)
            if ~isfield(store,field) || ~isgraphics(store.(field))
                h = fill(ax, x, y, 1, props{:});
            else
                h = store.(field);
                set(h, 'XData', x, 'YData', y, 'Visible','on');
            end
        end

        function store = hideIfExists(store, field)
            if isfield(store,field)
                h = store.(field);
                if iscell(h)
                    for k=1:numel(h)
                        if isgraphics(h{k}), set(h{k},'Visible','off'); end
                    end
                elseif isgraphics(h)
                    set(h,'Visible','off');
                end
            end
        end
    end
end
