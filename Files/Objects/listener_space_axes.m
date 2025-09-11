classdef listener_space_axes < handle
    %LISTENER_SPACE_AXES (optimized layering + object reuse)
    % - Stable z-order via 3 hggroups: bg (zones), mid (amp/taper), fg (lines/labels)
    % - Reuse graphics; no per-frame deletes/uistack
    % - Batched dashed zone rays in a single line (NaN-separated)
    % - Less transparency/edges; OpenGL renderer + childorder sort

    properties
        main_axes
        room
        receiver
        virtual_source_points
        loudspeaker_points
        binaural_source_points
        renderer_illustration

        parentFig
        zoomPanel
        btnMinus
        btnPlus
        zoomFactorIn  = 0.80;  % <1 zooms in
        zoomFactorOut = 1.25;  % >1 zooms out
        zoomFactor

        rendererSetting
        lastDrawTic = [];     % for throttling
        minDrawDt = 1/15;     % max 10 FPS for visuals

        simulatorHandle
    end

    methods
        function obj = listener_space_axes(varargin)
            if ~isempty(varargin) && isa(varargin{1},'matlab.graphics.axis.Axes')
                ax = varargin{1};
                obj.main_axes = ax;
            else
                obj.main_axes = axes('Units','normalized','Position',[0.08 0.18 0.88 0.78]);
                ax = obj.main_axes;
            end

            axis(obj.main_axes,'equal'); axis(obj.main_axes,'tight');
            obj.zoomFactor = 1;
            obj.lastDrawTic = tic;
            obj.renderer_illustration = struct();
            grid(obj.main_axes,'on');
            set(obj.main_axes,'ButtonDownFcn',[]);   % don’t rely on gca
            xlabel(obj.main_axes, 'x -> [m]');
            ylabel(obj.main_axes, 'y -> [m]');
            title('Main view')

            % Fast renderer defaults
            set(ancestor(ax,'figure'),'Renderer','opengl','GraphicsSmoothing','off');
            set(ax,'SortMethod','childorder','Clipping','on');

            if numel(varargin) > 1
                obj.room = varargin{2};
                hold(obj.main_axes,'on');
                obj.room.drawRoom(obj.main_axes);
                view(obj.main_axes, 2);
                xlim(obj.main_axes,[min(obj.room.vertices(:,2))-0.25 , max(obj.room.vertices(:,2))+0.25]);
                ylim(obj.main_axes,[min(obj.room.vertices(:,3))-0.25 , max(obj.room.vertices(:,3))+0.25]);
            else
                xlim(obj.main_axes,[-4,4]);
                ylim(obj.main_axes,[-4,4]);
            end

            % ---- Zoom UI (bottom-right of this axes) ----
            ax = obj.main_axes;
            obj.parentFig = ancestor(ax,'figure');

            % Put the floating controls into the SAME parent as the axes
            parentContainer = ax.Parent;
            bgc = ancestor(parentContainer,'figure').Color;

            obj.zoomPanel = uipanel( ...
                'Parent', parentContainer, ...
                'Units','pixels','BorderType','none', ...
                'BackgroundColor', bgc, ...
                'Tag','ZoomPanel');

            obj.btnMinus = uicontrol('Parent',obj.zoomPanel, ...
                'Style','pushbutton','String','–','Units','pixels', ...
                'Position',[0 0 28 28], 'TooltipString','Zoom out', ...
                'FontWeight','bold', ...
                'Callback',@(~,~)obj.zoomAxes(obj.zoomFactorOut), ...
                'BackgroundColor',[0.94 0.94 0.94]);

            obj.btnPlus = uicontrol('Parent',obj.zoomPanel, ...
                'Style','pushbutton','String','+','Units','pixels', ...
                'Position',[28 0 28 28], 'TooltipString','Zoom in', ...
                'FontWeight','bold', ...
                'Callback',@(~,~)obj.zoomAxes(obj.zoomFactorIn), ...
                'BackgroundColor',[0.94 0.94 0.94]);

            % Reposition on resize of BOTH the figure and the axes' parent panel
            prevFigSCF = obj.parentFig.SizeChangedFcn;
            obj.parentFig.SizeChangedFcn = @(src,evt) obj.onFigureResized(prevFigSCF, src, evt);

            try
                parentPrevSCF = parentContainer.SizeChangedFcn;
                parentContainer.SizeChangedFcn = @(src,evt) obj.onParentResized(parentPrevSCF, src, evt);
            catch
            end
            try
                prevAxSCF = obj.main_axes.SizeChangedFcn;
                obj.main_axes.SizeChangedFcn = @(src,evt) obj.onAxesResized(prevAxSCF, src, evt);
            catch
            end

            obj.positionZoomPanel();   % initial placement

            obj.rendererSetting.WFSrenderer = struct('showRef',true,'showAmp',true,'showTaper',true,'showZones',true);
        end

        function onFigureResized(obj, prevFcn, src, evt)
            if ~isempty(prevFcn)
                try, feval(prevFcn, src, evt); catch, end
            end
            obj.positionZoomPanel();
        end

        function onAxesResized(obj, prevFcn, src, evt)
            if ~isempty(prevFcn)
                try, feval(prevFcn, src, evt); catch, end
            end
            obj.positionZoomPanel();
        end

        function delete(obj)
            if ~isempty(obj.zoomPanel) && isvalid(obj.zoomPanel), delete(obj.zoomPanel); end
        end

        function hAx = get_main_axes(obj)
            hAx = obj.main_axes;
        end

        function receiver = draw_head(obj,pos,orientation,R)
            N = 20;
            fi = (0:2*pi/N:2*pi*(1-1/N));
            x_head = cos(fi)*R;
            y_head = sin(fi)*R;
            fi_nose =  [10 5   0    -5  -10 ];
            A_nose = R*[1  1.15 1.2 1.15   1 ];
            x_nose = cosd(linspace(fi_nose(1),fi_nose(end),N)).*interp1(fi_nose,A_nose,linspace(fi_nose(1),fi_nose(end),N),'spline');
            y_nose = sind(linspace(fi_nose(1),fi_nose(end),N)).*interp1(fi_nose,A_nose,linspace(fi_nose(1),fi_nose(end),N),'spline');
            fi_lear =  fliplr([18.5 18 16 0 -5 -15] + 90);
            A_lear  = fliplr(R*[1 1.08 1.1 1.04 1.06 1]);
            x_lear = cosd(linspace(fi_lear(1),fi_lear(end),N)).*interp1(fi_lear,A_lear,linspace(fi_lear(1),fi_lear(end),N));
            y_lear = sind(linspace(fi_lear(1),fi_lear(end),N)).*interp1(fi_lear,A_lear,linspace(fi_lear(1),fi_lear(end),N));
            fi_rear =  -[18.5 18 16 0 -5 -15]-90;
            A_rear =  R*[1 1.08 1.1 1.04 1.06 1];
            x_rear = cosd(linspace(fi_rear(1),fi_rear(end),N)).*interp1(fi_rear,A_rear,linspace(fi_rear(1),fi_rear(end),N));
            y_rear = sind(linspace(fi_rear(1),fi_rear(end),N)).*interp1(fi_rear,A_rear,linspace(fi_rear(1),fi_rear(end),N));
            x_torso = cos(fi)*R*0.7-R/7;
            y_torso = sin(fi)*R*1.7;
            x_rec = [x_torso;x_head;x_lear;x_rear;x_nose]';
            y_rec = [y_torso;y_head;y_lear;y_rear;y_nose]';
            x_rec = x_rec - mean(mean(x_rec));
            y_rec = y_rec - mean(mean(y_rec));
            c = [37 160 217;
                77  41  14;
                255 206 180;
                255 206 180;
                255 206 180]/255;
            receiver = patch(obj.main_axes, x_rec+pos(1) ,y_rec + pos(2),[0;1;1;1;1], 'Tag','Receiver');
            set(receiver,'FaceVertexCData',c,'FaceLighting','none');
            receiver.UserData = struct( 'Label', 1,...
                'Origin', [ mean(receiver.Vertices(:,1)), mean(receiver.Vertices(:,2))  ],...
                'Orientation', 0 );
            rotate(receiver, [0,0,1], asind(orientation(2)), [receiver.UserData.Origin(1),receiver.UserData.Origin(2),0]);
            receiver.UserData.Orientation = asind(orientation(2));
            obj.receiver = receiver;
        end

        function update_receiver( obj, receiver_position, receiver_orientation )
            if isempty(obj.receiver) || ~isvalid(obj.receiver), return; end
            v = [cosd(obj.receiver.UserData.Orientation), sind(obj.receiver.UserData.Orientation)];
            w = receiver_orientation(1:2);
            dfi = atan2(w(2)*v(1)-w(1)*v(2), w(1)*v(1)+w(2)*v(2));

            rotate(obj.receiver,[0 0 1], dfi*180/pi,...
                [obj.receiver.UserData.Origin(1),obj.receiver.UserData.Origin(2),0]);
            obj.receiver.UserData.Orientation = obj.receiver.UserData.Orientation + dfi*180/pi;

            dx = receiver_position(1:2) - obj.receiver.UserData.Origin;
            obj.receiver.UserData.Origin = obj.receiver.UserData.Origin + dx;
            obj.receiver.Vertices(:,1:2) = obj.receiver.Vertices(:,1:2) + dx;
        end

        function loudspeaker = draw_loudspeaker(obj,pos,R,orientation,idx)
            x1 = [-1.8 -1.8 -1 -1]'*R;
            y1 = [-1    1   1 -1]'*R;
            x2 = [-1 -1  0  0]'*R;
            y2 = [-1  1 1.5 -1.5]'*R;
            x = [x1,x2];
            y = [y1,y2];
            x = x - mean(mean(x));
            y = y - mean(mean(y));

            c = [0.2 0.2 0.2;
                0.5 0.5 0.5];
            loudspeaker = patch(obj.main_axes, x + pos(1), y + pos(2),[0;1], 'Tag',sprintf('loudspeaker_%i',idx));
            set(loudspeaker,'FaceVertexCData',c,'FaceLighting','none');
            loudspeaker.UserData = struct( 'Label', idx,...
                'Origin', [ mean(loudspeaker.Vertices(:,1)), mean(loudspeaker.Vertices(:,2))  ],...
                'Orientation', orientation );
            rotate(loudspeaker,[0 0 1], orientation,...
                [loudspeaker.UserData.Origin(1),loudspeaker.UserData.Origin(2),0]);
        end

        function virtual_source = draw_virtual_source(obj,pos,orientation,idx)
            N = 20;
            R = [180,40,45,50];
            A = [0.03,0.08,0.12,0.16]*1.5;
            w = 0.002;
            x_source = A(1)*cosd(linspace(-R(1),R(1),N));
            y_source = A(1)*sind(linspace(-R(1),R(1),N));
            x_w1 = [A(2)*cosd(linspace(-R(2),R(2),N/2)) (A(2)-w)*cosd(linspace(R(2),-R(2),N/2))];
            y_w1 = [A(2)*sind(linspace(-R(2),R(2),N/2)) (A(2)-w)*sind(linspace(R(2),-R(2),N/2))];
            x_w2 = [A(3)*cosd(linspace(-R(3),R(3),N/2)) (A(3)-w)*cosd(linspace(R(3),-R(3),N/2))];
            y_w2 = [A(3)*sind(linspace(-R(3),R(3),N/2)) (A(3)-w)*sind(linspace(R(3),-R(3),N/2))];
            x_w3 = [A(4)*cosd(linspace(-R(4),R(4),N/2)) (A(4)-w)*cosd(linspace(R(4),-R(4),N/2))];
            y_w3 = [A(4)*sind(linspace(-R(4),R(4),N/2)) (A(4)-w)*sind(linspace(R(4),-R(4),N/2))];

            x = [x_source;x_w1;x_w2;x_w3]';
            y = [y_source;y_w1;y_w2;y_w3]';

            x = x - mean(mean(x));
            y = y - mean(mean(y));

            c = [255 0 0;
                0 0 0;
                0 0 0;
                0 0 0]/255;
            virtual_source = patch(obj.main_axes, x + pos(1), y + pos(2),[1;0;0;0],'Tag',sprintf('s_%i',idx));
            set(virtual_source,'FaceVertexCData',c,'FaceLighting','none');

            virtual_source.UserData = struct( 'Label', idx,...
                'Origin', [ mean(virtual_source.Vertices(:,1)), mean(virtual_source.Vertices(:,2))  ],...
                'Orientation', orientation );

            rotate(virtual_source,[0 0 1], orientation,[pos(1),pos(2),0]);
            virtual_source.UserData.text = text(obj.main_axes, pos(1)+0.1,pos(2)+0.1,virtual_source.Tag);
            addlistener(virtual_source,'ObjectBeingDestroyed', @(src,~) cleanupPatch(src));

            function cleanupPatch(src)
                if isfield(src.UserData,'text')
                    h = src.UserData.text;
                    if isgraphics(h) && isvalid(h)
                        delete(h);
                    end
                end
            end
        end

        function positionZoomPanel(obj)
            if isempty(obj.zoomPanel) || ~isvalid(obj.zoomPanel) || ...
                    isempty(obj.main_axes) || ~isvalid(obj.main_axes)
                return;
            end

            % compute plot box in PIXELS (accounts for axis equal / aspect)
            pb = obj.getPlotBoxPixels();

            margin = 1; panelW = 56; panelH = 28;
            obj.zoomPanel.Position = [pb(1) + pb(3) - panelW - margin, ...
                pb(2) + margin, ...
                panelW, panelH];
            uistack(obj.zoomPanel,'top');
        end
        function pb = getPlotBoxPixels(obj)
            % Returns [x y w h] of the actual data plot box in PIXELS,
            % i.e., the square/rect that respects aspect ratio (axis equal).
            ax = obj.main_axes;

            % Axes inner rectangle (in pixels)
            oldU = ax.Units; ax.Units = 'pixels';
            pos  = ax.Position;                % axes box (including labels)
            ti   = ax.TightInset;              % label/tick margins
            ax.Units = oldU;

            inner = [pos(1)+ti(1), pos(2)+ti(2), ...
                pos(3)-ti(1)-ti(3), pos(4)-ti(2)-ti(4)];

            % If something is off, fall back to inner
            if any(~isfinite(inner)) || inner(3) <= 0 || inner(4) <= 0
                pb = inner; return;
            end

            % Data aspect ratio from current limits
            xl = xlim(ax); yl = ylim(ax);
            if any(~isfinite([xl yl])) || diff(yl) == 0
                pb = inner; return;
            end
            dataAR = diff(xl) / diff(yl);      % width / height ratio

            % Fit a rectangle of aspect dataAR inside the inner box
            w = inner(3); h = inner(4);
            if w/h > dataAR
                % too wide → letterbox horizontally
                w2 = h * dataAR;
                x  = inner(1) + (w - w2)/2;
                y  = inner(2);
                pb = [x, y, w2, h];
            else
                % too tall → letterbox vertically
                h2 = w / dataAR;
                x  = inner(1);
                y  = inner(2) + (h - h2)/2;
                pb = [x, y, w, h2];
            end
        end

        function zoomAxes(obj, factor)
            ax = obj.main_axes;
            xl = xlim(ax); yl = ylim(ax);
            cx = mean(xl); cy = mean(yl);
            hw = (diff(xl)/2)*factor;
            hh = (diff(yl)/2)*factor;
            xlim(ax, [cx - hw, cx + hw]);
            ylim(ax, [cy - hh, cy + hh]);
            obj.zoomFactor = obj.zoomFactor .* factor;
        end

        function request_draw_visuals(obj, renderer)
            if toc(obj.lastDrawTic) < obj.minDrawDt
                return; % skip this frame
            end
            obj.lastDrawTic = tic;
            obj.draw_renderer_properties(renderer);
            drawnow limitrate nocallbacks
        end

        function obj = draw_renderer_properties(obj, renderer)
            switch lower(class(renderer))
                case 'wfs_renderer'
                    obj.draw_wfs_properties(renderer);
                case 'vbap_renderer'
                case 'nfc_hoa_renderer'
                case 'direct_playback'
                case 'dolby'
            end
        end

        function obj = draw_wfs_properties(obj, renderer)
            rp = renderer.get_renderer_visuals();   % x0, n0, amp, ref_pos, isClosed, etc.
            ax = obj.main_axes;
            hold(ax,'on');

            % Ensure layer groups exist
            L = obj.ensureLayers();

            % --- Reference contour ---
            if obj.rendererSetting.WFSrenderer.showRef
                xRef = rp.ref_pos(:,1);  yRef = rp.ref_pos(:,2);
                if isfield(rp,'isClosed') && rp.isClosed
                    xRef = [xRef; xRef(1)]; yRef = [yRef; yRef(1)];
                end
                obj.renderer_illustration.wfsRefContour = ...
                    createOrUpdateLine(obj.renderer_illustration, 'wfsRefContour', L.fg, xRef, yRef, ...
                    {'LineStyle','--','LineWidth',1.25,'Color',[0 0 0],'HitTest','off','Tag','WFSRefContour'});
            end

            % --- Amplitude distribution fill (MID layer) ---
            if obj.rendererSetting.WFSrenderer.showAmp
                a = abs(rp.amp);
                ma = max(a(:));
                if isempty(ma) || ma==0, normA = zeros(size(a)); else, normA = a./ma; end

                scale = 0.5;
                win_x = rp.x0(:,1) + normA.*rp.n0(:,1)*scale;
                win_y = rp.x0(:,2) + normA.*rp.n0(:,2)*scale;

                if isfield(rp,'isClosed') && rp.isClosed
                    ix0 = 2*find(normA==0,1,'first');
                    Xpoly = [rp.x0(ix0+1:end,1); rp.x0(1:ix0,1); rp.x0(ix0+1,1); flipud(win_x(1:ix0+1)); flipud(win_x(ix0+1:end))];
                    Ypoly = [rp.x0(ix0+1:end,2); rp.x0(1:ix0,2); rp.x0(ix0+1,2); flipud(win_y(1:ix0+1)); flipud(win_y(ix0+1:end))];
                else
                    Xpoly = [rp.x0(:,1); flipud(win_x)];
                    Ypoly = [rp.x0(:,2); flipud(win_y)];
                end

                obj.renderer_illustration.amplitudeDistribution = ...
                    createOrUpdatePatch(obj.renderer_illustration, 'amplitudeDistribution', L.mid, Xpoly, Ypoly, ...
                    {'FaceColor',[1 1 1]*0.95,'FaceAlpha',0.90,'EdgeColor',[100,100,100]/255,'LineWidth',1.5,'Tag','WFSamplitudeDistribution'});
            end

            % --- Window/taper fill (MID layer) ---
            if obj.rendererSetting.WFSrenderer.showTaper
                a = abs(rp.win_taper);
                ma = max(a(:));
                if isempty(ma) || ma==0, normA = zeros(size(a)); else, normA = a./ma; end
                scale = 0.15;
                win_x2 = rp.x0(:,1) + normA.*rp.n0(:,1)*scale;
                win_y2 = rp.x0(:,2) + normA.*rp.n0(:,2)*scale;

                if isfield(rp,'isClosed') && rp.isClosed
                    ix0 = 2*find(normA==0,1,'first');
                    Xpoly2 = [rp.x0(ix0+1:end,1); rp.x0(1:ix0,1); rp.x0(ix0+1,1); flipud(win_x2(1:ix0+1)); flipud(win_x2(ix0+1:end))];
                    Ypoly2 = [rp.x0(ix0+1:end,2); rp.x0(1:ix0,2); rp.x0(ix0+1,2); flipud(win_y2(1:ix0+1)); flipud(win_y2(ix0+1:end))];
                else
                    Xpoly2 = [rp.x0(:,1); flipud(win_x2)];
                    Ypoly2 = [rp.x0(:,2); flipud(win_y2)];
                end

                obj.renderer_illustration.windowFunction = ...
                    createOrUpdatePatch(obj.renderer_illustration, 'windowFunction', L.mid, Xpoly2, Ypoly2, ...
                    {'FaceColor',[245,154,159]/255 ,'FaceAlpha',0.50,'EdgeColor',[209,151,137]/255,'LineWidth',1.5,'Tag','WFSwindowFunction'});
            end

            % --- Zones (BG layer) ---
            if obj.rendererSetting.WFSrenderer.showZones
                obj.draw_wfs_zones(rp, renderer);
            else
                % hide if present
                if isfield(obj.renderer_illustration,'zoneFills')
                    for kk = 1:obj.renderer_illustration.zoneFillsN
                        if isgraphics(obj.renderer_illustration.zoneFills(kk))
                            set(obj.renderer_illustration.zoneFills(kk),'Visible','off');
                        end
                    end
                end
                if isfield(obj.renderer_illustration,'zoneRays') && isgraphics(obj.renderer_illustration.zoneRays)
                    set(obj.renderer_illustration.zoneRays,'XData',NaN,'YData',NaN);
                end
            end

            if isprop(obj,'simulatorHandle') && ~isempty(obj.simulatorHandle) && isvalid(obj.simulatorHandle)
                obj.simulatorHandle.refreshRendererCopies();
            end

            % ---- nested creators (parent-aware, reuse) ----
            function h = createOrUpdateLine(store, field, parent, x, y, props)
                if ~isfield(store,field) || ~isvalid(store.(field))
                    h = line('Parent',parent,'XData',x,'YData',y,props{:});
                    obj.renderer_illustration.(field) = h;
                else
                    h = store.(field);
                    set(h, 'Parent',parent, 'XData', x, 'YData', y, 'Visible','on');
                end
            end
            function h = createOrUpdatePatch(store, field, parent, x, y, props)
                if ~isfield(store,field) || ~isvalid(store.(field))
                    h = patch('Parent', parent, 'XData', x, 'YData', y, 'CData', 1, props{:});
                    obj.renderer_illustration.(field) = h;
                else
                    h = store.(field);
                    set(h, 'Parent',parent, 'XData', x, 'YData', y, 'Visible','on');
                end
            end
        end

        function settings = get_settings(obj)
            settings = obj.rendererSetting;
        end

        function obj = set_settings(obj,setting)
            obj.rendererSetting = setting;
        end

        function obj = clearObject(obj)
            if ~isempty(obj.receiver) && isvalid(obj.receiver), delete(obj.receiver); end
            obj.receiver = [];

            if ~isempty(obj.virtual_source_points)
                for n = 1:numel(obj.virtual_source_points)
                    if isgraphics(obj.virtual_source_points{n}), delete(obj.virtual_source_points{n}); end
                end
            end
            obj.virtual_source_points = {};

            if ~isempty(obj.loudspeaker_points)
                for n = 1:numel(obj.loudspeaker_points)
                    if isgraphics(obj.loudspeaker_points{n}), delete(obj.loudspeaker_points{n}); end
                end
            end
            obj.loudspeaker_points = {};

            if ~isempty(obj.binaural_source_points)
                for n = 1:numel(obj.binaural_source_points)
                    if isgraphics(obj.binaural_source_points{n}), delete(obj.binaural_source_points{n}); end
                end
            end
            obj.binaural_source_points = {};

            if ~isempty(obj.renderer_illustration)
                fns = fieldnames(obj.renderer_illustration);
                for k = 1:numel(fns)
                    if strcmp(fns{k},'layers'), continue; end % keep layer groups alive
                    h = obj.renderer_illustration.(fns{k});
                    if isgraphics(h)
                        try
                            if isvalid(h), delete(h); end
                        catch
                        end
                    elseif isa(h,'matlab.graphics.primitive.Patch')
                        try, delete(h); catch, end
                    end
                end
                % remove non-layer fields
                keepLayers = struct();
                if isfield(obj.renderer_illustration,'layers')
                    keepLayers.layers = obj.renderer_illustration.layers;
                end
                obj.renderer_illustration = keepLayers;
            end
        end

        function onParentResized(obj, prevFcn, src, evt)
            if ~isempty(prevFcn)
                try, feval(prevFcn, src, evt); catch, end
            end
            obj.positionZoomPanel();
        end

        function reattachZoomUI(obj)
            % Call this after you move the axes to a new parent panel
            if isempty(obj.zoomPanel) || ~isvalid(obj.zoomPanel), return; end
            newParent = obj.main_axes.Parent;
            if obj.zoomPanel.Parent ~= newParent
                obj.zoomPanel.Parent = newParent;
            end
            % rehook parent panel resize callback
            try
                parentPrevSCF = newParent.SizeChangedFcn;
                newParent.SizeChangedFcn = @(src,evt) obj.onParentResized(parentPrevSCF, src, evt);
            catch
            end
            obj.positionZoomPanel();
            uistack(obj.zoomPanel,'top');
        end

        function draw_wfs_zones(obj, rp, renderer)
            % Optimized zone drawer:
            % - Preallocates/keeps a pool of patch handles for zone fills (BG layer)
            % - Uses ONE line object for all dashed boundary rays (FG layer)
            % - No per-frame deletes/uistack calls

            ax = obj.main_axes;
            L  = obj.ensureLayers();

            % ---- one-time setup (or reuse) ----
            if ~isfield(obj.renderer_illustration,'zonesInit') || isempty(obj.renderer_illustration.zonesInit)
                % Group for fills (kept at bottom)
                obj.renderer_illustration.zonesGroup = hggroup('Parent',L.bg, ...
                    'HitTest','off','PickableParts','none','Tag','ZonesGroup');

                % Preallocate a small pool; will grow if needed
                obj.renderer_illustration.zoneFills = gobjects(0);
                obj.renderer_illustration.zoneFillsN = 0;

                % One dashed line to hold ALL boundary rays
                obj.renderer_illustration.zoneRays = line('Parent',L.fg, ...
                    'XData',NaN,'YData',NaN,'LineStyle','--','Color',[0 0 0], ...
                    'LineWidth',0.75,'HitTest','off','Tag','ZonesRays');

                obj.renderer_illustration.zonesInit = true;
            end

            % ---- normalize zones ----
            zones = normalizeZones(rp);
            if isempty(zones)
                hideAllFills();
                set(obj.renderer_illustration.zoneRays,'XData',NaN,'YData',NaN);
                return;
            end

            % ---- ensure we have enough preallocated fills ----
            need = numel(zones);
            have = obj.renderer_illustration.zoneFillsN;
            if need > have
                growBy = max(3, need - have);
                newF = gobjects(growBy,1);
                for i = 1:growBy
                    newF(i) = patch('Parent', obj.renderer_illustration.zonesGroup, ...
                        'XData',NaN,'YData',NaN,'CData',1, ...
                        'FaceColor',[0.9 0.9 0.9],'FaceAlpha',0.50, ...
                        'EdgeColor','none','HitTest','off', ...
                        'Tag', sprintf('zoneFill_%d', have+i));
                end
                if have == 0
                    obj.renderer_illustration.zoneFills = newF;
                else
                    obj.renderer_illustration.zoneFills = [obj.renderer_illustration.zoneFills; newF];
                end
                obj.renderer_illustration.zoneFillsN = have + growBy;
            end

            % ---- geometry constants ----
            xl = xlim(ax); yl = ylim(ax);
            rayLen  = 50 * max(diff(xl), diff(yl));   % long-enough rays

            srcType = lower(renderer.virtual_source.source_type.Shape);  % 'point_source'|'plane_wave'
            F       = rp.F;

            % ---- build data for fills and rays (batched) ----
            raysX = []; raysY = [];

            for k = 1:numel(zones)
                z = zones(k);
                if isempty(z.boundaries) || numel(z.boundaries)~=2, continue; end
                i1 = z.boundaries(1); i2 = z.boundaries(2);
                if any(~isfinite([i1 i2])) || i1<1 || i2<1 || i1>size(rp.x0,1) || i2>size(rp.x0,1)
                    continue;
                end

                seg = betweenIndices(i1, i2, size(rp.x0,1), isfield(rp,'isClosed') && rp.isClosed);
                if isempty(seg), seg = i1; end

                % boundary rays (two per zone)
                [xb1,yb1] = boundaryRay(i1, srcType, F, rayLen, rp, renderer);
                [xb2,yb2] = boundaryRay(i2, srcType, F, rayLen, rp, renderer);

                raysX = [raysX, xb1, NaN, xb2, NaN]; %#ok<AGROW>
                raysY = [raysY, yb1, NaN, yb2, NaN]; %#ok<AGROW>

                % zone polygon
                [Xp,Yp] = zonePolygon(seg, i1, i2, srcType, F, rayLen, rp, renderer);

                % update k-th fill
                h = obj.renderer_illustration.zoneFills(k);
                set(h, 'XData', Xp, 'YData', Yp, ...
                    'FaceColor', getColor(z.type), 'FaceAlpha', 0.50, 'Visible','on');
            end

            % hide unused fills in the pool
            for k = (numel(zones)+1) : obj.renderer_illustration.zoneFillsN
                set(obj.renderer_illustration.zoneFills(k),'Visible','off');
            end

            % one shot update for all rays
            set(obj.renderer_illustration.zoneRays,'XData',raysX,'YData',raysY);

            % =================== nested helpers ===================
            function hideAllFills()
                for kk = 1:obj.renderer_illustration.zoneFillsN
                    if isgraphics(obj.renderer_illustration.zoneFills(kk))
                        set(obj.renderer_illustration.zoneFills(kk),'Visible','off');
                    end
                end
            end
            function zonesOut = normalizeZones(rpIn)
                if isfield(rpIn,'zone') && ~isempty(rpIn.zone), zonesOut = rpIn.zone; return; end
                if isfield(rpIn,'zones') && numel(rpIn.zones)==4
                    z = rpIn.zones(:).';
                    zonesOut = struct( ...
                        'boundaries',{[z(1) z(2)] , [z(2) z(3)] , [z(3) z(4)]}, ...
                        'type',     {'tapered'     , 'illuminated' , 'tapered'   });
                else
                    zonesOut = [];
                end
            end
            function idx = betweenIndices(iStart,iEnd,N,isClosed)
                if isClosed
                    if iEnd >= iStart, idx = iStart:iEnd; else, idx = [iStart:N, 1:iEnd]; end
                else
                    if iStart > iEnd, [iStart,iEnd] = deal(iEnd,iStart); end
                    idx = iStart:iEnd;
                end
                idx = unique(idx(:).','stable');
            end
            function [x,y] = boundaryRay(ii, srcType, Fval, len, rpIn, rend)
                dir = rpIn.kh(ii,:);   % unit
                switch srcType
                    case 'plane_wave'
                        x0i = rpIn.x0(ii,1); y0i = rpIn.x0(ii,2);
                        x = x0i + [-len len]*dir(1);
                        y = y0i + [-len len]*dir(2);
                    case 'point_source'
                        vs = rend.virtual_source.position;
                        if Fval == 1
                            x = [vs(1), vs(1) + len*dir(1)];
                            y = [vs(2), vs(2) + len*dir(2)];
                        else
                            x0i = rpIn.x0(ii,1); y0i = rpIn.x0(ii,2);
                            x = [x0i, x0i - len*dir(1)];
                            y = [y0i, y0i - len*dir(2)];
                        end
                end
            end
            function [Xp,Yp] = zonePolygon(segIdx,i1,i2,srcType,Fval,len,rpIn,rend)
                dir1 = rpIn.kh(i1,:); dir2 = rpIn.kh(i2,:);
                switch srcType
                    case 'plane_wave'
                        x0seg = rpIn.x0(segIdx,1); y0seg = rpIn.x0(segIdx,2);
                        p1 = rpIn.x0(i1,:) + len*dir1;
                        p2 = rpIn.x0(i2,:) + len*dir2;
                        Xp = [x0seg; p2(1); p1(1)];
                        Yp = [y0seg; p2(2); p1(2)];
                    case 'point_source'
                        vs = rend.virtual_source.position(:).';
                        if Fval == 1
                            p1 = vs + len*dir1;
                            p2 = vs + len*dir2;
                            x0seg = rpIn.x0(segIdx,1); y0seg = rpIn.x0(segIdx,2);
                            Xp = [x0seg; p2(1); p1(1)];
                            Yp = [y0seg; p2(2); p1(2)];
                        else
                            p1 = rpIn.x0(i1,:) - len*dir1;
                            p2 = rpIn.x0(i2,:) - len*dir2;
                            Xp = [vs(1); p2(1); p1(1)];
                            Yp = [vs(2); p2(2); p1(2)];
                        end
                end
                Xp = Xp(:); Yp = Yp(:);
            end
            function c = getColor(kind)
                switch lower(kind)
                    case 'illuminated', c = [255 248 174]/255;
                    case 'tapered',     c = [236 227 128]/255;
                    otherwise,          c = [0.9 0.9 0.9];
                end
            end
        end

        function L = ensureLayers(obj)
            % Create (or reuse) 3 groups with stable stacking:
            % bg (zones) < mid (amp/taper) < fg (lines/labels)
            if ~isfield(obj.renderer_illustration,'layers') || ...
               ~all(isfield(obj.renderer_illustration.layers,{'bg','mid','fg'})) || ...
               ~all(isgraphics(struct2array(obj.renderer_illustration.layers)))

                ax = obj.main_axes;
                L.bg  = hggroup('Parent',ax,'Tag','layer-bg','HitTest','off','PickableParts','none');
                L.mid = hggroup('Parent',ax,'Tag','layer-mid','HitTest','off','PickableParts','none');
                L.fg  = hggroup('Parent',ax,'Tag','layer-fg','HitTest','off','PickableParts','none');
                try
                    uistack(L.bg ,'bottom');  % run ONCE
                    uistack(L.mid,'top');
                    uistack(L.fg ,'top');
                catch
                end
                obj.renderer_illustration.layers = L;
                set(ax,'SortMethod','childorder');
                set(ancestor(ax,'figure'),'Renderer','opengl','GraphicsSmoothing','off');
            else
                L = obj.renderer_illustration.layers;
            end
        end

    end
end
