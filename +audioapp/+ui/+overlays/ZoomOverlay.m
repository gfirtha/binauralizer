classdef ZoomOverlay < handle
    % Floating + / – buttons anchored to bottom-right of a given axes.

    properties
        ax
        panel
        btnMinus
        btnPlus
        cbZoomOut
        cbZoomIn
    end

    methods
        function obj = ZoomOverlay(ax, cbZoomOut, cbZoomIn)
            obj.ax = ax;
            obj.cbZoomOut = cbZoomOut;
            obj.cbZoomIn  = cbZoomIn;

            parentContainer = ax.Parent;
            bgc = ancestor(parentContainer,'figure').Color;

            obj.panel = uipanel( ...
                'Parent', parentContainer, ...
                'Units','pixels','BorderType','none', ...
                'BackgroundColor', bgc, ...
                'Tag','ZoomPanel');

            obj.btnMinus = uicontrol('Parent',obj.panel, ...
                'Style','pushbutton','String','–','Units','pixels', ...
                'Position',[0 0 28 28], 'TooltipString','Zoom out', ...
                'FontWeight','bold', 'Callback',@(~,~)obj.cbZoomOut(), ...
                'BackgroundColor',[0.94 0.94 0.94]);

            obj.btnPlus = uicontrol('Parent',obj.panel, ...
                'Style','pushbutton','String','+','Units','pixels', ...
                'Position',[28 0 28 28], 'TooltipString','Zoom in', ...
                'FontWeight','bold', 'Callback',@(~,~)obj.cbZoomIn(), ...
                'BackgroundColor',[0.94 0.94 0.94]);

            % React to resizes
            fig = ancestor(parentContainer,'figure');
            prevFigSCF = fig.SizeChangedFcn;
            fig.SizeChangedFcn = @(src,evt) obj.onFigureResized(prevFigSCF, src, evt);

            try
                prevParSCF = parentContainer.SizeChangedFcn;
                parentContainer.SizeChangedFcn = @(src,evt) obj.onParentResized(prevParSCF, src, evt);
            catch
            end

            try
                prevAxSCF = ax.SizeChangedFcn;
                ax.SizeChangedFcn = @(src,evt) obj.onAxesResized(prevAxSCF, src, evt);
            catch
            end

            obj.reposition();
        end

        function delete(obj)
            if ~isempty(obj.panel) && isvalid(obj.panel), delete(obj.panel); end
        end

        function attachTo(obj, ax)
            % Call when the axes is reparented
            obj.ax = ax;
            newParent = ax.Parent;
            if obj.panel.Parent ~= newParent
                obj.panel.Parent = newParent;
            end
            obj.reposition();
            uistack(obj.panel,'top');
        end

        function onFigureResized(obj, prevFcn, src, evt)
            if ~isempty(prevFcn), try, feval(prevFcn, src, evt); catch, end, end
            obj.reposition();
        end
        function onParentResized(obj, prevFcn, src, evt)
            if ~isempty(prevFcn), try, feval(prevFcn, src, evt); catch, end, end
            obj.reposition();
        end
        function onAxesResized(obj, prevFcn, src, evt)
            if ~isempty(prevFcn), try, feval(prevFcn, src, evt); catch, end, end
            obj.reposition();
        end

        function reposition(obj)
            if isempty(obj.panel) || ~isvalid(obj.panel) || isempty(obj.ax) || ~isvalid(obj.ax)
                return;
            end

            ax = obj.ax;
            par = obj.panel.Parent;

            % Make sure layout is up to date before querying TightInset
            drawnow limitrate nocallbacks

            % Axes position in PIXELS relative to its parent
            oldAxUnits = ax.Units; ax.Units = 'pixels';
            p  = ax.Position;              % [x y w h]
            ti = get(ax,'TightInset');     % [l b r t]
            ax.Units = oldAxUnits;

            % Be defensive about TightInset
            if ~isnumeric(ti) || numel(ti)~=4 || any(~isfinite(ti))
                ti = [0 0 0 0];
            end

            % "Plot box" = axes pos minus tight insets
            pbX = p(1) + ti(1);
            pbY = p(2) + ti(2);
            pbW = max(0, p(3) - ti(1) - ti(3));
            pbH = max(0, p(4) - ti(2) - ti(4));

            % Desired panel size and margin
            margin = 6; panelW = 56; panelH = 28;

            % Bottom-right inside the plot box
            x = pbX + pbW - panelW - margin;
            y = pbY + margin;

            % Clamp to parent bounds (also in PIXELS)
            oldParUnits = par.Units; par.Units = 'pixels';
            parPos = par.Position;          % [x y w h]
            par.Units = oldParUnits;

            x = min(max(x, 1), parPos(3) - panelW - 1);
            y = min(max(y, 1), parPos(4) - panelH - 1);

            obj.panel.Position = [x y panelW panelH];
            uistack(obj.panel,'top');
        end

        function attachTo(obj, ax)
            % Call when the axes is reparented/moved
            obj.ax = ax;
            newParent = ax.Parent;
            if obj.panel.Parent ~= newParent
                obj.panel.Parent = newParent;
            end

            % Re-hook resize callbacks to the new parent/axes/figure
            fig = ancestor(newParent,'figure');
            prevFigSCF = fig.SizeChangedFcn;
            fig.SizeChangedFcn = @(src,evt) obj.onFigureResized(prevFigSCF, src, evt);

            try
                prevParSCF = newParent.SizeChangedFcn;
                newParent.SizeChangedFcn = @(src,evt) obj.onParentResized(prevParSCF, src, evt);
            catch
            end

            try
                prevAxSCF = ax.SizeChangedFcn;
                ax.SizeChangedFcn = @(src,evt) obj.onAxesResized(prevAxSCF, src, evt);
            catch
            end

            obj.reposition();
            uistack(obj.panel,'top');
        end

    end
end
