classdef SourceTable < handle
    properties
        C
        Table matlab.ui.control.Table
        ColRatios = [0.12 0.17 0.12 0.10 0.035 0.035]   % Source | Model | Renderer | Volume | S | M
        SetupLsnr event.listener
        Vol  = []           % user volume per source (the virtual source gets Vol x solo/mute)
        Solo = logical([])
        Mute = logical([])
    end

    methods
        function obj = SourceTable(parent, normPos, ctx)
            obj.C = ctx;
            obj.Table = uitable(parent,'Units','normalized','Position',normPos, ...
                'ColumnName',{'Source','Model','Renderer','Volume','S','M'},...
                'ColumnEditable',[true true false true true true],...
                'CellEditCallback',@(s,e)obj.onEdit(e));
            obj.refresh();
            fig = ancestor(parent,'figure');
            fig.SizeChangedFcn = @(~,~) obj.resize();
            obj.resize();

            obj.SetupLsnr = addlistener(obj.C, 'SetupChanged', @(~,~)obj.refresh());
        end

        function resize(obj)
            pad = 24;
            pos = getpixelposition(obj.Table,true);
            W   = max(120, pos(3)-pad);
            r   = obj.ColRatios(:)'/sum(obj.ColRatios);
            w   = max(60, floor(W * r));
            try
                obj.Table.ColumnWidth = w;            % uifigure path
            catch
                obj.Table.ColumnWidth = num2cell(w);  % classic figure path
            end
        end

        function refresh(obj)
            N = numel(obj.C.Scene.virtual_sources);
            if numel(obj.Vol) ~= N                    % new scene size: fresh solo / mute state
                obj.Vol = arrayfun(@(k) obj.getVSVolume(k), 1:N);
                obj.Solo = false(1, N); obj.Mute = false(1, N);
            end
            data = cell(N,6);
            for k = 1:N
                % Source name (as before)
                if isvalid(obj.C.SceneGUI.virtual_source_points{k})
                    data{k,1} = sprintf('%s', obj.C.SceneGUI.virtual_source_points{k}.UserData.text.String);
                else
                    data{k,1} = sprintf('s_%d', k);
                end
                % Model: read from VS (normalized)
                data{k,2} = obj.normalizeModel( obj.getVSModel(k) );
                % Renderer: current app-wide renderer (keep as your placeholder or per-VS if you add it)
                data{k,3} = obj.C.Setup.Rendering_mode;
                % Volume (the user's setting; the source itself gets Vol x solo/mute)
                data{k,4} = obj.Vol(k);
                data{k,5} = obj.Solo(k);
                data{k,6} = obj.Mute(k);
            end
            set(obj.Table,'Data',data);

            % Column formats: Source (free text), Model (popup), Renderer (popup), Volume (numeric), S, M (checkboxes)
            set(obj.Table,'ColumnFormat',{[], {'point_source','plane_wave'}, obj.allowedRenderers(), 'numeric', 'logical', 'logical'});
            obj.applyGains();      % a rebuilt scene has fresh sources: re-apply solo / mute
        end

        function onEdit(obj, evt)
            r = evt.Indices(1); c = evt.Indices(2);

            switch c
                case 1 % Source name (label)
                    newUI = evt.NewData;
                    obj.C.SceneGUI.virtual_source_points{r}.Tag = newUI;
                    obj.C.SceneGUI.virtual_source_points{r}.UserData.text.String = newUI;

                case 2 % Model (point_source / plane_wave)
                    newVal = obj.normalizeModel(evt.NewData);

                    % set on the virtual source
                    try
                        obj.C.Scene.virtual_sources{r}.source_type.Shape = newVal;
                    catch ME
                        warning('Failed to set source model on VS %d: %s', r, ME.message);
                        obj.Table.Data{r,c} = obj.normalizeModel(evt.PreviousData);
                        return;
                    end
                    % rebuild scene and refresh
                    obj.C.Scene.scene_renderer.SFS_renderer{r}.update_renderer; %% !!!!
                    obj.C.SceneGUI.request_draw_visuals(obj.C.Scene.scene_renderer.SFS_renderer{r});
%                    obj.C.Scene.scene_renderer.SFS_renderer{r}.update_renderer_settings;%%!!
                    obj.refresh();

                case 3 % Renderer change (if/when you implement per-VS renderer)
                    % (left as your placeholder; currently Renderer column is not editable)
                    % You can wire this to set a per-VS renderer_type and call obj.C.rebuildSoundScene().
                    return;

                case 4 % Volume
                    val = max(0,min(1,double(evt.NewData)));
                    obj.Vol(r) = val;
                    obj.applyGains();
                    % keep cell display consistent (clamped)
                    obj.Table.Data{r,c} = val;

                case 5 % Solo
                    obj.Solo(r) = logical(evt.NewData);
                    obj.applyGains();

                case 6 % Mute
                    obj.Mute(r) = logical(evt.NewData);
                    obj.applyGains();
            end
        end

        function applyGains(obj)
            % gain = volume, 0 when muted, 0 when other sources are soloed
            N = min(numel(obj.Vol), numel(obj.C.Scene.virtual_sources));
            anySolo = any(obj.Solo(1:N));
            for k = 1:N
                on = ~obj.Mute(k) && (~anySolo || obj.Solo(k));
                obj.setVSVolume(k, obj.Vol(k)*on);
            end
        end

        function g = getVSVolume(obj,idx)
            g = 1;
            try
                if ismethod(obj.C.Scene.virtual_sources{idx},'get_gain')
                    g = obj.C.Scene.virtual_sources{idx}.get_gain();
                end
            catch
            end
            if isempty(g) || ~isfinite(g), g = 1; end
        end

        function setVSVolume(obj,idx, val)
            try
                if ismethod(obj.C.Scene.virtual_sources{idx},'set_gain')
                    obj.C.Scene.virtual_sources{idx}.set_gain(val);
                end
            catch
                % no-op if the VS doesn’t expose an amp interface
            end
        end

        function list = allowedRenderers(obj)
            % keep your geometry-based mapping here or query helpers
            switch lower(obj.C.Setup.Loudspeaker_setup.Shape)
                case 'linear',   list = {'WFS','VBAP','Direct playback'};
                case 'circular', list = {'WFS','HOA','NFC_HOA','VBAP','Direct playback'};
                case 'stereo',   list = {'VBAP','Direct playback'};
                case '5.0',      list = {'Dolby Surround Decoder','HOA','NFC_HOA','VBAP','Direct playback'};
                otherwise,       list = {'VBAP','Direct playback'};
            end
        end
    end

    methods (Access = private)
        function s = getVSModel(obj, idx)
            s = 'point_source'; % default
            try
                s = obj.C.Scene.virtual_sources{idx}.source_type.Shape;
            catch
            end
            s = obj.normalizeModel(s);
        end

        function setVSModel(obj, idx, newVal)
            obj.C.Scene.virtual_sources{idx}.source_type.Shape = newVal;
            obj.C.rebuildSoundScene;
        end

        function s = normalizeModel(~, s)
            % Normalize to {'point_source','plane_wave'}
            s = char(string(s));
            s = lower(strtrim(s));
            s = strrep(s,'-','_');
            s = strrep(s,' ','_');
            if startsWith(s,'point'), s = 'point_source'; end
            if startsWith(s,'plane'), s = 'plane_wave'; end
            if ~ismember(s, {'point_source','plane_wave'})
                s = 'point_source';
            end
        end
    end
end
