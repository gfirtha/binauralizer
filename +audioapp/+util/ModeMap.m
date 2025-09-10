classdef ModeMap
    %MODEMAP Central mapping between UI labels and internal renderer IDs.
    %
    % Public API (all static):
    %   uiList()              -> cellstr of UI labels (order for popup)
    %   toInternal(name)      -> normalized internal id (e.g., 'nfc_hoa')
    %   toUI(nameOrId)        -> canonical UI label from id or any synonym
    %   uiIndex(uiList, name) -> index in a given UI list for a name/id

    methods (Static)
        function list = uiList()
            items = audioapp.util.ModeMap.registry();
            list  = {items.ui};
        end

        function id = toInternal(name)
            if isempty(name), id = 'wfs'; return; end
            key = audioapp.util.ModeMap.norm(name);
            items = audioapp.util.ModeMap.registry();
            % Try exact on IDs
            for k = 1:numel(items)
                if strcmp(key, audioapp.util.ModeMap.norm(items(k).id))
                    id = items(k).id; return;
                end
            end
            % Try UI
            for k = 1:numel(items)
                if strcmp(key, audioapp.util.ModeMap.norm(items(k).ui))
                    id = items(k).id; return;
                end
            end
            % Try synonyms
            for k = 1:numel(items)
                syn = items(k).syn;
                for j = 1:numel(syn)
                    if strcmp(key, audioapp.util.ModeMap.norm(syn{j}))
                        id = items(k).id; return;
                    end
                end
            end
            % Fallback
            id = 'wfs';
        end

        function ui = toUI(nameOrId)
            if isempty(nameOrId), ui = 'WFS'; return; end
            id = audioapp.util.ModeMap.toInternal(nameOrId);
            items = audioapp.util.ModeMap.registry();
            for k = 1:numel(items)
                if strcmp(id, items(k).id)
                    ui = items(k).ui; return;
                end
            end
            ui = 'WFS';
        end

        function idx = uiIndex(uiListCell, nameOrId)
            if isempty(uiListCell), idx = 1; return; end
            targetUI = audioapp.util.ModeMap.toUI(nameOrId);
            idx = find(strcmp(uiListCell, targetUI), 1);
            if isempty(idx), idx = 1; end
        end
    end

    methods (Static, Access = private)
        function items = registry()
            % One place to change labels/IDs/synonyms.
            %
            % ui  : label shown in popup
            % id  : internal, passed around the engine
            % syn : accepted inputs that normalize to 'id'
            items = struct( ...
                'ui', {}, 'id', {}, 'syn', {} );

            items(end+1) = struct( ...
                'ui','WFS', ...
                'id','wfs', ...
                'syn',{{'wave field synthesis','wave-field-synthesis'}});

            items(end+1) = struct( ...
                'ui','VBAP', ...
                'id','vbap', ...
                'syn',{{}});

            items(end+1) = struct( ...
                'ui','HOA', ...
                'id','hoa', ...
                'syn',{{'higher order ambisonics','ambisonics'}});

            items(end+1) = struct( ...
                'ui','NFC_HOA', ...
                'id','nfc_hoa', ...
                'syn',{{'nfc-hoa','near-field-compensated hoa','near field compensated hoa'}});

            % UI says "Encoder", your internal class says "Decoder".
            % We accept both and normalize to 'dolby_surround_decoder'.
            items(end+1) = struct( ...
                'ui','Dolby Surround Encoder', ...
                'id','dolby_surround_decoder', ...
                'syn',{{'dolby','dolby surround','dolby-surround', ...
                        'dolby surround decoder','dolby-surround-decoder', ...
                        'dolby surround encoder','dolby-surround-encoder'}});

            items(end+1) = struct( ...
                'ui','Direct playback', ...
                'id','direct_playback', ...
                'syn',{{'direct','direct-playback'}});

            % If you add new modes later, only touch here.
        end

        function s = norm(s0)
            if isstring(s0), s0 = char(s0); end
            if ~ischar(s0), s0 = ''; end
            s = lower(s0);
            s = strrep(s,'-','_');
            s = strrep(s,' ','_');
            s = regexprep(s,'__+','_');
            s = strtrim(s);
        end
    end
end
