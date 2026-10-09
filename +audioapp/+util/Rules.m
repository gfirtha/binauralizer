classdef Rules
    % Centralized renderer <-> geometry compatibility & defaults.
    %
    % Renderer-specific logic (which geometries are allowed, preferred
    % geometry, etc.) now lives in RendererRegistry descriptors.
    % This class delegates to the registry and keeps only:
    %   - geometry normalization & default N
    %   - the incompatibility warning dialog
    %   - legacy wrappers so existing callers keep working

    methods (Static)

        % ─── Delegates to RendererRegistry ───

        function tf = isutilible(renderer, geom)
            tf = audioapp.RendererRegistry.instance().supportsGeometry(renderer, geom);
        end

        function geom = preferredGeometryForRenderer(renderer, ~)
            geom = audioapp.RendererRegistry.instance().preferredGeometry(renderer);
        end

        function renderer = preferredRendererForGeometry(geom, ~)
            renderer = audioapp.RendererRegistry.instance().preferredRendererForGeometry(geom);
        end

        function [modeOut, geomOut, changed] = resolveOnRendererChange(rendererIn, geomCur)
            reg = audioapp.RendererRegistry.instance();
            [modeOut, geomOut, changed] = reg.resolveOnRendererChange(rendererIn, geomCur);
            if changed
                audioapp.util.Rules.showIncompatWarn(rendererIn, geomCur, ...
                    sprintf('Switched SSD geometry to "%s" to match renderer "%s".', geomOut, modeOut));
            end
        end

        function [modeOut, geomOut, changed] = resolveOnGeometryChange(geomIn, modeCur)
            reg = audioapp.RendererRegistry.instance();
            [modeOut, geomOut, changed] = reg.resolveOnGeometryChange(geomIn, modeCur);
            if changed
                audioapp.util.Rules.showIncompatWarn(modeCur, geomIn, ...
                    sprintf('Switched renderer to "%s" to match geometry "%s".', modeOut, geomOut));
            end
        end

        % ─── Naming helpers (legacy wrappers) ───

        function r = normR(r)
            r = lower(strrep(strrep(strtrim(r),'-','_'),' ','_'));
        end

        function r = canonR(r)
            % Returns the UI-friendly name via registry
            r = audioapp.RendererRegistry.instance().toUI(r);
        end

        function g = normG(g)
            g = lower(strtrim(g));
            switch g
                case {'5.0 (circular)','5_0','5.0'}, g = '5.0';
                case {'atmos 5.1.4'}, g = '5.1.4';
                case {'atmos 7.1.4'}, g = '7.1.4';
                case {'soundbar (16 ch)'}, g = 'soundbar';
                otherwise, g = g;
            end
        end

        % ─── Geometry defaults (not renderer-specific, stays here) ───

        function N = defaultNForGeometry(geom)
            switch audioapp.util.Rules.normG(geom)
                case 'stereo', N = 2;
                case '5.0',    N = 5;
                case '5.1.4',  N = 9;     % + LFE channel without a loudspeaker
                case '7.1.4',  N = 11;
                case 'soundbar', N = 16;  % + 4 when the discrete surrounds are on
                case 'linear', N = 64;
                case 'circular', N = 64;
                case 'general', N = 96;
                case 'genius_array', N = 8;
                otherwise,     N = 2;
            end
        end

        % ─── Warning dialog ───

        function showIncompatWarn(rendererUI, geomUI, varargin)
            if nargin >= 3 && ~isempty(varargin{1})
                autoAdjustNote = sprintf('\n\nAdjustment applied:\n%s', varargin{1});
            else
                autoAdjustNote = '';
            end

            reg = audioapp.RendererRegistry.instance();
            % Build allowed-combinations list dynamically from registry
            uiNames = reg.uiList();
            lines = {};
            for k = 1:numel(uiNames)
                d = reg.get(uiNames{k});
                if any(strcmp(d.geometries,'*'))
                    geomStr = 'any';
                else
                    geomStr = strjoin(d.geometries, ', ');
                end
                lines{end+1} = sprintf('  %s  ->  %s', d.uiName, geomStr); %#ok<AGROW>
            end

            msg = sprintf(['The selected renderer "%s" is incompatible with the current geometry "%s".\n\n' ...
                'Allowed combinations:\n%s\n'], ...
                rendererUI, geomUI, strjoin(lines, '\n'));
            warndlg([msg autoAdjustNote], 'Incompatible selection', 'modal');
        end
    end
end
