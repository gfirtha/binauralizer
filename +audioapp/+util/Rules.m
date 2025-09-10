classdef Rules
    % Centralized renderer <-> geometry utilibility & defaults
    methods (Static)
        function tf = isutilible(renderer, geom)
            r = audioapp.util.Rules.normR(renderer);
            g = audioapp.util.Rules.normG(geom);
            switch r
                case {'wfs', 'vbap'}
                    tf = any(strcmp(g, {'linear','circular','5.0','general'}));
                case {'hoa','nfc_hoa'}
                    tf = any(strcmp(g, {'circular','5.0'}));
                case 'dolby_surround_encoder'
                    tf = strcmp(g, 'stereo'); % 2.0 layout
                case {'td_stereo', 'ctc'}
                    tf = any(strcmp(g, {'circular','5.0','stereo'}));
                otherwise % vbap, direct_playback, etc.
                    tf = true;
            end
        end

        function geom = preferredGeometryForRenderer(renderer, currentGeom)
            r = audioapp.util.Rules.normR(renderer);
            candidates = {'linear','circular','stereo','5.0'};
            switch r
                case 'wfs',              pref = {'linear','circular','5.0'};
                case {'hoa','nfc_hoa'},  pref = {'circular','5.0'};
                case 'dolby_surround_encoder', pref = {'stereo'};
                case {'ctc','td_stereo'},  pref = {'5.0','stereo'};
                otherwise,               pref = {audioapp.util.Rules.normG(currentGeom), candidates{:}};
            end
            geom = pref{1};
        end

        function renderer = preferredRendererForGeometry(geom, currentRenderer)
            g = audioapp.util.Rules.normG(geom);
            switch g
                case 'linear',  pref = {'wfs','vbap','direct_playback'};
                case 'circular',pref = {'wfs','hoa','nfc_hoa','vbap','direct_playback'};
                case 'stereo',  pref = {'vbap','direct_playback','dolby_surround_encoder'};
                case '5.0',     pref = {'hoa','nfc_hoa','vbap','direct_playback'};
                otherwise,      pref = {'vbap','direct_playback'};
            end
            renderer = pref{1};
        end

        function N = defaultNForGeometry(geom)
            switch audioapp.util.Rules.normG(geom)
                case 'stereo', N = 2;
                case '5.0',    N = 5;
                case 'linear', N = 64;
                case 'circular', N = 64;
                case 'general', N = 96;
                otherwise,     N = 2;
            end
        end

        % Resolve when the USER changed the renderer
        function [modeOut, geomOut, changed] = resolveOnRendererChange(rendererIn, geomCur)
            modeOut = audioapp.util.Rules.canonR(rendererIn);
            geomOut = audioapp.util.Rules.normG(geomCur);
            if ~audioapp.util.Rules.isutilible(modeOut, geomOut)
                geomOut = audioapp.util.Rules.preferredGeometryForRenderer(modeOut, geomOut);
                changed = true;
                audioapp.util.Rules.showIncompatWarn(modeOut, geomCur,...
                    sprintf('Switched SSD geometry to "%s" to match renderer "%s".', geomOut, modeOut) );
            else
                changed = false;
            end
        end

        % Resolve when the USER changed the geometry
        function [modeOut, geomOut, changed] = resolveOnGeometryChange(geomIn, modeCur)
            geomOut = audioapp.util.Rules.normG(geomIn);
            modeOut = audioapp.util.Rules.normR(modeCur);
            if ~audioapp.util.Rules.isutilible(modeOut, geomOut)
                modeOut = audioapp.util.Rules.preferredRendererForGeometry(geomOut, modeCur);
                changed = true;
                audioapp.util.Rules.showIncompatWarn(modeCur, geomOut, ...
                    sprintf('Switched SSD geometry to "%s" to match renderer "%s".', geomOut, modeOut));

            else
                changed = false;
            end
            modeOut = audioapp.util.Rules.canonR(modeOut);
        end

        % ---------- naming helpers ----------
        function r = normR(r)
            r = lower(strrep(strrep(strtrim(r),'-','_'),' ','_'));
        end
        function r = canonR(r)
            r = audioapp.util.Rules.normR(r);
            switch r
                case 'dolby_surround_encoder'
                    r = 'Dolby Surround Encoder';
                case 'nfc_hoa'
                    r = 'NFC_HOA';
                case 'direct_playback'
                    r = 'Direct playback';
                case 'hoa'
                    r = 'HOA';
                case 'vbap'
                    r = 'VBAP';
                case 'wfs'
                    r = 'WFS';
                case 'TD stereo'
                    r = 'TD_stereo';
                otherwise
                    r = r;
            end
        end
        function g = normG(g)
            g = lower(strtrim(g));
            switch g
                case {'5.0 (circular)','5_0','5.0'}, g = '5.0';
                otherwise, g = g;
            end
        end

        function showIncompatWarn(rendererUI, geomUI, varargin)
            if nargin >= 3 && ~isempty(varargin{1})
                autoAdjustNote = sprintf('\n\nAdjustment applied:\n%s', varargin{1});
            else
                autoAdjustNote = '';
            end

            msg = sprintf(['The selected renderer "%s" is incompatible with the current geometry "%s".\n\n',...
                'Allowed combinations:\n',...
                '  • WFS                 → Linear or Circular or 5.0\n',...
                '  • HOA / NFC-HOA       → Circular or 5.0\n',...
                '  • Dolby Surround Enc. → 2.0 only\n',...
                '  • Direct playback     → Arbitrary\n'], ...
                rendererUI, geomUI);
            warndlg([msg autoAdjustNote], 'Incompatible selection', 'modal');
        end
    end
end
