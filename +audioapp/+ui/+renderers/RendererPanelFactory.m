classdef RendererPanelFactory
    methods (Static)
        function pnl = create(modeUI, parent, C)
            switch upper(modeUI)
                case 'WFS'
                    pnl = audioapp.ui.renderers.WFSPanel(C);    pnl.build(parent);
                case 'HOA'
                    pnl = audioapp.ui.renderers.HOAPanel(C);    pnl.build(parent);
                case 'NFC_HOA'
                    pnl = audioapp.ui.renderers.NFCHOAPanel(C); pnl.build(parent);
                case 'VBAP'
                    pnl = audioapp.ui.renderers.VBAPPanel(C);   pnl.build(parent);
                case 'DBAP'
                    pnl = audioapp.ui.renderers.DBAPPanel(C);   pnl.build(parent);
                case 'DOLBY SURROUND ENCODER'
                    pnl = audioapp.ui.renderers.DolbyPanel(C);  pnl.build(parent);
                case 'CTC'
                    pnl = audioapp.ui.renderers.CTCPanel(C); pnl.build(parent);
                case 'TD_STEREO'
                    pnl = audioapp.ui.renderers.TDPanel(C); pnl.build(parent);
                otherwise % Direct playback, etc.
                    pnl = audioapp.ui.renderers.DirectPanel(C); pnl.build(parent);
            end
        end
    end
end
