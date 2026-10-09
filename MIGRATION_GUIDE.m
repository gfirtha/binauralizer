%% ═══════════════════════════════════════════════════════════════════
%%  MIGRATION GUIDE — RendererRegistry bevezetése
%%
%%  ELŐTTE: 13 helyen kellett módosítani egy új renderer bevezetéséhez
%%  UTÁNA:  3 fájl (renderer + panel + 1 blokk a RendererRegistry-ben)
%% ═══════════════════════════════════════════════════════════════════

%% ─────────────────────────────────────────────────────────────────
%%  1.  sound_scene_renderer.make_renderer_for()
%%      ELŐTTE: 15-ágú switch-case
%%      UTÁNA:
%% ─────────────────────────────────────────────────────────────────

% function r = make_renderer_for(obj, vs, SSD, loudspeaker_array, setup)
%     rtype = obj.get_renderer_type(vs);
%     reg = RendererRegistry.instance();
%     r = reg.create(rtype, vs, SSD, loudspeaker_array, obj.receiver, setup);
% end

%% ─────────────────────────────────────────────────────────────────
%%  2.  AppContext.rendererDictionary()
%%      ELŐTTE: switch-case 10 ággal
%%      UTÁNA:  TÖRÖLHETŐ — helyette:
%% ─────────────────────────────────────────────────────────────────

% function modeOut = rendererDictionary(~, modeIn)
%     modeOut = RendererRegistry.instance().resolveId(modeIn);
% end

%% ─────────────────────────────────────────────────────────────────
%%  3.  ModeMap.m
%%      ELŐTTE: saját registry() + toInternal() + toUI() + uiList()
%%      UTÁNA:  TELJES FÁJL TÖRÖLHETŐ — RendererRegistry átveszi:
%%              - ModeMap.toInternal(x)  →  RendererRegistry.instance().resolveId(x)
%%              - ModeMap.toUI(x)        →  RendererRegistry.instance().toUI(x)
%%              - ModeMap.uiList()       →  RendererRegistry.instance().uiList()
%% ─────────────────────────────────────────────────────────────────

%% ─────────────────────────────────────────────────────────────────
%%  4–6.  Rules.isutilible / preferredGeometry / preferredRenderer
%%        ELŐTTE: 3 külön switch-case
%%        UTÁNA:  delegálás:
%% ─────────────────────────────────────────────────────────────────

% function tf = isutilible(renderer, geom)
%     tf = RendererRegistry.instance().supportsGeometry(renderer, geom);
% end
%
% function geom = preferredGeometryForRenderer(renderer, ~)
%     geom = RendererRegistry.instance().preferredGeometry(renderer);
% end
%
% function renderer = preferredRendererForGeometry(geom, ~)
%     renderer = RendererRegistry.instance().preferredRendererForGeometry(geom);
% end

%% ─────────────────────────────────────────────────────────────────
%%  7.  Rules.canonR()
%%      ELŐTTE: switch-case 8 ággal
%%      UTÁNA:  TÖRÖLHETŐ — toUI() pótolja:
%% ─────────────────────────────────────────────────────────────────

% function r = canonR(r)
%     r = RendererRegistry.instance().toUI(r);
% end

%% ─────────────────────────────────────────────────────────────────
%%  8.  RendererPanelFactory.create()
%%      ELŐTTE: switch-case 9 ággal
%%      UTÁNA:
%% ─────────────────────────────────────────────────────────────────

% classdef RendererPanelFactory
%     methods (Static)
%         function pnl = create(modeUI, parent, C)
%             pnl = RendererRegistry.instance().createPanel(modeUI, parent, C);
%         end
%     end
% end

%% ─────────────────────────────────────────────────────────────────
%%  9.  RendererTab constructor — popup lista
%%      ELŐTTE: modes = {'WFS','VBAP','DBAP',...} (kézzel)
%%      UTÁNA:
%% ─────────────────────────────────────────────────────────────────

% modes = RendererRegistry.instance().uiList();

%% ─────────────────────────────────────────────────────────────────
%%  10. RendererTab.refreshFromSetup() — class-name inference
%%      ELŐTTE: 15 soros if-elseif lánc (contains(cls,'wfs')...)
%%      UTÁNA:  RendererPanel kap egy 'Mode' property-t:
%% ─────────────────────────────────────────────────────────────────

% % A RendererDescriptor.createPanel() beállítja:
% %   pnl.Mode = descriptor.id;
% %
% % Ezután a refreshFromSetup()-ban:
% curId = RendererRegistry.instance().resolveId(obj.C.Setup.Rendering_mode);
% if isempty(obj.CurrentPanel) || ~isvalid(obj.CurrentPanel) ...
%         || ~strcmp(obj.CurrentPanel.Mode, curId)
%     obj.rebuildPanel();
% else
%     obj.CurrentPanel.syncFromSetup();
% end

%% ─────────────────────────────────────────────────────────────────
%%  11. AppContext.defaultSetup() — Renderer_setup default értékek
%%      ELŐTTE: egyetlen monolitikus struct literál
%%      UTÁNA:  a közös mezők maradnak, a renderer-specifikusak
%%              a registry-ből jönnek:
%% ─────────────────────────────────────────────────────────────────

% function S = defaultSetup(obj)
%     ...  % közös mezők (SampleRate, Block_size, HRTF, stb.)
%     reg = RendererRegistry.instance();
%     rendererDefaults = reg.defaultSetup(S.Rendering_mode);
%     fns = fieldnames(rendererDefaults);
%     for k = 1:numel(fns)
%         S.Renderer_setup.(fns{k}) = rendererDefaults.(fns{k});
%     end
% end

%% ─────────────────────────────────────────────────────────────────
%%  12–13. Új renderer + panel fájl létrehozása
%%         Ez természetesen megmarad — de CSAK ennyi a teendő.
%% ─────────────────────────────────────────────────────────────────


%% ═══════════════════════════════════════════════════════════════════
%%  ÖSSZEFOGLALÓ
%% ═══════════════════════════════════════════════════════════════════
%%
%%  ┌───────────────────────────────────────────────────────────────┐
%%  │  ÚJ RENDERER HOZZÁADÁSA — TELJES CHECKLIST                  │
%%  │                                                               │
%%  │  1. Írd meg:   my_renderer.m  (< base_renderer)             │
%%  │  2. Írd meg:   MyPanel.m      (< RendererPanel)             │
%%  │  3. Regisztráld a RendererRegistry.registerAll()-ban:        │
%%  │                                                               │
%%  │     obj.register(RendererDescriptor( ...                     │
%%  │         'id',           'my_renderer', ...                   │
%%  │         'uiName',       'My Renderer', ...                   │
%%  │         'synonyms',     {{'my-renderer'}}, ...               │
%%  │         'factory',      @(vs,SSD,ls,rec,s) ..., ...          │
%%  │         'panelClass',   'audioapp.ui.renderers.MyPanel', ... │
%%  │         'geometries',   {{'linear','circular'}}, ...         │
%%  │         'prefGeometry', 'circular', ...                      │
%%  │         'defaultSetup', struct('Param1', 42) ));             │
%%  │                                                               │
%%  │  KÉSZ. Nincs más teendő.                                     │
%%  └───────────────────────────────────────────────────────────────┘
%%
%%  Ami ELTŰNIK (törölhető vagy 1-soros wrapper marad):
%%    - AppContext.rendererDictionary()    → reg.resolveId()
%%    - ModeMap.m (teljes fájl)           → reg.resolveId/toUI/uiList
%%    - Rules.isutilible()                → reg.supportsGeometry()
%%    - Rules.preferredGeometry*()        → reg.preferredGeometry()
%%    - Rules.preferredRenderer*()        → reg.preferredRendererForGeometry()
%%    - Rules.canonR()                    → reg.toUI()
%%    - RendererPanelFactory.create()     → reg.createPanel()
%%    - RendererTab hardcoded modes list  → reg.uiList()
%%    - RendererTab class-name inference  → panel.Mode property
%%    - make_renderer_for() switch-case   → reg.create()
%%    - defaultSetup() renderer fields    → reg.defaultSetup()
