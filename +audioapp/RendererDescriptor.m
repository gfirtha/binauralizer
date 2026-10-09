classdef RendererDescriptor
    % RendererDescriptor  Single source of truth for one renderer mode.
    %
    % Every renderer registers ONE descriptor.  The rest of the app
    % (Rules, ModeMap, RendererPanelFactory, sound_scene_renderer,
    %  RendererTab, AppContext) reads from here — ZERO switch-cases needed.
    %
    % === Új renderer bevezetése: EGYETLEN fájl, EGYETLEN struct ===
    %
    %   1.  Írd meg a renderer osztályt  (pl.  my_renderer < base_renderer)
    %   2.  Írd meg a UI panelt          (pl.  MyPanel < RendererPanel)
    %   3.  Regisztráld a RendererRegistry.m-ben — egyetlen blokk:
    %
    %         d = RendererDescriptor( ...
    %             'id',            'my_renderer', ...
    %             'uiName',        'My Renderer', ...
    %             'synonyms',      {{'my-renderer','myrenderer'}}, ...
    %             'factory',       @(vs,SSD,ls,rec,s) my_renderer(vs,SSD,s.Renderer_setup), ...
    %             'panelClass',    'audioapp.ui.renderers.MyPanel', ...
    %             'geometries',    {{'linear','circular','5.0'}}, ...
    %             'prefGeometry',  'circular', ...
    %             'defaultSetup',  struct('MyParam1', 42, 'MyParam2', true) );
    %
    %   Kész. Nincs más teendő.

    properties
        % ─── Identity ───
        id            (1,:) char       % internal key, e.g. 'wfs'
        uiName        (1,:) char       % UI label, e.g. 'WFS'
        synonyms      (1,:) cell       % alternative strings that map to this id

        % ─── Engine factory ───
        %   Signature: @(vs, SSD, loudspeakers, receiver, setup) → renderer
        %   A factory csak azt adja tovább, amit a renderer tényleg használ.
        factory                               % function_handle (no type constraint — MATLAB cannot default-construct one)

        % ─── UI panel ───
        %   Fully qualified class name, e.g. 'audioapp.ui.renderers.WFSPanel'
        %   '' → uses DirectPanel (fallback)
        panelClass    (1,:) char

        % ─── Geometry compatibility ───
        %   Cell of geometry strings this renderer supports.
        %   '*' means "accepts anything".
        geometries    (1,:) cell       % e.g. {'linear','circular','5.0'}
        prefGeometry  (1,:) char       % preferred/default geometry

        % ─── Default Renderer_setup fields ───
        %   Struct of fields that AppContext.defaultSetup() should include
        %   when this mode is active.  Merged into Setup.Renderer_setup.
        defaultSetup  (1,1) struct

        % ─── Optional metadata ───
        sortOrder     (1,1) double = 100  % UI popup ordering (lower = higher)
        hidden        (1,1) logical = false % true → not shown in UI popup
    end

    methods
        function obj = RendererDescriptor(varargin)
            % Name-value constructor:
            %   RendererDescriptor('id','wfs', 'uiName','WFS', ...)
            p = inputParser;
            p.addParameter('id',            '',     @ischar);
            p.addParameter('uiName',        '',     @ischar);
            p.addParameter('synonyms',      {},     @iscell);
            p.addParameter('factory',       @(vs,SSD,ls,rec,s) direct_renderer(vs,SSD));
            p.addParameter('panelClass',    '',     @ischar);
            p.addParameter('geometries',    {'*'},  @iscell);
            p.addParameter('prefGeometry',  'circular', @ischar);
            p.addParameter('defaultSetup',  struct(), @isstruct);
            p.addParameter('sortOrder',     100,    @isnumeric);
            p.addParameter('hidden',        false,  @islogical);
            p.parse(varargin{:});
            r = p.Results;
            obj.id           = r.id;
            obj.uiName       = r.uiName;
            obj.synonyms     = r.synonyms;
            obj.factory      = r.factory;
            obj.panelClass   = r.panelClass;
            obj.geometries   = r.geometries;
            obj.prefGeometry = r.prefGeometry;
            obj.defaultSetup = r.defaultSetup;
            obj.sortOrder    = r.sortOrder;
            obj.hidden       = r.hidden;

            % Unwrap nested cells: {{'a','b'}} → {'a','b'}
            % This happens because MATLAB name-value parsing wraps cell values
            % in an extra layer when passed as  'key', {{'a','b'}}
            while iscell(obj.geometries) && numel(obj.geometries)==1 && iscell(obj.geometries{1})
                obj.geometries = obj.geometries{1};
            end
            while iscell(obj.synonyms) && numel(obj.synonyms)==1 && iscell(obj.synonyms{1})
                obj.synonyms = obj.synonyms{1};
            end
        end

        % ─── Queries ───
        function tf = supportsGeometry(obj, geom)
            geom = lower(strtrim(geom));
            % Normalize common aliases
            switch geom
                case {'5.0 (circular)', '5_0'}, geom = '5.0';
                case {'5.1 (circular)'},         geom = '5.1';
                case {'atmos 5.1.4'},            geom = '5.1.4';
                case {'atmos 7.1.4'},            geom = '7.1.4';
                case {'soundbar (16 ch)'},       geom = 'soundbar';
                case {'genius array'},           geom = 'genius_array';
            end
            if any(strcmpi(obj.geometries, '*'))
                tf = true;
            else
                tf = any(strcmpi(obj.geometries, geom));
            end
        end

        function r = createRenderer(obj, vs, SSD, loudspeakers, receiver, setup)
            r = obj.factory(vs, SSD, loudspeakers, receiver, setup);
        end

        function pnl = createPanel(obj, parent, C)
            if isempty(obj.panelClass)
                pnl = audioapp.ui.renderers.DirectPanel(C);
            else
                ctor = str2func(obj.panelClass);
                pnl  = ctor(C);
            end
            pnl.build(parent);
        end
    end
end
