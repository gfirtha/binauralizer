classdef AppContext < handle
    properties
        % Toolboxes
        HasDSP  (1,1) logical = ~isempty(ver('dsp'));
        HasAudio(1,1) logical = ~isempty(ver('audio'));

        % Engine & services (created on demand)
        Engine  audioapp.model.AudioEngine
        Video   audioapp.services.VideoService

        % Scene & settings
        Scene               % sound_scene instance
        SceneGUI            % listener_space_axes
        Setup               % sound_scene_setup struct
        InputDisplayName    % string for UI
        mainAx

        % Designed source trajectories (Trajectories tab); kept across
        % scene rebuilds, followed during playback with the media time
        Trajectories audioapp.model.TrajectorySet
    end

    properties (Access = private)
        WarnedMissingDSP   (1,1) logical = false;
        WarnedMissingAudio (1,1) logical = false;
    end

    events
        SetupChanged
        GeometryChanged
    end

    methods
        function obj = AppContext()
            % init defaults that don't depend on UI yet
            obj.Setup = obj.defaultSetup;
            obj.InputDisplayName = 'Input';
            obj.Trajectories = audioapp.model.TrajectorySet();

            % Engine & playback capability + user warning
            obj.ensureEngine();                           % creates engine
            obj.updateEnginePlaybackCapability();         % enable/disable playback
            obj.warnIfPlaybackDisabled();                 % show dialog(s) if needed

            % Safe audio load (won’t start playback here)
            try
                obj.Engine.loadAudio('SoundSamples/gitL.wav');
            catch
                % non-fatal; just continue
            end

            % Scene
            obj.initializeScene;
            try
                obj.SceneGUI.draw_renderer_properties(obj.Scene.scene_renderer.SFS_renderer{1});
            catch
            end
        end

        function initializeScene(obj)
            obj.SceneGUI = listener_space_axes();
            obj.Scene    = sound_scene(obj.SceneGUI, obj.Setup);
            obj.Scene.TestChannelFcn = @(ch) obj.onTestChannel(ch);
            obj.Scene.OnGeometryChangedFcn = @() notify(obj, 'GeometryChanged');
        end

        function rebuildSoundScene(obj)
            obj.SceneGUI.clearObject;
            obj.Scene = sound_scene(obj.SceneGUI, obj.Setup);
            obj.Scene.TestChannelFcn = @(ch) obj.onTestChannel(ch);
            obj.Scene.OnGeometryChangedFcn = @() notify(obj, 'GeometryChanged');
            try
                obj.SceneGUI.draw_renderer_properties(obj.Scene.scene_renderer.SFS_renderer{1});
            catch
            end
            notify(obj,'SetupChanged');
        end

        function ensureEngine(obj)
            if isempty(obj.Engine) || ~isvalid(obj.Engine)
                obj.Engine = audioapp.model.AudioEngine(obj);
            end
            % Keep playback policy in sync even if ensureEngine() is called later
            obj.updateEnginePlaybackCapability();
        end

        function ensureVideo(obj)
            if isempty(obj.Video) || ~isvalid(obj.Video)
                obj.Video = audioapp.services.VideoService(obj);
            end
        end

        function S = defaultSetup(obj)
            hrtf = SOFAload('BuK_ED_corr.sofa');
            fs   = hrtf.Data.SamplingRate;

            % Common / global defaults
            S = struct( ...
                'N_in',1, ...
                'SourceNames',{{}}, ...
                'N_out',2, ...
                'MasterVolume',1,...
                'SampleRate',fs, ...
                'Block_size',1024*2, ...
                'HRTF',hrtf, ...
                'Binauralization',true, ...
                'HRTF_extrapolation','linear', ...
                'Downmixing_enabled',false, ...
                'Decorrelation',1, ...
                'ChannelMap',[], ...   % [] = use layout default from get_default_layout
                'Loudspeaker_setup',struct('Shape','5.0','R',2,'N',5,'Height',1.7), ...
                'Rendering_mode','vbap', ...
                'Renderer_setup',struct('HRTF_database',hrtf), ...
                'Loudspeaker_type',struct('Shape','point_source','R',0.025), ...
                'Virtual_source_type',struct('Shape','point_source','R',0.01) );

            % Source separation / upmix (Source separation tab)
            S.Separation = struct('Method', 'stft_objects', 'K', 4, 'Width', 30, 'Base', 30, ...
                'AmbAz', 110, 'AmbGain', 1, 'Tau', 0.05);

            % Merge renderer-specific defaults from the registry
            reg = audioapp.RendererRegistry.instance();
            rendererDefaults = reg.defaultSetup(S.Rendering_mode);
            fns = fieldnames(rendererDefaults);
            for k = 1:numel(fns)
                S.Renderer_setup.(fns{k}) = rendererDefaults.(fns{k});
            end
        end

        function setRendererMode(obj, modeUI)
            reg = audioapp.RendererRegistry.instance();
            geomCur = obj.Setup.Loudspeaker_setup.Shape;
            [modeOut, geomOut, ~] = reg.resolveOnRendererChange(modeUI, geomCur);

            obj.Setup.Rendering_mode = modeOut;
            if ~strcmpi(geomOut, geomCur)
                obj.Setup.Loudspeaker_setup.Shape = geomOut;
                obj.Setup.Loudspeaker_setup.N     = audioapp.util.Rules.defaultNForGeometry(geomOut);
            end

            % Merge renderer-specific defaults into Renderer_setup
            rendererDefaults = reg.defaultSetup(modeOut);
            fns = fieldnames(rendererDefaults);
            for k = 1:numel(fns)
                if ~isfield(obj.Setup.Renderer_setup, fns{k})
                    obj.Setup.Renderer_setup.(fns{k}) = rendererDefaults.(fns{k});
                end
            end

            obj.rebuildSoundScene();
        end

        function modeOut = rendererDictionary(~, modeIn)
            % LEGACY wrapper — delegates to RendererRegistry
            modeOut = audioapp.RendererRegistry.instance().resolveId(modeIn);
        end

        function setGeometry(obj, geomAny, N)
            reg = audioapp.RendererRegistry.instance();
            modeCur = obj.Setup.Rendering_mode;
            [modeOut, geomOut, changed] = reg.resolveOnGeometryChange(geomAny, modeCur);

            obj.Setup.Loudspeaker_setup.Shape = geomOut;
            obj.Setup.Loudspeaker_setup.N     = audioapp.util.Rules.defaultNForGeometry(geomOut);
            if strcmpi(geomOut, 'soundbar')
                % home-cinema scale: fronts on a 1.5 m circle, the bar between
                % them; WFS amplitude referenced at the listener (bar depth)
                obj.Setup.Loudspeaker_setup.R = 1.5;
                obj.Setup.Renderer_setup.DeltaY = 1.5*cosd(30);
            end
            if ~strcmpi(modeOut, modeCur)
                obj.Setup.Rendering_mode = modeOut;
                % the auto-selected renderer's defaults, as setRendererMode does
                rendererDefaults = reg.defaultSetup(modeOut);
                fns = fieldnames(rendererDefaults);
                for k = 1:numel(fns)
                    if ~isfield(obj.Setup.Renderer_setup, fns{k})
                        obj.Setup.Renderer_setup.(fns{k}) = rendererDefaults.(fns{k});
                    end
                end
            end
            obj.rebuildSoundScene();
        end

        function touchSetup(obj)
            notify(obj,'SetupChanged');
        end
    end

    methods (Access = private)
        function updateEnginePlaybackCapability(obj)
            canPlay = obj.HasDSP && obj.HasAudio;
            if isempty(obj.Engine) || ~isvalid(obj.Engine), return; end
            try
                % Prefer method if your engine exposes it
                if ismethod(obj.Engine,'setPlaybackEnabled')
                    obj.Engine.setPlaybackEnabled(canPlay);
                % Or a property
                elseif isprop(obj.Engine,'PlaybackEnabled')
                    obj.Engine.PlaybackEnabled = canPlay;
                end
            catch
                % non-fatal: engine may not expose such control
            end
        end

        function warnIfPlaybackDisabled(obj)
            if ~obj.HasDSP && ~obj.WarnedMissingDSP
                obj.WarnedMissingDSP = true;
                warndlg( ...
                    ['DSP System Toolbox not found.' newline ...
                     'Streaming input and real-time playback are disabled.' newline ...
                     'You can still open files, inspect sources, and render fields.'], ...
                    'DSP toolbox missing','modal');
            end
            if ~obj.HasAudio && ~obj.WarnedMissingAudio
                obj.WarnedMissingAudio = true;
                warndlg( ...
                    ['Audio Toolbox not found.' newline ...
                     'Real-time audio output is disabled.'], ...
                    'Audio toolbox missing','modal');
            end
        end

        function onTestChannel(obj, chIdx)
            % Called when user double-clicks a loudspeaker
            try
                obj.ensureEngine();
                obj.Engine.testChannel(chIdx);
            catch ME
                fprintf('Speaker test failed: %s\n', ME.message);
            end
        end
    end
end
