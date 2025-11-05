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
    end

    properties (Access = private)
        WarnedMissingDSP   (1,1) logical = false;
        WarnedMissingAudio (1,1) logical = false;
    end

    events
        SetupChanged
    end

    methods
        function obj = AppContext()
            % init defaults that don't depend on UI yet
            obj.Setup = obj.defaultSetup;
            obj.InputDisplayName = 'Input';

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
        end

        function rebuildSoundScene(obj)
            obj.SceneGUI.clearObject;
            obj.Scene = sound_scene(obj.SceneGUI, obj.Setup);
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
            S = struct( ...
                'N_in',1, ...
                'N_out',2, ...
                'MasterVolume',1,...
                'SampleRate',fs, ...
                'Block_size',1024*2, ...
                'HRTF',hrtf, ...
                'Binauralization',true, ...
                'HRTF_extrapolation','linear', ...
                'Downmixing_enabled',false, ...
                'Decorrelation',1, ...
                'Loudspeaker_setup',struct('Shape','linear','R',2,'N',128,'Height',1.7), ...
                'Rendering_mode','WFS', ...
                'Renderer_setup',struct('Antialiasing',false,'Tapering',0.25,'ReferenceMode','concentric', ...
                'DeltaY',2,'Rref',1,'RefPoint',[0 0],'RefDistance',2, ...
                'EnergyPreserving',true,'Spread',0,'HOAorder',10,...
                'Plant_model','HRTF','VS_model','HRTF','HRTF_database',hrtf), ...
                'Loudspeaker_type',struct('Shape','point_source','R',0.025), ...
                'Virtual_source_type',struct('Shape','point_source','R',0.01) );
        end

        function setRendererMode(obj, modeUI)
            geomCur = obj.Setup.Loudspeaker_setup.Shape;
            [modeOut, geomOut, ~] = audioapp.util.Rules.resolveOnRendererChange(modeUI, geomCur);

            obj.Setup.Rendering_mode = obj.rendererDictionary(modeOut);
            if ~strcmpi(geomOut, geomCur)
                obj.Setup.Loudspeaker_setup.Shape = geomOut;
                obj.Setup.Loudspeaker_setup.N     = audioapp.util.Rules.defaultNForGeometry(geomOut);
            end
            obj.rebuildSoundScene();
        end

        function modeOut = rendererDictionary(obj, modeIn)
            switch lower(modeIn)
                case {'direct playback','direct_playback'}
                    modeOut = 'direct_playback';
                case 'wfs'
                    modeOut = 'wfs';
                case 'vbap'
                    modeOut = 'vbap';
                case 'dbap'
                    modeOut = 'dbap';
                case 'hoa'
                    modeOut = 'hoa';
                case 'nfc_hoa'
                    modeOut = 'nfc_hoa';
                case 'ctc'
                    modeOut = 'ctc';
                case 'td_stereo'
                    modeOut = 'td_stereo';
                case {'dolby surround encoder', 'dolby_surround_decoder'}
                    modeOut = 'dolby_surround_decoder';
            end
        end

        function setGeometry(obj, geomAny, N)
            modeCur = obj.Setup.Rendering_mode;
            [modeOut, geomOut, ~] = audioapp.util.Rules.resolveOnGeometryChange(geomAny, modeCur);

            obj.Setup.Loudspeaker_setup.Shape = geomOut;
            obj.Setup.Loudspeaker_setup.N     = audioapp.util.Rules.defaultNForGeometry(geomOut);
            if ~strcmpi(modeOut, modeCur)
                obj.Setup.Rendering_mode = modeOut;
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
    end
end
