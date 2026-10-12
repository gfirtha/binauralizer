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

        % Loudspeakers of the Loudspeaker designer (designed_loudspeaker by
        % name): assigned to SSD loudspeakers on the SSD setup tab. Filled from
        % Data/Loudspeakers (design files with an exported 'speaker') and
        % by the designer's "Send to Binauralizer"
        Speakers
    end

    properties (Access = private)
        WarnedMissingDSP   (1,1) logical = false;
        WarnedMissingAudio (1,1) logical = false;
    end

    events
        SetupChanged
        GeometryChanged
        SpeakersChanged     % designed loudspeakers: library or placements changed
    end

    methods
        function obj = AppContext()
            % init defaults that don't depend on UI yet
            obj.Setup = obj.defaultSetup;
            obj.InputDisplayName = 'Input';
            obj.Trajectories = audioapp.model.TrajectorySet();
            obj.Speakers = containers.Map('KeyType', 'char', 'ValueType', 'any');
            obj.scanSpeakerDir();

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
            if ~isempty(obj.Video) && isvalid(obj.Video) && obj.Video.hasVideo()
                obj.Video.arrangeScene();             % TV at the centre, centre under the TV
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

            % Designed loudspeakers in the scene (speaker_placement handles)
            S.DesignedSpeakers = {};

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

        % --- Designed loudspeakers ---
        function addSpeaker(obj, S)
            % Register designed_loudspeaker S (replaces one of the same
            % name; its placements take the new model, or are removed when
            % the number of drivers changed)
            obj.Speakers(S.Name) = S;
            L = obj.Setup.DesignedSpeakers;
            for q = numel(L):-1:1
                P = L{q};
                if strcmp(P.Model.Name, S.Name) && P.Model ~= S
                    if S.nChannels == numel(P.Elements)
                        obj.Scene.detach_speaker(P); P.Model = S; obj.Scene.attach_speaker(P);
                    else
                        obj.removeSpeaker(P);
                    end
                end
            end
            notify(obj, 'SpeakersChanged');
        end

        function n = speakerNames(obj)
            n = sort(obj.Speakers.keys);
        end

        function S = loadSpeakerFile(obj, file)
            S = designed_loudspeaker.load(file);
            obj.addSpeaker(S);
        end

        function scanSpeakerDir(obj)
            % design files with an exported loudspeaker in Data/Loudspeakers
            d = audioapp.AppContext.speakerDir();
            if ~isfolder(d), return; end
            F = dir(fullfile(d, '*.mat'));
            for k = 1:numel(F)
                file = fullfile(F(k).folder, F(k).name);
                if ~designed_loudspeaker.isSpeakerFile(file), continue; end
                try
                    S = designed_loudspeaker.load(file);
                    obj.Speakers(S.Name) = S;
                catch ME
                    warning('AppContext:speaker', 'Skipping %s: %s', F(k).name, ME.message);
                end
            end
        end

        function P = placeSpeaker(obj, S, elements)
            % SSD loudspeakers (loudspeaker_array indices, one per driver,
            % in driver order) become the drivers of designed loudspeaker S,
            % placed where they are (speaker_placement.onElements).
            % Placements sharing a loudspeaker are removed.
            if numel(elements) ~= S.nChannels
                error('AppContext:speaker', '"%s" has %d drivers: select %d loudspeakers.', S.Name, S.nChannels, S.nChannels);
            end
            L = obj.Setup.DesignedSpeakers;
            for q = numel(L):-1:1
                if any(ismember(L{q}.Elements, elements)), obj.removeSpeaker(L{q}); end
            end
            A = obj.Scene.loudspeaker_array(elements);
            pos = cell2mat(cellfun(@(l) l.position(1:3), A(:), 'UniformOutput', false));
            ori = cell2mat(cellfun(@(l) l.orientation(1:2), A(:), 'UniformOutput', false));
            P = speaker_placement.onElements(S, elements, pos, ori, obj.Scene.receiver.position);
            obj.Setup.DesignedSpeakers{end+1} = P;
            obj.Scene.attach_speaker(P);
            notify(obj, 'SpeakersChanged');
        end

        function removeSpeaker(obj, P)
            obj.Scene.detach_speaker(P);
            L = obj.Setup.DesignedSpeakers;
            obj.Setup.DesignedSpeakers = L(~cellfun(@(x) x == P, L));
            notify(obj, 'SpeakersChanged');
        end

        function P = placementOfElement(obj, e)
            % the placement whose driver SSD loudspeaker e is ([] if none)
            P = [];
            L = obj.Setup.DesignedSpeakers;
            for q = 1:numel(L)
                if any(L{q}.Elements == e), P = L{q}; return; end
            end
        end
    end

    methods (Static)
        function d = speakerDir()
            d = fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Data', 'Loudspeakers');
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
