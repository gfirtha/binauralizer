classdef AudioEngine < handle
    properties
        C                       % AppContext
        Input                   % dsp.AudioFileReader (or [])
        SRC                     % dsp.SampleRateConverter (or [])
        RB                      % dsp.AsyncBuffer
        Out                     % audioDeviceWriter (or [])
        stop_now  logical = false
        pause_now logical = false
        is_paused logical = false
        seek_to   double  = -1    % negative = no pending seek
        totalSamplesOut double = 0  % media position [samples at the processing rate]

        % ── Playback speed (tape-style varispeed: slower = lower pitch) ──
        Speed double = 1
        VS = []                  % audioapp.dsp.Varispeed

        MediaDuration double = 0   % total duration in seconds
        MediaFilename char = ''    % path to currently loaded file
        InputFs       double = 0   % input file sample rate

        MasterEQ   audioapp.dsp.MasterEQ

        % ── Object-based input (Dolby Atmos E-AC-3 JOC via Cavern) ──
        Objects  = []           % cavern_reader, or [] for channel-based files
        ObjectLevel = []        % smoothed per-object RMS (display only)
        BeamAnalyzer             % audioapp.model.DirectBeamAnalyzer
        BeamVizCounter double = 0
        BeamVizStride  double = 2

        % ── Loudspeaker test signal (double-click on a loudspeaker) ──
        Running   logical = false   % a playback loop is feeding the output
        TestTone  = []              % cached test signal at TestToneFs
        TestToneFs double = 0

        % ── Microphone input ──
        MicInput                % audioDeviceReader (or [])
        UseMic    logical = false
        MicDeviceName  char = ''
        MicNumChannels double = 1
    end

    methods
        function obj = AudioEngine(ctx)
            obj.C = ctx;
            obj.ensureOutput();
            if obj.C.HasDSP
                obj.RB = dsp.AsyncBuffer(20*obj.C.Setup.Block_size);
            else
                obj.RB = [];
            end
            obj.MasterEQ = audioapp.dsp.MasterEQ(ctx);  % disabled by default
            obj.BeamAnalyzer = audioapp.model.DirectBeamAnalyzer(2);

        end

        % ===== Device / IO =================================================
        function ensureOutput(obj)
            if ~obj.C.HasAudio
                obj.Out = [];
                return;
            end
            if isempty(obj.Out)
                obj.Out = audioDeviceWriter;  % system default
            end
            try, obj.Out.SampleRate = obj.C.Setup.SampleRate; catch, end
            try, obj.Out.BufferSize = obj.C.Setup.Block_size; catch, end
        end

        function initIO(obj, fs_proc, blk)
            obj.C.Setup.SampleRate = fs_proc;
            obj.C.Setup.Block_size = blk;

            obj.ensureOutput();
            if ~isempty(obj.Out)
                try, obj.Out.SampleRate = fs_proc;  catch, end
                try, obj.Out.BufferSize = blk;      catch, end
            end

            if obj.C.HasDSP
                obj.RB = dsp.AsyncBuffer(20*blk);       % room for a slowed-down input block
            else
                obj.RB = [];
            end

            if ~isempty(obj.MasterEQ) && isvalid(obj.MasterEQ)
                obj.MasterEQ.setFs(fs_proc);
            end
            if ~isempty(obj.BeamAnalyzer) && isvalid(obj.BeamAnalyzer)
                obj.BeamAnalyzer.reset();
            end

        end

        function refreshOutput(obj)
            obj.ensureOutput();
        end

        % ===== Input device listing ========================================
        function list = getInputDevices(obj)
            % Return cell array of available audio input device names
            if ~obj.C.HasAudio
                list = {'(no input devices found)'};
                return;
            end
            try
                list = getAudioDevices(audioDeviceReader);
                if isstring(list), list = cellstr(list); end
                if isempty(list), list = {'(no input devices found)'}; end
            catch
                try
                    S = audiodevinfo;
                    list = arrayfun(@(k) S.input(k).Name, 1:numel(S.input), 'uni', 0);
                    if isempty(list), list = {'(no input devices found)'}; end
                catch
                    list = {'(no input devices found)'};
                end
            end
        end

        % ===== Microphone input ============================================
        function loadMicInput(obj, deviceName, numChannels)
            % Set up microphone as the live input source.
            %   deviceName  - name of the audio input device
            %   numChannels - number of channels to capture (default: 1)

            if nargin < 3 || isempty(numChannels), numChannels = 1; end

            if ~obj.C.HasDSP || ~obj.C.HasAudio
                error('AudioEngine:noToolbox', ...
                    'DSP System Toolbox and Audio Toolbox are required for mic input.');
            end

            % Release any existing file reader
            if ~isempty(obj.Input)
                try release(obj.Input); catch, end
                obj.Input = [];
            end
            % Release any existing mic reader
            if ~isempty(obj.MicInput)
                try release(obj.MicInput); catch, end
                obj.MicInput = [];
            end

            % Create audioDeviceReader
            fs = obj.C.Setup.SampleRate;
            blk = obj.C.Setup.Block_size;

            obj.MicInput = audioDeviceReader;
            obj.MicInput.Device = deviceName;
            obj.MicInput.SampleRate = fs;
            obj.MicInput.SamplesPerFrame = blk;
            try
                obj.MicInput.NumChannels = numChannels;
            catch
                % Some versions use ChannelMapping instead
                try
                    obj.MicInput.ChannelMapping = 1:numChannels;
                catch
                end
            end

            obj.UseMic = true;
            obj.MicDeviceName = deviceName;
            obj.MicNumChannels = numChannels;
            obj.MediaFilename = '';
            obj.MediaDuration = 0;
            obj.InputFs = fs;
            obj.SRC = [];  % no SRC needed — mic runs at process sample rate

            % Update setup
            obj.C.Setup.N_in = numChannels;
            [~, obj.C.Setup.SourceNames] = separator_base.channel_layout(numChannels);
            obj.C.InputDisplayName = sprintf('[MIC] %s (%dch)', deviceName, numChannels);

            % Reset IO
            obj.initIO(fs, blk);

            fprintf('Mic input activated: %s, %d ch, %d Hz, block=%d\n', ...
                deviceName, numChannels, fs, blk);
        end

        function clearMicInput(obj)
            % Deactivate mic input mode
            if ~isempty(obj.MicInput)
                try release(obj.MicInput); catch, end
                obj.MicInput = [];
            end
            obj.UseMic = false;
            obj.MicDeviceName = '';
        end

        % ===== Unified media loader =======================================
        function kind = loadMedia(obj, filename)
            % kind = 'audio' or 'video'

            % Clear mic mode
            obj.clearMicInput();

            obj.C.InputDisplayName = filename;
            obj.C.Setup.SourceNames = {};
            obj.MediaFilename = filename;
            obj.MediaDuration = 0;
            obj.totalSamplesOut = 0;
            obj.seek_to = -1;

            % Determine if this is a video file by extension
            [~,~,ext] = fileparts(filename);
            videoExts = {'.mp4','.mkv','.avi','.webm','.mov','.flv','.wmv','.m4v','.ts','.mpg','.mpeg'};
            isVideoFile = any(strcmpi(ext, videoExts));

            % Open video (only for video files)
            hasVid = false;
            if isVideoFile
                try
                    obj.C.ensureVideo();
                    obj.C.Video.open(filename);
                    hasVid = obj.C.Video.HasVLC;  % only count as video if VLC works
                catch
                    if ~isempty(obj.C.Video), obj.C.Video.close(); end
                end
            end
            kind = ternary(hasVid, 'video', 'audio');

            % ── Object-based audio (Atmos) through Cavern, if it can ──
            obj.Objects = [];
            if obj.tryOpenObjects(filename)
                obj.Input = [];
                obj.initIO(obj.C.Setup.SampleRate, obj.C.Setup.Block_size);
                return;
            end

            % ── Audio reader ──
            % Strategy:
            %   1. Try dsp.AudioFileReader directly (works for wav, mp3, mp4, ogg, flac)
            %   2. Try MATLAB audioread (works for similar formats)
            %   3. Try ffmpeg to extract audio to temp WAV (works for mkv, avi, webm, etc.)
            if obj.C.HasDSP
                fs_in = [];
                audioOK = false;

                % Attempt 1: dsp.AudioFileReader directly
                try
                    obj.Input = dsp.AudioFileReader( ...
                        filename, ...
                        'SamplesPerFrame', obj.C.Setup.Block_size, ...
                        'PlayCount', 1);
                    fs_in = obj.Input.SampleRate;
                    fileInfo = info(obj.Input);
                    obj.C.Setup.N_in = fileInfo.NumChannels;
                    audioOK = true;
                catch
                end

                % Attempt 2: audioread → temp WAV
                if ~audioOK
                    try
                        [y, Fs] = audioread(filename);
                        tmp = fullfile(tempdir, ['binauralizer_audio_' char(java.util.UUID.randomUUID) '.wav']);
                        audiowrite(tmp, y, Fs);
                        obj.Input = dsp.AudioFileReader(tmp, ...
                            'SamplesPerFrame', obj.C.Setup.Block_size, 'PlayCount', 1);
                        fs_in = obj.Input.SampleRate;
                        fInfo = info(obj.Input);
                        obj.C.Setup.N_in = fInfo.NumChannels;
                        audioOK = true;
                    catch
                    end
                end

                % Attempt 3: ffmpeg extraction → temp WAV
                if ~audioOK
                    tmp = obj.extractAudioWithFFmpeg(filename);
                    if ~isempty(tmp)
                        try
                            obj.Input = dsp.AudioFileReader(tmp, ...
                                'SamplesPerFrame', obj.C.Setup.Block_size, 'PlayCount', 1);
                            fs_in = obj.Input.SampleRate;
                            fInfo = info(obj.Input);
                            obj.C.Setup.N_in = fInfo.NumChannels;
                            audioOK = true;
                        catch
                        end
                    end
                end

                if ~audioOK
                    obj.Input = [];
                    error('AudioEngine:openFailed', ...
                        ['Could not open audio from: %s\n\n' ...
                         'For mkv/webm/avi support, install ffmpeg and add it to your system PATH.\n' ...
                         'Download: https://ffmpeg.org/download.html'], filename);
                end

                [~, obj.C.Setup.SourceNames] = separator_base.channel_layout(obj.C.Setup.N_in);
                obj.initSRC(fs_in, obj.C.Setup.SampleRate);
                obj.InputFs = fs_in;

                % Detect duration using audioinfo (most reliable)
                try
                    ai = audioinfo(filename);
                    obj.MediaDuration = ai.Duration;
                catch
                    % Fallback: estimate from dsp.AudioFileReader
                    try
                        fInfo = info(obj.Input);
                        if isfield(fInfo,'NumSamples') && fInfo.NumSamples > 0
                            obj.MediaDuration = fInfo.NumSamples / fs_in;
                        end
                    catch
                    end
                end
            else
                obj.Input = [];
            end

            % Reset IO so ring buffer depth matches current block size
            obj.initIO(obj.C.Setup.SampleRate, obj.C.Setup.Block_size);
        end

        % tiny wrappers (optional)
        function loadAudio(obj, filename), obj.loadMedia(filename); end
        function loadVideo(obj, filename), obj.loadMedia(filename); end

        function initSRC(obj, fs_in, fs_proc)
            if isempty(fs_in) || isempty(fs_proc) || ~obj.C.HasDSP
                obj.SRC = [];
                return;
            end
            if fs_in ~= fs_proc
                obj.SRC = dsp.SampleRateConverter('InputSampleRate',fs_in,'OutputSampleRate',fs_proc);
                try, obj.SRC.StopbandAttenuation = 80; catch, end
            else
                obj.SRC = [];
            end
        end

        % ===== Device management API (used by OutputTab) ===================
        function list = getOutputDevices(obj)
            if ~obj.C.HasAudio, list = {'(no output devices found)'}; return; end
            try
                list = getAudioDevices(audioDeviceWriter);
                if isstring(list), list = cellstr(list); end
                if isempty(list), list = {'(no output devices found)'}; end
            catch
                try
                    S = audiodevinfo;
                    list = arrayfun(@(k) S.output(k).Name, 1:numel(S.output), 'uni', 0);
                    if isempty(list), list = {'(no output devices found)'}; end
                catch
                    list = {'(no output devices found)'};
                end
            end
        end

        function name = getCurrentDeviceName(obj)
            if isempty(obj.Out), name = '(no device)'; return; end
            try
                name = obj.Out.Device;
                if isempty(name), name = '(system default)'; end
            catch
                name = '(unknown)';
            end
        end

        function setOutputDevice(obj, name)
            if ~obj.C.HasAudio, return; end
            obj.ensureOutput();
            try, release(obj.Out); catch, end
            try, obj.Out.Device = name; catch, end
            obj.refreshOutput();
        end

        function info = getCurrentDeviceInfo(obj)
            info = struct('Name', obj.getCurrentDeviceName(), ...
                'MaxOutputChannels', NaN, ...
                'SampleRate', obj.C.Setup.SampleRate, ...
                'BufferSize', obj.C.Setup.Block_size);
            if isempty(obj.Out), return; end
            try
                inf = obj.Out.info;
                info.MaxOutputChannels = inf.MaximumOutputChannels;
            catch
            end
            try, info.SampleRate = obj.Out.SampleRate;  catch, end
            try, info.BufferSize = obj.Out.BufferSize;  catch, end
        end

        % ===== Playback (self-looping with A/V sync) =======================
        function play(obj)
            % ── Mic mode: live input ──
            if obj.UseMic
                obj.playMic();
                return;
            end

            % ── File mode (original behaviour) ──
            if ~(obj.C.HasDSP && obj.C.HasAudio && (~isempty(obj.Input) || ~isempty(obj.Objects)))
                errordlg('Playback disabled: missing DSP/Audio toolbox or no input.','Cannot play','modal');
                return;
            end

            blk = obj.C.Setup.Block_size;
            fs  = obj.C.Setup.SampleRate;
            obj.initIO(fs, blk);

            % Configure audioDeviceWriter
            if ~isempty(obj.Out)
                try release(obj.Out); catch, end
                try obj.Out.SampleRate = fs; catch, end
                try obj.Out.BufferSize = blk; catch, end
            end

            % Reset audio input
            if ~isempty(obj.Objects)
                obj.Objects.seek(0);
            else
                reset(obj.Input);
            end

            hasVid = ~isempty(obj.C.Video) && obj.C.Video.hasVideo();
            % VLC will be launched by playStart() after first audio block

            obj.stop_now  = false;
            obj.pause_now = false;
            obj.is_paused = false;
            obj.seek_to   = -1;
            obj.totalSamplesOut = 0;
            obj.VS = audioapp.dsp.Varispeed();
            drawnowCounter = 0;
            vlcStarted = false;
            obj.Running = true;
            runningGuard = onCleanup(@() obj.setRunning(false)); %#ok<NASGU>

            while ~obj.stop_now

                % ── Handle seek request ──
                if obj.seek_to >= 0
                    seekSec = obj.seek_to;
                    obj.seek_to = -1;
                    obj.totalSamplesOut = round(seekSec * fs);

                    % Seek audio: reset and skip forward
                    try
                        skipFrames = floor(seekSec * obj.InputFs / blk);
                        if ~isempty(obj.Objects)
                            % Container seek in CavernServe takes 10-20 s on a
                            % long film: wait for it without blocking the UI
                            % (the output underruns to silence meanwhile)
                            obj.Objects.seekAsync(skipFrames*blk);
                            while ~obj.Objects.seekDone() && ~obj.stop_now
                                drawnow limitrate
                                java.lang.Thread.sleep(30);
                            end
                            seekSec = obj.Objects.samplesRead / fs;   % where it landed
                            obj.totalSamplesOut = obj.Objects.samplesRead;
                        else
                            reset(obj.Input);
                            for sk = 1:skipFrames
                                if isDone(obj.Input), break; end
                                obj.Input();  % read and discard
                            end
                        end
                    catch
                    end

                    % Flush ring buffer
                    if ~isempty(obj.RB) && obj.RB.NumUnreadSamples > 0
                        read(obj.RB, obj.RB.NumUnreadSamples);
                    end
                    obj.VS.reset();

                    % Seek video
                    if hasVid
                        obj.C.Video.seekTo(seekSec);
                    end
                end

                % ── Handle pause ──
                if obj.pause_now && ~obj.is_paused
                    obj.is_paused = true;
                    obj.pause_now = false;
                    if hasVid, obj.C.Video.vlcPause(); end
                end
                if obj.is_paused
                    drawnow limitrate
                    java.lang.Thread.sleep(50);
                    continue;
                end

                if ~isempty(obj.Objects)
                    % Object-based input: samples + positions from Cavern
                    % (loops by itself at the end of the stream)
                    y = obj.readObjects(blk);
                else
                    % ── Loop the audio file ──
                    if isDone(obj.Input)
                        reset(obj.Input);
                        obj.totalSamplesOut = 0;
                        if hasVid
                            % Just rewind VLC, don't restart — it's --repeat anyway
                            obj.C.Video.rewind();
                        end
                    end

                    % Read from file into ring buffer
                    x = obj.Input();
                    if ~isempty(obj.SRC), y = obj.SRC(x); else, y = x; end
                end
                if obj.Speed ~= 1 || obj.VS.isActive()
                    y = obj.VS.process(y, obj.Speed);   % slow motion: stretched in time
                end
                write(obj.RB, y);
                obj.applyTrajectories();            % designed source paths (media time)

                % Process all available full blocks
                while obj.RB.NumUnreadSamples >= blk && ~obj.stop_now
                    inFrame = read(obj.RB, blk);

                    if ~isempty(obj.MasterEQ) && isvalid(obj.MasterEQ)
                        inFrame = obj.MasterEQ.apply(inFrame);
                    end

                    obj.updateBeamOverlay(inFrame);

                    outFrame = obj.C.Setup.MasterVolume * obj.C.Scene.render_sound_scene( ...
                        inFrame, obj.C.Setup.Binauralization, ...
                        obj.C.Setup.Downmixing_enabled, obj.C.Setup.Decorrelation);

                    obj.Out(outFrame);
                    % media time advances Speed samples per output sample
                    obj.totalSamplesOut = obj.totalSamplesOut + size(outFrame,1)*obj.Speed;

                    % Start VLC after first audio block
                    if hasVid && ~vlcStarted
                        obj.C.Video.playStart();
                        vlcStarted = true;
                    end
                end

                % Let UI callbacks run
                drawnowCounter = drawnowCounter + 1;
                % every block: rendering takes ~1 ms now, and 'limitrate'
                % caps the redraws; every 8th block made dragging lag ~0.3 s
                if drawnowCounter >= 1
                    drawnow limitrate
                    drawnowCounter = 0;
                end
            end

            % Drain partial block
            rem = obj.RB.NumUnreadSamples;
            if rem > 0
                pad = zeros(blk-rem, size(y,2), 'like', y);
                last = [read(obj.RB, rem); pad];
                obj.updateBeamOverlay(last);
                outFrame = obj.C.Scene.render_sound_scene( ...
                    last, obj.C.Setup.Binauralization, ...
                    obj.C.Setup.Downmixing_enabled, obj.C.Setup.Decorrelation);
                obj.Out(outFrame);
            end

            try release(obj.Out); catch, end
            try release(obj.Input); catch, end
            % Kill VLC entirely on stop — next Play will launch fresh
            if hasVid
                obj.C.Video.playStop();
            end
        end

        % ===== Source separation / upmix as object input ===================
        function openSeparation(obj, cfg)
            % Play the loaded media through a separator: its objects become
            % the virtual sources (same path as the Atmos objects)
            if isempty(obj.MediaFilename)
                error('AudioEngine:noMedia', 'Load a media file first.');
            end
            obj.stop();
            fs = obj.C.Setup.SampleRate;
            rd = separation_reader(obj.MediaFilename, make_separator(cfg), fs);
            obj.Objects = rd;
            obj.Input = [];
            obj.SRC = [];
            obj.InputFs = fs;
            obj.MediaDuration = rd.duration();
            obj.totalSamplesOut = 0;
            obj.C.Setup.N_in = rd.nObj;
            obj.C.Setup.SourceNames = rd.Names;
            obj.C.InputDisplayName = sprintf('[Separated: %s, %d objects] %s', cfg.Method, rd.nObj, obj.MediaFilename);
            obj.initIO(fs, obj.C.Setup.Block_size);
        end

        function tf = isSeparated(obj)
            tf = ~isempty(obj.Objects) && isa(obj.Objects, 'separation_reader');
        end

        function bypassSeparation(obj)
            % Back to the media's own channels (or Atmos objects)
            if ~obj.isSeparated(), return; end
            obj.stop();
            obj.loadMedia(obj.MediaFilename);
        end

        % ===== Loudspeaker test signal ========================================
        function testChannel(obj, k)
            % Play the test signal from loudspeaker k (loudspeaker_array
            % index). During playback the running loop mixes it in;
            % otherwise a short loop renders it alone.
            sc = obj.C.Scene;
            if isempty(sc) || k < 1 || k > numel(sc.loudspeaker_array), return; end
            sc.start_test_signal(k, obj.testTone(obj.C.Setup.SampleRate));
            if obj.Running, return; end
            if ~obj.C.HasAudio, return; end
            obj.runTestLoop();
        end

        function x = testTone(obj, fs)
            % Classic Windows chimes (C:\Windows\Media\chimes.wav), mono,
            % peak 0.5; a synthetic bell if the file is not there
            if isempty(obj.TestTone) || obj.TestToneFs ~= fs
                try
                    [x, fs0] = audioread(fullfile(getenv('WINDIR'), 'Media', 'chimes.wav'));
                    x = mean(x, 2);
                    if fs0 ~= fs
                        t0 = (0:numel(x)-1)'/fs0;
                        x = interp1(t0, x, (0:floor(t0(end)*fs))'/fs, 'pchip');
                    end
                catch
                    t = (0:round(1.2*fs)-1)'/fs;
                    x = sin(2*pi*880*t).*exp(-4*t) + 0.5*sin(2*pi*1320*t).*exp(-5*t) ...
                        + 0.3*sin(2*pi*1760*t).*exp(-6*t);
                end
                obj.TestTone = 0.5*x/max(abs(x));
                obj.TestToneFs = fs;
            end
            x = obj.TestTone;
        end

        function runTestLoop(obj)
            blk = obj.C.Setup.Block_size;
            sc  = obj.C.Scene;
            obj.ensureOutput();
            if isempty(obj.Out), return; end
            try release(obj.Out); catch, end
            try obj.Out.SampleRate = obj.C.Setup.SampleRate; catch, end
            try obj.Out.BufferSize = blk; catch, end
            obj.Running = true;
            runningGuard = onCleanup(@() obj.setRunning(false)); %#ok<NASGU>
            silence = zeros(blk, numel(sc.virtual_sources));
            while sc.test_active()
                outFrame = obj.C.Setup.MasterVolume * sc.render_sound_scene( ...
                    silence, obj.C.Setup.Binauralization, ...
                    obj.C.Setup.Downmixing_enabled, obj.C.Setup.Decorrelation);
                obj.Out(outFrame);
                drawnow limitrate          % another double-click restarts the tone
            end
            try release(obj.Out); catch, end
        end

        function updateBeamOverlay(obj, inFrame)
            if isempty(obj.BeamAnalyzer) || ~isvalid(obj.BeamAnalyzer)
                obj.BeamAnalyzer = audioapp.model.DirectBeamAnalyzer(2);
            end
            obj.BeamVizCounter = obj.BeamVizCounter + 1;
            if obj.BeamVizCounter < obj.BeamVizStride
                return;
            end
            obj.BeamVizCounter = 0;

            try
                names = {};
                if isfield(obj.C.Setup, 'SourceNames')
                    names = obj.C.Setup.SourceNames;
                end
                data = obj.BeamAnalyzer.update(inFrame, names);
                gui = obj.C.SceneGUI;
                if isempty(gui) || ~isvalid(gui)
                    return;
                end
                if data.active
                    gui.update_beam_polar(data.anglesDeg, data.pattern, data.musicPattern, data.phantomPattern);
                else
                    gui.set_beam_polar_visible(false);
                end
            catch ME
                fprintf(2, 'Beam overlay update failed: %s\n', ME.message);
            end
        end

        function applyTrajectories(obj, t)
            % Move the sources that follow a designed trajectory to their
            % positions at media time t (default: the current output time)
            T = obj.C.Trajectories;
            if isempty(T) || ~T.hasAssignments(), return; end
            if nargin < 2, t = obj.totalSamplesOut / max(obj.C.Setup.SampleRate, 1); end
            sc = obj.C.Scene;
            [idx, P] = T.positions(t, numel(sc.virtual_sources));
            if ~isempty(idx), sc.apply_source_positions(idx, P); end
        end

        function setRunning(obj, tf)
            if isvalid(obj), obj.Running = tf; end
        end

        function stop(obj)
            obj.stop_now = true;
            obj.pause_now = false;
            obj.is_paused = false;
        end

        function pausePlayback(obj)
            if ~obj.is_paused
                obj.pause_now = true;
            end
        end

        function resumePlayback(obj)
            if obj.is_paused
                obj.is_paused = false;
                obj.pause_now = false;
            end
        end

        function setSpeed(obj, s)
            % Playback speed (0 < s <= 1): audio slowed tape-style, the
            % video (VLC) at the same rate
            obj.Speed = s;
            if ~isempty(obj.C.Video) && obj.C.Video.hasVideo()
                obj.C.Video.setRate(s);
            end
        end

        function seekTo(obj, timeSec)
            % Request a seek to timeSec. The play loop will handle it.
            obj.seek_to = max(0, timeSec);
        end

        function t = getCurrentTime(obj)
            % Return current playback position in seconds
            t = obj.totalSamplesOut / max(obj.C.Setup.SampleRate, 1);
        end

        function d = getDuration(obj)
            if obj.UseMic
                d = 0;  % live input has no fixed duration
                return;
            end
            d = obj.MediaDuration;
            % Also check video duration (might be longer if extracted differently)
            if ~isempty(obj.C.Video) && obj.C.Video.hasVideo()
                vd = obj.C.Video.getDuration();
                if vd > 0, d = max(d, vd); end
            end
        end
    end

    % ═════════════════════════════════════════════════════
    %  Mic playback loop (private)
    % ═════════════════════════════════════════════════════
    methods (Access = private)

        function ok = tryOpenObjects(obj, filename)
            % Open the file with Cavern; accept it if it decodes and runs at
            % the processing sample rate. Formats Cavern cannot render
            % (AAC, TrueHD in this build, ...) fall back to the file reader.
            ok = false;
            [~,~,ext] = fileparts(filename);
            if ~any(strcmpi(ext, {'.ac3','.eac3','.ec3','.mp4','.m4a','.m4v','.mkv','.mka','.mov','.webm'}))
                return;
            end
            try
                cr = cavern_reader(filename);
                cr.read(obj.C.Setup.Block_size);             % some renderers fail on the first block
                if cr.failed || cr.fs ~= obj.C.Setup.SampleRate
                    if cr.fs ~= obj.C.Setup.SampleRate
                        warning('AudioEngine:objectsFs', 'Cavern stream at %g Hz, processing at %g Hz: using the file reader.', ...
                            cr.fs, obj.C.Setup.SampleRate);
                    end
                    delete(cr);
                    return;
                end
                cr.seek(0);                                  % back to the start
                obj.Objects = cr;
            catch err
                msg = strtok(err.message, newline);
                fprintf(2, 'Cavern: %s - using the file reader.\n', msg);
                if any(strcmp(err.identifier, {'cavern_reader:notBuilt','cavern_reader:protocol'}))
                    % set-up problem, not a format Cavern does not handle: say it out loud
                    warndlg(sprintf(['Object-based (Atmos) decoding is unavailable:\n\n%s\n\n' ...
                        'The file is played as plain channels.'], msg), 'Cavern', 'modal');
                end
                obj.Objects = [];
                return;
            end
            obj.C.Setup.N_in = obj.Objects.nObj;
            obj.C.Setup.SourceNames = arrayfun(@(k) sprintf('Obj %d', k), 1:obj.Objects.nObj, 'UniformOutput', false);
            obj.InputFs = obj.Objects.fs;
            obj.MediaDuration = obj.Objects.duration();
            obj.SRC = [];
            obj.C.InputDisplayName = sprintf('[Atmos/Cavern, %d objects] %s', obj.Objects.nObj, filename);
            ok = true;
        end

        function y = readObjects(obj, blk)
            % Next block of object samples; moves the scene's virtual sources
            if obj.Objects.atEnd()
                obj.Objects.seek(0);                             % loop from the start
                obj.totalSamplesOut = 0;
                if ~isempty(obj.C.Video) && obj.C.Video.hasVideo(), obj.C.Video.rewind(); end
            end
            [y, pos] = obj.Objects.read(blk);
            % object level for the display: fast attack, 0.5 s release
            lev = sqrt(mean(y.^2, 1))';
            rel = exp(-blk/obj.Objects.fs/0.5);
            if numel(obj.ObjectLevel) ~= numel(lev), obj.ObjectLevel = lev; end
            obj.ObjectLevel = max(lev, rel*obj.ObjectLevel + (1-rel)*lev);
            sc = obj.C.Scene;
            sc.ObjectRoomRadius = obj.C.Setup.Loudspeaker_setup.R + 1;  % = default source circle
            sc.apply_object_positions(pos, obj.ObjectLevel);
        end

        function playMic(obj)
            % Real-time mic input → render → output loop
            if ~(obj.C.HasDSP && obj.C.HasAudio && ~isempty(obj.MicInput))
                errordlg('Mic playback disabled: missing toolbox or no mic configured.', ...
                    'Cannot play','modal');
                return;
            end

            blk = obj.C.Setup.Block_size;
            fs  = obj.C.Setup.SampleRate;
            obj.initIO(fs, blk);

            % Configure output
            if ~isempty(obj.Out)
                try release(obj.Out); catch, end
                try obj.Out.SampleRate = fs; catch, end
                try obj.Out.BufferSize = blk; catch, end
            end

            % Setup mic reader to match
            try release(obj.MicInput); catch, end
            obj.MicInput.SampleRate = fs;
            obj.MicInput.SamplesPerFrame = blk;

            obj.stop_now  = false;
            obj.pause_now = false;
            obj.is_paused = false;
            obj.totalSamplesOut = 0;
            drawnowCounter = 0;

            nChExpected = obj.C.Setup.N_in;
            obj.Running = true;
            runningGuard = onCleanup(@() obj.setRunning(false)); %#ok<NASGU>

            fprintf('Mic playback started: %s\n', obj.MicDeviceName);

            while ~obj.stop_now

                % ── Handle pause ──
                if obj.pause_now && ~obj.is_paused
                    obj.is_paused = true;
                    obj.pause_now = false;
                end
                if obj.is_paused
                    % Still read from mic to keep the buffer flushed,
                    % but discard the data
                    try
                        obj.MicInput(); %#ok<VUNUS>
                    catch
                    end
                    drawnow limitrate
                    java.lang.Thread.sleep(20);
                    continue;
                end

                % ── Read from microphone ──
                inFrame = obj.MicInput();

                % Ensure correct number of channels
                nChGot = size(inFrame, 2);
                if nChGot > nChExpected
                    inFrame = inFrame(:, 1:nChExpected);
                elseif nChGot < nChExpected
                    inFrame = [inFrame, zeros(blk, nChExpected - nChGot)]; %#ok<AGROW>
                end

                % ── Master EQ ──
                if ~isempty(obj.MasterEQ) && isvalid(obj.MasterEQ)
                    inFrame = obj.MasterEQ.apply(inFrame);
                end
                obj.applyTrajectories();            % designed source paths
                obj.updateBeamOverlay(inFrame);

                % ── Render through sound scene ──
                outFrame = obj.C.Setup.MasterVolume * obj.C.Scene.render_sound_scene( ...
                    inFrame, obj.C.Setup.Binauralization, ...
                    obj.C.Setup.Downmixing_enabled, obj.C.Setup.Decorrelation);

                % ── Output ──
                obj.Out(outFrame);
                obj.totalSamplesOut = obj.totalSamplesOut + size(outFrame,1);

                % Let UI callbacks run
                drawnowCounter = drawnowCounter + 1;
                % every block: rendering takes ~1 ms now, and 'limitrate'
                % caps the redraws; every 8th block made dragging lag ~0.3 s
                if drawnowCounter >= 1
                    drawnow limitrate
                    drawnowCounter = 0;
                end
            end

            % Cleanup
            try release(obj.Out); catch, end
            % Don't release MicInput here — keep it configured for next play
            fprintf('Mic playback stopped.\n');
        end

        function wavPath = extractAudioWithFFmpeg(~, inputFile)
            % Extract audio track from any media file to a temp WAV
            % using ffmpeg. Returns '' if ffmpeg is not available.
            wavPath = '';

            % Find ffmpeg
            ffmpegExe = '';
            if ispc
                [st, res] = system('where ffmpeg 2>nul');
                if st == 0
                    lines = strsplit(strtrim(res), newline);
                    ffmpegExe = strtrim(lines{1});
                end
            else
                [st, res] = system('which ffmpeg 2>/dev/null');
                if st == 0
                    ffmpegExe = strtrim(res);
                end
            end

            % Also check common locations
            if isempty(ffmpegExe)
                commonPaths = { ...
                    'C:\ffmpeg\bin\ffmpeg.exe', ...
                    'C:\Program Files\ffmpeg\bin\ffmpeg.exe', ...
                    '/usr/local/bin/ffmpeg', ...
                    '/opt/homebrew/bin/ffmpeg'};
                for k = 1:numel(commonPaths)
                    if isfile(commonPaths{k})
                        ffmpegExe = commonPaths{k};
                        break;
                    end
                end
            end

            if isempty(ffmpegExe)
                return;  % ffmpeg not found
            end

            % Extract audio to temp WAV (PCM 16-bit, preserve sample rate & channels)
            wavPath = fullfile(tempdir, ['binauralizer_ffmpeg_' char(java.util.UUID.randomUUID) '.wav']);
            cmd = sprintf('"%s" -y -i "%s" -vn -acodec pcm_s16le "%s"', ...
                ffmpegExe, inputFile, wavPath);

            fprintf('Extracting audio with ffmpeg...\n');
            [status, cmdout] = system(cmd);

            if status ~= 0 || ~isfile(wavPath)
                fprintf('ffmpeg failed: %s\n', cmdout);
                wavPath = '';
            else
                fprintf('Audio extracted to: %s\n', wavPath);
            end
        end
    end
end

function y = ternary(cond, a, b); if cond, y=a; else, y=b; end; end
