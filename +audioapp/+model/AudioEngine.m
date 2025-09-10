classdef AudioEngine < handle
    properties
        C                       % AppContext
        Input                   % dsp.AudioFileReader (or [])
        SRC                     % dsp.SampleRateConverter (or [])
        RB                      % dsp.AsyncBuffer
        Out                     % audioDeviceWriter (or [])
        stop_now logical = false

        MasterEQ   audioapp.dsp.MasterEQ
    end

    methods
        function obj = AudioEngine(ctx)
            obj.C = ctx;
            obj.ensureOutput();
            if obj.C.HasDSP
                obj.RB = dsp.AsyncBuffer(10*obj.C.Setup.Block_size);
            else
                obj.RB = [];
            end
            obj.MasterEQ = audioapp.dsp.MasterEQ(ctx);  % disabled by default

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
                obj.RB = dsp.AsyncBuffer(10*blk);
            else
                obj.RB = [];
            end

            if ~isempty(obj.MasterEQ) && isvalid(obj.MasterEQ)
                obj.MasterEQ.setFs(fs_proc);
            end

        end

        function refreshOutput(obj)
            obj.ensureOutput();
        end

        % ===== Unified media loader =======================================
        function kind = loadMedia(obj, filename)
            % kind = 'audio' or 'video'
            obj.C.InputDisplayName = filename;

            % Try to open video (if it fails, we treat as audio-only)
            hasVid = false;
            try
                obj.C.ensureVideo();
                obj.C.Video.open(filename);
                hasVid = true;
            catch
                % No usable video stream -> ensure video is closed
                if ~isempty(obj.C.Video), obj.C.Video.close(); end
            end
            kind = ternary(hasVid,'video','audio');

            % Audio reader (from file; video files with embedded audio work too)
            if obj.C.HasDSP
                fs_in = [];
                try
                    % IMPORTANT: PlayCount = 1 so *we* control looping
                    obj.Input = dsp.AudioFileReader( ...
                        filename, ...
                        'SamplesPerFrame', obj.C.Setup.Block_size, ...
                        'PlayCount', 1);
                    fs_in = obj.Input.SampleRate;
                    obj.C.Setup.N_in = obj.Input.NumChannels;
                catch
                    % Fallback: decode entire audio, write temp WAV
                    try
                        [y,Fs] = audioread(filename);
                        tmp = fullfile(tempdir, ['binauralizer_audio_' char(java.util.UUID.randomUUID) '.wav']);
                        audiowrite(tmp, y, Fs);
                        obj.Input = dsp.AudioFileReader(tmp, ...
                            'SamplesPerFrame', obj.C.Setup.Block_size, 'PlayCount', 1);
                        fs_in = obj.Input.SampleRate;
                        obj.C.Setup.N_in = obj.Input.NumChannels;
                    catch ME
                        obj.Input = [];
                        error('Audio open failed: %s', ME.message);
                    end
                end
                obj.initSRC(fs_in, obj.C.Setup.SampleRate);
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
            if ~(obj.C.HasDSP && obj.C.HasAudio && ~isempty(obj.Input))
                errordlg('Playback disabled: missing DSP/Audio toolbox or no input.','Cannot play','modal');
                return;
            end

            obj.initIO(obj.C.Setup.SampleRate, obj.C.Setup.Block_size);
            blk = obj.C.Setup.Block_size;

            % Reset both at start for a clean sync point
            reset(obj.Input);
            if obj.C.Video.hasVideo()
                obj.C.Video.rewind();
                obj.C.Video.ensureWindow();
                obj.C.Video.showFirstFrame();
            end

            obj.stop_now = false;
            while ~obj.stop_now
                if isDone(obj.Input)
                    reset(obj.Input);
                    if obj.C.Video.hasVideo(), obj.C.Video.rewind(); obj.C.Video.showFirstFrame(); end
                end

                x = obj.Input();
                if ~isempty(obj.SRC), y = obj.SRC(x); else, y = x; end
                write(obj.RB, y);

                while obj.RB.NumUnreadSamples >= blk && ~obj.stop_now
                    inFrame = read(obj.RB, blk);
                    outFrame = obj.C.Setup.MasterVolume*obj.C.Scene.render_sound_scene( ...
                        inFrame, obj.C.Setup.Binauralization, obj.C.Setup.Downmixing_enabled, obj.C.Setup.Decorrelation);
                    % Apply master EQ (no-op if disabled)
                    if ~isempty(obj.MasterEQ) && isvalid(obj.MasterEQ)
                        outFrame = obj.MasterEQ.apply(outFrame);
                    end

                    obj.Out(outFrame);
                    if obj.C.Video.hasVideo()
                        obj.C.Video.onAudioSamplesProcessed(size(outFrame,1));
                    end
                end

                % <<< Important: let UI callbacks (Stop) run
                drawnow limitrate
            end


            % Optional tail on user stop: drain partial block
            rem = obj.RB.NumUnreadSamples;
            if rem > 0
                pad = zeros(blk-rem, size(y,2), 'like', y);
                last = [read(obj.RB, rem); pad];
                outFrame = obj.C.Scene.render_sound_scene( ...
                    last, obj.C.Setup.Binauralization, obj.C.Setup.Downmixing_enabled, obj.C.Setup.Decorrelation);
                obj.Out(outFrame);
                if obj.C.Video.hasVideo()
                    obj.C.Video.onAudioSamplesProcessed(size(outFrame,1));
                end
            end

            try, release(obj.Out);   catch, end
            try, release(obj.Input); catch, end
        end

        function stop(obj)
            obj.stop_now = true;
        end
    end
end

function y = ternary(cond, a, b); if cond, y=a; else, y=b; end; end
