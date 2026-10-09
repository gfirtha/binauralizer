classdef VideoService < handle
    % VIDEOSERVICE  Video playback via external VLC with RC control.
    %
    % Design:
    %   - Stop kills VLC entirely; Play launches a fresh VLC instance.
    %   - VLC keyboard IS enabled (space = pause toggle, f = fullscreen).
    %   - A periodic sync timer detects if VLC was paused/resumed externally
    %     and mirrors that to the AudioEngine.
    %   - The ONLY method that sends 'pause' to VLC is vlcToggle().
    %     It also flips VLCIsPlaying. External toggles (user pressing
    %     space in VLC) are detected by queryIsPlaying() and reconciled.

    properties
        C
        Filename  char = ''
        Duration  double = 0
        HasVLC    logical = false
        VLCPort   double = 0
        VLCSock
        FPS       double = 30
        VLCIsPlaying logical = false  % our tracking of VLC state
        Rate      double = 1          % playback rate (slow motion)
    end

    properties (Constant, Access = private)
        VLC_PATHS = { ...
            'C:\Program Files\VideoLAN\VLC\vlc.exe', ...
            'C:\Program Files (x86)\VideoLAN\VLC\vlc.exe', ...
            '/Applications/VLC.app/Contents/MacOS/VLC', ...
            '/usr/bin/vlc', ...
            '/snap/bin/vlc'};
    end

    methods
        function obj = VideoService(ctx)
            obj.C = ctx;
        end

        function tf = hasVideo(obj)
            tf = isscalar(obj) && isvalid(obj) && ~isempty(obj.Filename);
        end

        % ═══════════════════════════════════════════════
        %  Open / Close
        % ═══════════════════════════════════════════════
        function open(obj, filename)
            obj.killVLC();
            obj.Filename = filename;

            vlcExe = obj.findVLC();
            obj.HasVLC = ~isempty(vlcExe);

            if ~obj.HasVLC
                warndlg(['VLC media player is required for video playback.' newline newline ...
                    'Please install VLC from:' newline ...
                    'https://www.videolan.org/vlc/' newline newline ...
                    'Audio will still be processed.'], ...
                    'VLC not found', 'modal');
            end

            try
                vr = VideoReader(filename);
                obj.FPS = vr.FrameRate;
                obj.Duration = vr.Duration;
                delete(vr);
            catch
                obj.FPS = 30;
                obj.Duration = 0;
            end
            % Don't launch VLC yet — wait for Play button
        end

        function close(obj)
            obj.killVLC();
            obj.Filename = '';
            obj.FPS = 30;
            obj.Duration = 0;
            obj.VLCIsPlaying = false;
        end

        % ═══════════════════════════════════════════════
        %  Playback control
        % ═══════════════════════════════════════════════
        function playStart(obj)
            % Launch a fresh VLC and start playing
            if ~obj.HasVLC || isempty(obj.Filename), return; end
            obj.killVLC();           % kill any leftover
            obj.launchVLC(obj.Filename);  % VLC starts PLAYING
            % VLC is now playing from t=0 — no need to send any command
        end

        function playStop(obj)
            % Kill VLC entirely
            obj.killVLC();
        end

        function vlcPause(obj)
            if ~obj.HasVLC || isempty(obj.VLCSock), return; end
            if obj.VLCIsPlaying
                obj.vlcToggle();
            end
        end

        function vlcResume(obj)
            if ~obj.HasVLC || isempty(obj.VLCSock), return; end
            if ~obj.VLCIsPlaying
                obj.vlcToggle();
            end
        end

        function seekTo(obj, timeSec)
            if obj.HasVLC && ~isempty(obj.VLCSock)
                obj.vlcCmd(sprintf('seek %d', round(timeSec)));
            end
        end

        function rewind(obj)
            if obj.HasVLC && ~isempty(obj.VLCSock)
                obj.vlcCmd('seek 0');
                java.lang.Thread.sleep(50);
            end
        end

        function setRate(obj, r)
            obj.Rate = r;
            if obj.HasVLC && ~isempty(obj.VLCSock)
                obj.vlcCmd(sprintf('rate %g', r));
            end
        end

        function toggleFullscreen(obj)
            if obj.HasVLC && ~isempty(obj.VLCSock)
                obj.vlcCmd('fullscreen');
            end
        end

        function ensureWindow(~), end
        function showFirstFrame(~), end
        function onAudioSamplesProcessed(~, ~), end
        function d = getDuration(obj), d = obj.Duration; end

        % ═══════════════════════════════════════════════
        %  Sync: detect if user pressed Space in VLC window
        % ═══════════════════════════════════════════════
        function changed = syncFromVLC(obj)
            % Query VLC's actual state and update VLCIsPlaying.
            % Returns true if the state changed (user toggled in VLC).
            changed = false;
            if ~obj.HasVLC || isempty(obj.VLCSock), return; end
            actual = obj.queryIsPlaying();
            if actual ~= obj.VLCIsPlaying
                obj.VLCIsPlaying = actual;
                changed = true;
            end
        end
    end

    % ═══════════════════════════════════════════════════════
    %  Private
    % ═══════════════════════════════════════════════════════
    methods (Access = private)
        function vlcToggle(obj)
            % The ONLY place we send 'pause' and flip our state.
            obj.vlcCmd('pause');
            obj.VLCIsPlaying = ~obj.VLCIsPlaying;
        end

        function playing = queryIsPlaying(obj)
            % Ask VLC RC 'is_playing' — returns true/false.
            playing = obj.VLCIsPlaying;  % fallback on error
            if isempty(obj.VLCSock), return; end
            try
                % Flush stale data
                if obj.VLCSock.NumBytesAvailable > 0
                    read(obj.VLCSock);
                end
                writeline(obj.VLCSock, 'is_playing');
                java.lang.Thread.sleep(40);
                if obj.VLCSock.NumBytesAvailable > 0
                    raw = char(read(obj.VLCSock)');
                    % Response: "1\r\n> " or "0\r\n> "
                    % Find the first '0' or '1' character
                    idx1 = find(raw == '1', 1);
                    idx0 = find(raw == '0', 1);
                    if ~isempty(idx1) && (isempty(idx0) || idx1 < idx0)
                        playing = true;
                    elseif ~isempty(idx0)
                        playing = false;
                    end
                end
            catch
            end
        end

        function vlcExe = findVLC(~)
            vlcExe = '';
            paths = audioapp.services.VideoService.VLC_PATHS;
            for k = 1:numel(paths)
                if isfile(paths{k})
                    vlcExe = paths{k};
                    return;
                end
            end
            if ispc
                [st, ~] = system('where vlc 2>nul');
            else
                [st, ~] = system('which vlc 2>/dev/null');
            end
            if st == 0
                vlcExe = 'vlc';
            end
        end

        function launchVLC(obj, filename)
            vlcExe = obj.findVLC();
            if isempty(vlcExe), obj.HasVLC = false; return; end

            obj.VLCPort = 10000 + randi(50000);

            absFile = filename;
            if ~isfile(absFile), absFile = fullfile(pwd, filename); end
            if ~isfile(absFile)
                tmp = which(filename);
                if ~isempty(tmp), absFile = tmp; end
            end

            % VLC keyboard IS enabled — space = play/pause, f = fullscreen
            % --global-key-toggle-fullscreen "Ctrl+Return" for Ctrl+Enter
            cmd = sprintf(['"%s" --intf rc --rc-host 127.0.0.1:%d' ...
                ' --no-video-title-show --no-audio --repeat' ...
                ' --key-toggle-fullscreen "Ctrl+Return" --rate=%g' ...
                ' "%s"'], vlcExe, obj.VLCPort, obj.Rate, absFile);

            if ispc
                try
                    pb = java.lang.ProcessBuilder({'cmd', '/c', ['start /B /MIN "" ' cmd]});
                    pb.redirectErrorStream(true);
                    pb.start();
                catch
                    system(sprintf('start /B /MIN "" %s', cmd));
                end
            else
                system([cmd ' >/dev/null 2>&1 &']);
            end

            % Connect to RC
            obj.VLCSock = [];
            for attempt = 1:40
                java.lang.Thread.sleep(100);
                try
                    obj.VLCSock = tcpclient('127.0.0.1', obj.VLCPort, ...
                        'Timeout', 1, 'ConnectTimeout', 1);
                    break;
                catch
                end
            end

            if isempty(obj.VLCSock)
                warning('VideoService:vlcConnect', 'Could not connect to VLC.');
                obj.HasVLC = false;
                return;
            end

            % Drain greeting
            java.lang.Thread.sleep(500);
            try read(obj.VLCSock); catch, end

            % VLC just started with a file → it's playing
            obj.VLCIsPlaying = true;
        end

        function vlcCmd(obj, cmd)
            if isempty(obj.VLCSock), return; end
            try
                if obj.VLCSock.NumBytesAvailable > 0
                    read(obj.VLCSock);
                end
                writeline(obj.VLCSock, cmd);
                java.lang.Thread.sleep(10);
            catch
                obj.VLCSock = [];
            end
        end

        function killVLC(obj)
            if ~isempty(obj.VLCSock)
                try obj.vlcCmd('quit'); java.lang.Thread.sleep(200); catch, end
                try delete(obj.VLCSock); catch, end
                obj.VLCSock = [];
            end
            obj.VLCIsPlaying = false;
        end
    end
end
