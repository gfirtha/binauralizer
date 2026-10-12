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
    %   - Virtual (3D / first-person scene views): no VLC, the picture goes
    %     onto the scene's virtual TV (VirtualTV: ffmpeg frames on the
    %     audio clock). setVirtual switches over mid-film; the film is
    %     decoded by one of the two only.

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
        Virtual   logical = false     % picture on the virtual TV instead of VLC
        TV                            % audioapp.services.VirtualTV
        Playing   logical = false     % between playStart and playStop
        AVOffset  double = 0          % s, + = picture later (fine lip-sync trim)
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
            if isempty(obj.TV), obj.TV = audioapp.services.VirtualTV(); end
            obj.TV.open(filename);
            obj.TV.ShowFcn = @(img) obj.gui().set_tv_image(img);
            if ~isempty(obj.gui())
                obj.gui().show_tv(true);
                obj.Virtual = any(strcmp(obj.gui().ViewMode, {'3d', 'listener'}));
            end

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
            if ~isempty(obj.TV), obj.TV.close(); end
            if ~isempty(obj.gui()), obj.gui().show_tv(false); end
            obj.Filename = '';
            obj.FPS = 30;
            obj.Duration = 0;
            obj.VLCIsPlaying = false;
        end

        % ═══════════════════════════════════════════════
        %  Playback control
        % ═══════════════════════════════════════════════
        function playStart(obj)
            % Launch a fresh VLC and start playing (virtual: the TV's decoder)
            if isempty(obj.Filename), return; end
            obj.Playing = true;
            if ~isempty(obj.gui()), obj.gui().set_tv_locked(true); end
            if obj.Virtual
                obj.TV.start(obj.heardTime());
                return;
            end
            if ~obj.HasVLC, return; end
            obj.killVLC();           % kill any leftover
            obj.launchVLC(obj.Filename);  % VLC starts PLAYING
            % VLC is now playing from t=0 — no need to send any command
        end

        function playStop(obj)
            % Kill VLC entirely
            obj.killVLC();
            obj.Playing = false;
            if ~isempty(obj.TV), obj.TV.stop(); end
            if ~isempty(obj.gui()), obj.gui().set_tv_locked(false); end
        end

        function setVirtual(obj, tf)
            % picture on the virtual TV (3D / first-person views) or in VLC;
            % mid-film the other one takes over at the current time
            if tf == obj.Virtual, return; end
            obj.Virtual = tf;
            if ~obj.Playing || isempty(obj.Filename), return; end
            t = obj.heardTime();
            paused = ~isempty(obj.C.Engine) && obj.C.Engine.is_paused;
            if tf
                obj.killVLC();
                if paused
                    obj.gui().set_tv_image(obj.TV.still(t));
                else
                    obj.TV.start(t);
                end
            else
                obj.TV.stop();
                if obj.HasVLC
                    obj.launchVLC(obj.Filename);
                    obj.seekTo(t);
                    if paused, obj.vlcPause(); end
                end
            end
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
            if obj.Virtual
                if obj.Playing, obj.TV.start(timeSec); end
                return;
            end
            if obj.HasVLC && ~isempty(obj.VLCSock)
                obj.vlcCmd(sprintf('seek %d', round(timeSec)));
            end
        end

        function rewind(obj)
            if obj.Virtual
                if obj.Playing, obj.TV.start(0); end
                return;
            end
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

        function arrangeScene(obj)
            % Home-cinema layout with a film: the TV stands right in front of
            % the centre loudspeaker (the array's loudspeaker closest to the
            % listener's facing direction; nothing moves, the rendering
            % stays as it is); a 65" screen. Without loudspeakers: in front
            % of the listener at R.
            g = obj.gui(); S = obj.C.Scene;
            if isempty(g) || isempty(S) || ~obj.hasVideo(), return; end
            rec = S.receiver.position(1:2); ori = S.receiver.orientation(1:2);
            f = ori/max(norm(ori), eps); side = [f(2), -f(1)];   % front, right
            L = S.loudspeaker_array;
            P = zeros(0, 2);
            for n = 1:numel(L), P(end+1, :) = L{n}.position(1:2); end %#ok<AGROW>
            if isempty(P)
                p = rec + obj.C.Setup.Loudspeaker_setup.R*f;
            else
                d = P - rec;
                az = atan2d(d*side', d*f');                  % 0 = front, + = right
                [~, c] = min(abs(az));
                p = P(c, :);
            end
            diag = 65;
            u = (rec - p)/max(norm(rec - p), eps);           % from the centre towards the listener
            g.draw_tv(p + 0.25*u, atan2d(u(2), u(1)), diag); % just in front of the centre loudspeaker
            g.show_tv(true);
        end

        function showFirstFrame(obj)
            % TV picture while not playing: a frame from the first minutes
            if isempty(obj.TV) || isempty(obj.gui()), return; end
            img = obj.TV.still(min(10, 0.1*obj.TV.Duration));
            if ~isempty(img), obj.gui().set_tv_image(img); end
        end

        function onAudioSamplesProcessed(obj, ~)
            % after each audio block of the play loop: the virtual TV's frame
            if ~obj.Virtual || ~obj.Playing || isempty(obj.TV), return; end
            t = obj.heardTime() + obj.AVOffset;
            if isempty(obj.TV.Proc)
                % decoder ended (paused switch-over, end of the video stream)
                recent = ~isempty(obj.TV.Started) && toc(obj.TV.Started) < 1;   % no respawn storm
                if ~recent && (obj.TV.Duration == 0 || t < obj.TV.Duration - 0.5), obj.TV.start(t); end
                return;
            end
            obj.TV.update(t);                  % shown by the loop's drawnow limitrate (<= 20 Hz)
        end
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
        function g = gui(obj)
            g = [];
            try g = obj.C.SceneGUI; catch, end
            if ~isempty(g) && ~isvalid(g), g = []; end
        end

        function t = heardTime(obj)
            % media time the listener hears now: output count minus the
            % output latency (~ two blocks in the device queue)
            t = 0;
            E = obj.C.Engine;
            if isempty(E), return; end
            fs = max(obj.C.Setup.SampleRate, 1);
            t = max(0, (E.totalSamplesOut - 2*obj.C.Setup.Block_size*E.Speed)/fs);
        end

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
