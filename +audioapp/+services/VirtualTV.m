classdef VirtualTV < handle
    % VIRTUALTV  Video frames for the virtual TV of the 3D / first-person
    % scene views: an ffmpeg process decodes the film scaled down to a
    % texture (rawvideo rgb24 on its stdout, any container ffmpeg reads,
    % mkv included); update(t) shows the frame due at media time t (the
    % time the listener hears), so the picture follows the audio clock:
    % lip sync, pause (no updates), slow motion and seeks for free.
    %
    % Frames are read only when ffmpeg already writes one (no waiting for
    % the process start or a slow decode in the audio loop); late frames
    % are dropped; only a jump of the clock (backwards, or > 5 s ahead)
    % restarts the decoder at the right time (loop to the start).

    properties
        Filename = ''
        FfmpegExe = ''
        W = 480                 % texture size (upload + redraw cost on the audio thread)
        H = 270
        FPS = 25                % output frame rate (source, halved above 30)
        Duration = 0
        ShowFcn = []            % @(img) put the frame on the TV
        Lead = 0.015            % s: the frame goes up this much early (drawnow)
    end

    properties (SetAccess = private)
        Proc = []               % java.lang.Process
        Stream = []
        T0 = 0                  % media time of frame 0 of the running decoder
        K = 0                   % frames read since T0
        Started = []            % tic of the last decoder start
    end

    methods
        function ok = open(obj, filename)
            % probe the film (size, frame rate); false without ffmpeg
            obj.close();
            obj.Filename = filename;
            obj.Duration = 0;                             % unknown unless probed
            obj.FfmpegExe = audioapp.services.VirtualTV.findFfmpeg();
            ok = ~isempty(obj.FfmpegExe);
            if ~ok, return; end
            [d, ~, e] = fileparts(obj.FfmpegExe);
            probe = fullfile(d, ['ffprobe' e]);
            [st, out] = system(sprintf(['"%s" -v error -select_streams v:0 -show_entries ' ...
                'stream=width,height,avg_frame_rate,r_frame_rate:format=duration -of default=nw=1 "%s"'], probe, filename));
            w = 16; h = 9; fps = 25;
            if st == 0
                v = @(key) regexp(out, [key '=([^\r\n]+)'], 'tokens', 'once');
                t = v('width');  if ~isempty(t), w = str2double(t{1}); end
                t = v('height'); if ~isempty(t), h = str2double(t{1}); end
                % avg_frame_rate: r_frame_rate of an mkv is often a
                % timebase-like 1000/1, which the fps filter would turn
                % into hundreds of duplicate frames per second
                for key = {'r_frame_rate', 'avg_frame_rate'}
                    t = v(key{1});
                    if ~isempty(t), f = str2num(t{1}); if isscalar(f) && f >= 5 && f <= 120, fps = f; end, end %#ok<ST2NM>
                end
                t = v('duration'); if ~isempty(t), obj.Duration = str2double(t{1}); end
            end
            while fps > 30, fps = fps/2; end           % 50/60 fps films: every second frame
            obj.FPS = fps;
            obj.H = 2*round(obj.W*h/w/2);
            fprintf('Virtual TV: %dx%d -> %dx%d at %g fps, %.1f s (probe status %d)\n', w, h, obj.W, obj.H, fps, obj.Duration, st);
        end

        function close(obj)
            obj.stop();
            obj.Filename = '';
        end

        function start(obj, t)
            % (re)start the decoder at media time t
            obj.stop();
            if isempty(obj.FfmpegExe) || isempty(obj.Filename), return; end
            t = max(0, t);
            if obj.Duration > 0 && t >= obj.Duration, return; end   % past the end of the picture
            args = {obj.FfmpegExe, '-v', 'error', '-ss', sprintf('%.3f', t), ...
                '-i', obj.Filename, '-an', '-sn', '-dn', ...
                '-vf', sprintf('fps=%g,scale=%d:%d', obj.FPS, obj.W, obj.H), ...
                '-f', 'rawvideo', '-pix_fmt', 'rgb24', '-'};
            pb = java.lang.ProcessBuilder(args);
            pb.redirectError(java.io.File(audioapp.services.VirtualTV.nullDevice()));
            obj.Proc = pb.start();
            obj.Stream = obj.Proc.getInputStream();
            obj.T0 = t; obj.K = 0; obj.Started = tic;
        end

        function stop(obj)
            if ~isempty(obj.Proc)
                try obj.Stream.close(); catch, end
                try obj.Proc.destroy(); catch, end
            end
            obj.Proc = []; obj.Stream = [];
        end

        function shown = update(obj, t)
            % show the frame due at media time t; true if a new one went up
            shown = false;
            if isempty(obj.Proc), return; end
            tNext = obj.T0 + obj.K/obj.FPS;               % pts of the next frame in the pipe
            if t < obj.T0 - 0.25 || t > tNext + 5
                fprintf('Virtual TV: decoder restart at %.2f s (picture at %.2f s)\n', t, tNext);
                obj.start(t);                             % a jump: loop to the start, missed seek
                return;
            end
            img = [];
            for n = 1:12                                  % behind: the late frames are dropped
                if obj.T0 + obj.K/obj.FPS > t + obj.Lead, break; end
                try
                    if obj.Stream.available() <= 0, break; end   % not decoded yet: no waiting
                    b = org.apache.commons.io.IOUtils.toByteArray(obj.Stream, obj.W*obj.H*3);
                catch
                    obj.stop(); break;                    % end of the stream
                end
                obj.K = obj.K + 1;
                img = b;
            end
            if isempty(img), return; end
            img = permute(reshape(typecast(img, 'uint8'), 3, obj.W, obj.H), [3 2 1]);
            if ~isempty(obj.ShowFcn), obj.ShowFcn(img); end
            shown = true;
        end

        function img = still(obj, t)
            % one frame at time t (TV picture while not playing)
            img = [];
            if isempty(obj.FfmpegExe) || isempty(obj.Filename), return; end
            obj.start(t);
            t0 = tic;
            while toc(t0) < 5
                try
                    if obj.Stream.available() > 0
                        b = org.apache.commons.io.IOUtils.toByteArray(obj.Stream, obj.W*obj.H*3);
                        img = permute(reshape(typecast(b, 'uint8'), 3, obj.W, obj.H), [3 2 1]);
                        break;
                    end
                catch
                    break;
                end
                java.lang.Thread.sleep(20);
            end
            obj.stop();
        end

        function delete(obj)
            obj.stop();
        end
    end

    methods (Static)
        function exe = findFfmpeg()
            exe = '';
            if ispc
                [st, res] = system('where ffmpeg 2>nul');
            else
                [st, res] = system('which ffmpeg 2>/dev/null');
            end
            if st == 0
                lines = strsplit(strtrim(res), newline);
                exe = strtrim(lines{1});
                return;
            end
            cand = {'C:\ffmpeg\bin\ffmpeg.exe', 'C:\Program Files\ffmpeg\bin\ffmpeg.exe', ...
                '/usr/local/bin/ffmpeg', '/opt/homebrew/bin/ffmpeg'};
            for k = 1:numel(cand)
                if isfile(cand{k}), exe = cand{k}; return; end
            end
        end

        function d = nullDevice()
            if ispc, d = 'NUL'; else, d = '/dev/null'; end
        end
    end
end
