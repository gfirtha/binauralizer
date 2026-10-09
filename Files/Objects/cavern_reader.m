classdef cavern_reader < handle
    % CAVERN_READER  Object-based audio (Dolby Atmos E-AC-3 JOC, ...) via Cavern.
    %
    % Cavern (C#) runs in a helper process, External_software/cavern_build/
    % CavernServe.exe, and streams the decoded objects over a localhost TCP
    % socket. (Hosting it in MATLAB needs .NET Core, which the desktop
    % session cannot use: it loads .NET Framework at start-up.) Every read
    % returns the samples of all rendered objects (beds included, as static
    % objects) and their positions at the end of the block:
    %
    %   cr = cavern_reader(file);
    %   [X, pos] = cr.read(2048);   % X [n x Nobj], pos [Nobj x 3]
    %
    % pos is normalised to the Atmos room: [right, front, up] with right and
    % front in -1..1 (listener at the centre) and up in 0..1.
    %
    % Cavern licence: free, but no part may be sold - research use only.

    properties (SetAccess = private)
        file
        fs              % sample rate [Hz]
        nObj            % number of objects (beds included)
        isLFE           % [Nobj x 1] logical (not reported by the stream: all false)
        len             % stream length [samples per object]
        samplesRead = 0
        failed = false  % decoding stopped (Cavern exception); later reads give silence
        lastError = ''
        corrupt         % [Nobj x 1] logical: objects muted after an out-of-range decode
    end

    properties (Access = private)
        proc            % java.lang.Process of CavernServe
        chan            % java.nio.channels.SocketChannel (blocking)
        lastPos         % [Nobj x 3]
        seeking = false % an asynchronous seek is pending
        seekBuf         % partial seek answer
    end

    methods
        function obj = cavern_reader(file)
            obj.file = char(file);
            exe = cavern_reader.serverPath();
            pb = java.lang.ProcessBuilder({exe, obj.file});
            pb.redirectErrorStream(true);
            obj.proc = pb.start();
            out = java.io.BufferedReader(java.io.InputStreamReader(obj.proc.getInputStream()));
            line = char(out.readLine());                 % "PORT n" or "ERROR msg"
            if ~startsWith(line, 'PORT ')
                obj.proc.destroy();
                error('cavern_reader:open', 'Cavern cannot open %s: %s', obj.file, ...
                    strtrim(erase(line, 'ERROR ')));
            end
            port = str2double(extractAfter(line, 'PORT '));
            % Java NIO, blocking: tcpclient's typed/polled reads cost ~10-18 ms per block
            obj.chan = java.nio.channels.SocketChannel.open(java.net.InetSocketAddress('127.0.0.1', port));
            obj.chan.socket().setTcpNoDelay(true);
            obj.chan.configureBlocking(false);          % recv polls, with a timeout
            hdr = typecast(obj.recv(16), 'int32');
            if hdr(1) ~= int32(hex2dec('4E525643'))
                error('cavern_reader:protocol', 'Unexpected answer from CavernServe.');
            end
            obj.nObj = double(hdr(3));
            obj.fs = double(hdr(4));
            obj.len = double(typecast(obj.recv(8), 'int64'));
            obj.isLFE = false(obj.nObj, 1);
            obj.corrupt = false(obj.nObj, 1);
            obj.lastPos = zeros(obj.nObj, 3);
        end

        function [X, pos] = read(obj, n)
            % Next n samples of every object and the object positions.
            obj.send([uint8('R'), typecast(int32(n), 'uint8')]);
            raw = obj.recv(4*(1 + n*obj.nObj + 3*obj.nObj));
            status = typecast(raw(1:4), 'int32');
            v = double(typecast(raw(5:end), 'single'));
            X = reshape(v(1:n*obj.nObj), obj.nObj, n)';
            pos = reshape(v(n*obj.nObj+1:end), 3, obj.nObj)';
            obj.lastPos = pos;
            if status == 2 && ~obj.failed
                L = double(typecast(obj.recv(4), 'int32'));
                obj.lastError = native2unicode(obj.recv(L), 'UTF-8');
                obj.failed = true;
                % e.g. Cavern's OAMD decoder on some streams: keep the
                % session alive, report once, play silence from here on
                warning('cavern_reader:decode', 'Cavern stopped decoding: %s', obj.lastError);
            end
            if status == 0
                obj.samplesRead = obj.samplesRead + n;
            end
            % Guard: Cavern's LFE decoding (AHT path unfinished) can return
            % samples ~+50 dB over full scale. Mute an object for good once
            % it exceeds +18 dBFS - no legitimate decode gets there.
            bad = ~obj.corrupt & max(abs(X), [], 1)' > 8;
            if any(bad)
                obj.corrupt = obj.corrupt | bad;
                warning('cavern_reader:corrupt', ...
                    'Object(s) %s returned |x| > 8 at %.2f s (decoder fault) - muted.', ...
                    mat2str(find(bad)'-1), obj.samplesRead/obj.fs);
            end
            X(:, obj.corrupt) = 0;
        end

        function seek(obj, sample)
            % Restart and decode up to the given sample (Cavern decodes
            % sequentially: costs ~1/4..1/12 of the skipped duration)
            obj.send([uint8('S'), typecast(int64(sample), 'uint8')]);
            r = obj.recv(12, 120);
            status = typecast(r(1:4), 'int32');
            obj.samplesRead = double(typecast(r(5:12), 'int64'));   % container seeks land on a block
            obj.failed = status == 2;
            obj.corrupt(:) = false;
        end

        function seekAsync(obj, sample)
            % Start a seek without waiting (container seeks of long files take
            % 10-20 s in Cavern); poll with seekDone() and keep the UI alive.
            obj.send([uint8('S'), typecast(int64(sample), 'uint8')]);
            obj.seekBuf = uint8([]);
            obj.seeking = true;
        end

        function done = seekDone(obj)
            % Non-blocking: true once the pending seek has completed.
            done = ~obj.seeking;
            if done, return; end
            if ~obj.proc.isAlive()
                error('cavern_reader:closed', 'CavernServe closed the connection.');
            end
            bb = java.nio.ByteBuffer.allocate(12 - numel(obj.seekBuf));
            obj.chan.read(bb);
            got = bb.position();
            if got > 0
                a = typecast(bb.array(), 'uint8')';
                obj.seekBuf = [obj.seekBuf, a(1:got)];
            end
            if numel(obj.seekBuf) == 12
                status = typecast(obj.seekBuf(1:4), 'int32');
                obj.samplesRead = double(typecast(obj.seekBuf(5:12), 'int64'));
                obj.failed = status == 2;
                obj.corrupt(:) = false;
                obj.seeking = false;
                done = true;
            end
        end

        function tf = atEnd(obj)
            tf = obj.failed || obj.samplesRead >= obj.len;
        end

        function pos = positions(obj)
            pos = obj.lastPos;
        end

        function t = duration(obj)
            t = obj.len / obj.fs;
        end

        function delete(obj)
            try, obj.send(uint8('Q')); catch, end
            try, obj.chan.close(); catch, end
            try, obj.proc.destroy(); catch, end
        end
    end

    methods (Access = private)
        function send(obj, bytes)
            bb = java.nio.ByteBuffer.wrap(typecast(uint8(bytes(:)'), 'int8'));
            while bb.hasRemaining()
                obj.chan.write(bb);
            end
        end

        function raw = recv(obj, n, timeout)
            % exactly n bytes; errors after `timeout` s (default 10) so that a
            % stuck or dead server cannot freeze MATLAB
            if nargin < 3, timeout = 10; end
            bb = java.nio.ByteBuffer.allocate(n);
            t0 = tic;
            while bb.hasRemaining()
                got = obj.chan.read(bb);
                if got < 0 || ~obj.proc.isAlive()
                    error('cavern_reader:closed', 'CavernServe closed the connection.');
                elseif got == 0
                    if toc(t0) > timeout
                        error('cavern_reader:timeout', 'CavernServe did not answer within %g s.', timeout);
                    end
                    java.lang.Thread.sleep(1);
                end
            end
            raw = typecast(bb.array(), 'uint8')';
        end
    end

    methods (Static)
        function exe = serverPath()
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
            exe = fullfile(root, 'External_software', 'cavern_build', 'CavernServe.exe');
            if ~isfile(exe)
                error('cavern_reader:notBuilt', ['CavernServe is not built. Run (needs the .NET 8 SDK):\n' ...
                    '  dotnet build External_software/CavernServe/CavernServe.csproj -c Release ' ...
                    '-o External_software/cavern_build']);
            end
        end
    end
end
