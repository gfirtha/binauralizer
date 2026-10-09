classdef spatial_audio_out < handle
    %SPATIAL_AUDIO_OUT  Output through the Windows spatial audio API.
    %
    %   Drop-in for audioDeviceWriter in AudioEngine: out(Y) plays a block
    %   (frames x channels, WAVE channel order), release(out) stops,
    %   out.info gives MaximumOutputChannels. The channels are sent to
    %   SpatialAudioBridge.exe (External_software/spatial_bridge_build) over
    %   a localhost socket; it plays them as static objects of the
    %   endpoint's active spatial format, e.g. "Dolby Atmos for home
    %   theater" -> Dolby MAT over HDMI, height channels included.
    %   Layouts by channel count: 2, 6 (5.1), 8 (7.1), 10 (5.1.4), 12 (7.1.4).
    %   The bridge acknowledges every block once its buffer has room again,
    %   which paces the caller in real time like a sound card.

    properties
        Device = ''         % endpoint name (substring)
        SampleRate = 48000
        BufferSize = 2048
    end

    properties (SetAccess = private)
        nCh = 0
        Format = ''         % active spatial format reported by the bridge
    end

    properties (Access = private)
        proc = []
        chan = []
    end

    methods
        function obj = spatial_audio_out(device)
            if nargin > 0, obj.Device = device; end
        end

        function delete(obj)
            obj.release();
        end

        function varargout = subsref(obj, s)
            % out(Y) plays a block, like a System object
            if strcmp(s(1).type, '()') && numel(s) == 1
                obj.write(s(1).subs{1});
                varargout = {};
            else
                [varargout{1:nargout}] = builtin('subsref', obj, s);
            end
        end

        function write(obj, Y)
            if isempty(obj.chan) || size(Y,2) ~= obj.nCh
                obj.release();
                obj.open(size(Y,2), size(Y,1));
            end
            n = size(Y,1);
            if n ~= obj.BufferSize                      % the bridge works in fixed blocks
                Y = [Y; zeros(max(0, obj.BufferSize - n), obj.nCh)];
                Y = Y(1:obj.BufferSize, :);
            end
            try
                obj.send(typecast(single(reshape(Y.', [], 1)), 'uint8'));
                obj.recv(1, 5);                          % ack: paces in real time
            catch ME
                obj.release();                           % never leave a bridge holding the endpoint
                rethrow(ME);
            end
        end

        function release(obj)
            if ~isempty(obj.chan)
                try, obj.chan.close(); catch, end
                obj.chan = [];
            end
            if ~isempty(obj.proc)
                try
                    if ~obj.proc.waitFor(2, java.util.concurrent.TimeUnit.SECONDS), obj.proc.destroyForcibly(); end
                catch
                    try, obj.proc.destroyForcibly(); catch, end
                end
                obj.proc = [];
            end
            obj.nCh = 0;
        end

        function s = info(obj)
            s = struct('MaximumOutputChannels', 12, 'DeviceName', obj.Device, 'Format', obj.Format);
        end
    end

    methods (Access = private)
        function open(obj, nCh, blk)
            obj.BufferSize = blk;
            pb = java.lang.ProcessBuilder({spatial_audio_out.bridgePath(), obj.Device});
            pb.redirectErrorStream(true);
            obj.proc = pb.start();
            rd = java.io.BufferedReader(java.io.InputStreamReader(obj.proc.getInputStream()));
            port = [];
            while isempty(port)
                line = rd.readLine();
                if isempty(line)
                    error('spatial_audio_out:open', 'SpatialAudioBridge exited without answering.');
                end
                line = char(line);
                if startsWith(line, 'PORT ')
                    port = str2double(extractAfter(line, 'PORT '));
                elseif startsWith(line, 'INFO ')
                    obj.Format = extractAfter(line, 'format ');
                elseif startsWith(line, 'ERROR ')
                    obj.proc.destroy(); obj.proc = [];
                    error('spatial_audio_out:open', '%s', extractAfter(line, 'ERROR '));
                end
            end
            obj.chan = java.nio.channels.SocketChannel.open(java.net.InetSocketAddress('127.0.0.1', port));
            obj.chan.socket().setTcpNoDelay(true);
            obj.chan.configureBlocking(false);
            obj.send(typecast(int32([hex2dec('42415053'), obj.SampleRate, nCh, blk]), 'uint8'));
            st = typecast(obj.recv(4, 10), 'int32');
            if st ~= 0
                len = typecast(obj.recv(4, 2), 'int32');
                msg = native2unicode(obj.recv(double(len), 2), 'UTF-8');
                obj.release();
                error('spatial_audio_out:stream', 'Spatial audio stream: %s', msg);
            end
            obj.nCh = nCh;
        end

        function send(obj, bytes)
            bb = java.nio.ByteBuffer.wrap(typecast(uint8(bytes(:)'), 'int8'));
            while bb.hasRemaining()
                obj.chan.write(bb);
            end
        end

        function raw = recv(obj, n, timeout)
            bb = java.nio.ByteBuffer.allocate(n);
            t0 = tic;
            while bb.hasRemaining()
                got = obj.chan.read(bb);
                if got < 0 || ~obj.proc.isAlive()
                    error('spatial_audio_out:closed', 'SpatialAudioBridge closed the connection.');
                elseif got == 0
                    if toc(t0) > timeout
                        error('spatial_audio_out:timeout', 'SpatialAudioBridge did not answer within %g s.', timeout);
                    end
                    java.lang.Thread.sleep(1);
                end
            end
            raw = typecast(bb.array(), 'uint8')';
        end
    end

    methods (Static)
        function exe = bridgePath()
            root = fileparts(fileparts(fileparts(mfilename('fullpath'))));     % binauralizer/
            exe = fullfile(root, 'External_software', 'spatial_bridge_build', 'SpatialAudioBridge.exe');
            if ~isfile(exe)
                error('spatial_audio_out:notBuilt', ['SpatialAudioBridge is not built: dotnet build ' ...
                    'External_software/SpatialAudioBridge -c Release -o External_software/spatial_bridge_build']);
            end
        end

        function list = devices()
            % Endpoints with an active spatial format: {name, format}
            list = cell(0, 2);
            try
                [st, out] = system(['"' spatial_audio_out.bridgePath() '" --list']);
            catch
                return;
            end
            if st ~= 0, return; end
            for ln = strsplit(strtrim(out), newline)
                p = strsplit(strtrim(ln{1}), sprintf('\t'));
                if numel(p) == 2 && ~strcmp(p{2}, '(off)'), list(end+1, :) = p; end %#ok<AGROW>
            end
        end
    end
end
