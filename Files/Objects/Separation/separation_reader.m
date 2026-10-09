classdef separation_reader < handle
    %SEPARATION_READER  Object provider: media -> separator -> objects.
    %
    %   Same interface as cavern_reader, so AudioEngine plays it as
    %   object-based input:  [X, pos] = read(n), seek, seekAsync/seekDone,
    %   atEnd, duration, nObj, fs, samplesRead, failed.
    %   pos rows are normalised room coordinates (see separator_base).
    %   The media is read from a random-access WAV at the processing rate
    %   (the file itself, or an ffmpeg / audioread conversion).

    properties (SetAccess = private)
        file            % source media
        wav             % random-access WAV actually read
        sep             % separator_base
        fs
        nIn
        nObj
        total           % samples in the media
        samplesRead = 0 % media position of the next read
        failed = false
        Names = {}
    end

    properties (Access = private)
        tmp = ''        % converted WAV (cached for the session)
        buf = []        % read-ahead: samples buf0+1 .. buf0+size(buf,1)
        buf0 = 0
    end

    properties (Constant)
        ReadAhead = 0.5 % [s] per file access (audioread per block is slow)
    end

    methods
        function obj = separation_reader(file, sep, fs)
            obj.file = file; obj.sep = sep; obj.fs = fs;
            [obj.wav, obj.tmp] = separation_reader.as_wav(file, fs);
            ai = audioinfo(obj.wav);
            obj.nIn = ai.NumChannels;
            obj.total = ai.TotalSamples;
            sep.init(fs, obj.nIn);
            obj.nObj = sep.nOut();
            obj.Names = sep.Names;
        end

        function delete(~)
            % the converted WAV stays cached for this MATLAB session (as_wav)
        end

        function [Y, pos] = read(obj, n)
            a = obj.samplesRead + 1;
            b = min(obj.samplesRead + n, obj.total);
            if b >= a
                if a <= obj.buf0 || b > obj.buf0 + size(obj.buf,1)
                    e = min(obj.total, a - 1 + max(n, round(obj.ReadAhead*obj.fs)));
                    obj.buf = audioread(obj.wav, [a e]);
                    obj.buf0 = a - 1;
                end
                X = obj.buf(a - obj.buf0 : b - obj.buf0, :);
            else
                X = zeros(0, obj.nIn);
            end
            X = [X; zeros(n - size(X,1), obj.nIn)];
            obj.samplesRead = obj.samplesRead + n;
            Y = obj.sep.process(X);
            pos = obj.sep.Pos;
        end

        function seek(obj, sample)
            obj.samplesRead = max(0, min(round(sample), obj.total));
            obj.sep.reset();
        end

        function seekAsync(obj, sample), obj.seek(sample); end
        function done = seekDone(~), done = true; end

        function tf = atEnd(obj)
            tf = obj.samplesRead >= obj.total + obj.sep.Latency;
        end

        function p = positions(obj), p = obj.sep.Pos; end
        function t = duration(obj), t = obj.total/obj.fs; end
    end

    methods (Static)
        function [wav, tmp] = as_wav(file, fs)
            % WAV at fs that audioread can read by sample range. Conversions
            % are cached per (file, fs) for the session: a film is extracted
            % once, not for every analysis or re-open.
            persistent cache
            if isempty(cache), cache = containers.Map(); end
            key = sprintf('%s|%d', file, fs);
            if isKey(cache, key) && isfile(cache(key))
                wav = cache(key); tmp = wav; return;
            end
            [wav, tmp] = separation_reader.convert(file, fs);
            if ~isempty(tmp), cache(key) = tmp; end
        end

        function [wav, tmp] = convert(file, fs)
            tmp = '';
            if ~isfile(file) && ~isempty(which(file)), file = which(file); end   % relative to the MATLAB path (ffmpeg needs a real path)
            [~,~,ext] = fileparts(file);
            if strcmpi(ext, '.wav')
                try
                    ai = audioinfo(file);
                    if ai.SampleRate == fs, wav = file; return; end
                catch
                end
            end
            tmp = fullfile(tempdir, ['binauralizer_sep_' char(java.util.UUID.randomUUID) '.wav']);
            ff = separation_reader.ffmpeg();
            if ~isempty(ff)
                fprintf('Extracting audio with ffmpeg (%d Hz)...\n', fs);
                cmd = sprintf('"%s" -y -v error -i "%s" -vn -ar %d -acodec pcm_s16le "%s"', ff, file, fs, tmp);
                [st, out] = system(cmd);
                if st == 0 && isfile(tmp), wav = tmp; return; end
                fprintf(2, 'ffmpeg: %s\n', strtrim(out));
            end
            % no ffmpeg: whole file through audioread (short media only)
            [x, fs0] = audioread(file);
            if fs0 ~= fs
                t0 = (0:size(x,1)-1)'/fs0;
                x = interp1(t0, x, (0:floor(t0(end)*fs))'/fs, 'pchip');
            end
            audiowrite(tmp, x, fs, 'BitsPerSample', 24);
            wav = tmp;
        end

        function exe = ffmpeg()
            exe = '';
            if ispc
                [st, res] = system('where ffmpeg 2>nul');
            else
                [st, res] = system('which ffmpeg 2>/dev/null');
            end
            if st == 0
                lines = strsplit(strtrim(res), newline);
                exe = strtrim(lines{1});
            end
        end
    end
end
