classdef mbr_room_separator < separator_base
    %MBR_ROOM_SEPARATOR  Parametric reproduction of a measured room
    %(mbr_decompose_rir) as Binauralizer objects: the media, downmixed to
    %mono, is convolved with
    %   - the local impulse response of every reliable reflection
    %     (beamformed + time-focused), at its measured direction, and
    %   - the SDM impulse responses of the diffuse residual, at nDiff fixed
    %     directions.
    %   Used through the separator path of the AudioEngine (make_separator
    %   method 'measured_room'): its objects become the virtual sources.
    %   Objects sit on the unit sphere (direction only, the scene maps
    %   |Pos| = 1 to its object room radius); delays and levels are in the
    %   IRs (t = 0: direct sound), so no extra latency.

    properties
        Room                % struct of mbr_decompose_rir
        Gain = 1            % overall level
        DiffuseGain = 1     % level of the SDM residual objects
        UseDiffuse = true
    end

    properties (Access = private)
        IR = []
        conv = {}
        blk = 0
    end

    methods
        function obj = mbr_room_separator(room)
            if nargin > 0, obj.Room = room; end
        end

        function init(obj, fs, nIn)
            obj.fs = fs; obj.nIn = nIn;
            D = obj.Room;
            if D.fs ~= fs
                [pq, qq] = rat(fs/D.fs);
                ri = @(x) resample(x, pq, qq);
            else
                ri = @(x) x;
            end
            r = D.refl;
            obj.IR = ri(D.objIR);
            obj.Pos = separator_base.dir2pos([r.az], [r.el]);
            obj.Names = [{'direct'}, arrayfun(@(k) sprintf('refl %d (%.1f ms)', k-1, r(k).t*1000), ...
                2:numel(r), 'UniformOutput', false)];
            if obj.UseDiffuse && ~isempty(D.diffIR)
                obj.IR = [obj.IR, obj.DiffuseGain*ri(D.diffIR)];
                obj.Pos = [obj.Pos; separator_base.dir2pos(D.diffDirs(:,1), D.diffDirs(:,2))];
                obj.Names = [obj.Names, arrayfun(@(m) sprintf('diffuse %d', m), 1:size(D.diffIR,2), 'UniformOutput', false)];
            end
            % level: direct-sound IR peak -> Gain
            obj.IR = obj.Gain*obj.IR/max(abs(obj.IR(:,1)) + eps);
            obj.Latency = 0;
            obj.conv = {}; obj.blk = 0;
        end

        function reset(obj)
            obj.conv = {}; obj.blk = 0;              % rebuilt (state cleared) on the next block
        end

        function Y = process(obj, X)
            n = size(X, 1);
            if n ~= obj.blk || isempty(obj.conv)
                obj.blk = n;
                obj.conv = cell(1, size(obj.IR, 2));
                for k = 1:size(obj.IR, 2)
                    obj.conv{k} = OLS_convolver(obj.IR(:,k), n);
                end
            end
            x = mean(X, 2);
            Y = zeros(n, size(obj.IR, 2));
            for k = 1:size(obj.IR, 2)
                Y(:,k) = obj.conv{k}.convolve(x);
            end
        end
    end
end
