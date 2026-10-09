classdef Varispeed < handle
    % VARISPEED  Tape-style playback speed change of a block stream: the
    % input is read with step Rate (< 1 = slower, pitch drops with it) by
    % 4-point Hermite (Catmull-Rom) interpolation. Stateful across blocks,
    % so the output is continuous; its length per block varies
    % (about numel(x)/Rate). Rate = 1 at an integer read position passes
    % the samples through unchanged.

    properties (Access = private)
        buf = []        % unread input, with one sample of history in front
        pos = 2         % read position in buf (1-based, fractional)
    end

    methods
        function y = process(obj, x, rate)
            if isempty(obj.buf) || size(obj.buf, 2) ~= size(x, 2)
                obj.buf = zeros(1, size(x, 2), 'like', x);
                obj.pos = 2;
            end
            if rate == 1, obj.pos = round(obj.pos); end % back to 1x: exact pass-through
            obj.buf = [obj.buf; x];
            N = size(obj.buf, 1);
            last = N - 2;                               % needs i+2 <= N
            if obj.pos > last
                y = zeros(0, size(x, 2), 'like', x);
                return;
            end
            p = (obj.pos:rate:last)';
            i = floor(p); t = p - i;
            xm1 = obj.buf(i-1,:); x0 = obj.buf(i,:); x1 = obj.buf(i+1,:); x2 = obj.buf(i+2,:);
            c1 = 0.5*(x1 - xm1);
            c2 = xm1 - 2.5*x0 + 2*x1 - 0.5*x2;
            c3 = 0.5*(x2 - xm1) + 1.5*(x0 - x1);
            y = ((c3.*t + c2).*t + c1).*t + x0;
            % keep from the sample before the next read position on
            obj.pos = p(end) + rate;
            k = floor(obj.pos) - 2;
            if k > 0
                obj.buf = obj.buf(k+1:end, :);
                obj.pos = obj.pos - k;
            end
        end

        function tf = isActive(obj)
            % still holding samples from a slowed-down part
            tf = size(obj.buf, 1) > 1;
        end

        function reset(obj)
            obj.buf = [];
            obj.pos = 2;
        end
    end
end
