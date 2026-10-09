classdef passthrough_separator < separator_base
    %PASSTHROUGH_SEPARATOR  Reference: every input channel is an object at
    %its standard loudspeaker direction (stereo +-30, 5.1/7.1 ITU, LFE
    %unpositioned).

    methods
        function init(obj, fs, nIn)
            obj.fs = fs; obj.nIn = nIn;
            [az, obj.Names] = separator_base.channel_layout(nIn);
            obj.Pos = separator_base.dir2pos(az);
            obj.Pos(isnan(az), :) = 0;
            obj.Latency = 0;
        end

        function Y = process(~, X)
            Y = X;
        end

        function reset(~)
        end
    end
end
