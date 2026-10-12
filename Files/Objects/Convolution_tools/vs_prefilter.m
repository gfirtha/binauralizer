classdef vs_prefilter < handle
    %VS_PREFILTER  FIR filter on a virtual source's input (block by block,
    %zero latency, crossfade on a filter change): the response of a designed
    %loudspeaker's driver towards the listener (sound_scene.place_speaker).
    %   The partitioned convolver is (re)built for the block length it gets.

    properties (SetAccess = private)
        h = []                  % current FIR
        conv = []               % brir_convolver (second ear unused)
    end

    methods
        function set_filter(obj, h)
            obj.h = h(:);
            if ~isempty(obj.conv), obj.conv.set_filter([obj.h, zeros(size(obj.h))]); end
        end

        function y = process(obj, x)
            if isempty(obj.h), y = x; return; end
            if isempty(obj.conv) || obj.conv.block_size ~= numel(x)
                obj.conv = brir_convolver(numel(x));
                obj.conv.set_filter([obj.h, zeros(size(obj.h))]);
            end
            y = obj.conv.process(x(:));
            y = y(:, 1);
        end
    end
end
