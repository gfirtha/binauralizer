classdef (Abstract) separator_base < handle
    %SEPARATOR_BASE  Streaming source separator / upmixer interface.
    %
    %   init(fs, nIn) once, then Y = process(X) block by block (any block
    %   length; output has the same length, delayed by Latency samples),
    %   reset() after a seek. Each output column is one object; Names and
    %   Pos (normalised room coordinates [right front up], as Cavern
    %   objects: unit circle around the listener at ear height) describe it.
    %   [0 0 0] = unpositioned (e.g. LFE).

    properties (SetAccess = protected)
        fs = 48000
        nIn = 2
        Names = {}
        Pos = zeros(0,3)
        Latency = 0                 % [samples]
    end

    methods (Abstract)
        init(obj, fs, nIn)
        Y = process(obj, X)
        reset(obj)
    end

    methods
        function n = nOut(obj)
            n = size(obj.Pos, 1);
        end
    end

    methods (Static)
        function p = dir2pos(az, el)
            % azimuth (deg, 0 = front, + = left), elevation (deg) -> [x y z]
            if nargin < 2, el = 0*az; end
            az = az(:); el = el(:);
            p = [-sind(az).*cosd(el), cosd(az).*cosd(el), sind(el)];
        end

        function [az, names] = channel_layout(n)
            % Azimuths of the standard input layouts (NaN = LFE)
            switch n
                case 1, az = 0; names = {'C'};
                case 2, az = [30 -30]; names = {'L','R'};
                case 5, az = [30 -30 0 110 -110]; names = {'L','R','C','Ls','Rs'};   % 5.0 (5.1 without LFE)
                case 6, az = [30 -30 0 NaN 110 -110]; names = {'L','R','C','LFE','Ls','Rs'};
                case 8, az = [30 -30 0 NaN 150 -150 90 -90]; names = {'L','R','C','LFE','Lb','Rb','Ls','Rs'};
                otherwise
                    az = mod(360*(0:n-1)/n + 180, 360) - 180;
                    names = arrayfun(@(k) sprintf('ch%d', k), 1:n, 'UniformOutput', false);
            end
        end
    end
end
