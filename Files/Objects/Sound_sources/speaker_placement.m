classdef speaker_placement < handle
    %SPEAKER_PLACEMENT  A designed loudspeaker (designed_loudspeaker) in the
    %scene: which SSD loudspeakers are its drivers, where it stands, where
    %it faces.
    %
    %   Elements(k) is the SSD loudspeaker (loudspeaker_array index) that is
    %   driver k: it sits at the driver's position, along its axis, and its
    %   directivity is the model's response (DesignedDirectivity;
    %   sound_scene.attach_speaker). Kept in Setup.DesignedSpeakers across
    %   scene rebuilds.

    properties
        Model                   % designed_loudspeaker
        Elements = []           % [1 x nCh] SSD loudspeaker index per driver
        Position = [0 0 0]      % [m] model origin in the scene
        Yaw = 0                 % [deg] facing direction (0 = +x)
    end

    methods
        function obj = speaker_placement(model, elements, position, yaw)
            obj.Model = model;
            obj.Elements = elements(:)';
            if nargin > 2, obj.Position = [position(:)', zeros(1, 3 - numel(position))]; end
            if nargin > 3, obj.Yaw = yaw; end
        end

        function R = rotation(obj)
            R = [cosd(obj.Yaw) -sind(obj.Yaw) 0; sind(obj.Yaw) cosd(obj.Yaw) 0; 0 0 1];
        end

        function [P, A] = driverWorld(obj)
            % driver positions [nCh x 3] and axes [nCh x 3] in the scene
            R = obj.rotation();
            D = obj.Model.Drivers;
            P = reshape([D.pos], 3, [])'*R' + obj.Position;
            A = reshape([D.axis], 3, [])'*R';
        end

        function d = toLocal(obj, v)
            % scene direction(s) [n x 3] -> local frame of the model
            d = v*obj.rotation();
        end
    end

    methods (Static)
        function obj = onElements(model, elements, pos, ori, listener)
            % placement whose drivers land where the SSD loudspeakers are:
            % pos / ori [n x >=2] their positions / axes, listener [1 x >=2].
            % One driver: exactly at the loudspeaker, along its axis. An
            % array: its driver centroid at their centroid, its front (+x)
            % along their mean axis (towards the listener if that cancels).
            n = numel(elements);
            pos = [pos, zeros(n, 3 - size(pos, 2))]; pos = pos(:, 1:3);
            D = model.Drivers;
            dp = reshape([D.pos], 3, [])'; da = reshape([D.axis], 3, [])';
            o = mean(ori(:, 1:2), 1);
            c = mean(pos, 1);
            if n == 1
                yaw = atan2d(o(2), o(1)) - atan2d(da(1,2), da(1,1));
                ref = dp(1, :);
            else
                if norm(o) < 1e-6, o = listener(1:2) - c(1:2); end
                yaw = atan2d(o(2), o(1));
                ref = mean(dp, 1);
            end
            R = [cosd(yaw) -sind(yaw) 0; sind(yaw) cosd(yaw) 0; 0 0 1];
            obj = speaker_placement(model, elements, c - ref*R', yaw);
        end
    end
end
