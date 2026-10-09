classdef object_slots < handle
    %OBJECT_SLOTS  K output objects fed by tracks (soundfield_tracker):
    %which track drives which object, its position and fade gain, and the
    %soft assignment weights of a component to the objects.
    %
    %   A slot takes a confirmed track that is novel (an event) or MinAge
    %   old, fades in (FadeIn), follows it (PosTau if the track moves,
    %   PosTauStill if it stands), and after the track died fades out
    %   (FadeOut) at its last position before it can be reused - a
    %   sounding object never jumps. Two standing objects closer than
    %   DupDeg: the younger fades out, its track is blocked.
    %   weights(azS, azI): Gaussian (AssignSigma) x gain^2 within AssignDeg,
    %   by the smoothed direction azS for standing objects and the
    %   instantaneous azI for moving ones; normalised if they sum above 1
    %   (the rest belongs to the bed).

    properties
        K = 6
        MinAge = 0.1
        FadeIn = 0.02
        FadeOut = 0.1
        PosTau = 0.015
        PosTauStill = 0.5
        AssignSigma = 15
        AssignDeg = 40
        DupDeg = 20
        StableTime = 0.1
    end

    properties (SetAccess = private)
        id; az; g                       % track id (0 = free), azimuth [deg], fade gain
        moving; active; renewed         % per slot flags (renewed: got a new track now)
        blocked = zeros(0, 2)           % [blocked track id, blocking track id]
    end

    methods
        function obj = object_slots(varargin)
            for k = 1:2:numel(varargin), obj.(varargin{k}) = varargin{k+1}; end
            obj.reset();
        end

        function reset(obj)
            K = obj.K;
            obj.id = zeros(1, K); obj.g = zeros(1, K);
            obj.az = linspace(0, 360, K + 1); obj.az = obj.az(1:K);
            obj.moving = false(1, K); obj.active = false(1, K); obj.renewed = false(1, K);
            obj.blocked = zeros(0, 2);
        end

        function update(obj, tr, t, dt)
            % tr: live tracks of soundfield_tracker.step
            K = obj.K;
            ids = [tr.id]; azs = arrayfun(@(s) s.x(1), tr);
            % a block lasts only while the blocking track lives
            obj.blocked = obj.blocked(ismember(obj.blocked(:,2), ids), :);
            elig = [tr.conf] & ([tr.novel] | t - [tr.age] >= obj.MinAge) & ~ismember(ids, obj.blocked(:,1));
            target = zeros(1, K); obj.moving(:) = false; obj.renewed(:) = false;
            for k = 1:K
                if obj.id(k) == 0, continue; end
                j = find(ids == obj.id(k) & [tr.conf], 1);
                if isempty(j) || ismember(obj.id(k), obj.blocked(:,1))
                    if obj.g(k) < 0.01, obj.id(k) = 0; end          % faded out: free
                    continue;
                end
                target(k) = 1;
                da = mod(azs(j) - obj.az(k) + 180, 360) - 180;
                obj.moving(k) = tr(j).still < obj.StableTime;
                if obj.moving(k), tp = obj.PosTau; else, tp = obj.PosTauStill; end
                obj.az(k) = obj.az(k) + da*(1 - exp(-dt/tp));
            end
            for k = find(target & ~obj.moving)
                for q = find(target & ~obj.moving)
                    if q ~= k && obj.id(q) < obj.id(k) && ...
                            abs(mod(obj.az(q) - obj.az(k) + 180, 360) - 180) < obj.DupDeg
                        obj.blocked(end+1,:) = [obj.id(k), obj.id(q)]; target(k) = 0; break;
                    end
                end
            end
            if size(obj.blocked, 1) > 200, obj.blocked = obj.blocked(end-99:end, :); end
            for j = find(elig)
                if any(obj.id == ids(j)), continue; end
                k = find(obj.id == 0, 1);
                if isempty(k), break; end                   % all slots busy
                obj.id(k) = ids(j); obj.az(k) = azs(j); target(k) = 1; obj.renewed(k) = true;
            end
            obj.az = mod(obj.az + 180, 360) - 180;
            up = target > obj.g;
            tau = obj.FadeOut*ones(1, K); tau(up) = obj.FadeIn;
            obj.g = obj.g + (target - obj.g).*(1 - exp(-dt./tau));
            obj.active = obj.g > 0.01;
        end

        function [w, sw] = weights(obj, azS, azI)
            d = abs(mod(obj.az - azS + 180, 360) - 180);
            dm = abs(mod(obj.az - azI + 180, 360) - 180);
            d(obj.moving) = dm(obj.moving);
            w = obj.g.^2 .* exp(-d.^2/(2*obj.AssignSigma^2));
            w(~obj.active | d > obj.AssignDeg) = 0;
            sw = sum(w);
            if sw > 1, w = w/sw; sw = 1; end
        end
    end
end
