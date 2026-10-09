classdef motion_profile
    %MOTION_PROFILE  How a source travels along a trajectory.
    %
    %   Before T0 the source waits at arc length S0. From T0 on (tau = t-T0):
    %     - no Profile:  v(tau) = V0 + A*tau (stops when v reaches 0)
    %     - Profile [tau_k v_k] (rows): piecewise-linear speed keyframes,
    %       constant before the first / after the last key
    %   EndMode at the end of the path: 'stop' | 'loop' | 'pingpong'.

    properties
        T0 = 0          % start time [s] (media time)
        V0 = 1          % initial speed [m/s]
        A  = 0          % acceleration [m/s^2]
        Profile = zeros(0,2)
        EndMode = 'loop'
        S0 = 0          % start offset along the path [m]
    end

    methods
        function obj = motion_profile(varargin)
            for k = 1:2:numel(varargin)
                obj.(varargin{k}) = varargin{k+1};
            end
        end

        function d = distance(obj, t)
            % Path length travelled at media time(s) t (before EndMode)
            tau = max(t(:) - obj.T0, 0);
            if isempty(obj.Profile)
                v0 = max(obj.V0, 0); a = obj.A;
                if a < 0 && v0 > 0
                    tstop = v0/(-a);                % decelerates to standstill
                    tau = min(tau, tstop);
                elseif a < 0
                    tau = 0*tau;
                end
                d = v0*tau + a*tau.^2/2;
            else
                K = sortrows(obj.Profile, 1);
                tk = K(:,1); vk = max(K(:,2), 0);
                if tk(1) > 0, tk = [0; tk]; vk = [vk(1); vk]; end
                sk = [0; cumsum(diff(tk).*(vk(1:end-1) + vk(2:end))/2)];
                d = zeros(size(tau));
                for i = 1:numel(tau)
                    j = find(tk <= tau(i), 1, 'last');
                    dt = tau(i) - tk(j);
                    if j == numel(tk)
                        d(i) = sk(j) + vk(j)*dt;
                    else
                        g = (vk(j+1) - vk(j))/(tk(j+1) - tk(j));
                        d(i) = sk(j) + vk(j)*dt + g*dt^2/2;
                    end
                end
            end
        end

        function s = s_at(obj, t, L)
            % Arc length on a path of length L at media time(s) t
            s = obj.S0 + obj.distance(t);
            if L <= 0, s = 0*s; return; end
            switch lower(obj.EndMode)
                case 'loop'
                    s = mod(s, L);
                case 'pingpong'
                    s = mod(s, 2*L);
                    s(s > L) = 2*L - s(s > L);
                otherwise                            % 'stop'
                    s = min(s, L);
            end
        end

        function S = to_struct(obj)
            S = struct('T0', obj.T0, 'V0', obj.V0, 'A', obj.A, 'Profile', obj.Profile, ...
                'EndMode', obj.EndMode, 'S0', obj.S0);
        end
    end

    methods (Static)
        function obj = from_struct(S)
            obj = motion_profile();
            f = fieldnames(S);
            for k = 1:numel(f)
                v = S.(f{k});
                if strcmp(f{k}, 'Profile'), v = reshape(v, [], 2); end
                obj.(f{k}) = v;
            end
        end

        function K = parse_profile(str)
            % '0:1, 5:3, 10:0'  ->  [0 1; 5 3; 10 0]   ('' -> no profile)
            K = zeros(0,2);
            str = strtrim(char(str));
            if isempty(str), return; end
            parts = strsplit(str, {',', ';'});
            for k = 1:numel(parts)
                tv = sscanf(strrep(parts{k}, ':', ' '), '%f');
                if numel(tv) == 2, K(end+1,:) = tv'; end %#ok<AGROW>
            end
        end

        function str = format_profile(K)
            if isempty(K), str = ''; return; end
            str = strjoin(arrayfun(@(i) sprintf('%g:%g', K(i,1), K(i,2)), 1:size(K,1), 'UniformOutput', false), ', ');
        end
    end
end
