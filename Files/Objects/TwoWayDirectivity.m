classdef TwoWayDirectivity < Directivity
    % TWOWAYDIRECTIVITY  Two-way loudspeaker directivity model.
    %   Low-frequency driver (R(1)) and high-frequency driver (R(2)) with
    %   a 4th-order Linkwitz-Riley crossover at 400 Hz.
    %
    %   D(theta, f) = LP(f) * piston(R1, theta, f) + HP(f) * piston(R2, theta, f)

    properties
        R       (1,2) double = [0.05, 0.013]   % [woofer, tweeter] radii (m)
        fc      (1,1) double = 400              % crossover frequency (Hz)
        Nfilter (1,1) double = 4                % filter order
    end

    methods
        function obj = TwoWayDirectivity(R, fc, Nfilter)
            if nargin >= 1, obj.R = R; end
            if nargin >= 2, obj.fc = fc; end
            if nargin >= 3, obj.Nfilter = Nfilter; end
        end

        function D = evaluate(obj, theta, f)
            theta = theta(:);      % [M×1]
            f = f(:)';             % [1×Nf]
            c = audioapp.util.Physics.speedOfSound();

            wc = 2*pi * obj.fc;
            w  = 2*pi * f;        % [1×Nf]

            % Linkwitz-Riley crossover filters
            LP = (1 ./ sqrt(1 + (w/wc).^(2*obj.Nfilter))).^2;  % [1×Nf]
            HP = (1 ./ sqrt(1 + (wc./w).^(2*obj.Nfilter))).^2;  % [1×Nf]

            % Avoid division by zero at DC
            HP(abs(w) < 1e-12) = 0;
            LP(abs(w) < 1e-12) = 1;

            k = w / c;   % [1×Nf]

            [Theta, K] = meshgrid(theta, k);   % [Nf × M]

            % Woofer piston
            arg1 = obj.R(1) * K .* sin(Theta);
            D1 = 2 * besselj(1, arg1) ./ arg1;
            D1(abs(arg1) < 1e-12) = 1;

            % Tweeter piston
            arg2 = obj.R(2) * K .* sin(Theta);
            D2 = 2 * besselj(1, arg2) ./ arg2;
            D2(abs(arg2) < 1e-12) = 1;

            % Combine with crossover
            Dmx = bsxfun(@times, D1, LP(:)) + bsxfun(@times, D2, HP(:));

            % Group delay compensation
            fs_approx = max(f) * 2;  % approximate fs
            if fs_approx > 0
                delay = exp(-1i * w * obj.Nfilter / fs_approx / 2);
                Dmx = bsxfun(@times, Dmx, delay(:));
            end

            % DC and on-axis fixes
            Dmx(1,:) = 1;
            Dmx(:,1) = 1;  % theta=0 column

            D = Dmx.';  % → [M × Nf]
        end
    end
end
