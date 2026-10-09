classdef PistonDirectivity < Directivity
    % PISTONDIRECTIVITY  Circular rigid piston in an infinite baffle.
    %   D(theta, f) = 2 * J1(ka sin(theta)) / (ka sin(theta))
    %   where a = radius, k = 2*pi*f/c.
    %
    %   On-axis (theta=0) and DC (f=0): D = 1.

    properties
        radius   (1,1) double = 0.025   % piston radius (m)
    end

    methods
        function obj = PistonDirectivity(R)
            if nargin >= 1, obj.radius = R; end
        end

        function D = evaluate(obj, theta, f)
            theta = theta(:);     % [M×1]
            f = f(:)';            % [1×Nf]
            c = audioapp.util.Physics.speedOfSound();
            k = 2*pi * f / c;    % [1×Nf]

            [Theta, K] = meshgrid(theta, k);  % [Nf × M]
            arg = obj.radius * K .* sin(Theta);

            Dmx = 2 * besselj(1, arg) ./ arg;

            % Handle singularities: on-axis and DC
            Dmx(abs(arg) < 1e-12) = 1;
            Dmx(1,:) = 1;       % DC row

            % Hermitian symmetry at Nyquist (if even-length FFT)
            if mod(size(Dmx,1),2) == 0
                Dmx(end/2+1,:) = real(Dmx(end/2+1,:));
            end

            D = Dmx.';  % → [M × Nf]
        end
    end
end
