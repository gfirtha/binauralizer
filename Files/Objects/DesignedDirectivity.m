classdef DesignedDirectivity < Directivity
    % DESIGNEDDIRECTIVITY  An SSD loudspeaker that is driver Channel of a
    % designed loudspeaker (speaker_placement -> designed_loudspeaker).
    %
    %   The loudspeaker sits at the driver's position; the model's response
    %   is normalised with the propagation from the driver removed (H = 1 is
    %   a point source), so it multiplies the point source exactly like a
    %   directivity factor.
    %   responseAt(x, N_fft, fs): the full response at the point x (world
    %   frame), spherical-wave model incl. near field (binaural_renderer
    %   uses it at the listener). evaluate(theta, f): angle-only fallback
    %   for the field displays - the H model in the horizontal plane at the
    %   angle theta from the driver axis (left / right not distinguished).

    properties
        Placement           % speaker_placement
        Channel = 1         % driver index in the model
    end

    methods
        function obj = DesignedDirectivity(P, k)
            obj.Placement = P; obj.Channel = k;
        end

        function D = evaluate(obj, theta, f)
            % [M x Nf]
            S = obj.Placement.Model; f = f(:);
            a = S.Drivers(obj.Channel).axis; az0 = atan2(a(2), a(1));
            th = az0 + theta(:);
            Hm = S.directions(obj.Channel, [cos(th), sin(th), zeros(numel(th), 1)]);   % [M x nf stored]
            G = S.interpF(Hm.', abs(f));                                             % [Nf x M]
            neg = f < 0; G(neg, :) = conj(G(neg, :));
            D = G.';
        end

        function D = responseAt(obj, x, N_fft, fs)
            % [N_fft x 1] response at the world point x, FFT order of
            % Directivity.evaluate_freq_domain
            P = obj.Placement;
            x = [x(:)', zeros(1, 3 - numel(x))];
            xl = P.toLocal(x(1:3) - P.Position);
            fp = (0:N_fft/2)'/N_fft*fs;                     % the positive half; the rest Hermitian
            Dp = P.Model.response(obj.Channel, xl, fp);
            Dp(end) = real(Dp(end));
            D = [Dp; conj(Dp(end-1:-1:2))];
        end
    end
end
