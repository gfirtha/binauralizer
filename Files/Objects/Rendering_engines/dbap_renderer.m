classdef dbap_renderer < base_renderer
%DBAP_RENDERER  Distance-Based Amplitude Panning (DBAP)
%
% Gains: g_i = (1 / r_i^a) / sqrt(sum_j 1/r_j^(2a)), with r_i regularized
% by rs: r_i = sqrt(rs^2 + ||x_i - x_src||^2).
%
% setup fields used (all optional):
%   .DBAPExponent   -> a (default 1.0 ≈ 6 dB per doubling)
%   .DBAPRs         -> rs regularization in meters (default 0.1)

    properties
        setup
        G_vec  % column vector of per-loudspeaker gains
    end

    methods
        function obj = dbap_renderer(virtual_source, loudspeaker_array, setup)
            % Constructor matches base_renderer signature used elsewhere
            obj = obj@base_renderer(virtual_source, loudspeaker_array);

            % Small defaults; merge with provided setup
            def = struct('DBAPExponent', 1.0, ...   % 6 dB per doubling
                         'DBAPRs',       0.1);
            if nargin < 3 || isempty(setup), setup = struct(); end
            obj.setup = obj.mergeSetup(def, setup);

            % Allocate output_signal one-per-SSD element
            for n = 1:numel(obj.ssd.loudspeakers)
                obj.output_signal{n} = signal; %#ok<AGROW>
            end

            obj.update_renderer();
        end

        function obj = update_settings(obj, setup)
            % Accept a whole app.Setup struct or just DBAP fields
            obj.setup = obj.mergeSetup(obj.setup, setup);
            obj.update_renderer();
        end

        function obj = update_renderer(obj, ~)
            % Compute DBAP gains for current source & SSD
            X  = obj.ssd.x0;                 % N×2 loudspeaker positions
            N  = size(X,1);
            a  = obj.setup.DBAPExponent;
            rs = obj.setup.DBAPRs;

            % Source "position"
            switch lower(obj.virtual_source.source_type.Shape)
                case 'plane_wave'
                    % Treat as a very far point in the propagation direction
                    dir = ([0,0] - obj.virtual_source.position(:).');
                    dir = dir / max(norm(dir),1e-12);
                    xsrc = dir * 100;  % 100 m surrogate
                otherwise % 'point_source'
                    xsrc = obj.virtual_source.position(:).';
            end

            d2  = sum( (X - xsrc).^2, 2 );           % squared Euclidean
            r   = sqrt( rs^2 + d2 );                 % regularized distance

            inv_ra   = 1 ./ (r.^a);                  % 1/r^a
            norm_fac = 1 / sqrt( sum( inv_ra.^2 ) ); % energy preserving
            g        = norm_fac * inv_ra;

            % Store as column and keep length consistent with SSD
            obj.G_vec = reshape(g, [N 1]);
        end

        function render(obj)
            % Apply gains to time series and push to SSD bus
            x = obj.virtual_source.source_signal.time_series;
            for n = 1:numel(obj.output_signal)
                obj.output_signal{n}.set_signal( obj.G_vec(n) * x );
            end
            obj.add_output_to_ssd_signal;
        end

        % Optional: a concise info string (shows up in UI if you use it)
        function info = get_renderer_info(obj)
            lines = [
                "DBAP renderer — distance-based amplitude panning"
                sprintf("Exponent a = %.3g, regularization rs = %.3g m", ...
                        obj.setup.DBAPExponent, obj.setup.DBAPRs)
                "Gains normalize energy: ∑ g_i^2 = 1."
            ];
            info = strjoin(lines, newline);
        end
    end

    methods (Access = private)
        function out = mergeSetup(~, base, in)
            % Merge fields from 'in' (struct or nested Renderer_setup) into base
            out = base;
            if isstruct(in)
                % Accept either flat fields or Setup.Renderer_setup.* style
                cand = in;
                if isfield(in,'Renderer_setup') && isstruct(in.Renderer_setup)
                    cand = in.Renderer_setup;
                end
                f = fieldnames(cand);
                for k = 1:numel(f)
                    out.(f{k}) = cand.(f{k});
                end
            end
        end
    end
end
