classdef (Abstract) Directivity < handle
    % DIRECTIVITY  Abstract interface for source directivity patterns.
    %
    % Every sound source (loudspeaker or virtual source) can own a
    % Directivity object that knows how to compute its radiation pattern.
    %
    % Subclasses implement `evaluate(theta, f)` which returns the complex
    % directivity factor D(theta, f) for given angles and frequencies.
    %
    % evaluate_freq_domain caches its result so that repeated calls with the
    % same (N_fft, fs, theta) are free.  Sources with identical directivity
    % parameters can share the same Directivity object (flyweight pattern)
    % via Directivity.fromSourceType().

    properties (Access = protected)
        cache_key  = []    % [N_fft, fs, numel(theta)] — cache invalidation
        cache_D    = []    % [N_fft × Ntheta] cached matrix
    end

    methods (Abstract)
        D = evaluate(obj, theta, f)
    end

    methods
        function D = evaluate_freq_domain(obj, N_fft, fs, theta)
            % EVALUATE_FREQ_DOMAIN  Full FFT-length directivity matrix.
            %   Returns [N_fft × Ntheta].  Result is cached: repeated calls
            %   with the same (N_fft, fs, theta) return instantly.

            theta_v = theta(:)';
            key = [N_fft, fs, theta_v];
            if ~isempty(obj.cache_D) && numel(obj.cache_key) == numel(key) && all(obj.cache_key == key)
                D = obj.cache_D;
                return;
            end

            f = [(0:N_fft/2-1)'; (-N_fft/2:-1)'] / N_fft * fs;
            D = obj.evaluate(theta_v, f(:)).';  % [Ntheta×Nf].' → [N_fft×Ntheta]

            % Only cache large grids (small queries like single-angle are cheap)
            if numel(theta_v) >= 10
                obj.cache_key = key;
                obj.cache_D   = D;
            end
        end

        function invalidate_cache(obj)
            % INVALIDATE_CACHE  Force recomputation on next call.
            obj.cache_key = [];
            obj.cache_D   = [];
        end

        function Dval = evaluate_at(obj, theta, f)
            % EVALUATE_AT  Single angle, single frequency convenience.
            D = obj.evaluate(theta, f);
            Dval = D(1,1);
        end
    end

    methods (Static)
        function d = fromSourceType(source_type)
            % FROMSOURCETYPE  Flyweight factory: returns a shared Directivity
            % object for each unique (Shape, R) combination.
            %
            % Multiple loudspeakers with the same type get the SAME object,
            % so the cache is computed once and shared by all of them.

            persistent pool
            if isempty(pool)
                pool = containers.Map('KeyType','char','ValueType','any');
            end

            % Build a unique key from shape + radius
            switch source_type.Shape
                case 'point_source'
                    key = 'omni';
                case 'circular_piston'
                    key = sprintf('piston_%.6f', source_type.R(1));
                case 'two_way_speaker'
                    key = sprintf('twoway_%.6f_%.6f', source_type.R(1), source_type.R(2));
                otherwise
                    key = sprintf('omni_unknown_%s', source_type.Shape);
            end

            if pool.isKey(key)
                d = pool(key);
                return;
            end

            % Create new instance
            switch source_type.Shape
                case 'point_source'
                    d = OmnidirectionalDirectivity();
                case 'circular_piston'
                    d = PistonDirectivity(source_type.R(1));
                case 'two_way_speaker'
                    d = TwoWayDirectivity(source_type.R);
                otherwise
                    warning('Directivity:unknownType', ...
                        'Unknown source type "%s", falling back to omnidirectional.', source_type.Shape);
                    d = OmnidirectionalDirectivity();
            end

            pool(key) = d;
        end
    end
end
