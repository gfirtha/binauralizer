classdef binaural_renderer < handle
    %BINAURAL_RENDERER  Per-loudspeaker binaural rendering (HRTF + directivity).

    properties
        binaural_source
        receiver
        output_bus
        source_directivity      % legacy: directivity_table (kept for compat)
        output_signal
        extrapolator
        N_fft                   % FFT length for directivity evaluation
        fs                      % sample rate
        level_scale = 1         % calibrated playback: digital SSD signal -> ear pressure [Pa]
                                % (AppContext.applyLevelScales; 1 = uncalibrated, as before)
    end
    properties (SetAccess = protected)
        binaural_filter
    end

    methods
        function obj = binaural_renderer(source, receiver, output_bus, directivity, extrapolator)
            obj.binaural_source = source;
            obj.receiver = receiver;
            obj.output_bus = output_bus;
            obj.source_directivity = directivity;  % legacy fallback
            obj.extrapolator = extrapolator;

            obj.fs = source.source_signal.fs;
            obj.N_fft = length(source.source_signal.time_series);

            obj.binaural_filter  = OLS_convolver(obj.padHrtf(obj.extrapolator.extrapolate_hrtf(source,receiver)), ...
                                        obj.N_fft, 'spectrum');
            obj.output_signal = signal;
            obj.update_directivity;
        end

        function obj = update_hrtf(obj)
             obj.binaural_filter.update_coefficients( obj.padHrtf(obj.extrapolator.extrapolate_hrtf(obj.binaural_source,obj.receiver))  );
        end

        function hp = get_packed_ir(obj)
            % Two-ear impulse response (HRTF x directivity) of this
            % loudspeaker as h_L + 1i*h_R, exactly as the convolver applies
            % it. [] if the filter is partitioned (longer than a block).
            hp = [];
            if ~strcmp(obj.binaural_filter.mode, 'short_convolution'), return; end
            hp = ifft(obj.binaural_filter.coefficients{1});
        end

        function obj = update_directivity(obj)
            v_vec = (obj.receiver.position - obj.binaural_source.position);
            ori = obj.binaural_source.orientation;
            % Ensure same dimension for dot product
            D = min(numel(v_vec), numel(ori));
            theta0 = acos(max(-1, min(1, ...
                dot(v_vec(1:D), ori(1:D)) / (norm(v_vec) + eps) )));

            % The prefilter must match the convolver's internal FFT size
            Nfft_conv = obj.binaural_filter.N_fft;

            % Driver of a designed loudspeaker: its response at each ear
            % (spherical-wave model incl. near field, not only the angle; the
            % two ears differ in the near field), one prefilter per ear
            if isprop(obj.binaural_source, 'directivity') && isa(obj.binaural_source.directivity, 'DesignedDirectivity')
                [xL, xR] = obj.earPositions();
                Dd = obj.binaural_source.directivity;
                D_spectrum = [Dd.responseAt(xL, Nfft_conv, obj.fs), Dd.responseAt(xR, Nfft_conv, obj.fs)];
                obj.binaural_filter.update_prefilter(obj.level_scale*D_spectrum, 'frequency_domain');
            % Prefer the source's own Directivity object (new system)
            elseif isprop(obj.binaural_source, 'directivity') && ...
                    ~isempty(obj.binaural_source.directivity) && ...
                    isa(obj.binaural_source.directivity, 'Directivity')
                D_spectrum = obj.binaural_source.directivity.evaluate_freq_domain( ...
                    Nfft_conv, obj.fs, theta0);
                obj.binaural_filter.update_prefilter(obj.level_scale*D_spectrum(:,1), 'frequency_domain');
            elseif ~isempty(obj.source_directivity) && isprop(obj.source_directivity, 'theta')
                % Legacy fallback: use directivity_table
                [~,ind2] = min((obj.source_directivity.theta - theta0).^2);
                obj.binaural_filter.update_prefilter(obj.level_scale*obj.source_directivity.directivity_mx(:,ind2), 'frequency_domain');
            elseif obj.level_scale ~= 1
                obj.binaural_filter.update_prefilter(obj.level_scale*ones(Nfft_conv, 1), 'frequency_domain');
            end
        end

        function h = padHrtf(obj, h)
            % the HRTF zero-padded to the block length: the convolver's
            % filter length then also holds the directivity prefilter's
            % extension (HRTF * directivity; without it the tail was cut).
            % Free: the FFT size is 2 x block either way.
            if size(h, 1) < obj.N_fft, h = [h; zeros(obj.N_fft - size(h, 1), size(h, 2))]; end
        end

        function [xL, xR] = earPositions(obj)
            % the ears: +-ear_radius across the head, the nose along the
            % receiver orientation (hrtf_extrapolator: nose -> +x, left -> +y)
            p = obj.receiver.position(:)'; p = [p, zeros(1, 3 - numel(p))];
            o = obj.receiver.orientation(:)'; o = [o, zeros(1, 3 - numel(o))]; o = o/max(norm(o), eps);
            l = cross([0 0 1], o); if norm(l) < 1e-6, l = [0 1 0]; end
            l = l/norm(l); a = 0.0875;
            xL = p(1:3) + a*l; xR = p(1:3) - a*l;
        end

        function obj = update_renderer(obj,type)
            switch type
                case 'receiver_moved'
                    obj.update_hrtf;
                    obj.update_directivity;
                case 'receiver_rotated'
                    obj.update_hrtf;
                    if isprop(obj.binaural_source, 'directivity') && isa(obj.binaural_source.directivity, 'DesignedDirectivity')
                        obj.update_directivity;                 % the ears move with the head
                    end
                case 'loudspeaker_moved'
                    obj.update_hrtf;
                    obj.update_directivity;
                case 'loudspeaker_rotated'
                    obj.update_directivity;
            end
        end

        function obj = render(obj)
            [filter_out,ix] = obj.binaural_filter.convolve(obj.binaural_source.source_signal.time_series);
            obj.output_signal.set_spectrum(filter_out, ix);
            obj.output_bus.source_signal.add_spectra(obj.output_signal);
        end
    end
end

