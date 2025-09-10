classdef binaural_renderer < handle
    %BINAURAL_JOB Summary of this class goes here
    %   Detailed explanation goes here

    properties
        binaural_source
        receiver
        output_bus
        source_directivity
        output_signal
        extrapolator
    end
    properties (SetAccess = protected)
        binaural_filter
    end

    methods
        function obj = binaural_renderer(source, receiver, output_bus, directivity, extrapolator)
            obj.binaural_source = source;
            obj.receiver = receiver;
            obj.output_bus = output_bus;
            obj.source_directivity = directivity;
            obj.extrapolator = extrapolator;

            obj.binaural_filter  = OLS_convolver(obj.extrapolator.extrapolate_hrtf(source,receiver), ...
                                        length(obj.binaural_source.source_signal.time_series),'spectrum');
            obj.output_signal = signal;
            obj.update_directivity;
        end

        function obj = update_hrtf(obj)
             obj.binaural_filter.update_coefficients( obj.extrapolator.extrapolate_hrtf(obj.binaural_source,obj.receiver)  );
        end

        function obj = update_directivity(obj)
            v_vec = (obj.receiver.position-obj.binaural_source.position);
            theta0 = acos(v_vec*obj.binaural_source.orientation'/norm(v_vec));
            [~,ind2] = min((obj.source_directivity.theta - theta0).^2);
            obj.binaural_filter.update_prefilter(obj.source_directivity.directivity_mx(:,ind2),'frequency_domain');

        end

        function obj = update_renderer(obj,type)
            switch type
                case 'receiver_moved'
                    obj.update_hrtf;
                    obj.update_directivity;
                case 'receiver_rotated'
                    obj.update_hrtf;
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

