classdef binaural_renderer < handle
    %BINAURAL_JOB Summary of this class goes here
    %   Detailed explanation goes here
    
    % Input: (set of?) binaural_sources
    % Output: 2-channel binaural signals
    properties
        binaural_source
        receiver
        hrtf_coefficients
        source_directivity
        output_signal
    end
    properties (SetAccess = protected)
        binaural_filter
    end
    
    methods
        function obj = binaural_renderer(source, receiver, directivity)
            obj.binaural_source = source;
            obj.receiver = receiver;
            obj.source_directivity = directivity;
            v_vec = obj.binaural_source.position-obj.receiver.position;
            v_vec = v_vec/norm(v_vec);
            
            [theta_1,~] = cart2pol(v_vec(1), v_vec(2));
            [theta_2,~] = cart2pol(obj.receiver.orientation(1),obj.receiver.orientation(2));
            [~,ind] = min( sum(  (bsxfun(@minus, obj.binaural_source.hrtf.SourcePosition(:,[1,2]) ...
                                ,mod([theta_1 - theta_2,0]*180/pi, 360))).^2,2  ) );
            obj.hrtf_coefficients = [squeeze(obj.binaural_source.hrtf.Data.IR(ind,1,:)),...
                                     squeeze(obj.binaural_source.hrtf.Data.IR(ind,2,:))];
            obj.binaural_filter  = OLS_convolver(obj.hrtf_coefficients, length(obj.binaural_source.source_signal.time_series),'spectrum');
            obj.output_signal = signal;

            theta_dir = acos(-v_vec*obj.binaural_source.orientation');
            [~,ind2] = min((obj.source_directivity.theta - theta_dir).^2);
            obj.binaural_filter.update_prefilter(obj.source_directivity.directivity_mx(:,ind2),'frequency_domain');
        end

        function obj = update_hrtf(obj)
            v_vec = obj.binaural_source.position-obj.receiver.position;
            v_vec = v_vec/norm(v_vec);
            [theta_1,~] = cart2pol(v_vec(1), v_vec(2));
            [theta_2,~] = cart2pol(obj.receiver.orientation(1),obj.receiver.orientation(2));
            [~,ind] = min( sum(  (bsxfun(@minus, obj.binaural_source.hrtf.SourcePosition(:,[1,2]) ...
                                ,mod([theta_1 - theta_2,0]*180/pi, 360))).^2,2  ) );
            obj.hrtf_coefficients = [squeeze(obj.binaural_source.hrtf.Data.IR(ind,1,:)),...
                                     squeeze(obj.binaural_source.hrtf.Data.IR(ind,2,:))];
            obj.binaural_filter.update_coefficients(obj.hrtf_coefficients);
        end

        function obj = update_directivity(obj)
            v_vec = (obj.receiver.position-obj.binaural_source.position); 
            v_vec = v_vec/norm(v_vec);
            theta0 = acos(v_vec*obj.binaural_source.orientation');
            [~,ind2] = min((obj.source_directivity.theta - theta0).^2);
            obj.binaural_filter.update_prefilter(obj.source_directivity.directivity_mx(:,ind2),'frequency_domain');
        end
         
        function obj = render(obj)
            [filter_out,ix] = obj.binaural_filter.convolve(obj.binaural_source.source_signal.time_series);
            obj.output_signal.set_spectrum(1/norm(obj.binaural_source.position-obj.receiver.position(1:2)) * filter_out, ix);
        end
    end
end

