classdef binaural_renderer < handle
    %BINAURAL_JOB Summary of this class goes here
    %   Detailed explanation goes here
    
    % Input: (set of?) binaural_sources
    % Output: 2-channel binaural signals
    properties
        binaural_source
        receiver
        hrtf_coefficients
        output_signal
    end
    properties (SetAccess = protected)
        binaural_filter
    end
    
    methods
        function obj = binaural_renderer(source, receiver)
            obj.binaural_source = source;
            obj.receiver = receiver;
            
            [theta,~] = cart2pol(obj.binaural_source.position(1)-obj.receiver.position(1),...
                                 obj.binaural_source.position(2)-obj.receiver.position(2));
            [~,ind] = min( sum(  (bsxfun(@minus, obj.binaural_source.hrtf.SourcePosition(:,[1,2]) ...
                                ,mod([theta,0]*180/pi - [obj.receiver.orientation,0] ,360))).^2,2  ) );
            obj.hrtf_coefficients = [squeeze(obj.binaural_source.hrtf.Data.IR(ind,1,:)),...
                                     squeeze(obj.binaural_source.hrtf.Data.IR(ind,2,:))];
            obj.binaural_filter  = OLS_convolver(obj.hrtf_coefficients, length(obj.binaural_source.source_signal.time_series),'spectrum');
            obj.output_signal = signal;
        end
        
        function obj = update_hrtf(obj)
            [theta,~] = cart2pol(obj.binaural_source.position(1)-obj.receiver.position(1),...
                obj.binaural_source.position(2)-obj.receiver.position(2));
            [~,ind] = min( sum(  (bsxfun(@minus, obj.binaural_source.hrtf.SourcePosition(:,[1,2]) ...
                           ,mod([theta,0]*180/pi - [obj.receiver.orientation,0],360)  )).^2, 2  ) );
            obj.hrtf_coefficients = [squeeze(obj.binaural_source.hrtf.Data.IR(ind,1,:)),...
                                     squeeze(obj.binaural_source.hrtf.Data.IR(ind,2,:))];
            obj.binaural_filter.update_coefficients(obj.hrtf_coefficients);
        end
        
        function obj = render(obj)
            [filter_out,ix] = obj.binaural_filter.convolve(obj.binaural_source.source_signal.time_series);
            obj.output_signal.set_spectrum(1/norm(obj.binaural_source.position-obj.receiver.position(1:2)) * filter_out, ix);
        end
    end
end

