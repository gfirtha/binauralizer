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
            
            [theta,~] = cart2pol(obj.binaural_source.position(1)-obj.receiver(1),...
                                 obj.binaural_source.position(2)-obj.receiver(2));
            [~,ind] = min( sum(  (bsxfun(@minus, obj.binaural_source.hrtf.SourcePosition(:,[1,2]) ...
                ,mod([theta,0]*180/pi,360))).^2,2  ) );
            obj.hrtf_coefficients = [squeeze(obj.binaural_source.hrtf.Data.IR(ind,1,:)),...
                                     squeeze(obj.binaural_source.hrtf.Data.IR(ind,2,:))];
            obj.binaural_filter  = OLS_convolver(obj.hrtf_coefficients, length(obj.binaural_source.source_signal));
            obj.output_signal = zeros(size(obj.binaural_source.source_signal,1), 2);
        end
        
        function obj = update_hrtf(obj)
            [theta,~] = cart2pol(obj.binaural_source.position(1)-obj.receiver(1),...
                obj.binaural_source.position(2)-obj.receiver(2));
            [~,ind] = min( sum(  (bsxfun(@minus, obj.binaural_source.hrtf.SourcePosition(:,[1,2]) ...
                ,mod([theta,0]*180/pi,360))).^2,2  ) );
            obj.hrtf_coefficients = [squeeze(obj.binaural_source.hrtf.Data.IR(ind,1,:)),...
                                     squeeze(obj.binaural_source.hrtf.Data.IR(ind,2,:))];
            obj.binaural_filter.update_coefficients(obj.hrtf_coefficients);
        end
        
        function obj = render(obj)
            obj.output_signal = (1/norm(obj.binaural_source.position-obj.receiver) *...
                obj.binaural_filter.convolve(obj.binaural_source.source_signal));
            if size( obj.output_signal , 2 ) == 1
                obj.output_signal = repmat(obj.output_signal , 1, 2);
            end
        end
    end
end

