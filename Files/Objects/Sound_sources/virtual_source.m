classdef virtual_source < handle
    %VIRTUAL_SOURCE Summary of this class goes here
    %   Detailed explanation goes here
    events
        position_changed
    end
    properties
        gain
        source_index
        source_signal
        position
        orientation
        source_type
        renderer_type
    end

    methods
        function obj = virtual_source(idx, position, orientation, source_type, renderer_type)
            %VIRTUAL_SOURCE Construct an instance of this class
            %   Detailed explanation goes here
            obj.source_index = idx;
            obj.position = position;
            obj.source_signal = signal;
            obj.orientation = orientation;
            obj.source_type = source_type;
            obj.renderer_type = renderer_type;
            obj.gain = 1;
        end

        function obj = set_position(obj, pos_in)
            obj.position = pos_in;
            notify(obj, 'position_changed');
        end
        function obj = set_input(obj,varargin)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            if length(varargin) == 1
                obj.source_signal.set_signal(obj.gain*varargin{1});
            else
                obj.source_signal.set_signal(obj.gain*varargin{1},varargin{2});
            end
        end

        function obj = set_gain(obj, gain)
            obj.gain = gain;
        end
        function gain = get_gain(obj)
            gain = obj.gain;
        end

        function field = get_sound_field(obj,X,Y,f0,t0)
            switch obj.source_type.Shape
                case 'point_source'
                    R = sqrt( (X-obj.position(1)).^2 + (Y-obj.position(2)).^2 );
                    field = obj.gain*1/4/pi*exp(1i*2*pi*f0*(t0-R/343.1))./R;
                case 'plane_wave'
                    k = -obj.position/norm(obj.position);
                    field = 1/4/pi*exp(1i*2*pi*f0*(t0-(X*k(1)+Y*k(2))/343.1));
            end
        end

    end
end

