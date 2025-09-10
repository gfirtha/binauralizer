classdef headphone < handle
    %VIRTUAL_SOURCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        source_index
        source_signal
        position
    end
    
    methods
        function obj = headphone(idx, position)
            obj.source_index = idx;
            obj.position = position;
            obj.source_signal = signal;

        end
        function obj = set_output(obj,varargin)
            if length(varargin) == 1
                obj.source_signal.set_signal(varargin{1});
            else
                obj.source_signal.set_signal(varargin{1},varargin{2});
            end
        end
    end
end

