classdef reverberant_source < handle
    %REVERBERANT_SOURCE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        source_index
        source_signal
        position
        orientation
        image_sources
    end
    
    methods
        function obj = reverberant_source(idx, position, orientation, room, N)
            %REVERBERANT_SOURCE Construct an instance of this class
            %   Detailed explanation goes here
            obj.source_index = idx;
            obj.position = position;
            obj.orientation = orientation;
            obj.source_signal = signal;
            obj.image_sources = [];
            obj = room.generate_mirror_source_tree(obj, N);
        end
        
        function D = get_positions(obj)
            %REVERBERANT_SOURCE Construct an instance of this class
            %   Detailed explanation goes here
            D = [obj.position cell2mat(cellfun( @(x) x.position,obj.image_sources,'UniformOutput',false))];

        end
        function obj = set_input(obj,varargin)
            if length(varargin) == 1
                obj.source_signal.set_signal(varargin{1});
            else
                obj.source_signal.set_signal(varargin{1},varargin{2});
            end
        end
    end
end

