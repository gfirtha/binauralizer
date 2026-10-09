classdef loudspeaker < handle
    %LOUDSPEAKER  A single loudspeaker element in the SSD.
    %
    %   position    : 1×3  [x, y, z]  — full 3D position
    %   orientation : 1×2 or 1×3       — normal direction (2D or 3D)

    events
        position_changed
    end
    properties
        source_index      % SSD contour order (1..N, geometric)
        channelIndex      % physical output channel (default = source_index)
        source_signal
        position          % 1×3 [x, y, z]
        orientation
        source_type
        directivity       % Directivity object (lazy-init from source_type)
    end
    
    methods
        function obj = loudspeaker(idx, position, orientation, type, height_or_z)
            % LOUDSPEAKER  Constructor.
            %   loudspeaker(idx, pos, ori, type, z)
            %
            %   pos can be 1×2 [x,y] or 1×3 [x,y,z].
            %   If 1×2, the 5th argument is used as z (default 0).
            %   If 1×3, the 5th argument is ignored.

            obj.source_index = idx;
            obj.channelIndex = idx;
            obj.orientation = orientation;
            obj.source_signal = signal;
            obj.source_type = type;
            obj.directivity = Directivity.fromSourceType(type);

            % Ensure position is 1×3
            if numel(position) >= 3
                obj.position = position(1:3);
            else
                if nargin >= 5 && ~isempty(height_or_z)
                    obj.position = [position(:)', height_or_z];
                else
                    obj.position = [position(:)', 0];
                end
            end
        end

        function obj = set_position(obj, pos_in)
            if numel(pos_in) >= 3
                obj.position = pos_in(1:3);
            else
                obj.position(1:2) = pos_in(1:2);
            end
            notify(obj, 'position_changed');
        end

        function obj = set_orientation(obj, orientation_in)
            obj.orientation = orientation_in;
            notify(obj, 'position_changed');
        end
        
        function obj = set_output(obj,varargin)
            if length(varargin) == 1
                obj.source_signal.set_signal(varargin{1});
            else
                obj.source_signal.set_signal(varargin{1},varargin{2});
            end
        end

        % Legacy compatibility: .height reads/writes position(3)
        function h = get_height(obj)
            h = obj.position(3);
        end
        function obj = set_height(obj, h)
            obj.position(3) = h;
            notify(obj, 'position_changed');
        end
    end
end
