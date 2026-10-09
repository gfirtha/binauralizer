classdef room_model < handle
    %ROOM_MODEL  Room geometry + materials for the Room builder module.
    %
    %   PLACEHOLDER (module 3). Implemented: geometry import (.obj, .stl)
    %   and bounds. Planned:
    %     - surface grouping and materials (absorption / scattering per
    %       surface, APD files as in EntelFEM)
    %     - simulation backends: image sources, ray tracing (phased beam
    %       tracer), low-frequency FEM / FDTD - the EntelFEM engines are the
    %       reference; their classes clash with binauralizer names
    %       (receiver, signal), so they are to be called through a bridge,
    %       not put on the same path
    %     - export to the Binauralizer: BRIRs / image-source sets per
    %       loudspeaker for the room simulator renderer
    %
    %   m = room_model.from_file(file);  m.V [Nv x 3], m.F [Nf x 3|4]

    properties
        Name = ''
        File = ''
        V = zeros(0,3)          % vertices [m]
        F = zeros(0,3)          % faces (vertex indices, triangles or quads)
        SurfaceGroup = []       % [Nf x 1] group index per face (placeholder)
        Materials = struct('name', {'default'}, 'alpha', {0.1}, 'scattering', {0.1})
    end

    methods (Static)
        function m = from_file(file)
            m = room_model();
            [~, nm, ext] = fileparts(file);
            m.Name = nm; m.File = file;
            switch lower(ext)
                case '.stl'
                    T = stlread(file);
                    m.V = T.Points; m.F = T.ConnectivityList;
                case '.obj'
                    [m.V, m.F] = room_model.read_obj(file);
                otherwise
                    error('room_model:format', 'Unsupported geometry format "%s" (.obj, .stl).', ext);
            end
            m.SurfaceGroup = ones(size(m.F,1), 1);
        end

        function [V, F] = read_obj(file)
            % vertices ('v x y z') and faces ('f a b c [d]', a/b/c index forms)
            L = readlines(file);
            V = zeros(0,3); F = zeros(0,4);
            for k = 1:numel(L)
                s = strtrim(L(k));
                if startsWith(s, "v ")
                    V(end+1,:) = sscanf(extractAfter(s, 2), '%f', 3)'; %#ok<AGROW>
                elseif startsWith(s, "f ")
                    tok = split(strtrim(extractAfter(s, 2)));
                    idx = cellfun(@(t) sscanf(t, '%d', 1), cellstr(tok))';
                    idx(idx < 0) = size(V,1) + 1 + idx(idx < 0);
                    for j = 2:numel(idx)-1          % fan triangulation
                        F(end+1,:) = [idx(1) idx(j) idx(j+1) NaN]; %#ok<AGROW>
                    end
                end
            end
            F = F(:,1:3);
        end
    end

    methods
        function b = bounds(obj)
            % [min; max] corner of the bounding box
            b = [min(obj.V, [], 1); max(obj.V, [], 1)];
        end

        function n = nFaces(obj), n = size(obj.F, 1); end

        function simulate(~, varargin) %#ok<STOUT>
            error('room_model:notImplemented', ...
                'Room simulation is not implemented yet (planned: image sources, ray tracing, FEM/FDTD).');
        end

        function export_to_binauralizer(~, varargin)
            error('room_model:notImplemented', ...
                'Export to the Binauralizer is not implemented yet (planned: room simulator renderer input).');
        end
    end
end
