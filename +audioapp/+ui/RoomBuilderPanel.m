classdef RoomBuilderPanel < handle
    % RoomBuilderPanel  Top-level "Room builder" tab (module 3).
    %
    %   PLACEHOLDER. Working: import a room geometry (.obj / .stl) and view
    %   it (room_model). Planned (disabled controls): surface groups and
    %   materials, simulation (image sources / ray tracing / FEM-FDTD) and
    %   sending the result to the Binauralizer.

    properties
        C
        Root
        TxtFile
        TxtInfo
        LstSurfaces
        PopMaterial
        PopMethod
        Ax
        Model = []          % room_model
    end

    properties (Constant)
        SimMethods = {'Image sources', 'Ray tracing (phased beams)', 'Low-frequency FEM / FDTD'}
    end

    methods
        function obj = RoomBuilderPanel(parent, C)
            obj.C = C;
            obj.build(parent);
        end

        function build(obj, parent)
            u = {'Units', 'normalized'};
            bg = get(ancestor(parent, 'figure'), 'Color');
            obj.Root = uipanel('Parent', parent, u{:}, 'Position', [0 0 1 1], 'BorderType', 'none');
            L = uipanel('Parent', obj.Root, u{:}, 'Position', [0.01 0.02 0.28 0.96], 'Title', 'Room builder');
            txt = @(s, y, h) uicontrol(L, 'Style', 'text', 'String', s, u{:}, 'Position', [0.04 y 0.92 h], ...
                'HorizontalAlignment', 'left', 'BackgroundColor', bg);

            % --- geometry ---
            txt('Geometry:', 0.94, 0.035);
            uicontrol(L, 'Style', 'pushbutton', 'String', 'Import geometry (.obj / .stl)...', u{:}, ...
                'Position', [0.04 0.885 0.92 0.05], 'Callback', @(~,~) obj.onImport());
            obj.TxtFile = txt('-', 0.84, 0.04); set(obj.TxtFile, 'FontAngle', 'italic');
            obj.TxtInfo = txt('', 0.80, 0.04);

            % --- surfaces / materials (planned) ---
            txt('Surfaces / materials (planned):', 0.75, 0.035);
            obj.LstSurfaces = uicontrol(L, 'Style', 'listbox', 'String', {}, u{:}, ...
                'Position', [0.04 0.55 0.92 0.19], 'Enable', 'off');
            obj.PopMaterial = uicontrol(L, 'Style', 'popupmenu', 'String', {'default (alpha 0.1)'}, u{:}, ...
                'Position', [0.04 0.50 0.92 0.04], 'Enable', 'off');

            % --- simulation (planned) ---
            txt('Simulation (planned):', 0.43, 0.035);
            obj.PopMethod = uicontrol(L, 'Style', 'popupmenu', 'String', obj.SimMethods, u{:}, ...
                'Position', [0.04 0.385 0.92 0.04], 'Enable', 'off');
            uicontrol(L, 'Style', 'pushbutton', 'String', 'Run simulation', u{:}, ...
                'Position', [0.04 0.32 0.92 0.05], 'Enable', 'off');
            uicontrol(L, 'Style', 'pushbutton', 'String', 'Send to Binauralizer', u{:}, ...
                'Position', [0.04 0.17 0.92 0.055], 'FontWeight', 'bold', 'Enable', 'off');
            txt(['Planned: materials per surface group, image-source / ray-tracing / ' ...
                 'FEM-FDTD simulation, export of BRIRs to the room simulator renderer.'], 0.02, 0.12);

            obj.Ax = axes('Parent', obj.Root, u{:}, 'Position', [0.35 0.08 0.62 0.86]);
            title(obj.Ax, 'Room geometry'); axis(obj.Ax, 'equal'); grid(obj.Ax, 'on'); view(obj.Ax, 3);
            xlabel(obj.Ax, 'x [m]'); ylabel(obj.Ax, 'y [m]'); zlabel(obj.Ax, 'z [m]');
        end
    end

    methods (Access = private)
        function onImport(obj)
            [f, p] = uigetfile({'*.obj;*.stl', 'Room geometry (*.obj, *.stl)'}, 'Import room geometry');
            if isequal(f, 0), return; end
            try
                obj.Model = room_model.from_file(fullfile(p, f));
            catch ME
                set(obj.TxtInfo, 'String', ['Failed: ' ME.message]); return;
            end
            m = obj.Model; b = m.bounds();
            set(obj.TxtFile, 'String', f);
            set(obj.TxtInfo, 'String', sprintf('%d vertices, %d faces, %.1f x %.1f x %.1f m', ...
                size(m.V,1), m.nFaces(), b(2,1)-b(1,1), b(2,2)-b(1,2), b(2,3)-b(1,3)));
            set(obj.LstSurfaces, 'String', {sprintf('all faces (%d)', m.nFaces())});
            cla(obj.Ax);
            patch(obj.Ax, 'Vertices', m.V, 'Faces', m.F, 'FaceColor', [0.85 0.88 0.92], ...
                'FaceAlpha', 0.35, 'EdgeColor', [0.3 0.3 0.4]);
            axis(obj.Ax, 'equal'); grid(obj.Ax, 'on'); view(obj.Ax, 3);
            title(obj.Ax, ['Room geometry: ' m.Name], 'Interpreter', 'none');
        end
    end
end
