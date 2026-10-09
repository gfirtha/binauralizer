classdef sound_scene < handle
    %SOUND_SCENE Summary of this class goes here
    %   Should contain description for the objects, present in the sound
    %   scene and interfaces positions with the gui

    properties
        virtual_sources % Virtual source: WFS or HOA or stereo etc source (should be a class with input signal,)
        loudspeaker_array  % SSD for SFS
        headphone
        receiver
        env
        scene_renderer
        downmixing_matrix
        channelMap        % permutation vector: SSD index → physical channel
        % [] or 1:N = identity (SSD order == channel order)
        TestChannelFcn    % @(chIdx) callback — set externally by AppContext
        OnGeometryChangedFcn  % @() callback — fired when receiver/VS moves
    end

    properties (Access = private)
        % Drag coalescing: mouse-motion callbacks only move the model and
        % record what changed; the renderer updates run at most every
        % UpdateDelay seconds (and once more after the last motion).
        pendingUpdates = struct('receiver','', 'vs',[], 'vsAll',false, 'vsRot',false, ...
                                'ls',[], 'lsRot',[])
        updateTimer = []
        UpdateDelay = 0.03
        UpdateBudget = 0.015        % [s] renderer + BRIR update time per flush for moved sources
        styleCounter = 0            % object markers are restyled every 3rd block
        gui_ref                     % listener_space_axes, for the overlays

        % Loudspeaker test signal (double-click on a loudspeaker): mixed
        % in at loudspeaker-signal level, i.e. into that output channel, or
        % in binaural mode through that loudspeaker's two-ear response
        testSig  = []               % samples still to play (column)
        testLs   = 0                % loudspeaker_array index
        testH    = []               % [Nh x 2] two-ear IR of the loudspeaker
        testTail = []               % binaural overlap-add tail
    end

    properties
        ObjectRoomRadius = 3        % [m] half-width of the Atmos room in the scene
        ObjectRoomHeight = 1.5      % [m] Atmos ceiling (up = 1) above ear height (up = 0)
        ObjectPlayback = false      % object-based input drives the sources
    end

    methods (Static, Access = private)
        function P = emptyPending()
            P = struct('receiver','', 'vs',[], 'vsAll',false, 'vsRot',false, 'ls',[], 'lsRot',[]);
        end

        function [position, normal] = default_input_layout(nIn, R)
            % Media channels are in file/channel order, not SSD contour order.
            % For standard layouts, invert get_default_layout's channel map
            % so row k describes input channel k.
            switch nIn
                case 5
                    [p0, n0, chMap] = get_default_layout(5, R, '5.0');
                case 6
                    [p0, n0, chMap] = get_default_layout(6, R, '5.1');
                otherwise
                    [position, normal] = get_default_layout(nIn, R, 'circular');
                    return;
            end

            position = zeros(nIn, size(p0, 2));
            normal = zeros(nIn, size(n0, 2));
            for k = 1:min(size(p0, 1), numel(chMap))
                ch = chMap(k);
                if ch >= 1 && ch <= nIn
                    position(ch,:) = p0(k,:);
                    normal(ch,:) = n0(k,:);
                end
            end
        end
    end

    methods
        function obj = sound_scene(gui, setup, varargin)
            obj.env = environment;
            R0 = setup.Loudspeaker_setup.R;


            if isempty(gui.room)

                if isempty(varargin)
                    rec_pos = [0,0];
                    rec_orientation = [0,1];
                    [pos_vs,n_vs] = sound_scene.default_input_layout(setup.N_in, R0 + 1);
                else
                    rec_pos = varargin{1}.rec_position;
                    rec_orientation = varargin{1}.rec_orientation;
                    pos_vs = varargin{1}.vs_position;
                    n_vs = varargin{1}.vs_orientation;
                end
                [pos_ssd, n_ssd, ssdChannelMap] = get_default_layout(setup.Loudspeaker_setup.N, R0, setup.Loudspeaker_setup.Shape, setup.Loudspeaker_setup);
            else
                if isempty(varargin)
                    rec_pos = [1,1];
                    rec_orientation = [0,1];
                else
                    rec_pos = varargin{1}.rec_position;
                    rec_orientation = varargin{1}.rec_orientation;
                end
                R0 = setup.Loudspeaker_setup.R;
                [pos_vs,n_vs] = sound_scene.default_input_layout(setup.N_in, R0 + 1);
                pos_vs(:,1:2) = pos_vs(:,1:2) + [1,1];
                [pos_ssd, n_ssd, ssdChannelMap] = get_default_layout(setup.Loudspeaker_setup.N, R0*0.2, 'circular');
                pos_ssd(:,1:2) = pos_ssd(:,1:2) + [1,1];
            end

            % Store channel map: prefer explicit Setup.ChannelMap, else use layout default
            if isfield(setup, 'ChannelMap') && ~isempty(setup.ChannelMap)
                ssdChannelMap = setup.ChannelMap;
            end

            % z=0 is ear height reference. All positions default to z=0.
            % 3D layouts (e.g. Atmos) specify z offsets relative to this.
            if numel(rec_pos) < 3
                rec_pos = [rec_pos(:)', 0];
            end

            sourceNames = {};
            if isfield(setup, 'SourceNames') && ~isempty(setup.SourceNames)
                if isstring(setup.SourceNames)
                    sourceNames = cellstr(setup.SourceNames);
                elseif ischar(setup.SourceNames)
                    sourceNames = {setup.SourceNames};
                elseif iscell(setup.SourceNames)
                    sourceNames = setup.SourceNames(:).';
                end
            end

            obj.create_receiver(gui, rec_pos, rec_orientation);
            gui.Move3DFcn = @(kind, idx, xyz) obj.move_from_view(kind, idx, xyz);
            for n = 1 : setup.N_in
                sourceName = [];
                if n <= numel(sourceNames)
                    sourceName = sourceNames{n};
                end
                obj.create_virtual_source(pos_vs(n,:),n_vs(n,:), gui, setup.Virtual_source_type, setup.Rendering_mode, sourceName);
                obj.virtual_sources{n}.set_input(zeros(setup.Block_size,1),setup.SampleRate);
            end
            switch setup.Rendering_mode
                case 'direct_playback'
                    for n = 1 : length(obj.virtual_sources)
                        obj.create_loudspeaker( obj.virtual_sources{n}.position,...
                            obj.virtual_sources{n}.orientation, ...
                            gui, setup.Loudspeaker_type);
                        obj.loudspeaker_array{n}.set_output(zeros(setup.Block_size,1),setup.SampleRate);
                        if isvalid(gui.virtual_source_points{n}.UserData.text)
                            gui.loudspeaker_points{n}.UserData.text = gui.virtual_source_points{n}.UserData.text;
                        end
                        delete(gui.virtual_source_points{n});
                    end
                    ssdChannelMap = 1:setup.N_in;
                otherwise
                    for n = 1 : size(pos_ssd,1)   % the layout decides (e.g. soundbar + surrounds)
                        obj.create_loudspeaker( pos_ssd(n,:), n_ssd(n,:), gui, setup.Loudspeaker_type);
                        obj.loudspeaker_array{n}.set_output(zeros(setup.Block_size,1),setup.SampleRate);
                    end
                    if isempty(gui.room)
                        gui.main_axes.XLim = 1.5*(R0+0.5)*[-1,1];
                        gui.main_axes.YLim = 1.5*(R0+0.5)*[-1,1];
                    end
            end
            obj.headphone = headphone(1,obj.receiver.position);
            obj.scene_renderer = sound_scene_renderer(obj.virtual_sources,obj.loudspeaker_array, obj.headphone, obj.receiver, setup);

            obj.downmixing_matrix = get_downmixing_mx(obj.loudspeaker_array,setup.N_out);

            % Channel mapping: SSD contour order → physical output channel
            % The SSD internally reorders its loudspeakers by angle for
            % geometry operations, but sound_scene.loudspeaker_array stays
            % in the original order. Since render_sound_scene reads signals
            % from loudspeaker_array (original order), channelMap must match
            % that original order — no permutation needed.
            obj.channelMap = ssdChannelMap;
            for k = 1:numel(obj.loudspeaker_array)
                if k <= numel(obj.channelMap)
                    obj.loudspeaker_array{k}.channelIndex = obj.channelMap(k);
                end
            end
            gui.refresh_view();     % markers of the current (side/3D) view
        end

        function start_test_signal(obj, k, sig)
            % Play sig from loudspeaker k (loudspeaker_array index) on top of
            % whatever the scene renders; restarts a running test signal.
            obj.testSig = sig(:);
            obj.testLs = k;
            obj.testH = obj.loudspeaker_ir(k);
            obj.testTail = [];
        end

        function tf = test_active(obj)
            tf = ~isempty(obj.testSig) || any(obj.testTail(:));
        end

        function h = loudspeaker_ir(obj, k)
            % Two-ear impulse response (HRTF x directivity) of loudspeaker k,
            % exactly as the per-loudspeaker binaural chain applies it
            h = [];
            R = obj.scene_renderer;
            if isempty(R), return; end
            idx = obj.loudspeaker_array{k}.source_index;
            m = find(cellfun(@(l) l.source_index == idx, R.ssd.loudspeakers), 1);
            if isempty(m), return; end
            hp = R.binaural_renderer{m}.get_packed_ir();
            if isempty(hp), return; end
            Nh = find(hp ~= 0, 1, 'last');
            h = [real(hp(1:Nh)), imag(hp(1:Nh))];
        end

        function move_from_view(obj, kind, idx, xyz)
            % Drag in a side view of the scene axes (x and height)
            switch kind
                case 'vs'
                    obj.virtual_sources{idx}.set_position(xyz);
                    obj.face_origin(idx);
                    obj.schedule_update('virtual_source_moved', idx);
                case 'ls'
                    obj.loudspeaker_array{idx}.set_position(xyz);
                    obj.schedule_update('loudspeaker_moved', idx);
            end
        end

        function obj = change_renderer(obj, setup)
            obj.scene_renderer = sound_scene_renderer(obj.virtual_sources,obj.loudspeaker_array, obj.headphone, obj.receiver, setup);
        end

        function obj = create_receiver(obj, gui, pos,orientation)
            R = 0.0875;  % average head radius (m)
            obj.receiver = receiver(pos,orientation);
            gui.receiver = gui.draw_head(pos,orientation,R);

            obj.gui_ref = gui;
            draggable(gui.receiver,@update_receiver_position, @update_receiver_orientation);
            function update_receiver_position(receiver)
                obj.receiver.set_position( gui.receiver.UserData.Origin );
                obj.schedule_update('receiver_moved');
            end
            function update_receiver_orientation(receiver)
                obj.receiver.set_orientation(  [cosd(gui.receiver.UserData.Orientation),...
                    sind(gui.receiver.UserData.Orientation)] );
                obj.schedule_update('receiver_rotated');
            end

        end

        function obj = delete_receiver(obj, gui)
            obj.receiver = {};
            delete(gui.receiver);
        end

        function obj = create_virtual_source(obj, position, orientation, gui, source_type, rendering, source_name)
            if nargin < 7
                source_name = [];
            end
            idx = length(obj.virtual_sources) + 1;
            obj.virtual_sources{idx} = virtual_source(idx, position, orientation, source_type, rendering);
            gui.virtual_source_points{idx} = gui.draw_virtual_source(obj.virtual_sources{idx}.position,...
                cart2pol(orientation(1),orientation(2))*180/pi,idx,source_name);
            draggable(gui.virtual_source_points{idx},@update_virtual_position, @update_virtual_orientation);
            obj.face_origin(idx, gui);
            function update_virtual_position(virtual_source)
                k = virtual_source.UserData.Label;
                obj.virtual_sources{k}.set_position(gui.virtual_source_points{k}.UserData.Origin);
                obj.face_origin(k, gui);
                obj.schedule_update('virtual_source_moved', k);
            end
            function update_virtual_orientation(virtual_source)
                % sources always face the origin: a manual turn snaps back
                obj.face_origin(virtual_source.UserData.Label, gui);
            end
        end

        function face_origin(obj, k, gui)
            % Virtual sources always face the origin (model orientation in
            % 3D, the icon in the horizontal plane)
            vs = obj.virtual_sources{k};
            p = vs.position;
            if norm(p) < 1e-9, return; end
            vs.orientation = -p/norm(p);
            if nargin < 3, gui = obj.gui_ref; end
            if ~isempty(gui) && isvalid(gui) && norm(p(1:2)) > 1e-9
                gui.orient_virtual_source(k, atan2d(-p(2), -p(1)));
            end
        end

        function obj = delete_virtual_source(obj, virt_source_idx, gui)
            obj.virtual_sources{virt_source_idx} = {};
            delete(gui.virtual_source_points{virt_source_idx});
        end

        function obj = create_loudspeaker(obj, position, orientation, gui, type, height)
            idx = length(obj.loudspeaker_array) + 1;
            % position may be 1×2; height becomes z-coordinate
            if numel(position) < 3 && nargin >= 6 && ~isempty(height)
                position = [position(:)', height];
            end
            obj.loudspeaker_array{idx} = loudspeaker(idx, position, orientation, type);
            gui.loudspeaker_points{idx} = ...
                gui.draw_loudspeaker(obj.loudspeaker_array{idx}.position,type.R(1),cart2pol(orientation(1),orientation(2))*180/pi,idx);

            draggable(gui.loudspeaker_points{idx},@update_loudspeaker_position, @update_loudspeaker_orientation, @on_loudspeaker_dblclick);
            function update_loudspeaker_position(loudspeaker)
                obj.loudspeaker_array{loudspeaker.UserData.Label}.set_position...
                    (gui.loudspeaker_points{loudspeaker.UserData.Label}.UserData.Origin);
                obj.schedule_update('loudspeaker_moved', loudspeaker.UserData.Label);
            end
            function update_loudspeaker_orientation(loudspeaker)
                obj.loudspeaker_array{loudspeaker.UserData.Label}.set_orientation...
                    ( [ cosd(gui.loudspeaker_points{loudspeaker.UserData.Label}.UserData.Orientation),...
                    sind(gui.loudspeaker_points{loudspeaker.UserData.Label}.UserData.Orientation)]);
                obj.schedule_update('loudspeaker_rotated', loudspeaker.UserData.Label);
            end
            function on_loudspeaker_dblclick(lsObj)
                chIdx = lsObj.UserData.Label;
                if ~isempty(obj.TestChannelFcn)
                    obj.TestChannelFcn(chIdx);
                else
                    fprintf('Double-click on speaker %d (no test callback set)\n', chIdx);
                end
            end
        end

        function obj = delete_loudspeaker(obj, bin_source_idx, gui)
            obj.loudspeaker_array(bin_source_idx) = [];
            delete(gui.loudspeaker_points{bin_source_idx});
        end


        % Reserved for multiple sound scene scenario (for comparing approaches)
        function obj = save_sound_scene(obj)
        end
        function obj = load_sound_scene(obj)
        end

        function duplum = duplicate_scene(obj)
            duplum = obj;
        end

        function setup = save_setup(obj)
            setup.vs_position = cell2mat(cellfun(@(x) x.position, obj.virtual_sources, 'UniformOutput',0)');
            setup.vs_orientation = cell2mat(cellfun(@(x) x.orientation, obj.virtual_sources, 'UniformOutput',0)');
            setup.rec_position = obj.receiver.position;
            setup.rec_orientation = obj.receiver.orientation;
            setup.ls_position = cell2mat(cellfun(@(x) x.position, obj.loudspeaker_array, 'UniformOutput',0)');
            setup.ls_orientation = cell2mat(cellfun(@(x) x.orientation, obj.loudspeaker_array, 'UniformOutput',0)');



        end

        function obj = delete(obj,gui)
            if ~isempty(obj.updateTimer) && isvalid(obj.updateTimer)
                stop(obj.updateTimer); delete(obj.updateTimer);
            end
            obj.delete_receiver(gui);
            N = length(obj.virtual_sources);
            for n = 1 : N
                obj.delete_virtual_source(N-n+1,gui);
            end
            N = length(obj.loudspeaker_array);
            for n = 1 : N
                obj.delete_loudspeaker(N-n+1,gui);
            end
            for n = 1 : length(gui.main_axes.Children)
                delete(gui.main_axes.Children(1))
            end
            clear obj
        end

        function schedule_update(obj, kind, label)
            % Record a geometry change; the renderer update runs in
            % flush_updates within UpdateDelay (one update per burst of
            % mouse-motion events instead of one per event).
            % label: virtual source / loudspeaker index (vs: optional, all if omitted)
            P = obj.pendingUpdates;
            switch kind
                case 'receiver_moved',          P.receiver = 'receiver_moved';
                case 'receiver_rotated'
                    if isempty(P.receiver), P.receiver = 'receiver_rotated'; end
                case 'virtual_source_moved'
                    if nargin < 3 || isempty(label), P.vsAll = true;
                    else, P.vs = union(P.vs, label); end
                case 'virtual_source_rotated',  P.vsRot = true;
                case 'loudspeaker_moved',       P.ls = union(P.ls, label);
                case 'loudspeaker_rotated',     P.lsRot = union(P.lsRot, label);
            end
            obj.pendingUpdates = P;
            obj.arm_update_timer();
        end

        function flush_updates(obj)
            % Apply the recorded geometry changes. Moved virtual sources
            % update only their own renderer, within UpdateBudget seconds
            % per flush (round-robin; the rest is re-armed), so that many
            % moving objects (Atmos playback) do not starve the audio loop.
            P = obj.pendingUpdates;
            obj.pendingUpdates = obj.emptyPending();
            R = obj.scene_renderer;
            if isempty(R), return; end
            t0 = tic;
            if ~isempty(P.receiver)
                R.update_binaural_renderers(P.receiver, []);
                R.update_SFS_renderers(P.receiver);
            end
            targeted = ~any(cellfun(@(r) isa(r, 'beamformer_renderer'), R.SFS_renderer));
            if P.vsAll || (~isempty(P.vs) && ~targeted)
                R.update_SFS_renderers('virtual_source_moved');
            elseif ~isempty(P.vs)
                done = 0;
                for l = P.vs(:)'
                    R.update_SFS_renderer(l, 'virtual_source_moved');
                    R.prepare_fast_path();      % its BRIR too: the budget covers both
                    done = done + 1;
                    if toc(t0) > obj.UpdateBudget, break; end
                end
                rest = P.vs(done+1:end);
                if ~isempty(rest)
                    obj.pendingUpdates.vs = union(obj.pendingUpdates.vs, rest);
                end
            end
            if P.vsRot, R.update_SFS_renderers('virtual_source_rotated'); end
            if ~isempty(P.ls)
                R.update_SFS_renderers('loudspeaker_moved');
                for l = P.ls(:)', R.update_binaural_renderers('loudspeaker_moved', l); end
            end
            if ~isempty(P.lsRot)
                R.update_SFS_renderers('loudspeaker_rotated');
                for l = P.lsRot(:)', R.update_binaural_renderers('loudspeaker_rotated', l); end
            end
            % rebuild the binaural BRIRs here, not inside the next audio block
            R.prepare_fast_path();
            if ~isempty(obj.gui_ref) && isvalid(obj.gui_ref) && ~obj.ObjectPlayback
                for kk = 1:numel(R.SFS_renderer)
                    obj.gui_ref.request_draw_visuals(R.SFS_renderer{kk});
                end
            end
            if ~isempty(obj.OnGeometryChangedFcn) && ~obj.ObjectPlayback
                try obj.OnGeometryChangedFcn(); catch, end
            end
            if ~isempty(obj.pendingUpdates.vs)
                obj.arm_update_timer();             % budget left sources over
            end
        end

        function apply_object_positions(obj, pos, level)
            % Object-based input (cavern_reader): move the virtual sources.
            %   pos   [Nobj x 3] normalised [right, front, up] (-1..1, -1..1, 0..1)
            %   level [Nobj x 1] block RMS, for the marker styling
            % The Atmos room is mapped to a box around the listener: half-width
            % ObjectRoomRadius, ear height (up = 0) to ObjectRoomHeight
            % (up = 1). The sources move in 3D; each renderer handles the
            % height its own way (binaural/VBAP: direction, 2.5D: distance).
            obj.ObjectPlayback = true;
            Rr = obj.ObjectRoomRadius;
            rec = obj.receiver.position(1:3);
            n = min(size(pos,1), numel(obj.virtual_sources));
            moved = [];
            for k = 1:n
                if ~any(pos(k,:))
                    continue;    % Cavern reports [0 0 0] until the object's first metadata
                end
                xyz = rec + [Rr*pos(k,1:2), obj.ObjectRoomHeight*pos(k,3)];
                vs = obj.virtual_sources{k};
                if norm(vs.position(1:3) - xyz) > 0.01      % 1 cm
                    vs.set_position(xyz);
                    moved(end+1) = k; %#ok<AGROW>
                    if ~isempty(obj.gui_ref) && isvalid(obj.gui_ref)
                        obj.gui_ref.move_virtual_source(k, xyz);
                    end
                    obj.face_origin(k);
                end
            end
            obj.styleCounter = obj.styleCounter + 1;
            if ~isempty(obj.gui_ref) && isvalid(obj.gui_ref) && nargin > 2 && mod(obj.styleCounter, 3) == 1
                obj.gui_ref.style_virtual_sources(pos(1:n,3), level(1:n));
            end
            if ~isempty(moved)
                obj.schedule_update('virtual_source_moved', moved);
            end
        end

        function apply_source_positions(obj, idx, P)
            % Move virtual sources idx to rows of P [x y z] (trajectories);
            % same targeted, time-budgeted update as the object positions
            moved = [];
            hasGui = ~isempty(obj.gui_ref) && isvalid(obj.gui_ref);
            for j = 1:numel(idx)
                k = idx(j);
                if k > numel(obj.virtual_sources), continue; end
                vs = obj.virtual_sources{k};
                if norm(vs.position(1:3) - P(j,:)) > 0.005      % 5 mm
                    vs.set_position(P(j,:));
                    if hasGui, obj.gui_ref.move_virtual_source(k, P(j,:)); end
                    obj.face_origin(k);
                    moved(end+1) = k; %#ok<AGROW>
                end
            end
            if ~isempty(moved)
                obj.schedule_update('virtual_source_moved', moved);
            end
        end

        function arm_update_timer(obj)
            if isempty(obj.updateTimer) || ~isvalid(obj.updateTimer)
                obj.updateTimer = timer('ExecutionMode','singleShot', ...
                    'StartDelay', obj.UpdateDelay, 'Name','SceneUpdate', ...
                    'TimerFcn', @(~,~) obj.flush_updates());
            end
            if strcmp(obj.updateTimer.Running, 'off')
                start(obj.updateTimer);
            end
        end

        function output = render_sound_scene(obj,input, Binauralization, Downmixing_enabled, Decorrelation)

            % ── Direct playback bypass ──
            % In direct mode, the input goes straight to output without
            % SSD/renderer processing. This avoids issues with mismatched
            % SSD size vs. soundcard channel count.
            isDirect = false;
            try
                isDirect = strcmpi(obj.scene_renderer.Setup.Rendering_mode, 'direct_playback');
            catch
            end
            test = obj.next_test_block(size(input,1));   % [] when no test signal

            if isDirect && ~Binauralization
                output = input;
                if ~isempty(test) && obj.testLs <= size(output,2)
                    output(:,obj.testLs) = output(:,obj.testLs) + test;
                end
                return;
            end

            % ── Normal rendering path ──
            % Headphone bus starts from zero; the loudspeaker signals are
            % zeroed by the scene renderer, and only when it computes them
            % (the fast binaural path skips them)
            obj.headphone.source_signal.clear_signal;

            % Render (all renderers work in SSD contour order)
            obj.scene_renderer.render(input, Binauralization);
            if ~Binauralization
                % Collect signals in SSD order
                raw = cell2mat(cellfun( @(x) x.source_signal.get_signal, obj.loudspeaker_array, 'UniformOutput',false));
                if ~isempty(test)
                    raw(:,obj.testLs) = raw(:,obj.testLs) + test;
                end

                % Remap SSD order → physical channel order
                if ~isempty(obj.channelMap)
                    Nch = size(raw, 2);
                    Nmap = numel(obj.channelMap);
                    Nout = max(Nmap, max(obj.channelMap));
                    remapped = zeros(size(raw,1), Nout);
                    for k = 1:min(Nch, Nmap)
                        remapped(:, obj.channelMap(k)) = raw(:, k);
                    end
                    raw = remapped;
                end

                if Downmixing_enabled
                    output = raw * obj.downmixing_matrix';
                else
                    output = raw;
                end
            elseif Binauralization
                output = obj.headphone.source_signal.get_signal;
                if ~isempty(test) || any(obj.testTail(:))
                    output = output + obj.binaural_test_block(test, size(output,1));
                end
            end

        end

        function x = next_test_block(obj, n)
            % Next n samples of the loudspeaker test signal ([] if none)
            x = [];
            if isempty(obj.testSig), return; end
            m = min(n, numel(obj.testSig));
            x = [obj.testSig(1:m); zeros(n-m, 1)];
            obj.testSig = obj.testSig(m+1:end);
        end

        function y = binaural_test_block(obj, x, n)
            % Test block through the loudspeaker's two-ear IR (overlap-add)
            if isempty(x), x = zeros(n,1); end
            if isempty(obj.testH)
                full = 0.5*[x x];                 % no IR available: centre
            else
                full = [conv(x, obj.testH(:,1)), conv(x, obj.testH(:,2))];
            end
            t = obj.testTail;
            L = max(size(full,1), size(t,1));
            full(end+1:L, :) = 0;
            if ~isempty(t), full(1:size(t,1), :) = full(1:size(t,1), :) + t; end
            y = full(1:n, :);
            obj.testTail = full(n+1:end, :);
            if ~isempty(obj.testTail) && max(abs(obj.testTail(:))) < 1e-7
                obj.testTail = [];
            end
        end

    end
end
