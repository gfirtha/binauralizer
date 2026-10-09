classdef sound_scene_renderer < handle
    % sound_scene_renderer
    % Keeps renderer stack for each virtual source and supports switching
    % renderer "modes" using either UI names (e.g., 'NFC_HOA') or internal
    % ids (e.g., 'nfc_hoa') via set_mode(...).

    properties
        % Renderers & processing
        SFS_renderer            % cell array of per-source renderers
        decorrelator
        binaural_renderer       % cell array, one per loudspeaker

        % HRTF
        HRTF_extrapolators

        % ---- Cached scene references (so we can rebuild on mode change)
        virtual_sources         % cell array of virtual source objects
        ssd                     % secondary_source_distribution (owns loudspeakers + geometry)
        headphone
        receiver
        Setup                   % copy of setup struct used to build

        % ---- Fast binaural path
        % For LTI renderers the chain VS -> driving filters -> loudspeakers
        % -> HRTFs collapses to one two-ear "system BRIR" per virtual
        % source, rebuilt only when a filter changes. Loudspeaker signals
        % are then NOT computed; set KeepLoudspeakerSignals (e.g. live
        % sound-field monitor) or FastBinaural = false to force the
        % per-loudspeaker chain.
        FastBinaural           = true
        KeepLoudspeakerSignals = false
    end

    properties (Access = private)
        fast_conv     = {}      % cell{Nvs} of brir_convolver
        fast_versions = []      % SFS filter_version snapshot the BRIRs were built from
        fast_bin_version = -1   % binaural_version the cached spectra belong to
        fast_B        = []      % cached loudspeaker two-ear spectra
        fast_ok       = false   % last build succeeded
        fast_active   = false   % fast path produced the previous block
        binaural_version = 0    % bumped by update_binaural_renderers
    end

    methods
        function obj = sound_scene_renderer(virtual_sources, loudspeaker_array, headphone, receiver, setup)

            % Cache scene references for future rebuilds (set_mode, etc.)
            obj.virtual_sources   = virtual_sources;
            obj.headphone         = headphone;
            obj.receiver          = receiver;
            obj.Setup             = setup;

            % Build SSD — the single owner of loudspeaker geometry + objects
            obj.ssd = secondary_source_distribution(loudspeaker_array, [0 0]);

            % Instantiate SFS renderers per virtual source
            obj.SFS_renderer = cell(1, numel(virtual_sources));
            for n = 1:numel(virtual_sources)
                obj.SFS_renderer{n} = obj.make_renderer_for(virtual_sources{n}, setup);
            end

            % HRTF extrapolator(s)
            if ndims(setup.HRTF.Data.IR) > 3
                obj.HRTF_extrapolators = brir_extrapolator(setup.HRTF);
            else
                obj.HRTF_extrapolators = hrtf_extrapolator(setup.HRTF, setup.HRTF_extrapolation);
            end

            % Wire CTC support into beamformer renderers.
            % The first beamformer instance is marked as "primary" — it
            % iterates over all scene VSs in CTC mode. The others are
            % no-ops in CTC mode to avoid double-summing into the SSD.
            firstBF = true;
            for n = 1:numel(obj.SFS_renderer)
                if isa(obj.SFS_renderer{n}, 'beamformer_renderer')
                    obj.SFS_renderer{n}.scene_virtual_sources = virtual_sources;
                    obj.SFS_renderer{n}.hrtf_extrap           = obj.HRTF_extrapolators;
                    obj.SFS_renderer{n}.ctc_primary           = firstBF;
                    firstBF = false;
                end
            end

            % Binaural renderers per loudspeaker
            % Each loudspeaker now owns its own Directivity object.
            obj.binaural_renderer = cell(1, obj.ssd.N);
            for m = 1:obj.ssd.N
                obj.binaural_renderer{m} = binaural_renderer( ...
                    obj.ssd.loudspeakers{m}, receiver, headphone, ...
                    [], obj.HRTF_extrapolators);
                % Note: directivity_table argument is [] — binaural_renderer
                % uses source.directivity (new path) automatically.
            end

            % Decorrelator
            obj.decorrelator = decorrelator(setup.Decorrelation, obj.ssd.loudspeakers, 1);
        end

        % ----- Public API -------------------------------------------------

        function set_mode(obj, modeInt)
            % Switch all virtual sources to the requested renderer mode.
            obj.Setup.Rendering_mode = modeInt;

            % Rebuild SSD geometry (in case positions changed)
            obj.ssd.update_geometry(obj.ssd.center);

            for n = 1:numel(obj.virtual_sources)
                vs = obj.virtual_sources{n};
                if isprop(vs,'renderer_type') || isfield(vs,'renderer_type')
                    try, vs.renderer_type = modeInt; catch, end
                end
                obj.SFS_renderer{n} = obj.make_renderer_for(vs, obj.Setup);
            end
            obj.fast_conv = {};          % new renderers: rebuild from scratch
            obj.fast_versions = [];
            obj.fast_B = [];
        end

        function update_SFS_renderers(obj, type)
            % Refresh SSD geometry from loudspeaker handle positions (only
            % loudspeaker edits change it; classifying it costs an SVD +
            % boundary() per call, i.e. per drag update)
            if nargin < 2 || startsWith(string(type), 'loudspeaker')
                obj.ssd.update_geometry();
            end
            % Then update all renderers
            cellfun(@(x) x.update_renderer(type), obj.SFS_renderer);
        end

        function update_SFS_renderer(obj, idx, type)
            % Update one source's renderer (a moved virtual source affects
            % only its own driving functions)
            obj.SFS_renderer{idx}.update_renderer(type);
        end

        function update_renderer_settings(obj, setup)
            % Update settings on all renderers and store the latest copy
            obj.Setup = setup;
            cellfun(@(x) x.update_settings(setup.Renderer_setup), obj.SFS_renderer);
        end

        function update_binaural_renderers(obj, type, N)
            if isempty(N)
                cellfun(@(x) x.update_renderer(type), obj.binaural_renderer);
            else
                obj.binaural_renderer{N}.update_renderer(type);
            end
            obj.binaural_version = obj.binaural_version + 1;
        end

        function render(obj, input, binaural_mode)
            % First push this block's signal into every virtual source.
            % This must happen before any renderer.render() runs, so that
            % multi-VS renderers (e.g. beamformer in CTC mode) see all
            % current-block signals when they iterate over scene VSs.
            for m = 1:numel(obj.SFS_renderer)
                obj.SFS_renderer{m}.virtual_source.source_signal.set_signal( ...
                    obj.SFS_renderer{m}.virtual_source.gain * input(:,m));
            end

            if binaural_mode && obj.use_fast_path()
                y = 0;
                for m = 1:numel(obj.SFS_renderer)
                    y = y + obj.fast_conv{m}.process( ...
                        obj.SFS_renderer{m}.virtual_source.source_signal.time_series);
                end
                obj.headphone.source_signal.set_signal(y);
                if ~obj.fast_active
                    obj.ssd.clear_signals;   % no stale loudspeaker signals
                    obj.fast_active = true;
                end
                return;
            end
            if obj.fast_active
                % back on the per-loudspeaker chain: the headphone bus is
                % summed as spectra again, drop the fast path's time signal
                obj.headphone.source_signal.time_series = [];
                obj.fast_active = false;
            end
            obj.ssd.clear_signals;

            % Then render each renderer in turn.
            for m = 1:numel(obj.SFS_renderer)
                obj.SFS_renderer{m}.render;
            end
            if (obj.Setup.Decorrelation ~= 1)
                obj.decorrelator.render;
            end
            if binaural_mode
                cellfun(@(x) x.render, obj.binaural_renderer);
            end
        end

        function prepare_fast_path(obj)
            % Rebuild stale system BRIRs now (e.g. right after a scene
            % update) so the next audio block does not pay for it.
            if obj.FastBinaural && isfield(obj.Setup, 'Binauralization') && obj.Setup.Binauralization
                obj.use_fast_path();
            end
        end

        function info = get_renderer_info(obj)
            % Convenience: ask the first renderer for info if it supports it
            info = '(no renderer info available)';
            if ~isempty(obj.SFS_renderer) && ismethod(obj.SFS_renderer{1}, 'get_renderer_info')
                try
                    info = obj.SFS_renderer{1}.get_renderer_info();
                catch ME
                    info = sprintf('Error in get_renderer_info(): %s', ME.message);
                end
            end
        end
    end

    methods (Access = private)
        function tf = use_fast_path(obj)
            tf = false;
            if ~obj.FastBinaural || obj.KeepLoudspeakerSignals || obj.Setup.Decorrelation ~= 1
                return;
            end
            % binaural filters only change via update_binaural_renderers,
            % which bumps binaural_version (polling 128 objects costs ~1 ms)
            sfs = cellfun(@(r) r.filter_version, obj.SFS_renderer);
            binChanged = obj.binaural_version ~= obj.fast_bin_version;
            if numel(sfs) ~= numel(obj.fast_versions)
                stale = true(size(sfs));
            else
                stale = binChanged | sfs ~= obj.fast_versions;
            end
            if any(stale)
                % only the sources whose filters changed are rebuilt
                obj.fast_ok = obj.build_system_brirs(find(stale), binChanged);
                obj.fast_versions = sfs;
                obj.fast_bin_version = obj.binaural_version;
            end
            tf = obj.fast_ok;
        end

        function ok = build_system_brirs(obj, which, binChanged)
            % H_m(w) = sum_n D_mn(w) * B_n(w): driving filter of loudspeaker n
            % for source m times its two-ear transfer (HRTF x directivity).
            ok = false;
            Nvs = numel(obj.SFS_renderer);
            Nls = numel(obj.binaural_renderer);
            d = cell(1, Nvs);
            for m = which(:)'
                d{m} = obj.SFS_renderer{m}.get_driving_filters();
                if isempty(d{m}) || size(d{m},2) ~= Nls, return; end
                d{m} = obj.finite_or_silent(d{m}, m);
            end
            Nd = max(cellfun(@(x) size(x,1), d(which)));

            % Loudspeaker two-ear spectra, cached until the HRTFs change
            if binChanged || isempty(obj.fast_B) || Nd ~= obj.fast_B.Nd
                hp = obj.binaural_renderer{1}.get_packed_ir();
                if isempty(hp), return; end
                HP = zeros(numel(hp), Nls);
                for n = 1:Nls
                    HP(:,n) = obj.binaural_renderer{n}.get_packed_ir();
                end
                % support over all loudspeakers (a directivity prefilter can
                % spread the IR over the convolver's whole FFT length)
                Nh = find(any(HP ~= 0, 2), 1, 'last');
                Nc = 2^nextpow2(Nd + Nh - 1);      % no circular wrap
                k1 = 1:Nc/2+1;                     % real IRs: one-sided spectra suffice
                BL = fft(real(HP(1:Nh,:)), Nc, 1); BR = fft(imag(HP(1:Nh,:)), Nc, 1);
                obj.fast_B = struct('Nd', Nd, 'Nh', Nh, 'Nc', Nc, 'BL', BL(k1,:), 'BR', BR(k1,:));
                which = 1:Nvs;                     % every source needs the new HRTFs
                for m = which
                    if isempty(d{m})
                        d{m} = obj.SFS_renderer{m}.get_driving_filters();
                        if isempty(d{m}) || size(d{m},2) ~= Nls, return; end
                        d{m} = obj.finite_or_silent(d{m}, m);
                    end
                end
            end
            B = obj.fast_B;

            blk = numel(obj.SFS_renderer{1}.virtual_source.source_signal.time_series);
            if numel(obj.fast_conv) ~= Nvs
                obj.fast_conv = cell(1, Nvs);
                for m = 1:Nvs, obj.fast_conv{m} = brir_convolver(blk); end
            end
            for m = which(:)'
                Fd = fft(d{m}, B.Nc, 1);
                Fd = Fd(1:B.Nc/2+1, :);
                h = ifft([sum(Fd.*B.BL, 2), sum(Fd.*B.BR, 2)], B.Nc, 1, 'symmetric');
                obj.fast_conv{m}.set_filter(h(1:size(d{m},1)+B.Nh-1, :));
            end
            ok = true;
        end

        function d = finite_or_silent(~, d, m)
            % One non-finite driving filter would turn every source's
            % output into NaN once summed: silence that source instead.
            if ~all(isfinite(d(:)))
                warning('sound_scene_renderer:nonFinite', ...
                    'Virtual source %d has non-finite driving functions - muted.', m);
                d(:) = 0;
            end
        end

        function r = make_renderer_for(obj, vs, setup)
            % Create a concrete renderer via the central RendererRegistry.
            rtype = 'direct_playback';
            try
                if isprop(vs,'renderer_type'); rtype = vs.renderer_type; end
            catch
            end
            if isstruct(vs) && isfield(vs,'renderer_type')
                rtype = vs.renderer_type;
            end

            reg = audioapp.RendererRegistry.instance();
            r = reg.create(rtype, vs, obj.ssd, obj.ssd.loudspeakers, obj.receiver, setup);
            if isa(r, 'wfs_renderer')
                r.LazyFilterBank = true;   % fast path needs the spectra only
            end
        end

  
    end
end
