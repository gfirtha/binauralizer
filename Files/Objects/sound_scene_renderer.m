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

        % HRTF / directivity
        directivity_tables
        HRTF_extrapolators

        % ---- Cached scene references (so we can rebuild on mode change)
        virtual_sources         % cell array of virtual source objects
        loudspeaker_array       % array/cell describing loudspeakers
        headphone
        receiver
        Setup                   % copy of setup struct used to build
    end

    methods
        function obj = sound_scene_renderer(virtual_sources, loudspeaker_array, headphone, receiver, setup)

            % Cache scene references for future rebuilds (set_mode, etc.)
            obj.virtual_sources   = virtual_sources;
            obj.loudspeaker_array = loudspeaker_array;
            obj.headphone         = headphone;
            obj.receiver          = receiver;
            obj.Setup             = setup;

            % Build SSD (used by several renderer types)
            SSD = secondary_source_distribution(loudspeaker_array, [0 0]);

            % Build directivity tables (unique per source type)
            N_fft = 2^nextpow2( min(setup.Block_size + size(setup.HRTF.Data.IR,3), 2*setup.Block_size) - 1 );
            cnt = 0;
            for n = 1:numel(virtual_sources)
                if isempty(get_dirtable_idx(obj.directivity_tables, virtual_sources{n}))
                    cnt = cnt + 1;
                    obj.directivity_tables{cnt} = directivity_table(virtual_sources{n}.source_type, N_fft, setup.SampleRate);
                end
            end
            for m = 1:numel(loudspeaker_array)
                if isempty(get_dirtable_idx(obj.directivity_tables, loudspeaker_array{m}))
                    cnt = cnt + 1;
                    obj.directivity_tables{cnt} = directivity_table(loudspeaker_array{m}.source_type, N_fft, setup.SampleRate);
                end
            end

            % Instantiate SFS renderers per virtual source
            obj.SFS_renderer = cell(1, numel(virtual_sources));
            for n = 1:numel(virtual_sources)
                obj.SFS_renderer{n} = obj.make_renderer_for(virtual_sources{n}, SSD, loudspeaker_array, setup);
            end

            % HRTF extrapolator(s)
            if ndims(setup.HRTF.Data.IR) > 3
                % BRIRs
                obj.HRTF_extrapolators = brir_extrapolator(setup.HRTF);
            else
                obj.HRTF_extrapolators = hrtf_extrapolator(setup.HRTF, setup.HRTF_extrapolation);
            end

            % Binaural renderers per loudspeaker
            obj.binaural_renderer = cell(1, numel(loudspeaker_array));
            for m = 1:numel(loudspeaker_array)
                idx = get_dirtable_idx(obj.directivity_tables, loudspeaker_array{m});
                obj.binaural_renderer{m} = binaural_renderer( ...
                    loudspeaker_array{m}, receiver, headphone, ...
                    obj.directivity_tables{idx}, obj.HRTF_extrapolators);
            end

            % Decorrelator
            obj.decorrelator = decorrelator(setup.Decorrelation, loudspeaker_array, 1);
        end

        % ----- Public API -------------------------------------------------

        function set_mode(obj, modeInt)
            % Switch all virtual sources to the requested renderer mode.
            % Accepts 'WFS' | 'VBAP' | 'HOA' | 'NFC_HOA' | 'Direct playback' |
            % 'Dolby Surround Decoder' or internal ids like 'wfs','vbap',...

            % Update cached Setup (optional, so the value is traceable)
            obj.Setup.Rendering_mode = modeInt;

            % Rebuild SSD (in case geometry changed elsewhere)
            SSD = secondary_source_distribution(obj.loudspeaker_array, [0 0]);

            for n = 1:numel(obj.virtual_sources)
                vs = obj.virtual_sources{n};
                % Update VS renderer_type if present
                if isprop(vs,'renderer_type') || isfield(vs,'renderer_type')
                    try, vs.renderer_type = modeInt; catch, end
                end
                % Rebuild concrete renderer
                obj.SFS_renderer{n} = obj.make_renderer_for(vs, SSD, obj.loudspeaker_array, obj.Setup);
            end
        end

        function update_SFS_renderers(obj, type)
            % For backward-compatibility with your earlier code
            cellfun(@(x) x.update_renderer(type), obj.SFS_renderer);
        end

        function update_renderer_settings(obj, setup)
            % Update settings on all renderers and store the latest copy
            obj.Setup = setup;
            cellfun(@(x) x.update_settings(setup), obj.SFS_renderer);
        end

        function update_binaural_renderers(obj, type, N)
            if isempty(N)
                cellfun(@(x) x.update_renderer(type), obj.binaural_renderer);
            else
                obj.binaural_renderer{N}.update_renderer(type);
            end
        end

        function render(obj, input, binaural_mode)
            % Push one block through all renderers; mix/route assumed external
            for m = 1:numel(obj.SFS_renderer)
                %TODO!
                
                obj.SFS_renderer{m}.virtual_source.source_signal.set_signal(obj.SFS_renderer{m}.virtual_source.gain*input(:,m));
                obj.SFS_renderer{m}.render;
            end
            if obj.decorrelator.mode
                obj.decorrelator.render;
            end
            if binaural_mode
                cellfun(@(x) x.render, obj.binaural_renderer);
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
        function r = make_renderer_for(obj, vs, SSD, loudspeaker_array, setup)
            % Create a concrete renderer instance for a virtual source
            % based on its (normalized) renderer_type.
            % Accepts VS with possible field/property 'renderer_type'.
            rtype = 'direct_playback';
            try
                if isprop(vs,'renderer_type'); rtype = vs.renderer_type; end
            catch
            end
            if isstruct(vs) && isfield(vs,'renderer_type')
                rtype = vs.renderer_type;
            end

            switch lower(rtype)
                case 'direct_playback'
                    r = direct_renderer(vs, SSD);
                case 'vbap'
                    r = vbap_renderer(vs, SSD);
                case 'dbap'
                    r = dbap_renderer(vs, SSD);
                case 'wfs'
                    r = wfs_renderer(vs, SSD, setup.Renderer_setup);
                case 'lwfs_foc'
                    r = lwfs_foc_renderer(vs, loudspeaker_array, setup.SampleRate);
                case 'lwfs_nara'
                    r = lwfs_renderer_nara(vs, loudspeaker_array, setup.SampleRate);
                case 'nara'
                    r = nara_renderer(vs, loudspeaker_array, setup.SampleRate);
                case 'hoa'
                    r = hoa_renderer(vs, SSD, setup.Renderer_setup.HOAorder);
                case 'nfc_hoa'
                    r = nfc_hoa_renderer(vs, SSD, setup.Renderer_setup.HOAorder);
                case 'vbap_wfs'
                    r = vbap_wfs_renderer(vs, loudspeaker_array, setup.SampleRate);
                case 'td_stereo'
                    r = time_delay_renderer(vs, SSD, setup);
                case 'ctc'
                    rs = setup.Renderer_setup;
                    r  = ctc_renderer(vs, SSD, obj.receiver, rs);
                case 'dolby_surround_decoder'
                    r = dolby_surround_renderer(vs, SSD);
                case 'room_simulator'
                    r = room_renderer(vs, obj.receiver, loudspeaker_array, setup.Renderer_setup);
                otherwise
                    % Fallback: direct playback to avoid hard failure
                    warning('Unknown renderer type "%s". Falling back to direct_playback.', rtype);
                    r = direct_renderer(vs, SSD);
            end
        end

  
    end
end
