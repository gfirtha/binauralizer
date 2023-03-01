classdef sound_scene_renderer < handle
    %RENDERER Summary of this class goes here
    %   Detailed explanation goes here
    % TODO: make input, rendeder out and binauarl buses
    % and "wire" them up
    properties
        SFS_renderer
        binaural_renderer
        directivity_tables
    end

    methods
        function obj = sound_scene_renderer(virtual_sources,loudspeaker_array,headphone, receiver, setup)

            % Get all required directivity characteristics
            N_fft = 2^nextpow2( min(setup.Block_size + size(setup.HRTF.Data.IR,3), 2*setup.Block_size) - 1 );
            cnt = 0;
            for n = 1 : length(virtual_sources)
                if isempty(get_dirtable_idx( obj.directivity_tables, virtual_sources{n}))
                    cnt = cnt + 1;
                    obj.directivity_tables{cnt} = directivity_table(virtual_sources{n}.source_type, N_fft, setup.SampleRate);
                end
            end
            for m = 1 : length(loudspeaker_array)
                if isempty(get_dirtable_idx( obj.directivity_tables, loudspeaker_array{m}))
                    cnt = cnt + 1;
                    obj.directivity_tables{cnt} = directivity_table(loudspeaker_array{m}.source_type, N_fft, setup.SampleRate);
                end
            end

            for n = 1 : length(virtual_sources)
                idx = get_dirtable_idx(obj.directivity_tables,virtual_sources{n});
                switch virtual_sources{n}.renderer_type
                    case 'Direct_playback'
                        obj.SFS_renderer{n} = direct_renderer(virtual_sources{n}, loudspeaker_array);
                    case 'VBAP'
                        obj.SFS_renderer{n} = vbap_renderer(virtual_sources{n}, loudspeaker_array);
                    case 'DBAP'
                        obj.SFS_renderer{n} = dbap_renderer(virtual_sources{n}, loudspeaker_array);
                    case 'WFS'
                        obj.SFS_renderer{n} = wfs_renderer(virtual_sources{n}, loudspeaker_array, setup.SampleRate,obj.directivity_tables{idx}, setup.Renderer_setup.Antialiasing);
                    case 'VBAP_WFS'
                        obj.SFS_renderer{n} = vbap_wfs_renderer(virtual_sources{n}, loudspeaker_array, setup.SampleRate );
                    case 'TD_stereo'
                        obj.SFS_renderer{n} = time_delay_renderer(virtual_sources{n}, loudspeaker_array, setup.SampleRate);
                    case 'CTC'
                        obj.SFS_renderer{n} = ctc_renderer(virtual_sources{n}, loudspeaker_array, receiver, setup.SampleRate, setup.Renderer_setup.Plant_model, setup.Renderer_setup.VS_model, setup.Renderer_setup.HRTF_database,setup.Renderer_setup.N_filt);
                end
                for n = 1 : length(loudspeaker_array)
                    idx = get_dirtable_idx(obj.directivity_tables,loudspeaker_array{n});
                    obj.binaural_renderer{n} = binaural_renderer(loudspeaker_array{n}, receiver,headphone, obj.directivity_tables{idx});
                end
            end
        end

        function update_SFS_renderers(obj, type)
            cellfun( @(x) x.update_renderer(type), obj.SFS_renderer);
        end

        function update_binaural_renderers(obj, type)
            cellfun( @(x) x.update_renderer(type), obj.binaural_renderer);
        end

        function render(obj, input, binaural_mode)
            for m = 1 : length(obj.SFS_renderer)
                obj.SFS_renderer{m}.virtual_source.source_signal.set_signal(input(:,m));
                obj.SFS_renderer{m}.render;
            end
            if binaural_mode
                cellfun( @(x) x.render, obj.binaural_renderer);
            end
        end
    end
end