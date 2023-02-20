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
        function obj = sound_scene_renderer(virtual_sources,binaural_sources,receiver, setup)
            
            % Get all required directivity characteristics
            N_fft = 2^nextpow2( min(setup.Block_size + size(setup.HRTF.Data.IR,3), 2*setup.Block_size) - 1 );
            cnt = 0;
            for n = 1 : length(virtual_sources)
                if isempty(get_dirtable_idx( obj.directivity_tables, virtual_sources{n}))
                    cnt = cnt + 1;
                    obj.directivity_tables{cnt} = directivity_table(virtual_sources{n}.source_type, N_fft, setup.Input_stream.SampleRate);
                end
            end
            for m = 1 : length(binaural_sources)
                if isempty(get_dirtable_idx( obj.directivity_tables, binaural_sources{m}))
                    cnt = cnt + 1;
                    obj.directivity_tables{cnt} = directivity_table(binaural_sources{m}.source_type, N_fft, setup.Input_stream.SampleRate);
                end
            end
            
            if isempty(virtual_sources) % direct binauralization scenario
                for n = 1 : length(binaural_sources)
                    idx = get_dirtable_idx(obj.directivity_tables,binaural_sources{n});
                    obj.binaural_renderer{n} = binaural_renderer(binaural_sources{n}, receiver,obj.directivity_tables{idx});
                end
            else % virtual sound field synthesis scenario
                for n = 1 : length(virtual_sources)
                    idx = get_dirtable_idx(obj.directivity_tables,virtual_sources{n});
                    switch virtual_sources{n}.renderer_type
                        case 'VBAP'
                            obj.SFS_renderer{n} = vbap_renderer(virtual_sources{n}, binaural_sources);
                        case 'DBAP'                       
                            obj.SFS_renderer{n} = dbap_renderer(virtual_sources{n}, binaural_sources);
                        case 'WFS'
                            obj.SFS_renderer{n} = wfs_renderer(virtual_sources{n}, binaural_sources, setup.Input_stream.SampleRate,obj.directivity_tables{idx}, setup.renderer_setup.Antialiasing);
                        case 'VBAP_WFS'
                            obj.SFS_renderer{n} = vbap_wfs_renderer(virtual_sources{n}, binaural_sources, setup.Input_stream.SampleRate );
                        case 'TD_stereo'
                            obj.SFS_renderer{n} = time_delay_renderer(virtual_sources{n}, binaural_sources, setup.Input_stream.SampleRate);
                    end
                end
                for n = 1 : length(binaural_sources)
                    idx = get_dirtable_idx(obj.directivity_tables,binaural_sources{n});
                    obj.binaural_renderer{n} = binaural_renderer(binaural_sources{n}, receiver,obj.directivity_tables{idx});
                end
            end
        end
        
        function update_SFS_renderers(obj, idx)
            obj.SFS_renderer{idx}.update_renderer;
        end
        
        function update_binaural_renderers(obj, idx, type)
            switch type
                case 'receiver_moved'
                    obj.binaural_renderer{idx}.update_hrtf;
                    obj.binaural_renderer{idx}.update_directivity;
                case 'receiver_rotated'
                    obj.binaural_renderer{idx}.update_hrtf;
                case 'source_moved'
                    obj.binaural_renderer{idx}.update_hrtf;
                    obj.binaural_renderer{idx}.update_directivity;
                case 'source_rotated'
                    obj.binaural_renderer{idx}.update_directivity;
            end
        end
        
        function output = render(obj, input)
            %% Only binauralization job
            if (isempty(obj.SFS_renderer))
                output_signal = signal;
                for n = 1 : length(obj.binaural_renderer)
                    obj.binaural_renderer{n}.binaural_source.source_signal.set_signal(input(:,n));
                    obj.binaural_renderer{n}.render;
                    output_signal.add_spectra(obj.binaural_renderer{n}.output_signal);
                end
                output = output_signal.get_signal;
                %% Sound field synthesis job
            else
                SFS_output = 0;
                for m = 1 : length(obj.SFS_renderer)
                    obj.SFS_renderer{m}.virtual_source.source_signal.set_signal(input(:,m));
                    obj.SFS_renderer{m}.render;
                    
                    SFS_output = SFS_output + cell2mat(cellfun( @(x) x.time_series,...
                        obj.SFS_renderer{m}.output_signal , 'UniformOutput', false));
                end
                output_signal = signal;
                for n = 1 : length(obj.binaural_renderer)
                    obj.binaural_renderer{n}.binaural_source.source_signal.set_signal(SFS_output(:,n));
                    obj.binaural_renderer{n}.render;
                    output_signal.add_spectra(obj.binaural_renderer{n}.output_signal);
                end
                output = output_signal.get_signal;
            end
            
        end
    end
end