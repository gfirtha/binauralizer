function sep = make_separator(cfg)
%MAKE_SEPARATOR  Separator from a config struct (Setup.Separation).
%   cfg.Method: 'passthrough' | 'matrix' | 'stft' (fixed slots) |
%               'stft_objects' (tracked, moving objects) |
%               'avendano' (Avendano-Jot coherence / panning index, slots) |
%               'avendano_ring' (5.1 four-pair ring, C excluded from pairs) |
%               'avendano_5ch' (5.1 five-pair bank, C included) |
%               'soundfield' (sound-field sectors + tracking: moving objects + bed) |
%               'wavenumber' (local wavenumber + wavefront curvature at the listener,
%                             stereo first Lt/Rt-decoded; objects with distance) |
%               'phantom_objects' (phantom/ICA objects + fast tracks + music bed) |
%               'phantom_traj' (coherent phantom-pan transient trajectories + music bed) |
%               'ltrt' (Lt/Rt matrix decoder only: L C R Ls Rs) |
%               'cluster' (TF bins clustered on their wavenumber direction ->
%                          sources, each with its own wavenumber vector)
%   stft, avendano: cfg.K, cfg.Width, cfg.Base, cfg.AmbAz, cfg.AmbGain, cfg.Tau

switch lower(cfg.Method)
    case 'passthrough'
        sep = passthrough_separator();
    case 'matrix'
        sep = matrix_separator();
    case {'stft', 'stft_objects'}
        sep = stft_separator();
        if strcmpi(cfg.Method, 'stft_objects'), sep.Mode = 'objects'; end
        f = {'K', 'Width', 'Base', 'AmbAz', 'AmbGain', 'Tau'};
        for k = 1:numel(f)
            if isfield(cfg, f{k}), sep.(f{k}) = cfg.(f{k}); end
        end
    case 'avendano'
        sep = avendano_separator();
        f = {'K', 'Width', 'Base', 'AmbAz', 'AmbGain', 'Tau'};
        for k = 1:numel(f)
            if isfield(cfg, f{k}), sep.(f{k}) = cfg.(f{k}); end
        end
    case 'avendano_ring'
        sep = avendano_ring_separator();
        f = {'K', 'Width', 'Base', 'AmbAz', 'AmbGain', 'Tau'};
        for k = 1:numel(f)
            if isfield(cfg, f{k}), sep.(f{k}) = cfg.(f{k}); end
        end
    case 'avendano_5ch'
        sep = avendano_5ch_separator();
        f = {'K', 'Width', 'Base', 'AmbAz', 'AmbGain', 'Tau'};
        for k = 1:numel(f)
            if isfield(cfg, f{k}), sep.(f{k}) = cfg.(f{k}); end
        end
    case 'wavenumber'
        sep = wavenumber_separator();
        f = {'K', 'AmbGain'};
        for k = 1:numel(f)
            if isfield(cfg, f{k}), sep.(f{k}) = cfg.(f{k}); end
        end
    case 'phantom_objects'
        sep = phantom_object_separator();
        f = {'K', 'AmbGain', 'Tau'};
        for k = 1:numel(f)
            if isfield(cfg, f{k}), sep.(f{k}) = cfg.(f{k}); end
        end
    case 'phantom_traj'
        sep = phantom_trajectory_separator();
        f = {'K', 'AmbGain', 'Tau'};
        for k = 1:numel(f)
            if isfield(cfg, f{k}), sep.(f{k}) = cfg.(f{k}); end
        end
    case 'cluster'
        sep = cluster_separator();
        f = {'K', 'AmbGain'};
        for k = 1:numel(f)
            if isfield(cfg, f{k}), sep.(f{k}) = cfg.(f{k}); end
        end
    case 'ltrt'
        sep = ltrt_decoder();
    case 'soundfield'
        sep = soundfield_separator();
        f = {'K', 'AmbGain', 'AmbAz'};
        for k = 1:numel(f)
            if isfield(cfg, f{k}), sep.(f{k}) = cfg.(f{k}); end
        end
    otherwise
        error('make_separator:method', 'Unknown separation method "%s".', cfg.Method);
end
end
