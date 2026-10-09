function sep = make_separator(cfg)
%MAKE_SEPARATOR  Separator from a config struct (Setup.Separation).
%   cfg.Method: 'passthrough' | 'ltrt' (Lt/Rt matrix decoder only: L C R Ls Rs) |
%               'stft' (PCA primary/ambient, fixed slots) |
%               'stft_objects' (tracked, moving objects) |
%               'avendano' (Avendano-Jot coherence / panning index, slots) |
%               'soundfield' (sound-field sectors + tracking: moving objects + bed) |
%               'wavenumber' (local wavenumber + wavefront curvature at the listener,
%                             stereo first Lt/Rt-decoded; objects with distance) |
%               'cluster' (TF bins clustered on their wavenumber direction ->
%                          sources, each with its own wavenumber vector)
%   Legacy (Files/Objects/Separation/legacy, not offered in the UI, see
%   legacy/README.md): 'matrix', 'avendano_ring', 'avendano_5ch',
%   'phantom_objects', 'phantom_traj'.
%   stft, avendano: cfg.K, cfg.Width, cfg.Base, cfg.AmbAz, cfg.AmbGain, cfg.Tau
%   A description of every method: separator_info(method).

switch lower(cfg.Method)
    case 'passthrough'
        sep = passthrough_separator();
    case 'ltrt'
        sep = ltrt_decoder();
    case {'stft', 'stft_objects'}
        sep = stft_separator();
        if strcmpi(cfg.Method, 'stft_objects'), sep.Mode = 'objects'; end
        sep = set_fields(sep, cfg, {'K', 'Width', 'Base', 'AmbAz', 'AmbGain', 'Tau'});
    case 'avendano'
        sep = set_fields(avendano_separator(), cfg, {'K', 'Width', 'Base', 'AmbAz', 'AmbGain', 'Tau'});
    case 'soundfield'
        sep = set_fields(soundfield_separator(), cfg, {'K', 'AmbGain', 'AmbAz'});
    case 'wavenumber'
        sep = set_fields(wavenumber_separator(), cfg, {'K', 'AmbGain'});
    case 'cluster'
        sep = set_fields(cluster_separator(), cfg, {'K', 'AmbGain'});
    % ---- legacy ----
    case 'matrix'
        sep = matrix_separator();
    case 'avendano_ring'
        sep = set_fields(avendano_ring_separator(), cfg, {'K', 'Width', 'Base', 'AmbAz', 'AmbGain', 'Tau'});
    case 'avendano_5ch'
        sep = set_fields(avendano_5ch_separator(), cfg, {'K', 'Width', 'Base', 'AmbAz', 'AmbGain', 'Tau'});
    case 'phantom_objects'
        sep = set_fields(phantom_object_separator(), cfg, {'K', 'AmbGain', 'Tau'});
    case 'phantom_traj'
        sep = set_fields(phantom_trajectory_separator(), cfg, {'K', 'AmbGain', 'Tau'});
    otherwise
        error('make_separator:method', 'Unknown separation method "%s".', cfg.Method);
end
end

function sep = set_fields(sep, cfg, f)
for k = 1:numel(f)
    if isfield(cfg, f{k}), sep.(f{k}) = cfg.(f{k}); end
end
end
