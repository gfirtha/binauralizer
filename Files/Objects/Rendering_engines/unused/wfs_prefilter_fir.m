function [h, delay] = wfs_prefilter_fir(dim, N, fl, fu, fs, c)
    % Create pre-equalization filter for WFS.
    % Rising slope with 3dB/oct ('2.5D') or 6dB/oct ('2D' and '3D').
    % Constant magnitude below fl and above fu.
    % Type 1 linear phase FIR filter of order N.
    % Simple design via "frequency sampling method".
    
    if mod(N, 2)
        error('N must be an even int.');
    end

    bins = floor(N / 2 + 1);
    delta_f = fs / (2 * bins - 1);
    f = (0:(bins-1)) * delta_f;
    if strcmp(dim, '2D') || strcmp(dim, '3D')
        alpha = 1;
    elseif strcmp(dim, '2.5D')
        alpha = 0.5;
    else
        error('Invalid dim value. Must be ''2D'', ''2.5D'' or ''3D''.');
    end
    desired = (2 * pi * f / c) .^ alpha;
    low_shelf = (2 * pi * fl / c) ^ alpha;
    high_shelf = (2 * pi * fu / c) ^ alpha;
    desired = max(desired, low_shelf);
    desired = min(desired, high_shelf);
    desired = desired / desired(end);

    h = ifft(desired, 2*bins - 1, 'symmetric');
    h = circshift(h, [0, bins - 1]);

    E = pi * fs / (2 * c);
    h = h * sqrt(E) / norm(h);

    delay = (bins - 1) / fs;
end
