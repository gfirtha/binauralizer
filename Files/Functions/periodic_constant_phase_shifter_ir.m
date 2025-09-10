function h = periodic_constant_phase_shifter_ir(Ndft, phase_angle)
    n = 0:(Ndft-1);
    h = zeros(1, Ndft);

    if mod(Ndft, 2) == 0
        n_odd = n(mod(n, 2) == 1);
        h(mod(n, 2) == 1) = 2 / Ndft ./ tan(pi * n_odd / Ndft);
    elseif mod(Ndft, 2) == 1
        n_odd = n(mod(n, 2) == 1);
        n_even_nonzero = n(mod(n, 2) == 0 & n ~= 0);
        h(mod(n, 2) == 1) = 1 / Ndft ./ tan(pi * n_odd / 2 / Ndft);
        h(mod(n, 2) == 0 & n ~= 0) = ...
            1 / Ndft ./ tan(pi * (n_even_nonzero + Ndft) / 2 / Ndft);
    end
    h = -sin(phase_angle) * h;
    h(1) = cos(phase_angle);
end
