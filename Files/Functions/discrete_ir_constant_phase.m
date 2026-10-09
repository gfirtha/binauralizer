function h = discrete_ir_constant_phase(n, phase_angle)
    idx_zero = (n == 0);
    idx_odd = ~idx_zero & mod(n, 2) == 1;
    h = zeros(size(n));
    h(idx_zero) = cos(phase_angle);
    h(idx_odd) = -2 / pi ./ n(idx_odd) .* sin(phase_angle);
end
