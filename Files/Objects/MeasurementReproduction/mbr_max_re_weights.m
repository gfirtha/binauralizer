function a_n = mbr_max_re_weights(N)
%MBR_MAX_RE_WEIGHTS  max-rE weights per ACN channel [1 x (N+1)^2]
%(port of get_max_rE_weights, M. Bubrik's thesis code).
n = repelem(0:N, 2*(0:N) + 1);
a_n = zeros(size(n));
x = cos(2.4062/(N + 1.51));
for k = 0:N
    a_n(n == k) = legendreP(k, x);
end
end
