function y = righthalf_kaiser(n, beta)
    y = kaiser(2*n-1, beta);
    y = y(n:end);
end