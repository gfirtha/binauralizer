function h = lagrange_filter(filter_order, dfrac)
    % Design Lagrange interpolation filter.
    
    h = ones(1, filter_order + 1);
    for k = 1:(filter_order + 1)
        for m = 1:(filter_order + 1)
            if m == k
                continue;
            end
            h(k) = h(k) * (dfrac - m) / (k - m);
        end
    end
end
