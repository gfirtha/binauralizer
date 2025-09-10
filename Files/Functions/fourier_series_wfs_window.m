function [orders, Am] = fourier_series_wfs_window(max_order, phin)
% Eq (50)
    orders = (-max_order:max_order)';
    idx_ones = abs(orders) == 1;
    idx_even = mod(orders, 2) == 0;
    
    Am = zeros(size(orders));
    Am(idx_ones) = 1 / 4;
    Am(idx_even) = (-1).^(orders(idx_even)/2) ./ (pi * (1 - orders(idx_even).^2));
    Am = Am .* exp(-1i * orders * phin);
end