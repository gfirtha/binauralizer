
function max_order = max_order_circular_harmonics(N, max_order)
    % Compute order of 2D HOA.
    if isempty(max_order)
        max_order = fix(N / 2);
    end
end
