function iout = get_circular_index(i,v)
    N = length(v);
    if i<0
        iout = N+i+1;
    else 
        iout = i + 1;
    end
end

