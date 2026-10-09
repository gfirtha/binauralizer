function [seg, i0, i1] = mbr_window(h, fs, win)
%MBR_WINDOW  Samples of h in the time window win = [t0 t1] (ms), as the
%thesis scripts select them: i = round(t/1000*fs) + 1, clamped.
i0 = max(1, round(win(1)/1000*fs) + 1);
i1 = min(size(h,1), round(win(2)/1000*fs) + 1);
if i1 - i0 < 10
    error('mbr_window:short', 'Window too short (%d samples).', i1 - i0 + 1);
end
seg = h(i0:i1, :);
end
