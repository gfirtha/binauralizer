% MAIN_SOUNDFIELD_MAP  Direction-over-time map of a channel mix from the
% virtual sound field at the listener (soundfield_analyzer: intensity
% direction + diffuseness, transient-adaptive block switching): sector
% analysis (several simultaneous directions per band) vs one intensity
% vector, both above 1.5 kHz. Saves soundfield_map_<name>_<t0>s.png.
addpath(genpath('Files'));

file = 'Data/SoundSamples/Eli_50.wav';
segs = [11.45 0.6; 23.0 0.6; 11 3; 22.9 3];   % excerpts [start length] [s]: bullets at 11 / 23 s, zoom + tail
be = logspace(log10(180), log10(16000), 41);   % 40 bands

info = audioinfo(file); fs = info.SampleRate;
[~, nm] = fileparts(file);
for s = 1:size(segs, 1)
t0 = segs(s,1); len = segs(s,2);
x = audioread(file, [round(t0*fs) + 1, min(info.TotalSamples, round((t0 + len)*fs))]);

As = soundfield_analyzer('BandEdges', be, 'Mode', 'sectors', 'HistBand', [1500 Inf]).analyze(x, fs);
Al = soundfield_analyzer('BandEdges', be, 'HistBand', [1500 Inf]).analyze(x, fs);
Tr = soundfield_tracker('Band', [1500 Inf]).track(As);
Tr = Tr(arrayfun(@(s) s.t(end) - s.t(1) >= 0.03, Tr));      % drop very short tracks

fig = figure('Position', [50 50 1400 950], 'Color', 'w');
tl = tiledlayout(fig, 4, 1, 'TileSpacing', 'compact');
az = (As.AzEdges(1:end-1) + As.AzEdges(2:end))/2;
show = @(ax, A, ttl, norm) map(ax, A, az, t0, ttl, norm);

ax1 = nexttile(tl); show(ax1, As, 'Sectors: direct power [dB], > 1.5 kHz, with tracks', false);
hold(ax1, 'on'); cols = lines(max(1, numel(Tr)));
for k = 1:numel(Tr)
    a = Tr(k).az; a(~Tr(k).hit) = NaN;
    a(abs(diff([a(1) a])) > 180) = NaN;                       % do not draw the wrap-around
    plot(ax1, t0 + Tr(k).t, a, '-', 'Color', cols(k,:), 'LineWidth', 1.5);
    text(ax1, t0 + Tr(k).t(1), Tr(k).az(1), sprintf(' %d', Tr(k).id), 'Color', 'w', 'FontSize', 7);
end
hold(ax1, 'off');
ax2 = nexttile(tl); show(ax2, As, 'Sectors: per-column normalised (red ticks: short blocks)', true);
g = find(As.gShort); hold(ax2, 'on');
plot(ax2, t0 + (g - 0.5)*1024/fs, 175*ones(size(g)), 'r|', 'MarkerSize', 6); hold(ax2, 'off');
ax3 = nexttile(tl); show(ax3, Al, 'One intensity vector: per-column normalised', true);
ax4 = nexttile(tl);
fr = round(0.01*fs); nf = floor(size(x,1)/fr);
lev = 10*log10(squeeze(mean(reshape(x(1:nf*fr,:).^2, fr, nf, []), 1)) + 1e-10);
plot(ax4, t0 + ((1:nf) - 0.5)*fr/fs, lev); hold(ax4, 'on');
yyaxis(ax4, 'right'); plot(ax4, t0 + As.tc, As.psiCol, 'k'); ylabel(ax4, 'diffuseness'); ylim(ax4, [0 1]);
yyaxis(ax4, 'left'); ylabel(ax4, 'channel level [dB]');
[~, names] = separator_base.channel_layout(size(x,2));
legend(ax4, [names, {'diffuseness'}], 'Location', 'eastoutside'); grid(ax4, 'on');
xlabel(ax4, 't [s]'); xlim(ax4, t0 + [0 len]);
linkaxes([ax1 ax2 ax3 ax4], 'x');
exportgraphics(fig, sprintf('soundfield_map_%s_%gs.png', nm, t0), 'Resolution', 110);
fprintf('%g s: %d tracks >= 30 ms\n', t0, numel(Tr));
for k = 1:numel(Tr)
    h = Tr(k).hit;
    fprintf('   id %3d  %.3f-%.3f s  az %+4.0f -> %+4.0f  (min %+4.0f, max %+4.0f)  peak vel %5.0f deg/s  lev %5.1f dB\n', ...
        Tr(k).id, t0 + Tr(k).t(1), t0 + Tr(k).t(end), Tr(k).az(1), Tr(k).az(end), min(Tr(k).az(h)), max(Tr(k).az(h)), ...
        max(abs(Tr(k).vel(h))), max(Tr(k).lev));
end
fprintf('%g s: short frames in %d of %d granules\n', t0, sum(As.gShort), numel(As.gShort));
close(fig);
end

function map(ax, A, az, t0, ttl, norm)
    M = conv2(A.Hist, ones(3,1)/3, 'same');         % 15 deg azimuth smoothing
    if norm, M = M ./ max(max(M, [], 1), eps); M = 10*log10(M + 1e-6); lim = [-30 0];
    else, M = 10*log10(M + 1e-12); lim = max(M(:)) + [-50 0]; end
    imagesc(ax, t0 + A.tc, az, M); axis(ax, 'xy'); clim(ax, lim); colormap(ax, turbo);
    yticks(ax, -180:45:180); ylabel(ax, 'az [deg] (+ = left)'); title(ax, ttl);
    hold(ax, 'on');
    for a = [30 -30 0 110 -110], yline(ax, a, ':', 'Color', [1 1 1 0.5]); end
    hold(ax, 'off');
end
