function I = separator_info(method)
%SEPARATOR_INFO  Description of a separation method (Setup.Separation.Method)
%for the "Separator info" button of the Source separation tab.
%
%   I = separator_info(method) -> struct:
%     Title   one-line name
%     Status  'reference' | 'working' | 'experimental' | 'legacy'
%     Params  UI parameters the method reads (names of SeparationPanel.Params)
%     Text    detailed description (char, newline separated)
%   The text ends with the class header of the implementing class.
%
%   The status / measured notes come from the synthetic benchmark of
%   2026-10-09 (static panned sources, a fly-by over a centre vocal, a 5.1
%   scene with a C vocal, an L-Ls phantom and a circling source).

m = lower(method);
switch m
    case 'passthrough'
        cls = 'passthrough_separator';
        I.Title = 'Passthrough (channels as objects)';
        I.Status = 'reference';
        I.Params = {};
        txt = {
            'Every input channel becomes one object at its standard loudspeaker'
            'direction: stereo +-30 deg, 5.0/5.1/7.1 ITU layout, LFE unpositioned.'
            'No analysis, no latency. The reference that every other method should'
            'beat: it is what the binauralizer plays without separation.'
            ''
            'Outputs: one object per input channel.'};
    case 'ltrt'
        cls = 'ltrt_decoder';
        I.Title = 'Lt/Rt matrix decoder (stereo -> L C R Ls Rs)';
        I.Status = 'working';
        I.Params = {};
        txt = {
            'Frequency-domain active 4-2-4 matrix decoder (in the spirit of the'
            'classic film surround decoders). Per critical band, from smoothed'
            'powers: left/right dominance lr and centre/surround dominance cs of'
            'mid M = (L+R)/sqrt2 and side D = (L-R)/sqrt2.'
            '  C   = max(cs,0) * M,      Sur = max(-cs,0) * D'
            '  L'' = L - (C+Sur)/sqrt2,   R'' = R - (C-Sur)/sqrt2'
            '  Ls, Rs = Sur steered by lr, IN PHASE (no anti-phase pair, no delay).'
            'Also used inside the wavenumber and cluster separators (same STFT).'
            ''
            'Outputs (stereo): L C R Ls Rs at 30 0 -30 110 -110 deg. Other input'
            'widths pass through unchanged. Latency N/2 (N = 1024).'
            ''
            'Benchmark: better than the broadband matrix decoder on every scene'
            '(centre/side sources land in C / L / R). Not power preserving for'
            'mixed content (about -1 dB).'};
    case 'stft'
        cls = 'stft_separator';
        I.Title = 'STFT fixed directions (slots)';
        I.Status = 'working';
        I.Params = {'K', 'Width', 'Base', 'AmbAz', 'AmbGain', 'Tau'};
        txt = {
            'Primary / ambient decomposition of stereo after Goodwin-Jot and'
            'He-Tan-Gan (PCA). Per STFT bin, from the 2x2 covariance smoothed over'
            'Tau: principal eigenvector u = panning direction of the dominant source,'
            'primary P = sqrt(g) u''x with g = (l1-l2)/l1, ambience = the power'
            'complementary rest. Direction from the stereo tangent law of |u|'
            '(content base +-Base), stretched to +-Width.'
            'Each bin''s primary part is panned (constant power) between the two'
            'neighbouring of K FIXED front slots; the ambience pair goes to +-AmbAz'
            '(scaled by AmbGain). Objects stand still, the direction per bin stays'
            'continuous.'
            '5.1/7.1: L/R are decomposed, C goes to the centre slot, surrounds to'
            'the ambience objects, LFE dropped.'
            ''
            'Outputs: K front slots + 2 ambience objects. Latency N/2 (N = 2048).'
            ''
            'Benchmark: power preserving (0.97-1.0), robust, but fixed slots: a'
            'source between two slots is split between them.'
            'Known issue: mono input is counted twice (centre +7.7 dB).'};
    case 'stft_objects'
        cls = 'stft_separator';
        I.Title = 'STFT tracked objects (moving sources)';
        I.Status = 'working';
        I.Params = {'K', 'Width', 'Base', 'AmbAz', 'AmbGain', 'Tau'};
        txt = {
            'Same PCA primary/ambient decomposition as "STFT fixed directions", but'
            'instead of fixed slots up to K MOVING objects follow the sources:'
            'the primary energy''s direction histogram (smoothed in time and'
            'azimuth) is peak-picked every frame, peaks are tracked (nearest track'
            'within a gate, smoothed azimuth, attack/release activity), and each'
            'bin''s primary part goes to the active tracks near its direction'
            '(Gaussian masks, power complementary). Object positions follow the'
            'tracks; the direction map shows the tracks.'
            'Informed demixing (checkbox, experimental): where a bin lies between'
            'two active sources, the 2x2 mixing of their panning vectors is'
            'inverted (regularised LS). Marked in the code as making results worse.'
            ''
            'Outputs: K objects + 2 ambience objects. Latency N/2 (N = 2048).'
            ''
            'Benchmark: the BEST stereo method for moving sources (fly-by captured'
            'in one object, r2 = 0.92, ~4 deg tracking error) and good on static'
            'panned sources. Recommended default for stereo.'
            'Known issues: when no track is active the primary energy is dropped;'
            'on 5.1 input the C channel is added to object 3, which can move.'};
    case 'avendano'
        cls = 'avendano_separator';
        I.Title = 'Avendano-Jot coherence / panning index (slots)';
        I.Status = 'working';
        I.Params = {'K', 'Width', 'Base', 'AmbAz', 'AmbGain', 'Tau'};
        txt = {
            'Stereo upmix after Avendano & Jot (ICASSP/AES 2002, JAES 2004,'
            'WASPAA 2003). Per STFT bin: inter-channel coherence (smoothed over'
            'Tau) through a tanh mask -> ambience share; instantaneous similarity'
            'and panning index -> amplitude-panning vector a. The direct part is'
            'downmixed onto a and panned between the two neighbouring of K fixed'
            'front slots; the downmix residual joins the ambience pair at +-AmbAz.'
            'HardPanGuard: a hard-panned source (coherence ~ 0) is not taken as'
            'ambience (mask scaled by the channel balance).'
            ''
            'Outputs: K front slots + 2 ambience objects. Latency N/2 (N = 2048).'
            ''
            'Benchmark: slightly behind "STFT fixed directions", which covers the'
            'same use. Known issues: an exactly silent channel makes sign(0) = 0 ->'
            'a hard-panned source is put in the centre; with inter-channel phase'
            'differences the output is up to +3 dB (the header''s "power kept"'
            'claim does not hold).'
            'Literature: Documentation/Separation/AvendanoJot.'};
    case 'soundfield'
        cls = 'soundfield_separator';
        I.Title = 'Sound field sectors + tracking (moving objects, 360 deg)';
        I.Status = 'experimental';
        I.Params = {'K', 'AmbAz', 'AmbGain'};
        txt = {
            'Virtual DirAC analysis of the field the channels would create at the'
            'listener (every channel = plane wave from its loudspeaker). Per band:'
            'max-rE sector beams (several simultaneous directions), coherent'
            'neighbour sectors merged into groups, per group intensity direction,'
            'energy and diffuseness. Direct power -> direction histogram -> peaks'
            '-> Kalman tracker ([azimuth; angular velocity], novelty births) ->'
            'K object slots (fade in/out, no jumps). Each group''s direct part goes'
            'to the objects nearby; diffuse and unassigned parts stay in the bed.'
            ''
            'Outputs: K objects + bed (stereo: 2 beds at +-AmbAz; 5.1: 5 beds at'
            'the loudspeaker directions, LFE dropped). Latency N/2 (N = 1024).'
            ''
            'Benchmark: power preserving, but on static stereo a panned guitar ended'
            'in the rear bed (untracked front content plays from +-110 deg) and'
            'the run was ~1x real time in MATLAB (the slowest method). The'
            'offline analyzer (main_soundfield_map) and this streaming path have'
            'separate detection code.'};
    case 'wavenumber'
        cls = 'wavenumber_separator';
        I.Title = 'Local wavenumber + curvature (objects with distance)';
        I.Status = 'experimental';
        I.Params = {'K', 'AmbGain'};
        txt = {
            'Ray-theory description of the local wavefront at the listener (after'
            'Firtha, "Limitations of WFS"): channels = point sources at R m, closed-'
            'form pressure, gradient and Hessian at the origin -> local wavenumber'
            'vector k^ and wavefront curvature per bin. Stereo is Lt/Rt-decoded'
            'first (L C R Ls Rs).'
            'Per critical band: one direction (power-weighted k^) and a direct share'
            '(directional consistency^2); transient-adaptive block switching for'
            'the tracker; Kalman tracker -> K object slots; band content (omni'
            'pressure) to the objects, the rest to the bed. The wavefront curvature'
            'can give each object a distance (UseDistance, off by default).'
            ''
            'Outputs: K objects + bed (the decoded / input channels, LFE dropped).'
            'Latency N/2.'
            ''
            'Benchmark: good tracking of a fly-by (~5 deg) and of a centre vocal;'
            'physics of the gradient / Hessian formulas re-derived and correct.'
            'Limit: ONE direction per band, so two simultaneous sources in a band'
            'are not separated. Objects sit at radius 0.67, beds at 1.'};
    case 'cluster'
        cls = 'cluster_separator';
        I.Title = 'Wavenumber clustering (TF bins -> sources)';
        I.Status = 'experimental';
        I.Params = {'K', 'AmbGain'};
        txt = {
            'Time-frequency bins are clustered on the direction of their local'
            'wavenumber vector (wavenumber analysis, stereo Lt/Rt-decoded first):'
            'online EM of a von Mises mixture + one uniform "diffuse" component'
            '(the bed), statistics with forgetting, so a source is what keeps'
            'explaining energy from one direction. Births from the diffuse'
            'residual histogram or onsets, deaths and merges. Per source a Kalman'
            'filter predicts its direction for the next E-step (fly-bys stay in'
            'their component). Soft masks = posteriors; each source is an object'
            '(fade in/out), the diffuse posterior stays in the bed.'
            ''
            'Outputs: K objects + bed (decoded / input channels, LFE dropped).'
            'Latency N/2.'
            ''
            'Benchmark: the BEST method on 5.1 (L-Ls phantom found at 71 deg for'
            'a true 70 deg; best of all on the circling source, still ~30 deg'
            'error). Power preserving per band. Weaker than stft_objects on'
            'stereo. A hybrid EM: the mean comes from the Kalman state.'};
    case 'measured_room'
        cls = 'mbr_room_separator';
        I.Title = 'Measured room: reflections + SDM diffuse part';
        I.Status = 'experimental';
        I.Params = {};
        txt = {
            'Parametric reproduction of a spherical-array RIR (Measurement-based'
            'reproduction tab, mbr_decompose_rir): every reliably detected'
            'reflection is a point source with its own local impulse response'
            '(max-rE beam towards it, time-focused; omni above the aliasing'
            'frequency); the residual - checked for diffuseness by DirAC'
            'diffuseness and Avendano-Jot coherence - is rendered by the Spatial'
            'Decomposition Method onto fixed directions. The media (mono) is'
            'convolved with every IR. Sent from that tab, not offered in the'
            'separation method list.'};
    % ---- legacy ----
    case 'matrix'
        cls = 'matrix_separator';
        I.Title = 'Matrix (broadband Lt/Rt) [legacy]';
        I.Status = 'legacy';
        I.Params = {};
        txt = {
            'Broadband passive Lt/Rt matrix with correlation steering, anti-phase'
            'delayed surround pair. Superseded by the Lt/Rt matrix decoder (per'
            'band, better on every benchmark scene); output power 0.5.'};
    case 'avendano_ring'
        cls = 'avendano_ring_separator';
        I.Title = 'Avendano-Jot 5.1 ring pair bank [legacy]';
        I.Status = 'legacy';
        I.Params = {'K', 'Width', 'Base', 'AmbAz', 'AmbGain', 'Tau'};
        txt = {
            'Avendano-Jot on the four non-centre 5.1 pairs. Drops the C channel'
            'entirely (dialogue lost: centre vocal r2 = 0.00 in the benchmark).'};
    case 'avendano_5ch'
        cls = 'avendano_5ch_separator';
        I.Title = 'Avendano-Jot 5.1 five-pair bank [legacy]';
        I.Status = 'legacy';
        I.Params = {'K', 'Width', 'Base', 'AmbAz', 'AmbGain', 'Tau'};
        txt = {
            'Avendano-Jot on five adjacent 5.1 pairs; each channel split by 1/sqrt2'
            'and re-summed coherently: direct part +3 dB (output power 1.6).'};
    case 'phantom_objects'
        cls = 'phantom_object_separator';
        I.Title = 'Phantom objects + music bed [legacy]';
        I.Status = 'legacy';
        I.Params = {'K', 'AmbGain', 'Tau'};
        txt = {
            'Broadband covariance + ICA + phantom-gain dictionary. Broken: output'
            'power x3.5 (overlap-add gain), synthesis uses the first dictionary'
            'template instead of the matched gains (5.1 C lost), mono bed.'};
    case 'phantom_traj'
        cls = 'phantom_trajectory_separator';
        I.Title = 'Phantom trajectory [legacy]';
        I.Status = 'legacy';
        I.Params = {'K', 'AmbGain', 'Tau'};
        txt = {
            'Subset of "Phantom objects" (one eigenvector, no ICA) with the same'
            'defects, plus slot signals subtracted twice. Output power x4.'};
    otherwise
        error('separator_info:method', 'Unknown separation method "%s".', method);
end
if isempty(I.Params), p = 'none'; else, p = strjoin(I.Params, ', '); end
I.Text = strjoin([{sprintf('%s   [%s]', I.Title, upper(I.Status)), ...
    sprintf('Method ''%s'', class %s', m, cls), ''}, txt(:)', ...
    {'', ['UI parameters used: ' p]}], newline);
hdr = class_header(cls);
if ~isempty(hdr)
    I.Text = [I.Text newline '---- Implementation notes (class header) ----' newline hdr];
end
end

function h = class_header(cls)
% leading comment block of the class file, '%' stripped
h = '';
f = which(cls);
if isempty(f), return; end
L = readlines(f);
k = find(startsWith(strtrim(L), '%'), 1);
if isempty(k), return; end
out = strings(0, 1);
while k <= numel(L) && startsWith(strtrim(L(k)), '%')
    out(end+1) = regexprep(L(k), '^\s*%', ''); %#ok<AGROW>
    k = k + 1;
end
h = char(strjoin(out, newline));
end
