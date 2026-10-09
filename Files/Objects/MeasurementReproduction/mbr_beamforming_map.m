function R = mbr_beamforming_map(h, arr, fs, win, opt)
%MBR_BEAMFORMING_MAP  Steered-response power map of a time window of a
%spherical-array RIR: max-rE SH beamformer (port of
%offl_conventional_window_beamforming.m, M. Bubrik).
%
%   R = mbr_beamforming_map(h, arr, fs, [t0 t1] ms, opt)
%     h   [T x Nch] RIR, arr mbr_mic_array
%     opt.N (3), opt.Nf (4096), opt.fcLP (3000 Hz)
%   R.azim [1 x Na] (rad), R.zen [Nz x 1] (rad), R.map [Nz x Na] power,
%   normalised to max 1; R.polarAz, R.polar (horizontal plane);
%   R.peak = [azimuth zenith] (deg)
%
%   Per Nf-sample block: SHT (pinv), radial equalisation H = Hrad Hreg Hlp
%   (linear convolution, block output length kept, as the OLS of the
%   original), then sum over time of |Y(dir) diag(a_n) psi|^2.

if nargin < 5, opt = struct(); end
N    = getopt(opt, 'N', min(3, arr.Nmax));
Nf   = getopt(opt, 'Nf', 4096);
fcLP = getopt(opt, 'fcLP', 3000);

seg = mbr_window(h, fs, win);
F = mbr_encoder_filters(arr.R, N, fs, Nf, fcLP);
hf = ifft(F.H, [], 1, 'symmetric');                 % [Nf x N+1] (as the original)

nB = ceil(size(seg,1)/Nf);
x = zeros(nB*Nf, size(seg,2)); x(1:size(seg,1),:) = seg;
psi = arr.sht(x, N);                                 % [nB*Nf x nSH]
L = 2^nextpow2(size(psi,1) + Nf - 1);
Psi = fft(psi, L, 1);
Hf = fft(hf, L, 1);
psi = real(ifft(Psi.*Hf(:, F.ell + 1), [], 1));
psi = psi(1:nB*Nf, :);

a_n = mbr_max_re_weights(N);
Na = 90; Nz = 46;
R.azim = (0:Na-1)/Na*2*pi;
R.zen  = (0:Nz-1)'/(Nz-1)*pi;
[Az, Ze] = meshgrid(R.azim, R.zen);
B = getSpherHarmMx(Ze(:), Az(:), N, 'real').*a_n;
R.map = reshape(sum(abs(B*psi.').^2, 2), size(Az));
R.map = R.map/max(R.map(:));
R.polarAz = (0:359)/360*2*pi;
Bp = getSpherHarmMx(pi/2*ones(360,1), R.polarAz(:), N, 'real').*a_n;
R.polar = sum(abs(Bp*psi.').^2, 2); R.polar = R.polar/max(R.polar);
[~, i] = max(R.map(:));
R.peak = [Az(i) Ze(i)]*180/pi;
R.method = 'beamforming';
R.window = win;
end

function v = getopt(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
