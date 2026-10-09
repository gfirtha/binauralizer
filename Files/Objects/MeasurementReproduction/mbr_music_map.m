function R = mbr_music_map(h, arr, fs, win, opt)
%MBR_MUSIC_MAP  SH-MUSIC direction map of a time window of a
%spherical-array RIR (port of offline_MUSIC_windowed.m V9, M. Bubrik).
%
%   R = mbr_music_map(h, arr, fs, [t0 t1] ms, opt)
%   opt (defaults of the thesis): N 3, L 6 (signal subspace dimension),
%       fBand [800 1600] Hz, Nfft 4096, alphaTik 0.01 (Tikhonov of the
%       modal equaliser), whitening true, diagLoad 1e-4, gridDeg 2,
%       minWinMs 0.5 (shorter windows are widened around their centre)
%   R.phi [1 x Nphi], R.theta [1 x Nth] (deg), R.map [Nth x Nphi] dB
%   (max 0), R.peak = [azimuth polar] (deg), R.eig (descending),
%   R.eigRatio = eig1/eig2 [dB]
%
%   Steps: Tukey(0.2) window, FFT, SHT (complex SH, pinv), modal
%   equalisation b_n*/(|b_n|^2 + alpha) per bin, cross-spectral matrix
%   averaged over the band, whitening by the equalised sensor-noise
%   covariance, eigen-decomposition, MUSIC pseudo-spectrum of the
%   noise subspace on the equalised steering vectors.

if nargin < 5, opt = struct(); end
N        = getopt(opt, 'N', min(3, arr.Nmax));
L        = getopt(opt, 'L', 6);
fBand    = getopt(opt, 'fBand', [800 1600]);
Nfft     = getopt(opt, 'Nfft', 4096);
alphaTik = getopt(opt, 'alphaTik', 0.01);
whiten   = getopt(opt, 'whitening', true);
diagLoad = getopt(opt, 'diagLoad', 1e-4);
gridDeg  = getopt(opt, 'gridDeg', 2);
minWinMs = getopt(opt, 'minWinMs', 0.5);
c = 343.1;
M = (N+1)^2;
L = max(1, min(L, M - 2));

if diff(win) < minWinMs
    tc = mean(win); win = tc + [-1 1]*minWinMs/2;
end

% modal equaliser per bin
SHT_op = pinv(arr.sh_matrix(N, 'complex'));          % [M x Nch]
f_vec = (0:Nfft-1)'/Nfft*fs;
f_idx = find(f_vec >= fBand(1) & f_vec <= fBand(2) & f_vec > 0);
f_idx = f_idx(f_idx <= floor(Nfft/2) + 1);
nn = repelem((0:N)', 2*(0:N)' + 1);
Q = numel(f_idx);
Binv = zeros(M, Q); Davg = zeros(M, 1);
for q = 1:Q
    kr = 2*pi*f_vec(f_idx(q))/c*arr.R;
    bn = 4*pi*1i.^(nn + 1)./(kr^2.*dh2(kr, nn));
    Binv(:,q) = conj(bn)./(abs(bn).^2 + alphaTik + 1e-30);
    Davg = Davg + abs(bn).^2./(abs(bn).^2 + alphaTik + 1e-30);
end
Davg = Davg/Q;

% whitening operator (equalised sensor noise)
if whiten
    Sn = zeros(M);
    for q = 1:Q
        E = Binv(:,q).*SHT_op;
        Sn = Sn + E*E';
    end
    Sn = Sn/Q;
    [U, Dd] = eig(Sn + diagLoad*trace(Sn)*eye(M));
    Wh = U*diag(1./sqrt(real(diag(Dd))))*U';
else
    Wh = eye(M);
end

% cross-spectral matrix of the window
[seg] = mbr_window(h, fs, win);
P = fft(seg.*tukeywin(size(seg,1), 0.2), Nfft, 1);
S = zeros(M);
for q = 1:Q
    a = Binv(:,q).*(SHT_op*P(f_idx(q),:).');
    S = S + a*a';
end
S = Wh*(S/Q)*Wh;
S = (S + S')/2;
[V, Dd] = eig(S);
[lam, ix] = sort(real(diag(Dd)), 'descend');
V = V(:, ix);

% MUSIC pseudo-spectrum
theta = deg2rad(0:gridDeg:180);
phi   = deg2rad(0:gridDeg:360-gridDeg);
[TH, PH] = meshgrid(theta, phi);                    % [Nphi x Nth]
Yg = conj(getSpherHarmMx(TH(:), PH(:), N, 'complex')).';
Yg = Wh*(Davg.*Yg);
Un = V(:, L+1:end);
W = Un'*Yg;
Pm = 1./real(sum(conj(W).*W, 1));
Pm = reshape(Pm, size(TH));
R.map = (10*log10(Pm/max(Pm(:)))).';               % [Nth x Nphi]
[~, i] = max(Pm(:));
R.peak = real([PH(i) TH(i)])*180/pi;
R.phi = rad2deg(phi); R.theta = rad2deg(theta);
R.eig = lam;
R.eigRatio = real(10*log10(lam(1)/max(lam(2), eps)));
R.L = L; R.fBand = fBand; R.bins = Q;
R.method = 'music';
R.window = win;
end

function d = dh2(x, n)
% derivative of the spherical Hankel function of the 2nd kind (thesis form)
jn = @(x, n) sqrt(pi./(2*x)).*besselj(n + 0.5, x);
yn = @(x, n) sqrt(pi./(2*x)).*bessely(n + 0.5, x);
h2 = @(x, n) jn(x, n) - 1i*yn(x, n);
d = n./x.*h2(x, n) - h2(x, n + 1);
end

function v = getopt(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
