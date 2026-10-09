function R = mbr_herglotz_echoes(h, arr, fs, win, opt)
%MBR_HERGLOTZ_ECHOES  Echo detection and localisation in a time window of
%a spherical-array RIR: Herglotz localisation function with Tikhonov
%regularisation, Gaussian (von Mises-like) mixture fit and model-order
%selection, after Massé et al. (2026). Port of
%offl_herglotz_echo_detection.m (M. Bubrik).
%
%   R = mbr_herglotz_echoes(h, arr, fs, [t0 t1] ms, opt)
%   opt (thesis defaults): N 3, Nf 4096, fcLP 2100 Hz, D 1024 (Fibonacci
%       quadrature nodes), muRel 1e-4 (Tikhonov mu^2 = muRel ||D D^H||),
%       fBand [] (auto: max-directivity band up to kr = N, edge diTol
%       0.1 dB), Pmax 10, lam0 30, kappa 0.75 (C = ln L + kappa P)
%   R.echoes struct array: az, zen, el (deg), relDb (re strongest),
%       lambda (sharpness); R.map [Nz x Na] (dB, max 0) on R.azim/R.zen
%       (rad); R.polarAz/R.polar (horizontal cut) and R.polarModel;
%       R.C, R.L model-selection trace; R.band [fLo fHi] Hz.

if nargin < 5, opt = struct(); end
N     = getopt(opt, 'N', min(3, arr.Nmax));
Nf    = getopt(opt, 'Nf', 4096);
fcLP  = getopt(opt, 'fcLP', 2100);
D     = getopt(opt, 'D', 1024);
muRel = getopt(opt, 'muRel', 1e-4);
fBand = getopt(opt, 'fBand', []);
diTol = getopt(opt, 'diTol', 0.1);
Pmax  = getopt(opt, 'Pmax', 10);
lam0  = getopt(opt, 'lam0', 30);
kappa = getopt(opt, 'kappa', 0.75);
c0 = 343.1;
nSH = (N+1)^2;

F = mbr_encoder_filters(arr.R, N, fs, Nf, fcLP);
ell = F.ell;

% usable band: maximum directivity index (Massé II-B / III-B)
Nh = Nf/2 + 1; fh = F.freq(1:Nh);
wl = 2*(0:N)' + 1;
DI = abs(F.c(1:Nh,:)*wl).^2./(abs(F.c(1:Nh,:)).^2*wl);
if isempty(fBand)
    fAlias = N*c0/(2*pi*arr.R);
    [~, ia] = min(abs(fh - fAlias));
    bad = find(10*log10(DI(1:ia)) < 10*log10(DI(ia)) - diTol, 1, 'last');
    if isempty(bad), ilo = 2; else, ilo = max(bad + 1, 2); end
    kb = (ilo:ia)';
else
    kb = find(fh >= fBand(1) & fh <= fBand(2));
end
R.band = [fh(kb(1)) fh(kb(end))];

% quadrature grid and Tikhonov inverses
[thF, phF, xyzF] = fibonacci_sphere(D);
Yfit = getSpherHarmMx(thF, phF, N, 'real');
beta = 4*pi/D;
Kb = numel(kb);
c_ell = F.c(kb, ell + 1);
Ainv = zeros(nSH, nSH, Kb);
for k = 1:Kb
    Dk = beta*(c_ell(k,:).').*Yfit.';
    Gk = Dk*Dk';
    Ainv(:,:,k) = inv(Gk + muRel*norm(Gk)*eye(nSH));
end

% window -> encoded spectra -> regularised plane-wave amplitudes
seg = mbr_window(h, fs, win);
if size(seg,1) > Nf, seg = seg(1:Nf,:); end
Psi = fft(arr.sht(seg, N), Nf, 1);
U = (F.H(kb, ell + 1).*Psi(kb,:)).';
W = zeros(nSH, Kb);
for k = 1:Kb, W(:,k) = Ainv(:,:,k)*U(:,k); end
B = beta*conj(c_ell).'.*W;

% frequency-averaged localisation function (Eq. 11)
Mfit = abs(Yfit*B);
gam = max(Mfit, [], 1).';
Gtot = sum(gam);
mbar = (Mfit*gam)/Gtot;
scale = max(mbar);
mbar = mbar/scale;

Na = 90; Nz = 46;
R.azim = (0:Na-1)/Na*2*pi;
R.zen = (0:Nz-1)'/(Nz-1)*pi;
[Az, Ze] = meshgrid(R.azim, R.zen);
Mmap = (abs(getSpherHarmMx(Ze(:), Az(:), N, 'real')*B)*gam)/Gtot/scale;
R.map = 20*log10(reshape(Mmap, size(Az)) + eps);
R.polarAz = (0:359)/360*2*pi;
Ypol = getSpherHarmMx(pi/2*ones(360,1), R.polarAz(:), N, 'real');
R.polar = (abs(Ypol*B)*gam)/Gtot/scale;
xyz_pol = ang2xyz(pi/2*ones(360,1), R.polarAz(:));

% Gaussian mixture fit + echo count (Eqs. 12-15)
[X, R.C, R.L] = estimate_echoes(xyzF, mbar, Pmax, lam0, kappa);
R.polarModel = rbf_model(X, xyz_pol);
P = size(X, 2);
az = mod(X(2,:), 2*pi)*180/pi;
zen = mod(X(1,:), 2*pi)*180/pi; zen(zen > 180) = 360 - zen(zen > 180);
hrel = X(3,:)/max(X(3,:));
R.echoes = struct('az', num2cell(az), 'zen', num2cell(zen), 'el', num2cell(90 - zen), ...
    'relDb', num2cell(20*log10(hrel)), 'lambda', num2cell(exp(X(4,:))));
R.P = P;
R.method = 'herglotz';
R.window = win;
end

% ---------------- model fit (unchanged from the thesis code) ----------------
function [X, Cs, Ls] = estimate_echoes(xyz, mbar, Pmax, lam0, kappa)
w = mbar; Mtot = sum(mbar);
X = zeros(4, 0);
Cprev = inf; Cs = []; Ls = [];
for P = 1:Pmax
    if isempty(X), mhat = zeros(size(mbar)); else, mhat = rbf_model(X, xyz); end
    [rm, j] = max(mbar - mhat);
    new = [acos(xyz(j,3)); atan2(xyz(j,2), xyz(j,1)); max(rm, 1e-3); log(lam0)];
    [Xn, Ln] = bb_fit([X, new], xyz, mbar, w, Mtot);
    C = log(max(Ln, eps)) + kappa*P;
    Cs(P) = C; Ls(P) = Ln; %#ok<AGROW>
    if C > Cprev, break; end
    X = Xn; Cprev = C;
end
end

function [X, Lb] = bb_fit(X0, xyz, mbar, w, Mtot)
% gradient descent with Barzilai-Borwein steps on Eq. (14)
P = size(X0, 2);
x = X0(:);
lo = repmat([-inf; -inf; 1e-6; log(0.5)], P, 1);
hi = repmat([ inf;  inf;  inf; log(300)], P, 1);
fun = @(x) rbf_loss(x, xyz, mbar, w, Mtot, P);
[L, g] = fun(x);
xb = x; Lb = L;
stepCap = 0.2;
alpha = 0.05/(max(abs(g)) + eps);
for it = 1:400
    xn = min(max(x - alpha*g, lo), hi);
    [Ln, gn] = fun(xn);
    s = xn - x; y = gn - g; sy = s.'*y;
    if sy > eps, alpha = (s.'*s)/sy; else, alpha = 0.05/(max(abs(gn)) + eps); end
    alpha = min(alpha, stepCap/(max(abs(gn)) + eps));
    x = xn; g = gn; L = Ln;
    if L < Lb, Lb = L; xb = x; end
    if norm(g, inf) < 1e-9, break; end
end
X = reshape(xb, 4, P);
end

function [L, g, mhat] = rbf_loss(x, xyz, mbar, w, Mtot, P)
% mhat_d = sum_p |h_p| exp(lambda_p (d_p.r_d - 1)), lambda = exp(x4)
X = reshape(x, 4, P);
th = X(1,:); ph = X(2,:); h = X(3,:); lam = exp(X(4,:));
dp = [sin(th).*cos(ph); sin(th).*sin(ph); cos(th)];
C = xyz*dp;
E = exp(lam.*(C - 1));
T = abs(h).*E;
mhat = sum(T, 2);
r = mbar - mhat;
L = sum(w.*r.^2)/Mtot;
if nargout > 1
    dLdm = -2*(w.*r)/Mtot;
    dT_dC = T.*lam;
    ddth = [cos(th).*cos(ph); cos(th).*sin(ph); -sin(th)];
    ddph = [-sin(th).*sin(ph); sin(th).*cos(ph); zeros(1,P)];
    gth = sum(dLdm.*dT_dC.*(xyz*ddth), 1);
    gph = sum(dLdm.*dT_dC.*(xyz*ddph), 1);
    gh  = sum(dLdm.*sign(h).*E, 1);
    gl  = sum(dLdm.*T.*(C - 1).*lam, 1);
    g = [gth; gph; gh; gl];
    g = g(:);
end
end

function mhat = rbf_model(X, xyz)
P = size(X, 2);
if P == 0, mhat = zeros(size(xyz,1),1); return; end
[~, ~, mhat] = rbf_loss(X(:), xyz, zeros(size(xyz,1),1), zeros(size(xyz,1),1), 1, P);
end

function xyz = ang2xyz(th, ph)
xyz = [sin(th).*cos(ph), sin(th).*sin(ph), cos(th)];
end

function [th, ph, xyz] = fibonacci_sphere(D)
i = (0:D-1)' + 0.5;
z = 1 - 2*i/D;
th = acos(z);
ph = mod(pi*(1 + sqrt(5))*i, 2*pi);
xyz = ang2xyz(th, ph);
end

function v = getopt(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
