function D = mbr_decompose_rir(h, arr, fs, opt)
%MBR_DECOMPOSE_RIR  Spherical-array RIR -> reliable early reflections as
%point sources with their own local impulse responses + a residual that
%is checked for diffuseness and rendered by the Spatial Decomposition
%Method (SDM). Input of the parametric reproduction (mbr_room_separator).
%
%   D = mbr_decompose_rir(h, arr, fs, opt)
%
%   1. Encoding: SHT + rigid-sphere radial equaliser + crossover
%      regularisation (mbr_encoder_filters, no low-pass), the Nf/2 filter
%      delay removed: psi [T x (N+1)^2]. A unit plane wave from d gives
%      psi_nm = Hreg_n * Y_nm(d) s(t). Omni pressure p0: psi_00/Y_00 below
%      fx, above it the microphone mean equalised for a diffuse field (the
%      SH order 0 is spatially aliased there).
%   2. Detection: half-overlapping windows (winMs) from the direct sound
%      over rangeMs, echoes by mbr_herglotz_echoes; an echo is reliable if
%      its level re the window's strongest >= confDb, its sharpness lambda
%      is in [lamMin lamMax] (a plane wave seen at order N = 3 fits lambda
%      ~5.5; ~0.5 = no peak, >> 6 = an over-fitted side lobe), the window
%      energy >= floorDb re the strongest window, and its arrival (step 3)
%      is not within edgeMs of the window edges (a tail leaking in from
%      the neighbour window peaks at the edge; a real echo is inside one of
%      the overlapping windows). The direct sound (strongest echo of the
%      first window) is always kept; echoes closer than mergeMs (1.5) and
%      mergeDeg (45, ~ the order-3 beam width) are merged.
%   3. Local IR of a reflection ("beamforming and focus"): below fx the
%      max-rE beam towards d over the whole RIR, normalised per frequency so
%      that a unit plane wave from d gives a unit output (G(f) = sum_n a_n
%      (2n+1)/(4 pi) Hreg_n(f)); above fx (spatial aliasing, default
%      N c/(2 pi R)) the least-squares deconvolution of the 3 microphones
%      nearest to d with their rigid-sphere responses. Complementary
%      raised-cosine crossover (1/3 octave). Arrival = peak of the result
%      inside the window; time focus per band (Tukey [t - preMs, t +
%      postMs]: < fFocus 1/5 ms, fFocus..fx 0.5/2 ms, > fx 0.3/1 ms - the
%      low band cannot be localised in 1 ms).
%   4. Residual at the microphones: every reflection re-synthesised with the
%      rigid-sphere model (order 20) and subtracted, then encoded like the
%      RIR: psi_res, p_res.
%   5. Diffuseness check (fFocus..fx) of the original and of the residual: DirAC
%      diffuseness 1 - |<I>|/<E> (first-order SH) and the Avendano-Jot
%      short-time coherence of opposite microphones (low = diffuse), per
%      5 ms frame.
%   6. SDM of the residual: DOA per sample = smoothed pseudo-intensity
%      p v (below fx); p_res is split onto nDiff fixed directions with
%      energy-preserving weights w_m ~ exp(kappaSdm (u.d_m - 1)).
%
%   opt: rangeMs 50, winMs 2.5, confDb -6, lamMin 3, lamMax 9, edgeMs 0.25,
%        floorDb -30, maxRefl 24, mergeMs 1.5, mergeDeg 45,
%        preMs [1 0.5 0.3], postMs [5 2 1] (per band), fFocus 500, fx [] (auto),
%        nDiff 16, kappaSdm 6,
%        sdmSmoothMs 1, irLenMs 400 (length of the object IRs after the
%        direct sound), herglotz struct (passed on)
%   D: refl (t [s re direct], tAbs, az, el, relDb, lambda, gain), objIR
%      [L x K], diffIR [L x M], diffDirs [M x 2] (az el, deg), fs, tDirect,
%      p0, pres (from the direct sound on), diag (t, psiOrig, psiRes,
%      cohOrig, cohRes), energy (refl, res, total; from the direct sound on)

if nargin < 4, opt = struct(); end
rangeMs = getopt(opt, 'rangeMs', 50);
winMs   = getopt(opt, 'winMs', 2.5);
confDb  = getopt(opt, 'confDb', -6);
lamMin  = getopt(opt, 'lamMin', 3);
lamMax  = getopt(opt, 'lamMax', 9);
edgeMs  = getopt(opt, 'edgeMs', 0.25);
mergeMs = getopt(opt, 'mergeMs', 1.5);
mergeDeg = getopt(opt, 'mergeDeg', 45);       % ~ the order-3 beam width: side-lobe fits of one wave
floorDb = getopt(opt, 'floorDb', -30);
maxRefl = getopt(opt, 'maxRefl', 24);
preMs   = getopt(opt, 'preMs', [1.0 0.5 0.3]);     % per band: < fBand(1), fBand(1)..fx, > fx
postMs  = getopt(opt, 'postMs', [5.0 2.0 1.0]);
fBand1  = getopt(opt, 'fFocus', 500);              % lower band edge of the focus [Hz]
nDiff   = getopt(opt, 'nDiff', 16);
kSdm    = getopt(opt, 'kappaSdm', 6);
smMs    = getopt(opt, 'sdmSmoothMs', 1);
irLenMs = getopt(opt, 'irLenMs', 400);
hopt    = getopt(opt, 'herglotz', struct());
N = min(3, arr.Nmax);
c0 = 343.1;
fx = getopt(opt, 'fx', N*c0/(2*pi*arr.R));

T = size(h, 1);
nSH = (N+1)^2;

%% 1. encoding of the whole RIR + rigid-sphere model
[psi, Lf] = encode_sh(h, arr, N, fs);
Y00 = 1/sqrt(4*pi);
fL = (0:Lf-1)'/Lf*fs; fL = min(fL, fs - fL);       % |f| on the long grid
LP = crossover_lp(fL, fx);
F = mbr_encoder_filters(arr.R, N, fs, 4096, Inf);
U = arr.x0./sqrt(sum(arr.x0.^2, 2));
Lsim = 40;                                         % rigid-sphere series: valid to kR ~ 35 (19 kHz at R = 10 cm)
f1L = (0:Lf/2)'/Lf*fs;
BwL = sphere_bw(f1L, arr.R, Lsim, c0).*crossover_lp(f1L, 19000);   % rigid-sphere modal weights (one-sided), faded out where the series stops converging
Gam = diffuse_mean_gain(arr.R, U, Lsim, c0, fs);   % power of the mic mean in a diffuse field
GamL = interp1(Gam.f, Gam.g, fL, 'linear', 'extrap');
HPomni = (1 - LP)./sqrt(max(GamL, 1e-6));
% omni pressure: SH order 0 below fx, diffuse-equalised mic mean above
p0 = lin_filter(psi(:,1)/Y00, LP) + lin_filter(mean(h, 2), HPomni);

% first-order "velocity" axes of the real SH (sign of the ACN channels)
Yax = getSpherHarmMx([pi/2; pi/2; 0], [0; pi/2; 0], 1, 'real');   % +x, +y, +z
vIdx = [4 2 3]; vSgn = sign([Yax(1,4) Yax(2,2) Yax(3,3)]);
c1 = sqrt(3/(4*pi));
a_n = mbr_max_re_weights(N);
Gbeam = zeros(Lf, 1);                              % beam gain for a unit plane wave
for n = 0:N
    Hn = interp1((0:2048)'/4096*fs, real(F.Hreg(1:2049, n+1)), fL, 'linear', 'extrap');
    Gbeam = Gbeam + a_n(n^2+1)*(2*n+1)/(4*pi)*Hn;
end
Gbeam = max(Gbeam, 1e-3*max(Gbeam));

%% direct sound
env = sqrt(mean(h.^2, 2));
iD = find(env >= max(env)*10^(-20/20), 1);
t0ms = max(0, (iD - 1)/fs*1000 - 0.5);

%% 2. detection
a_n = mbr_max_re_weights(N);
cand = struct('tw', {}, 'az', {}, 'zen', {}, 'el', {}, 'relDb', {}, 'lambda', {}, 'hr', {}, 'win', {}, 'E', {});
hop = winMs/2;
nW = max(1, floor((rangeMs - winMs)/hop) + 1);
for w = 1:nW
    win = t0ms + (w-1)*hop + [0 winMs];
    if win(2)/1000*fs > T, break; end
    try
        R = mbr_herglotz_echoes(h, arr, fs, win, hopt);
    catch
        continue;
    end
    seg = mbr_window(p0, fs, win);
    hr = 10.^([R.echoes.relDb]/20);
    for k = 1:R.P
        e = R.echoes(k);
        cand(end+1) = struct('tw', w, 'az', wrap180(e.az), 'zen', e.zen, 'el', e.el, 'relDb', e.relDb, ...
            'lambda', e.lambda, 'hr', hr(k)/sum(hr), 'win', win, 'E', sum(seg.^2)); %#ok<AGROW>
    end
end
if isempty(cand), error('mbr_decompose_rir:none', 'No echoes detected.'); end
Emax = max([cand.E]);
first = find([cand.tw] == cand(1).tw); [~, j] = max([cand(first).relDb]); dIdx = first(j);
ok = [cand.relDb] >= confDb & [cand.lambda] >= lamMin & [cand.lambda] <= lamMax & [cand.E] >= Emax*10^(floorDb/10);
ok(dIdx) = true;
keep = find(ok);

%% 3. local IRs: normalised beam below fx, LS deconvolution of the nearest
%     microphones (rigid-sphere model) above; arrival + time focus
refl = struct('t', {}, 'tAbs', {}, 'az', {}, 'el', {}, 'relDb', {}, 'lambda', {}, 'gain', {});
S = {}; Dk = zeros(0, 3);
pre = round(preMs/1000*fs); post = round(postMs/1000*fs); edge = round(edgeMs/1000*fs);
for q = keep
    c = cand(q);
    d = [cosd(c.el)*cosd(c.az), cosd(c.el)*sind(c.az), sind(c.el)];
    Yd = getSpherHarmMx(acos(d(3)), atan2(d(2), d(1)), N, 'real');
    beam = lin_filter(psi*(a_n.*Yd).', LP./Gbeam);  % full length: the Hreg kernels are long
    i0 = max(1, round((c.win(1) - 2)/1000*fs) + 1); i1 = min(T, round((c.win(2) + 2 + postMs(1))/1000*fs) + 1);
    nS = i1 - i0 + 1; Ls = 2^nextpow2(nS + 64);
    f1 = (0:Ls/2)'/Ls*fs;
    Hm = sphere_bw(f1, arr.R, Lsim, c0)*legendre_all(U*d', Lsim);   % [nb x Q] mic responses
    [~, near] = maxk(U*d', 3);
    Pm = fft(h(i0:i1, near), Ls); Pm = Pm(1:Ls/2+1, :);
    H3 = Hm(:, near);
    S1 = sum(conj(H3).*Pm, 2)./(sum(abs(H3).^2, 2) + 1e-6);
    fS = min((0:Ls-1)'/Ls*fs, fs - (0:Ls-1)'/Ls*fs);
    hpS = 1 - crossover_lp(fS, fx);
    s_hf = real(ifft([S1; conj(S1(end-1:-1:2))].*hpS));
    s = beam(i0:i1) + s_hf(1:nS);
    j0 = max(1, round(c.win(1)/1000*fs) + 1 - i0 + 1); j1 = min(nS, round(c.win(2)/1000*fs) + 1 - i0 + 1);
    [~, jp] = max(abs(s(j0:j1))); jp = jp + j0 - 1;
    if q ~= dIdx && (jp - j0 < edge || j1 - jp < edge), continue; end   % leaked tail of a neighbour
    % frequency-dependent time focus: longer windows for the low band
    % (its content cannot be localised in 1 ms), complementary 3-band split
    i0w = max(1, i0 - post(1)); i1w = min(T, i1 + post(1));
    Lw = 2^nextpow2(i1w - i0w + 1 + 64);
    fW = min((0:Lw-1)'/Lw*fs, fs - (0:Lw-1)'/Lw*fs);
    Bl = crossover_lp(fW, fBand1); Bh = 1 - crossover_lp(fW, fx); Bm = 1 - Bl - Bh;
    sx = zeros(i1w - i0w + 1, 1); sx(i0 - i0w + (1:nS)) = s;
    jpw = jp + i0 - i0w;
    sk = zeros(T, 1);
    Bands = {Bl, Bm, Bh};
    for bnd = 1:3
        sb = real(ifft(fft(sx, Lw).*Bands{bnd})); sb = sb(1:numel(sx));
        w = zeros(numel(sx), 1);
        a = max(1, jpw - pre(bnd)); b = min(numel(sx), jpw + post(bnd));
        w(a:b) = tukeywin(b - a + 1, 0.5);
        sk(i0w:i1w) = sk(i0w:i1w) + sb.*w;
    end
    tAbs = (i0 - 1 + jp - 1)/fs;
    dup = false;
    for m = 1:numel(refl)
        if abs(refl(m).tAbs - tAbs) < mergeMs/1000 && acosd(min(1, Dk(m,:)*d')) < mergeDeg, dup = true; break; end
    end
    if dup, continue; end
    refl(end+1) = struct('t', 0, 'tAbs', tAbs, 'az', c.az, 'el', c.el, 'relDb', c.relDb, ...
        'lambda', c.lambda, 'gain', sqrt(sum(sk.^2))); %#ok<AGROW>
    S{end+1} = sk; Dk(end+1,:) = d; %#ok<AGROW>
end
% direct first, then the strongest up to maxRefl, in time order
g = [refl.gain];
rest = 2:numel(refl);
[~, o] = sort(g(rest), 'descend'); rest = rest(o(1:min(numel(o), maxRefl - 1)));
[~, o] = sort([refl(rest).tAbs]); sel = [1 rest(o)];
refl = refl(sel); S = S(sel); Dk = Dk(sel,:);
tD = refl(1).tAbs;
g0 = refl(1).gain;
for k = 1:numel(refl), refl(k).t = refl(k).tAbs - tD; refl(k).gain = refl(k).gain/g0; end

%% 4. residual at the microphones: the reflections re-synthesised with the
%     rigid-sphere model and subtracted, then encoded like the RIR
hres = h;
for k = 1:numel(refl)
    Sk = fft(S{k}, Lf); Sk = Sk(1:Lf/2+1);
    Mk = Sk.*(BwL*legendre_all(U*Dk(k,:)', Lsim));
    x = real(ifft([Mk; conj(Mk(end-1:-1:2,:))]));
    hres = hres - x(1:T, :);
end
psiRes = encode_sh(hres, arr, N, fs);
psiRes1 = psiRes(:,1:4);
pres = lin_filter(psiRes(:,1)/Y00, LP) + lin_filter(mean(hres, 2), HPomni);
Ssum = zeros(T, 1);
for k = 1:numel(refl), Ssum = Ssum + S{k}; end

%% 5. diffuseness diagnostics (below fx)
lpT = @(x) lowpass_fft(x, fs, fx);
iS = round(tD*fs) + 1;
iO = max(1, iS - round(1.5e-3*fs));                % object IR start: 1.5 ms before the direct peak
[p_o, v_o] = pv(lpT(psi(:,1:4)), vIdx, vSgn, Y00, c1);
[p_r, v_r] = pv(lpT(psiRes1), vIdx, vSgn, Y00, c1);
bpT = @(x) lowpass_fft(x, fs, fx) - lowpass_fft(x, fs, fBand1);
[p_o, v_o] = pv(bpT(psi(:,1:4)), vIdx, vSgn, Y00, c1);
[p_rd, v_rd] = pv(bpT(psiRes1), vIdx, vSgn, Y00, c1);
D.diag = diffuseness(p_o, v_o, p_rd, v_rd, h, hres, U, fs, iS, fx, fBand1);
D.diag.band = [fBand1 fx];

%% 6. SDM of the residual (objects start preAlignMs before the direct sound)
u = smoothdata(p_r.*v_r, 1, 'movmean', max(1, round(smMs/1000*fs)));
u = u./max(sqrt(sum(u.^2, 2)), eps);
[dirs, xyzD] = fib_dirs(nDiff);
W = exp(kSdm*(u*xyzD' - 1));
W = W./max(sqrt(sum(W.^2, 2)), eps);               % energy preserving
L = min(T - iO + 1, round(irLenMs/1000*fs));
D.diffIR = pres(iO:iO+L-1).*W(iO:iO+L-1, :);
D.diffDirs = dirs;

% object IRs aligned to the direct sound (t = 0 at sample 1)
D.objIR = zeros(L, numel(refl));
for k = 1:numel(refl)
    x = S{k}(iO:end);
    D.objIR(:,k) = x(1:L);
end
D.refl = refl;
D.preDelay = (iS - iO)/fs;                           % IR time of the direct sound [s]
D.hres = hres;
D.fs = fs; D.tDirect = tD; D.fx = fx;
D.p0 = p0; D.pres = pres;
D.energy = struct('refl', sum(Ssum(iS:end).^2), 'res', sum(pres(iS:end).^2), 'total', sum(p0(iS:end).^2));
D.nCandidates = numel(cand); D.nReliable = numel(keep); D.candidates = cand;
end

% ---------------------------------------------------------------------
function [p, v] = pv(psi1, vIdx, vSgn, Y00, c1)
% omni pressure and particle-velocity-like vector (|v| = |p| for a plane
% wave, pointing to the source) from first-order real SH
p = psi1(:,1)/Y00;
v = psi1(:, vIdx).*vSgn/c1;
end

function Dg = diffuseness(p_o, v_o, p_r, v_r, h, hres, U, fs, iS, fx, fLo)
fr = round(0.005*fs); nF = floor((numel(p_o) - iS + 1)/fr);
Dg.t = ((0:nF-1)' + 0.5)*fr/fs*1000;                % ms after the direct sound
Dg.psiOrig = dirac_psi(p_o, v_o, iS, fr, nF);
Dg.psiRes  = dirac_psi(p_r, v_r, iS, fr, nF);
% opposite microphone pairs (most antipodal partner of each microphone)
[~, j] = min(U*U', [], 2);
pairs = unique(sort([(1:size(U,1))' j], 2), 'rows');
Dg.cohOrig = pair_coherence(h, pairs, fs, iS, fr, nF, fLo, fx);
Dg.cohRes  = pair_coherence(hres, pairs, fs, iS, fr, nF, fLo, fx);
end

function psiD = dirac_psi(p, v, iS, fr, nF)
% DirAC diffuseness over 3 frames (15 ms)
I = p.*v; E = 0.5*(p.^2 + sum(v.^2, 2));
psiD = zeros(nF, 1);
for f = 1:nF
    r = iS + (max(0, f-2)*fr : min(nF, f+1)*fr - 1);
    psiD(f) = 1 - norm(sum(I(r,:), 1))/max(sum(E(r)), eps);
end
end

function coh = pair_coherence(x, pairs, fs, iS, fr, nF, fLo, fHi)
% Avendano-Jot inter-channel short-time coherence |<X Y*>| /
% sqrt(<|X|^2><|Y|^2>) of opposite microphones (STFT 256, hop fr,
% recursive smoothing tau 10 ms), averaged over fLo..fHi and the pairs:
% ~1 for a single wave, low for a diffuse field
Nw = 256; win = hann(Nw, 'periodic');
k = (0:Nw/2)'/Nw*fs; kb = k >= fLo & k <= fHi;
a = exp(-fr/fs/0.01);
coh = zeros(nF, 1);
Sab = 0; Saa = 0; Sbb = 0;
for f = 1:nF
    c0 = iS + (f - 1)*fr + fr/2;
    idx = round(c0 - Nw/2) + (0:Nw-1);
    if idx(1) < 1 || idx(end) > size(x, 1), continue; end
    X = fft(x(idx,:).*win); X = X(1:Nw/2+1,:);
    XA = X(:, pairs(:,1)); XB = X(:, pairs(:,2));
    Sab = a*Sab + (1-a)*XA.*conj(XB); Saa = a*Saa + (1-a)*abs(XA).^2; Sbb = a*Sbb + (1-a)*abs(XB).^2;
    c = abs(Sab)./sqrt(Saa.*Sbb + eps);
    coh(f) = mean(mean(c(kb,:)));
end
end

function [psi, Lf] = encode_sh(h, arr, N, fs)
% SHT + radial equaliser + crossover regularisation, filter delay removed
T = size(h, 1); Nf = 4096;
F = mbr_encoder_filters(arr.R, N, fs, Nf, Inf);
hrad = ifft(F.H, [], 1, 'symmetric');
Lf = 2^nextpow2(T + Nf);
Hl = fft(hrad, Lf, 1);
psi = real(ifft(fft(arr.sht(h, N), Lf, 1).*Hl(:, F.ell + 1), [], 1));
psi = psi(Nf/2 + (1:T), :);
end

function y = lin_filter(x, H)
% zero-phase filtering with the full-spectrum response H (length >= rows of x)
Lf = numel(H);
y = real(ifft(fft(x, Lf, 1).*H, [], 1));
y = y(1:size(x,1), :);
end

function Bw = sphere_bw(f, R, L, c0)
% (2n+1)/(4 pi) B_n(kR), n = 0..L: pressure on a rigid sphere for a unit
% plane wave = Bw * P_n(cos gamma) (mbr_simulate_sma_rir convention)
kr = 2*pi*f(:)/c0*R;
Bw = zeros(numel(kr), L+1);
ok = kr > 1e-6;
x = kr(ok);
jn = @(n) sqrt(pi./(2*x)).*besselj(n + 0.5, x);
yn = @(n) sqrt(pi./(2*x)).*bessely(n + 0.5, x);
h2 = @(n) jn(n) - 1i*yn(n);
for n = 0:L
    dh = n./x.*h2(n) - h2(n + 1);
    Bw(ok, n+1) = (2*n + 1)*1i^n./(1i*x.^2.*dh);
end
Bw(~ok, 1) = 1;
end

function Pn = legendre_all(x, L)
% [L+1 x numel(x)] Legendre polynomials P_0..P_L
x = x(:).';
Pn = zeros(L+1, numel(x)); Pn(1,:) = 1;
if L >= 1, Pn(2,:) = x; end
for n = 2:L
    Pn(n+1,:) = ((2*n - 1)*x.*Pn(n,:) - (n - 1)*Pn(n-1,:))/n;
end
end

function G = diffuse_mean_gain(R, U, L, c0, fs)
% power of the microphone mean for a unit diffuse field (average over
% directions of |mean_q H_q(f; d)|^2), on a coarse grid
G.f = linspace(0, fs/2, 1025)';
Bw = sphere_bw(G.f, R, L, c0);
[~, X] = fib_dirs(400);
g = zeros(numel(G.f), 1);
for k = 1:size(X, 1)
    g = g + abs(mean(Bw*legendre_all(U*X(k,:)', L), 2)).^2;
end
G.g = g/size(X, 1);
end

function L = crossover_lp(f, fx)
% zero-phase raised-cosine low-pass, 1/3 octave wide around fx (|f| input)
lo = fx*2^(-1/6); hi = fx*2^(1/6);
L = ones(size(f));
L(f >= hi) = 0;
m = f > lo & f < hi;
L(m) = 0.5 + 0.5*cos(pi*log2(f(m)/lo)/log2(hi/lo));
end

function y = lowpass_fft(x, fs, fx)
Lx = 2^nextpow2(size(x,1));
f = (0:Lx-1)'/Lx*fs; f = min(f, fs - f);
y = real(ifft(fft(x, Lx, 1).*crossover_lp(f, fx), [], 1));
y = y(1:size(x,1), :);
end

function [dirs, xyz] = fib_dirs(M)
i = (0:M-1)' + 0.5;
z = 1 - 2*i/M; th = acos(z); ph = mod(pi*(1 + sqrt(5))*i, 2*pi);
xyz = [sin(th).*cos(ph), sin(th).*sin(ph), cos(th)];
dirs = [wrap180(ph*180/pi), 90 - th*180/pi];
end

function a = wrap180(a), a = mod(a + 180, 360) - 180; end

function v = getopt(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
