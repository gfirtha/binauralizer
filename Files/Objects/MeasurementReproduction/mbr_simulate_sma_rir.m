function [h, truth] = mbr_simulate_sma_rir(arr, fs, echoes, opt)
%MBR_SIMULATE_SMA_RIR  Synthetic rigid-sphere array RIR of plane-wave
%impulse echoes (Massé et al. 2026, Sec. IV-A, Eq. 16; port of
%masse_chapter4_simulation_v4.m, M. Bubrik).
%
%   [h, truth] = mbr_simulate_sma_rir(arr, fs, echoes, opt)
%     echoes [P x 3]: azimuth (deg), elevation (deg), level (dB)
%     opt: Nsamples 256, Nf 4096, Lsim 20 (synthesis order), bpHz
%          [125 16000] (band, 1/12-octave cosine edges), delaySpanMs
%          [0.4 2.5] (arrivals spread uniformly), marginMs 1.0
%   h [Nsamples x Nch]; truth: az, el, relDb, delays (s)

if nargin < 4, opt = struct(); end
Nsamples = getopt(opt, 'Nsamples', 256);
Nf       = getopt(opt, 'Nf', 4096);
Lsim     = getopt(opt, 'Lsim', 20);
bpHz     = getopt(opt, 'bpHz', [125 16000]);
span     = getopt(opt, 'delaySpanMs', [0.4 2.5]);
margin   = getopt(opt, 'marginMs', 1.0);
c0 = 343.1;

Q = arr.nCh();
u = arr.x0./sqrt(sum(arr.x0.^2, 2));
P = size(echoes, 1);
amp = 10.^(echoes(:,3)/20);
delays = (linspace(span(1), span(2), P)' + margin)/1000;

nh = Nf/2 + 1;
fh = (0:nh-1)'/Nf*fs;
kr = 2*pi*fh/c0*arr.R; kr(1) = 1e-12;
jn = @(x, n) sqrt(pi./(2*x)).*besselj(n + 0.5, x);
yn = @(x, n) sqrt(pi./(2*x)).*bessely(n + 0.5, x);
h2 = @(x, n) jn(x, n) - 1i*yn(x, n);
dh2 = @(x, n) n./x.*h2(x, n) - h2(x, n + 1);
Bn = zeros(nh, Lsim + 1);
for n = 0:Lsim
    Bn(:,n+1) = 4*pi*1i^n./(1i*kr.^2.*dh2(kr, n));
end

t1 = bpHz(1)*(2^(1/12) - 1); t2 = bpHz(2)*(2^(1/12) - 1);
mask = zeros(nh, 1);
ix = fh >= bpHz(1) + t1 & fh <= bpHz(2) - t2; mask(ix) = 1;
ix = fh > bpHz(1) & fh < bpHz(1) + t1;       mask(ix) = 0.5 - 0.5*cos(pi*(fh(ix) - bpHz(1))/t1);
ix = fh > bpHz(2) - t2 & fh < bpHz(2);       mask(ix) = 0.5 + 0.5*cos(pi*(bpHz(2) - fh(ix))/t2);
mask(1) = 0;

half = zeros(nh, Q);
for p = 1:P
    d = [cosd(echoes(p,2))*cosd(echoes(p,1)); cosd(echoes(p,2))*sind(echoes(p,1)); sind(echoes(p,2))];
    xg = u*d;
    K = zeros(nh, Q);
    Pm1 = ones(Q, 1); P0 = xg;                 % Legendre recursion
    for n = 0:Lsim
        if n == 0, Pn = Pm1;
        elseif n == 1, Pn = P0;
        else
            Pn = ((2*n - 1)*xg.*P0 - (n - 1)*Pm1)/n;
            Pm1 = P0; P0 = Pn;
        end
        K = K + ((2*n + 1)/(4*pi))*(Bn(:,n+1)*Pn.');
    end
    half = half + amp(p).*(K.*exp(-1i*2*pi*fh*delays(p))).*mask;
end
full = [half; conj(half(end-1:-1:2,:))];
hall = real(ifft(full, [], 1));
h = hall(1:Nsamples, :);
truth = struct('az', echoes(:,1), 'el', echoes(:,2), 'relDb', echoes(:,3), 'delays', delays, 'fs', fs);
end

function v = getopt(s, f, d)
if isfield(s, f) && ~isempty(s.(f)), v = s.(f); else, v = d; end
end
