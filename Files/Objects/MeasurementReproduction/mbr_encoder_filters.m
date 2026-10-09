function F = mbr_encoder_filters(R, N, fs, Nf, fcLP, mode)
%MBR_ENCODER_FILTERS  SH-domain encoding filters of a rigid-sphere array
%(port of get_beamforming_filters / get_regularization_filters of M.
%Bubrik's thesis code; renamed, the binauralizer has functions of the same
%names with other signatures).
%
%   F = mbr_encoder_filters(R, N, fs, Nf, fcLP, mode)
%     R     sphere radius [m], N SH order, fs [Hz], Nf FFT length
%     fcLP  low-pass corner of the 6th order roll-off [Hz] (Inf: none)
%     mode  'modal' (default, radial equaliser 1/b_n) | 'DAS'
%   F.Hrad [Nf x N+1]  radial filters (delayed by Nf/2 samples, as before)
%   F.Hreg [Nf x N+1]  crossover regularisation (order-n bands opened at
%                      fc = [90 680 1650 2600]*0.042/R Hz)
%   F.Hlp  [Nf x 1]    1/sqrt(1 + (f/fcLP)^6)
%   F.H    = Hrad.*Hreg.*Hlp, F.c = Hreg.*Hlp (the c_l(f) of Massé et al.)
%   F.freq [Nf x 1], F.ell [1 x (N+1)^2] degree of each ACN channel.

if nargin < 5 || isempty(fcLP), fcLP = Inf; end
if nargin < 6, mode = 'modal'; end
c0 = 343.1;
freq = (0:Nf-1)'/Nf*fs;
k = 2*pi*freq/c0;
x = k*R;

% radial filters
H = zeros(Nf, N+1);
for n = 0:N
    switch lower(mode)
        case 'modal'
            H(:,n+1) = 1/(4*pi*1i^n)*(1i*x.^2.*dsph_hankel2(n, x));
        case 'das'
            jn = sqrt(pi./(2*x)).*besselj(n + 0.5, x);
            H(:,n+1) = 4*pi*(-1i)^(n-1)*x.^2.*dsph_hankel2(n, x).*abs(jn).^2;
        otherwise
            error('mbr_encoder_filters:mode', 'Unknown radial filter mode "%s".', mode);
    end
end
H(~isfinite(H)) = 0;
h = fftshift(ifft(H, [], 1, 'symmetric'), 1);
F.Hrad = fft(h, [], 1);

% crossover regularisation (Linkwitz-like band split, order n opens above fc(n))
fc = [90 680 1650 2600]'*(0.042/R);
fc(fc > 20e3) = 20e3;
f = freq.';
nn = [2 3 4 5];
B = zeros(5, Nf);
B(1,:) = 1./(1 + (f/fc(1)).^nn(1));
for b = 2:4
    B(b,:) = (f/fc(b-1)).^nn(b-1)./(1 + (f/fc(b-1)).^nn(b-1)) .* 1./(1 + (f/fc(b)).^nn(b));
end
B(5,:) = (f/fc(4)).^nn(4)./(1 + (f/fc(4)).^nn(4));
B = B./sum(B, 1);
Hr = flipud(cumsum(flipud(B), 1));            % H_reg_n = sum of the bands >= n
Hr = Hr(1:N+1, :);
hr = ifft(Hr.', Nf, 1, 'symmetric');
F.Hreg = fft(hr);

F.Hlp = 1./sqrt(1 + (freq/fcLP).^6);
F.H = F.Hrad.*F.Hreg.*F.Hlp;
F.c = F.Hreg.*F.Hlp;
F.freq = freq;
F.ell = repelem(0:N, 2*(0:N) + 1);
F.N = N; F.fs = fs; F.Nf = Nf; F.R = R;
end

function dh = dsph_hankel2(n, x)
% derivative of the spherical Hankel function of the 2nd kind
h = @(m) sqrt(pi./(2*x)).*besselh(m + 0.5, 2, x);
dh = -h(n + 1) + (n./x).*h(n);
end
