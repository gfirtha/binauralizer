function D = hrtf_df_average(sofa)
%HRTF_DF_AVERAGE  Diffuse-field average of an HRTF set: the power mean of
%its magnitude responses over all directions and both ears - what the set
%gives at its microphones in a diffuse field. It tells where the set was
%measured / how it was equalised, and is the reference of the headphone EQ
%of binaural playback (headphone_db.eq mode 'hrtf').
%
%   D = hrtf_df_average(sofa)        SOFA struct (SOFAload), HRIR data
%   D: f [Hz], dB (1/3-octave smoothed, 0 dB at 1 kHz), gain (the peak
%      2 ... 5 kHz over 1 kHz [dB]), fPeak, flat (max |dB| 200 Hz ... 8 kHz),
%      type 'eardrum' | 'diffuse-field equalised' | 'blocked ear canal',
%      name (the set's database / file name), text (one line)
%   Classification: an eardrum / ear-simulator set has the ear-canal
%   resonance in its diffuse-field average (>= 9 dB at 2-5 kHz); a
%   diffuse-field equalised set (DTF) is flat (within 5 dB); in between:
%   blocked ear canal (the concha only). Directions are weighted equally
%   (dense, roughly uniform grids).

IR = sofa.Data.IR;                                  % [M x R x N]
if ndims(IR) > 3, error('hrtf_df_average:brir', 'BRIR sets (room responses) have no diffuse-field average.'); end
fs = sofa.Data.SamplingRate;
N = 2^nextpow2(max(size(IR, 3), 512));
H = fft(permute(IR, [3 1 2]), N);                  % [N x M x R]
P = mean(mean(abs(H(1:N/2+1, :, :)).^2, 2), 3);    % power mean: directions, ears
f = (0:N/2)'*fs/N;
L = 10*log10(P + eps);
% 1/3-octave smoothing on the power
Ls = L;
for k = 2:numel(f)
    in = abs(log2(f/f(k))) <= 1/6;
    Ls(k) = 10*log10(mean(10.^(L(in)/10)));
end
Ls(1) = Ls(2);
Ls = Ls - interp1(f, Ls, 1000);
b = f >= 2000 & f <= 5000; [g, i] = max(Ls(b)); fb = f(b);
w = f >= 200 & f <= 8000; flat = max(abs(Ls(w)));
if g >= 9, type = 'eardrum'; elseif flat <= 5, type = 'diffuse-field equalised'; else, type = 'blocked ear canal'; end
nm = '';
try, nm = sofa.GLOBAL_DatabaseName; catch, end
if isempty(nm), try, [~, nm] = fileparts(sofa.GLOBAL_Title); catch, end, end
D = struct('f', f, 'dB', Ls, 'gain', g, 'fPeak', fb(i), 'flat', flat, 'type', type, 'name', nm, ...
    'text', sprintf('%s type (diffuse-field average: %+.0f dB at %.1f kHz)', type, g, fb(i)/1e3));
end
