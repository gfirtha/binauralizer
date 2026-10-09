addpath(genpath('Files/'));
addpath(genpath('+audioapp/'));

fs = 44100;
block_len = 1024;
Nls = 8;
dx = 0.05;
dist = 0.5;

y_s = linspace(-(Nls-1)/2*dx, (Nls-1)/2*dx, Nls);
x_s = ones(1, Nls) * dist;
z_s = zeros(1, Nls);

ls = cell(1, Nls);
ls_type = struct('Shape', 'point_source', 'R', 0.025);
for n = 1:Nls
    ls{n} = loudspeaker(n, [x_s(n), y_s(n), z_s(n)], [-1, 0, 0], ls_type);
    ls{n}.set_output(zeros(block_len, 1), fs);
end

ssd = secondary_source_distribution(ls, [dist, 0]);
rec = receiver([0, 0, 0], [1, 0, 0]);

vs_type = struct('Shape', 'point_source', 'R', 0.025);
vs = virtual_source(1, [dist, 0.15, 0], [-1, 0, 0], vs_type, 'beamformer');
vs.set_input(zeros(block_len, 1), fs);

rs = struct( ...
    'BeamformerMode', 'focus', ...
    'BeamformerNumTaps', 512, ...
    'Tapering', 0.3, ...
    'Antialiasing', true, ...
    'LowFreqCut', false, ...
    'BeamLPos', [dist, -0.15, 0], ...
    'BeamRPos', [dist,  0.15, 0]);

r = beamformer_renderer(vs, ssd, rs, rec);

[f0, HLL0] = local_get_hll(r);

r_base_off = r.wfsL.driving_spectrum;
r.LowFreqCut = true;
r.update_renderer();
[f1, HLL1] = local_get_hll(r);
r_base_on = r.wfsL.driving_spectrum;

low_mask = (f0 >= 50) & (f0 <= 1000);
delta_lowcut_db = local_max_finite_diff(local_db(r_base_on(low_mask, :)), local_db(r_base_off(low_mask, :)));
fprintf('LOWCUT_DRIVING_MAX_DB_DIFF=%.6f\n', delta_lowcut_db);

eq = zeros(256, Nls);
b = fir1(64, 1200/(fs/2), 'low');
for n = 1:Nls
    eq(1:numel(b), n) = b(:);
end
tmp_eq_file = fullfile(tempdir, 'tmp_beamformer_eq.mat');
hinv = eq; %#ok<NASGU>
save(tmp_eq_file, 'hinv', 'fs');

r.loadSpeakerEQ(tmp_eq_file);
r.SpeakerEQEnabled = false;
[f2, HLL2] = local_get_hll(r);
r_eq_off = local_effective_df(r);
r.SpeakerEQEnabled = true;
[f3, HLL3] = local_get_hll(r);
r_eq_on = local_effective_df(r);

delta_eq_db = local_max_finite_diff(local_db(r_eq_on(low_mask, :)), local_db(r_eq_off(low_mask, :)));
fprintf('SPK_EQ_DRIVING_MAX_DB_DIFF=%.6f\n', delta_eq_db);

fprintf('ANALYTIC_HLL_FINITE_COUNT=%d\n', sum(isfinite(HLL3)));
fprintf('BASE_DF_FINITE_COUNT=%d\n', sum(isfinite(r_base_on(:))));
fprintf('EQ_DF_FINITE_COUNT=%d\n', sum(isfinite(r_eq_on(:))));

if ~(delta_lowcut_db > 0.05)
    error('ValidationFailed:LowCut', ...
        'Low-freq pre-EQ had no measurable effect in ear-transfer driving spectra.');
end

if ~(delta_eq_db > 0.05)
    error('ValidationFailed:SpeakerEQ', ...
        'Speaker EQ had no measurable effect in ear-transfer driving spectra.');
end

fprintf('EAR_TRANSFER_VALIDATION_OK\n');

function [f, HLL] = local_get_hll(r)
    [f, HLL, ~, ~, ~] = r.compute_ear_transfer('analytic');
end

function y = local_db(x)
    y = 20 * log10(abs(x) + 1e-15);
end

function Df = local_effective_df(r)
    Df = r.wfsL.driving_spectrum;
    if r.SpeakerEQEnabled && ~isempty(r.eqImpulseResponses)
        Nf = size(Df, 1);
        Neq = min(size(r.eqImpulseResponses, 2), size(Df, 2));
        Nfft = 2 * (Nf - 1);
        H_eq = fft(r.eqImpulseResponses(:, 1:Neq), Nfft, 1);
        H_eq = H_eq(1:Nf, :);
        Df(:, 1:Neq) = Df(:, 1:Neq) .* H_eq;
    end
end

function out = local_max_finite_diff(a, b)
    d = abs(a - b);
    d = d(isfinite(d));
    if isempty(d)
        out = NaN;
    else
        out = max(d);
    end
end
