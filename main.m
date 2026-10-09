%% =========================================================
%  main.m  —  Binaurális Szintézis Analízis
%
%  Futtatási sorrend:
%    1. Config betöltés
%    2. Path beállítás
%    3. Validáció
%    4. Geometria
%    5. Renderer létrehozása
%    6. Fülszámítás frekvenciatartományban
%    7. Vizualizáció
%    8. Optimalizáció (opcionális)
%
%  Kapcsolók a config.m-ben:
%    cfg.drive_method  — 'beamforming', 'wfs', 'lsfi', 'ctc'
%    cfg.ear_method    — 'sphere', 'hrtf', 'freefield'
%    cfg.compare_mode  — true: sphere + freefield egyszerre
%    cfg.optimize      — true: grid search
% =========================================================

clear; close all; clc;
clear classes;


%% --- 1. CONFIG ---
config;


%% --- 2. PATH ---
addpath(genpath('Files/'));
addpath(genpath('+audioapp/'));


%% --- 3. VALIDÁCIÓ ---
assert(ismember(cfg.drive_method, {'beamforming','wfs','lsfi','ctc'}), ...
    'Ismeretlen drive_method: %s', cfg.drive_method);
assert(ismember(cfg.ear_method, {'sphere','hrtf','freefield'}), ...
    'Ismeretlen ear_method: %s', cfg.ear_method);
assert(geo.N >= 2,          'Legalább 2 hangszóró szükséges.');
assert(geo.dx > 0,          'geo.dx pozitív kell legyen.');
assert(geo.head_radius > 0, 'geo.head_radius pozitív kell legyen.');

if ~isfield(cfg, 'compare_mode'), cfg.compare_mode = false; end
if ~isfield(cfg, 'optimize'),     cfg.optimize     = false; end

disp('Validáció: OK');


%% --- 4. GEOMETRIA ---
% loudspeaker + SSD + receiver
% Kiváltja: setup_geometry()

block_len = 2 * cfg.N_taps;

switch geo.type
    case 'linear'
        % Hangszórók a hallgató előtt X=geo.dist-ben, Y tengelyen elosztva
        y_s = linspace(-(geo.N-1)/2*geo.dx, (geo.N-1)/2*geo.dx, geo.N);
        x_s = ones(1, geo.N) * geo.dist;
        z_s = zeros(1, geo.N);
    otherwise
        error('geo.type = ''%s'' nincs implementálva.', geo.type);
end

ls = cell(1, geo.N);
for n = 1:geo.N
    ls{n} = loudspeaker(n, [x_s(n), y_s(n), z_s(n)], [-1,0,0], geo.ls_type);
    ls{n}.set_output(zeros(block_len, 1), cfg.fs);
end

ssd = secondary_source_distribution(ls, [geo.dist, 0]);
rec = receiver(geo.rec_pos, geo.rec_ori);

fprintf('Geometria: %d hangszóró, %s, dx=%.0f mm, dist=%.0f cm\n', ...
    geo.N, geo.type, geo.dx*1000, geo.dist*100);
fprintf('Hallgató:  [%.2f, %.2f, %.2f] m, ori: [%.0f, %.0f, %.0f]\n', ...
    geo.rec_pos, geo.rec_ori);


%% --- 5. RENDERER ---
% beamformer_renderer / ctc_renderer
% Kiváltja: drive_beamforming, drive_wfs, drive_ctc

vs_type = struct('Shape', 'point_source', 'R', 0.025);

switch cfg.drive_method

    case 'beamforming'
        % DAS: síkhullám irányítás cfg.theta_steer szögbe
        % A beamR_pos a jobb fül iránya, beamL_pos a bal fülé
        % A szög a hangszórósor tengelyétől mérve (Y irány)
        R_steer = 2.0;   % irányítási vektor hossza [m]
        center  = [geo.dist, 0, 0];
        beamR   = center + R_steer * [0,  sind(cfg.theta_steer), 0];
        beamL   = center + R_steer * [0, -sind(cfg.theta_steer), 0];

        vs = virtual_source(1, beamR, [-1,0,0], vs_type, 'beamformer');
        vs.set_input(zeros(block_len,1), cfg.fs);

        rs = struct(...
            'BeamformerMode',    'das', ...
            'BeamformerNumTaps', cfg.N_taps, ...
            'Tapering',          cfg.tapering, ...
            'Antialiasing',      cfg.antialiasing, ...
            'BeamLPos',          beamL, ...
            'BeamRPos',          beamR);

        renderer = beamformer_renderer(vs, ssd, rs, rec);
        fprintf('Driving: DAS beamforming (theta=%.0f fok)\n', cfg.theta_steer);
        fprintf('  BeamL: [%.3f, %.3f] m\n', beamL(1), beamL(2));
        fprintf('  BeamR: [%.3f, %.3f] m\n', beamR(1), beamR(2));

    case 'wfs'
        % WFS: pontforrás fókuszálás (cfg.xs, cfg.ys) pozícióban
        vs = virtual_source(1, [cfg.xs, cfg.ys, 0], [-1,0,0], vs_type, 'beamformer');
        vs.set_input(zeros(block_len,1), cfg.fs);

        rs = struct(...
            'BeamformerMode',    'focus', ...
            'BeamformerNumTaps', cfg.N_taps, ...
            'Tapering',          cfg.tapering, ...
            'Antialiasing',      cfg.antialiasing);

        renderer = beamformer_renderer(vs, ssd, rs, rec);
        fprintf('Driving: WFS focus (VS: [%.2f, %.2f] m)\n', cfg.xs, cfg.ys);

    case 'lsfi'
        % LSFI: optimalizált beamforming (Hohnerlein & Ahrens 2017)
        % Explicit null a kontralaterális fülön, WNG megszorítással
        R_steer = 2.0;
        center  = [geo.dist, 0, 0];
        beamR   = center + R_steer * [0,  sind(cfg.theta_steer), 0];
        beamL   = center + R_steer * [0, -sind(cfg.theta_steer), 0];

        vs = virtual_source(1, beamR, [-1,0,0], vs_type, 'beamformer');
        vs.set_input(zeros(block_len,1), cfg.fs);

        rs = struct(...
            'BeamformerMode',    'lsfi', ...
            'BeamformerNumTaps', cfg.N_taps, ...
            'Tapering',          cfg.tapering, ...
            'Antialiasing',      cfg.antialiasing, ...
            'BeamLPos',          beamL, ...
            'BeamRPos',          beamR);

        renderer = beamformer_renderer(vs, ssd, rs, rec);

        % LSFI specifikus beállítások
        renderer.lsfi_null_width  = cfg.lsfi_null_width;
        renderer.lsfi_null_depth  = cfg.lsfi_null_depth;
        renderer.lsfi_fir_length  = cfg.lsfi_fir_length;
        renderer.update_renderer();

        fprintf('Driving: LSFI beamforming\n');
        fprintf('  Null szélesség: %.0f fok, mélység: %.0f dB\n', ...
            cfg.lsfi_null_width, cfg.lsfi_null_depth);

    case 'ctc'
        % CTC renderer
        vs = virtual_source(1, [cfg.xs, cfg.ys, 0], [-1,0,0], vs_type, 'ctc');
        vs.set_input(zeros(block_len,1), cfg.fs);

        switch cfg.ctc_model
            case 'sphere',    plant = 'rigid_sphere';
            case 'hrtf',      plant = 'HRTF';
            case 'freefield', plant = 'point_source';
        end

        ctc_setup = struct('Renderer_setup', struct(...
            'CTCPlantModel', plant, ...
            'CTCVSModel',    plant, ...
            'CTCLambda',     cfg.beta, ...
            'CTCRhead',      geo.head_radius));

        if strcmp(cfg.ctc_model, 'hrtf')
            SOFAstart;
            ctc_setup.HRTF_database = SOFAload(cfg.sofa_file);
        end

        renderer = ctc_renderer(vs, ls, rec, ctc_setup);
        fprintf('Driving: CTC (%s modell)\n', cfg.ctc_model);
end


%% --- 6. FÜLSZÁMÍTÁS FREKVENCIATARTOMÁNYBAN ---
% Az új beamformer_renderer compute_ear_transfer 4 értéket ad vissza:
%   H_LL: L nyaláb → bal fül
%   H_LR: L nyaláb → jobb fül  (crossTalk)
%   H_RL: R nyaláb → bal fül   (crossTalk)
%   H_RR: R nyaláb → jobb fül

fprintf('\nFülszámítás...\n');

% --- Sphere modell ---
if strcmp(cfg.ear_method,'sphere') || cfg.compare_mode
    if isa(renderer, 'beamformer_renderer')
        [f, H_LL_s, H_LR_s, H_RL_s, H_RR_s] = renderer.compute_ear_transfer('analytic');
    else
        [f, H_L_s, H_R_s] = renderer.compute_ear_transfer('analytic');
        H_LL_s = H_L_s; H_RR_s = H_R_s;
        H_LR_s = zeros(size(H_L_s)); H_RL_s = zeros(size(H_R_s));
    end
    fprintf('  Sphere: kész (%d frekvenciapont)\n', numel(f));
end

% --- HRTF modell ---
if strcmp(cfg.ear_method,'hrtf') || cfg.compare_mode
    try
        SOFAstart;
        hrtf_extrap = hrtf_extrapolator(SOFAload(cfg.sofa_file), 'nearest');
        if isa(renderer, 'beamformer_renderer')
            [f, H_LL_h, H_LR_h, H_RL_h, H_RR_h] = renderer.compute_ear_transfer('hrtf', hrtf_extrap);
        else
            [f, H_L_h, H_R_h] = renderer.compute_ear_transfer('hrtf', hrtf_extrap);
            H_LL_h = H_L_h; H_RR_h = H_R_h;
            H_LR_h = zeros(size(H_L_h)); H_RL_h = zeros(size(H_R_h));
        end
        fprintf('  HRTF: kész\n');
    catch ME
        warning('HRTF sikertelen: %s - kihagyva.', ME.message);
    end
end

% --- Freefield modell ---
if strcmp(cfg.ear_method,'freefield') || cfg.compare_mode
    if ~exist('f','var')
        [f, H_LL_s, H_LR_s, H_RL_s, H_RR_s] = renderer.compute_ear_transfer('analytic');
    end
    Nf  = numel(f);
    Nls = geo.N;

    % Fülpontok
    n_hat = geo.rec_ori(:)' / norm(geo.rec_ori);
    up    = [0, 0, 1];
    perp  = cross(n_hat, up);
    if norm(perp) < 1e-10, perp = cross(n_hat, [1,0,0]); end
    perp  = perp / norm(perp);
    x_ear_L = geo.rec_pos + geo.head_radius * perp;
    x_ear_R = geo.rec_pos - geo.head_radius * perp;

    % Beamformer szűrők kinyerése a filter bankból
    if isa(renderer, 'beamformer_renderer') && strcmpi(renderer.wfs_mode, 'lsfi')
        % LSFI: FIR mátrixból
        DfL = fft(renderer.lsfi_fir_L, 2*(Nf-1), 1);
        DfL = DfL(1:Nf, :);
        DfR = fft(renderer.lsfi_fir_R, 2*(Nf-1), 1);
        DfR = DfR(1:Nf, :);
    elseif isa(renderer, 'beamformer_renderer')
        % WFS: filter bankból
        DfL = zeros(Nf, Nls);
        DfR = zeros(Nf, Nls);
        for n = 1:Nls
            h_L = renderer.wfsL.WFS_filter_bank{n}.unfiltered_coefficients{1};
            h_R = renderer.wfsR.WFS_filter_bank{n}.unfiltered_coefficients{1};
            DfL(:,n) = fft(h_L, 2*(Nf-1));
            DfR(:,n) = fft(h_R, 2*(Nf-1));
            DfL(:,n) = DfL(1:Nf, n);
            DfR(:,n) = DfR(1:Nf, n);
        end
    end

    diff_obj = Diffractor();
    diff_obj.scatterer_pos    = geo.rec_pos;
    diff_obj.scatterer_radius = geo.head_radius;
    diff_obj.source_type      = 'point';

    H_LL_ff = zeros(Nf,1); H_LR_ff = zeros(Nf,1);
    H_RL_ff = zeros(Nf,1); H_RR_ff = zeros(Nf,1);

    for n = 1:Nls
        diff_obj.source_pos = ssd.x0(n,:);
        p0 = diff_obj.evaluate_freefield(f', [x_ear_L; x_ear_R]);
        H_LL_ff = H_LL_ff + DfL(:,n) .* p0(1,:).';
        H_LR_ff = H_LR_ff + DfL(:,n) .* p0(2,:).';
        H_RL_ff = H_RL_ff + DfR(:,n) .* p0(1,:).';
        H_RR_ff = H_RR_ff + DfR(:,n) .* p0(2,:).';
    end
    fprintf('  Freefield: kész\n');
end


%% --- 7. VIZUALIZÁCIÓ ---
% results struct összerakása -> plot_ear_vs_frequency

results  = struct();
cfg.f_vec = f';

if exist('H_LL_s','var')
    % ILD = ipsilaterális / kontralaterális arány
    % L nyaláb: bal fül hangos (H_LL), jobb fül halk (H_LR) kellene
    % R nyaláb: jobb fül hangos (H_RR), bal fül halk (H_RL) kellene
    ILD_L = 20*log10(abs(H_LL_s) ./ (abs(H_LR_s) + eps));
    ILD_R = 20*log10(abs(H_RR_s) ./ (abs(H_RL_s) + eps));
    results.ear_sphere.p_ears  = [H_LL_s'; H_RR_s'];
    results.ear_sphere.ILD     = (ILD_L + ILD_R) / 2;
    results.ear_sphere.H_LL    = H_LL_s;
    results.ear_sphere.H_LR    = H_LR_s;
    results.ear_sphere.H_RL    = H_RL_s;
    results.ear_sphere.H_RR    = H_RR_s;
end

if exist('H_LL_h','var')
    ILD_L = 20*log10(abs(H_LL_h) ./ (abs(H_LR_h) + eps));
    ILD_R = 20*log10(abs(H_RR_h) ./ (abs(H_RL_h) + eps));
    results.ear_hrtf.p_ears = [H_LL_h'; H_RR_h'];
    results.ear_hrtf.ILD    = (ILD_L + ILD_R) / 2;
    results.ear_hrtf.H_LL   = H_LL_h;
    results.ear_hrtf.H_LR   = H_LR_h;
end

if exist('H_LL_ff','var')
    ILD_L = 20*log10(abs(H_LL_ff) ./ (abs(H_LR_ff) + eps));
    ILD_R = 20*log10(abs(H_RR_ff) ./ (abs(H_RL_ff) + eps));
    results.ear_freefield.p_ears = [H_LL_ff'; H_RR_ff'];
    results.ear_freefield.ILD    = (ILD_L + ILD_R) / 2;
end

plot_ear_vs_frequency(results, cfg);


%% --- 8. OPTIMALIZÁCIÓ (opcionális) ---
if cfg.optimize
    if ~isa(renderer, 'beamformer_renderer')
        warning('Optimalizáció csak beamformer_renderer esetén elérhető.');
    else
        fprintf('\nGrid search: %dx%d, %d-%d Hz...\n', ...
            cfg.grid_Nx, cfg.grid_Ny, cfg.f_opt_low, cfg.f_opt_high);

        result = renderer.run_optimal_search(...
            cfg.f_opt_low, cfg.f_opt_high, 'analytic', [], ...
            cfg.grid_Nx, cfg.grid_Ny);

        fprintf('Optimális pozíció: [%.3f, %.3f] m\n', result.opt_pos);
        fprintf('Érték: %.1f dB\n', result.opt_val);

        % Heatmap
        figure('Name','Grid search','Color','w');
        imagesc(result.X(1,:), result.Y(:,1), result.suppression_map);
        set(gca,'YDir','normal'); colorbar; colormap('jet');
        hold on;
        plot(ssd.x0(:,1), ssd.x0(:,2), 'ws', 'MarkerSize',8, ...
             'MarkerFaceColor','w', 'DisplayName','Hangszórók');
        plot(geo.rec_pos(1), geo.rec_pos(2), 'k^', 'MarkerSize',10, ...
             'MarkerFaceColor','k', 'DisplayName','Hallgató');
        plot(result.opt_pos(1), result.opt_pos(2), 'w*', ...
             'MarkerSize',15, 'LineWidth',2, 'DisplayName','Optimum');
        xlabel('x [m]'); ylabel('y [m]');
        title(sprintf('Optimum: [%.2f, %.2f] m — %.1f dB', ...
            result.opt_pos(1), result.opt_pos(2), result.opt_val));
        legend('Location','best');

        % ILD összehasonlítás eredeti vs optimális
        if exist('H_LL_s','var')
            renderer.beamR_pos = [result.opt_pos, 0];
            renderer.beamL_pos = [-result.opt_pos(1), result.opt_pos(2), 0];
            renderer.update_renderer();
            [~, H_LL_opt, H_LR_opt, H_RL_opt, H_RR_opt] = ...
                renderer.compute_ear_transfer('analytic');
            ILD_opt = (20*log10(abs(H_LL_opt)./(abs(H_LR_opt)+eps)) + ...
                       20*log10(abs(H_RR_opt)./(abs(H_RL_opt)+eps))) / 2;

            figure('Name','Eredeti vs Optimális','Color','w');
            hold on; grid on;
            semilogx(f, results.ear_sphere.ILD, 'b-', 'LineWidth',2, ...
                     'DisplayName','Eredeti pozíció');
            semilogx(f, ILD_opt, 'r--', 'LineWidth',2, ...
                     'DisplayName','Optimális pozíció');
            yline(0,'k-','Alpha',0.3);
            xlabel('Frekvencia [Hz]'); ylabel('ILD [dB]');
            title('ILD összehasonlítás'); legend('Location','best');
            set(gca,'XScale','log');
            xlim(cfg.f_plot_range);
        end
    end
end

disp('Kész.');
