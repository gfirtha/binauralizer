%% =========================================================
%  config.m  —  Binaurális Szintézis Konfiguráció
%
%  A binaulizerből importált osztályokat használja:
%    loudspeaker, secondary_source_distribution, receiver
%    beamformer_renderer, ctc_renderer, Diffractor
% =========================================================

%% --- FIZIKAI ÁLLANDÓK ---
cfg.c    = 343;
cfg.rho0 = 1.21;


%% --- HANGSZÓRÓSOR GEOMETRIA ---
geo.N    = 8;
geo.dx   = 0.05;        % Hangszórók közötti távolság [m]
geo.dist = 0.5;         % Sor távolsága a hallgatótól X irányban [m]
geo.type = 'linear';

geo.ls_type = struct('Shape', 'point_source', 'R', 0.025);
%   Shape: 'point_source', 'circular_piston', 'two_way_speaker'


%% --- FEJMODELL ---
geo.head_radius = 0.0875;   % Fejsugár [m]

% Hallgató az origóban, X irányba néz (a hangszórósor felé)
geo.rec_pos = [0, 0, 0];
geo.rec_ori = [1, 0, 0];


%% --- FREKVENCIA BEÁLLÍTÁSOK ---
cfg.fs     = 44100;
cfg.N_taps = 512;


%% --- MÓDSZER VÁLASZTÁS ---
cfg.drive_method = 'beamforming';
%   'beamforming' — DAS síkhullám irányítás
%   'wfs'         — WFS pontforrás fókuszálás
%   'lsfi'        — LSFI optimalizált beamforming (Hohnerlein & Ahrens 2017)
%   'ctc'         — Crosstalk Cancellation

cfg.ear_method = 'sphere';
%   'sphere'    — merev gömb modell (Diffractor)
%   'hrtf'      — mért HRTF adatbázis
%   'freefield' — szabad tér 1/r referencia

cfg.compare_mode = false;
%   false — csak cfg.ear_method fut
%   true  — sphere + freefield egyszerre


%% --- NYALÁB POZÍCIÓK ---
% A beamformer_renderer kettős nyalábot használ:
%   beamL_pos → bal fül irányába sugároz
%   beamR_pos → jobb fül irányába sugároz
%
% Ha [0,0,0]-ban hagyod, a renderer automatikusan
% ±15 cm-es szimmetrikus elhelyezést használ a VS körül.
cfg.beamL_pos = [0, 0, 0];   % [0,0,0] = automatikus
cfg.beamR_pos = [0, 0, 0];   % [0,0,0] = automatikus

% DAS mód: irányítási szög [fok]
% 0 = egyenesen előre (X irány), + = Y+ oldal
cfg.theta_steer = 15;

% WFS / LSFI virtuális forrás pozíciója
cfg.xs = 0;      % [m]
cfg.ys = 0.3;    % [m]

% Közös beállítások
cfg.tapering     = 0.3;
cfg.antialiasing = false;


%% --- LSFI BEÁLLÍTÁSOK ---
% Csak cfg.drive_method = 'lsfi' esetén aktív
cfg.lsfi_null_width = 9;     % Null szélesség [fok] a kontralaterális fül körül
cfg.lsfi_null_depth = -40;   % Elvárt null mélység [dB]
cfg.lsfi_fir_length = 1024;  % LSFI FIR szűrő hossz [minták]
%
% Fizikai értelmezés:
%   lsfi_null_width  — mennyire széles legyen a "néma zóna" a stop irány körül
%                      kisebb = pontosabb null, de érzékenyebb
%   lsfi_null_depth  — milyen mélyre nyomja a stop irányt
%                      -40 dB ≈ 100-szoros csillapítás
%   lsfi_fir_length  — a szűrő frekvencia-felbontását határozza meg
%                      hosszabb = pontosabb, de lassabb számítás


%% --- CTC BEÁLLÍTÁSOK ---
% Csak cfg.drive_method = 'ctc' esetén aktív
cfg.ctc_model = 'sphere';   % 'sphere', 'hrtf', 'freefield'
cfg.beta      = 1e-5;       % Tikhonov regularizáció


%% --- HRTF ---
cfg.sofa_file = 'HRTFs/BuK_ED_corr.sofa';


%% --- VIZUALIZÁCIÓ ---
cfg.f_plot_range = [100, 6000];   % Megjelenített frekvencia tartomány [Hz]
cfg.dB_floor     = -40;
cfg.dB_ceiling   =  0;


%% --- OPTIMALIZÁCIÓ ---
cfg.optimize   = false;   % true: grid search fut az optimális pozícióhoz
cfg.f_opt_low  = 500;     % Sáv alja [Hz]
cfg.f_opt_high = 5000;    % Sáv teteje [Hz]
cfg.grid_Nx    = 40;
cfg.grid_Ny    = 40;
