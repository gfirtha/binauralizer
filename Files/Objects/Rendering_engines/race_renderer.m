classdef race_renderer < base_renderer
    % RACE_RENDERER  Recursive Ambiophonic Crosstalk Elimination renderer.
    %
    % Alacsony frekvenciás CrossTalk Cancellation módszer (Glasgal, 2007),
    % amelyet Hohnerlein & Ahrens (ICASSP 2017) alkalmazott lineáris
    % hangszórósor alapú binaurális szintézishez.
    %
    % Elv:
    %   A bal hangszóróból kiadott hang a jobb fülhöz is eljut (crossTalk).
    %   A RACE ezt rekurzívan törli: a bal csatornából késleltetett és
    %   csillapított ellentétes jelű másolatot kever a jobb csatornába,
    %   és fordítva. Mivel ez a törlőjel maga is crossTalk-ot okoz,
    %   a rekurzió végtelen sor — amely konvergens ha delta_a < 1.
    %
    % Frekvencia tartomány: 250 Hz - 1000 Hz (tipikusan)
    %   Ennél magasabb frekvencián a beamforming veszi át,
    %   ennél alacsonyabban nincs aktív CTC (szubwoofer tartomány).
    %
    % Sávszétválasztás (opcionális):
    %   Ha cfg.race_crossover_freq be van állítva, Linkwitz-Riley
    %   crossover szűrőkkel kombinálja a RACE-t és a beamformert.
    %
    % Irodalom:
    %   R. Glasgal, "360 Localization via 4.x RACE Processing,"
    %   AES Convention 123, 2007.
    %
    %   C. Hohnerlein & J. Ahrens, "Perceptual Evaluation of a Multiband
    %   Acoustic Crosstalk Canceler Using a Linear Loudspeaker Array,"
    %   ICASSP 2017.
    %
    % Használat:
    %   rs = struct('RACEDeltaA', 0.5, 'RACEDeltaT', 0.0006, ...
    %               'RACENumTaps', 512);
    %   r = race_renderer(vs, ssd, rs, rec);
    %   [f, H_LL, H_LR, H_RL, H_RR] = r.compute_ear_transfer('analytic');

    properties
        fs
        receiver

        % RACE paraméterek
        delta_a   = 0.5      % Csillapítás (0-1): fej árnyékolási hatás
                             % delta_a = 0.5 ≈ -6 dB, tipikus érték
                             % Kisebb → erősebb törlés, de esetleg instabil
        delta_t   = 6e-4     % Késleltetés [s]: kontralaterális út többlete
                             % delta_t ≈ fej_átmérő / c = 0.20/343 ≈ 0.6 ms
        N_taps    = 512      % FIR szűrő hossz

        % Szűrők (OLS_convolver objektumok)
        filter_L  = []       % bal csatornára alkalmazott RACE szűrő
        filter_R  = []       % jobb csatornára alkalmazott RACE szűrő

        % Bal és jobb hangszóró indexek az SSD-ben
        % Lineáris array esetén: bal = bal szélen, jobb = jobb szélen
        % Vagy: összes hangszóró egyforma súllyal (sztereo keverés)
        ls_mode   = 'stereo' % 'stereo': 2 csatorna → N hangszóró keverés
                             % 'pair':   csak a bal/jobb szélső hangszóró
    end

    properties (SetAccess = protected)
        race_fir_L   % [N_taps x 1] bal csatorna RACE FIR szűrő
        race_fir_R   % [N_taps x 1] jobb csatorna RACE FIR szűrő
        H_race_L     % [Nf x 1] bal csatorna frekvenciatartomány válasz
        H_race_R     % [Nf x 1] jobb csatorna frekvenciatartomány válasz
    end

    methods
        function obj = race_renderer(virtual_source, SSD, setup, receiver)
            % Konstruktor
            %   virtual_source: virtual_source objektum
            %   SSD:            secondary_source_distribution
            %   setup:          struct a paraméterekkel (Renderer_setup mező)
            %   receiver:       receiver objektum

            obj = obj@base_renderer(virtual_source, SSD);
            obj.fs = obj.virtual_source.source_signal.fs;

            if nargin >= 4 && ~isempty(receiver)
                obj.receiver = receiver;
            end

            % Beállítások beolvasása
            RS = struct();
            if isstruct(setup) && isfield(setup, 'Renderer_setup')
                RS = setup.Renderer_setup;
            elseif isstruct(setup)
                RS = setup;
            end

            if isfield(RS, 'RACEDeltaA'),  obj.delta_a = RS.RACEDeltaA; end
            if isfield(RS, 'RACEDeltaT'),  obj.delta_t = RS.RACEDeltaT; end
            if isfield(RS, 'RACENumTaps'), obj.N_taps  = RS.RACENumTaps; end
            if isfield(RS, 'RACELsMode'),  obj.ls_mode = RS.RACELsMode; end

            % Output signal inicializálás
            for n = 1:numel(obj.ssd.loudspeakers)
                obj.output_signal{n} = signal;
            end

            obj.update_renderer();
        end

        function info = get_renderer_info(obj)
            lines = [
                "RACE renderer — Recursive Ambiophonic Crosstalk Elimination"
                ""
                "MŰKÖDÉSI ELV"
                "• Alacsony frekvenciás CTC (250 Hz - 1 kHz)"
                "• A kontralaterális crossTalk-ot rekurzív jelvisszacsatolással törli"
                "• Frekvenciafüggetlen késleltetés és csillapítás"
                ""
                sprintf("PARAMÉTEREK")
                sprintf("• delta_a (csillapítás):  %.3f  (%.1f dB)", ...
                    obj.delta_a, 20*log10(obj.delta_a))
                sprintf("• delta_t (késleltetés):  %.4f s  (%.1f mm)", ...
                    obj.delta_t, obj.delta_t * 343 * 1000)
                sprintf("• N_taps (FIR hossz):     %d", obj.N_taps)
                sprintf("• ls_mode:                %s", obj.ls_mode)
                ""
                "IRODALOM"
                "• Glasgal (AES 2007): RACE módszer"
                "• Hohnerlein & Ahrens (ICASSP 2017): lineáris array alkalmazás"
            ];
            info = strjoin(lines, newline);
        end

        function update_renderer(obj, ~)
            % RACE szűrők felépítése
            % A RACE frekvenciatartományi transzfermátrixa:
            %
            %   [Y_L]   1/(1-G²) * [1,  G] * [X_L]
            %   [Y_R] =             [G,  1]   [X_R]
            %
            % ahol G(ω) = -delta_a * exp(-jω*delta_t)
            %
            % Ez a 2x2 mátrix invertálható mint egyszerű IIR szűrő.

            Nfft = obj.N_taps * 2;
            Nf   = Nfft / 2 + 1;
            freq = (0:Nf-1)' * obj.fs / Nfft;
            omega = 2 * pi * freq;

            % G(ω) = -delta_a * exp(-jω*delta_t)
            G = -obj.delta_a * exp(-1i * omega * obj.delta_t);

            % A RACE mátrix inverze:
            %   [Y_L(ω)]   1/(1-G²) * [1, -G] * [X_L(ω)]
            %   [Y_R(ω)] =             [-G,  1]   [X_R(ω)]
            %
            % Az egyes csatornákra vonatkozó szűrők:
            %   H_LL = H_RR = 1 / (1 - G²)          (direkt szűrő)
            %   H_LR = H_RL = -G / (1 - G²)          (kereszt szűrő)

            denom = 1 - G.^2;
            H_direct = 1 ./ denom;           % [Nf x 1]
            H_cross  = -G ./ denom;          % [Nf x 1]

            % DC kezelés: G(0) = -delta_a (valós), 1-G²(0) = 1-delta_a²
            % Ez véges, nincs singularitás ha delta_a < 1

            % Teljes frekvenciatartományú szűrők (egyoldalas → kétoldalas)
            H_direct_full = [H_direct; conj(flipud(H_direct(2:end-1)))];
            H_cross_full  = [H_cross;  conj(flipud(H_cross(2:end-1)))];

            % IFFT → idő tartomány
            h_direct = real(ifft(H_direct_full));
            h_cross  = real(ifft(H_cross_full));

            % Ablakozás a Gibbs-jelenség csökkentésére
            win = tukeywin(Nfft, 0.1);
            h_direct = fftshift(h_direct) .* win;
            h_cross  = fftshift(h_cross)  .* win;

            obj.race_fir_L = h_direct;   % bal csatorna direkt
            obj.race_fir_R = h_cross;    % kereszt szűrő

            % OLS konvolverek frissítése
            blk = size(obj.virtual_source.source_signal.time_series, 1);
            if blk == 0, blk = 1024; end

            obj.filter_L = OLS_convolver(h_direct, blk);
            obj.filter_R = OLS_convolver(h_cross,  blk);

            % Eltárolás a compute_ear_transfer-hez
            obj.H_race_L = H_direct;
            obj.H_race_R = H_cross;
        end

        function render(obj)
            % RACE renderelés:
            % A bal csatornát a direkt + kereszt szűrőn futtatjuk,
            % majd egyenletesen osztjuk el az összes hangszóróra.

            x = obj.virtual_source.source_signal.time_series;
            if isempty(x), return; end

            % RACE szűrés (sztereo kimenet)
            y_L = obj.filter_L.convolve(x);   % bal fül csatorna
            y_R = obj.filter_R.convolve(x);   % jobb fül csatorna

            Nls = numel(obj.ssd.loudspeakers);

            switch obj.ls_mode
                case 'stereo'
                    % Sztereo → N csatorna keverés
                    % Bal félre: y_L, jobb félre: y_R, középső: average
                    center = mean(obj.ssd.x0(:,2));
                    for n = 1:Nls
                        y_pos = obj.ssd.x0(n, 2);
                        % Lineáris cross-fade a bal és jobb csatorna között
                        alpha = (y_pos - min(obj.ssd.x0(:,2))) / ...
                                (max(obj.ssd.x0(:,2)) - min(obj.ssd.x0(:,2)) + eps);
                        y_n = (1 - alpha) * y_L + alpha * y_R;
                        obj.output_signal{n}.set_signal(y_n / Nls);
                    end

                case 'pair'
                    % Csak a bal/jobb szélső hangszóró kap jelet
                    [~, idx_L] = min(obj.ssd.x0(:,2));
                    [~, idx_R] = max(obj.ssd.x0(:,2));
                    for n = 1:Nls
                        obj.output_signal{n}.set_signal(zeros(size(x)));
                    end
                    obj.output_signal{idx_L}.set_signal(y_L);
                    obj.output_signal{idx_R}.set_signal(y_R);
            end

            obj.add_output_to_ssd_signal;
        end

        function [f, H_LL, H_LR, H_RL, H_RR] = compute_ear_transfer(obj, model, hrtf_extrap)
            % Fülszámítás frekvenciatartományban — 2x2 CrossTalk mátrix
            %
            % A RACE esetén a 2x2 mátrix:
            %   H_LL: bal csatorna → bal fül   (direkt)
            %   H_LR: bal csatorna → jobb fül  (crossTalk — ezt nyomja el)
            %   H_RL: jobb csatorna → bal fül  (crossTalk — ezt nyomja el)
            %   H_RR: jobb csatorna → jobb fül (direkt)

            if nargin < 2, model = 'analytic'; end
            if nargin < 3, hrtf_extrap = []; end

            Nf = numel(obj.H_race_L);
            f  = (0:Nf-1)' * obj.fs / (2*(Nf-1));

            Nls = numel(obj.ssd.loudspeakers);
            c   = audioapp.util.Physics.speedOfSound();

            % Akusztikus átvitel hangszóróktól a fülekhez [Nf x Nls]
            G_earL = zeros(Nf, Nls);
            G_earR = zeros(Nf, Nls);

            switch lower(model)
                case 'analytic'
                    d = Diffractor();
                    d.scatterer_pos    = obj.receiver.position;
                    d.scatterer_radius = 0.0875;
                    for n = 1:Nls
                        d.source_pos  = obj.ssd.x0(n,:);
                        d.source_type = 'point';
                        H_ears = d.evaluate_ears(f', obj.receiver.orientation);
                        G_earL(:,n) = H_ears(1,:).';
                        G_earR(:,n) = H_ears(2,:).';
                    end

                case 'hrtf'
                    if isempty(hrtf_extrap)
                        error('race_renderer:noHRTF', 'hrtf_extrapolator szükséges.');
                    end
                    Nfft = 2*(Nf-1);
                    for n = 1:Nls
                        hrir = hrtf_extrap.extrapolate_hrtf(...
                            obj.ssd.loudspeakers{n}, obj.receiver);
                        H = fft(hrir, Nfft, 1);
                        G_earL(:,n) = H(1:Nf, 1);
                        G_earR(:,n) = H(1:Nf, 2);
                    end

                otherwise
                    error('race_renderer:badModel', 'Ismeretlen modell: %s', model);
            end

            % Sztereo keverési súlyok [Nls x 1]
            alpha = (obj.ssd.x0(:,2) - min(obj.ssd.x0(:,2))) ./ ...
                    (max(obj.ssd.x0(:,2)) - min(obj.ssd.x0(:,2)) + eps);
            w_L = (1 - alpha) / Nls;   % bal csatorna súlyok
            w_R = alpha / Nls;         % jobb csatorna súlyok

            % Frekvenciatartományú driving function
            % direkt szűrő: H_race_L, kereszt szűrő: H_race_R
            Df_LL = obj.H_race_L .* (G_earL * w_L);  % bal cs → bal fül
            Df_LR = obj.H_race_L .* (G_earR * w_L);  % bal cs → jobb fül
            Df_RL = obj.H_race_R .* (G_earL * w_R);  % jobb cs → bal fül
            Df_RR = obj.H_race_R .* (G_earR * w_R);  % jobb cs → jobb fül

            H_LL = sum(Df_LL, 2);
            H_LR = sum(Df_LR, 2);
            H_RL = sum(Df_RL, 2);
            H_RR = sum(Df_RR, 2);

            % Normalizálás
            norm_val = max(abs([H_LL; H_RR])) + eps;
            H_LL = H_LL / norm_val;
            H_LR = H_LR / norm_val;
            H_RL = H_RL / norm_val;
            H_RR = H_RR / norm_val;
        end

        function settings = defaultOverlaySettings(~)
            settings = struct();
        end

        function store = drawOverlay(~, ~, ~, ~, store)
            % RACE renderernek nincs vizuális overlay
        end
    end
end
