classdef sound_scene_simulator < handle
    %SCENE_SIMULATOR Summary of this class goes here
    %   Detailed explanation goes here

    properties
        temporary_sound_scene
        gui
        simulation_type
        harmonic_type
        f0
        t
        t0
        pulse_width
        x
        y
        fig
        simulation_axes
        simulation_result
        p_field

        % --- UI / state for visualization
        viewMode            % 'real' | 'mag' | 'amperr' | 'ref'
        btnGroup            % uibuttongroup handle
        timeFreqSlider      % handle to the bottom slider
        lastSliderVal = 0.5 % remember slider position to refresh view

        % --- Colorbar + caxis UI
        cbarHandle           % colorbar handle
        cminEdit             % uicontrol handle (min)
        cmaxEdit             % uicontrol handle (max)
        cLabel               % uicontrol text label "cLim:"
        % Slider readout + checkbox
        sliderValueLabel    % uicontrol text showing Hz or ms
        rendererCheckbox    % uicontrol checkbox
        rendererCopies = [] % array of copied handles into simulation axes


        % --- Monitor update throttle (timer-based)
        monitorTimer          % timer handle
        monitorRate (1,1) double = 20   % target FPS (refreshes per second)

        masterEQ   % audioapp.dsp.MasterEQ (opcionális)

        % --- Head diffraction
        diffractor          % Diffractor object (lazy-init)
        includeHeadEffect   = false   % checkbox state
        headCheckbox        % uicontrol handle

        % --- Fast field evaluation + animation
        dx_sim     = 0.02   % grid step [m]; display is bilinear-interpolated
        engine              % field_engine (cached geometry)
        Zharm      = []     % complex harmonic field, Re{Zharm} = field at t0
        playButton          % uicontrol toggle (harmonic / impulse)
        animTimer           % timer driving the animation
        animFps    = 30
        animPhase  = 0      % harmonic: animation phase [rad]
        animFrame  = 0      % frame counter (sparse UI refresh)
        animSignals = []    % impulse: loudspeaker signals of the running animation
        tOnset     = 0      % impulse: time the first loudspeaker sounds [s]
        impulseClim = []    % impulse: fixed colour scale (set on first plot)
    end

    methods
        function obj = sound_scene_simulator(varargin)
            %SCENE_SIMULATOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.temporary_sound_scene = varargin{1}.duplicate_scene;
            obj.gui = varargin{2};
            obj.gui.simulatorHandle = obj;

            obj.simulation_type = varargin{3};
            obj.harmonic_type = 'real';
            fs = obj.temporary_sound_scene.virtual_sources{1}.source_signal.fs;
            N = length(obj.temporary_sound_scene.virtual_sources{1}.source_signal.time_series);

            t0_default = N/2/fs;

            if nargin >= 5 && ~isempty(varargin{5})
                obj.masterEQ = varargin{5};
            end

            switch obj.simulation_type
                case 'harmonic'
                    obj.t0 = t0_default;
                    obj.f0 = varargin{4};
                case 'impulse'
                    obj.t0 = t0_default;
                    obj.pulse_width = varargin{4};
                case 'monitor'
                    obj.t0 = 0;
            end
            dx = obj.dx_sim;
            obj.x = [obj.gui.main_axes.XLim(1):dx:obj.gui.main_axes.XLim(2)]';
            obj.y = [obj.gui.main_axes.YLim(1):dx:obj.gui.main_axes.YLim(2)]';
            xr = [0,0];
            obj.engine = field_engine(audioapp.util.Physics.speedOfSound(), fs);
            obj.engine.set_grid(obj.x, obj.y);
            obj.fig = figure('units','normalized','position',[0.15 0.2 0.45 0.6], ...
                'Color','w', 'CloseRequestFcn', @sound_scene_simulator.closeFigure);
            % NB: not a method of obj - the figure's callbacks may hold the
            % last references to this object, and destroying it from inside
            % one of its own methods crashed R2023a on window close.

            set(obj.fig,'defaulttextinterpreter','latex')

            obj.simulation_axes = axes('Units','normalized','Position',[0.0,0.2,0.7,0.7]);
            % image object: CData updates are far cheaper than pcolor +
            % shading interp; bilinear display smooths the 2 cm grid
            obj.simulation_result = imagesc(obj.simulation_axes, obj.x, obj.y, ...
                zeros(numel(obj.y), numel(obj.x)));
            try, obj.simulation_result.Interpolation = 'bilinear'; catch, end
            set(obj.simulation_axes, 'YDir', 'normal', 'Layer', 'top');
            colormap(obj.simulation_axes, sound_scene_simulator.fieldColormap('diverging'));
            axis equal tight
            ftsize = 12;  %case1
            %ftsize = 16; %case2
            copyobj(obj.gui.receiver, obj.simulation_axes);
            for n = 1 : length(obj.gui.virtual_source_points)
                if isvalid(obj.gui.virtual_source_points{n})
                    copyobj(obj.gui.virtual_source_points{n}, obj.simulation_axes);
                end
            end
            for n = 1 : length(obj.gui.loudspeaker_points)
                copyobj(obj.gui.loudspeaker_points{n}, obj.simulation_axes);
            end
            xlim(obj.gui.main_axes.XLim)
            ylim(obj.gui.main_axes.YLim)
            clim([-1,1])
            xlabel( '$x$ [m]', 'FontSize', ftsize );
            ylabel( '$y$ [m]', 'FontSize', ftsize );
            set(gca,'FontName','Times New Roman');
            allAxesInFigure = findall(obj.fig,'type','axes');
            b = get(gca,'XTickLabel');
            set(allAxesInFigure,'XTickLabel',b,'FontSize',ftsize);
            set(gca,'XTickLabelMode','auto')
            b = get(gca,'YTickLabel');
            set(allAxesInFigure,'YTickLabel',b,'FontSize',ftsize);
            set(gca,'YTickLabelMode','auto')
            hold on
            plot(xr(1),xr(2), 'k+', 'MarkerSize', 20);

            % --- Colorbar
            obj.cbarHandle = colorbar(obj.simulation_axes);
            obj.cbarHandle.Location = 'eastoutside';
            obj.cbarHandle.Box = 'on';

            % Position helpers (normalized to the figure)
            axPos   = obj.simulation_axes.Position;   % [x y w h]
            cbWidth = 0.03;                           % slim colorbar width
            gap     = -0.1;                          % gap between axes and colorbar

            % Resize axes a tad to make room for cbar + edits
            axRightEdge = axPos(1) + axPos(3);
            extraRight  = cbWidth + 0.02;             % space for cbar + 2 edit boxes
            obj.simulation_axes.Position = [axPos(1), axPos(2), 1 - axPos(1) - extraRight, axPos(4)];

            % Recompute after axes shrink
            axPos = obj.simulation_axes.Position;

            % Park the colorbar
            obj.cbarHandle.Position(1) = axPos(1) + axPos(3) + gap;
            obj.cbarHandle.Position(3) = cbWidth;

            % --- cLim edit controls stacked to the right of the colorbar
            editPanelX = obj.cbarHandle.Position(1) + obj.cbarHandle.Position(3) + gap +0.105;
            editW      = 0.06;  % width of the edit boxes
            editH      = 0.05; % height of the edit boxes

            obj.cLabel = uicontrol('Style','text','Parent',obj.fig,'Units','normalized', ...
                'Position',[editPanelX, axPos(2) + axPos(4) - editH, editW, editH], ...
                'String','cLim:', 'HorizontalAlignment','left');

            obj.cminEdit = uicontrol('Style','edit','Parent',obj.fig,'Units','normalized', ...
                'Position',[editPanelX, obj.cbarHandle.Position(2) - 0.025, editW, editH], ...
                'String','', 'BackgroundColor','white', ...
                'Callback', @(src,~) obj.onClimEdited(src,'min'));

            obj.cmaxEdit = uicontrol('Style','edit','Parent',obj.fig,'Units','normalized', ...
                'Position',[editPanelX, obj.cbarHandle.Position(2) + obj.cbarHandle.Position(4) - editH + 0.025, editW, editH], ...
                'String','', 'BackgroundColor','white', ...
                'Callback', @(src,~) obj.onClimEdited(src,'max'));


            %% Add interactive slider for time (impulse) or frequency (harmonic)
            switch obj.simulation_type
                case 'impulse'
                    uicontrol('Style','text','Parent',obj.fig,'Units','normalized',...
                        'Position',[0.05 0.05 0.2 0.03],'String','Time (s)');
                    obj.timeFreqSlider = uicontrol('Style','slider','Parent',obj.fig,'Units','normalized',...
                        'Position',[0.275 0.05 0.4 0.03],'Min',-2e-2+t0_default,'Max',2e-2+t0_default,'Value',t0_default,...
                        'Callback', @(src,~) obj.update_time(src.Value));

                    % -> readout to the right of slider
                    obj.sliderValueLabel = uicontrol('Style','text','Parent',obj.fig,'Units','normalized',...
                        'Position',[0.68 0.05 0.12 0.03], 'HorizontalAlignment','left', ...
                        'String','');  % filled by updateSliderLabel below

                case 'harmonic'
                    uicontrol('Style','text','Parent',obj.fig,'Units','normalized',...
                        'Position',[0.05 0.05 0.2 0.03],'String','Frequency (Hz)');
                    obj.timeFreqSlider = uicontrol('Style','slider','Parent',obj.fig,'Units','normalized',...
                        'Position',[0.275 0.05 0.4 0.03],'Min',50,'Max',5000,'Value',obj.f0,...
                        'Callback', @(src,~) obj.update_frequency(src.Value));

                    % -> readout to the right of slider
                    obj.sliderValueLabel = uicontrol('Style','text','Parent',obj.fig,'Units','normalized',...
                        'Position',[0.68 0.05 0.12 0.03], 'HorizontalAlignment','left', ...
                        'String','');  % filled by updateSliderLabel below

                    % ---- Horizontal visualization mode selector (top row) ----
                    obj.viewMode = 'real';
                    obj.btnGroup = uibuttongroup('Parent',obj.fig, ...
                        'Units','normalized','Position',[0.11 0.955 0.60 0.04], ...
                        'BorderType','none', ...
                        'SelectionChangedFcn', @(~,evt) obj.onViewChanged(evt));

                    uicontrol('Parent',obj.btnGroup,'Style','radiobutton', ...
                        'Units','normalized','Position',[0.00 -0.25 0.25 1], ...
                        'String','Re(P_{synth})','Tag','real');

                    uicontrol('Parent',obj.btnGroup,'Style','radiobutton', ...
                        'Units','normalized','Position',[0.25 -0.25 0.25 1], ...
                        'String','|P_{synth}|','Tag','mag');

                    uicontrol('Parent',obj.btnGroup,'Style','radiobutton', ...
                        'Units','normalized','Position',[0.5 -0.25 0.25 1], ...
                        'String','Amp. error','Tag','amperr');

                    % the virtual sources' own field (designed loudspeakers:
                    % their radiated field incl. near field and directivity)
                    uicontrol('Parent',obj.btnGroup,'Style','radiobutton', ...
                        'Units','normalized','Position',[0.75 -0.25 0.25 1], ...
                        'String','Re(P_{ref})','Tag','ref');

                    obj.btnGroup.SelectedObject = findobj(obj.btnGroup.Children,'Tag','real');
                case 'monitor'
            end

            % ---- Top-right checkbox (independent of mode) ----
            obj.rendererCheckbox = uicontrol('Style','checkbox','Parent',obj.fig,'Units','normalized', ...
                'Position',[0.75 0.955 0.23 0.04], 'String','Show renderer visuals', ...
                'Value',0, 'Callback', @(src,~) obj.onToggleRendererVisuals(src.Value));

            % ---- Play / pause animation (harmonic: phase, impulse: time) ----
            if any(strcmpi(obj.simulation_type, {'harmonic','impulse'}))
                obj.playButton = uicontrol('Style','togglebutton','Parent',obj.fig, ...
                    'Units','normalized','Position',[0.82 0.045 0.13 0.045], ...
                    'String',[char(9654) ' Play'], ...
                    'Callback', @(src,~) obj.onPlayToggled(src.Value));
            end

            % ---- Head diffraction checkbox (harmonic mode only) ----
            if strcmpi(obj.simulation_type, 'harmonic')
                obj.headCheckbox = uicontrol('Style','checkbox','Parent',obj.fig,'Units','normalized', ...
                    'Position',[0.50 0.955 0.23 0.04], 'String','Include head effect', ...
                    'Value',0, 'Callback', @(src,~) obj.onToggleHeadEffect(src.Value));
            end

            obj.set_input(obj.simulation_type);
            if strcmpi(obj.simulation_type,'impulse')
                obj.setTimeAxisFromOnset();
            end
            field = obj.get_field;
            obj.plot_field(field);
            obj.updateSliderLabel();

            if strcmpi(obj.simulation_type,'monitor')
                obj.startMonitorLoop();    % run at ~10 FPS
            end

        end

        function set_input(obj,type)
            input = repmat( [zeros(length(obj.temporary_sound_scene.virtual_sources{1}.source_signal.time_series),1)] ...
                ,1 , length(obj.temporary_sound_scene.virtual_sources));
            obj.temporary_sound_scene.render_sound_scene(input,false, false);

            fs = obj.temporary_sound_scene.virtual_sources{1}.source_signal.fs;
            N = length(obj.temporary_sound_scene.virtual_sources{1}.source_signal.time_series);
            obj.t = (1:N)/fs;
            switch type
                case 'harmonic'
                    for ti = 1 : 2
                        input = repmat( cos(2*pi*(obj.t + (ti-1)*N/fs)*obj.f0)' ...
                            ,1 , length(obj.temporary_sound_scene.virtual_sources));
                        obj.temporary_sound_scene.render_sound_scene(input,false,false);
                    end
                case 'harmonic_sin'
                    for ti = 1 : 2
                        input = repmat( sin(2*pi*(obj.t + (ti-1)*N/fs)*obj.f0)' ...
                            ,1 , length(obj.temporary_sound_scene.virtual_sources));
                        obj.temporary_sound_scene.render_sound_scene(input,false,false);
                    end
                case 'impulse'
                    input = repmat([hann(obj.pulse_width); zeros(N-obj.pulse_width,1)], ...
                        1, length(obj.temporary_sound_scene.virtual_sources));

                    % --- Master EQ alkalmazása, ha elérhető és engedélyezett ---
                    if ~isempty(obj.masterEQ) && isvalid(obj.masterEQ) && obj.masterEQ.Enabled
                        try
                            obj.masterEQ.setFs(fs);           % biztos, ami biztos: FS szinkron
                        catch
                        end
                        % MasterEQ.apply bemenet: [Nsamp x Nchan]
                        input = obj.masterEQ.apply(input);
                    end

                    obj.temporary_sound_scene.render_sound_scene(input,false,false);

                case 'monitor'
            end
        end


        function field = evaluate_field(obj, t_slider)
            c = audioapp.util.Physics.speedOfSound();
            % t_slider = (t_slider-0.5)/10e1;
            [X,Y] = meshgrid(obj.x,obj.y);
            field = zeros(size(X));
            switch obj.simulation_type
                case 'harmonic'
                    t0 = obj.t(end);
                    %title(obj.simulation_axes , sprintf('Reproduced field with harmonic excitation at f_0 = %d Hz',obj.f0));
                case 'impulse'
                    if isempty(obj.temporary_sound_scene.virtual_sources)
                        mean_d =  mean(sqrt(sum((bsxfun(@minus, obj.temporary_sound_scene.receiver.position,...
                            (cell2mat(cellfun( @(x) x.position', obj.temporary_sound_scene.loudspeaker_array, 'UniformOutput', false)))')).^2,2)));
                        t0 = mean_d/c  + t_slider;
                    else
                        mean_a =  mean(sqrt(sum((bsxfun(@minus, obj.temporary_sound_scene.receiver.position,...
                            (cell2mat(cellfun( @(x) x.position', obj.temporary_sound_scene.loudspeaker_array, 'UniformOutput', false)))')).^2,2)));
                        mean_b =  mean(sqrt(sum((bsxfun(@minus, obj.temporary_sound_scene.receiver.position,...
                            (cell2mat(cellfun( @(x) x.position', obj.temporary_sound_scene.virtual_sources, 'UniformOutput', false)))')).^2,2)));
                        t0 = t_slider;
                    end
                case 'monitor'
                    if isempty(obj.temporary_sound_scene.virtual_sources)
                        mean_d =  mean(sqrt(sum((bsxfun(@minus, obj.temporary_sound_scene.receiver.position,...
                            (cell2mat(cellfun( @(x) x.position', obj.temporary_sound_scene.loudspeaker_array, 'UniformOutput', false)))')).^2,2)));
                        t0 = mean_d/c;
                    else
                        mean_a =  mean(sqrt(sum((bsxfun(@minus, obj.temporary_sound_scene.receiver.position,...
                            (cell2mat(cellfun( @(x) x.position', obj.temporary_sound_scene.loudspeaker_array, 'UniformOutput', false)))')).^2,2)));
                        mean_b =  mean(sqrt(sum((bsxfun(@minus, obj.temporary_sound_scene.receiver.position,...
                            (cell2mat(cellfun( @(x) x.position', obj.temporary_sound_scene.virtual_sources, 'UniformOutput', false)))')).^2,2)));
                        t0 = 0.0213;
                    end

            end

            % Omnidirectional loudspeakers: one vectorised gather on the
            % cached geometry (field_engine) instead of an interp2 per speaker
            ls_all = obj.temporary_sound_scene.loudspeaker_array;
            if obj.allOmni(ls_all)
                [S, pos, ori] = obj.loudspeakerData(ls_all);
                if ~isempty(S)
                    obj.engine.prepare(pos, ori);
                    field = obj.engine.snapshot(obj.displayBandLimit(S), t0);
                end
                return;
            end

            for n = 1 : length(obj.temporary_sound_scene.loudspeaker_array)
                ls = obj.temporary_sound_scene.loudspeaker_array{n};
                ls_pos = ls.position;
                sig = ls.source_signal.time_series;
                sig = sig(:);  % force column
                Nsig = numel(sig);
                if Nsig == 0, continue; end

                theta_grid = (0:180)/180*pi;
                Dir_mx = ls.directivity.evaluate_freq_domain(Nsig, ls.source_signal.fs, theta_grid);
                % Dir_mx: [Nsig × 181]  (for omni: all ones → no change)

                in = fft(sig, Nsig);  % [Nsig × 1]
                Directive_signal = real(ifft(bsxfun(@times, Dir_mx, in), [], 1));
                % Directive_signal: [Nsig × 181]

                t_local = (1:Nsig) / ls.source_signal.fs;
                [Theta, T_grid] = meshgrid(theta_grid, t_local(:));

                xr = [X(:) - ls_pos(1), Y(:) - ls_pos(2)];
                R = reshape(sqrt(sum(xr.^2, 2)), size(X,1), size(X,2));
                v = bsxfun(@times, xr, 1./sqrt(sum(xr.^2, 2)));
                ori2 = ls.orientation(1:2);
                vn = v * ori2(:);
                Fi = reshape(real(acos(vn)), size(X,1), size(X,2));
                T_ret = t0 - R/c;
                field0 = 1./(4*pi*R) .* interp2(Theta, T_grid, Directive_signal, Fi, T_ret, 'linear');
                field0(isnan(field0)) = 0;

                field = field + field0;
            end
        end

        function field = get_field(obj)
            switch obj.simulation_type
                case 'harmonic'
                    % One render pass + closed-form evaluation: the steady
                    % state loudspeaker signals are s_n(t) = Re{A_n e^{jwt}},
                    % so the complex field Z gives Re (cos input), Im (sin
                    % input) and |P| at once.
                    obj.set_input('harmonic');
                    obj.Zharm = obj.harmonicField();
                    obj.animPhase = 0;
                    field = obj.harmonicView(obj.Zharm);

                case 'impulse'
                    obj.set_input('impulse');
                    field = obj.evaluate_field(obj.t0);
                case 'monitor'
                    field = obj.evaluate_field(0);
            end

        end
        function plot_field(obj,field)
            % explicit axes everywhere: this also runs from timers, where
            % gca/gcf may point at another figure
            ax = obj.simulation_axes;
            cmap = 'diverging';
            switch obj.simulation_type
                case 'harmonic'
                    switch obj.viewMode
                        case {'real', 'ref'}
                            set(obj.simulation_result,'CData',field)
                            if median(median(abs(field))) > 0
                                clim(ax, [-1,1]*5*median(median(abs(field))));
                            end
                        case 'mag'
                            set(obj.simulation_result,'CData',20*log10(field))
                            clim(ax, 20*log10([0.1,1]*5*median(median(abs(field)))))
                            cmap = 'magnitude';
                        case 'amperr'    % <-- was 'abserr' earlier
                            % If you want absolute error in dB, ensure 'field' is already a magnitude
                            set(obj.simulation_result,'CData',20*log10(abs(field)))
                            % Choose a sensible default range for error in dB:
                            clim(ax, [-80,-10])
                            cmap = 'magnitude';
                    end
                case 'impulse'
                    set(obj.simulation_result,'CData',field)
                    if isempty(obj.impulseClim)
                        obj.impulseClim = obj.listenerScale(field);
                    end
                    clim(ax, [-1,1]*obj.impulseClim)
                case 'monitor'
                    set(obj.simulation_result,'CData',field)
                    clim(ax, [-1,1]*1e-2)
            end
            colormap(ax, sound_scene_simulator.fieldColormap(cmap));
            set(obj.fig,'PaperPositionMode','auto');

            % Keep the colorbar + edit boxes in sync
            if ~isempty(obj.cbarHandle) && isvalid(obj.cbarHandle)
                obj.syncClimUI();
            end
        end

        function plot_transfer(obj)
            delete(obj.fig);
            N = length(obj.temporary_sound_scene.virtual_sources{1}.source_signal.time_series);
            transfer = zeros(N,1);
            fs = obj.temporary_sound_scene.virtual_sources{1}.source_signal.fs;

            for n = 1 : length(obj.temporary_sound_scene.loudspeaker_array)
                ls = obj.temporary_sound_scene.loudspeaker_array{n};
                s0 = ls.source_signal.get_spectrum;

                xr = obj.temporary_sound_scene.receiver.position - ls.position;
                R = sqrt(sum(xr.^2, 2));
                v = xr / R;
                ori = ls.orientation;
                D_dim = min(numel(v), numel(ori));
                vn = dot(v(1:D_dim), ori(1:D_dim));
                fi0 = real(acos(max(-1, min(1, vn))));

                % Get directivity at the receiver angle using the source's Directivity object
                freq_vec = (0:N-1)'/N*fs;
                D_at_angle = ls.directivity.evaluate(fi0, freq_vec');  % [1 × N]
                Dir0 = D_at_angle(:);  % [N × 1]

                c = audioapp.util.Physics.speedOfSound();
                transfer = transfer + 1/(4*pi) * s0 .* exp(-1i*2*pi*freq_vec/c*R) .* Dir0 / R;
            end
            types = cellfun(@(x) x.source_type.Shape,obj.temporary_sound_scene.virtual_sources,'UniformOutput',0);
            pw_source = cell2mat(cellfun(@(x) strcmp(x,'plane_wave'), types, 'UniformOutput',0));
            ps_source = cell2mat(cellfun(@(x) strcmp(x,'point_source'), types, 'UniformOutput',0));

            Nps = length( find(ps_source) );
            Npw = length( find(pw_source) );
            if Nps>0
                xs = cell2mat(cellfun(@(x) x.position,obj.temporary_sound_scene.virtual_sources(ps_source),'UniformOutput',0)') ;
                xr = obj.temporary_sound_scene.receiver.position;
                A0 = 1/4/pi/mean(sqrt(sum((xs-xr).^2,2)),1)*length( find(ps_source) ) + 1/4/pi*length( find(pw_source) ) ;
            else
                A0 = Npw/4/pi;
            end
            fig = figure('units','normalized','outerposition',[0.2 0.3 0.35 0.4]);
            fs = obj.temporary_sound_scene.virtual_sources{1}.source_signal.fs;
            freq = (0:length(transfer)-1)/N*fs;
            semilogx(freq,20*log10(abs(transfer) / A0),'LineWidth',1.5);
            grid on
            xlim([10 ,20e3]);
            ylim([-20,10]);
            xlabel('f [Hz]')
            ylabel('P_{rec} [dB]')
            title(sprintf('Synthesized field transfer at x = %.1f, y = %.1f',...
                obj.temporary_sound_scene.receiver.position(1),...
                obj.temporary_sound_scene.receiver.position(2)));
            p = [];
            for ri = 1 : length(obj.temporary_sound_scene.scene_renderer.SFS_renderer)
                if isa((obj.temporary_sound_scene.scene_renderer.SFS_renderer{ri}),'wfs_renderer')
                    f(ri,:) = obj.temporary_sound_scene.scene_renderer.SFS_renderer{ri}.get_lower_cutoff_frequency(obj.temporary_sound_scene.receiver);
                    yl = get(gca,'yLim');
                    hold on
                    p(ri) = plot([1,1]*max(f(ri,:)),yl,'--k','LineWidth',1.25,'Color',[0 0 0]+(0.15*ri));
                    tag{ri} = sprintf('f_c^{%s}',obj.gui.virtual_source_points{ri}.UserData.text.String)
                    %%
                end
            end
            if ~isempty(p)
                legend(p, tag);
            end




        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%% CALLBACKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function update_frequency(obj, freq)
            obj.f0 = freq;
            field = obj.get_field;      % renders the new excitation itself
            obj.plot_field(field);
            obj.updateSliderLabel();
        end

        function update_time(obj, t0)
            obj.t0 = t0;                  % <- store t0 so readout reflects it
            field = obj.evaluate_field(t0);
            obj.plot_field(field);
            obj.updateSliderLabel();
        end

        function onViewChanged(obj, evt)
            obj.viewMode = get(evt.NewValue,'Tag');   % 'real' | 'mag' | 'amperr'
            if strcmpi(obj.simulation_type,'harmonic') && ~isempty(obj.Zharm)
                field = obj.harmonicView(obj.Zharm);  % same f0: no re-render
            else
                field = obj.get_field;
            end
            obj.plot_field(field);
        end

        function field = harmonicView(obj, Z)
            % Real part / magnitude / amplitude error of the complex field
            switch obj.viewMode
                case 'real'
                    field = real(Z);
                case 'mag'
                    field = abs(Z);
                case 'amperr'
                    field = abs(Z) - abs(obj.referenceField(obj.t0));
                case 'ref'
                    field = real(obj.referenceField(obj.t(end)));
                    return;                 % no head scattering on the target field
            end
            % ── Add head diffraction if enabled ──
            if obj.includeHeadEffect
                field = obj.addHeadScattering(field);
            end
        end

        function P = referenceField(obj, t0)
            % sum of the virtual sources' own fields at f0 (the target)
            [X,Y] = meshgrid(obj.x,obj.y);
            P = zeros(size(X));
            for n = 1 : length(obj.temporary_sound_scene.virtual_sources)
                P = P + obj.temporary_sound_scene.virtual_sources{n}.get_sound_field(X,Y,obj.f0,t0);
            end
        end

        function Z = harmonicField(obj)
            % Complex harmonic field from the rendered steady-state block:
            % least-squares fit s_n(t) = a_n cos(wt) + b_n sin(wt) on the
            % block's own time axis gives A_n = a_n - j b_n, and
            % Re{Z} reproduces the time-domain snapshot at t = obj.t(end).
            ls_all = obj.temporary_sound_scene.loudspeaker_array;
            [S, pos, ori] = obj.loudspeakerData(ls_all);
            fs = obj.temporary_sound_scene.virtual_sources{1}.source_signal.fs;
            w = 2*pi*obj.f0;
            t = (1:size(S,1))'/fs;
            ab = [cos(w*t), sin(w*t)] \ S;
            A = (ab(1,:) - 1i*ab(2,:)).';
            dirs = {};
            if ~obj.allOmni(ls_all)
                dirs = cellfun(@(l) l.directivity, ls_all, 'UniformOutput', false);
            end
            obj.engine.prepare(pos, ori);
            Z = exp(1i*w*obj.t(end)) * obj.engine.harmonic(A, obj.f0, dirs);
        end

        % ── Animation (Play button) ──
        function onPlayToggled(obj, isOn)
            if ~isOn, obj.stopAnimation(); return; end
            if strcmpi(obj.simulation_type,'harmonic') && ~strcmpi(obj.viewMode,'real')
                % only the instantaneous field moves; switch to Re{P}
                obj.viewMode = 'real';
                if ~isempty(obj.btnGroup) && isvalid(obj.btnGroup)
                    obj.btnGroup.SelectedObject = findobj(obj.btnGroup.Children,'Tag','real');
                end
                obj.plot_field(obj.harmonicView(obj.Zharm));
            end
            if strcmpi(obj.simulation_type,'impulse') && obj.t0 >= obj.timeFreqSlider.Max
                obj.t0 = obj.timeFreqSlider.Min;      % replay from the first onset
            end
            set(obj.playButton, 'String', [char(10074) char(10074) ' Pause']);
            obj.stopAnimationTimer();
            obj.animTimer = timer('ExecutionMode','fixedRate', ...
                'Period', round(1000/obj.animFps)/1000, 'BusyMode','drop', ...
                'TimerFcn', @(~,~) obj.onAnimTick(), 'Name','SoundFieldAnimation', ...
                'UserData', struct('fig', obj.fig, 'scene_renderer', []));  % found by closeFigure
            start(obj.animTimer);
        end

        function stopAnimation(obj)
            obj.stopAnimationTimer();
            obj.animSignals = [];
            if ~isempty(obj.playButton) && isvalid(obj.playButton)
                set(obj.playButton, 'Value', 0, 'String', [char(9654) ' Play']);
            end
            % leave the widgets consistent with the last frame
            if strcmpi(obj.simulation_type,'impulse') && ~isempty(obj.timeFreqSlider) ...
                    && isvalid(obj.timeFreqSlider)
                obj.timeFreqSlider.Value = obj.t0;
                obj.updateSliderLabel();
            end
        end

        function stopAnimationTimer(obj)
            if ~isempty(obj.animTimer) && isvalid(obj.animTimer)
                try, stop(obj.animTimer); catch, end
                try, delete(obj.animTimer); catch, end
            end
            obj.animTimer = [];
        end

        function onAnimTick(obj)
            if isempty(obj.fig) || ~ishandle(obj.fig)
                obj.stopAnimationTimer();
                return;
            end
            try
                switch obj.simulation_type
                    case 'harmonic'
                        % one period in 1.5 s
                        obj.animPhase = mod(obj.animPhase + 2*pi/(1.5*obj.animFps), 2*pi);
                        % free field only: the head scattering is a full
                        % series evaluation, too slow to redo every frame
                        set(obj.simulation_result, 'CData', single(real(obj.Zharm*exp(1i*obj.animPhase))));
                    case 'impulse'
                        % 1 m of propagation per ~2 s, looping over the slider range
                        s = obj.timeFreqSlider;
                        t = obj.t0 + 1/(2*obj.animFps)/audioapp.util.Physics.speedOfSound();
                        if t > s.Max, t = s.Min; end
                        obj.t0 = t;
                        ls_all = obj.temporary_sound_scene.loudspeaker_array;
                        if obj.allOmni(ls_all)
                            % the rendered signals are fixed while playing:
                            % collect them once per run
                            if isempty(obj.animSignals)
                                [S, pos, ori] = obj.loudspeakerData(ls_all);
                                obj.animSignals = obj.displayBandLimit(S);
                                obj.engine.prepare(pos, ori);
                            end
                            F = obj.engine.snapshot(obj.animSignals, t);
                        else
                            F = obj.evaluate_field(t);
                        end
                        set(obj.simulation_result, 'CData', single(F));
                        % UI widgets cost more than the field: refresh them sparsely
                        obj.animFrame = obj.animFrame + 1;
                        if mod(obj.animFrame, 5) == 0
                            s.Value = t;
                            obj.updateSliderLabel();
                        end
                end
                drawnow limitrate
            catch
                obj.stopAnimation();
            end
        end

        function setTimeAxisFromOnset(obj)
            % Time axis of the impulse view: 0 = the moment the first
            % loudspeaker actually sounds (the rendering latency differs
            % from method to method), end = the last loudspeaker's signal
            % has left the displayed area. Opens on the main wavefront at
            % the listener.
            % "Sounds" = within 40 dB of the loudest sample: lower
            % thresholds catch the long pre-ringing of the sqrt(jw)
            % prefilter (-60 dB starts ~9 ms before the actual front).
            ls_all = obj.temporary_sound_scene.loudspeaker_array;
            [S, pos] = obj.loudspeakerData(ls_all);
            if isempty(S), return; end
            fs = obj.temporary_sound_scene.virtual_sources{1}.source_signal.fs;
            c = audioapp.util.Physics.speedOfSound();
            env = max(abs(S), [], 2);
            on = find(env > 1e-2*max(env), 1, 'first');
            off = find(env > 1e-2*max(env), 1, 'last');
            if isempty(on), return; end
            [Xc, Yc] = meshgrid(obj.x([1 end]), obj.y([1 end]));
            far = max(sqrt((Xc(:)' - pos(:,1)).^2 + (Yc(:)' - pos(:,2)).^2), [], 'all');
            obj.tOnset = on/fs;
            tEnd = off/fs + far/c;
            % main wavefront at the listener: peak of the displayed
            % (band-limited) pressure there
            rec = obj.temporary_sound_scene.receiver.position(1:2);
            R = sqrt(sum((pos - rec).^2, 2));
            Sd = obj.displayBandLimit(S);
            t = obj.tOnset : 1/fs : tEnd;
            p = zeros(size(t));
            for n = 1:size(Sd,2)
                p = p + interp1((1:size(Sd,1))/fs, Sd(:,n), t - R(n)/c, 'linear', 0) / (4*pi*R(n));
            end
            [~, ip] = max(abs(p));
            obj.t0 = t(ip);
            obj.impulseClim = [];                             % rescale on the next plot
            s = obj.timeFreqSlider;
            if ~isempty(s) && isvalid(s)
                span = tEnd - obj.tOnset;
                set(s, 'Min', obj.tOnset, 'Max', tEnd, 'Value', obj.t0, ...
                    'SliderStep', min(1, [1e-4, 2e-3]/span));
            end
        end

        function a = listenerScale(obj, field)
            % Colour scale of the impulse view: 1.2 x the peak within 0.5 m
            % of the listener (the frame opens with the wavefront there).
            % A global max would be set by the 1/R peaks right at the
            % loudspeakers and leave the reproduced wavefront invisible.
            [X, Y] = meshgrid(obj.x, obj.y);
            rec = obj.temporary_sound_scene.receiver.position(1:2);
            near = (X - rec(1)).^2 + (Y - rec(2)).^2 < 0.5^2;
            a = 1.2*max(abs(field(near)), [], 'all');
            if isempty(a) || ~(a > 0)
                a = prctile(abs(field(:)), 99.5);
            end
            if ~(a > 0), a = 1e-2; end
        end

        function S = displayBandLimit(obj, S)
            % Zero-phase low-pass of the loudspeaker signals to what the
            % grid can show (spatial Nyquist c/(2 dx)): pulses thinner than
            % a pixel would otherwise alias into speckle. Flat up to 70 %
            % of that frequency, raised-cosine roll-off to zero at it.
            fs = obj.temporary_sound_scene.virtual_sources{1}.source_signal.fs;
            fN = obj.displayBandwidth();
            if fN >= fs/2, return; end
            N = size(S,1);
            f = abs([0:ceil(N/2)-1, -floor(N/2):-1]') * fs / N;
            H = 0.5*(1 + cos(pi*min(max((f - 0.7*fN)/(0.3*fN), 0), 1)));
            S = real(ifft(fft(S) .* H));
        end

        function fN = displayBandwidth(obj)
            fN = audioapp.util.Physics.speedOfSound() / (2*obj.dx_sim);
        end


        function syncClimUI(obj)
            % Reflect current caxis into the edit boxes
            cl = caxis(obj.simulation_axes);
            set(obj.cminEdit,'String',sprintf('%g',cl(1)));
            set(obj.cmaxEdit,'String',sprintf('%g',cl(2)));

            % Optional: label units depending on view
            switch obj.simulation_type
                case 'harmonic'
                    switch obj.viewMode
                        case 'real'
                            obj.cbarHandle.Label.String = 'Re\{P_{synth}\}';
                        case 'ref'
                            obj.cbarHandle.Label.String = 'Re\{P_{ref}\}';
                        case 'mag'
                            obj.cbarHandle.Label.String = '|P_{synth}| (dB)';
                        case 'amperr'
                            obj.cbarHandle.Label.String = 'Amplitude error (dB)'; % or whatever you prefer
                    end
                case 'impulse'
                    obj.cbarHandle.Label.String = sprintf('p(t), shown up to %.1f kHz', ...
                        obj.displayBandwidth()/1e3);
            end
        end

        function onClimEdited(obj, src, whichSide)
            % Read both edits; replace the edited side with the user's value
            cl = caxis(obj.simulation_axes);
            try
                val = str2double(get(src,'String'));
                if isnan(val) || ~isfinite(val)
                    % reset text if invalid
                    obj.syncClimUI();
                    return;
                end
                if strcmp(whichSide,'min')
                    cl(1) = val;
                else
                    cl(2) = val;
                end
                % Ensure proper ordering
                if cl(1) == cl(2)
                    % Avoid singular caxis; nudge upper
                    cl(2) = cl(2) + eps;
                elseif cl(1) > cl(2)
                    % Swap if out of order
                    cl = fliplr(cl);
                end
                caxis(obj.simulation_axes, cl);
                obj.syncClimUI();
            catch
                obj.syncClimUI();
            end
        end

        function updateSliderLabel(obj)
            switch obj.simulation_type
                case 'harmonic'
                    % show kHz with 3 sig figs
                    txt = sprintf('%0.3f kHz', obj.f0/1000);
                case 'impulse'
                    % time since the first loudspeaker starts to sound
                    txt = sprintf('%0.1f ms', (obj.t0 - obj.tOnset)*1e3);
            end
            if ~isempty(obj.sliderValueLabel) && isvalid(obj.sliderValueLabel)
                set(obj.sliderValueLabel,'String',txt);
            end
        end
        function onToggleRendererVisuals(obj, isOn)
            if isOn
                % If we haven’t copied yet (or copies were deleted), copy now
                if isempty(obj.rendererCopies) || any(~isvalid(obj.rendererCopies))
                    obj.copyRendererVisuals();
                end
                obj.setRendererCopiesVisible('on');
            else
                delete(obj.rendererCopies);
                obj.rendererCopies = [];
            end
        end
        function refreshRendererCopies(obj)
            % If checkbox is on, recopy from gui.renderer_illustration
            if ~isempty(obj.rendererCheckbox) && isvalid(obj.rendererCheckbox) ...
                    && get(obj.rendererCheckbox,'Value') == 1
                % Delete old copies
                if ~isempty(obj.rendererCopies)
                    delete(obj.rendererCopies(ishghandle(obj.rendererCopies)));
                end
                obj.rendererCopies = [];
                % Recopy fresh
                obj.copyRendererVisuals();
            end
        end

        function copyRendererVisuals(obj)
            src = obj.gui.renderer_illustration;
            if isempty(src), return; end
            flds = fieldnames(src);
            copies = gobjects(0);
            for k = 1:numel(flds)
                h = src.(flds{k});
                if ishghandle(h)
                    c = copyobj(h, obj.simulation_axes);
                    set(c, 'HitTest','off','PickableParts','none', ...
                        'Visible', get(h,'Visible'));
                    c.UserData = struct('src', h, 'field', flds{k}); % <-- store link
                    copies(end+1,1) = c; %#ok<AGROW>
                end
            end
            obj.rendererCopies = copies;
        end


        function setRendererCopiesVisible(obj, mode)
            for n = 1:numel(obj.rendererCopies)
                c = obj.rendererCopies(n);
                if ~isvalid(c), continue; end
                ud = get(c,'UserData');
                srcVis = 'off';
                if isstruct(ud) && isfield(ud,'src') && ishghandle(ud.src)
                    srcVis = get(ud.src,'Visible');
                end
                switch mode
                    case 'off'
                        set(c,'Visible','off');
                    case {'on','sync'}
                        set(c,'Visible', srcVis);  % mirror source
                end
            end
        end

        function startMonitorLoop(obj)
            % Stop any existing timer first
            obj.stopMonitorLoop();
            % the monitor reads the live loudspeaker signals, which the
            % fast binaural path does not compute
            obj.setKeepLoudspeakerSignals(true);

            obj.monitorTimer = timer( ...
                'ExecutionMode','fixedSpacing', ...
                'Period', 1/obj.monitorRate, ...
                'TimerFcn', @(~,~) obj.onMonitorTick(), ...
                'Name','SoundFieldMonitorTimer', ...
                'UserData', struct('fig', obj.fig, ...             % found by closeFigure
                    'scene_renderer', obj.temporary_sound_scene.scene_renderer));

            try
                start(obj.monitorTimer);
            catch
                % If timers are disabled in this context, just ignore
            end
        end

        function stopMonitorLoop(obj)
            if ~isempty(obj.monitorTimer) && isvalid(obj.monitorTimer)
                try, stop(obj.monitorTimer); catch, end
                try, delete(obj.monitorTimer); catch, end
                obj.setKeepLoudspeakerSignals(false);
            end
            obj.monitorTimer = [];
        end

        function setKeepLoudspeakerSignals(obj, tf)
            try
                obj.temporary_sound_scene.scene_renderer.KeepLoudspeakerSignals = tf;
            catch
            end
        end

        function onMonitorTick(obj)
            % Stop if the figure/axes are gone
            if isempty(obj.fig) || ~ishandle(obj.fig) || ...
                    isempty(obj.simulation_axes) || ~ishandle(obj.simulation_axes)
                obj.stopMonitorLoop();
                return;
            end

            % Compute + plot one frame
            try
                fld = obj.evaluate_field(0);       % use current signals
                obj.plot_field(fld);
                drawnow limitrate nocallbacks      % yield without flooding
            catch
                % On errors during live updates, stop gracefully
                obj.stopMonitorLoop();
            end
        end

        function delete(obj)
            % Make sure timer is cleaned up if user closes the window
            obj.stopAnimationTimer();
            obj.stopMonitorLoop();
            if strcmpi(obj.simulation_type, 'monitor')
                % closeFigure removed the timer already: re-enable the
                % fast binaural path here
                obj.setKeepLoudspeakerSignals(false);
            end
        end

        % ── Head diffraction toggle ──
        function onToggleHeadEffect(obj, isOn)
            obj.includeHeadEffect = logical(isOn);

            % Lazy-init the Diffractor with receiver position and radius
            if isOn && isempty(obj.diffractor)
                obj.diffractor = Diffractor();
                obj.diffractor.scatterer_pos = obj.temporary_sound_scene.receiver.position;
                obj.diffractor.scatterer_radius = 0.0875;  % ~8.75 cm head radius
            end

            % Refresh the field
            field = obj.get_field();
            obj.plot_field(field);
        end

        % ── Compute and superpose scattered field from all loudspeakers ──
        function field_out = addHeadScattering(obj, field_free)
            % field_free: [Ny x Nx] complex free-field pressure at f0
            % field_out:  same size, with scattering added

            if isempty(obj.diffractor)
                obj.diffractor = Diffractor();
                obj.diffractor.scatterer_pos = obj.temporary_sound_scene.receiver.position;
                obj.diffractor.scatterer_radius = 0.0875;
            end

            d = obj.diffractor;
            d.scatterer_pos = obj.temporary_sound_scene.receiver.position;

            [X, Y] = meshgrid(obj.x, obj.y);
            x_eval = [X(:), Y(:)];    % [Npts x 2]
            Npts = size(x_eval, 1);

            % For each loudspeaker, compute (total - incident) = scattered field
            % and add it to the free-field result
            p_scattered_total = zeros(Npts, 1);

            ls_array = obj.temporary_sound_scene.loudspeaker_array;
            for n = 1:numel(ls_array)
                src_pos = ls_array{n}.position;

                % Get the driving amplitude for this loudspeaker at f0
                % (already included in field_free via evaluate_field,
                %  so we need the complex amplitude to weight the scattering)
                sig = ls_array{n}.source_signal.time_series;
                if isempty(sig), continue; end

                % Estimate complex amplitude at f0 from the signal
                % (the signal was set by set_input('harmonic') before evaluate_field)
                Nfft = numel(sig);
                fs = ls_array{n}.source_signal.fs;
                S = fft(sig);
                f_bin = round(obj.f0 / fs * Nfft) + 1;
                f_bin = max(1, min(f_bin, Nfft));
                A_complex = S(f_bin) * 2 / Nfft;   % complex amplitude at f0

                if abs(A_complex) < 1e-15, continue; end

                % Set diffractor source
                d.source_pos  = src_pos;
                d.source_type = 'point';

                % Total field with scattering
                p_total = d.evaluate(obj.f0, x_eval);       % [Npts x 1]

                % Free-field from this source alone
                p_ff = d.evaluate_freefield(obj.f0, x_eval); % [Npts x 1]

                % Scattered = total - incident
                p_scat = p_total - p_ff;

                % Weight by the loudspeaker's complex driving amplitude
                p_scattered_total = p_scattered_total + A_complex * p_scat;
            end

            % Add scattered field to the full free-field result
            field_out = field_free + reshape(real(p_scattered_total), size(field_free));
        end
    end

    methods (Static)
        function closeFigure(fig, ~)
            % Stop + delete this window's timers (they hold the simulator),
            % then the window; no simulator method is on the stack while
            % the object is released.
            ts = timerfindall;
            for k = 1:numel(ts)
                ud = ts(k).UserData;
                if isstruct(ud) && isfield(ud, 'fig') && isequal(ud.fig, fig)
                    try, stop(ts(k)); catch, end
                    try, delete(ts(k)); catch, end
                    % a closed monitor no longer needs the loudspeaker
                    % signals: give the fast binaural path back
                    if ~isempty(ud.scene_renderer) && isvalid(ud.scene_renderer)
                        ud.scene_renderer.KeepLoudspeakerSignals = false;
                    end
                end
            end
            delete(fig);
        end

        function tf = allOmni(ls_all)
            tf = all(cellfun(@(l) isa(l.directivity, 'OmnidirectionalDirectivity'), ls_all));
        end

        function [S, pos, ori] = loudspeakerData(ls_all)
            % S(:,n): current loudspeaker signals; pos/ori: [Nls x 2]
            sig = cellfun(@(l) l.source_signal.time_series(:), ls_all, 'UniformOutput', false);
            if any(cellfun(@isempty, sig)), S = []; else, S = [sig{:}]; end
            pos = cell2mat(cellfun(@(l) l.position(1:2), ls_all(:), 'UniformOutput', false));
            ori = cell2mat(cellfun(@(l) l.orientation(1:2), ls_all(:), 'UniformOutput', false));
        end

        function cmap = fieldColormap(kind)
            % 'diverging': blue-white-red for signed pressure (zero = white)
            % 'magnitude': perceptually uniform sequential map for dB views
            switch kind
                case 'diverging'
                    k = [ 5  48  97;  33 102 172;  67 147 195; 146 197 222; 209 229 240; ...
                         247 247 247; 253 219 199; 244 165 130; 214  96  77; 178  24  43; 103   0  31]/255;
                    cmap = interp1(linspace(0,1,size(k,1)), k, linspace(0,1,256), 'pchip');
                    cmap = min(max(cmap, 0), 1);
                otherwise
                    cmap = parula(256);
            end
        end
    end
end

