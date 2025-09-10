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
        viewMode            % 'real' | 'mag' | 'amperr'
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
            dx = 1e-2;
            obj.x = [obj.gui.main_axes.XLim(1):dx:obj.gui.main_axes.XLim(2)]';
            obj.y = [obj.gui.main_axes.YLim(1):dx:obj.gui.main_axes.YLim(2)]';
            [X,Y] = meshgrid(obj.x,obj.y);
            xr = [0,0];
            obj.fig = figure('units','normalized','position',[0.15 0.2 0.45 0.6]);

            set(obj.fig,'defaulttextinterpreter','latex')

            obj.simulation_axes = axes('Units','normalized','Position',[0.0,0.2,0.7,0.7]);
            obj.simulation_result  = pcolor(obj.x,obj.y,zeros(size(X)));
            %   obj.simulation_axes.Position  =  [0.2,0.12,0.7,0.7]; % case 1
            % obj.simulation_axes.Position = [0.14,0.12,0.85,0.875]; %case 2
            shading interp
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
                        'Position',[0.275 0.05 0.4 0.03],'Min',-0.5e-2+t0_default,'Max',2e-2+t0_default,'Value',t0_default,...
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
                        'Units','normalized','Position',[0.00 -0.25 0.34 1], ...
                        'String','Re(P_{synth})','Tag','real');

                    uicontrol('Parent',obj.btnGroup,'Style','radiobutton', ...
                        'Units','normalized','Position',[0.25 -0.25 0.33 1], ...
                        'String','|P_{synth}|','Tag','mag');

                    uicontrol('Parent',obj.btnGroup,'Style','radiobutton', ...
                        'Units','normalized','Position',[0.5 -0.25 0.33 1], ...
                        'String','Amp. error','Tag','amperr');

                    obj.btnGroup.SelectedObject = findobj(obj.btnGroup.Children,'Tag','real');
                case 'monitor'
            end

            % ---- Top-right checkbox (independent of mode) ----
            obj.rendererCheckbox = uicontrol('Style','checkbox','Parent',obj.fig,'Units','normalized', ...
                'Position',[0.75 0.955 0.23 0.04], 'String','Show renderer visuals', ...
                'Value',0, 'Callback', @(src,~) obj.onToggleRendererVisuals(src.Value));

            obj.set_input(obj.simulation_type);
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
                    input = repmat( [hann(obj.pulse_width);zeros(N-obj.pulse_width,1)] ...
                        ,1 , length(obj.temporary_sound_scene.virtual_sources));
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

            for n = 1 : length(obj.temporary_sound_scene.loudspeaker_array)
                switch obj.temporary_sound_scene.scene_renderer.binaural_renderer{n}.binaural_source.source_type.Shape
                    case 'point_source'
                        R = sqrt( (X-obj.temporary_sound_scene.scene_renderer.binaural_renderer{n}.binaural_source.position(1)).^2 +...
                            (Y-obj.temporary_sound_scene.scene_renderer.binaural_renderer{n}.binaural_source.position(2)).^2 );
                        T = t0 - R/c;
                        field0 = 1./(4*pi*R).*reshape(...
                            interp1(obj.t, obj.temporary_sound_scene.scene_renderer.binaural_renderer{n}.binaural_source.source_signal.time_series, T(:),'linear'),...
                            size(X,1),size(X,2));
                        field0(isnan(field0)) = 0;
                    otherwise
                        idx = get_dirtable_idx(obj.temporary_sound_scene.scene_renderer.directivity_tables,...
                            obj.temporary_sound_scene.loudspeaker_array{n});
                        Dir0 = obj.temporary_sound_scene.scene_renderer.directivity_tables{idx}.directivity_mx;
                        in =  fft(obj.temporary_sound_scene.loudspeaker_array{n}.source_signal.time_series,size(Dir0,1));
                        Directive_signal = ifft(bsxfun( @times, Dir0, in ),[],1);
                        Directive_signal = Directive_signal(1:end/2,:);
                        theta = obj.temporary_sound_scene.scene_renderer.directivity_tables{idx}.theta;
                        [Theta, T] = meshgrid(theta,obj.t');

                        xr = [X(:)-obj.temporary_sound_scene.loudspeaker_array{n}.position(1),...
                            Y(:)-obj.temporary_sound_scene.loudspeaker_array{n}.position(2)];
                        R = reshape(sqrt( sum(xr.^2,2) ) , size(X,1),size(X,2));
                        v = bsxfun(@times,xr,1./(sqrt( sum(xr.^2,2) )));
                        vn = v*obj.temporary_sound_scene.loudspeaker_array{n}.orientation';
                        Fi = reshape(acos(vn), size(X,1),size(X,2));
                        Fi = real(Fi);
                        T_ret = t0 - R/c;
                        field0 = 1./(4*pi*R).*interp2(Theta,T,Directive_signal, Fi, T_ret, 'linear');
                        field0(isnan(field0)) = 0;
                end
                field = field + field0;
            end
        end

        function field = get_field(obj)
            switch obj.simulation_type
                case 'harmonic'
                    switch obj.viewMode
                        case 'real'
                            obj.set_input('harmonic');
                            ReP = obj.evaluate_field(obj.t0);
                            field = ReP;
                        case 'mag'
                            obj.set_input('harmonic');
                            ReP = obj.evaluate_field(obj.t0);
                            obj.set_input('harmonic_sin');
                            ImP = obj.evaluate_field(obj.t0);
                            field = abs(ReP + 1i*ImP);
                        case 'amperr'
                            obj.set_input('harmonic');
                            ReP = obj.evaluate_field(obj.t0);
                            obj.set_input('harmonic_sin');
                            ImP = obj.evaluate_field(obj.t0);
                            field = ReP + 1i*ImP;

                            % Calculate reference field
                            [X,Y] = meshgrid(obj.x,obj.y);
                            field_ref = zeros(size(X));

                            for n = 1 : length(obj.temporary_sound_scene.virtual_sources)
                                field_ref = field_ref + obj.temporary_sound_scene.virtual_sources{n}.get_sound_field(X,Y,obj.f0,obj.t0);
                            end
                            field = abs(field) - abs(field_ref);

                    end
                case 'impulse'
                    obj.set_input('impulse');
                    field = obj.evaluate_field(obj.t0);
                case 'monitor'
                    field = obj.evaluate_field(0);
            end

        end
        function plot_field(obj,field)
            switch obj.simulation_type
                case 'harmonic'
                    switch obj.viewMode
                        case 'real'
                            set(obj.simulation_result,'CData',field)
                            caxis([-1,1]*5*median(median(abs(field))))
                        case 'mag'
                            set(obj.simulation_result,'CData',20*log10(field))
                            caxis(20*log10([0.1,1]*5*median(median(abs(field)))))
                        case 'amperr'    % <-- was 'abserr' earlier
                            % If you want absolute error in dB, ensure 'field' is already a magnitude
                            set(obj.simulation_result,'CData',20*log10(abs(field)))
                            % Choose a sensible default range for error in dB:
                            caxis([-80,-10])
                    end
                case 'impulse'
                    set(obj.simulation_result,'CData',field)
                    caxis([-1,1]*1e-2)
                case 'monitor'
                    set(obj.simulation_result,'CData',field)
                    caxis([-1,1]*1e-2)
            end
            set(gcf,'PaperPositionMode','auto');

            % Keep the colorbar + edit boxes in sync
            if ~isempty(obj.cbarHandle) && isvalid(obj.cbarHandle)
                obj.syncClimUI();
            end
        end

        function plot_transfer(obj)
            delete(obj.fig);
            N = length(obj.temporary_sound_scene.virtual_sources{1}.source_signal.time_series);
            transfer = zeros(N,1);
            for n = 1 : length(obj.temporary_sound_scene.loudspeaker_array)
                s0 = obj.temporary_sound_scene.loudspeaker_array{n}.source_signal.get_spectrum;

                idx = get_dirtable_idx(obj.temporary_sound_scene.scene_renderer.directivity_tables,...
                    obj.temporary_sound_scene.loudspeaker_array{n});
                Dir_mx = obj.temporary_sound_scene.scene_renderer.directivity_tables{idx}.directivity_mx;

                xr = [obj.temporary_sound_scene.receiver.position(1)-obj.temporary_sound_scene.loudspeaker_array{n}.position(1),...
                    obj.temporary_sound_scene.receiver.position(2)-obj.temporary_sound_scene.loudspeaker_array{n}.position(2)];
                R = sqrt( sum(xr.^2,2) );
                v = bsxfun(@times,xr,1./(sqrt( sum(xr.^2,2) )));
                vn = v*obj.temporary_sound_scene.loudspeaker_array{n}.orientation';
                fi0 = real(acos(vn));

                theta = obj.temporary_sound_scene.scene_renderer.directivity_tables{idx}.theta;
                Dir0 = fft(ifft(interp1(theta, Dir_mx.', fi0, 'linear'),length(s0)));

                fs = obj.temporary_sound_scene.virtual_sources{1}.source_signal.fs;
                freq_vec = (0:N-1)'/N*fs;
                c = audioapp.util.Physics.speedOfSound();
                transfer = transfer + 1/(4*pi)*s0.*exp(-1i*2*pi*freq_vec/c*R).*Dir0.'/R;
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
            obj.set_input(obj.simulation_type);
            field = obj.get_field;
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
            field = obj.get_field;
            obj.plot_field(field);
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
                        case 'mag'
                            obj.cbarHandle.Label.String = '|P_{synth}| (dB)';
                        case 'amperr'
                            obj.cbarHandle.Label.String = 'Amplitude error (dB)'; % or whatever you prefer
                    end
                case 'impulse'
                    obj.cbarHandle.Label.String = 'p(t)';
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
                    fs = obj.temporary_sound_scene.virtual_sources{1}.source_signal.fs;
                    N = length(obj.temporary_sound_scene.virtual_sources{1}.source_signal.time_series);


                    % the slider value is in seconds (per your evaluate_field),
                    % show milliseconds (relative t0 selection)
                    txt = sprintf('%0.1f ms', ( obj.t0 - N/2/fs)*1e3);
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

            obj.monitorTimer = timer( ...
                'ExecutionMode','fixedSpacing', ...
                'Period', 1/obj.monitorRate, ...
                'TimerFcn', @(~,~) obj.onMonitorTick(), ...
                'Name','SoundFieldMonitorTimer');

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
            end
            obj.monitorTimer = [];
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
            obj.stopMonitorLoop();
        end
    end
end

