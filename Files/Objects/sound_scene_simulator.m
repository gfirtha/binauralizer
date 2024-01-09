classdef sound_scene_simulator < handle
    %SCENE_SIMULATOR Summary of this class goes here
    %   Detailed explanation goes here

    properties
        temporary_sound_scene
        excitation_type
        f0
        t
        x
        y
        fig
        simulation_axes
        simulation_result
    end

    methods
        function obj = sound_scene_simulator(varargin)
            %SCENE_SIMULATOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.temporary_sound_scene = varargin{1}.duplicate_scene;
            gui = varargin{2};
            obj.excitation_type = varargin{3};
            dx = 1e-2;
            obj.x = [gui.main_axes.XLim(1):dx:gui.main_axes.XLim(2)]';
            obj.y = [gui.main_axes.YLim(1):dx:gui.main_axes.YLim(2)]';
            [X,Y] = meshgrid(obj.x,obj.y);
            xr = [0,0];
            obj.fig = figure('units','normalized','outerposition',[0.2 0.3 0.35 0.7]);

            set(obj.fig,'defaulttextinterpreter','latex')

            obj.simulation_axes = axes;
            obj.simulation_result  = pcolor(obj.x,obj.y,zeros(size(X)));
            obj.simulation_axes.Position  =  [0.11,0.12,0.85,0.875]; % case 1
           % obj.simulation_axes.Position = [0.14,0.12,0.85,0.875]; %case 2
            shading interp
            axis equal tight
            ftsize = 12;  %case1
            %ftsize = 16; %case2
            copyobj(gui.receiver, obj.simulation_axes);
            for n = 1 : length(gui.virtual_source_points)
                if isvalid(gui.virtual_source_points{n})
                copyobj(gui.virtual_source_points{n}, obj.simulation_axes);
                end
            end
            for n = 1 : length(gui.loudspeaker_points)
                copyobj(gui.loudspeaker_points{n}, obj.simulation_axes);
            end
            xlim(gui.main_axes.XLim)
            ylim(gui.main_axes.YLim)
            caxis([-1,1])
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

            source_pos = obj.temporary_sound_scene.virtual_sources{1}.position;
            x0 = cell2mat(cellfun( @(x) x.position, obj.temporary_sound_scene.loudspeaker_array, 'UniformOutput', false)');
            n0 = cell2mat(cellfun( @(x) x.orientation,  obj.temporary_sound_scene.loudspeaker_array, 'UniformOutput', false)');

            kR = xr-x0;
            kR = kR./sqrt(sum(kR.^2,2));
            kP = [0 0]-source_pos;
            kP = kP./sqrt(sum(kP.^2,2));
            [~,ix0] = min(sum((kR-kP).^2,2));
            k0 = kP;

           % LineLength = 0.65;headWidth = 5*1.5;
           % headLength = 5*1.5;
           % ah = annotation('arrow','headStyle','cback2','HeadLength',headLength,'HeadWidth',headWidth,'LineWidth',1,'Color','Black','LineStyle','--');
           % set(ah,'parent',gca);
           % set(ah,'position',[x0(ix0,1),x0(ix0,2), LineLength*k0(1,1), LineLength*k0(1,2)]);    
          % plot(x0(ix0,1),x0(ix0,2),'.w','MarkerSize',20)
            
            input = repmat( [zeros(length(obj.temporary_sound_scene.virtual_sources{1}.source_signal.time_series),1)] ...
                ,1 , length(obj.temporary_sound_scene.virtual_sources));
            obj.temporary_sound_scene.render_sound_scene(input,false, false);

            fs = obj.temporary_sound_scene.virtual_sources{1}.source_signal.fs;
            N = length(obj.temporary_sound_scene.virtual_sources{1}.source_signal.time_series);
            obj.t = (1:N)/fs;
            switch obj.excitation_type
                case 'harmonic'
                    obj.f0 = varargin{4};
                    input = repmat( cos(2*pi*obj.t*obj.f0)' ...
                        ,1 , length(obj.temporary_sound_scene.virtual_sources));
                    obj.temporary_sound_scene.render_sound_scene(input,false,false);
                case 'impulse'
                    W = varargin{4};
                    input = repmat( [hann(W);zeros(N-W,1)] ...
                        ,1 , length(obj.temporary_sound_scene.virtual_sources));
                    obj.temporary_sound_scene.render_sound_scene(input,false,false);
            end
        end

        function simulate(obj, t_slider)
            c = 343.1;
            t_slider = (t_slider-0.5)/10e1;
            [X,Y] = meshgrid(obj.x,obj.y);
            field = zeros(size(X));
            switch obj.excitation_type
                case 'harmonic'
                    t0 = obj.t(end) + t_slider;
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
                        t0 = max(mean_a,mean_b)/c + t_slider;
                    end
                    %title(obj.simulation_axes, sprintf('Reproduced field with impulse excitation at t_0 = %0.1f ms',1e3*t0));

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
            field_ref = 0;
            set(obj.simulation_result,'CData',field)
            switch obj.excitation_type
                case 'harmonic'
                    caxis([-1,1]*5*median(median(abs(field))))
                case 'impulse'
                    caxis([-1,1]*0.02)
            end
           set(gcf,'PaperPositionMode','auto');
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
                c = 343.1;
                transfer = transfer + 1/(4*pi)*s0.*exp(-1i*2*pi*freq_vec/c*R).*Dir0.'/R;
            end
            fig = figure('units','normalized','outerposition',[0.2 0.3 0.35 0.4]);
            fs = obj.temporary_sound_scene.virtual_sources{1}.source_signal.fs;
            freq = (0:length(transfer)-1)/N*fs;
            semilogx(freq,20*log10(abs(transfer)),'LineWidth',1.5);
            grid on
            xlim([40,20e3]);
            ylim([-10,10]+mean(20*log10(abs(transfer))));
            xlabel('f [Hz]')
            ylabel('P_{rec} [dB]')
            title(sprintf('Synthesized field transfer at x = %d, y = %d',...
                obj.temporary_sound_scene.receiver.position(1),...
                obj.temporary_sound_scene.receiver.position(2)));
       %     save('nara_left')
        end
    end
end

