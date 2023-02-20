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
            obj.temporary_sound_scene = varargin{1};
            gui = varargin{2};
            obj.excitation_type = varargin{3};
            dx = 1e-2;
            obj.x = [gui.main_axes.XLim(1):dx:gui.main_axes.XLim(2)]';
            obj.y = [gui.main_axes.YLim(1):dx:gui.main_axes.YLim(2)]';
            [X,Y] = meshgrid(obj.x,obj.y);
            obj.fig = figure('units','normalized','outerposition',[0 0.3 0.4 0.7]);
            obj.simulation_axes = axes;
            obj.simulation_result  = pcolor(obj.x,obj.y,zeros(size(X)));
            shading interp
            
            axis equal tight
            copyobj(gui.receiver, obj.simulation_axes);
            for n = 1 : length(gui.virtual_source_points)
                copyobj(gui.virtual_source_points{n}, obj.simulation_axes);
            end
            for n = 1 : length(gui.binaural_source_points)
                copyobj(gui.binaural_source_points{n}, obj.simulation_axes);
            end
            xlim(gui.main_axes.XLim)
            ylim(gui.main_axes.YLim)
            caxis([-1,1])
            xlabel('x -> [m]');
            ylabel('y -> [m]');
            
            input = repmat( [zeros(length(obj.temporary_sound_scene.binaural_sources{1}.source_signal.time_series),1)] ...
                ,1 , length(obj.temporary_sound_scene.binaural_sources));
            obj.temporary_sound_scene.scene_renderer.render(input);
            
            fs = obj.temporary_sound_scene.binaural_sources{1}.source_signal.fs;
            N = length(obj.temporary_sound_scene.binaural_sources{1}.source_signal.time_series);
            obj.t = (1:N)/fs;
            switch obj.excitation_type
                case 'harmonic'
                    obj.f0 = varargin{4};
                    input = repmat( cos(2*pi*obj.t*obj.f0)' ...
                        ,1 , length(obj.temporary_sound_scene.binaural_sources));
                    obj.temporary_sound_scene.scene_renderer.render(input);
                case 'impulse'
                    W = varargin{4};
                    W = 5;
                    input = repmat( [hann(W);zeros(N-W,1)] ...
                        ,1 , length(obj.temporary_sound_scene.binaural_sources));
                    obj.temporary_sound_scene.scene_renderer.render(input);
            end            
        end
        
        function simulate(obj, t_slider)
            c = 343.1;
            t_slider = (t_slider-0.5)/2e1;
            [X,Y] = meshgrid(obj.x,obj.y);
            field = zeros(size(X));
            switch obj.excitation_type
                case 'harmonic'
                    t0 = obj.t(end) + t_slider;
                    title(obj.simulation_axes , sprintf('Reproduced field with harmonic excitation at f_0 = %d Hz',obj.f0));
                case 'impulse'
                    if isempty(obj.temporary_sound_scene.virtual_sources)
                        mean_d =  mean(sqrt(sum((bsxfun(@minus, obj.temporary_sound_scene.receiver.position,...
                            (cell2mat(cellfun( @(x) x.position', obj.temporary_sound_scene.binaural_sources, 'UniformOutput', false)))')).^2,2)));
                        t0 = mean_d/c  + t_slider;
                    else
                        mean_a =  mean(sqrt(sum((bsxfun(@minus, obj.temporary_sound_scene.receiver.position,...
                            (cell2mat(cellfun( @(x) x.position', obj.temporary_sound_scene.binaural_sources, 'UniformOutput', false)))')).^2,2)));
                        mean_b =  mean(sqrt(sum((bsxfun(@minus, obj.temporary_sound_scene.receiver.position,...
                            (cell2mat(cellfun( @(x) x.position', obj.temporary_sound_scene.virtual_sources, 'UniformOutput', false)))')).^2,2)));
                        t0 = max(mean_a,mean_b)/c + t_slider;
                    end
                    title(obj.simulation_axes, sprintf('Reproduced field with impulse excitation at t_0 = %0.1f ms',1e3*t0));
                    
            end
            
            for n = 1 : length(obj.temporary_sound_scene.binaural_sources)
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
                            obj.temporary_sound_scene.binaural_sources{n});
                        Dir0 = obj.temporary_sound_scene.scene_renderer.directivity_tables{idx}.directivity_mx;
                        in =  fft(obj.temporary_sound_scene.binaural_sources{n}.source_signal.time_series,size(Dir0,1));
                        Directive_signal = ifft(bsxfun( @times, Dir0, in ),[],1);
                        Directive_signal = Directive_signal(1:end/2,:);
                        theta = obj.temporary_sound_scene.scene_renderer.directivity_tables{idx}.theta;
                        [Theta, T] = meshgrid(theta,obj.t');
                        
                        xr = [X(:)-obj.temporary_sound_scene.scene_renderer.binaural_renderer{n}.binaural_source.position(1),...
                            Y(:)-obj.temporary_sound_scene.scene_renderer.binaural_renderer{n}.binaural_source.position(2)];
                        R = reshape(sqrt( sum(xr.^2,2) ) , size(X,1),size(X,2));
                        v = bsxfun(@times,xr,1./(sqrt( sum(xr.^2,2) )));
                        vn = v*obj.temporary_sound_scene.scene_renderer.binaural_renderer{n}.binaural_source.orientation';
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
                    caxis([-1,1]*abs(field(round(end/2),round(end/2))))
                case 'impulse'
                    caxis([-1,1]/30)
            end
        end
    end
end

