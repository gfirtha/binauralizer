classdef hrtf_extrapolator < handle
    %HRTF_EXTRAPOLATOR Summary of this class goes here
    %   Detailed explanation goes here

    properties
        mode
        measurement_grid
        HRTF_database
        HRTF_SHT
        HRTF_spectrum
        HRTF_amp
        HRTF_phase
    end

    methods
        function obj = hrtf_extrapolator(hrtf, mode)
            %HRTF_EXTRAPOLATOR Construct an instance of this class
            %   Detailed explanation goes here
            obj.mode = mode;
            obj.HRTF_database = hrtf;

            name = strcat(obj.HRTF_database.GLOBAL_DatabaseName,'_',obj.HRTF_database.GLOBAL_ListenerShortName);
            if all(abs(hrtf.SourcePosition(:,2))<1e-9)
                obj.measurement_grid = circular_design(obj.HRTF_database.SourcePosition,name);
            else
                obj.measurement_grid = spherical_design('Direct',obj.HRTF_database.SourcePosition,name);
            end

            obj.HRTF_spectrum = fft( obj.HRTF_database.Data.IR,[],3 );
            obj.HRTF_amp = abs(obj.HRTF_spectrum);
            obj.HRTF_phase = unwrap(angle(obj.HRTF_spectrum),[],3);
            Nmax = max(get_maxN(obj.measurement_grid),15);
            fname = fullfile('Data/Spherical_designs/',[name,'_SHT.mat']);
            if exist(fname)
                temp = open(fname);
                obj.HRTF_SHT = temp.SHT;
            else
                SHT = obj.measurement_grid.SHT_on_sphere(obj.HRTF_database.Data.IR,Nmax);
                obj.HRTF_SHT = SHT;
                save(fname,'SHT');
            end
        end

        function hrtf_out = extrapolate_hrtf(obj,source, receiver)
            xr = [receiver.position, 1.7];
            xs = [source.position, source.height];
            o = [receiver.orientation,0];
            xe = (get_rotation_mx(o',[1,0,0]')*(xs-xr)')';
            switch obj.mode
                case 'nearest'
                    hrtf_out = obj.measurement_grid.interpolate(xe', obj.HRTF_database.Data.IR,'nearest')*...
                        mean(obj.measurement_grid.Phi0(:,end))/norm(xe);
                case 'linear'
                    out = obj.measurement_grid.interpolate(xe', obj.HRTF_spectrum,'linear')*...
                        mean(obj.measurement_grid.Phi0(:,end))/norm(xe);
                    hrtf_out = ifft(out,[],1);
                case'SHT'
                    [azim_out,elev_out,r_out] = cart2sph(xe(1),xe(2),xe(3));
                    Nmax = sqrt(size(obj.HRTF_SHT,1))-1;
                    H_rad = get_radial_filters(r_out, mean(obj.measurement_grid.Phi0(:,3)), size(obj.HRTF_SHT,3), obj.HRTF_database.Data.SamplingRate, Nmax);
                    HRTF_f = 2*obj.HRTF_SHT.*reshape(H_rad,[size(H_rad,1),1,size(H_rad,2)]);
                    hrtf_out = squeeze(ifft( ISHT( HRTF_f, sqrt(size(obj.HRTF_SHT,1))-1, pi/2-elev_out, azim_out, 'real') ,[],3,'symmetric'))';
                case 'WFS'
                    mode = 'focused';
                    switch class(obj.measurement_grid)
                        case 'spherical_design'
                            x0 = obj.measurement_grid.x0;
                            n_in = -obj.measurement_grid.n_out;
                            dS = obj.measurement_grid.voronoi_cells.dS;
                            c = 343;
                            vkp = x0-xe;

                            R0 = sqrt(sum((x0).^2,2));
                            R = sqrt(sum(vkp.^2,2));
                            vkp_u = vkp./R;
                            vkp_n = sum( vkp_u.*n_in ,2 );
                            win = sum(n_in.*xe/norm(xe),2); win(win<=0) = 0; win = win.^2;
                            amp = -2/(c*4*pi).*win.*vkp_n.*dS ./R;
                            delay = -R/c;
                            delay = delay - min(delay);

                            fs = obj.HRTF_database.Data.SamplingRate;
                            Nt = size(obj.HRTF_database.Data.IR,3);
                            omega = 2*pi*(0:Nt-1)'/Nt*fs;
                            vkt = vkp_u-vkp_n.*n_in;
                            vkt_length = sqrt(sum( vkt.^2,2));
                            dx_dir = obj.measurement_grid.get_delauney_dx(vkt);

                            wc =  c*pi./dx_dir./abs(vkt_length);
                            wc(isinf(wc)) = 1e20;

                            [Wc,W] = meshgrid(wc,omega);
                            Nbut = 4;
                            AA = 1./sqrt(1+(W./Wc).^(2*Nbut));
                            D_abs = reshape(  1i*amp.*omega'.*AA.', [size(amp,1),1,size(omega,1)]);
                            D_phase = reshape( -delay*omega', [size(amp,1),1,size(omega,1)]);
                            HRTF_out = squeeze( sum( obj.HRTF_amp.*D_abs.*exp(1i*(obj.HRTF_phase+D_phase)), 1) ).';
                            hrtf_out = ifft( HRTF_out,[],1,'symmetric');
                        case 'circular_design'
                            x0 = obj.measurement_grid.x0;
                            n_in = -obj.measurement_grid.n_out;
                            edges = obj.measurement_grid.edges;
                            dS = obj.measurement_grid.get_dS(x0,edges);
                            %                         [x0,n_in,dS,edges] = obj.measurement_grid.append_point(xe);
                            c = 340;
                            switch mode
                                case 'normal'
                                    xe = xe(1:2);
                                case 'focused'
                                    xe = -xe(1:2);
                            end
                            vkp = x0-xe;
                            R0 = sqrt(sum((x0).^2,2));
                            R = sqrt(sum(vkp.^2,2));
                            vkp_u = vkp./R;
                            vkp_n = sum( vkp_u.*n_in ,2 );

                            switch mode
                                case 'normal'
                                    if ~all(vkp_n<0)
                                        dref = sqrt(R.*R0./(R+R0));
                                        win = vkp_n>=0;
                                        amp = -sqrt(1/(c*2*pi)).*win.*dref.*vkp_n.*dS ./R;
                                        delay = R/c;
                                    else
                                        dref = sqrt(R.*R0./(R-R0));
                                        win = sum(n_in.*xe/norm(xe),2); win(win>=0) = 0; win = win.^4;
                                        amp = sqrt(-1/(c*2*pi)).*win.*dref.*vkp_n.*dS ./R;
                                        delay = -R/c;
                                    end
                                case 'focused'
                                    dref = sqrt(R.*R0./(R-R0));
                                    dref(isinf(dref)) = 1;
                                    win = sum(n_in.*xe/norm(xe),2); win(win<=0) = 0; win = win.^2;
                                    amp = sqrt(1/(c*2*pi)).*win.*dref.*vkp_n.*dS ./R;
                                    delay = -R/c;
                            end
                            delay = delay - min(delay);

                            fs = obj.HRTF_database.Data.SamplingRate;
                            Nt = size(obj.HRTF_database.Data.IR,3);
                            omega = 2*pi*(0:Nt-1)'/Nt*fs;

                            vkt = vkp_u-vkp_n.*n_in;
                            vkt_length = sqrt(sum( vkt.^2,2));
                            dx_dir = obj.measurement_grid.get_directive_dx(vkt,x0,edges);
                            [~,ixmin] = sort(sum((x0+xe/norm(xe)*mean(R0)).^2,2));
                            wc =  c*pi./dx_dir./abs(vkt_length);
                            wc(ixmin(1)) = inf;
                            wc(isinf(wc)) = 1e20;

                            [Wc,W] = meshgrid(wc,omega);
                            Nbut = 4;
                            AA = 1./sqrt(1+(W./Wc).^(2*Nbut));
                            D_abs = reshape(  amp.*sqrt(1i*omega').*AA.', [size(amp,1),1,size(omega,1)]);
                            D_phase = reshape( -delay*omega', [size(amp,1),1,size(omega,1)]);
                            HRTF_out =  squeeze( sum( obj.HRTF_amp.*D_abs.*exp(1i*(obj.HRTF_phase+D_phase)), 1) ).';
                            hrtf_out = ifft( HRTF_out,[],1,'symmetric');

                    end
            end
        end
    end
end

