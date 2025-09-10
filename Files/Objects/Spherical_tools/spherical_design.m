classdef spherical_design < handle
    %SPHERICAL_DESIGN Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name    % Design name: name of HRTF-set, or the design type
        N0      % number of sampling points
        x0      % sampling positions, Descartes
        Phi0    % sampling positions, Spherical
        delauney_trias

        delauney_circles
        voronoi_cells
        min_vol_ellipses
        quadrature_points
    end

    methods
        function obj = spherical_design(varargin)
            %SPHERICAL_DESIGN Construct an instance of this class
            %   Detailed explanation goes here

            % Definition mode:
            %   Direct:  sample grid is given directly
            %   Design:  design mode is given with N0 being the no. of
            %            sampling points
            def_mode = varargin{1};
            switch def_mode
                case 'Direct'
                    obj.Phi0 = varargin{2};
                    if max(abs(obj.Phi0(:,1)))>2*pi
                        obj.Phi0(:,1:2) = obj.Phi0(:,1:2)*pi/180;
                        [x0,y0,z0] = sph2cart(obj.Phi0(:,1),obj.Phi0(:,2),obj.Phi0(:,3));
                        obj.x0 = [x0 y0 z0]*(mean(obj.Phi0(:,3))).^2./(sum( [x0 y0 z0].^2,2));
                    end
                    obj.name = varargin{3};
                    obj.N0 = length(obj.x0);
                case 'Design'

                    obj.N0 = varargin{2};
                    R = varargin{4};
                    obj.name = [varargin{3},'_',mat2str(varargin{2}),'_R',mat2str(R)];
                    switch varargin{3}
                        case 'T'
                            files = dir('Data/Spherical_designs/t_design_points');
                    end
                    for n = 3 : length(files)
                        temp = split(files(n).name,'.');
                        N_available(n-2) = str2num(temp{2});
                    end
                    [~,file_ix] = min(abs(N_available-obj.N0));
                    fid = fopen(fullfile(files(file_ix+2).folder,files(file_ix+2).name), 'r');
                    x0 = fscanf(fid,'%f');
                    obj.x0 = R*reshape(x0,[3,size(x0,1)/3])';
                    [azim0,elev0,r0] = cart2sph(obj.x0(:,1),obj.x0(:,2),obj.x0(:,3));
                    obj.Phi0 = [azim0,elev0,r0];
                    fclose(fid);

            end

            fname = fullfile('Data','Spherical_designs',[obj.name,'.mat']);
            if isfile(fname)
                in = open(fname);
                obj.voronoi_cells = in.voronoi_cells;
                obj.min_vol_ellipses = in.min_vol_ellipses;
                obj.delauney_trias = in.delauney_trias;
                obj.delauney_circles = in.delauney_circles;
            else
                R0 = mean(obj.Phi0(:,3));
                [ face_num, tria_faces ] = sphere_delaunay( size(obj.x0,1) , obj.x0'/R0 );
                obj.delauney_trias = tria_faces;
                obj.get_delauney_circles;
                [v_vor, K, n_in, dS] = get_voronoi_polygons(obj.delauney_trias, obj.x0);
                %[v_vor, K, n_in, dS] = get_voronoi_polygons(obj.delauney_trias, obj.x0);

                obj.voronoi_cells = struct('vertices',[],'connectivity',[],'adjacency',{},'dS',[], 'n_in',[]);
                obj.voronoi_cells(1).vertices = v_vor;
                obj.voronoi_cells(1).connectivity = K';
                obj.voronoi_cells(1).dS = dS';
                obj.voronoi_cells(1).n_in = n_in';
                obj.voronoi_cells(1).adjacency = obj.get_adjacency(size(obj.x0),obj.delauney_trias);
                obj.get_mv_ellipses;
                voronoi_cells = obj.voronoi_cells;
                min_vol_ellipses = obj.min_vol_ellipses;
                delauney_trias = obj.delauney_trias;
                delauney_circles = obj.delauney_circles;
                save(fname,'voronoi_cells','min_vol_ellipses','delauney_trias','delauney_circles');
            end
        end

        function obj = translate_sphere(obj,dx)
            obj.x0 = obj.x0 + dx;
            [azim0,elev0,r0] = cart2sph(obj.x0(:,1),obj.x0(:,2),obj.x0(:,3));
            obj.Phi0 = [azim0,elev0,r0];

        end

        function [f_out, x_interp] = interpolate(obj, x_target, f_in,mode)
            % interpolates input function to poistions over the sphere into
            % the direction x_target
            switch mode
                case 'linear'
                    [ixout,x_interp] = obj.find_triangle2point(x_target);
                    ix = obj.delauney_trias(:,ixout(1));
                    vertices = obj.x0(ix,:);
                    w = vertices'\x_interp;
                    f_out = squeeze(sum(f_in(ix,:,:).*w,1)).';
                case 'nearest'
                    [~,ix] = min(sum(( obj.x0- x_target'./norm(x_target).*mean(obj.Phi0(:,3)) ).^2,2));
                    f_out = squeeze(f_in(ix,:,:))';
            end
        end


        function adjacency = get_adjacency(obj,N,trias)
            for n = 1 : N
                ixs = (trias(1,:)==n)|(trias(2,:)==n)|(trias(3,:)==n);
                adjacency{n} = setdiff(unique(trias(:,ixs)),n);
            end
        end

        function get_delauney_circles(obj)
            xc = zeros(size(obj.delauney_trias,2),3);
            r = zeros(size(obj.delauney_trias,2),1);
            fi0 = zeros(size(obj.delauney_trias,2),1);
            for n = 1 : size(obj.delauney_trias,2)
                tria = obj.delauney_trias(:,n);
                xv = obj.x0(tria,:);
                [xc(n,:),r(n),v1n,v2nb] = circlefit3d(xv(1,:),xv(2,:),xv(3,:));
                fi0(n) = atan(r(n)/norm(xc(n,:)));
            end
            xc = xc./sqrt(sum(xc.^2,2));
            obj.delauney_circles = struct('center',xc,'fi',fi0);
        end

        function obj = get_mv_ellipses(obj)
            K = obj.voronoi_cells.connectivity;
            vertices = obj.voronoi_cells.vertices;
            Tr_mx = zeros(2,3,length(K));
            R = zeros(2,length(K));

            for n = 1 : length(K)
                p1 = vertices(:,K{n}(1));
                p2 = vertices(:,K{n}(2));
                p3 = vertices(:,K{n}(end));
                v1 = (p2-p1)/norm(p2-p1);
                v2 = (p3-p1)/norm(p3-p1);
                Tmx = [v1 ,v2 ];
                points = [obj.x0(n,:);(obj.x0(obj.voronoi_cells.adjacency{n},:))]';
                tr_points = Tmx\(points-p1);
                [A, c0] = mvee_fit(tr_points);
                [~, D, V] = svd(A);
                R(:,n) = 1./sqrt(diag(D));
                Tr_mx(:,:,n) = pinv(Tmx*V);
                center(:,n) = Tmx*c0 + p1;

            end
            obj.min_vol_ellipses = struct('Tmx',Tr_mx,'R',R,'center',center);
        end

        function [r,kr0] = get_dx(obj,k0)
            k0 = k0./sqrt(sum(k0.^2,1));
            kr = k0'-sum(obj.voronoi_cells.n_in.*k0',2).*obj.voronoi_cells.n_in;
            kr0 = sqrt(sum(kr.^2,2));
            kr = kr./kr0;
            for n = 1 : size(kr,1)
                wtr(:,n) = squeeze(obj.min_vol_ellipses.Tmx(:,:,n))*(kr(n,:))';
            end
            r = ( prod(obj.min_vol_ellipses.R,1)./sqrt(sum((obj.min_vol_ellipses.R.*flipud(wtr)).^2,1)) )';
            r(isnan(r)) = 1e-10;
        end

        function dS = get_tria_areas(obj,varargin)
            if isempty(varargin)
                dS = [];
                for n = 1 : length(obj.x0)

                    [row,col] = find(obj.delauney_trias == n);
                    trias = obj.delauney_trias(:,col);
                    area = [];
                    for m = 1 : length(col)
                        tria = trias(:,m);
                        vertices = obj.x0(tria,:)';
                        v1 = vertices(:,2) - vertices(:,1);
                        v2 = vertices(:,3) - vertices(:,1);
                        area(m) = norm(cross(v1, v2, 1)) / 2;
                    end
                    dS(n) = sum(area)/3;
                end
            else
                ix0 = varargin{1};
                [row,col] = find(obj.delauney_trias == ix0);
                trias = obj.delauney_trias(:,col);
                area = [];
                for m = 1 : length(col)
                    tria = trias(:,m);
                    vertices = obj.x0(tria,:)';
                    v1 = vertices(:,2) - vertices(:,1);
                    v2 = vertices(:,3) - vertices(:,1);
                    area(m) = norm(cross(v1, v2, 1)) / 2;
                end
                dS = sum(area)/3;
            end
        end

        function [ix_tria,x_intersection] = find_triangle2point(obj,xp)
            % Calculates the Delauney triangle that is intersected by a
            % direction, given by position vector xp
            for n = 1 : size(obj.delauney_trias ,2)
                tria = obj.delauney_trias(:,n);
                p0 = obj.x0(tria(1),:);
                v1 = obj.x0(tria(2),:) - p0;
                v2 = obj.x0(tria(3),:) - p0;
                normal = cross(v1',v2');
                normal = normal ./ sqrt( sum( normal.^2 ) );
                d(n) = dot(p0,normal)/dot(xp/norm(xp),normal);
                xp_ = d(n)*xp/norm(xp);
                w(:,n) = [v1',v2']\(xp_-p0');
            end
            w = (round(1e8*w))/1e8;
            ix_tria = find(w(1,:)>=0 & w(2,:)>=0 & sum(w,1)<=1&d>0);
            x_intersection = xp/norm(xp)*d(ix_tria(1));
        end

        function [ix_pos,weights,x_extrap] = get_interpolation_weights(obj,x_target)
            [ixout,x_extrap] = obj.find_triangle2point(x_target);
            ix_pos = obj.delauney_trias(:,ixout(1));
            vertices = obj.x0(ix_pos,:);
            weights = vertices'\x_extrap;
            ix_pos(weights<1e-8) = [];
            weights(weights<1e-8) = [];
        end
        function [f_out, x_interp] = interpolate2position(obj, x_target, f_in)
            [ix,w,x_interp] = get_interpolation_weights(obj,x_target);
            amp = abs(f_in(ix,:,:));
            phase = unwrap(angle(f_in(ix,:,:)),[],3);
            amp_interp = squeeze(sum(amp.*w,1)).';
            phase_interp = squeeze(sum(phase.*w,1)).';
            f_out = amp_interp.*exp(1i*phase_interp);
        end

        function [obj] = get_quadrature_points(obj)
            order = 10;
            [xx, ww] = gaussquad2(order, 3);
            tria_set = ShapeSet.LinearTria;
            w = tria_set.eval(xx);

            x_quad = [];
            w_quad = [];
            w_interp = [];
            ind_interp = [];
            Ssphere = 0;
            for n = 1 : size(obj.delauney_trias,2)
                tria = obj.delauney_trias(:,n)';

                p0 = obj.x0(tria(1),:);
                v1 = obj.x0(tria(2),:) - p0;
                v2 = obj.x0(tria(3),:) - p0;
                normal = cross(v1',v2');
                S = norm(cross(v1',v2'))/2;
                Ssphere = Ssphere + S;
                normal = normal ./ sqrt( sum( normal.^2 ) );
                W = [v1',v2',normal];
                x_quad  = [x_quad ;  (W*[xx,zeros(size(xx,1),1)]')'+p0 ];
                w_quad = [w_quad; ww*S/0.5];
                w_interp = [w_interp; w];
                ind_interp = [ind_interp; repmat(tria,[size(xx,1),1])];
            end
            obj.quadrature_points = struct('quad_pos',x_quad,'quad_w',w_quad,'interp_w',w_interp,'interp_pos',ind_interp);
        end

        function [obj] = get_adaptive_quadrature_points(obj)
            edges = unique(sort([obj.delauney_trias([1,2],:)';obj.delauney_trias([1,3],:)';obj.delauney_trias([2,3],:)'],2),'rows');
            Le = sqrt(sum((obj.x0(edges(:,1),:)-obj.x0(edges(:,2),:)).^2,2));
            Lmin = min(Le);
            Lmax = max(Le);
            min_order = 6;
            max_order = 20;

            x_quad = [];
            w_quad = [];
            w_interp = [];
            ind_interp = [];
            for n = 1 : size(obj.delauney_trias,2)
                tria = obj.delauney_trias(:,n)';
                p0 = obj.x0(tria(1),:);
                v1 = obj.x0(tria(2),:) - p0;
                v2 = obj.x0(tria(3),:) - p0;
                normal = cross(v1',v2');
                Stria = norm(cross(v1',v2'))/2;

                L = max([norm(v1),norm(v2),norm(v1-v2)]);
                s0 = (L-Lmin)/(Lmax-Lmin);
                order = round(((max_order-min_order)*s0+min_order));
                %                [xx, ww] = gaussquad2(order, 3);
                Q = quadtriangle(order,'Domain',[0,0;1 0;0 1]);
                xx = Q.Points;
                ww = Q.Weights;
                tria_set = ShapeSet.LinearTria;
                w = tria_set.eval(xx);

                normal = normal ./ sqrt( sum( normal.^2 ) );
                W = [v1',v2',normal];
                x_quad  = [x_quad ;  (W*[xx,zeros(size(xx,1),1)]')'+p0 ];
                w_quad = [w_quad; ww*Stria/0.5];
                w_interp = [w_interp; w];
                ind_interp = [ind_interp; repmat(tria,[size(xx,1),1])];
            end
            obj.quadrature_points = struct('quad_pos',x_quad,'quad_w',w_quad,'interp_w',w_interp,'interp_pos',ind_interp);
        end

        function [adjacency] = get_delanuey_adjacency(obj)
            for n = 1 : size(obj.x0,1)
                [~,cols] = find(obj.delauney_trias-n == 0);
                adjacency{n} = unique(obj.delauney_trias(:,cols))';
            end
        end
        function [phase_out] = unwrap_phase(obj,phase_in,ix_start)
            adjacency = obj.get_delanuey_adjacency;
            figure
            N0 = 100;
            pl = obj.plot_data(phase_in(:,N0));
            axis equal tight
            clim([min(phase_in(:,100)),max(phase_in(:,N0))])
            R = get_rotation_mx(obj.x0(ix_start,:)',[0,0,1]');
            Xout = (R*obj.x0')';
            [~,elev] = cart2sph(Xout(:,1),Xout(:,2),Xout(:,3));
            zenith = pi/2-elev;

            ixs_to_unwrap = setdiff((1:size(obj.x0,1)),ix_start);
            ixs_unwrapped = ix_start;
            i = 0;
            while ~isempty(ixs_to_unwrap)
                i = i+1
                % Choose neighbours not unwrapped
                neighbours_to_unwrap = setdiff(adjacency{ix_start},ixs_unwrapped);
                % Here comes unwrapping
                phase_ref = phase_in(ix_start,:);
                phases = phase_in(neighbours_to_unwrap,:);
                phase_in(neighbours_to_unwrap,:) = phase_ref + mod( phases-phase_ref+pi, 2*pi )-pi;
                for m = 1 : length(neighbours_to_unwrap)
                    ixs_ref = find(zenith(ixs_unwrapped)<zenith(neighbours_to_unwrap(m)) & ...
                        zenith(ixs_unwrapped)>zenith(neighbours_to_unwrap(m))-pi/8 );
                    ixs_ref = ixs_unwrapped(ixs_ref);
                    ix_to_mod = 0;
                    while ~isempty(ix_to_mod)
                        ix_to_mod = find(phase_in(neighbours_to_unwrap(m),:) > mean(phase_in(ixs_ref,:),1));
                        phase_in(neighbours_to_unwrap(m),ix_to_mod) = phase_in(neighbours_to_unwrap(m),ix_to_mod)-2*pi;
                    end
                end
                % Now they are unwrapped
                ixs_unwrapped = [ixs_unwrapped, neighbours_to_unwrap];
                ixs_to_unwrap = setdiff(ixs_to_unwrap,neighbours_to_unwrap);

                % Possible next start:
                % Choose from unwrapped data, that has wrapped
                % neighbour
                N_wrapped_neigbours = cell2mat(cellfun(@(x) numel(intersect(x,ixs_to_unwrap)), adjacency(ixs_unwrapped),'UniformOutput',0));
                % possible starting points
                ixs_possible = ixs_unwrapped(N_wrapped_neigbours>0);
                % Check phase difference from unwrapped neighbours of possible starting points
                unwrapped_adj = cellfun(@(x) intersect(x,ixs_unwrapped), adjacency(ixs_possible),'UniformOutput',0);
                phase_adj = cellfun(@(x) phase_in(x,:), unwrapped_adj,'UniformOutput',0);

                phase_diff = cell2mat( cellfun( @(x,y) sum(sum(abs(x-y)<pi,2),1),phase_adj, mat2cell(phase_in(ixs_possible,:),ones(length(ixs_possible),1))','UniformOutput',0) );
                [~,ix0] = max(phase_diff);
                ix_start = ixs_possible(ix0);
                set(pl,'CData',phase_in(:,100))
                drawnow
                obj;
            end
            obj;
            phase_out = phase_in;

        end

        function [output,amp,phase] = interpolate2quadrature(obj,input)
            if ~isreal(input)
                A_in = abs(input);
                phase_in = unwrap(angle(input),[],3);
            end
            tic
            obj.get_adaptive_quadrature_points;
            toc
            sizein = size(input);
            output = zeros([size(obj.quadrature_points.quad_pos,1),sizein(2:end)]);
            amp = zeros([size(obj.quadrature_points.quad_pos,1),sizein(2:end)]);
            phase = zeros([size(obj.quadrature_points.quad_pos,1),sizein(2:end)]);
            for n = 1 : size(obj.quadrature_points.quad_pos,1)
                if  isreal(input)
                    output(n,:,:) =  sum(obj.quadrature_points.interp_w(n,:)'.*input(obj.quadrature_points.interp_pos(n,:)',:,:),1);
                else
                    A_interp = sum(obj.quadrature_points.interp_w(n,:)'.*A_in(obj.quadrature_points.interp_pos(n,:)',:,:),1);
                    phase_interp = sum(obj.quadrature_points.interp_w(n,:)'.*phase_in(obj.quadrature_points.interp_pos(n,:)',:,:),1);
                    output(n,:,:) =  A_interp.*exp(1i*phase_interp);
                    amp(n,:,:) =  A_interp;
                    phase(n,:,:) =  phase_interp;
                end
            end
        end

        function [S] = get_daluney_areas(obj)
            p0 = obj.x0(obj.delauney_trias(1,:)',:);
            v1 = obj.x0(obj.delauney_trias(2,:)',:)-p0;
            v2 = obj.x0(obj.delauney_trias(3,:)',:)-p0;
            S = sqrt(sum((cross(v1,v2,2)).^2,2))/2;
        end

        function [obj] = append_point(obj,xp)
            if ~(min(sum((obj.x0./sqrt(sum(obj.x0.^2,2)) - xp'/norm(xp)).^2,2)) < 1e-9)
                % Calculate new delauney triangles
                % First get triangles to be deleted,
                % d_p = index of neighbouring points
                ix_p = size(obj.x0,1) + 1;
                obj.x0 = [obj.x0;  xp'/norm(xp)*mean(obj.Phi0(:,3))  ];
                cents = obj.delauney_circles.center./sqrt(sum((obj.delauney_circles.center).^2,2));
                near_trias = find(acos(cents*xp/norm(xp)) < 1*obj.delauney_circles.fi);
                neighbour_x = [unique(obj.delauney_trias(:,near_trias));ix_p];

                % Calculate new Delauney triangles
                vx = get_rotation_mx(xp,[0,0,-1])*obj.x0(neighbour_x,:)';
                x = vx(1,:)./(1-vx(3,:));
                y = vx(2,:)./(1-vx(3,:));
                DT = delaunay(x,y);
                new_trias = neighbour_x(DT)';
                [~,col] = find(new_trias == ix_p);
                new_trias = new_trias(:,col);

                voronoi_vxs = zeros(3,size(new_trias,2));
                for n = 1 : size(new_trias,2)
                    voronoi_vxs(:,n) = get_circumcircle(new_trias(:,n), obj.x0',mean(sqrt(sum(obj.x0.^2,2))) );
                end
                obj.delauney_trias = [obj.delauney_trias, new_trias];
                obj.voronoi_cells.vertices = [obj.voronoi_cells.vertices, voronoi_vxs];

                for n = 1 : length(neighbour_x)
                    [rows,cols] = find(obj.delauney_trias' == neighbour_x(n));
                    rows = setdiff(rows, near_trias);
                    V = obj.voronoi_cells.vertices(:,rows') - mean(obj.voronoi_cells.vertices(:,rows'),2);
                    [order, nin] = sort_vertices(V./sqrt(sum(V.^2,1)));
                    obj.voronoi_cells.n_in(neighbour_x(n),:)  = nin';
                    obj.voronoi_cells.connectivity{neighbour_x(n)} = rows(order);
                end
                obj.voronoi_cells.dS = get_spherical_polygon_area(obj.voronoi_cells.connectivity,obj.voronoi_cells.vertices);
                %   obj.delauney_trias(:,near_trias) = [];
                %   obj.voronoi_cells.vertices(:,near_trias) = [];

            end
        end

        function [SHT_out] = SHT_on_sphere(obj,in,Nmax)
            input = fft(in,[],3);
            output = obj.interpolate2quadrature(input);

            [azim, elev]= cart2sph(obj.quadrature_points.quad_pos(:,1),obj.quadrature_points.quad_pos(:,2),obj.quadrature_points.quad_pos(:,3));
            zenith = pi/2 - elev;
            S = 4*pi*mean(obj.Phi0(:,3))^2;
            Y = getSpherHarmMx( zenith, azim, Nmax, 'real' );
            SHT_out = zeros(size(Y,2),size(in,2),size(in,3));
            for m = 1 : size(in,3)
                SHT_out(:,:,m) = (4*pi*obj.quadrature_points.quad_w.*Y./S)'*squeeze(output(:,:,m));
            end

        end

        % Plotting functions
        function plot_ellipse(obj,n)
            u = linspace(0, 2*pi, 90);
            X_ell = pinv(obj.min_vol_ellipses.Tmx(:,:,n)) * [ obj.min_vol_ellipses.R(1,n)*cos(u);...
                obj.min_vol_ellipses.R(2,n)*sin(u) ] + obj.min_vol_ellipses.center(:,n);
            hold on
            plot3(X_ell(1,:)',X_ell(2,:)',X_ell(3,:)','k','LineWidth',1.5);

        end
        function plot_dx(obj,k0,n)
            xo = obj.x0(n,:);
            quiver3(xo(1),xo(2),xo(3),k0(1),k0(2),k0(3),'b');
            k0 = k0 / norm(k0);
            kt = k0'-(obj.voronoi_cells.n_in(n,:)*k0).*obj.voronoi_cells.n_in(n,:);
            kt = kt / norm(kt);
            r = obj.get_dx(k0);
            v = r(n)*kt + obj.min_vol_ellipses.center(:,n)';
            scatter3(v(1),v(2),v(3),20,'b','filled');
        end
        function scatter_sphere(obj)
            scatter3(obj.x0(:,1),obj.x0(:,2),obj.x0(:,3),10,'k','filled');
            axis equal tight
            view(3)
        end
        function plot_voronoi(obj)
            K = obj.voronoi_cells.connectivity;
            vertices = obj.voronoi_cells.vertices;
            for n = 1 : size(K,1)
                patch(vertices(1,K{n}),vertices(2,K{n}),vertices(3,K{n}),'red')
                hold on
            end
            axis equal tight
            view(3)
        end
        function plot_poly(obj,K,vertices)
            for n = 1 : size(K,2)
                patch(vertices(1,K{n}),vertices(2,K{n}),vertices(3,K{n}),'red')
                hold on
            end
            axis equal tight
            view(3)
        end

        function s = plot_data(obj,in)
            s = trisurf(obj.delauney_trias.',obj.x0(:,1),obj.x0(:,2),obj.x0(:,3),in);
            axis equal tight

        end
        function plot_delauney_patch(varargin)
            obj = varargin{1};
            if length(varargin)>1
                fig = varargin{2};
                hold on
                for n = 1 : size(obj.delauney_trias,2)
                    patch(fig,obj.x0(obj.delauney_trias(:,n),1)',obj.x0(obj.delauney_trias(:,n),2)',obj.x0(obj.delauney_trias(:,n),3)','red','LineStyle','none')
                end
            else
                hold on
                for n = 1 : size(obj.delauney_trias,2)
                    patch(obj.x0(obj.delauney_trias(:,n),1)',obj.x0(obj.delauney_trias(:,n),2)',obj.x0(obj.delauney_trias(:,n),3)','red')
                end
                axis equal tight
                view(3)

            end
        end
        function plot_delauney(obj, varargin)
            if isempty(varargin)
                all_edges = unique(sort([obj.delauney_trias(1:2,:),obj.delauney_trias(2:3,:),obj.delauney_trias([1,3],:)]',2),'rows');
                hold on
                for n = 1 : size(all_edges,1)
                    startpt = obj.x0(all_edges(n,1),:);
                    endpt = obj.x0(all_edges(n,2),:);
                    plot3([startpt(1),endpt(1)],[startpt(2),endpt(2)],[startpt(3),endpt(3)],'b')
                end
                axis equal tight
                view(3)
            else
                ix = varargin{1};
                all_edges = unique(sort([obj.delauney_trias(1:2,ix),obj.delauney_trias(2:3,ix),obj.delauney_trias([1,3],ix)]',2),'rows');
                for n = 1 : size(all_edges,1)
                    startpt = obj.x0(all_edges(n,1),:);
                    endpt = obj.x0(all_edges(n,2),:);
                    plot3([startpt(1),endpt(1)],[startpt(2),endpt(2)],[startpt(3),endpt(3)],'k','LineWidth',1.5)
                end
            end
        end

    end
end

