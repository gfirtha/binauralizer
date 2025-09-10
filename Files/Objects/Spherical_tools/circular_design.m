classdef circular_design < handle
    %SPHERICAL_DESIGN Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name    % Design name: name of HRTF-set, or the design type
        N0      % number of sampling points
        x0      % sampling positions, Descartes
        Phi0    % sampling positions, Spherical
        n_out   % outward normal
        edges
        dS
        
    end

    methods
        function obj = circular_design(varargin)
            obj.Phi0 = varargin{1}(:,[1,3]);
            if max(abs(obj.Phi0(:,1)))>2*pi
                obj.Phi0(:,1) = obj.Phi0(:,1)*pi/180;
            end
            [x0,y0] = pol2cart(obj.Phi0(:,1),obj.Phi0(:,2));
            obj.x0 = [x0 y0]*(mean(obj.Phi0(:,2))).^2./(sum( [x0 y0].^2,2));
            obj.name = varargin{2};
            obj.N0 = length(obj.x0);
            obj.n_out = obj.x0./sqrt(sum(obj.x0.^2,2));
            obj.edges = convhulln(obj.x0);
            obj.dS = obj.get_dS(obj.x0,obj.edges);
        end

        function [ix_edge,x_intersection,w0] = find_edge(obj,xp)
            for n = 1 : size(obj.edges,1)
                w(:,n) = obj.x0(obj.edges(n,:),:)'\xp;
            end
            ix_edge = find(w(1,:)>=-1e-10&w(2,:)>=-1e-10);
            ix_edge = ix_edge(1);
            w0 = w(:,ix_edge)/norm(w(:,ix_edge),1);
            x_intersection = obj.x0(obj.edges(ix_edge,:),:)'*w0;
        end

        function [f_out, x_interp] = interpolate(obj, x_target, f_in,mode)
            % interpolates input function to poistions over the sphere into
            % the direction x_target
            switch mode
                case 'linear'
                    [ixout,x_interp,w0] = find_edge(obj,x_target(1:2));
                    ix0 = obj.edges(ixout,:);
                    f_out = squeeze(sum(f_in(ix0,:,:).*w0,1)).';
                case 'nearest'
                    [~,ix] = min(sum(( obj.x0- x_target(1:2)'./norm(x_target(1:2)).*mean(obj.Phi0(:,2)) ).^2,2));
                    f_out = squeeze(f_in(ix,:,:))';
            end
        end

        function [dS] = get_dS(~,x0,edges)
            dS = zeros(size(x0,1),1);
            for n = 1 : size(x0,1)
                [r,c] = find(edges==n);
                dS(n) = sum(sqrt(sum((x0(edges(r,1),:)-x0(edges(r,2),:)).^2,2)))/2;
            end
        end

        function [dx_dir] = get_directive_dx(obj,kt,x0,edges)
            dx_dir = zeros(size(kt,1),1);
            for n = 1 : size(x0,1)
                [r,c] = find(edges==n);
                v = [ (x0(edges(r(1),3-c(1)),:) - x0(edges(r(1),c(1)),:))',...
                (x0(edges(r(2),3-c(2)),:) - x0(edges(r(2),c(2)),:))'];                    
                ix0 = find(v'*kt(n,:)'>0);
                if isempty(ix0)
                    dx_dir(n) = 1e10;
                else
                    dx_dir(n) = norm(v(:,ix0));
                end
            end
        end

        function [x0,n_in,dS, edges] = append_point(obj,xp)
            xp = xp(1:2)/norm(xp)*mean(obj.Phi0(:,2));
            x0 = obj.x0;
            edges = obj.edges;
            dS = obj.get_dS(x0,edges);

            [ix_edge,~,w0] = obj.find_edge(xp');
            x0 = [x0;xp];
            n_in = -x0./sqrt(sum(x0.^2,2));

            edge0 = edges(ix_edge,:);
            edges(ix_edge,:) = [];
            edges = [edges; [edge0(1),size(x0,1);edge0(2),size(x0,1)]  ];
%            dS = obj.get_dS(x0,edges);
            dS(end+1) = w0'*dS(edge0);
            dS(edge0) = (1-w0).*dS(edge0);
        end

        function [SHT_out] = SHT_on_circle(obj,in,Nmax)

        end
    
        function plot_design(obj)
            figure;
            scatter(obj.x0(:,1),obj.x0(:,2),10,'k','filled');
            axis equal tight;
            hold on;
            quiver(obj.x0(:,1),obj.x0(:,2),obj.n_out(:,1),obj.n_out(:,2))
            for n = 1 : size(obj.edges,1)
                line(obj.x0(obj.edges(n,:),1)',obj.x0(obj.edges(n,:),2)','Color','black')
            end
        end

    end
end

