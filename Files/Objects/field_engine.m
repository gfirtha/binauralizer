classdef field_engine < handle
    % FIELD_ENGINE  Free-field pressure of a 2-D loudspeaker layout on a grid.
    %
    % Geometry (delay in samples and 1/(4*pi*R) for every pixel/loudspeaker
    % pair) is cached and reused until the grid or the loudspeakers move, so
    %   - a harmonic field is one complex matrix-vector product,
    %   - a time-domain snapshot is one vectorised gather,
    % instead of a per-loudspeaker interp2 for every frame.
    %
    %   E = field_engine(c, fs); E.set_grid(x, y); E.prepare(pos, ori);
    %   Z = E.harmonic(A, f0, dirs)   % complex, s_n(t) = Re{A_n e^{jwt}}
    %   p = E.snapshot(S, t0)         % S(i,n) = s_n(i/fs), p(x) at time t0

    properties
        c
        fs
        max_bytes = 4e8     % cache budget; larger grids are evaluated in chunks
        use_mex = true      % snapshot via field_snapshot_mex when it is compiled
    end

    properties (SetAccess = private)
        x                   % grid vectors
        y
        D = []              % [Npix x Nls] single: R/c*fs (propagation delay, samples)
        G = []              % [Npix x Nls] single: 1/(4*pi*R)
        cosTheta = []       % [Npix x Nls] single: angle to loudspeaker axis (lazy)
        pos = []
        ori = []
    end

    methods
        function obj = field_engine(c, fs)
            obj.c = c;
            obj.fs = fs;
        end

        function set_grid(obj, x, y)
            if isequal(x(:), obj.x) && isequal(y(:), obj.y), return; end
            obj.x = x(:); obj.y = y(:);
            obj.clear_cache();
        end

        function prepare(obj, pos, ori)
            % pos, ori: [Nls x >=2] loudspeaker positions / axis directions
            pos = pos(:,1:2); ori = ori(:,1:2);
            if isequal(pos, obj.pos) && isequal(ori, obj.ori) && ~isempty(obj.D), return; end
            obj.pos = pos; obj.ori = ori;
            obj.clear_cache();
            if obj.cache_fits()
                [obj.D, obj.G] = obj.geometry(1:numel(obj.y)*numel(obj.x));
            end
        end

        function Z = harmonic(obj, A, f0, dirs)
            % Z(x) = sum_n A_n D_n(theta,f0) exp(-jkR)/(4*pi*R), [Ny x Nx] complex
            if nargin < 4, dirs = {}; end
            w = 2*pi*f0/obj.fs;                 % phase per sample of delay
            Z = zeros(numel(obj.y)*numel(obj.x), 1, 'single');
            Z = complex(Z, Z);
            C = obj.row_chunks();
            for k = 1:numel(C)
                rows = C{k};
                [D_, G_, cT] = obj.chunk(rows, ~isempty(dirs));
                E = G_ .* exp(-1i*w*D_);
                if ~isempty(dirs)
                    E = E .* obj.directivity_at(dirs, cT, f0);
                end
                Z(rows) = E * single(A(:));
            end
            Z = reshape(double(Z), numel(obj.y), numel(obj.x));
        end

        function p = snapshot(obj, S, t0)
            % Time-domain field at time t0 [s]; S(i,n) is loudspeaker n at i/fs
            % (linear interpolation in time, zero outside the signal).
            [Nt, Nls] = size(S);
            S = single(S);
            off = single((0:Nls-1)*Nt);
            p = zeros(numel(obj.y)*numel(obj.x), 1, 'single');
            useMex = obj.use_mex && exist('field_snapshot_mex', 'file') == 3;
            C = obj.row_chunks();
            for k = 1:numel(C)
                rows = C{k};
                [D_, G_] = obj.chunk(rows, false);
                if useMex   % same maths, compiled + OpenMP (mex/build_field_mex.m)
                    p(rows) = field_snapshot_mex(D_, G_, S, t0*obj.fs);
                    continue;
                end
                q = single(t0*obj.fs) - D_;
                i0 = floor(q);
                a = q - i0;
                ok = i0 >= 1 & i0 < Nt;
                idx = i0 + off;
                idx(~ok) = 1;
                v = (1-a).*S(idx) + a.*S(idx+1);
                v(~ok) = 0;
                p(rows) = sum(G_ .* v, 2);
            end
            p = reshape(double(p), numel(obj.y), numel(obj.x));
        end
    end

    methods (Access = private)
        function clear_cache(obj)
            obj.D = []; obj.G = []; obj.cosTheta = [];
        end

        function tf = cache_fits(obj)
            tf = 3*4*numel(obj.x)*numel(obj.y)*size(obj.pos,1) <= obj.max_bytes;
        end

        function C = row_chunks(obj)
            % pixel index blocks: one block when the geometry is cached,
            % otherwise blocks that keep the temporaries within max_bytes
            Npix = numel(obj.x)*numel(obj.y);
            if ~isempty(obj.D)
                step = Npix;
            else
                step = max(1, floor(obj.max_bytes/(6*4*size(obj.pos,1))));
            end
            n = ceil(Npix/step);
            C = cell(1, n);
            for k = 1:n, C{k} = ((k-1)*step+1 : min(k*step, Npix))'; end
        end

        function [D_, G_, cT] = chunk(obj, rows, needAngle)
            cT = [];
            if ~isempty(obj.D)
                D_ = obj.D; G_ = obj.G;
                if needAngle
                    if isempty(obj.cosTheta), obj.cosTheta = obj.angle_cos(1:numel(obj.D(:,1))); end
                    cT = obj.cosTheta;
                end
            else
                [D_, G_] = obj.geometry(rows);
                if needAngle, cT = obj.angle_cos(rows); end
            end
        end

        function [D_, G_] = geometry(obj, rows)
            [X, Y] = meshgrid(obj.x, obj.y);
            dx = X(rows(:)) - obj.pos(:,1)';
            dy = Y(rows(:)) - obj.pos(:,2)';
            R = sqrt(dx.^2 + dy.^2);
            R = max(R, 1e-3);                       % keep the 1/R singularity finite
            D_ = single(R*(obj.fs/obj.c));
            G_ = single(1./(4*pi*R));
        end

        function cT = angle_cos(obj, rows)
            [X, Y] = meshgrid(obj.x, obj.y);
            dx = X(rows(:)) - obj.pos(:,1)';
            dy = Y(rows(:)) - obj.pos(:,2)';
            on = obj.ori ./ max(sqrt(sum(obj.ori.^2, 2)), eps);
            cT = single((dx.*on(:,1)' + dy.*on(:,2)') ./ max(sqrt(dx.^2 + dy.^2), 1e-3));
        end

        function Dm = directivity_at(~, dirs, cT, f0)
            % directivity of every loudspeaker at f0, tabulated over angle
            th = linspace(0, pi, 361);
            Dm = complex(zeros(size(cT), 'single'));
            for n = 1:numel(dirs)
                tab = dirs{n}.evaluate(th(:), f0);         % [361 x 1]
                Dm(:,n) = single(interp1(th, tab, real(acos(max(-1, min(1, double(cT(:,n))))))));
            end
        end
    end
end
