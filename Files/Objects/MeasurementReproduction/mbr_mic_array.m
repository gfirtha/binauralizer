classdef mbr_mic_array < handle
    %MBR_MIC_ARRAY  Spherical microphone array geometry for the
    %Measurement-based reproduction module (port of the array handling of
    %M. Bubrik's thesis code, spherical_design('File', ...)).
    %
    %   a = mbr_mic_array.firthenmike();        % 24-ch rigid sphere, r = 10 cm
    %   a = mbr_mic_array.from_file(file, R, frontIdx, name)
    %
    %   x0 [Nch x 3] microphone positions [m]; Phi0 [Nch x 3] = [azimuth,
    %   elevation, radius] (cart2sph); the array is rotated so that
    %   microphone frontIdx points to +x (front), as in the thesis code.
    %   Rotation: exact Rodrigues rotation. The thesis' get_rotation_mx had
    %   a(3)^3 instead of a(3)^2 and an unnormalised axis, i.e. a
    %   non-orthogonal matrix that distorted the microphone positions.
    %   Other layouts: from_file with any unit-vector list (planned: GUI).

    properties
        Name = ''
        R = 0.1             % sphere radius [m]
        x0 = zeros(0,3)     % microphone positions [m]
        Phi0 = zeros(0,3)   % [azimuth elevation radius]
        Nmax = 3            % highest SH order the layout supports
        Rigid = true        % rigid (scattering) sphere
    end

    methods (Static)
        function a = firthenmike(R)
            if nargin < 1, R = 0.1; end
            f = fullfile(fileparts(mfilename('fullpath')), 'arrays', 'firthenmike.dat');
            a = mbr_mic_array.from_file(f, R, 13, 'Firthenmike (24 ch, rigid sphere)');
            a.Nmax = 3;
        end

        function a = from_file(file, R, frontIdx, name)
            % file: whitespace separated unit vectors, 3 per microphone
            if nargin < 3, frontIdx = []; end
            if nargin < 4, [~, name] = fileparts(file); end
            fid = fopen(file, 'r');
            if fid < 0, error('mbr_mic_array:file', 'Cannot open %s.', file); end
            v = fscanf(fid, '%f'); fclose(fid);
            x0 = R*reshape(v, 3, [])';
            if ~isempty(frontIdx)
                Rm = mbr_mic_array.rotation(x0(frontIdx,:), [1 0 0]);
                x0 = x0*Rm';
            end
            a = mbr_mic_array();
            a.Name = name; a.R = R; a.x0 = x0;
            [az, el, r] = cart2sph(x0(:,1), x0(:,2), x0(:,3));
            a.Phi0 = [az el r];
            a.Nmax = floor(sqrt(size(x0,1))) - 1;
        end

        function Rm = rotation(u, v)
            % rotation matrix with Rm*u = v (Rodrigues, exact)
            u = u(:)/norm(u); v = v(:)/norm(v);
            a = cross(u, v); c = dot(u, v);
            if norm(a) < 1e-12
                if c > 0, Rm = eye(3); return; end
                p = null(u'); p = p(:,1);              % 180 deg about any normal axis
                Rm = 2*(p*p') - eye(3); return;
            end
            K = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
            Rm = eye(3) + K + K*K/(1 + c);
        end
    end

    methods
        function n = nCh(obj), n = size(obj.x0, 1); end

        function Y = sh_matrix(obj, N, type)
            % [Nch x (N+1)^2] SH matrix at the microphone directions (ACN)
            if nargin < 3, type = 'real'; end
            Y = getSpherHarmMx(pi/2 - obj.Phi0(:,2), obj.Phi0(:,1), N, type);
        end

        function psi = sht(obj, p, N, type)
            % least-squares SHT of microphone signals p [T x Nch] -> [T x (N+1)^2]
            if nargin < 4, type = 'real'; end
            psi = p*pinv(obj.sh_matrix(N, type)).';
        end
    end
end
