classdef hrtf_db
    % HRTF_DB  HRTF sets ("ear type") for the Binauralizer: the local SOFA
    % files (Data/HRTFs) and the online catalogue of sofacoustics.org
    % (~60 databases: ARI, CIPIC, HUTUBS, SADIE, MIT KEMAR, 3D3A, OlHeaD, ...;
    % dummy heads and individual listeners), downloaded on demand into
    % Data/HRTFs/<database>/.
    %
    %   L = hrtf_db.local()              struct array name, file
    %   D = hrtf_db.databases()          names of the online databases
    %   F = hrtf_db.files(db)            SOFA files of a database (names)
    %   file = hrtf_db.download(db, name)   the local copy (downloaded once)
    %   A free-field / eardrum HRTF (e.g. a head-and-torso simulator with an
    %   ear simulator) goes with the 'eardrum' headphone equalisation;
    %   diffuse-field equalised sets (DTF) with the 'diffuse' one.

    properties (Constant)
        Base = 'https://sofacoustics.org/data/database/'
    end

    methods (Static)
        function d = dataDir()
            d = fullfile(fileparts(fileparts(fileparts(fileparts(mfilename('fullpath'))))), 'Data', 'HRTFs');
        end

        function L = local()
            S = dir(fullfile(hrtf_db.dataDir(), '**', '*.sofa'));
            L = struct('name', {}, 'file', {});
            for k = 1:numel(S)
                rel = erase(fullfile(S(k).folder, S(k).name), [hrtf_db.dataDir() filesep]);
                L(end+1) = struct('name', rel, 'file', fullfile(S(k).folder, S(k).name)); %#ok<AGROW>
            end
        end

        function D = databases()
            D = hrtf_db.links(hrtf_db.Base, true);
        end

        function F = files(db)
            F = hrtf_db.links([hrtf_db.Base hrtf_db.enc(db) '/'], false);
        end

        function file = download(db, name)
            d = fullfile(hrtf_db.dataDir(), db);
            if ~isfolder(d), mkdir(d); end
            file = fullfile(d, name);
            if ~isfile(file)
                websave(file, [hrtf_db.Base hrtf_db.enc(db) '/' hrtf_db.enc(name)], weboptions('Timeout', 120));
            end
        end
    end

    methods (Static, Access = private)
        function N = links(url, dirs)
            % the directory listing's entries: sub-directories or .sofa files
            html = webread(url, weboptions('Timeout', 30, 'ContentType', 'text'));
            h = regexp(html, 'href="([^"?#]+)"', 'tokens'); h = cellfun(@(c) c{1}, h, 'UniformOutput', false);
            h = h(~startsWith(h, {'/', 'http', '..'}));
            if dirs, h = h(endsWith(h, '/')); h = erase(h, '/'); else, h = h(endsWith(lower(h), '.sofa')); end
            N = unique(cellfun(@(s) char(java.net.URLDecoder.decode(s, 'UTF-8')), h, 'UniformOutput', false), 'stable');
        end

        function e = enc(s)
            e = strrep(char(java.net.URLEncoder.encode(s, 'UTF-8')), '+', '%20');
        end
    end
end
