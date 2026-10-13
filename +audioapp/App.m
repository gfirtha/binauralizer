classdef App < handle
    properties
        C   % audioapp.AppContext
        UI  % audioapp.ui.MainWindow
    end

    methods
        function obj = App()
            addpath(pwd);              % your paths, install(), etc. if needed
            % the app's own Files/install.m by its full path: another
            % install.m on the path (e.g. NiHu's, which also calls savepath
            % and opens its demo in the help browser) must not run instead.
            % Its content runs here (not run(), which changes into Files/
            % and breaks its relative paths).
            eval(fileread(fullfile(fileparts(fileparts(mfilename('fullpath'))), 'Files', 'install.m')));
            % Object-based (Atmos) input runs Cavern in a helper process
            try
                cavern_reader.serverPath();
            catch err
                fprintf(2, 'Cavern not available (object-based Atmos input disabled):\n  %s\n', ...
                    strtok(err.message, newline));
            end
            obj.C = audioapp.AppContext();
            obj.UI = audioapp.ui.MainWindow(obj.C);   % builds figure, tabs, transport
        end
    end
end
