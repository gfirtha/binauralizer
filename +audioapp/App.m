classdef App < handle
    properties
        C   % audioapp.AppContext
        UI  % audioapp.ui.MainWindow
    end

    methods
        function obj = App()
            addpath(pwd);              % your paths, install(), etc. if needed
            install;
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
