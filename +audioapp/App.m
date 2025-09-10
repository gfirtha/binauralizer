classdef App < handle
    properties
        C   % audioapp.AppContext
        UI  % audioapp.ui.MainWindow
    end

    methods
        function obj = App()
            addpath(pwd);              % your paths, install(), etc. if needed
            install;
            obj.C = audioapp.AppContext();
            obj.UI = audioapp.ui.MainWindow(obj.C);   % builds figure, tabs, transport
        end
    end
end
