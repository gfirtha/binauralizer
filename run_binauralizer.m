function app = run_binauralizer()
%RUN_BINAURALIZER  
addpath(genpath(fileparts(mfilename('fullpath'))));
try
    app = audioapp.App();
catch ME
    fprintf(2, 'Failed to start app:\n%s\n', ME.getReport('extended'));
    rethrow(ME);
end
end
