function app = run_binauralizer()
%RUN_BINAURALIZER  
addpath(genpath(fileparts(mfilename('fullpath'))));
% Reference copies of external work (e.g. External_software/Bubrik) carry
% functions with binauralizer names but other versions: keep them off the path
root = fileparts(mfilename('fullpath'));
bub = fullfile(root, 'External_software', 'Bubrik');
if isfolder(bub), rmpath(genpath(bub)); end
try
    app = audioapp.App();
catch ME
    fprintf(2, 'Failed to start app:\n%s\n', ME.getReport('extended'));
    rethrow(ME);
end
end 