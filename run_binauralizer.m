function app = run_binauralizer()
%RUN_BINAURALIZER  
addpath(genpath(fileparts(mfilename('fullpath'))));
% Reference copies of external work (e.g. External_software/Bubrik) carry
% functions with binauralizer names but other versions: keep them off the path
root = fileparts(mfilename('fullpath'));
for ext = {'Bubrik', 'Loudspeaker_model_bem', 'nihu'}        % (data only: its functions shadow ours, e.g. get_default_layout;
                                                              %  nihu: the NiHu sources - only the built mex in nihu_bem is used,
                                                             %  and its src/matlab/install.m shadows Files/install.m)
    d = fullfile(root, 'External_software', ext{1});
    if isfolder(d), rmpath(genpath(d)); end
end
try
    a = audioapp.App();
catch ME
    fprintf(2, 'Failed to start app:\n%s\n', ME.getReport('extended'));
    rethrow(ME);
end
setappdata(groot, 'binauralizerApp', a);    % kept alive without an output variable
if nargout > 0, app = a; end                % 'run_binauralizer' alone: no "ans = App" print
end 