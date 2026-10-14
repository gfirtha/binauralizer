function root = app_paths()
%APP_PATHS  The MATLAB path of the binauralizer tree, shared by its apps
%(run_binauralizer, run_loudspeaker_designer): the whole tree, minus the
%reference copies of external work that carry functions with our names but
%other versions (Bubrik; Loudspeaker_model_bem shadows get_default_layout;
%nihu: only the built mex in nihu_bem is used, its src/matlab/install.m
%shadows Files/install.m). Session only: never savepath (a shadowed
%install.m once rewrote the global pathdef).
%
%   root = app_paths()   root: the binauralizer folder

root = fileparts(mfilename('fullpath'));
addpath(genpath(root));
for ext = {'Bubrik', 'Loudspeaker_model_bem', 'nihu'}
    d = fullfile(root, 'External_software', ext{1});
    if isfolder(d), rmpath(genpath(d)); end
end
end
