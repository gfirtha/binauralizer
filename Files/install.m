addpath(genpath('Files'))
addpath(genpath('Data'))
addpath(genpath('temp'))
addpath(genpath('SoundSamples'))
addpath(genpath('HRTF_extrapolation'))
SOFAstart;
addpath(genpath('SourceFiles'));
addpath(genpath('ExternalSoftware'));
warning('off','all');

%addpath(pwd);                      % parent of +bz
%addpath(genpath(fullfile(pwd,'+bz')));
%rehash; clear classes
%which bz.util.make_input_src -all