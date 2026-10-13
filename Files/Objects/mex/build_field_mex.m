function build_field_mex()
%BUILD_FIELD_MEX  Compile field_snapshot_mex (OpenMP) next to this file: the
%gateway field_snapshot_mex.c and the MATLAB-independent kernel
%field_snapshot_kernel.c.
%   Needs a configured C compiler (mex -setup C); tested with MinGW64.
%   field_engine falls back to its MATLAB path when the MEX is missing.
here = fileparts(mfilename('fullpath'));
src  = fullfile(here, 'field_snapshot_mex.c');
mex('-O', '-outdir', here, src, fullfile(here, 'field_snapshot_kernel.c'), ...
    'CFLAGS=$CFLAGS -fopenmp -O3', 'LDFLAGS=$LDFLAGS -fopenmp');
end
