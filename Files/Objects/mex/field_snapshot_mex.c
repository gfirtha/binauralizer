/*
 * field_snapshot_mex.c  -  MATLAB gateway of field_snapshot_kernel
 * (time-domain free-field snapshot, field_engine)
 *
 *   p = field_snapshot_mex(D, G, S, t0fs)
 *
 *   D    [Npix x Nls] single : propagation delay in samples (R/c*fs)
 *   G    [Npix x Nls] single : 1/(4*pi*R)
 *   S    [Nt   x Nls] single : loudspeaker signals, S(i,n) = s_n(i/fs)
 *   t0fs double scalar       : evaluation time in samples (t0*fs)
 *   p    [Npix x 1]   single : sum_n G .* s_n(t0 - R/c), linear interpolation
 *                              in time, zero outside 1 <= i0 < Nt
 *
 * Build: build_field_mex.m (MinGW64 with -fopenmp).
 */
#include "mex.h"
#include "field_snapshot_kernel.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    if (nrhs != 4)
        mexErrMsgIdAndTxt("field_snapshot_mex:nrhs", "Usage: p = field_snapshot_mex(D, G, S, t0fs)");
    for (int k = 0; k < 3; ++k)
        if (!mxIsSingle(prhs[k]) || mxIsComplex(prhs[k]))
            mexErrMsgIdAndTxt("field_snapshot_mex:type", "D, G and S must be real single.");

    const mwSize Npix = mxGetM(prhs[0]);
    const mwSize Nls  = mxGetN(prhs[0]);
    const mwSize Nt   = mxGetM(prhs[2]);
    if (mxGetM(prhs[1]) != Npix || mxGetN(prhs[1]) != Nls || mxGetN(prhs[2]) != Nls)
        mexErrMsgIdAndTxt("field_snapshot_mex:size", "Inconsistent sizes of D, G and S.");

    plhs[0] = mxCreateNumericMatrix(Npix, 1, mxSINGLE_CLASS, mxREAL);
    field_snapshot((const float *) mxGetData(prhs[0]), (const float *) mxGetData(prhs[1]),
                   (const float *) mxGetData(prhs[2]), Npix, Nls, Nt, mxGetScalar(prhs[3]),
                   (float *) mxGetData(plhs[0]));
}
