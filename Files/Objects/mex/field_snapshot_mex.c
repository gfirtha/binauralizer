/*
 * field_snapshot_mex.c  -  time-domain free-field snapshot (field_engine)
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
 * Same maths as the MATLAB path of field_engine.snapshot. Pixels are
 * processed in blocks (OpenMP); within a block the loudspeaker loop walks
 * contiguous column segments of D and G.
 *
 * Build: build_field_mex.m (MinGW64 with -fopenmp).
 */
#include "mex.h"
#include <math.h>
#include <string.h>

#define BLOCK 1024

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

    const float *D = (const float *) mxGetData(prhs[0]);
    const float *G = (const float *) mxGetData(prhs[1]);
    const float *S = (const float *) mxGetData(prhs[2]);
    const double t = mxGetScalar(prhs[3]);

    plhs[0] = mxCreateNumericMatrix(Npix, 1, mxSINGLE_CLASS, mxREAL);
    float *P = (float *) mxGetData(plhs[0]);

    const long long nBlocks = (long long)((Npix + BLOCK - 1) / BLOCK);

    #pragma omp parallel for schedule(dynamic, 4)
    for (long long b = 0; b < nBlocks; ++b) {
        const mwSize p0 = (mwSize) b * BLOCK;
        const mwSize p1 = (p0 + BLOCK < Npix) ? p0 + BLOCK : Npix;
        double acc[BLOCK];
        memset(acc, 0, sizeof(acc));
        for (mwSize n = 0; n < Nls; ++n) {
            const float *Dn = D + n * Npix;
            const float *Gn = G + n * Npix;
            const float *Sn = S + n * Nt;            /* Sn[i-1] = S(i,n) */
            for (mwSize p = p0; p < p1; ++p) {
                const double q = t - (double) Dn[p];
                const double fi = floor(q);
                if (fi < 1.0 || fi >= (double) Nt) continue;
                const mwSize i0 = (mwSize) fi;      /* 1-based sample index */
                const double a = q - fi;
                const double v = (1.0 - a) * Sn[i0 - 1] + a * Sn[i0];
                acc[p - p0] += (double) Gn[p] * v;
            }
        }
        for (mwSize p = p0; p < p1; ++p) P[p] = (float) acc[p - p0];
    }
}
