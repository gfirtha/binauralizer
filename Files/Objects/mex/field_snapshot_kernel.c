/* field_snapshot_kernel.c  -  see field_snapshot_kernel.h. Same maths as the
 * MATLAB path of field_engine.snapshot. */
#include "field_snapshot_kernel.h"

#include <math.h>
#include <string.h>

#define BLOCK 1024

void field_snapshot(const float *D, const float *G, const float *S,
                    size_t Npix, size_t Nls, size_t Nt, double t, float *P)
{
    const long long nBlocks = (long long)((Npix + BLOCK - 1) / BLOCK);

    #pragma omp parallel for schedule(dynamic, 4)
    for (long long b = 0; b < nBlocks; ++b) {
        const size_t p0 = (size_t) b * BLOCK;
        const size_t p1 = (p0 + BLOCK < Npix) ? p0 + BLOCK : Npix;
        double acc[BLOCK];
        memset(acc, 0, sizeof(acc));
        for (size_t n = 0; n < Nls; ++n) {
            const float *Dn = D + n * Npix;
            const float *Gn = G + n * Npix;
            const float *Sn = S + n * Nt;            /* Sn[i-1] = S(i,n) */
            for (size_t p = p0; p < p1; ++p) {
                const double q = t - (double) Dn[p];
                const double fi = floor(q);
                if (fi < 1.0 || fi >= (double) Nt) continue;
                const size_t i0 = (size_t) fi;      /* 1-based sample index */
                const double a = q - fi;
                const double v = (1.0 - a) * Sn[i0 - 1] + a * Sn[i0];
                acc[p - p0] += (double) Gn[p] * v;
            }
        }
        for (size_t p = p0; p < p1; ++p) P[p] = (float) acc[p - p0];
    }
}
