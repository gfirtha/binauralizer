/*
 * field_snapshot_kernel.h  -  time-domain free-field snapshot (field_engine).
 *
 * Plain C interface, no MATLAB dependency (gateway: field_snapshot_mex.c).
 * Matrices are column-major.
 *
 *   D  [npix x nls] propagation delay in samples (R/c*fs)
 *   G  [npix x nls] 1/(4*pi*R)
 *   S  [nt x nls]   loudspeaker signals, S(i,n) = s_n(i/fs)
 *   t               evaluation time in samples (t0*fs)
 *   P  [npix]       sum_n G .* s_n(t0 - R/c), linear interpolation in time,
 *                   zero outside 1 <= i0 < nt
 *
 * Pixels are processed in blocks (OpenMP); within a block the loudspeaker
 * loop walks contiguous column segments of D and G.
 */
#ifndef FIELD_SNAPSHOT_KERNEL_H
#define FIELD_SNAPSHOT_KERNEL_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

void field_snapshot(const float *D, const float *G, const float *S,
                    size_t npix, size_t nls, size_t nt, double t, float *P);

#ifdef __cplusplus
}
#endif

#endif
