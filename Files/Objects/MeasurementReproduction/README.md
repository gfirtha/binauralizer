# Measurement-based reproduction

This module analyses spherical-array room impulse responses (24-channel Firthenmike rigid sphere, r = 10 cm, SH order N = 3). The goal is to reproduce the measured rooms in the Binauralizer.

**Credit.** The analysis methods are from Márton Bubrik's thesis work, ported from his scripts. The original code is kept for reference in `External_software/Bubrik`. That folder is not on the app's path (see `run_binauralizer.m`), because its `lib/` holds other versions of functions that the binauralizer also has: `get_regularization_filters`, `get_rotation_mx`, `getSphH` and `spherical_design`.

## Contents

| File | What it does | Origin |
|---|---|---|
| `mbr_mic_array` | Array geometry from a unit-vector file, with front microphone #13 rotated to +x; SHT by pseudo-inverse | `spherical_design('File', ...)` |
| `mbr_encoder_filters` | Rigid-sphere radial equaliser (modal / DAS), crossover regularisation, low-pass | `get_beamforming_filters`, `get_regularization_filters` |
| `mbr_beamforming_map` | max-rE SH beamformer power map of a time window | `offl_conventional_window_beamforming.m` |
| `mbr_music_map` | SH-MUSIC: modal equalisation, whitening, noise-subspace pseudo-spectrum | `offline_MUSIC_windowed.m` (V9) |
| `mbr_herglotz_echoes` | Herglotz localisation function, Tikhonov regularisation, Gaussian mixture fit with model-order selection (Massé et al.) | `offl_herglotz_echo_detection.m` |
| `mbr_simulate_sma_rir` | Synthetic rigid-sphere RIR of plane-wave echoes (Massé IV-A, Eq. 16) | `masse_chapter4_simulation_v4.m` |
| `mbr_max_re_weights`, `mbr_window` | Helpers | |

The example data are in `Data/Measurements/firthenmike/`: `rir_front_trunc.mat` (measured) and `rir_mass4_sim.mat` (the thesis simulation).

The GUI is `audioapp.ui.MeasurementPanel`, the "Measurement-based reproduction" tab.

## Validation (2026-10-09)

The port was compared with the original pipelines on the same array geometry. The originals were run through their own callbacks.

| Pipeline | Result |
|---|---|
| Herglotz, two windows | Identical: echo count, directions and levels |
| SH-MUSIC | Identical: DOA and eigenvalue ratio |
| Beamforming map | Equal to 1e-13 (see the display fix below) |
| Ground truth | 4 echoes simulated with the exact geometry and Legendre polynomials (N = 3): all 4 detected, angle errors 0.6–6.3°, levels within 1.7 dB |

## Fixes relative to the thesis code

1. **Array rotation.**
   - **Problem:** the thesis `get_rotation_mx` had `a(3)^3` instead of `a(3)^2` and an unnormalised axis, so the matrix was not orthogonal. Microphone radii ranged from 0.094 to 0.100 m, and directions were off by 1.4° on average (at most 2.1°).
   - **Fix:** the port uses the exact Rodrigues rotation.
   - **Effect on measured data:** in a weakly structured window (0.5–4.5 ms of `rir_front_trunc`) the Herglotz model order changes from 3 to 2. Elsewhere the changes are about 2°.
2. **Legendre recursion in the simulation.**
   - **Problem:** P_n was computed with the coefficients of P_{n+1}, i.e. ((2n+1)x P_{n-1} − n P_{n-2})/(n+1), giving errors up to 0.17 for n ≤ 6. The simulated field was therefore not exactly the intended plane waves.
   - **Consequence:** the thesis' `rir_mass4_sim.mat` yields a 5th, spurious echo.
3. **Beamforming map display.** The thesis showed the power map as 20·log10(power), i.e. in doubled dB. The panel shows 10·log10.

Not ported:
- the `'EQ'` radial filter mode (it calls undefined functions);
- `realtime_MUSIC.m` (live input through `mcdmex_msvc`; planned as live analysis).

## Planned

- **Microphone layouts:** selectable layouts via `mbr_mic_array.from_file` plus GUI.
- **Parametric reproduction:** direct sound plus the detected reflections (time, direction, level) as virtual sources in the Binauralizer.
- **Late part / diffuse residual:** SH/Ambisonic binaural decoding (see `hoa_renderer`).
