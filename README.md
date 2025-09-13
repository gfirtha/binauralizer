# Binauralizer

Interactive MATLAB framework for sound field synthesis, spatial rendering, and binaural monitoring.  
Load audio or video, pick a renderer (WFS, DBAP, CTC, time-delay stereo, direct playback), visualize the sound field, and listen in real time via HRTFs.

---

## Features

- **Real-time audio engine**
  - Block-based processing with ring buffer
  - Audio/video loader with A/V sync
  - Output device picker and device info
  - Binauralization via SOFA HRTFs
  - Optional **Master EQ** (graphic + parametric + simple low-pass) with hard bypass

- **Renderers**
  - **WFS** (generalized 2.5D):
    - Referencing: *concentric/parallel*, *fixed point*, *fixed distance*
    - Per-element anti-aliasing (optional) and Tukey tapering
    - Exposes illuminated/tapered **zones** for visualization
  - **DBAP** (Distance-Based Amplitude Panning), updated to the current engine
  - **CTC** (Cross-Talk Cancellation):
    - Plant models: HRTF / point source / rigid sphere
    - Source models: HRTF / point source / rigid sphere
  - **Time-delay stereo** (pairwise delay filters)
  - **Direct playback** (no spatialization)

- **Visualization**
  - **Sound-field simulator**:
    - Harmonic , Impulse, **Monitor (live)** modes
    - Amplitude error and magnitude views (harmonic)
    - Transfer function (dB) at the receiver
  - **Renderer overlays** in the main axes:
    - Reference contour, amplitude distribution, tapering window
    - Illuminated/tapered **zones** (filled polygons)
  - Floating **zoom** controls anchored to the axes

---

## Requirements

- MATLAB (recent release recommended)
- **Audio Toolbox** and **DSP System Toolbox**
- **SOFA API for MATLAB** (for HRTFs), e.g. `SOFAload` (Included in the code)
- At least one SOFA HRTF file (e.g., `HRTFs/BuK_ED_corr.sofa`)
- Optional: OS video codecs supported by MATLAB if you load videos

The app detects missing toolboxes and will disable playback with a clear message.

---

## Installation

```matlab
% Clone the repository then in MATLAB:
addpath(genpath('<path-to>/binauralizer'));
savepath; % optional
```

---
## Citing Binauralizer

If you use this software in academic work, please cite the repository:

```bibtex
@misc{Binauralizer,
  title        = {Binauralizer: MATLAB framework for sound field synthesis, spatial rendering, and binaural monitoring},
  author       = {Firtha, Gergely},
  howpublished = {\url{https://github.com/gfirtha/binauralizer}},
  year         = {2025},
  note         = {Version retrieved from GitHub}
}
```

For Wave Field Synthesis theory and implementation details, also consider:

- G. Firtha, P. Fiala, F. Schultz, S. Spors, “Improved Referencing Schemes for 2.5D Wave Field Synthesis Driving Functions,” IEEE/ACM TASLP, 25(5):1117–1127, 2017.
- G. Firtha, “A generalized Wave Field Synthesis framework with application for moving virtual sources,” PhD dissertation, BME, 2019.

---
