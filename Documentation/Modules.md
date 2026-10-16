# Binauralizer modules and roadmap

The app has four top-level modules, one tab each in `audioapp.ui.MainWindow`.

| # | Module | GUI | Model | Status |
|---|---|---|---|---|
| 1 | Binauralizer | `MainWindow` (scene view and side tabs) | `sound_scene`, renderers | working |
| 2 | Source separation | `audioapp.ui.SeparationPanel` | `Files/Objects/Separation` | working, cleanup in progress |
| 3 | Measurement-based reproduction | `audioapp.ui.MeasurementPanel` | `Files/Objects/MeasurementReproduction` (`mbr_*`) | analysis (M. Bubrik's thesis work) + parametric reproduction (reflections + SDM) |
| 4 | Room builder | `audioapp.ui.RoomBuilderPanel` | `Files/Objects/RoomBuilder/room_model` | placeholder (only geometry import works) |

Methods that are not implemented yet raise `<class>:notImplemented`. Controls that are planned but not working are disabled in the GUI.

## 1. Source separation: next steps

- **Status of the methods:** `Files/Objects/Separation/separator_info.m` (also shown by the Info button). The retired methods and the reasons are in `legacy/README.md`.
- **Common streaming bugs:**
  - `process()` pads the output with zeros, so the actual latency depends on the block size. The STFT framing is copy-pasted, so it should be moved into one shared helper and fixed there.
  - `Pos` runs ahead of the audio by the latency.
- **Method-specific bugs:**
  - Mono input is counted twice in `stft` and `avendano`.
  - `avendano`: `sign(0)` on a silent channel.
  - `stft_objects`: the 5.1 C channel is added to a moving object.
- **Open problem:** no method tracks a source that circles around the listener in a 5.1 mix.

## 3. Measurement-based reproduction

- **Status:** spherical-array RIR analysis per time window, with three methods: max-rE beamforming, SH-MUSIC and Herglotz echo detection. The port was validated against the thesis code; see `Files/Objects/MeasurementReproduction/README.md`, which also lists the bugs fixed in the port.
- **Reproduction (working, experimental):** "Decompose" splits the RIR into reliable reflections and a residual (`mbr_decompose_rir`).
  - Each reflection becomes a point source with its own local, beamformed and time-focused IR.
  - The residual is checked for diffuseness (DirAC diffuseness and Avendano–Jot coherence) and rendered by SDM onto fixed directions.
  - "Send to Binauralizer" plays the loaded media through these objects (`mbr_room_separator`); "Bypass" returns to the media's own channels.
- **Later:** selectable microphone layouts and live analysis (`realtime_MUSIC`).

## 4. Room builder

- **Geometry import:** `.obj` / `.stl` (done). Next: surface grouping (coplanar merge, as in EntelFEM's preprocessing).
- **Materials:** per surface group, absorption and scattering, APD files.
- **Simulation:** image sources, ray tracing (the phased beam tracer), and low-frequency FEM / FDTD.
  - The EntelFEM engines are the reference implementation.
  - Their class names clash with binauralizer classes (`receiver`, `signal`), so they must be called through a bridge (separate MATLAB path or process), not put on the same path.
- **Export:** per-loudspeaker BRIRs / image-source sets for the room simulator renderer.

The Loudspeaker designer is a separate app; the binauralizer reads the loudspeakers it
exports (`designed_loudspeaker`, `Data/Loudspeakers`) and its Pre-EQ files (`speaker_eq`).
