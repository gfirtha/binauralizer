# Legacy separators

These separators have been retired from the Source separation tab. They are kept here, still on the path, so that `make_separator` can construct them for old setups and comparisons.

The decision is based on a code review and a synthetic benchmark run on 2026-10-09. The benchmark used three scenes:
- static panned stereo sources;
- a fly-by over a centre vocal;
- a 5.1 scene with a vocal in C, a phantom between L and Ls, and a circling source.

| Method | Class | Why retired |
|---|---|---|
| `matrix` | `matrix_separator` | A broadband passive Lt/Rt decoder. The per-band `ltrt_decoder` is better on every scene. Output power is 0.5. |
| `avendano_ring` | `avendano_ring_separator` | A 5.1 four-pair bank that drops the C channel entirely, so dialogue is lost (centre vocal r² = 0.00). |
| `avendano_5ch` | `avendano_5ch_separator` | A 5.1 five-pair bank. Each channel is split by 1/√2 and then re-summed coherently, so the direct part comes out +3 dB (output power 1.6). |
| `phantom_objects` | `phantom_object_separator` | Broken, for three reasons: the overlap-add gain gives output power ×3.5; synthesis uses the first dictionary template `G(1,:)` instead of the matched gains, so the 5.1 C channel is lost; and the bed is mono. |
| `phantom_traj` | `phantom_trajectory_separator` | A subset of `phantom_objects` (one eigenvector, no ICA). It has the same defects, and it also subtracts the slot signals twice. Output power ×4. |

The pair-bank code stays inside `avendano_separator.m` (`EnablePairBank`, `PairBankMode`). The two subclasses above only switch it on.

Descriptions of all methods are in `separator_info.m`.
