# Avendano & Jot frequency-domain stereo upmix: reading list

This folder collects the literature on Avendano and Jot's STFT-domain upmix method. The method has two parts:

- **Ambience extraction** uses the inter-channel short-time coherence with a nonlinear (tanh) mask.
- **Source identification** uses the similarity and panning index, with panning-index windows for centre extraction, un-mix and re-panning.

The collection was assembled on 2026-10-08. All PDFs come from legal sources: an author-hosted copy (Columbia, D. Ellis), an author preprint (SigPort), and Google Patents. Each file was checked: it starts with `%PDF`, and its title matches the work.

**Recommended reading order:**
1. WASPAA 2003 paper (free, the cleanest statement of ψ and Ψ and of the window).
2. US7567845 or US7315624: the same text as the AES 2002 / JAES 2004 papers in patent form. It covers the coherence formulas, the tanh mask with its parameter table, the 2-to-3 and 2-to-5 upmix, and un-mix.
3. JAES 2004 (paywalled), if you have AES access.
4. For comparison: He/Tan/Gan 2014 (PCA and least-squares primary-ambient extraction).

---

## 1. ICASSP 2002 (ambience): paywalled

C. Avendano, J.-M. Jot, "Ambience extraction and synthesis from stereo signals for multi-channel audio up-mix," *Proc. IEEE ICASSP 2002*, Orlando, FL, May 2002, vol. II, pp. 1957-1960.
The page numbers come from memory and are NOT verified. Look it up on IEEE Xplore: https://ieeexplore.ieee.org (search the title). Patents cite the title variably, as "...Multichannel Audio Upmix" or "...Mixup".

- **Local file:** none (IEEE Xplore, paywalled). No legal free copy was found.
- **Content:** this is the ambience half of the method. The paper computes the inter-channel short-time coherence Φ(m,k) from recursively smoothed auto- and cross-spectra and passes it through a nonlinear mapping M[Φ]. It scales both channel STFTs by M[Φ] to get the ambience signals, then decorrelates them with all-pass filters to feed the surround channels. The same formulas appear verbatim in **US7567845** (see item 5), which you can read instead.

## 2. AES 22nd Int. Conf. 2002 (full upmix framework): paywalled

C. Avendano, J.-M. Jot, "Frequency Domain Techniques for Stereo to Multichannel Upmix," *Proc. AES 22nd Int. Conf. on Virtual, Synthetic and Entertainment Audio*, Espoo, Finland, June 2002, paper 000251, pp. 121-130. The page numbers are as cited in the WASPAA 2003 paper.

- **Link:** https://aes2.org/publications/elibrary-page/?id=11131
- **Local file:** none (AES E-Library, paywalled).
- **Content:** this paper first introduces the combined system: coherence-based ambience extraction, the similarity ψ and panning index, re-panning to L/C/R using the panning law and vector-based localization, and a 2-to-5 upmix. The WASPAA 2003 paper notes that this paper used the non-energy-preserving panning law (g_L = 1-α, g_R = α). The patent text in items 5a and 5c matches this description closely.

## 3. WASPAA 2003 (panning index, extraction window): downloaded

C. Avendano, "Frequency-domain source identification and manipulation in stereo mixes for enhancement, suppression and re-panning applications," *Proc. IEEE WASPAA 2003*, New Paltz, NY, Oct. 19-22, 2003, pp. 55-58.

- **Local file:** `2003_Avendano_WASPAA_source_identification_stereo_mixes.pdf` (author copy hosted at https://www.ee.columbia.edu/~dpwe/papers/Aven03-unmix.pdf)
- **Content:** this paper uses the energy-preserving panning model x1 = √(1-α²)·s, x2 = α·s, so ψ = 2α√(1-α²). It defines the following (verified from the PDF):
  - Similarity: ψ(m,k) = 2·|X1 X2*| / (|X1|² + |X2|²)  (eq. 2). It equals 1 for centre-panned sources and 0 for sources panned hard to one side. These are *instantaneous* STFT values, with no smoothing.
  - Partial similarities: ψ_i(m,k) = |X_i X_j*| / |X_i|², for i ≠ j  (eq. 3). Their difference is Δ = ψ1 - ψ2  (eq. 4), and Δ̂ = sign(Δ) ∈ {1, 0, -1}  (eq. 5).
  - Panning index: **Ψ(m,k) = [1 - ψ(m,k)]·Δ̂(m,k)**  (eq. 6), with Ψ ∈ [-1, 1] and Ψ = 0 at the centre.
  - Panning-index window: **Θ(m,k) = ν + (1-ν)·exp(-(Ψ(m,k)-Ψ0)² / (2ξ))**  (eq. 8). Here ν is a floor that avoids musical noise. The width is set by **ξ = -(Ψc-Ψ0)² / (2 log A)**  (eq. 9), where A is the rejection level (e.g. -60 dB) at the point Ψc. Ψc is the maximum panning-index error for a chosen SIR (signal-to-interference ratio). The example value is ξ = 0.006 for vocal suppression, corresponding to SIR = -60 dB.
  - Applications:
    - Suppression: Y_i = X_i(1-Θ).
    - Enhancement: Y_i = X_i(1+βΘ).
    - Re-panning: Y_i = X_i[1 + Θ(ρ_i-1)], with ρ_i = γ_i/α_i.
    - Resynthesis uses a least-squares inverse STFT (Griffin-Lim).
- **Sign-convention caveat:** see the note at the end of this file.

## 4. JAES 2004 (journal version): paywalled

C. Avendano, J.-M. Jot, "A Frequency-Domain Approach to Multichannel Upmix," *J. Audio Eng. Soc.*, vol. 52, no. 7/8, pp. 740-749, July/Aug. 2004. The page numbers come from memory and are NOT verified.

- **Link:** https://aes2.org/publications/elibrary-page/?id=13016
- **Local file:** none (AES E-Library, paywalled).
- **Content:** according to the abstract, this is the journal consolidation of items 1 and 2. A similarity measure estimates the panning of sources in the time-frequency plane, heuristic mapping functions extract the sources, and perceptual (localization-based) functions re-pan them, giving a 2-to-5 upmix. US8107631 cites it for the coherence-based soft ambience mask. Its exact parameter values could not be checked here. Treat the patent values below as the same system, but this is unverified.

## 5. Creative Technology patents (Google Patents PDFs, public)

All of the following are assigned to Creative Technology Ltd. Status is as shown on Google Patents on 2026-10-08.

| File | Patent | Title | Inventors | Filed / granted | Status (Google Patents) |
|---|---|---|---|---|---|
| `US7567845_Avendano_Jot_ambience_generation.pdf` | US7567845B1 | Ambience generation for stereo signals | Avendano, Jot | 2002-06-04 / 2009-07-28 | Expired - Lifetime (adj. exp. 2025-08-22) |
| `US7315624_Avendano_Jot_stream_segregation.pdf` | US7315624B2 | Stream segregation for stereo signals | Avendano, Jot | prio. 2002-06-04, appl. 2006-10-27 / 2008-01-01 | Expired - Lifetime (anticipated exp. 2022-06-04) |
| `US8019093_Avendano_Jot_stream_segregation.pdf` | US8019093B2 | Stream segregation for stereo signals | Avendano, Jot | prio. 2002-06-04, appl. 2007-10-31 / 2011-09-13 | Expired - Lifetime (adj. exp. 2024-11-05) |
| `US8280077_Avendano_Jot_stream_segregation.pdf` | US8280077B2 | Stream segregation for stereo signals | Avendano, Jot | prio. 2002-06-04, appl. 2007-10-31 / 2012-10-02 | Expired - Fee Related (adj. exp. 2024-08-24) |
| `US7412380_Avendano_Goodwin_Jot_ambience_extraction_modification.pdf` | US7412380B1 | Ambience extraction and modification for enhancement and upmix of audio signals | Avendano, Goodwin, Sridharan, Wolters, Jot | 2003-12-17 / 2008-08-12 | Expired - Lifetime (adj. exp. 2025-11-28) |
| `US8107631_Merimaa_Goodwin_Jot_correlation_ambience.pdf` | US8107631B2 | Correlation-based method for ambience extraction from two-channel audio signals | Merimaa, Goodwin, Jot (not Avendano) | prio. 2007-10-04 / 2012-01-31 | **ACTIVE** (adj. exp. 2030-12-01) |

**What each patent covers:**
- **US7567845 (5a)** and the US7315624 / US8019093 / US8280077 family (5b) share essentially the same description. That description is the JAES/AES method in full: coherence ambience extraction, the similarity ψ, the panning index Γ, 2-to-3 re-panning, un-mix, the energy histogram for identifying sources, 2-to-5 upmix, and a parameter table. They differ only in their claims:
  - US7567845 claims ambience extraction (a low inter-channel correlation measure, a nonlinear mapping, a forgetting factor).
  - US7315624 claims source separation by classifying time-frequency bins with similar panning coefficients.
  - US8019093 claims generating a centre signal from centre-panned bins.
  - US8280077 claims deriving a per-band localization position.
- **US7412380:** a user-controlled ratio of ambience to total signal per channel, for enhancement and upmix.
- **US8107631:** recursive correlation estimates with τ = 1/(f_c ln(1/λ)) and f_c = f_s/hop. It adds a bias compensation of the short-time correlation coefficient, |φ̂_LR| = max{0, 1 - (1-|φ_LR|)/λ} (eq. 44), and ambience masks derived under "equal levels" or "equal ratios" assumptions for the ambience. **This patent is still in force until about 2030.** Avoid copying its compensation and mask derivation in a product without a freedom-to-operate (FTO) review.

### Formulas verified in US7567845 (identical text in US7315624, cols. 5-12)

- Short-time statistics (1a-c): Φ_LL = Σ_n S_L S_L*, Φ_RR = Σ_n S_R S_R*, Φ_LR = Σ_n S_L S_R*.
- In practice they are computed recursively (3): **Φ_ij(m,k) = λ·Φ_ij(m-1,k) + (1-λ)·S_i(m,k)·S_j*(m,k)**.
- Inter-channel short-time coherence (2): **Φ(m,k) = |Φ_LR(m,k)| · [Φ_LL(m,k)·Φ_RR(m,k)]^(-1/2)**. It is near 1 where the direct path dominates, even for a source panned to one side, and near 0 in reverberant tails.
- Simplest ambience estimate: multiply both STFTs by 1 - Φ.
- General form (4a,b): **A_L = S_L·M[Φ], A_R = S_R·M[Φ]**, with
  **M[Φ] = 0.5·(μ_max - μ_min)·tanh{σπ(Φ_0 - Φ)} + 0.5·(μ_max + μ_min)**  (5).
  Low coherence passes and high coherence is attenuated. μ_max is normally 1, and μ_min > 0 avoids spectral-subtraction artifacts.
- **Defaults (parameter table, 44.1 kHz):**

  | Parameter | Value | Meaning |
  |---|---|---|
  | N | 1024 | STFT window |
  | K | 2048 | FFT size |
  | L | 256 | hop |
  | λ | 0.90 | forgetting factor |
  | σ | 8.00 | slope of M |
  | Φ_0 | 0.15 | breakpoint |
  | μ_min | 0.05 | floor |
  | Δ | 256 | rear-channel delay (samples) |
  | N_p | 15 | complex-conjugate pole pairs of the decorrelation all-pass filter G(z) |

  The residual direct path in the ambience is about -25 dB. Raising σ, or lowering Φ_0 and μ_min, makes the mask more aggressive.
- Panning model (6): s_L = Σ(1-α_i)s_i, s_R = Σα_i s_i. This is the linear, non-energy-preserving law, so α = 0.5 is the centre.
- Similarity (7): **Ψ(m,k) = 2|S_L S_R*| / (|S_L|² + |S_R|²)**, computed with *instantaneous* values ("instantaneous correlation"). For a single source, Ψ = 2(α-α²)/(α²+(1-α)²).
- Partial similarities (7a,b): **Ψ_L = |S_L S_R*| / |S_L|²** and **Ψ_R = |S_R S_L*| / |S_R|²**. Their difference is D = Ψ_L - Ψ_R (8). D' = +1 if D > 0, and -1 if D ≤ 0 (9).
- **Panning index (10): Γ(m,k) = [1 - Ψ(m,k)]·D'(m,k)**. Note that the patent names the similarity Ψ and the panning index Γ, while WASPAA uses ψ and Ψ.
- 2-to-3 re-panning: Γ_L = 1 for Γ < 0 (else 0), and Γ_R = 1 for Γ ≥ 0 (else 0). This gives the regions S_LL = S_L Γ_L, S_LR = S_R Γ_L, S_RL = S_L Γ_R, S_RR = S_R Γ_R. The outputs are s_LF = ISTFT{S_LL M_LF}, s_LC = ISTFT{S_LR M_LC}, s_RC = ISTFT{S_RL M_RC}, s_RF = ISTFT{S_RR M_RF}, with C = s_LC + s_RC. The modification functions M are the ratio of the desired to the input panning gains. They come from vector-based localization (velocity vectors below about 700 Hz, energy/intensity vectors above), with a loudspeaker angle of θ = 60° or 30° (patent Figs. 7-8, plotted only).
- Un-mix: a narrow window is centred on the panning index of the target. A hard selection (Θ = 1 where Γ = 0) is replaced by a smooth window. The example uses windows at Γ = 0 (speech), Γ = -0.8 and Γ = 0.27.
- Source identification: an energy histogram of Γ over [-1, 1] in steps of 0.01. Its peaks give the prominent sources.
- Front channels of the 2-to-5 system (Fig. 12): extract the centre from centre-panned bins with a nonlinear mapping of Γ, and send the residual to L and R.

## 6. Comparison: PCA / linear-estimation primary-ambient extraction (PAE)

**6a.** M. M. Goodwin, J.-M. Jot, "Primary-ambient signal decomposition and vector-based localization for spatial audio coding and enhancement," *Proc. IEEE ICASSP 2007*, Honolulu, HI, Apr. 2007, vol. I, pp. I-9 to I-12. The pages are as cited in patents.
- **Local file:** none (IEEE Xplore, paywalled). No legal free copy was found.
- **Content:** this paper does per-band PCA of the 2×2 STFT covariance. The primary component is the projection onto the principal eigenvector, and the ambience is the residual, followed by vector-based localization of the primary. Freely readable patent counterparts are on Google Patents (not downloaded): US8204237 "Adaptive primary-ambient decomposition", US9088855 "Vector-space methods for primary-ambient decomposition", and US8103005 "...using a complex similarity index".

**6b.** J. He, E.-L. Tan, W.-S. Gan, "Linear estimation based primary-ambient extraction for stereo audio signals," *IEEE/ACM Trans. Audio, Speech, Lang. Process.*, vol. 22, no. 2, pp. 505-517, Feb. 2014, DOI 10.1109/TASLP.2013.2297015.
- **Local file:** `2014_He_Tan_Gan_TASLP_linear_estimation_PAE_preprint.pdf` (author manuscript from SigPort, 37 pp.; IEEE "personal use" notice included)
- **Content:** this paper puts PCA and least-squares PAE into one linear-estimation framework, proposes error measures (leakage, distortion, interference), and compares PCA, LS and three LS variants. It is a good reference for what the coherence-mask method trades away.
- **For a broader review** (not downloaded, paywalled): J. He, *Spatial Audio Reproduction with Primary Ambient Extraction*, SpringerBriefs in Signal Processing, 2017, DOI 10.1007/978-981-10-1551-9.

---

## Note: sign convention of the panning index (important for implementation)

Both the WASPAA 2003 text and the patent text say that **positive D (or Δ) means the source is panned left**. The formulas as printed say the opposite.

Take ψ_L = |S_L S_R*|/|S_L|², which equals |S_R|/|S_L|. Then:

D = |S_R|/|S_L| - |S_L|/|S_R|

This is **negative** when |S_L| > |S_R|, so a left-panned source gives a negative D. The WASPAA Fig. 1b agrees with the formula: Ψ rises toward +1 as α, the gain of channel 2 (right), goes to 1. The patent's own numerical examples are mixed:
- the un-mix example labels a right-dominant trumpet (0.1/0.9) as Γ = -0.8 and a left-dominant guitar (0.7/0.3) as Γ = +0.27, which follows the text;
- its histogram text maps α = 0.1 to Γ = -0.8 and α = 0.7 to Γ = +0.275, which follows the formula.

**Recommendation:** define the sign explicitly in code, for example
Ψ = (1-ψ)·sign(|X_L|² - |X_R|²) for "positive means left",
and choose the window centres Ψ0 to match. The magnitude 1-ψ is the same either way.

## Not verified / open points

- The exact content of the ICASSP 2002, AES 2002 and JAES 2004 papers could not be read because they are paywalled. Whether the journal version uses smoothed Φ_ij inside ψ, or the same tanh parameters as the patent, is not confirmed.
- Page numbers for ICASSP 2002 and JAES 2004 come from memory.
- Patent expiry dates come from Google Patents' computed fields. They are not legal advice.
