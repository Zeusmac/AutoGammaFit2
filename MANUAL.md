# AutoGammaFit 2.0 — User Manual

## Contents
1. [Overview](#1-overview)
2. [Quick Start](#2-quick-start)
3. [Keyboard Shortcuts](#3-keyboard-shortcuts)
4. [AutoFit Tab](#4-autofit-tab)
5. [Manual Fit Tab](#5-manual-fit-tab)
6. [FWHM Tab](#6-fwhm-tab)
7. [Source Tab](#7-source-tab)
8. [Fit Results Tab](#8-fit-results-tab)
9. [Isotopes Tab](#9-isotopes-tab)
10. [Peak Table Tab](#10-peak-table-tab)
11. [Decay Tab](#11-decay-tab)
    - [11.1 Cuts Sub-tab](#111-cuts-sub-tab)
    - [11.2 Fitter Sub-tab](#112-fitter-sub-tab)
    - [11.3 Decay Models](#113-decay-models)
    - [11.4 Background Models](#114-background-models)
    - [11.5 Auto Model and T½ Seeding](#115-auto-model-and-t-seeding)
    - [11.6 Peak Counts vs Time](#116-peak-counts-vs-time)
    - [11.7 Fit Results and Peak Population](#117-fit-results-and-peak-population)
12. [Fitting Models and Mathematics](#12-fitting-models-and-mathematics)
    - [12.1 Standard N-Gaussian + Linear Background](#121-standard-n-gaussian--linear-background)
    - [12.2 Quadratic Background](#122-quadratic-background)
    - [12.3 Compton Step Model](#123-compton-step-model)
    - [12.4 Tied-Width Constraint](#124-tied-width-constraint)
    - [12.5 Same-Sigma Constraint](#125-same-sigma-constraint)
    - [12.6 Peak Area Calculation](#126-peak-area-calculation)
    - [12.7 Signal-to-Noise Ratio (post-fit)](#127-signal-to-noise-ratio-post-fit)
    - [12.8 Detector Resolution Model](#128-detector-resolution-model)
    - [12.9 Efficiency Model](#129-efficiency-model)
    - [12.10 NNDC Uncertainty Notation](#1210-nndc-uncertainty-notation)
    - [12.11 Chi-squared and Goodness of Fit](#1211-chi-squared-and-goodness-of-fit)
    - [12.12 Peak Significance](#1212-peak-significance)
    - [12.13 Decay Chain Mathematics](#1213-decay-chain-mathematics)
    - [**12.14 Complete Fit Statistics Reference**](#1214-complete-fit-statistics-reference)
13. [Algorithms](#13-algorithms)
14. [Cache System](#14-cache-system)
15. [References](#15-references)

---

## 1. Overview

AutoGammaFit is an interactive ROOT-based GUI for gamma-ray spectroscopy analysis.
It automates the full analysis chain — background removal, peak finding, multi-peak
fitting, isotope identification, efficiency correction, and decay curve analysis —
while giving the analyst fine-grained control at every step through a manual fitting
interface.

**Key capabilities:**
- Automated N-Gaussian fitting with AIC model selection (AutoFit)
- Manual fitting with MIGRAD, full covariance, residual display
- Multiple background models: linear, flat, quadratic, BG-anchor seeded
- Compton step correction (Erfc shelf model)
- Tied-width doublet fitting constrained to the detector resolution model
- S/N pre-screening to suppress spurious TSpectrum candidates
- Efficiency correction via log-polynomial model
- Channel→keV polynomial calibration
- Decay curve fitting (Bateman equations, A→B→C chains)
- Isotope identification, classification, and labeling
- Export to gnuScope binary format (.spe, .sqr)

**Launch:**
```bash
cd /path/to/AutoGammaFit && bin/gamma_gui
```

---

## 2. Quick Start

1. **Open ROOT file** — AutoFit tab → "Open ROOT File…"
2. **Open Isotope DB** — "Open Isotope DB…" (loads `Isotope_energys.txt` by default)
3. **Select histogram** from the list and classify it as "Gamma Spectrum"
4. **Run AutoFit (selected)** — fits TSpectrum-found peaks with N-Gaussian + linear BG
5. View results in the **Fit Results** tab; navigate peaks with Prev/Next
6. Assign isotope labels in the **Isotopes** tab

---

## 3. Keyboard Shortcuts

All shortcuts are global — they fire regardless of which panel or text field currently
has keyboard focus.

| Default | Action |
|---------|--------|
| **Ctrl+P** | Toggle **Choose Peaks** mode (place / remove seed markers on the canvas) |
| **Ctrl+F** | **Run Fit** — execute Manual Fit on the current peak seeds |
| **Ctrl+R** | **Run AutoFit (selected)** — fit the selected histogram |
| **Ctrl+L** | **Load Cache (selected)** — overlay stored fits on the canvas |
| **Ctrl+K** | **Clear Peaks** — remove all manual peak seed markers |
| **Ctrl+O** | **Open ROOT File** |
| **Ctrl+←** | Navigate to the **previous** fitted peak group |
| **Ctrl+→** | Navigate to the **next** fitted peak group |
| **Ctrl+B** | Toggle **Show BG curve** overlay |
| **Ctrl+E** | Toggle **Error Bars** on the histogram |
| **Ctrl+I** | Toggle **Isotope label** overlays on fitted peaks |
| **Esc** | Exit **Choose Peaks** mode (does nothing if already inactive) |

### Changing shortcuts

Click **Shortcuts…** (button row above the log panel) to open the shortcut editor.
Each row shows the action name and the current binding.

- **Rebind** — click the binding button, then press the desired key combination.
  Ctrl, Shift, and Alt modifiers are captured automatically.
- **Clear** — removes the binding for that action (shortcut is disabled).
- **Reset Defaults** — restores all bindings to the factory defaults shown above.
- **Save** — writes the current bindings to `.agf_shortcuts.conf` in the working
  directory and registers them immediately. Changes take effect without restarting.
- **Close** — closes the dialog; unsaved changes are discarded.

---

## 4. AutoFit Tab

AutoFit runs the full automated pipeline: background subtraction → peak finding →
peak grouping → adaptive multi-peak fitting → cache storage. It operates on one or
all Gamma Spectrum histograms in the loaded ROOT file.

### 3.1 Fit Options

| Option | Description |
|--------|-------------|
| Use Cached Seeds | Warm-start MIGRAD from the best previously converged parameters. Faster convergence and more consistent results on re-runs. Uncheck to always start from scratch. |
| Log-likelihood (L) | Minimise the Poisson log-likelihood instead of chi². Preferred for low-count spectra (any bin with < ~10 counts) where chi² is biased. |
| IMPROVE (M) | Run IMPROVE after MIGRAD to search for a better global minimum. Slower but helpful when the chi² surface has multiple local minima. |

### 3.2 Background Subtraction

The Compton continuum is estimated and subtracted before peak finding using
TSpectrum's iterative averaging (SNIP) algorithm [[Morhac et al. 1997](#ref-8)]. The subtracted
histogram is used only for peak finding and seeding — the actual Gaussian fit is
performed on the original histogram (background included as B₀ + B₁·x in the model).

- **Iterations**: number of averaging passes. Default 14. Increase to 20–30 for
  spectra with very broad features or overlapping Compton edges.
- **Preview / Apply to All**: verify the background shape on the canvas before
  committing to a full batch run.

### 3.3 Peak Finding

TSpectrum::Search scans the background-subtracted spectrum for local maxima
[[Morhac et al. 1997](#ref-8)].

- **Sigma (bins)**: the expected peak width in bins. Must roughly match the actual
  detector FWHM expressed in bins. Too low → misses broad peaks; too high → merges
  close doublets.
- **Threshold**: minimum peak height relative to the tallest peak in the spectrum
  (range 0–1). Default 0.02 means only peaks ≥ 2% of the maximum are found.
  Lower to find weak transitions; raise to suppress noise hits.

### 3.4 S/N Pre-screen

**What it does:** filters out TSpectrum candidates that are too weak to fit reliably,
*before* MIGRAD is called. This prevents wasted fit attempts on statistical
fluctuations and suppresses spurious cache entries in noisy regions of the spectrum.

**Why it matters:** TSpectrum is sensitive to noise, especially at low threshold
settings. Without a pre-screen, every noise bump gets a MIGRAD call, often producing
fits with parameters stuck at their bounds and chi²/ndf >> 1 that pollute the cache.

**How it works:** For each peak group, the algorithm sums the counts in a ±2σ signal
window on the background-subtracted working histogram, then estimates the noise from
the absolute counts in ±(4–7)σ sidebands:

```
signal = Σ max(h[b], 0)   for b in [E−2σ, E+2σ]

noise  = √( mean(|h[b]|) in sidebands  ×  n_peak_bins )

SN = signal / noise
```

If the best SN across all peaks in the group is below `Min S/N`, the entire group
is skipped and logged. Real gamma lines from a statistically significant source
typically yield SN > 10–100. Statistical noise bumps yield SN ≈ 1–2.

**Recommended settings:**
- `0` (default) — disabled; all groups are fitted. Use when you want to fit
  everything and filter manually in the Fit Results tab.
- `3` — suppresses fluctuations at 3σ significance or less. Good starting point
  for moderately noisy spectra.
- `5` — conservative; only clearly visible peaks are fitted.

*Technique inspired by the `GetSumFF` pre-screen in gnuScope [[Pavan & Tabor, FSU](#ref-14)].
Signal significance criterion: see [Bevington & Robinson (2003)](#ref-7), Ch. 4.*

### 3.5 Efficiency Correction

Stores a log-polynomial detector efficiency model. Once set, the Manual Fit peak
statistics panel reports an efficiency-corrected peak area alongside the raw area.

**Model** (see Section 12.8 for derivation):
```
ln(ε(E)) = a − b·ln(E) + c·(ln(E))² − d/E²
```

Enter the four parameters `a, b, c, d` determined from a source calibration
measurement (Source tab → Efficiency vs Energy plot), then click **Apply**. The
corrected area will appear in peak statistics after the next manual fit.

*Reference: [Knoll (2000)](#ref-6), Eq. 12.8; log-polynomial form also used in gnuScope
efficiency module [[Pavan & Tabor, FSU](#ref-14)]*

### 3.6 Rebin Histogram

Merges N adjacent bins into one, improving statistics at the cost of energy
resolution. Useful for weak sources or coarse ADC digitisation.

Poisson statistics are correctly preserved: `Sumw2()` is called before `Rebin(N)`
so uncertainties add in quadrature:
```
σ_merged = √(σ₁² + σ₂² + … + σ_N²) = √(c₁ + c₂ + … + c_N)
```
(the last equality holds for Poisson-distributed counts where σᵢ² = cᵢ).

The rebin factor is stored per-histogram in the cache metadata file so it persists
across sessions.

*Reference: Poisson error propagation — see [Bevington & Robinson (2003)](#ref-7), Ch. 3*

### 3.7 Custom Projections / Background Histogram Subtraction

- **Custom Projection**: slice a TH2 along X or Y with a user-defined cut range on
  the complementary axis. Creates a virtual 1D histogram in the session.
- **Background Histogram Subtraction**: `result = source − scale × background`.
  Sumw2() is applied so errors propagate as σ_result² = σ_source² + scale²·σ_bg².

### 3.8 Debug Sections

Toggle per-module diagnostic output to the log strip at the bottom of the window:

| Module | What it logs |
|--------|-------------|
| FITTER | AdaptiveFitter: seed values, retry stages, AIC scores |
| GROUPER | PeakGrouper: merge/split decisions, group boundaries |
| TRACKER | PeakTracker: FWHM fit updates, sigma bounds |
| DB | FitDatabase: cache load/save/seed operations |
| GAMMADB | GammaDB: isotope loading, match results |
| RESMODEL | ResolutionModel: parameter updates |
| FILEIO | ROOT file writes |
| PEAKFITTER | Per-peak fit results, SNR, matched isotopes |

---

## 5. Manual Fit Tab

Manual fitting gives full analyst control over peak placement, fit range, background
model, and parameter bounds. It uses MINUIT2 MIGRAD for precise parameter estimation
with a full covariance matrix, parameter uncertainty display, and optional MINOS
asymmetric errors.

### 4.1 Workflow

1. Load a histogram with the **Load** button (or select from the dropdown)
2. Choose a **View** mode: Raw, Background Subtracted, or BG Sub + AutoFit peaks
3. Click on the spectrum canvas to place peak seed markers (cyan crosshairs)
4. Set the **Fit Range** with "Set from Canvas" (two clicks) or type Lo/Hi
5. Set the background model and seeds (Section 4.2)
6. Click **Run Fit** — MIGRAD minimises chi² or log-likelihood
7. Inspect the **Peak Statistics** panel and residuals
8. Assign a label and classification, then **Accept & Save** to write to cache

### 4.2 Background Fitting

The background represents the Compton continuum and detector noise floor that
underlies the gamma peaks. Accurate background modelling directly affects the
fitted peak areas and therefore all downstream quantities (efficiency, activity).

#### Simple linear background (default)

```
BG(x) = B₀ + B₁·x
```

Seed `bg0` and `bg1` by typing in the Fit Parameters popup, or use
**Fit Background**: select a region with no gamma lines, click the two endpoints
on the canvas, then "Fit Background" estimates a straight line through the
off-peak region.

#### BG Anchor Regions

Specify two distinct off-peak spectral regions — one on each side of the peak
of interest. AutoGammaFit computes the mean count rate in each region and draws
a straight line through the two region centres:

```
B₁ = (ȳ₂ − ȳ₁) / (x̄₂ − x̄₁)
B₀ = ȳ₁ − B₁·x̄₁
```

This is more robust than a single-region fit when the continuum has a strong
slope, because it samples the background on *both* sides of the peak.

**When to use:** any time the background level changes noticeably across the fit
window (e.g. near a Compton edge, in a dense spectrum, or when the peak sits on
the shoulder of a much larger transition).

**Usage:**
1. Enter `Left Lo/Hi` — a clean off-peak region to the left of the peak
2. Enter `Right Lo/Hi` — a clean off-peak region to the right
3. Click **Apply Anchors → seed bg0, bg1**

*Reference: [Debertin & Helmer (1988)](#ref-5), Sec. 4.3 — background interpolation methods*

#### Flat background

Check **Flat background** to fix `B₁ = 0`. Use this after background subtraction
(View mode 2 or 3) where the continuum has already been removed and a residual
near-constant offset is all that remains.

#### Quadratic background

Check **Quadratic background** to add a curvature term:
```
BG(x) = B₀ + B₁·x + B₂·x²
```

**When to use:** when the residuals of a linear-background fit show a systematic
bowl or arch shape across the fit window. This indicates the Compton continuum
has curvature that a straight line cannot capture.

Common situations:
- Fitting a peak on the rising edge of a Compton peak from a higher-energy line
- Very wide fit windows (>20 FWHM) where the background curves measurably
- High-energy peaks (>2 MeV) where the Compton continuum has strong curvature

The B₂ parameter is seeded at zero and its limits are set conservatively to
prevent the quadratic term from dominating and absorbing peak area.

*Reference: [Debertin & Helmer (1988)](#ref-5), Ch. 3; standard practice in HPGe analysis*

### 4.3 Compton Step

**Physical background:** when a gamma photon enters a germanium detector, it can
undergo Compton scattering before depositing its full energy. Each Compton scatter
transfers part of the photon's energy to a recoil electron, which is detected. The
surviving photon then either escapes (partial energy deposition) or undergoes a
subsequent interaction. Photons that scatter once and then escape produce counts
at energies between zero and the Compton edge. When many such events accumulate,
they form a "shelf" — a step-like rise in counts on the low-energy side of the
full-energy peak.

This shelf is almost always present to some degree but only becomes significant
for:
- **High-energy peaks** (≥500 keV) where the Compton cross-section is large
- **Thick or high-Z surrounding materials** that scatter photons back into the
  detector
- **Tight-geometry sources** or detectors with significant back-scatter

**Model:** a complementary error function (erfc) is added for each Gaussian peak:
```
f_step,i(x) = Sᵢ · erfc( (x − Eᵢ) / (σᵢ · √2) )
```

The erfc function transitions from ~2Sᵢ (far left of peak) to ~0 (far right),
creating the characteristic step shape. The step amplitude Sᵢ is a free parameter,
seeded at 5% of the peak amplitude and constrained to [0, Aᵢ].

**When to use:**
- Visually obvious shelf on the low-energy side of the peak
- Residuals after a pure Gaussian fit show a systematic positive excess to the
  left of the peak centroid and a systematic negative residual at the peak itself
- Fitting the 511 keV annihilation peak, which almost always has a pronounced shelf
- Peaks at E > 500 keV in germanium

**When NOT to use:**
- Low-energy peaks (< 200 keV) where the shelf is negligible and adding Sᵢ will
  be underconstrained
- Low-statistics spectra where the step amplitude cannot be reliably determined —
  MIGRAD will push Sᵢ to a bound
- Background-subtracted spectra where the TSpectrum removal has already absorbed
  most of the shelf

*Reference: [Helmer & McCullagh (1979)](#ref-3) —
the original derivation and justification of the erfc step model for Ge detectors.
Also: [Debertin & Helmer (1988)](#ref-5), Sec. 4.4*

### 4.4 Tied Widths

**Physical background:** for a given detector, the Gaussian sigma of a full-energy
peak is determined by the detector resolution at that energy, described by the FWHM
model (Section 5). Every peak from the same detector in the same spectrum *must*
have a sigma consistent with this model.

In practice, when fitting a doublet (two peaks separated by less than ~2 FWHM),
MIGRAD has insufficient information to independently determine both sigmas. The
optimizer will often converge to an unphysical solution where one peak is very
narrow and the other very broad — a configuration that can fit the data just as
well statistically but has no physical meaning.

**What tied widths does:** fixes each sigma to the model prediction:
```
σ_fixed,i = FWHM(Eᵢ) / (2√(2·ln2)) = FWHM(Eᵢ) / 2.3548
```

This reduces the number of free parameters by N (one per peak), making the fit
more constrained and physically meaningful. The energies and amplitudes remain free.

**When to use:**
- **Doublets or triplets** where the peaks are within ~2 FWHM of each other
- When one peak is much weaker than the other (amplitude ratio > ~5:1) — the
  weak peak's sigma has too few constraining counts to converge independently
- When the free fit returns a sigma far outside the model prediction (e.g. 2× or
  0.5× σ_model) — indicates the fitter is compensating for something else
- Any time you have a well-calibrated FWHM model and want to enforce physical
  consistency

**When NOT to use:**
- Well-separated, well-resolved single peaks with good statistics — in this case
  the free sigma is a valuable cross-check on the FWHM model
- When your FWHM model is poorly calibrated — tying widths to a bad model will
  produce systematically wrong areas
- When an anomalously broad peak is physically expected (e.g. Doppler-broadened
  511 keV in flight, recoil-broadened lines from short-lived states)

**Requires:** a valid FWHM model loaded in the FWHM tab. If the model parameters
are all zero, tied widths will produce zero-sigma peaks.

*Reference: gnuScope tied-width fitting mode; [Debertin & Helmer (1988)](#ref-5), Sec. 4.2;
[Radford (1995)](#ref-4) — tied-width multiplet fitting*

### 4.5 Tied Widths — Same Sigma

An additional option under the width-tying controls forces **all Gaussians in the
current fit to share a single free sigma parameter**:

```
f(x) = Σᵢ Aᵢ · exp(−½·((x − Eᵢ)/σ)²)  +  BG(x)
```

All peak amplitudes and centroids remain free; only one sigma is fitted.

**When to use:**
- Doublets or multiplets where the peaks must have physically identical widths
  (e.g., an unresolved doublet from two transitions of the same nucleus at the
  same energy resolution)
- When the individual per-peak sigmas are not independently constrained (very
  close peaks with low statistics)

**When NOT to use:**
- Peaks separated by more than ~3 FWHM — real detectors show measurable width
  variation across energy, and a shared sigma will be wrong for at least one peak
- When the peaks are from transitions with genuinely different widths (e.g.,
  Doppler-broadened vs. narrow lines in the same window)

**Implementation note:** after fitting, the shared sigma value is copied back to
all peak entries in the cache so that the standard per-peak display and FWHM model
comparison work correctly.

### 4.6 Fit Options

| Option | Description |
|--------|-------------|
| Log-likelihood (L) | Minimise the Poisson log-likelihood. Preferred when any bin has fewer than ~10 counts. chi² is biased for low counts because the assumption of Gaussian bin errors breaks down. |
| IMPROVE (M) | Run IMPROVE after MIGRAD to search for a better global minimum. Adds ~2× fit time. Use when chi²/ndf is suspiciously high or parameters have hit bounds. |
| MINOS errors (E) | Compute asymmetric (profile-likelihood) parameter uncertainties. More accurate than symmetric MIGRAD errors near parameter boundaries. Significantly slower. |
| Show fit components | Overlay the background (green dashed) and individual Gaussian components (blue dashed) on the spectrum, in addition to the total fit (red). |
| Show isotope matches | Display the matched isotope name from the database above each peak energy label. Uncheck to show energies only. |

### 4.7 Fit Parameter Bounds

Click **Fit Parameters…** to edit the search bounds for each parameter class:

| Parameter | Default | Effect |
|-----------|---------|--------|
| Energy window ±keV | ±8 keV | How far the centroid can move from the clicked position |
| Amplitude lower fraction | 0.01× | Prevents amplitude collapsing to near-zero |
| Amplitude upper fraction | 20× | Prevents runaway amplitude |
| Sigma lower fraction | 0.2× model | Prevents unphysically narrow peaks |
| Sigma upper fraction | 4.0× model | Prevents runaway broadening |
| Range (σ) | window half-width | Auto fit range = peak ± Range × σ(E) |

Bounds are applied via MIGRAD's parameter limits. Note that parameters at their
bounds will have zero or underestimated uncertainties — the **WARNING: hit bound**
message in the log indicates this.

*Reference: [James & Roos (1975)](#ref-1), MINUIT documentation — behaviour of parameters
at limits*

### 4.8 Peak Statistics Panel

After each fit, the upper-right panel shows a block of statistics for every
Gaussian component. If a cached fit already exists for the same peak set, the
manual fit results appear first under `=== Manual Fit ===` and the cached fit
appears below under `=== Cached Fit ===` for direct comparison.

A complete reference with derivations and physical interpretation of every term
is given in [Section 12.14](#1214-complete-fit-statistics-reference). The quick
summary table is below; read §12.14 before making claims about weak peaks.

| Quantity | Formula | Quick interpretation |
|----------|---------|---------------------|
| Area ratio (fit/obs) | N_fit / N_obs\_above\_BG | 0.90–1.10 = GOOD; outside = model problem |
| Amplitude | A ± σ_A | Gaussian peak height in counts/bin |
| Centroid | E ± σ_E (NNDC) | Full-energy peak position in keV |
| Sigma | σ ± σ_σ | Gaussian width (not FWHM) |
| FWHM | 2.3548·σ ± 2.3548·σ_σ | Standard detector resolution measure |
| Peak area | A·σ·√(2π) / Δx | Total counts in peak above background |
| Area uncertainty | area·√((σ_A/A)² + (σ_σ/σ)²) | Propagated from fit covariance |
| SNR | N_peak / √N_BG | > 5 = clearly visible; < 3 = questionable |
| Significance (fit) | N_peak / σ(N_peak) | Fit-based significance; requires valid covariance |
| Significance (stat) | N_peak / √(N_peak + 2·N_BG) | Counting-based significance (p-value shown) |
| p-value (stat) | ½·erfc(S/√2) | Probability peak is background fluctuation |
| P/T (window) | N_peak / N_window | Fraction of window counts in peak |
| P/T (all) | N_peak / N_spectrum | Fraction of all spectrum counts in peak |
| Eff-corr area | N_peak / ε(E) | Emitted gammas; needs efficiency model loaded |
| bg0 | B₀ ± σ_B₀ | Background intercept (counts/bin at x = 0) |
| bg1 | B₁ ± σ_B₁ | Background slope (counts/bin per keV) |
| bg2 | B₂ ± σ_B₂ | Background curvature (quadratic term) |
| chi²/ndf | Σ(pull²) / (bins − params) | ≈ 1 = good fit; >> 1 = bad model |
| p-value (chi²) | TMath::Prob(χ², ndf) | > 0.05 = fit consistent with data |

### 4.9 Residuals

Check **Show residuals** to split the canvas vertically: spectrum + fit on top,
normalised pull residuals (data − fit) / σ_data on the bottom.

**Interpreting residuals:**
- Residuals uniformly scattered in [−2, +2] → good fit
- Systematic arch or bowl shape → wrong background model (try quadratic BG)
- Positive excess to the left, negative at peak centre → Compton step present
- Narrow positive spike between two peaks → missed peak (add another Gaussian)
- All residuals shifted by a constant → bg0 is wrong

### 4.10 Peak Navigation

Use **Prev / Next** or the **Select fit** dropdown to cycle through all cached fits
for the loaded histogram. The canvas zooms to each fit region automatically.
Use **+ / −** or the X range boxes to adjust the zoom level.

---

## 6. FWHM Tab

Measures and models the detector energy resolution FWHM(E).

**Why the FWHM model matters:** the resolution model is used throughout AutoGammaFit
to set sigma bounds in AutoFit, to seed sigma values in Manual Fit, and to enforce
physically consistent widths when Tie Widths is enabled. A well-calibrated model
directly improves fit quality and convergence.

**Physical origin of the resolution model:** the energy resolution of a semiconductor
detector has three contributions [[Knoll (2000)](#ref-6), Ch. 11]:

```
FWHM(E)² = a + b·E + c·E²
```

- **a** (constant term): electronic noise floor — readout amplifier and digitiser noise
- **b·E** (statistical term): Fano-factor-limited charge-carrier statistics.
  For Ge: b ≈ 8F·ε where F ≈ 0.1 is the Fano factor and ε ≈ 2.96 eV/pair
- **c·E²** (quadratic term): incomplete charge collection and ballistic deficit.
  Usually small for well-maintained detectors; non-zero for degraded detectors.

**Workflow:**
1. Run AutoFit on several calibrated spectra with well-known gamma lines
2. In the FWHM tab, select a histogram and click **Load** — extracts FWHM values
   from the cached fits
3. Click **Fit FWHM** — MIGRAD fits a, b, c to the FWHM vs E data
4. Inspect residuals; exclude outliers with Remove Mode
5. Click **Accept FWHM** — stores the model; it is now used by all fitting operations

*Reference: [Knoll (2000)](#ref-6), Ch. 11; [Debertin & Helmer (1988)](#ref-5), Sec. 2.4*

---

## 7. Source Tab

Calibration and efficiency analysis using a known radioactive source.

### 6.1 Load Source Description

A plain-text file describing the source:
```
isotope   Co-60
activity  37000        # Bq at caldate
halflife  1925.5       # days
caldate   2020-01-15
measdate  2025-03-10
livetime  3600.0       # seconds
# energy_keV  branching_ratio
1173.228  0.9985
1332.492  0.9998
```
Branching ratios (intensities) are per-decay values in the range (0, 1].
Reference values can be obtained from [NuDat 3.0](https://www.nndc.bnl.gov/nudat3/)
or the IAEA Nuclear Data Services.

### 6.2 Peak Assignment

**Auto Identify** matches each source line to the nearest cached fitted peak
within 3 FWHM of the reference energy.

**Manual Assign** lets you type an energy and assign it to the selected source
line — useful when auto-identification fails for a weak line or when the fitted
energy is slightly outside the tolerance.

### 6.3 Energy Calibration Plot

Plots (E_fit − E_ref) vs E_ref for all assigned lines. A good calibration gives
residuals scattered around zero with no systematic trend. A slope indicates a
gain error; a quadratic trend indicates nonlinearity in the ADC or shaping
amplifier.

Read off the polynomial coefficients from this plot to use as input for the
Channel→keV calibration (Section 6.6).

### 6.4 Efficiency vs Energy

Computes the absolute full-energy peak efficiency for each assigned line:
```
ε(E) = N_peak / (A_measured × BR × T_live)
```
where the source activity at measurement time is:
```
A_measured = A_cal × exp(−ln2 × Δt / T½)
```
and Δt is the time from `caldate` to `measdate`.

Plot the result on a log-log scale — the efficiency curve should be smooth and
monotonically decreasing (for Ge, above ~100 keV). Fit a log-polynomial to the
data (see Section 12.8) and enter the coefficients in the Efficiency Correction
group of the AutoFit tab.

*Reference: [Knoll (2000)](#ref-6), Sec. 12-VI; [Debertin & Helmer (1988)](#ref-5), Ch. 6*

### 6.5 FWHM vs Energy

Loads FWHM values from the source spectrum cache and plots FWHM vs Energy for
input to the FWHM tab model fitting.

### 6.6 Channel → keV Calibration

Applies a polynomial energy calibration to a spectrum stored in raw channel
numbers. This transforms the x-axis from ADC channels to keV.

**Model:**
```
E(ch) = a + b·ch + c·ch²
```

**Workflow:**
1. Fit known gamma lines in a channel-number spectrum using Manual Fit to get
   precise channel centroids
2. Use the Energy Calibration plot to fit the polynomial coefficients
3. Enter `a`, `b`, `c` in the calibration fields
4. Click **Apply to Selected Histogram**

A new histogram `<name>_cal` is created with a non-uniform energy axis derived
from the polynomial. Bin contents (counts) are unchanged; only the x-axis edges
are rescaled. The quadratic term `c` corrects for ADC nonlinearity and is often
small but non-negligible for 12-bit ADCs.

**Note:** for a purely linear calibration, set `c = 0`. The result is equivalent
to `TH1::GetXaxis()->SetLimits(E_min, E_max)`.

*Reference: [Knoll (2000)](#ref-6), Sec. 11-VII;
[Debertin & Helmer (1988)](#ref-5), Sec. 2.2*

---

## 8. Fit Results Tab

Lists all histograms with cached fits. Selecting an entry redraws the histogram
with all cached Gaussian overlays, energy labels, and isotope names.

### 7.1 Canvas Annotations

Add a custom title, axis labels, or free-text annotations (positioned in NDC
coordinates) to the displayed canvas before saving.

### 7.2 Export

- **Save Selected / Save All**: write ROOT TF1 objects and fit parameters to files
- **Export Cache CSV**: dump all fit entries for all histograms to a single CSV
  file suitable for import into spreadsheets or further analysis scripts
- **Export to gnuScope** / **Export All to gnuScope…**: write binary files
  readable by the gnuScope analysis program

### 7.3 gnuScope Binary Format

Files use the Fortran sequential unformatted record convention (4-byte record
length before and after each data block):

**1D (.spe):** single data record
```
[N×4]  [float₁ float₂ … floatN]  [N×4]
```
where N is the number of histogram bins and each float is the bin content.

**2D (.sqr):** header record followed by data record
```
Header: [12]  [maxXY]  [0]  [0]  [12]
Data:   [maxXY²×4]  [float matrix, row-major]  [maxXY²×4]
```
where maxXY = max(nX, nY). Non-square histograms are zero-padded to a
maxXY × maxXY matrix.

---

## 9. Isotopes Tab

Browse the gamma database and assign isotope labels to fitted peaks.

- **Auto Match All**: match each cached peak to the nearest database line within
  the configurable threshold (keV). Uses the loaded FWHM model as the matching
  window.
- **Apply Auto Matches**: commit all auto-matched labels to the cache files
- **Draw Schematic**: display a decay-chain diagram for the loaded parent nucleus
  using AME2020 and NuBase mass/half-life data (load via "Load AME Table" /
  "Load NuBase Table")

**Classification system:**
Parent, Daughter, Granddaughter, Beta-n Daughter, Beta-2n Daughter,
Beta-n Granddaughter, Beta-2n Granddaughter, Background, X-ray, Custom

Classifications are stored per peak label and used in the Decay tab and CSV exports.

*Reference: [NNDC NuDat database](#ref-10); [AME2020 — Wang et al. (2021)](#ref-11)*

---

## 10. Peak Table Tab

Aggregates peak data from multiple fit caches into a unified sortable/filterable
table. Use it to compare energies, areas, and intensities across all loaded
histograms.

### 10.1 Loading Caches

- **Scan All**: automatically discovers all `.dat` cache files in the current
  cache directory and adds them to the cache list.
- **Add Cache** / **Remove Cache**: manually select which cache files contribute
  to the table.
- **Rebuild Table**: recollects all entries from the selected caches and
  populates the table. Columns include: label, energy (NNDC notation),
  area ± error, chi²/ndf, sigma, classification, and optionally absolute
  intensity.

### 10.2 Intensity Calculation

Absolute emission probabilities (intensities) are computed as:
```
I = Area / (ε(E) × N)
```
where ε(E) is the detector efficiency at the peak energy (from the loaded
efficiency cache) and N is the **Total Population** — the total number of
decays of the parent nucleus during the measurement.

**Total Population field:** enter the number of decays directly, or click
**Populate from Decay Fit** to auto-fill from a decay fit result:
```
N = A₀ × T½ / ln(2)
```
where A₀ is the fitted initial activity (counts/ms at t=0) and T½ is the
parent half-life from the decay fit. This represents the total number of
parent decays from t=0 to t=∞.

### 10.3 Efficiency

Select a named efficiency cache from the dropdown (populated via **Scan Eff.
Caches**) or type the efficiency manually. The efficiency is evaluated at the
peak energy using the log-polynomial model (see §12.9).

### 10.4 Export

**Export CSV** writes all visible rows to a comma-separated file with all
columns, suitable for further analysis in ROOT macros or Python.

---

## 11. Decay Tab

Fits time-dependent gamma peak intensities extracted from a TH2 histogram
(gamma energy on one axis, time-since-beam-stop on the other). The tab is
divided into two sub-tabs: **Cuts** and **Fitter**.

### 11.1 Cuts Sub-tab

Controls which peak slice from the TH2 is used for the decay fit.

**2D Histogram group:**
- Select the TH2 from the dropdown and specify which axis carries gamma energy
  (X or Y). Changing this clears the peak list.

**Fitted Peaks group:**
- **Peaks cache**: select which gamma fit cache provides the peak list (or let
  it default to the projection of the TH2). Click **Scan** to populate the
  dropdown.
- **sigma window**: the ±N·σ integration window around the peak centroid used
  when projecting the TH2 onto the time axis. Determines which time-projection
  counts are attributed to this gamma line.
- **Refresh**: reloads the peak list from the selected cache, sorted by energy.
  Peaks are displayed as `label  E keV  (sig=σ)  [class]`.
- **Preview** (decay projection): projects the TH2 for the selected peak and
  sigma window and displays the raw time spectrum without fitting.
- **Rebin**: merge N adjacent time bins before fitting or preview. Click Apply or
  Reset.
- **Preview Gamma Peak**: shows the *gamma* spectrum zoomed to ±6σ around the
  selected peak, with red dashed lines marking the ±N·σ integration window.
  Use this to verify the cut before fitting.
- **Label / Apply Label**: set or update the isotope label stored in the gamma
  cache for the selected peak.
- **Load Cache / Save Cache**: load or save the decay fit results (half-lives,
  chi², parameters) to a per-TH2 `.dat` file in the cache directory.

**Peak Counts vs Time group** (see §11.6):
- Slice the TH2 into time bins and fit the gamma peak in each slice independently,
  producing a graph of peak counts vs time.

### 11.2 Fitter Sub-tab

Controls the decay model and runs the fit. After selecting a peak in the Cuts
sub-tab, switch to the Fitter sub-tab to configure and run the decay fit.

When a peak with a stored decay fit is selected, the Fitter sub-tab automatically
populates with the stored parameters and replays the fit on the canvas.

**Decay Model group:**
- **Model**: select the decay chain model (see §11.3).
- **BG type**: select the background functional form (see §11.4).
- **Equation**: shows the selected model formula for reference (auto-updates
  when the model or BG type changes).
- **T½ rows**: seed values and Fix checkboxes for the parent, daughter, and
  granddaughter half-lives. The row labels update automatically for β⁻n/β⁻2n
  models to show which chain level each T½ corresponds to.
- **T½ Exp BG**: seed for the exponential background T½ (only enabled when
  BG type is "Flat+Exp" or "Exp only").
- **Auto Model**: auto-selects the model ID from the selected peak's
  classification and calls Seed T½ (see §11.5).
- **Seed T½**: looks up the peak label in the nuclear database and fills T½
  from the chain (see §11.5).
- **Fit from / to / Full**: fit range on the time axis. Check Full to use the
  full histogram range.
- **Fit Decay**: runs the MIGRAD fit and displays results in the Fit Results
  panel below.

**Fit Results panel**: shows the fit equation, chi²/ndf, fitted parameters with
uncertainties, integrated peak counts in the fit window, extrapolated total peak
counts (0→∞), and an informational model classification.

### 11.3 Decay Models

| ID | Name | Signal formula | Notes |
|----|------|---------------|-------|
| 1 | Parent | A·exp(−ln2·t/T_P) | Single exponential decay |
| 2 | Daughter (β⁻) | Bateman 2-step | In-growth from parent β⁻ decay |
| 3 | Granddaughter (β⁻) | Bateman 3-step | A→B→C Bateman chain |
| 4 | Background only | (no signal) | Pure BG model |
| 5 | Daughter (β⁻n) | Bateman 2-step | Parent is a β⁻n product of higher nucleus |
| 6 | Daughter (β⁻2n) | Bateman 2-step | Parent is a β⁻2n product |
| 7 | Granddaughter (β⁻n) | Bateman 3-step | Second level is a β⁻n product |
| 8 | Granddaughter (β⁻2n) | Bateman 3-step | Second level is a β⁻2n product |

Models 5–8 use the same Bateman equations as 2–3 but with different T½ row labels
to clarify which decay branch each half-life belongs to. The mathematics and
parameter layout are identical.

For the Bateman equations and parameter layouts, see §12.13.

**Classification → model auto-mapping** (used by Auto Model):

| Classification | Model |
|---------------|-------|
| Parent | 1 |
| Daughter | 2 |
| Granddaughter | 3 |
| Background | 4 |
| bn Daughter | 5 |
| b2n Daughter | 6 |
| bn Granddaughter | 7 |
| b2n Granddaughter | 8 |

### 11.4 Background Models

The background model is selected independently of the signal model via the
**BG type** combo:

| BG type | Formula | When to use |
|---------|---------|-------------|
| Flat | BG | Constant beam-off background; pile-up pedestal |
| Flat + Exp | BG_flat + A_bg·exp(−ln2·t/T_bg) | Beam-related background that decays with its own T½ alongside a constant pedestal |
| Exp only | A_bg·exp(−ln2·t/T_bg) | Background with no constant component; all background follows an exponential decay |

The exponential background T½ (T_bg) is a free parameter seeded from the
**T½ Exp BG** entry. It can be fixed with the adjacent Fix checkbox.

**Note:** for the "Exp only" BG type applied to a "Background only" signal model
(model 4), the fit reduces to a single decaying exponential — useful for fitting a
beam-correlated background component without any in-growth signal.

### 11.5 Auto Model and T½ Seeding

**Auto Model button:**
1. Reads the classification of the selected peak from its gamma cache entry (or
   from the label→class map set in the Isotopes tab)
2. Maps the classification string to a model ID (table in §11.3)
3. Selects that model in the combo and updates the equation label
4. Calls Seed T½ automatically

**Seed T½ button:**
1. Reads the isotope label from the Label field
2. Looks up `nuclearDB_` (populated from NNDC via the Nuclear tab) using the key
   format `ASymbol` (e.g. `94Rb`)
3. Converts the half-life from seconds to milliseconds
4. Assigns T½ to the appropriate entry:
   - *Parent model*: fills T½ Parent
   - *Daughter models*: fills T½ Daughter; traverses `nucChainIsotopes_` to find
     the parent one level up and fills T½ Parent
   - *Granddaughter models*: fills T½ GDaughter; traverses two levels up for T½
     Daughter and T½ Parent

**Requires:** the nuclear chain must have been loaded and traced in the Nuclear tab.
Without chain data, only the direct T½ lookup for the labeled isotope is performed.

### 11.6 Peak Counts vs Time

An alternative approach that avoids the global Bateman fit. The TH2 is sliced into
time bins of width **Slice width (ms)**. In each slice, the gamma axis is projected
and the selected peak is fitted with a Gaussian + linear BG. The fitted peak counts
are plotted as a graph of counts vs time bin centre.

This is useful for:
- Visual inspection of the decay profile before fitting
- Cases where the background is not well-described by the Bateman assumption
- Identifying outlier time bins due to beam artifacts or pile-up

The fitted Gaussian amplitudes, sigmas, and centroids from each slice are stored
in a ROOT file in the cache directory:
```
fit_caches/<rootfile>/decay_slices_<TH2name>_<E>keV.root
```

### 11.7 Fit Results and Peak Population

After a successful decay fit, the **Fit Results** panel shows:

- Peak energy, sigma window, and integration cut [eMin, eMax]
- Chi²/NDF and fit status
- All fitted parameters with uncertainties (par names reflect the model)
- **Counts in fit window**: `(∫f(t)dt − BG·(thi−tlo)) / binWidth` — the
  signal integral over the fitted time range, corrected for background
- **Total peak counts (0→∞)**: `A × τ_P / binWidth` — extrapolated total
  signal counts assuming the decay started at t=0 and ran to infinity.
  For daughter and granddaughter models, the integral from 0→∞ still gives
  `A × τ_P` regardless of daughter half-lives (the Bateman equations conserve
  total production)
- Informational model classification (does not write to gamma cache)

The **Total Population** field in the Peak Table tab can be populated from
a decay fit by clicking **Populate from Decay Fit** in the Isotopes tab or
Peak Table tab. This uses `N = A₀ × T½_P / ln(2)` to give the total number of
parent decays.

*Reference: Bateman equations — [Bateman (1910)](#ref-12); ingrowth formalism —
[Debertin & Helmer (1988)](#ref-5), Appendix 1*

---

## 12. Fitting Models and Mathematics

### 12.1 Standard N-Gaussian + Linear Background

The core fit model for N peaks in a window [x_lo, x_hi]:
```
f(x) = Σᵢ Aᵢ · exp(−½·((x − Eᵢ)/σᵢ)²)  +  B₀  +  B₁·x
```
Parameters: Aᵢ (amplitude), Eᵢ (centroid), σᵢ (Gaussian sigma), B₀, B₁.
Total: **3N + 2** free parameters.

The Gaussian models the detector response to a mono-energetic gamma line, which
is broadened by charge-carrier statistics (Fano factor) and electronic noise.

### 12.2 Quadratic Background

Adds a curvature term to the background:
```
BG(x) = B₀ + B₁·x + B₂·x²
```
Total: **3N + 3** parameters. B₂ is seeded at zero.

Use when the Compton continuum has visible curvature across the fit window.

### 12.3 Compton Step Model

The Compton step arises because some photons deposit only a fraction of their
energy via single or multiple Compton scatters, creating a shelf of counts at
energies below the full-energy peak. The shelf has the shape of an error function:

```
f(x) = Σᵢ [Aᵢ·exp(−½·((x−Eᵢ)/σᵢ)²)  +  Sᵢ·erfc((x−Eᵢ)/(σᵢ·√2))]  +  BG(x)
```

The erfc function (complementary error function) transitions from ~2Sᵢ far to
the left of Eᵢ to ~0 far to the right. The step amplitude Sᵢ is constrained
to [0, Aᵢ].

Total parameters: **4N + NBG** where NBG = 2 (linear) or 3 (quadratic).

*Reference: [Helmer & McCullagh (1979)](#ref-3); also [Debertin & Helmer (1988)](#ref-5), Sec. 4.4*

### 12.4 Tied-Width Constraint

When enabled, each sigma is fixed to the resolution model prediction before fitting:
```
σ_fixed,i = FWHM(Eᵢ) / 2.3548
```
This removes N free parameters from the fit (one per peak), improving convergence
for close doublets and ensuring physical consistency.

*Reference: [Radford (1995)](#ref-4); gnuScope implementation*

### 12.5 Same-Sigma Constraint

When the same-sigma option is enabled, all N Gaussians in the fit share a single
free sigma parameter:
```
f(x) = Σᵢ Aᵢ · exp(−½·((x − Eᵢ)/σ)²)  +  BG(x)
```

The shared σ is estimated simultaneously with all amplitudes and centroids. This
reduces the total free-parameter count from 3N + N_BG to 2N + 1 + N_BG, saving
N − 1 degrees of freedom relative to an independent-sigma fit.

σ is constrained to [0.2 × σ_model, 4 × σ_model], where
`σ_model = FWHM(Ē) / 2.3548` is evaluated at the mean energy of the peak group.
After fitting, the shared sigma value is copied to all N peak entries in the cache
so that the resolution-model display and FWHM comparisons remain consistent.

**Use case:** doublets or multiplets where the peaks must have physically identical
widths (e.g., two transitions of the same nucleus falling at nearly the same energy,
or a case where the per-peak sigmas are individually underconstrained).

*Reference: shared-width constraint follows the same logical basis as the per-peak
tied-width mode — see [Radford (1995)](#ref-4) and [Debertin & Helmer (1988)](#ref-5), Sec. 4.2.*

### 12.6 Peak Area Calculation

ROOT stores histograms as counts-per-bin. The fitted amplitude A has units of
counts/bin, so the true peak area (total counts) is:

```
N_peak = Aᵢ · σᵢ · √(2π) / Δx
```

where Δx is the bin width in keV at the peak centroid. The factor √(2π) is the
integral of the unit Gaussian; dividing by Δx converts from the per-bin convention.

The area uncertainty is propagated in quadrature:
```
σ(N_peak) = N_peak · √((σ_A/A)² + (σ_σ/σ)²)
```

### 12.7 Signal-to-Noise Ratio (post-fit)

The post-fit SNR reported in peak statistics:
```
SNR = N_peak / √(N_BG)
```
where N_BG = |BG(Eᵢ)| · 5σᵢ / Δx is the estimated background count under
±2.5σ of the peak. This is the standard figure of merit for peak detectability.

*Reference: [Bevington & Robinson (2003)](#ref-7), Ch. 4; [Debertin & Helmer (1988)](#ref-5), Sec. 4.5*

### 12.8 Detector Resolution Model

```
FWHM(E) = √(a + b·E + c·E²)   [keV]
σ(E) = FWHM(E) / 2.3548
```

Physical interpretation of each term:
- **a** (noise): σ_noise² ∝ electronic noise variance
- **b·E** (statistical): FWHM_stat = 2.3548·√(F·ε·E), where F is the Fano factor
  (~0.1 for Ge) and ε is the energy per electron-hole pair (~2.96 eV for Ge)
- **c·E²** (incomplete collection): σ ∝ E, i.e. constant relative resolution;
  this term grows with detector degradation

For a pristine HPGe detector: a ≈ 0.5, b ≈ 1×10⁻³, c ≈ 0.

*Reference: [Knoll (2000)](#ref-6), Sec. 11-III; [Debertin & Helmer (1988)](#ref-5), Sec. 2.4*

### 12.9 Efficiency Model

The absolute full-energy peak efficiency ε(E) is the probability that a gamma
photon emitted by the source produces a count in the full-energy peak. It depends
on source-detector geometry, detector size, and self-absorption.

The log-polynomial model accurately describes ε(E) over a broad energy range
(typically 80–3000 keV):
```
ln(ε(E)) = a − b·ln(E) + c·(ln(E))² − d/E²
```

The 1/E² term accounts for photoelectric self-absorption at low energies.
The ln(E) and (ln(E))² terms describe the smooth decrease and curvature of
the efficiency at intermediate and high energies.

Efficiency-corrected peak area:
```
N_corrected = N_peak / ε(E)
```

This gives the number of gamma quanta emitted into the solid angle of the
detector, after correcting for detection probability.

*Reference: [Knoll (2000)](#ref-6), Eq. 12.8; [Debertin & Helmer (1988)](#ref-5), Sec. 6.2;
log-polynomial form also used in gnuScope [[Pavan & Tabor, FSU](#ref-14)]*

### 12.10 NNDC Uncertainty Notation

Fitted centroids are displayed in the standard nuclear data notation:
```
E(δ)
```
where δ is the uncertainty in units of the last displayed decimal place.

Examples:
- `1332.492(4)` → 1332.492 ± 0.004 keV
- `661.6(5)`    → 661.6 ± 0.5 keV
- `1460(2)`     → 1460 ± 2 keV

*Reference: [Browne & Tuli — NNDC notation convention](#ref-10);
[NuDat 3.0](https://www.nndc.bnl.gov/nudat3/)*

### 12.11 Chi-squared and Goodness of Fit

AutoGammaFit computes chi²/ndf directly from the pull residuals rather than
from ROOT's `r->Chi2()` (which returns −2 ln L for log-likelihood fits):
```
χ² = Σ_b [(y_b − f(x_b)) / σ_b]²
ndf = N_bins_in_window − N_free_parameters
chi²/ndf = χ² / ndf
```

**Interpretation:**
- chi²/ndf ≈ 1 → good fit; residuals are consistent with statistical noise
- chi²/ndf >> 1 → poor fit; systematic residuals or wrong model
- chi²/ndf << 1 → overestimated uncertainties or over-parameterised model

The p-value (probability that chi² ≥ observed, given correct model) is also
reported; p < 0.01 typically indicates a poor fit.

*Reference: [Bevington & Robinson (2003)](#ref-7), Ch. 6; [James & Roos (1975)](#ref-1)*

### 12.12 Peak Significance

Two independent significance metrics are reported for each fitted peak.

#### Significance (fit)

```
S_fit = N_peak / σ(N_peak)
```

This is the number of standard deviations the fitted peak area is above zero,
using the uncertainty propagated from MINUIT's covariance matrix:

```
σ(N_peak) = N_peak · √((σ_A/A)² + (σ_σ/σ)²)
```

**Interpretation:** S_fit answers the question *"how well-constrained is this
Gaussian by the fit?"* A high S_fit means the amplitude and width are well
determined; a low S_fit means one or both parameters are poorly constrained.

**Caveat:** if MINUIT converged on a boundary or did not produce a valid
covariance matrix, σ(N_peak) = 0 and S_fit is not reported. Always check
the **WARNING: hit bound** messages in the log.

#### Significance (stat)

```
S_stat = N_net / √(N_net + 2·N_BG)
```

where N_net is the fitted peak area and N_BG is the estimated background count
under ±2.5σ of the peak.

**Derivation:** the gross counts in the peak region are N_gross = N_net + N_BG.
Both terms follow Poisson statistics, so the variance of N_net is:

```
Var(N_net) = Var(N_gross) + Var(N_BG) = N_gross + N_BG = N_net + 2·N_BG
```

The factor of 2 on N_BG arises because the background must be *estimated* from
the fit rather than measured independently; estimating it carries its own
Poisson uncertainty equal to N_BG.

**Interpretation:** S_stat answers the question *"given the raw counts, how
likely is this peak to be a background fluctuation?"* It does not depend on
MINUIT convergence and remains valid even when fit errors are unreliable.

The reported one-sided p-value is:

```
p = ½ · erfc(S_stat / √2)
```

This is the probability that a pure-background measurement would produce at
least this many net counts by chance.

#### Threshold conventions

| S (σ) | p-value | Interpretation in nuclear counting |
|------:|--------:|-------------------------------------|
| < 2   | > 2.3 % | Not significant — likely background fluctuation |
| 2–3   | 2.3–0.13 % | Tentative — worth investigating |
| ≥ 3   | < 0.13 % | Evidence of a peak — above the Currie detection limit |
| ≥ 5   | < 3×10⁻⁷ | Significant — unambiguous peak |

In routine HPGe spectroscopy, 3σ is the standard minimum detection threshold
(the Currie critical level L_c). The 5σ convention is borrowed from particle
physics and is appropriate when claim of a new line must be very robust.

*Reference: [Currie (1968)](#ref-16) defines L_c (critical level) and L_D (detection
limit) from Poisson statistics. The N/√(N+2B) formula follows from standard
error propagation — [Bevington & Robinson (2003)](#ref-7), Ch. 3–4;
[Knoll (2000)](#ref-6), Sec. 3-V. One-sided p-value via the complementary
error function: [Bevington & Robinson (2003)](#ref-7), App. A.*

### 12.13 Decay Chain Mathematics

The decay fitter implements analytical solutions to the Bateman equations
[[Bateman (1910)](#ref-12)] for linear radioactive decay chains. All eight models
are built from three signal types (Parent, Daughter, Granddaughter) combined with
three background types (Flat, Flat+Exp, Exp-only).

#### Parent decay (Model 1)

```
S_P(t) = A · exp(−λ_P · t)
```

λ_P = ln2 / T½_P. Free parameters: A (initial count rate in counts/ms), T½_P (ms).

#### Daughter in-growth (Models 2, 5, 6)

The Bateman solution for a two-step linear chain A → B → [stable]:
```
         λ_D
S_D(t) = A · ————————————— · [exp(−λ_P·t) − exp(−λ_D·t)]
              λ_D − λ_P
```

where λ_P = ln2 / T½_P and λ_D = ln2 / T½_D. Free parameters: A, T½_P, T½_D.

At t = 0 the daughter count rate is zero; it grows as the parent populates it,
reaches a maximum at `t_max = ln(λ_D / λ_P) / (λ_D − λ_P)`, then decays as
the parent is exhausted.

**Degenerate case** (λ_P = λ_D): the formula reduces to:
```
S_D(t) = A · λ_P · t · exp(−λ_P · t)
```

Models 5 (β⁻n Daughter) and 6 (β⁻2n Daughter) use the same Bateman formula;
the half-life row labels are updated to clarify which chain level each T½
belongs to in the presence of delayed-neutron branching
[[Pfeiffer et al. (2002)](#ref-18)].

#### Granddaughter in-growth (Models 3, 7, 8)

The Bateman solution for a three-step linear chain A → B → C → [stable]:
```
S_G(t) = A · λ_D · λ_G · [           exp(−λ_P·t)
                             ——————————————————————————————
                             (λ_D − λ_P)(λ_G − λ_P)

                                         exp(−λ_D·t)
                           +  ——————————————————————————————
                              (λ_P − λ_D)(λ_G − λ_D)

                                         exp(−λ_G·t)
                           +  ——————————————————————————————  ]
                              (λ_P − λ_G)(λ_D − λ_G)
```

Free parameters: A, T½_P, T½_D, T½_G.

This expression assumes all three decay constants are distinct. When two decay
constants are equal (degenerate case), the limit must be taken by L'Hôpital's
rule; the implementation uses a small numerical separation ε = 10⁻⁶ ms⁻¹ to
avoid division-by-zero without specialising the code path.

Models 7 (β⁻n Granddaughter) and 8 (β⁻2n Granddaughter) use the same formula,
with T½ row labels updated for the β-delayed neutron chain topology
[[Pfeiffer et al. (2002)](#ref-18)].

#### Background models

The background accounts for beam-correlated backgrounds, scattered neutron
activation, or any time-varying background component that does not follow the
signal chain:

| BG type | Formula | Free parameters |
|---------|---------|----------------|
| Flat | BG | BG |
| Flat+Exp | BG_flat + A_bg · exp(−ln2 · t / T_bg) | BG_flat, A_bg, T_bg |
| Exp only | A_bg · exp(−ln2 · t / T_bg) | A_bg, T_bg |

The exponential component models a beam-correlated contamination that decays
with its own characteristic time T_bg independent of the signal chain.

#### Parameter layout

Parameters are ordered: [signal block] followed by [background block].

| Signal type | Signal parameters | N_sig |
|-------------|-----------------|-------|
| Parent | A, T½_P | 2 |
| Daughter | A, T½_P, T½_D | 3 |
| Granddaughter | A, T½_P, T½_D, T½_G | 4 |
| BG only | (none) | 0 |

| BG type | BG parameters | N_bg |
|---------|--------------|------|
| Flat | BG | 1 |
| Flat+Exp | BG_flat, A_bg, T_bg | 3 |
| Exp only | A_bg, T_bg | 2 |

Total parameter count = N_sig + N_bg. Signal parameters occupy indices [0, N_sig − 1];
background parameters occupy indices [N_sig, N_sig + N_bg − 1].

#### Conservation of total integral

For all signal types (Parent, Daughter, Granddaughter), the integral of the
signal function from 0 to ∞ satisfies:
```
∫₀^∞ S(t) dt = A / λ_P = A · T½_P / ln2
```

This holds because each step in the Bateman chain conserves total production:
every parent decay eventually produces exactly one daughter, which in turn
produces exactly one granddaughter. The total number of events in the peak
(0 → ∞) is therefore always set by the initial parent activity A and its mean
lifetime τ_P = 1/λ_P = T½_P / ln2, regardless of how many intermediate
steps exist. The reported **Total peak counts (0→∞)** uses this identity:
```
N_total = A · T½_P / (ln2 · Δt_bin)
```
where Δt_bin is the time bin width in ms.

*Reference: Bateman chain solutions — [Bateman (1910)](#ref-12).
Ingrowth formalism for detector efficiency measurements —
[Debertin & Helmer (1988)](#ref-5), Appendix 1.
β-delayed neutron emission topology —
[Pfeiffer et al. (2002)](#ref-18); half-life and decay-mode data —
[NUBASE2020 — Kondev et al. (2021)](#ref-19).*

### 12.14 Complete Fit Statistics Reference

This section defines and contextualises every quantity shown in the Peak
Statistics panel. The goal is that you can look at any line in that panel
and know exactly what it means, how it was computed, and whether the value
is reasonable for a real gamma-ray peak.

---

#### Area Ratio (fit/obs)

```
ratio = N_fit_total / N_obs_above_BG
```

Computed **before** the per-peak breakdown. `N_fit_total` is the sum of all
Gaussian integrals in the fit window. `N_obs_above_BG` is the sum of
`(bin_content − background_model)` for every bin in the window.

**What it means in a gamma spectrum:** if the fit model correctly describes
all peaks and the background model is accurate, these two quantities should
agree to within a few percent. The ratio is a global sanity check independent
of MINUIT convergence.

| Ratio | Grade | What to look for |
|-------|-------|-----------------|
| 0.90–1.10 | **GOOD** | Model accounts for essentially all the signal |
| 0.70–1.30 | **fair** | Mild mismatch — check background model |
| outside | **POOR** | Fit is missing a peak, background is wrong, or amplitude hit a bound |

A ratio > 1.1 usually means the background is underestimated (too low), so
the "observed above BG" is artificially small while the Gaussians picked up
extra area. A ratio < 0.9 usually means a peak was missed or the background
is overestimated.

*Reference: [Debertin & Helmer (1988)](#ref-5), Sec. 4.2*

---

#### Amplitude  (A ± σ_A)

The height of the fitted Gaussian in **counts per bin** at the centroid.
It is a raw MINUIT fit parameter, not directly interpretable as a physical
quantity — the same physical peak will have a different amplitude if you
change the bin width or rebin the histogram.

**In a gamma spectrum:** amplitude scales with count rate, bin width, and
detector efficiency. A large amplitude simply means many counts landed near
the peak centroid. To get a physically meaningful quantity, use the
**peak area** instead (which accounts for the bin width and peak width).

The quoted uncertainty σ_A comes directly from the diagonal element of
MINUIT's covariance matrix. It is reliable only when the fit has converged
and no parameter is at a bound. A zero or very small σ_A always means
a parameter hit a bound — check the log.

---

#### Centroid  (E ± σ_E keV,  NNDC notation)

The mean of the fitted Gaussian, in keV. This is the best estimate of the
**full-energy peak position** — the energy at which the detector most
frequently registers gamma photons of that transition.

**In a gamma spectrum:** the centroid should match the tabulated gamma-ray
energy to within your calibration uncertainty. If it doesn't, either the
calibration is wrong, the peak is a different transition, or a nearby line
is pulling the centroid (doublet).

The uncertainty σ_E reflects both the statistical precision of the
measurement (more counts → smaller σ_E) and the stiffness of the χ²
surface in the energy direction (correlated with σ_A). For a well-resolved
peak with ~1000 counts, σ_E ≈ FWHM/(2.35·√N) ≈ a fraction of a keV.

The energy is displayed in NNDC parenthesis notation: `1332.492(4)` means
`1332.492 ± 0.004 keV`. See [Section 12.10](#1210-nndc-uncertainty-notation).

*Reference: uncertainty on centroid — [Bevington & Robinson (2003)](#ref-7), Sec. 6.3;
NNDC notation — [Browne & Tuli](#ref-10)*

---

#### Sigma  (σ ± σ_σ keV)

The standard deviation of the fitted Gaussian, in keV. Along with the
centroid and amplitude, sigma is one of the three core fit parameters.

**In a gamma spectrum:** sigma is determined by the detector energy resolution.
For an HPGe detector, a typical value at 1.33 MeV is σ ≈ 0.5–0.6 keV
(FWHM ≈ 1.2–1.4 keV). Sigma increases with energy roughly as √E (statistical
term) plus a constant noise floor. A sigma much larger than the model
prediction indicates a doublet or a detector problem. A sigma much smaller
indicates the peak is undersampled (too few counts) or has hit a lower bound.

When **Tie Widths** is enabled, sigma is fixed to the FWHM model prediction
and σ_σ = 0.

---

#### FWHM  (keV ± keV)

```
FWHM = 2√(2·ln 2) · σ = 2.3548 · σ
σ(FWHM) = 2.3548 · σ_σ
```

The full width at half maximum of the fitted peak. FWHM is the conventional
measure of detector energy resolution and is directly comparable across
detectors and published specifications.

**In a gamma spectrum:** FWHM is quoted at specific calibration energies.
For HPGe detectors, a commonly cited benchmark is the FWHM at the 1332 keV
Co-60 line. Typical values:
- Excellent HPGe: ≤ 1.8 keV FWHM at 1332 keV
- Good HPGe: 1.8–2.5 keV
- NaI(Tl): ~60 keV (about 4.5% at 1332 keV)

FWHM in keV is equivalent to energy resolution R:
```
R (%) = 100 × FWHM(E) / E
```

For the 1332 keV Co-60 line at 1.8 keV FWHM: R = 0.135%.

The uncertainty on FWHM propagates directly from σ_σ. It is NOT the
uncertainty on the peak position — that is σ_E (the centroid uncertainty).

*Reference: resolution definition — [Knoll (2000)](#ref-6), Sec. 11-III;
[Debertin & Helmer (1988)](#ref-5), Sec. 2.4; IEEE Std 325 [[Leo (1994)](#ref-17)]*

---

#### Peak Area  (N_peak ± σ_N counts)

```
N_peak = A · σ · √(2π) / Δx
σ(N_peak) = N_peak · √( (σ_A/A)² + (σ_σ/σ)² )
```

The total number of counts in the Gaussian peak, above the background.
Δx is the bin width in keV at the peak centroid. This formula converts
from ROOT's per-bin convention (amplitude in counts/bin) to the true
count integral over all energies.

**In a gamma spectrum:** peak area is the physically fundamental quantity.
It is proportional to:
```
N_peak = A_source · BR · ε(E) · Ω · T_live
```
where A_source is the source activity (decays/s), BR is the branching ratio
(photons per decay), ε(E) is the full-energy peak efficiency, Ω is the
solid-angle factor, and T_live is the live measurement time.

All downstream quantities — efficiency, activity, isotope ratios — derive
from the peak area. The amplitude alone is **not** the right quantity to use.

**Warnings displayed in the panel:**
- `WARNING: area > spectrum total` — the Gaussian integral exceeds all counts
  in the entire spectrum. The fit is almost certainly wrong (amplitude at
  upper bound, very broad sigma, or wrong model).
- `WARNING: negative area` — A or σ is negative or zero. Check bounds.

The area uncertainty propagates A and σ in quadrature. It does **not** include
uncertainty from the background model; for that, use the BG anchor approach
(§4.2) to estimate BG systematic uncertainty separately.

*Reference: area formula — [Debertin & Helmer (1988)](#ref-5), Sec. 4.2;
error propagation — [Bevington & Robinson (2003)](#ref-7), Ch. 3;
[Knoll (2000)](#ref-6), Sec. 12-IV*

---

#### Signal-to-Noise Ratio  (SNR)

```
SNR = N_peak / √(N_BG)
```

where `N_BG = |BG(E)| · 5σ / Δx` is the estimated number of background
counts under ±2.5σ of the peak. BG(E) is the fitted background polynomial
evaluated at the peak centroid, and Δx is the bin width.

**What it means:** SNR compares the signal (peak counts) to the statistical
noise of the underlying background (Poisson fluctuations ≈ √N_BG). A peak
with SNR = 10 means the signal is 10 standard deviations above the expected
random fluctuation of the background in the peak region.

**In a gamma spectrum:**
- SNR >> 10: dominant, clearly visible peak
- SNR ≈ 5–10: well-defined peak
- SNR ≈ 3–5: marginal; visible by eye but requires careful background treatment
- SNR < 3: indistinguishable from background noise

Note: SNR uses the *fitted* background level. If the background model is wrong,
SNR will be wrong too. Compare SNR with `Significance (stat)`, which uses
the same N_BG estimate but from a different angle.

*Reference: [Bevington & Robinson (2003)](#ref-7), Ch. 4;
[Debertin & Helmer (1988)](#ref-5), Sec. 4.5*

---

#### Significance (fit)

```
S_fit = N_peak / σ(N_peak)
```

The number of standard deviations the fitted peak area is above zero,
using the uncertainty from MINUIT's covariance matrix.

**What it means:** S_fit answers *"how well-constrained is the Gaussian
amplitude by the fitter?"* It is a measure of fit confidence, not of
counting statistics. A high S_fit means the likelihood surface has a
clear minimum with respect to the peak amplitude — the fitter is certain
this Gaussian is needed. A low S_fit means the amplitude could shrink to
zero without significantly worsening chi².

**In a gamma spectrum:** S_fit can be large even for a weak physical peak
if the fit is well-behaved (good convergence, low background). Conversely,
it can be small for a real peak if the fitter is poorly constrained
(e.g., doublet with correlated parameters).

**Caveat:** S_fit is only trustworthy when:
1. MIGRAD converged (status = 0 or 1)
2. No parameters hit their bounds (`WARNING: hit bound` absent in log)
3. EDM is small (< 0.01)

If any of these conditions fail, S_fit may be artificially inflated or
suppressed. In that case, rely on `Significance (stat)` instead.

*Reference: parameter significance from covariance —
[James & Roos (1975)](#ref-1); [Bevington & Robinson (2003)](#ref-7), Sec. 6.2*

---

#### Significance (stat)  and  p-value

```
S_stat = N_peak / √(N_peak + 2·N_BG)

p = ½ · erfc( S_stat / √2 )
```

The Poisson counting significance, entirely independent of MINUIT convergence.
N_BG is the same estimated background count used in SNR (background under
±2.5σ of the peak).

**Derivation:** The net peak count is `N_net = N_gross − N_BG`. Both are
Poisson random variables, so:
```
Var(N_net) = Var(N_gross) + Var(N_BG) = N_gross + N_BG = N_net + 2·N_BG
```

The factor of 2 appears because the background must be *estimated* from the
fit rather than measured in an independent blank measurement — this estimation
carries its own Poisson uncertainty of N_BG. Under the null hypothesis
(background only), Var(N_net) ≈ 2·N_BG.

**What it means:** S_stat answers *"if only background were present, how
unlikely is it to observe this many net counts by chance?"* The p-value is
that probability.

**In a gamma spectrum:** S_stat is your primary detection threshold.
The standard critical level L_c in nuclear counting (after Currie, 1968)
corresponds to roughly S_stat ≈ 3σ (p ≈ 0.0013). Below this level, a peak
cannot be claimed with confidence.

| S_stat | p-value | Interpretation |
|-------:|--------:|----------------|
| < 2σ   | > 0.023 | Not significant — likely background fluctuation |
| 2–3σ   | 0.023–0.0013 | Tentative — worth further investigation |
| ≥ 3σ   | < 0.0013 | Evidence — above the Currie critical level |
| ≥ 5σ   | < 3×10⁻⁷ | Significant — unambiguous full-energy peak |

**What `p < 1e-10` means:** the p-value is smaller than 1 in 10 billion —
the peak is essentially certain to be real from a counting statistics
perspective. This typically occurs for any strong gamma line (>5σ).

**When S_stat and S_fit disagree:**
- S_stat >> S_fit: the peak is clearly present in the counts but the fitter
  is poorly constrained (e.g., doublet, parameter correlated with background).
  Trust S_stat; consider using tied widths or BG anchors.
- S_fit >> S_stat: the fitter found a clean minimum but the raw counts are
  marginal. This happens for narrow peaks in a high-background region.
  Both quantities are reporting real information; the peak is formally weak.

*Reference: [Currie (1968)](#ref-16) — critical level L_c and detection
limit L_D; [Knoll (2000)](#ref-6), Sec. 3-V;
[Bevington & Robinson (2003)](#ref-7), Ch. 3–4*

---

#### Peak-to-Total Ratio  P/T (window) and P/T (all)

```
P/T (window) = N_peak / N_window
P/T (all)    = N_peak / N_spectrum
```

`N_window` = total counts in the fit window (signal + background + all peaks).
`N_spectrum` = total counts in the entire loaded histogram.

**What P/T (window) means:** the fraction of counts in the local fit window
that belong to this peak. A ratio of 0.80 means 80% of the counts in the
window are under this Gaussian — the background and any other peaks account
for the remaining 20%. Useful for checking whether the window is dominated
by the peak or by background/neighbours.

**What P/T (all) means:** the peak area as a fraction of the entire spectrum.
This is the classical **peak-to-total ratio**, a standard figure of merit for
detector quality. A detector with high P/T has most of its counts in
full-energy peaks rather than Compton continuum.

**Typical values for HPGe at 1332 keV:**
- P/T (all) ≈ 0.50–0.65 for a 100% relative efficiency HPGe in close geometry
- P/T (all) ≈ 0.20–0.35 for a small-volume HPGe
- P/T (all) < 0.10 for NaI(Tl)

In practice, P/T (all) from a single peak in a mixed-isotope spectrum is not
a good detector figure of merit — it depends on the source composition. Use
it for single-line sources (e.g., Co-60 at 1332 keV in a clean spectrum).

*Reference: peak-to-total ratio — [Knoll (2000)](#ref-6), Sec. 12-II;
[Debertin & Helmer (1988)](#ref-5), Sec. 3.2*

---

#### Efficiency-Corrected Area  (Eff-corr area)

```
N_corrected = N_peak / ε(E)
```

Only displayed when the log-polynomial efficiency model has been loaded
(AutoFit tab → Efficiency Correction). ε(E) is evaluated from:
```
ln ε(E) = a − b·ln(E) + c·(ln E)² − d/E²
```

**What it means:** N_corrected is the number of gamma photons that were
emitted (into the detector solid angle) at energy E, after correcting for
the probability that any one photon was actually detected as a full-energy
event. It removes the detector's energy-dependent response.

**In a gamma spectrum:** combining N_corrected from multiple lines of the
same isotope should give activity × branching ratio × solid angle, which
is constant across lines if the efficiency model is correct. Inconsistency
across lines indicates either a wrong efficiency model or an interfering peak.
For absolute activity determination:
```
Activity (Bq) = N_corrected / (BR × T_live × Ω / 4π)
```

*Reference: [Knoll (2000)](#ref-6), Sec. 12-VI;
[Debertin & Helmer (1988)](#ref-5), Sec. 6.2*

---

#### Background Parameters  (bg0, bg1, bg2)

The fitted polynomial background under the peaks:
```
BG(x) = bg0 + bg1·x              (linear, default)
BG(x) = bg0 + bg1·x + bg2·x²     (quadratic, when enabled)
```

**What they mean in a gamma spectrum:** the background represents the Compton
continuum — scattered photons from higher-energy transitions that deposit
only part of their energy in the detector. bg1 (slope) reflects the gradient
of the Compton edge region at that energy; bg0 (intercept) reflects the
absolute level of the continuum.

**Interpreting bg0 and bg1:**
- bg0 units: counts/bin (this is the background at x = 0 keV, which is
  usually not physically meaningful; what matters is bg0 + bg1·E at the
  peak energy)
- bg1 units: (counts/bin) per keV — the slope of the continuum
- A negative bg1 means the continuum decreases with energy (normal for
  Compton continuum, which drops off above the Compton edge)
- A positive bg1 means the continuum rises with energy (unusual; check
  whether the fit window crosses a Compton edge from a higher-energy line)

**When bg2 is non-zero:** the quadratic term indicates curvature in the
continuum across the fit window. This is physical when the window sits on
the rising or falling edge of a broad Compton feature.

*Reference: background interpolation — [Debertin & Helmer (1988)](#ref-5),
Sec. 4.3; [Knoll (2000)](#ref-6), Sec. 12-IV*

---

#### chi²/ndf  (Goodness of Fit)

```
χ² = Σ_b [ (y_b − f(x_b)) / σ_b ]²
ndf = N_bins_in_window − N_free_parameters
chi²/ndf = χ² / ndf
```

Computed from the normalised pull residuals at every bin in the fit window.
Each term `(y_b − f(x_b)) / σ_b` is called a **pull** — it is the
discrepancy between data and model, measured in units of the bin uncertainty.

AutoGammaFit computes chi²/ndf directly from the pulls rather than using
ROOT's `r->Chi2()`, because the ROOT value is unreliable for log-likelihood
fits on background-subtracted histograms.

**What it means:**
- Each pull should be approximately a standard normal random variable if the
  model is correct.
- chi²/ndf = 1 means the average squared pull is 1 — the fit agrees with
  data at the level of statistical noise. This is ideal.
- chi²/ndf > 1 means the model is over-predicting scatter, or there are
  systematic discrepancies (wrong model).
- chi²/ndf < 1 usually means errors are overestimated, or the model has
  too many parameters for the available data.

**Typical diagnostic values:**

| chi²/ndf | Interpretation | Action |
|---------|---------------|--------|
| 0.8–1.2 | Good fit | Accept |
| 1.2–2.0 | Marginal | Check residual plot — may be acceptable |
| 2–5 | Poor | Try quadratic BG, Compton step, or add a peak |
| > 5 | Very poor | Wrong model — inspect residuals |
| < 0.5 | Over-fit / errors inflated | Check rebin factor; too many parameters |

**Note on degrees of freedom:** if the fit window has fewer bins than
parameters (e.g., very narrow window with many peaks), ndf ≤ 0 and chi²/ndf
falls back to chi²/N_bins. This is labelled in the display.

*Reference: chi² goodness of fit — [Bevington & Robinson (2003)](#ref-7),
Ch. 6; [James & Roos (1975)](#ref-1)*

---

#### p-value  (chi²)

```
p = P(χ² ≥ χ²_observed | ndf)  = TMath::Prob(chi2, ndf)
```

The probability, assuming the fit model is **correct**, of observing chi²
this large or larger by pure statistical chance.

**What it means:** p is a one-tailed tail probability of the chi² distribution
with ndf degrees of freedom. A small p means it is unlikely that the data
would look this bad if the model were right — i.e., the model is probably wrong.
A large p means the data are consistent with the model.

| p-value | Interpretation |
|---------|---------------|
| > 0.05 | Fit is consistent with data |
| 0.01–0.05 | Marginal — worth checking residuals |
| < 0.01 | Poor fit — model likely wrong |
| < 0.001 | Very poor fit |

**Common confusions:**
- A large p-value does **not** prove the model is correct — only that the
  data do not significantly contradict it.
- The p-value on chi² (fit quality) is completely separate from the p-value
  on `Significance (stat)` (peak detection). Both are shown; they answer
  different questions.

*Reference: [Bevington & Robinson (2003)](#ref-7), Ch. 6;
[Barlow & Beeston (1993)](#ref-9)*

---

#### Quick Diagnostic Checklist

When evaluating a fit in the statistics panel, check these in order:

1. **Area ratio** — should be 0.90–1.10. If not, the model is missing
   something before looking at individual peak quantities.
2. **chi²/ndf and p-value** — should be ≈ 1 and > 0.05. If chi²/ndf >> 1,
   the individual peak stats are unreliable.
3. **Significance (stat)** — must be ≥ 3σ to claim a peak. The p-value
   directly quantifies false-positive risk.
4. **Significance (fit)** — corroborates (stat) when fit errors are valid.
   Check the log for `WARNING: hit bound`.
5. **FWHM** — should match the resolution model (FWHM tab). A value 2× the
   model suggests a doublet; 0.5× suggests a parameter bound was hit.
6. **Peak area and uncertainty** — the ratio N_peak/σ(N_peak) should agree
   with `Significance (fit)` within rounding.

*Reference: full pipeline guidance — [Debertin & Helmer (1988)](#ref-5),
Sec. 4.5–4.7; [Leo (1994)](#ref-17), Ch. 7*

---

## 13. Algorithms

### 13.1 AutoFit Pipeline

1. **Background subtraction**: TSpectrum::Background uses the SNIP (Statistics-
   sensitive Non-linear Iterative Peak-clipping) algorithm. The background estimate
   is subtracted from a working copy of the histogram.
   *Reference: [Morhac et al. (1997)](#ref-8)*

2. **Peak search**: TSpectrum::Search scans the background-subtracted histogram
   for local maxima above the threshold.
   *Reference: [Morhac et al. (1997)](#ref-8)*

3. **S/N pre-screen** (if enabled): each group's signal-to-noise ratio is estimated
   from the ratio of peak-window signal to sideband noise on the BG-subtracted
   spectrum. Groups below the threshold are skipped.
   *Technique: gnuScope GetSumFF [[Pavan & Tabor, FSU](#ref-14)]; significance criterion: [Bevington & Robinson (2003)](#ref-7)*

4. **Peak grouping**: PeakGrouper merges peaks whose fit windows overlap (separation
   < 3σ of the lower-energy peak) into a single fit group.

5. **Adaptive fitting** (AdaptiveFitter::FitGroup): fits models with 1, 2, and up
   to maxPeaksPerGrp Gaussians + linear BG. Selects the model with the lowest
   AIC-penalised chi²/ndf score:
   ```
   score = chi²/ndf + kPenalty × (N_peaks − 1)
   ```
   Default kPenalty = 0.05; higher values favour simpler models.
   *Reference: [Akaike (1974)](#ref-2)*

6. **Double-Gaussian override**: the 511 keV annihilation peak is automatically
   modelled with a narrow + broad Gaussian (Doppler broadening). The model with
   the lower AIC score is kept.

7. **Residual-guided expansion**: if the best N-peak model leaves systematic
   residuals suggesting a missed peak, an additional Gaussian is seeded at the
   residual maximum and the group is re-fitted.

### 13.2 MIGRAD Minimisation

Both AutoFit and Manual Fit use MINUIT2 MIGRAD, a variable-metric gradient-descent
method that builds an approximation to the inverse Hessian of the log-likelihood or
chi² surface. MIGRAD converges when the estimated distance to the minimum (EDM) falls
below a tolerance.

The covariance matrix (and hence parameter uncertainties) is estimated from the
inverse Hessian at the minimum. Parameters at or near their bounds will have
underestimated uncertainties — check the **WARNING: hit bound** messages.

Optional: **MINOS** (option E) computes asymmetric uncertainties by scanning the
log-likelihood profile for each parameter independently. More accurate near bounds
but significantly slower.

*Reference: [James & Roos (1975)](#ref-1); MINUIT2 manual, CERN Program Library D506*

### 13.3 Chi² vs Log-likelihood

**Chi²:** minimises Σ (y_i − f(x_i))² / σᵢ² where σᵢ = √y_i (Poisson).
Valid when y_i >> 1 (typically y_i > 5–10 per bin). Biased for low counts because
the Poisson distribution is not symmetric.

**Poisson log-likelihood (option L):** minimises −2 ln L = 2 Σ [f(x_i) − y_i·ln f(x_i)].
Valid for any count level including zero. Preferred for weak peaks and low-statistics
spectra. The chi²/ndf interpretation of the fit quality is less direct.

*Reference: [Barlow & Beeston (1993)](#ref-9); [Bevington & Robinson (2003)](#ref-7), Ch. 8*

### 13.4 Background Anchor Algorithm

Given two anchor regions [lo₁, hi₁] and [lo₂, hi₂]:
1. Compute the mean bin content in each region: (x̄₁, ȳ₁) and (x̄₂, ȳ₂)
2. Fit a straight line through the two region centres:
   ```
   B₁ = (ȳ₂ − ȳ₁) / (x̄₂ − x̄₁)
   B₀ = ȳ₁ − B₁ · x̄₁
   ```
3. These values seed bg0 and bg1 as starting parameters for MIGRAD

The method samples the background on both sides of the peak, making it more robust
than a single off-peak region estimate when the continuum has a significant slope.

*Reference: [Debertin & Helmer (1988)](#ref-5), Sec. 4.3; [Knoll (2000)](#ref-6), Sec. 12-IV*

### 13.5 AIC Model Selection

When comparing models with different numbers of parameters, the Akaike Information
Criterion penalises added complexity:
```
AIC = 2k − 2 ln(L_max)
```
where k is the number of free parameters. For Gaussian noise this reduces to:
```
AIC ≈ χ² + 2k
```
AutoGammaFit uses a simplified per-peak penalty:
```
score = chi²/ndf + kPenalty × (N_peaks − 1)
```
with default kPenalty = 0.05. This means adding a peak must reduce chi²/ndf by
at least 0.05 to be favoured.

*Reference: [Akaike (1974)](#ref-2); [Burnham & Anderson (2002)](#ref-2), Ch. 2*

---

## 14. Cache System

Fit results are stored as plain-text cache files in `fit_caches/<rootfile>/`:
```
fit_caches/<rootfile>/fit_cache_<histogram_name>.dat
```

Each cache file contains one entry per fit region. Each entry records:
- Gaussian parameters (A, E, σ) with MIGRAD uncertainties
- Background parameters (B₀, B₁)
- Chi²/ndf, fit status, fit range [xlo, xhi]
- Peak label (isotope name) and classification
- Fit method (chi² or log-likelihood)

Special entries:
- `__RESOLUTION__`: three-parameter FWHM model (a, b, c)
- `__EXCLUDED_FWHM__`: list of energies excluded from the FWHM fit
- `rebin <name> <factor>`: stored rebin factors per histogram

Archived caches are stored in `fit_caches/<rootfile>/archive/` with a timestamp
suffix and can be restored via the **Restore Archived Cache** button.

### 14.1 Gamma Fit Cache Format

Each non-special line in a gamma cache file has the form:
```
<key>  <npar>  <chi2ndf>  <p0> <e0> ... <pN-1> <eN-1>  <residRMS>  <maxPull>
       <fitMethod>  <xlo>  <xhi>
```

- `key`: isotope label string (or `__RESOLUTION__`, etc.)
- `npar`: number of fit parameters
- `chi2ndf`: chi²/ndf
- `p0 e0 … pN-1 eN-1`: alternating parameter values and MIGRAD uncertainties
- `residRMS`, `maxPull`: residual diagnostics
- `fitMethod`: `0` = chi², `1` = log-likelihood
- `xlo xhi`: fit range in keV

META header lines (before any fit entries):
```
# META bg_subtracted 0|1
# META bg_iterations N
# META root_file <path>
```

### 14.2 Decay Fit Cache Format

Decay fit results are stored in:
```
fit_caches/<rootfile>/decay_cache_<TH2name>.dat
```

Each entry records one decay fit result for one peak:
```
<peakE>  <model>  <bgType>  <chi2ndf>  <status>  <eMin>  <eMax>  <Nsig>
<label>  <classification>
<npar>  <p0> <e0>  <p1> <e1>  …  <pN-1> <eN-1>
```

- `peakE`: peak centroid in keV
- `model`: decay model ID 1–8 (see §11.3)
- `bgType`: background type 1 = Flat, 2 = Flat+Exp, 3 = Exp only (see §11.4).
  Old cache files that predate the adaptable BG fitter omit this field and
  default to `bgType = 1` when loaded.
- `chi2ndf`, `status`: MIGRAD fit quality
- `eMin`, `eMax`, `Nsig`: TH2 cut range and sigma window used to generate the
  time projection
- `label`, `classification`: isotope label and classification string from the
  gamma cache
- `npar`, `p0 e0 …`: parameter values and uncertainties. Parameter ordering
  follows the layout table in §12.13 (signal block first, background block last).

Cache files are plain text — they can be inspected, edited, or merged with a text
editor if needed.

---

## 15. References

<a id="ref-1"></a>
1. **MIGRAD minimisation and MINOS errors:**
   F. James & M. Roos, "MINUIT — A System for Function Minimisation and Analysis of
   the Parameter Errors and Correlations," *Comput. Phys. Commun.* **10** (1975) 343–367.
   MINUIT2 manual: CERN Program Library D506.

<a id="ref-2"></a>
2. **AIC model selection:**
   H. Akaike, "A new look at the statistical model identification,"
   *IEEE Trans. Autom. Control* **19** (1974) 716–723.
   Practical guide: K.P. Burnham & D.R. Anderson, *Model Selection and Multimodel
   Inference*, 2nd ed., Springer, 2002.

<a id="ref-3"></a>
3. **Compton step (erfc shelf model):**
   R.G. Helmer & C.W. McCullagh, "Analytical functions for fitting peaks from Ge(Li)
   detectors," *Nucl. Instrum. Methods* **168** (1979) 593–599.

<a id="ref-4"></a>
4. **Tied-width multiplet fitting:**
   D.C. Radford, "ESCL8R and LEVIT8R: Software for interactive graphical analysis of
   HPGe coincidence data sets," *Nucl. Instrum. Methods* A **361** (1995) 297–305.

<a id="ref-5"></a>
5. **HPGe spectroscopy — BG methods, FWHM models, efficiency curves, peak areas:**
   K. Debertin & R.G. Helmer, *Gamma- and X-ray Spectrometry with Semiconductor
   Detectors*, North-Holland, Amsterdam, 1988.

<a id="ref-6"></a>
6. **Radiation detection — efficiency, resolution, statistics:**
   G.F. Knoll, *Radiation Detection and Measurement*, 4th ed., Wiley, 2010.

<a id="ref-7"></a>
7. **Statistical methods for experimental physics:**
   P.R. Bevington & D.K. Robinson, *Data Reduction and Error Analysis for the
   Physical Sciences*, 3rd ed., McGraw-Hill, 2003.

<a id="ref-8"></a>
8. **TSpectrum background estimation and peak search:**
   M. Morhac et al., "Background elimination methods for multidimensional coincidence
   γ-ray spectra," *Nucl. Instrum. Methods* A **401** (1997) 113–132.
   ROOT class documentation: `TSpectrum`.

<a id="ref-9"></a>
9. **Poisson log-likelihood for binned data:**
   R. Barlow & C. Beeston, "Fitting using finite Monte Carlo samples,"
   *Comput. Phys. Commun.* **77** (1993) 219–228.
   Also: W.T. Eadie et al., *Statistical Methods in Experimental Physics*,
   North-Holland, 1971.

<a id="ref-10"></a>
10. **NNDC notation and nuclear data:**
    E. Browne & J.K. Tuli, *Table of Radioactive Isotopes*, Wiley, 1986.
    NuDat 3.0: [https://www.nndc.bnl.gov/nudat3/](https://www.nndc.bnl.gov/nudat3/)

<a id="ref-11"></a>
11. **AME2020 atomic mass evaluation:**
    M. Wang et al., "The AME2020 atomic mass evaluation,"
    *Chinese Physics C* **45** (2021) 030003.

<a id="ref-12"></a>
12. **Bateman decay equations:**
    H. Bateman, "The solution of a system of differential equations occurring in the
    theory of radio-active transformations," *Proc. Cambridge Phil. Soc.* **15** (1910) 423.

<a id="ref-13"></a>
13. **Energy calibration polynomials:**
    G.F. Knoll, *Radiation Detection and Measurement* (2010), Sec. 11-VII.
    K. Debertin & R.G. Helmer (1988), Sec. 2.2.

<a id="ref-14"></a>
14. **S/N pre-screen technique and efficiency model form:**
    J. Pavan & K. Tabor, gnuScope — interactive gamma spectroscopy analysis program,
    Florida State University Nuclear Structure Group. Internal program documentation.

<a id="ref-15"></a>
15. **Fano factor and detector resolution:**
    U. Fano, "Ionization yield of radiations. II. The fluctuations of the number of
    ions," *Physical Review* **72** (1947) 26–29.

<a id="ref-16"></a>
16. **Detection limits and significance in counting experiments:**
    L.A. Currie, "Limits for qualitative detection and quantitative determination:
    application to radiochemistry," *Analytical Chemistry* **40** (1968) 586–593.
    Defines the critical level L_c (false-positive threshold) and detection limit L_D
    from Poisson counting statistics. The N/√(N+2B) significance formula follows
    directly from the Gaussian approximation to paired Poisson counts.

<a id="ref-17"></a>
17. **Experimental nuclear and particle physics techniques:**
    W.R. Leo, *Techniques for Nuclear and Particle Physics Experiments*, 2nd ed.,
    Springer-Verlag, 1994. Chapters 6–7 cover detector resolution, FWHM, peak
    fitting methodology, and signal significance in counting experiments. The
    peak-to-total ratio and detector figures of merit are discussed in Ch. 7.

<a id="ref-18"></a>
18. **β-delayed neutron emission — structure and systematics:**
    B. Pfeiffer, K.-L. Kratz & P. Möller, "Status of delayed-neutron precursor
    data: half-lives and neutron emission probabilities," *Progress in Nuclear
    Energy* **41** (2002) 39–69.
    Provides half-lives and P_n values for β⁻n and β⁻2n emitters relevant to
    fission-fragment decay spectroscopy. The daughter/granddaughter distinction
    for β-delayed neutron emitters follows the chain topology described therein.

<a id="ref-19"></a>
19. **NUBASE2020 nuclear properties evaluation:**
    F.G. Kondev, M. Wang, W.J. Huang, S. Naimi & G. Audi, "The NUBASE2020
    evaluation of nuclear physics properties," *Chinese Physics C* **45** (2021)
    030001.
    Tabulates half-lives, decay modes, excitation energies, and spin-parities
    for all known nuclides. Used by the nuclear database loader for T½ seeding
    in the Decay tab.
