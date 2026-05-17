# AutoGammaFit 2.0 — User Manual

## Contents
1. [Overview](#1-overview)
2. [Quick Start](#2-quick-start)
3. [AutoFit Tab](#3-autofit-tab)
4. [Manual Fit Tab](#4-manual-fit-tab)
5. [FWHM Tab](#5-fwhm-tab)
6. [Source Tab](#6-source-tab)
7. [Fit Results Tab](#7-fit-results-tab)
8. [Isotopes Tab](#8-isotopes-tab)
9. [Decay Tab](#9-decay-tab)
10. [Fitting Models and Mathematics](#10-fitting-models-and-mathematics)
11. [Algorithms](#11-algorithms)
12. [Cache System](#12-cache-system)
13. [References](#13-references)

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

## 3. AutoFit Tab

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

**Model** (see Section 10.7 for derivation):
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

## 4. Manual Fit Tab

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

### 4.5 Fit Options

| Option | Description |
|--------|-------------|
| Log-likelihood (L) | Minimise the Poisson log-likelihood. Preferred when any bin has fewer than ~10 counts. chi² is biased for low counts because the assumption of Gaussian bin errors breaks down. |
| IMPROVE (M) | Run IMPROVE after MIGRAD to search for a better global minimum. Adds ~2× fit time. Use when chi²/ndf is suspiciously high or parameters have hit bounds. |
| MINOS errors (E) | Compute asymmetric (profile-likelihood) parameter uncertainties. More accurate than symmetric MIGRAD errors near parameter boundaries. Significantly slower. |
| Show fit components | Overlay the background (green dashed) and individual Gaussian components (blue dashed) on the spectrum, in addition to the total fit (red). |
| Show isotope matches | Display the matched isotope name from the database above each peak energy label. Uncheck to show energies only. |

### 4.6 Fit Parameter Bounds

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

### 4.7 Peak Statistics Panel

After each fit, the panel displays for every Gaussian component:

| Quantity | Formula | Notes |
|----------|---------|-------|
| Amplitude | A ± σ_A | Raw fit parameter |
| Centroid | E ± σ_E in NNDC notation | e.g. 1332.492(4) = 1332.492 ± 0.004 keV |
| Sigma | σ ± σ_σ | Gaussian width |
| FWHM | 2.3548·σ ± 2.3548·σ_σ | Full width at half maximum |
| Peak area | A·σ·√(2π) / Δx | Δx = bin width; converts counts/bin to counts |
| Area uncertainty | area·√((σ_A/A)² + (σ_σ/σ)²) | Propagated in quadrature |
| SNR | N_peak / √N_BG | Signal-to-noise; N_BG = BG under ±2.5σ |
| P/T (window) | N_peak / N_window | Peak prominence in local context |
| P/T (all) | N_peak / N_spectrum | Detector figure of merit |
| Eff-corr area | N_peak / ε(E) | Only if efficiency model is loaded |
| chi²/ndf | Σ(pull²) / DOF | Computed from residuals, not from ROOT r->Chi2() |
| p-value | TMath::Prob(chi², ndf) | Probability of chi² ≥ observed under H₀ |

The area ratio (fit/obs) compares the total fitted Gaussian area to the observed
counts above the background. A ratio of 0.9–1.1 indicates a good fit.

### 4.8 Residuals

Check **Show residuals** to split the canvas vertically: spectrum + fit on top,
normalised pull residuals (data − fit) / σ_data on the bottom.

**Interpreting residuals:**
- Residuals uniformly scattered in [−2, +2] → good fit
- Systematic arch or bowl shape → wrong background model (try quadratic BG)
- Positive excess to the left, negative at peak centre → Compton step present
- Narrow positive spike between two peaks → missed peak (add another Gaussian)
- All residuals shifted by a constant → bg0 is wrong

### 4.9 Peak Navigation

Use **Prev / Next** or the **Select fit** dropdown to cycle through all cached fits
for the loaded histogram. The canvas zooms to each fit region automatically.
Use **+ / −** or the X range boxes to adjust the zoom level.

---

## 5. FWHM Tab

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

## 6. Source Tab

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
data (see Section 10.7) and enter the coefficients in the Efficiency Correction
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

## 7. Fit Results Tab

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

## 8. Isotopes Tab

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

## 9. Decay Tab

Extracts time-dependent gamma peak intensities from a TH2 histogram (e.g. gamma
energy on one axis, time-since-beam-stop on the other).

**Supported decay models:**

| Model | Equation |
|-------|----------|
| Parent | N(t) = N₀·e^(−λt) |
| Daughter | Bateman: N(t) = [λ_P·N₀_P / (λ_D−λ_P)]·(e^(−λ_P·t) − e^(−λ_D·t)) + N₀_D·e^(−λ_D·t) |
| Granddaughter | Full Bateman A→B→C chain |
| Background | const + A·e^(−λt) |

Half-lives are fitted with uncertainties and stored per peak in a separate
decay cache file. The peak intensity at each time slice is extracted by integrating
the fitted Gaussian within ±N·σ of the centroid.

*Reference: Bateman equations — [Bateman (1910)](#ref-12); ingrowth formalism —
[Debertin & Helmer (1988)](#ref-5), Appendix 1*

---

## 10. Fitting Models and Mathematics

### 10.1 Standard N-Gaussian + Linear Background

The core fit model for N peaks in a window [x_lo, x_hi]:
```
f(x) = Σᵢ Aᵢ · exp(−½·((x − Eᵢ)/σᵢ)²)  +  B₀  +  B₁·x
```
Parameters: Aᵢ (amplitude), Eᵢ (centroid), σᵢ (Gaussian sigma), B₀, B₁.
Total: **3N + 2** free parameters.

The Gaussian models the detector response to a mono-energetic gamma line, which
is broadened by charge-carrier statistics (Fano factor) and electronic noise.

### 10.2 Quadratic Background

Adds a curvature term to the background:
```
BG(x) = B₀ + B₁·x + B₂·x²
```
Total: **3N + 3** parameters. B₂ is seeded at zero.

Use when the Compton continuum has visible curvature across the fit window.

### 10.3 Compton Step Model

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

### 10.4 Tied-Width Constraint

When enabled, each sigma is fixed to the resolution model prediction before fitting:
```
σ_fixed,i = FWHM(Eᵢ) / 2.3548
```
This removes N free parameters from the fit (one per peak), improving convergence
for close doublets and ensuring physical consistency.

*Reference: [Radford (1995)](#ref-4); gnuScope implementation*

### 10.5 Peak Area Calculation

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

### 10.6 Signal-to-Noise Ratio (post-fit)

The post-fit SNR reported in peak statistics:
```
SNR = N_peak / √(N_BG)
```
where N_BG = |BG(Eᵢ)| · 5σᵢ / Δx is the estimated background count under
±2.5σ of the peak. This is the standard figure of merit for peak detectability.

*Reference: [Bevington & Robinson (2003)](#ref-7), Ch. 4; [Debertin & Helmer (1988)](#ref-5), Sec. 4.5*

### 10.7 Detector Resolution Model

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

### 10.8 Efficiency Model

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

### 10.9 NNDC Uncertainty Notation

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

### 10.10 Chi-squared and Goodness of Fit

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

---

## 11. Algorithms

### 11.1 AutoFit Pipeline

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

### 11.2 MIGRAD Minimisation

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

### 11.3 Chi² vs Log-likelihood

**Chi²:** minimises Σ (y_i − f(x_i))² / σᵢ² where σᵢ = √y_i (Poisson).
Valid when y_i >> 1 (typically y_i > 5–10 per bin). Biased for low counts because
the Poisson distribution is not symmetric.

**Poisson log-likelihood (option L):** minimises −2 ln L = 2 Σ [f(x_i) − y_i·ln f(x_i)].
Valid for any count level including zero. Preferred for weak peaks and low-statistics
spectra. The chi²/ndf interpretation of the fit quality is less direct.

*Reference: [Barlow & Beeston (1993)](#ref-9); [Bevington & Robinson (2003)](#ref-7), Ch. 8*

### 11.4 Background Anchor Algorithm

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

### 11.5 AIC Model Selection

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

## 12. Cache System

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

Cache files are plain text — they can be inspected, edited, or merged with a text
editor if needed.

---

## 13. References

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
