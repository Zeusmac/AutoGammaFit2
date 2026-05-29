"""
generate_spectra.py — Synthetic gamma spectrum generator for ML training.

Produces sliding-window datasets for:
  - Background MLP    (101-bin windows, label = true background at centre bin)
  - Peak detector     (51-bin windows,  label = 1 if centre bin is within a peak)
  - Centroid regressor(51-bin windows,  label = centroid offset in bins ∈ [-0.5, 0.5])

Output: data/synthetic.npz  containing:
    X_bg        (N_windows, 101) float32
    y_bg        (N_windows,)     float32
    X_peak      (M_windows,  51) float32
    y_peak      (M_windows,)     float32
    X_centroid  (K_windows,  51) float32
    y_centroid  (K_windows,)     float32   centroid offset in bins

Usage:
    python train/generate_spectra.py [--n-spectra 20000] [--bins 4096]
"""

import argparse
import numpy as np
import os

# ── constants ─────────────────────────────────────────────────────────────────
E_MIN   = 100.0    # keV
E_MAX   = 3000.0   # keV

BG_WIN  = 101
BG_HALF = BG_WIN // 2
PK_WIN  = 51
PK_HALF = PK_WIN // 2
CEN_WIN  = 51
CEN_HALF = CEN_WIN // 2


# ── FWHM model: FWHM² = a + b·E ───────────────────────────────────────────────
def make_fwhm_params(rng):
    """Return (a, b) for a typical HPGe-like FWHM² = a + b·E.
    Vary ±30% per spectrum to simulate different detectors."""
    a0, b0 = 1.0, 0.002   # baseline (keV²)
    scale  = rng.uniform(0.7, 1.3)
    return a0 * scale, b0 * scale


def fwhm_at(E_keV, a, b):
    return np.sqrt(np.maximum(a + b * E_keV, 0.01))


def sigma_at(E_keV, a, b):
    return fwhm_at(E_keV, a, b) / 2.3548200450309493


# ── background shapes ─────────────────────────────────────────────────────────
def make_background(E, rng):
    kind = rng.integers(0, 3)          # 0=flat, 1=linear, 2=exponential
    if kind == 0:
        c0 = rng.uniform(10, 2000)
        return np.full_like(E, c0)
    if kind == 1:
        c0 = rng.uniform(50, 2000)
        c1 = rng.uniform(-0.3, 0.3)   # counts/keV
        bg = c0 + c1 * (E - E_MIN)
        return np.maximum(bg, 1.0)
    # exponential
    A   = rng.uniform(500, 5000)
    tau = rng.uniform(300, 1500)       # keV
    return A * np.exp(-E / tau) + 5.0


# ── peak sum ──────────────────────────────────────────────────────────────────
def make_peaks(E, bg, rng, a_fwhm, b_fwhm):
    """Return (peak_spectrum, peak_centres, peak_sigmas)."""
    npeaks   = rng.poisson(5)
    spectrum = np.zeros_like(E)
    centres  = []
    sigmas   = []
    for _ in range(npeaks):
        Ei  = rng.uniform(E_MIN + 50, E_MAX - 50)
        sig = sigma_at(Ei, a_fwhm, b_fwhm)
        bg_at_peak = np.interp(Ei, E, bg)
        Ai  = rng.uniform(2, 500) * bg_at_peak
        spectrum += Ai * np.exp(-0.5 * ((E - Ei) / sig) ** 2)
        centres.append(Ei)
        sigmas.append(sig)
    return spectrum, centres, sigmas


# ── window extraction ─────────────────────────────────────────────────────────
def extract_bg_windows(raw_norm, bg_norm, nBins):
    """Yield (window_101, bg_label) for every bin."""
    padded = np.pad(raw_norm, BG_HALF, mode='constant')
    for b in range(nBins):
        win   = padded[b: b + BG_WIN].astype(np.float32)
        label = bg_norm[b]
        yield win, np.float32(label)


def extract_peak_windows(sub_norm, peak_labels, nBins):
    """Yield (window_51, label) for every bin."""
    padded = np.pad(sub_norm, PK_HALF, mode='constant')
    for b in range(nBins):
        win   = padded[b: b + PK_WIN].astype(np.float32)
        label = peak_labels[b]
        yield win, np.float32(label)


def extract_centroid_windows(sub_norm, centres, sigmas, E_axis, nBins):
    """Yield (window_51, offset_in_bins, sigma_in_bins) for each peak.

    offset_in_bins: fractional offset of true centroid from bin centre, in [-0.5, 0.5].
    sigma_in_bins:  true peak sigma in bins (positive float).
    Windows are normalised so the peak amplitude is 1.
    """
    if len(E_axis) < 2:
        return
    bin_w  = E_axis[1] - E_axis[0]
    padded = np.pad(sub_norm, CEN_HALF, mode='constant')
    for Ei, sig in zip(centres, sigmas):
        idx = int(np.searchsorted(E_axis, Ei))
        idx = np.clip(idx, 0, nBins - 1)
        offset     = float(np.clip((Ei - E_axis[idx]) / bin_w, -0.5, 0.5))
        sigma_bins = float(sig / bin_w)
        win = padded[idx: idx + CEN_WIN].astype(np.float32)
        win_max = win.max()
        if win_max > 0:
            win = win / win_max
        yield win, np.float32(offset), np.float32(sigma_bins)


# ── main ──────────────────────────────────────────────────────────────────────
def generate(n_spectra: int, bins: int, seed: int = 42, win_per_spec: int = 100):
    rng    = np.random.default_rng(seed)
    E      = np.linspace(E_MIN, E_MAX, bins)

    bg_X, bg_y        = [], []
    peak_X, peak_y    = [], []
    cen_X,  cen_y     = [], []
    sig_X,  sig_y     = [], []

    for i in range(n_spectra):
        if i % 500 == 0:
            print(f"  spectrum {i}/{n_spectra}", flush=True)

        a_fwhm, b_fwhm = make_fwhm_params(rng)
        bg              = make_background(E, rng)
        peaks_spec, centres, sigmas = make_peaks(E, bg, rng, a_fwhm, b_fwhm)

        observed = rng.poisson(bg + peaks_spec).astype(np.float64)

        raw_max  = np.maximum(observed.max(), 1.0)
        raw_norm = observed / raw_max
        bg_norm  = bg / raw_max

        # ── BG windows — random sample of win_per_spec centre bins ───────────
        half   = BG_HALF
        valid  = np.arange(half, bins - half)
        chosen = rng.choice(valid, size=min(win_per_spec, len(valid)), replace=False)
        padded = np.pad(raw_norm, half, mode='constant')
        for b in chosen:
            win = padded[b : b + BG_WIN].astype(np.float32)
            bg_X.append(win)
            bg_y.append(np.float32(bg_norm[b - half]))

        # ── Peak windows — random sample of win_per_spec centre bins ─────────
        sub = np.maximum(observed - bg, 0.0)
        sub_max  = np.maximum(sub.max(), 1.0)
        sub_norm = sub / sub_max

        peak_labels = np.zeros(bins, dtype=np.float32)
        for Ei in centres:
            idx = np.searchsorted(E, Ei)
            if 0 <= idx < bins:
                peak_labels[idx] = 1.0

        half   = PK_HALF
        valid  = np.arange(half, bins - half)
        chosen = rng.choice(valid, size=min(win_per_spec, len(valid)), replace=False)
        padded = np.pad(sub_norm, half, mode='constant')
        for b in chosen:
            win = padded[b : b + PK_WIN].astype(np.float32)
            peak_X.append(win)
            peak_y.append(np.float32(peak_labels[b - half]))

        # ── Centroid + sigma windows ──────────────────────────────────────────
        for win, offset, sigma_bins in extract_centroid_windows(sub_norm, centres, sigmas, E, bins):
            cen_X.append(win)
            cen_y.append(offset)
            sig_X.append(win.copy())   # same window, different label
            sig_y.append(sigma_bins)

    return (np.array(bg_X,   dtype=np.float32),
            np.array(bg_y,   dtype=np.float32),
            np.array(peak_X, dtype=np.float32),
            np.array(peak_y, dtype=np.float32),
            np.array(cen_X,  dtype=np.float32),
            np.array(cen_y,  dtype=np.float32),
            np.array(sig_X,  dtype=np.float32),
            np.array(sig_y,  dtype=np.float32))


def main():
    parser = argparse.ArgumentParser(description="Generate synthetic gamma spectra for ML training")
    parser.add_argument("--n-spectra", type=int, default=5000,
                        help="Number of synthetic spectra to generate (default 5000)")
    parser.add_argument("--bins",      type=int, default=4096,
                        help="Number of bins per spectrum (default 4096)")
    parser.add_argument("--seed",      type=int, default=42)
    parser.add_argument("--windows-per-spectrum", type=int, default=100,
                        help="Random BG/peak windows sampled per spectrum (default 100)")
    args = parser.parse_args()

    os.makedirs("data", exist_ok=True)
    outpath = "data/synthetic.npz"

    print(f"Generating {args.n_spectra} spectra × {args.bins} bins, "
          f"{args.windows_per_spectrum} windows/spectrum …", flush=True)
    X_bg, y_bg, X_peak, y_peak, X_cen, y_cen, X_sig, y_sig = generate(
        args.n_spectra, args.bins, args.seed, args.windows_per_spectrum)

    print(f"BG windows:       {len(X_bg):>10,}  (mean BG frac: {y_bg.mean():.4f})")
    print(f"Peak windows:     {len(X_peak):>10,}  (peak frac: {y_peak.mean():.4f})")
    print(f"Centroid windows: {len(X_cen):>10,}  (offset rms: {y_cen.std():.4f} bins)")
    print(f"Sigma windows:    {len(X_sig):>10,}  (sigma mean: {y_sig.mean():.4f} bins)")

    np.savez_compressed(outpath,
                        X_bg=X_bg,       y_bg=y_bg,
                        X_peak=X_peak,   y_peak=y_peak,
                        X_centroid=X_cen, y_centroid=y_cen,
                        X_sigma=X_sig,   y_sigma=y_sig)
    print(f"Saved → {outpath}")


if __name__ == "__main__":
    main()
