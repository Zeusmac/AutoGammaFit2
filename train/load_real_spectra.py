"""
load_real_spectra.py — Extract ML training windows from user-approved ROOT spectra.

Reads ml/approved_spectra.txt, loads each histogram via PyROOT, reconstructs the
true background from fit_caches/ fitted parameters, and extracts sliding windows.

Windows are appended to data/synthetic.npz (which must already exist from
generate_spectra.py).  The combined file is saved as data/combined.npz.

Auto-skip rules (configurable at top of file):
  - More than MAX_NEEDS_REFIT_FRAC of cache entries have needsRefit=1
  - Any fit entry has chi2/ndf > MAX_CHI2NDF
  - Fewer than MIN_BINS non-zero bins after background subtraction

Usage:
    python train/load_real_spectra.py [--approved ml/approved_spectra.txt]
"""

import argparse
import os
import sys
import numpy as np

# ── quality gates ─────────────────────────────────────────────────────────────
MAX_NEEDS_REFIT_FRAC = 0.20   # skip spectrum if > 20 % of entries need refit
MAX_CHI2NDF          = 5.0    # skip spectrum if any chi2/ndf exceeds this
MIN_BINS             = 100    # skip spectrum if too few non-zero bins after BG sub

BG_WIN  = 101
BG_HALF = BG_WIN // 2
PK_WIN  = 51
PK_HALF = PK_WIN // 2


# ── fit cache parser ──────────────────────────────────────────────────────────
def parse_fit_cache(cache_path):
    """
    Return list of dicts with keys:
      energy, sigma, amplitude, chi2ndf, needsRefit
    Ignores entries where needsRefit==1 for BG truth (they have unreliable peak params).
    Returns (entries, needs_refit_frac).
    """
    entries = []
    try:
        with open(cache_path) as f:
            content = f.read()
    except FileNotFoundError:
        return entries, 0.0

    for block in content.split("---"):
        block = block.strip()
        if not block:
            continue
        d = {}
        for line in block.splitlines():
            line = line.strip()
            if "=" in line:
                k, _, v = line.partition("=")
                d[k.strip()] = v.strip()
        try:
            d["energy"]     = float(d.get("energy", "nan"))
            d["sigma"]      = float(d.get("sigma",  "nan"))
            d["amplitude"]  = float(d.get("amplitude", "nan"))
            d["chi2ndf"]    = float(d.get("chi2ndf", "0"))
            d["needsRefit"] = int(d.get("needsRefit", "0"))
            entries.append(d)
        except (ValueError, KeyError):
            pass

    n_refit = sum(1 for e in entries if e["needsRefit"])
    frac    = n_refit / len(entries) if entries else 0.0
    return entries, frac


# ── window extraction (same logic as generate_spectra.py) ─────────────────────
def extract_bg_windows(raw_norm, bg_norm, nBins, exclude_mask):
    """Yield (window_101, bg_label) skipping bins covered by a known peak."""
    padded = np.pad(raw_norm, BG_HALF, mode='constant')
    for b in range(nBins):
        if exclude_mask[b]:
            continue
        win   = padded[b: b + BG_WIN].astype(np.float32)
        label = np.float32(bg_norm[b])
        yield win, label


def extract_peak_windows(sub_norm, peak_labels, nBins):
    padded = np.pad(sub_norm, PK_HALF, mode='constant')
    for b in range(nBins):
        win   = padded[b: b + PK_WIN].astype(np.float32)
        label = peak_labels[b]
        yield win, label


# ── process one histogram ─────────────────────────────────────────────────────
def process_histogram(h, cache_path):
    """
    Returns (bg_windows, bg_labels, peak_windows, peak_labels) or None on skip.
    h is a ROOT TH1 object.
    """
    nBins = h.GetNbinsX()
    E     = np.array([h.GetBinCenter(b) for b in range(1, nBins + 1)])

    entries, refit_frac = parse_fit_cache(cache_path)
    if refit_frac > MAX_NEEDS_REFIT_FRAC:
        print(f"    SKIP: refit frac {refit_frac:.2f} > {MAX_NEEDS_REFIT_FRAC}")
        return None
    if not entries:
        print("    SKIP: no cache entries")
        return None

    bad_chi2 = [e for e in entries if e["chi2ndf"] > MAX_CHI2NDF]
    if bad_chi2:
        print(f"    SKIP: {len(bad_chi2)} entries with chi2/ndf > {MAX_CHI2NDF}")
        return None

    # Good entries only (needsRefit==0)
    good = [e for e in entries if e["needsRefit"] == 0]

    # ── Reconstruct background truth ──────────────────────────────────────────
    # bg[b] = spectrum[b] − Σ_i A_i · Gauss(E[b]; E_i, σ_i)
    observed = np.array([h.GetBinContent(b) for b in range(1, nBins + 1)],
                        dtype=np.float64)
    peaks_sum = np.zeros(nBins, dtype=np.float64)
    for e in good:
        Ei  = e["energy"]
        sig = e["sigma"]
        Ai  = e["amplitude"]
        if sig <= 0:
            continue
        peaks_sum += Ai * np.exp(-0.5 * ((E - Ei) / sig) ** 2)

    bg = np.maximum(observed - peaks_sum, 0.0)

    # Exclude bins within 3σ of any peak (unreliable BG truth there)
    exclude = np.zeros(nBins, dtype=bool)
    for e in good:
        Ei  = e["energy"]
        sig = e["sigma"]
        if sig <= 0:
            continue
        exclude |= np.abs(E - Ei) < 3.0 * sig

    # ── Quality check ─────────────────────────────────────────────────────────
    sub = np.maximum(observed - bg, 0.0)
    if sub.sum() < MIN_BINS:
        print(f"    SKIP: too few non-zero bins after BG sub ({sub.sum():.0f})")
        return None

    # ── Normalise ─────────────────────────────────────────────────────────────
    raw_max  = max(observed.max(), 1.0)
    raw_norm = observed / raw_max
    bg_norm  = bg       / raw_max

    sub_max  = max(sub.max(), 1.0)
    sub_norm = sub / sub_max

    # Peak labels — bin is 1 if it holds a known peak centre
    peak_labels = np.zeros(nBins, dtype=np.float32)
    for e in good:
        idx = np.searchsorted(E, e["energy"])
        if 0 <= idx < nBins:
            peak_labels[idx] = 1.0

    bg_X, bg_y     = [], []
    peak_X, peak_y = [], []

    for win, lbl in extract_bg_windows(raw_norm, bg_norm, nBins, exclude):
        bg_X.append(win); bg_y.append(lbl)
    for win, lbl in extract_peak_windows(sub_norm, peak_labels, nBins):
        peak_X.append(win); peak_y.append(lbl)

    return (np.array(bg_X,   dtype=np.float32),
            np.array(bg_y,   dtype=np.float32),
            np.array(peak_X, dtype=np.float32),
            np.array(peak_y, dtype=np.float32))


# ── main ──────────────────────────────────────────────────────────────────────
def main():
    parser = argparse.ArgumentParser(
        description="Append real-spectrum windows to data/synthetic.npz")
    parser.add_argument("--approved",  default="ml/approved_spectra.txt")
    parser.add_argument("--caches",    default="fit_caches",
                        help="Directory containing fit cache .dat files")
    parser.add_argument("--synthetic", default="data/synthetic.npz",
                        help="Existing synthetic dataset to extend")
    parser.add_argument("--out",       default="data/combined.npz")
    args = parser.parse_args()

    # Load approved list
    try:
        with open(args.approved) as f:
            lines = [l.strip() for l in f if l.strip() and not l.startswith("#")]
    except FileNotFoundError:
        print(f"ERROR: {args.approved} not found. Approve spectra via the GUI first.")
        sys.exit(1)

    if not lines:
        print("No approved spectra. Nothing to do.")
        sys.exit(0)

    # Import ROOT (PyROOT from system installation)
    try:
        import ROOT          # noqa: F401
        ROOT.gROOT.SetBatch(True)
        ROOT.gErrorIgnoreLevel = ROOT.kError
    except ImportError:
        print("ERROR: PyROOT not available. Source thisroot.sh before running.")
        sys.exit(1)

    # Load existing synthetic data
    if not os.path.exists(args.synthetic):
        print(f"ERROR: {args.synthetic} not found — run generate_spectra.py first.")
        sys.exit(1)
    syn = np.load(args.synthetic)
    bg_X    = [syn["X_bg"]]
    bg_y    = [syn["y_bg"]]
    peak_X  = [syn["X_peak"]]
    peak_y  = [syn["y_peak"]]

    n_ok = 0
    open_files = {}   # cache open TFiles

    for line in lines:
        parts = line.split()
        if len(parts) < 2:
            print(f"  WARN: malformed line: {line!r}")
            continue
        hname   = parts[0]
        fpath   = parts[1]

        print(f"  Processing {hname}  ({fpath})")

        if fpath not in open_files:
            tf = ROOT.TFile.Open(fpath, "READ")
            if not tf or tf.IsZombie():
                print(f"    SKIP: cannot open {fpath}")
                continue
            open_files[fpath] = tf

        tf = open_files[fpath]
        h  = tf.Get(hname)
        if not h:
            print(f"    SKIP: histogram {hname!r} not found in {fpath}")
            continue

        # Derive cache path: fit_caches/<hname>.dat
        cache_path = os.path.join(args.caches, hname + ".dat")
        result = process_histogram(h, cache_path)
        if result is None:
            continue

        bX, by, pX, py = result
        bg_X.append(bX);   bg_y.append(by)
        peak_X.append(pX); peak_y.append(py)
        print(f"    OK  {len(bX):,} BG windows, {len(pX):,} peak windows "
              f"(peak frac {py.mean():.3f})")
        n_ok += 1

    for tf in open_files.values():
        tf.Close()

    # Concatenate and save
    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    X_bg    = np.concatenate(bg_X,   axis=0)
    y_bg    = np.concatenate(bg_y,   axis=0)
    X_peak  = np.concatenate(peak_X, axis=0)
    y_peak  = np.concatenate(peak_y, axis=0)

    np.savez_compressed(args.out,
                        X_bg=X_bg, y_bg=y_bg,
                        X_peak=X_peak, y_peak=y_peak)
    print(f"\nDone: {n_ok} real spectra added → {args.out}")
    print(f"  BG windows:   {len(X_bg):>10,}")
    print(f"  Peak windows: {len(X_peak):>10,}  (peak frac {y_peak.mean():.4f})")


if __name__ == "__main__":
    main()
