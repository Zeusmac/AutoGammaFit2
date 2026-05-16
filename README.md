I used Claude Ai to vibe code this. I also told it to make a manual:

# AutoGammaFit 2.0

Interactive ROOT-based GUI for fitting and analyzing gamma-ray spectra. Automated peak finding, manual Gaussian fitting, isotope identification, resolution modeling, energy calibration, efficiency correction, and decay chain visualization in a single window.

---

## Requirements

- ROOT 6.x with GUI support (`-lGui`)
- C++17-capable compiler
- Linux or macOS

---

## Build & Launch

```bash
# Build everything
make

# Build the GUI only
make gui

# Launch (run from project root)
./bin/gamma_gui
```

---

## Features

- **AutoFit** — TSpectrum peak finding + MIGRAD Gaussian fitting with AIC model selection, S/N pre-screening, and batch processing
- **Manual Fit** — click-to-place peaks, quadratic background, Compton step (Erfc), tied widths, BG anchor seeding
- **Source Tab** — calibration source analysis, energy calibration (channel→keV), efficiency vs energy, FWHM from source
- **FWHM Tab** — three-parameter resolution model `FWHM² = a + b·E + c·E²` with Fano limit overlay
- **Decay Tab** — half-life extraction from decay-time histograms
- **Fit Results** — batch export to ROOT file or CSV, canvas annotations
- **Isotopes Tab** — automatic isotope matching, decay schematic with AME/NUBASE data
- **Cache System** — plain-text fit caches with archive/restore support

---

## Documentation

See [MANUAL.md](MANUAL.md) for full usage instructions, algorithm descriptions, and references.

---

## Cache Files

Fit results are stored in `fit_caches/` relative to the working directory. Archived caches live in `fit_caches/archive/`. Cache files survive `make clean`.
