#ifndef GAMMAFITGUI_SHARED_H
#define GAMMAFITGUI_SHARED_H

#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "TH1.h"
#include "TF1.h"
#include "TMath.h"
#include "TSystem.h"
#include "TGClient.h"
#include "TGFileDialog.h"

// ─────────────────────────────────────────────────────────────────────────────
// File-scope constants shared across GammaFitGUI translation units
// ─────────────────────────────────────────────────────────────────────────────

static constexpr const char* kCacheDir         = "fit_caches";
static constexpr const char* kIsotopeDBDefault = "../Isotope_energys.txt";
static constexpr const char* kResolutionKey    = "__RESOLUTION__";
static constexpr const char* kExcludedFwhmKey  = "__EXCLUDED_FWHM__";
static constexpr const char* kFwhmInclKey      = "__FWHM_INCL__";
static constexpr const char* kEnergyCalKey     = "__ENERGY_CAL__";
static constexpr const char* kFwhmPrefix       = "FWHM: ";
static constexpr const char* kSettingsFile     = "gamma_gui.conf";
static constexpr const char* kResSourceKey     = "__RES_SOURCE__";

// ─────────────────────────────────────────────────────────────────────────────
// Widget IDs  -  one place to own all numeric IDs so collisions surface here
// ─────────────────────────────────────────────────────────────────────────────
enum WidgetID : Int_t {
    kWID_HistList         = 100,
    kWID_HistClassCombo   = 101,
    kWID_RecentCombo      = 102,
    kWID_ManualCombo      = 300,
    kWID_ResidualFitCombo = 301,
    kWID_HistViewCombo    = 700,
    kWID_FWHMCombo        = 800,
    kWID_PeakLabelCombo   = 921,
    kWID_BgSubSrcCombo    = 930,
    kWID_BgSubBgCombo     = 931,
};

// ─────────────────────────────────────────────────────────────────────────────
// FitLayout  -  describes the parameter layout of a fitted TF1.
//
// Supported layouts (where N = number of Gaussians):
//   Standard:                      3N+2   params  (bg0, bg1)
//   QuadBG:                        3N+3   params  (bg0, bg1, bg2)
//   ComptonStep:                   4N+2   params  (bg0, bg1, step_0..N-1)
//   QuadBG + ComptonStep:          4N+3   params  (bg0, bg1, bg2, step_0..N-1)
//   ExpoTail:                      5N+2   params  (T_0..N-1, beta_0..N-1)
//   QuadBG + ExpoTail:             5N+3   params
//   ComptonStep + ExpoTail:        6N+2   params  (step_0..N-1, T_0..N-1, beta_0..N-1)
//   QuadBG + ComptonStep + Expo:   6N+3   params
//
// Special case: Double-Gaussian (DG) model from AdaptiveFitter has 7 params
// but a different formula  -  detected by TryDetectDG(TF1*) below.
// ─────────────────────────────────────────────────────────────────────────────
struct FitLayout {
    int  n           = 0;      // number of Gaussians; 0 = unrecognised
    bool quadBG      = false;
    bool comptonStep = false;
    bool expoTail    = false;  // low-energy exponential tail (HYPERMET)
    bool dg          = false;  // true only when confirmed as DG model via TryDetectDG

    bool valid()    const { return n > 0; }

    // Index of bg0 in the parameter array.
    int  bgBase()   const { return dg ? 5 : 3*n; }

    // Number of background parameters.
    int  nBGPars()  const { return dg ? 2 : (quadBG ? 3 : 2); }

    // Total parameter count (matches NTotalPars for non-DG).
    int  totalPars() const {
        if (dg) return 7;
        int nStep = comptonStep ? n : 0;
        int nTail = expoTail   ? 2*n : 0;
        return 3*n + nBGPars() + nStep + nTail;
    }
};

// Detect layout from parameter count alone.
// Pass 1 covers standard + comptonStep (backward compatible).
// Pass 2 covers expoTail models (higher npar values, tried only if pass 1 fails).
// If the TF1 is available and might be the DG model, call TryDetectDG first.
inline FitLayout DetectLayout(int npar) {
    FitLayout lay;
    for (int nc = 1; nc <= 12; nc++) {
        if      (npar == 3*nc + 2) { lay.n = nc; return lay; }
        else if (npar == 3*nc + 3) { lay.n = nc; lay.quadBG      = true; return lay; }
        else if (npar == 4*nc + 2) { lay.n = nc; lay.comptonStep = true; return lay; }
        else if (npar == 4*nc + 3) { lay.n = nc; lay.quadBG = true; lay.comptonStep = true; return lay; }
    }
    for (int nc = 1; nc <= 12; nc++) {
        if      (npar == 5*nc + 2) { lay.n = nc; lay.expoTail = true; return lay; }
        else if (npar == 5*nc + 3) { lay.n = nc; lay.quadBG = true; lay.expoTail = true; return lay; }
        else if (npar == 6*nc + 2) { lay.n = nc; lay.comptonStep = true; lay.expoTail = true; return lay; }
        else if (npar == 6*nc + 3) { lay.n = nc; lay.quadBG = true; lay.comptonStep = true; lay.expoTail = true; return lay; }
    }
    return lay;  // invalid: valid() == false
}

// Confirm whether f is the AdaptiveFitter Double-Gaussian model.
// If so, populates lay with dg=true, n=1 and returns true.
// Call this before DetectLayout when you have a TF1 in hand.
inline bool TryDetectDG(TF1* f, FitLayout& lay) {
    if (!f || f->GetNpar() != 7) return false;
    std::string formula = f->GetExpFormula().Data();
    if (formula.find("TMath::Erfc") != std::string::npos) return false;
    // DG has exactly 2 exp( terms
    size_t cnt = 0, pos = 0;
    while ((pos = formula.find("exp(", pos)) != std::string::npos) { ++cnt; ++pos; }
    if (cnt != 2) return false;
    lay = FitLayout{};
    lay.n  = 1;
    lay.dg = true;
    return true;
}

// ─────────────────────────────────────────────────────────────────────────────
// FitModel  -  self-describing model for building TF1 formulas.
// ─────────────────────────────────────────────────────────────────────────────
struct FitModel {
    int  n           = 1;
    bool quadBG      = false;
    bool comptonStep = false;
    bool expoTail    = false;

    int bgBase()    const { return 3*n; }
    int nBGPars()   const { return quadBG ? 3 : 2; }
    int stepBase()  const { return 3*n + nBGPars(); }
    int tailBase()  const { return 3*n + nBGPars() + (comptonStep ? n : 0); }
    int totalPars() const {
        int nStep = comptonStep ? n : 0;
        int nTail = expoTail   ? 2*n : 0;
        return 3*n + nBGPars() + nStep + nTail;
    }

    std::string Formula() const;  // implemented via BuildNGaussFormulaEx below

    static FitModel FromLayout(const FitLayout& lay) {
        return {lay.n, lay.quadBG, lay.comptonStep, lay.expoTail};
    }
};

// ─────────────────────────────────────────────────────────────────────────────
// Shared helper functions
// ─────────────────────────────────────────────────────────────────────────────

inline void ApplyHistStyle(TH1* h, const char* title = nullptr) {
    if (title && title[0]) h->SetTitle(title);
    else if (!h->GetTitle() || h->GetTitle()[0] == '\0') h->SetTitle(h->GetName());
    h->GetXaxis()->SetTitle("Energy (keV)");
    h->GetYaxis()->SetTitle("Counts");
    h->SetMinimum(0);
}

// Set histogram y-range: minimum always 0, maximum = highest visible bin x margin.
inline void SetYMaxFromVisible(TH1* h, double margin = 1.30) {
    int first = h->GetXaxis()->GetFirst();
    int last  = h->GetXaxis()->GetLast();
    double ymax = 0.0;
    for (int b = first; b <= last; b++) {
        double v = h->GetBinContent(b);
        if (v > ymax) ymax = v;
    }
    h->SetMinimum(0.0);
    h->SetMaximum(ymax > 0.0 ? ymax * margin : 1.0);
}

// Set Y-axis title based on actual bin width  -  works for both rebinned and
// original histograms with non-unit bins.
inline void SetHistYTitle(TH1* h, bool isDecay = false) {
    double bw = h->GetBinWidth(1);
    if (bw > 1.01) {
        const char* unit = isDecay ? "ms" : "keV";
        h->GetYaxis()->SetTitle(Form("Counts / (%.4g %s)", bw, unit));
    } else {
        h->GetYaxis()->SetTitle("Counts");
    }
}

inline std::string Fmt(double v, int n = 3) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(n) << v;
    return ss.str();
}

// Format energy with uncertainty in NNDC parenthesis notation.
// e.g. E=1332.492, err=0.004 -> "1332.492(4)"
inline std::string NNDCFormat(double E, double err) {
    if (err <= 0.0 || !std::isfinite(err)) return Fmt(E, 1);
    double mag  = std::floor(std::log10(err));
    int    ndec = std::min(4, std::max(0, (int)(-mag)));
    double scale = std::pow(10.0, (double)ndec);
    int    unc   = std::max(1, (int)std::round(err * scale));
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(ndec) << E << "(" << unc << ")";
    return ss.str();
}

// ─── Formula builders ─────────────────────────────────────────────────────────

inline std::string BuildNGaussFormula(int n) {
    std::string f;
    for (int i = 0; i < n; i++) {
        if (i) f += "+";
        int p = 3 * i;
        f += "[" + std::to_string(p)   + "]*exp(-0.5*((x-["
           + std::to_string(p+1) + "])/["
           + std::to_string(p+2) + "])^2)";
    }
    f += "+[" + std::to_string(3*n) + "]+[" + std::to_string(3*n+1) + "]*x";
    return f;
}

// Extended formula builder  -  single source of truth for all TF1 strings.
// Parameter layout:
//   [3i], [3i+1], [3i+2]         = A_i, E_i, sig_i          for i in [0,n)
//   [3n], [3n+1]                  = bg0, bg1                 (always)
//   [3n+2]                        = bg2                      (quadBG only)
//   [3n+nBG+i]                    = step_i                   (comptonStep only)
//   [3n+nBG+nStep+i]              = T_i   (tail amplitude)   (expoTail only)
//   [3n+nBG+nStep+n+i]            = beta_i (tail slope)      (expoTail only)
inline int NBgPars(bool quadBG) { return quadBG ? 3 : 2; }
inline int StepParIdx(int n, bool quadBG, int peakI) { return 3*n + NBgPars(quadBG) + peakI; }
inline int TailAmplParIdx(int n, bool quadBG, bool comptonStep, int peakI) {
    return 3*n + NBgPars(quadBG) + (comptonStep ? n : 0) + peakI;
}
inline int TailSlopeParIdx(int n, bool quadBG, bool comptonStep, int peakI) {
    return 3*n + NBgPars(quadBG) + (comptonStep ? n : 0) + n + peakI;
}
inline int NTotalPars(int n, bool quadBG, bool comptonStep, bool expoTail = false) {
    int nStep = comptonStep ? n : 0;
    int nTail = expoTail   ? 2*n : 0;
    return 3*n + NBgPars(quadBG) + nStep + nTail;
}

inline std::string BuildNGaussFormulaEx(int n, bool quadBG = false,
                                        bool comptonStep = false, bool expoTail = false) {
    std::string f;
    for (int i = 0; i < n; i++) {
        if (i) f += "+";
        int p = 3*i;
        f += "[" + std::to_string(p)   + "]*exp(-0.5*((x-["
           + std::to_string(p+1) + "])/["
           + std::to_string(p+2) + "])^2)";
    }
    int bg = 3*n;
    int nbg = NBgPars(quadBG);
    f += "+[" + std::to_string(bg) + "]+[" + std::to_string(bg+1) + "]*x";
    if (quadBG)
        f += "+[" + std::to_string(bg+2) + "]*x*x";
    if (comptonStep) {
        for (int i = 0; i < n; i++) {
            int sIdx = 3*n + nbg + i;
            // Compton step: erfc models scattered photons below the full-energy peak.
            f += "+[" + std::to_string(sIdx) + "]*TMath::Erfc((x-["
               + std::to_string(3*i+1) + "])/(["
               + std::to_string(3*i+2) + "]*1.41421356))";
        }
    }
    if (expoTail) {
        int nStep = comptonStep ? n : 0;
        for (int i = 0; i < n; i++) {
            int tIdx = 3*n + nbg + nStep + i;
            int bIdx = 3*n + nbg + nStep + n + i;
            // HYPERMET low-energy exponential tail:
            // T/2 * exp((x-E)/β + s²/(2β²)) * erfc((x-E)/(s√2) + s/(β√2))
            // Naturally decays to zero on both sides of the peak centroid.
            auto s = [](int v){ return "[" + std::to_string(v) + "]"; };
            f += "+" + s(tIdx) + "/2*exp((x-" + s(3*i+1) + ")/" + s(bIdx)
               + "+" + s(3*i+2) + "*" + s(3*i+2)
               + "/(2*" + s(bIdx) + "*" + s(bIdx) + "))"
               + "*TMath::Erfc((x-" + s(3*i+1) + ")/(" + s(3*i+2) + "*1.41421356)"
               + "+" + s(3*i+2) + "/(" + s(bIdx) + "*1.41421356))";
        }
    }
    return f;
}

// FitModel::Formula() delegates to BuildNGaussFormulaEx
inline std::string FitModel::Formula() const {
    return BuildNGaussFormulaEx(n, quadBG, comptonStep, expoTail);
}

// ─────────────────────────────────────────────────────────────────────────────
// Classification helpers  -  used by Isotopes, Decay, and Results tabs
// ─────────────────────────────────────────────────────────────────────────────

inline std::string ClassToString(int sel, const std::string& custom)
{
    switch (sel) {
        case 2:  return "Parent";
        case 3:  return "Daughter";
        case 4:  return "Granddaughter";
        case 5:  return "Beta-n Daughter";
        case 6:  return "Beta-2n Daughter";
        case 7:  return "Beta-n Granddaughter";
        case 8:  return "Beta-2n Granddaughter";
        case 9:  return "Background";
        case 10: return custom.empty() ? "Custom" : ("Custom:" + custom);
        case 11: return "X-ray";
        default: return "";
    }
}

inline int ClassToComboIndex(const std::string& cls)
{
    if (cls.empty() || cls == "(none)")       return 1;
    if (cls == "Parent")                       return 2;
    if (cls == "Daughter")                     return 3;
    if (cls == "Granddaughter")                return 4;
    if (cls == "Beta-n Daughter" || cls == "Beta-n")   return 5;
    if (cls == "Beta-2n Daughter" || cls == "Beta-2n") return 6;
    if (cls == "Beta-n Granddaughter")         return 7;
    if (cls == "Beta-2n Granddaughter")        return 8;
    if (cls == "Background")                   return 9;
    if (cls.size() >= 6 && cls.substr(0,6) == "Custom") return 10;
    if (cls == "X-ray")                        return 11;
    return 1;
}

// ─────────────────────────────────────────────────────────────────────────────
// OpenFileDialog  -  drop-in for `new TGFileDialog(...)` that preserves CWD.
//
// ROOT's TGFileDialog calls gSystem->cd() as the user navigates, permanently
// changing the process working directory.  All relative cache paths built by
// CacheFileFor() break after the dialog closes.  This wrapper saves and
// restores the CWD around every dialog so callers are unaffected.
// ─────────────────────────────────────────────────────────────────────────────
inline void OpenFileDialog(TGWindow* main, EFileDialogMode mode, TGFileInfo* fi) {
    std::string cwd = gSystem->WorkingDirectory();
    new TGFileDialog(gClient->GetRoot(), main, mode, fi);
    gSystem->cd(cwd.c_str());
}

#endif // GAMMAFITGUI_SHARED_H
