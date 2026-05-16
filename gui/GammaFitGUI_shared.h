#ifndef GAMMAFITGUI_SHARED_H
#define GAMMAFITGUI_SHARED_H

#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include "TH1.h"
#include "TMath.h"

// ─────────────────────────────────────────────────────────────────────────────
// File-scope constants shared across GammaFitGUI translation units
// ─────────────────────────────────────────────────────────────────────────────

static constexpr const char* kCacheDir         = "fit_caches";
static constexpr const char* kIsotopeDBDefault = "../Isotope_energys.txt";
static constexpr const char* kResolutionKey    = "__RESOLUTION__";
static constexpr const char* kExcludedFwhmKey  = "__EXCLUDED_FWHM__";
static constexpr const char* kFwhmPrefix       = "FWHM: ";

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

// Set histogram y-range: minimum always 0, maximum = highest visible bin × margin.
inline void SetYMaxFromVisible(TH1* h, double margin = 1.10) {
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

inline std::string Fmt(double v, int n = 3) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(n) << v;
    return ss.str();
}

// Format energy with uncertainty in NNDC parenthesis notation.
// e.g. E=1332.492, err=0.004 → "1332.492(4)"
//      E=661.6,    err=0.5   → "661.6(5)"
//      E=1460.0,   err=2.0   → "1460(2)"
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

// Extended formula builder with optional quadratic BG and Compton step.
// Parameter layout:
//   [3i], [3i+1], [3i+2]  = A_i, E_i, sig_i  for i in [0,n)
//   [3n], [3n+1]           = bg0, bg1          (linear BG, always present)
//   [3n+2]                 = bg2               (only when quadBG=true)
//   [3n+nbg+i]             = step_i            (one Erfc amplitude per peak, when comptonStep=true)
inline int NBgPars(bool quadBG) { return quadBG ? 3 : 2; }
inline int StepParIdx(int n, bool quadBG, int peakI) { return 3*n + NBgPars(quadBG) + peakI; }
inline int NTotalPars(int n, bool quadBG, bool comptonStep) {
    return 3*n + NBgPars(quadBG) + (comptonStep ? n : 0);
}

inline std::string BuildNGaussFormulaEx(int n, bool quadBG = false, bool comptonStep = false) {
    std::string f;
    for (int i = 0; i < n; i++) {
        if (i) f += "+";
        int p = 3*i;
        f += "[" + std::to_string(p)   + "]*exp(-0.5*((x-["
           + std::to_string(p+1) + "])/["
           + std::to_string(p+2) + "])^2)";
    }
    int bg = 3*n;
    f += "+[" + std::to_string(bg) + "]+[" + std::to_string(bg+1) + "]*x";
    if (quadBG)
        f += "+[" + std::to_string(bg+2) + "]*x*x";
    if (comptonStep) {
        int nbg = NBgPars(quadBG);
        for (int i = 0; i < n; i++) {
            int sIdx = 3*n + nbg + i;
            // Compton step: erfc rises steeply to the left of the peak centroid;
            // models photons scattered into the continuum below the full-energy peak.
            f += "+[" + std::to_string(sIdx) + "]*TMath::Erfc((x-["
               + std::to_string(3*i+1) + "])/(["
               + std::to_string(3*i+2) + "]*1.41421356))";
        }
    }
    return f;
}

// ─────────────────────────────────────────────────────────────────────────────
// Classification helpers — used by Isotopes, Decay, and Results tabs
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

#endif // GAMMAFITGUI_SHARED_H
