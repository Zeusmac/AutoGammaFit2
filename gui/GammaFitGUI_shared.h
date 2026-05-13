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

// Set histogram y-max to the highest visible bin × margin, leaving room for labels.
inline void SetYMaxFromVisible(TH1* h, double margin = 1.30) {
    int first = h->GetXaxis()->GetFirst();
    int last  = h->GetXaxis()->GetLast();
    double ymax = 0.0;
    for (int b = first; b <= last; b++) {
        double v = h->GetBinContent(b);
        if (v > ymax) ymax = v;
    }
    h->SetMaximum(ymax > 0.0 ? ymax * margin : 1.0);
}

inline std::string Fmt(double v, int n = 3) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(n) << v;
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
