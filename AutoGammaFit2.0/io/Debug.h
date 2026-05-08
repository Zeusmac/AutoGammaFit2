#ifndef DEBUG_H
#define DEBUG_H

#include <string>
#include <vector>

class TH1;
class TF1;
class TFile;

namespace Debug {

// ─────────────────────────────────────────────────────────────────────────────
// Per-section debug toggles — all OFF by default.
// Enable only the sections you care about to avoid drowning in output.
//
//  FITTER     — AdaptiveFitter: param seeding, AIC selection, retry stages
//  GROUPER    — PeakGrouper: which peaks are merged or split
//  TRACKER    — PeakTracker: point filtering, FWHM fit, model update
//  DB         — FitDatabase: cache load / save / seed / store
//  GAMMADB    — GammaDB: line loading, match results, isotope scoring
//  RESMODEL   — ResolutionModel: parameter updates and exponential smoothing
//  FILEIO     — RootFileManager: every object written to disk
//  PEAKFITTER — PeakFitter main loop: per-peak results, SNR, background
//
// CLI: ./gamma_fit input.root --debug=FITTER,DB
//      ./gamma_fit input.root --debug-all
// ─────────────────────────────────────────────────────────────────────────────
enum Section {
    FITTER     = 0,
    GROUPER    = 1,
    TRACKER    = 2,
    DB         = 3,
    GAMMADB    = 4,
    RESMODEL   = 5,
    FILEIO     = 6,
    PEAKFITTER = 7,
    N_SECTIONS
};

// ─── Toggle control ──────────────────────────────────────────────────────────
void Set(Section s, bool enable);
void SetAll(bool enable);
bool IsEnabled(Section s);

// Enable sections from a comma-separated string, e.g. "FITTER,DB,TRACKER"
void EnableFromString(const std::string& csv);

// Print current toggle state to stdout
void PrintConfig();

// Legacy shims that operate on ALL sections at once
void SetEnabled(bool enable);
bool IsEnabled();

// Section name utilities
const char* SectionName(Section s);
Section     SectionFromName(const std::string& name);  // returns N_SECTIONS if unknown

// ─── Logging ─────────────────────────────────────────────────────────────────
void Log(Section s, const std::string& msg);
void Log(const std::string& msg);  // legacy — routes to FITTER

// ─── Specialised helpers (each guarded by its own section) ───────────────────
void LogGroup(int g, const std::vector<double>& peaks);

void LogPeak(int g, int i, double E, double sigma);

void LogFitAttempt(int nPeaks, double chi2,
                   double improvement, double resolutionRatio, int status);

void LogFitFailure(int nPeaks, int status, double chi2ndf);

void LogResolutionDecision(double groupWidth,
                           double expectedFWHM, int maxPeaks);

void DrawFitComponents(TFile* fout, TH1* h, TF1* model,
                       const std::string& name);

void DumpHistogram(TFile* fout, TH1* h,
                   const std::string& dir = "debug");

} // namespace Debug

#endif
