#include "Debug.h"

#include <iostream>
#include <sstream>
#include <algorithm>
#include <cctype>

#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TF1.h"

namespace Debug {

// ─── State ───────────────────────────────────────────────────────────────────
static bool sections_[N_SECTIONS] = {};

static const char* kNames_[N_SECTIONS] = {
    "FITTER", "GROUPER", "TRACKER", "DB",
    "GAMMADB", "RESMODEL", "FILEIO", "PEAKFITTER",
    "GUI", "CACHE", "MANUAL"
};

// ─── Toggle control ──────────────────────────────────────────────────────────
void Set(Section s, bool enable) {
    if (s < N_SECTIONS) sections_[s] = enable;
}

void SetAll(bool enable) {
    for (int i = 0; i < N_SECTIONS; ++i) sections_[i] = enable;
}

bool IsEnabled(Section s) {
    return s < N_SECTIONS && sections_[s];
}

bool IsEnabled() {
    for (int i = 0; i < N_SECTIONS; ++i)
        if (sections_[i]) return true;
    return false;
}

void SetEnabled(bool enable) { SetAll(enable); }

const char* SectionName(Section s) {
    return (s < N_SECTIONS) ? kNames_[s] : "UNKNOWN";
}

Section SectionFromName(const std::string& raw) {
    std::string name = raw;
    std::transform(name.begin(), name.end(), name.begin(), ::toupper);
    while (!name.empty() && std::isspace((unsigned char)name.front())) name.erase(name.begin());
    while (!name.empty() && std::isspace((unsigned char)name.back()))  name.pop_back();
    for (int i = 0; i < N_SECTIONS; ++i)
        if (name == kNames_[i]) return static_cast<Section>(i);
    return N_SECTIONS;
}

void EnableFromString(const std::string& csv) {
    std::istringstream ss(csv);
    std::string token;
    while (std::getline(ss, token, ',')) {
        Section s = SectionFromName(token);
        if (s < N_SECTIONS)
            Set(s, true);
        else
            std::cerr << "[Debug] Unknown section name: \"" << token
                      << "\"  valid: FITTER GROUPER TRACKER DB GAMMADB RESMODEL FILEIO PEAKFITTER\n";
    }
}

void PrintConfig() {
    std::cout << "[Debug] Enabled sections:";
    bool any = false;
    for (int i = 0; i < N_SECTIONS; ++i) {
        if (sections_[i]) { std::cout << " " << kNames_[i]; any = true; }
    }
    if (!any) std::cout << " (none — all debug output suppressed)";
    std::cout << "\n";
}

// ─── Logging ─────────────────────────────────────────────────────────────────
void Log(Section s, const std::string& msg) {
    if (!IsEnabled(s)) return;
    std::cout << "[" << kNames_[s] << "] " << msg << "\n";
}

void Log(const std::string& msg) {
    Log(FITTER, msg);
}

// ─── Specialised helpers ─────────────────────────────────────────────────────
void LogGroup(int g, const std::vector<double>& peaks) {
    if (!IsEnabled(GROUPER)) return;
    std::cout << "[GROUPER] Group " << g << "  peaks:";
    for (double p : peaks) std::cout << " " << p;
    std::cout << "\n";
}

void LogPeak(int g, int i, double E, double sigma) {
    if (!IsEnabled(PEAKFITTER)) return;
    std::cout << "[PEAKFITTER] Group " << g << "  Peak " << i
              << "  E=" << E << "  sigma=" << sigma << "\n";
}

void LogFitAttempt(int nPeaks, double chi2, double improvement,
                   double resolutionRatio, int status) {
    if (!IsEnabled(FITTER)) return;
    std::cout << "[FITTER] " << nPeaks << "-peak"
              << "  chi2=" << chi2
              << "  improve=" << improvement
              << "  FWHMratio=" << resolutionRatio
              << "  status=" << status << "\n";
}

void LogFitFailure(int nPeaks, int status, double chi2ndf) {
    if (!IsEnabled(FITTER)) return;
    std::cerr << "[FITTER] WARNING  nPeaks=" << nPeaks
              << "  status=" << status
              << "  chi2/ndf=" << chi2ndf << "\n";
}

void LogResolutionDecision(double groupWidth, double expectedFWHM, int maxPeaks) {
    if (!IsEnabled(GROUPER)) return;
    std::cout << "[GROUPER] groupWidth=" << groupWidth
              << "  expectedFWHM=" << expectedFWHM
              << "  maxPeaks=" << maxPeaks << "\n";
}

void DumpHistogram(TFile* fout, TH1* h, const std::string& dir) {
    if (!IsEnabled(FILEIO) || !fout || !h) return;
    TDirectory* d = fout->GetDirectory(dir.c_str());
    if (!d) d = fout->mkdir(dir.c_str());
    d->cd();
    d->WriteObject(h, h->GetName());
    fout->cd();
    std::cout << "[FILEIO] DumpHistogram: " << h->GetName()
              << "  dir=" << dir << "\n";
}

void DrawFitComponents(TFile* fout, TH1* h, TF1* model, const std::string& name) {
    if (!IsEnabled(FILEIO) || !fout || !model) return;
    TDirectory* d = fout->GetDirectory("debugFits");
    if (!d) d = fout->mkdir("debugFits");
    d->cd();
    model->SetName(name.c_str());
    d->WriteObject(model, name.c_str());
    fout->cd();
    std::cout << "[FILEIO] DrawFitComponents: " << name << "\n";
}

} // namespace Debug
