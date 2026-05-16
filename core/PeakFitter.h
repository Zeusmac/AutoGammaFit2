#ifndef PEAKFITTER_H
#define PEAKFITTER_H

#include "GammaDB.h"
#include "PeakTracker.h"
#include "ResolutionModel.h"
#include "FitStorage.h"
#include "FitDatabase.h"

#include "TH1.h"
#include "TFile.h"

class TCanvas;

struct PeakFitterBgOptions {
    bool   subtractBg      = true;
    int    iterations      = 14;
    double tspecSigma      = 2.0;    // TSpectrum::Search sigma (bins)
    double tspecThresh     = 0.02;   // TSpectrum::Search threshold (fraction of max)
    bool   useLogLikelihood = true;  // AdaptiveFitter: "L" option
    bool   useImprove       = false; // AdaptiveFitter: "M" option (IMPROVE)
    double snMinRatio      = 0.0;   // S/N pre-screen: skip groups below threshold (0 = disabled)
};

// Flat resolution model used as a fallback when the energy-calibrated model
// gives nonsensical sigma for non-energy axes (e.g. time projections).
struct ConstantResModel {
    double fwhm = 5.0;
    double FWHM(double /*E*/) const { return fwhm; }
    double Sigma(double /*E*/) const { return fwhm / 2.3548200450309493; }
};

class PeakFitter {
public:
    using BgOptions = PeakFitterBgOptions;

    PeakFitter(GammaDB& db,
               PeakTracker* tracker,
               ResolutionModel& res,
               FitStorage& storage,
               FitDatabase* fitdb = nullptr);

    void FitHistogram(TH1* h, TFile* fout, bool enableTracking,
                      TCanvas* extCanvas = nullptr,
                      BgOptions bg = BgOptions{},
                      const std::vector<double>& forcedSeeds = {});

    void SetFitDatabase(FitDatabase* db) { fitdb = db; }

private:
    GammaDB& db;
    PeakTracker* tracker;
    ResolutionModel& res;
    FitStorage& storage;
    FitDatabase* fitdb;

    bool IsTooWide(double sigma, double E);
    TF1* BuildModel(int nPeaks, double x);
};

#endif