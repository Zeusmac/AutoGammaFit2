#ifndef PEAKFITTER_H
#define PEAKFITTER_H

#include "GammaDB.h"
#include "PeakTracker.h"
#include "ResolutionModel.h"
#include "FitStorage.h"
#include "FitDatabase.h"

#include "TH1.h"
#include "TFile.h"

#include <string>
#include <vector>

class TCanvas;

struct PeakFitterBgOptions {
    bool   subtractBg      = true;
    int    iterations      = 14;
    double tspecSigma      = 2.0;    // TSpectrum::Search sigma (bins)
    double tspecThresh     = 0.02;   // TSpectrum::Search threshold (fraction of max)
    bool   useLogLikelihood = true;  // AdaptiveFitter: "L" option
    bool   useImprove       = false; // AdaptiveFitter: "M" option (IMPROVE)
    double snMinRatio      = 0.0;   // S/N pre-screen: skip groups below threshold (0 = disabled)
    // ML peak finder (requires HAS_ONNX build flag + trained models)
    bool        useML               = false;
    float       mlPeakThresh        = 0.5f;
    std::string mlBgModelPath       = "ml/bg_model.onnx";
    std::string mlPeakModelPath     = "ml/peak_model.onnx";
    std::string mlCentroidModelPath = "ml/centroid_model.onnx";
    // ML sigma constraint — constrain MIGRAD sigma to ±8% of model prediction
    bool        useMLSigma          = false;
    std::string mlSigmaModelPath    = "ml/sigma_model.onnx";
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

    // Clear cached ML model instances so the next FitHistogram() call re-loads
    // them from disk.  Call from the GUI "Reload ML Models" button.
    static void ReloadMLModels();

    // Run just the background MLP and return a background TH1* owned by the
    // caller (nullptr if ONNX not compiled in or model fails to load).
    static TH1* GetMLBackground(TH1* raw,
                                const std::string& modelPath = "ml/bg_model.onnx");

    // Run just the peak-detector MLP on a background-subtracted histogram and
    // return a TH1* where bin content = peak probability [0,1].  Caller owns.
    static TH1* GetMLPeakProbabilities(TH1* bgSub,
                                       const std::string& modelPath = "ml/peak_model.onnx");

    // Run the centroid-regression MLP on a 51-bin window around E and return
    // the sub-bin centroid offset in keV.  Returns NaN when unavailable.
    static double GetMLCentroidOffset(TH1* h, double E, double sig,
                                      const std::string& modelPath = "ml/centroid_model.onnx");

    // Run the sigma-regression MLP on a 51-bin window around E and return
    // the predicted sigma in keV.  Returns NaN when unavailable.
    static double GetMLSigmaWidth(TH1* h, double E,
                                  const std::string& modelPath = "ml/sigma_model.onnx");

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