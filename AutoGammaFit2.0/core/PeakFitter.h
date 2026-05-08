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

class PeakFitter {
public:
    PeakFitter(GammaDB& db,
               PeakTracker* tracker,
               ResolutionModel& res,
               FitStorage& storage,
               FitDatabase* fitdb = nullptr);

    void FitHistogram(TH1* h, TFile* fout, bool enableTracking, TCanvas* extCanvas = nullptr);

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