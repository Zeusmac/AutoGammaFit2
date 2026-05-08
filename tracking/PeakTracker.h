#ifndef PEAKTRACKER_H
#define PEAKTRACKER_H

#include <map>
#include <vector>
#include <string>
#include "ResolutionModel.h"

class TFile;

class ResolutionModel;   // forward declare

class PeakTracker {
public:

    // -----------------------------
    // CORE INTERFACE
    // -----------------------------
    void Add(double energy,
             double time,
             double counts_val,
             double error,
             double sigma,
             double sigma_err,
             int    fit_status = 0,
             double fit_edm    = 0.0);

    // Resolution feedback loop
    void FitResolutionModel(ResolutionModel& model, TFile* fout, const std::string& tag);
    
    // -----------------------------
    // INTERNAL STORAGE
    // -----------------------------
    struct PeakPoint {
        double energy;
        double time;
        double counts;
        double error;
        double sigma;
        double sigma_err;
        int    fit_status = 0;
        double fit_edm    = 0.0;
    };

    std::map<double, std::vector<PeakPoint>> data;
    // -----------------------------
    // ACCESSORS
    // -----------------------------
   const std::map<double, std::vector<PeakPoint>>& GetData() const {
    return data;
    }
    // -----------------------------
    // CONFIGURATION
    // -----------------------------
    void SetEnergyTolerance(double tol);

private:


    // -----------------------------
    // MATCHING / GROUPING
    // -----------------------------
    double GetKey(double energy, double tol_override = -1);
    double energy_tol = 1.0;  // keV grouping tolerance
};

#endif