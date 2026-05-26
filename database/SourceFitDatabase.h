#ifndef SOURCEFITDATABASE_H
#define SOURCEFITDATABASE_H

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <limits>
#include <memory>

#include "FitDatabase.h"
#include "Debug.h"

// ─────────────────────────────────────────────────────────────────────────────
// SourceFitDatabase - Separate cache for source spectrum fits
// ─────────────────────────────────────────────────────────────────────────────
// Stores fitted peaks from calibration/source spectra with channel and energy info
// allowing reconstruction of energy vs channel plots.

struct SourcePeakEntry {
    std::string key;  // unique identifier for the source histogram
    
    // Peak energy reference (from database or literature value)
    double energy = 0.0;  // keV
    
    // Fitted parameters
    std::vector<double> params;  // fit parameters (same as FitEntry)
    std::vector<double> paramErrors;
    
    // Fit quality
    double chi2ndf = std::numeric_limits<double>::max();
    double residualRMS = std::numeric_limits<double>::max();
    double maxPull = std::numeric_limits<double>::max();
    
    // Gaussian parameters extracted from fit
    double centroidChannel = 0.0;  // channel position of peak
    double centroidChannelErr = 0.0;
    double sigma = 0.0;            // peak width in channels
    double sigmaErr = 0.0;
    double area = 0.0;             // integrated area
    double areaErr = 0.0;
    
    // Source info
    std::string sourceHistogram;   // which histogram this came from
    std::string isotope;           // isotope ID (e.g., "Co-60")
    double activity = 0.0;         // Bq at measurement time
    std::string fitDate;           // when the fit was performed
    
    // Fit method and window
    std::string fitMethod;
    double xlo = 0.0, xhi = 0.0;   // fit window in keV
};

class SourceFitDatabase {
public:
    SourceFitDatabase() = default;
    ~SourceFitDatabase() = default;
    
    // Enable/disable use of cached seeding
    bool useCachedSeeds = true;
    
    // Load from file: fit_caches/source_fits.txt
    bool Load(const std::string& cachePath = "fit_caches/source_fits.txt");
    
    // Save to file
    bool Save(const std::string& cachePath = "fit_caches/source_fits.txt") const;
    
    // Add a peak entry
    void AddEntry(const SourcePeakEntry& entry);
    
    // Get entry by source histogram name and energy
    const SourcePeakEntry* GetEntry(const std::string& sourceHist, double energy, 
                                    double toleranceKeV = 1.0) const;
    
    // Get all entries for a particular source histogram
    std::vector<const SourcePeakEntry*> GetEntriesForSource(
        const std::string& sourceHist) const;
    
    // Get all entries for an isotope
    std::vector<const SourcePeakEntry*> GetEntriesForIsotope(
        const std::string& isotope) const;
    
    // Clear database
    void Clear() { entries_.clear(); }
    
    // Access all entries
    const std::vector<SourcePeakEntry>& GetAllEntries() const { return entries_; }
    std::vector<SourcePeakEntry>& GetAllEntries() { return entries_; }
    
    // Get collection of unique source histograms
    std::vector<std::string> GetUniqueSourceHistograms() const;
    
    // Get collection of unique isotopes
    std::vector<std::string> GetUniqueIsotopes() const;
    
private:
    std::vector<SourcePeakEntry> entries_;
    
    // Parsing utilities
    bool ParseLine(const std::string& line, SourcePeakEntry& entry);
    std::string SerializeEntry(const SourcePeakEntry& entry) const;
};

#endif // SOURCEFITDATABASE_H
