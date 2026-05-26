#ifndef ENERGYCALIBRATION_H
#define ENERGYCALIBRATION_H

#include <map>
#include <string>
#include <vector>
#include <memory>

#include "TH1.h"
#include "TGraph.h"
#include "CalibrationModel.h"

// ─────────────────────────────────────────────────────────────────────────────
// EnergyCalibration - Manager for energy calibrations and their application
// ─────────────────────────────────────────────────────────────────────────────
class EnergyCalibration {
public:
    EnergyCalibration() = default;
    ~EnergyCalibration();
    
    // Store a named calibration
    void AddCalibration(const std::string& name, 
                       const LinearCalibrationModel& model);
    
    // Get named calibration
    const LinearCalibrationModel* GetCalibration(const std::string& name) const;
    LinearCalibrationModel* GetCalibration(const std::string& name);
    
    // Remove calibration
    void RemoveCalibration(const std::string& name) {
        calibrations_.erase(name);
    }
    
    // Clear all calibrations
    void Clear() { calibrations_.clear(); }
    
    // Get all calibration names
    std::vector<std::string> GetCalibrationNames() const;
    
    // Apply calibration to histogram: creates new histogram with energy bins
    // Returns new histogram (caller owns it) or nullptr on failure
    TH1* ApplyCalibrationToHistogram(TH1* hist, 
                                     const std::string& calibName,
                                     const std::string& newHistName = "");
    
    // Apply calibration to histogram in-place (rebins with energy axis)
    // Returns true on success
    bool ApplyCalibrationInPlace(TH1*& hist, 
                                 const std::string& calibName);
    
    // Convert channel array to energy array using calibration
    std::vector<double> ChannelsToEnergies(
        const std::vector<double>& channels,
        const std::string& calibName) const;
    
    // Load all calibrations from file
    bool LoadFromFile(const std::string& dirPath);
    
    // Save all calibrations to file
    bool SaveToFile(const std::string& dirPath) const;
    
    // Get calibration statistics
    int GetCalibrationCount() const { return calibrations_.size(); }
    
private:
    std::map<std::string, LinearCalibrationModel> calibrations_;
    
    // Helper: create energy bin edges from channel bin edges
    std::vector<double> ComputeEnergyBinEdges(
        const TH1* hist,
        const LinearCalibrationModel& calib) const;
};

#endif // ENERGYCALIBRATION_H
