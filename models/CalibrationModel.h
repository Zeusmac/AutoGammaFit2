#ifndef CALIBRATIONMODEL_H
#define CALIBRATIONMODEL_H

#include <vector>
#include <string>
#include <memory>
#include <cmath>
#include <limits>

#include "TF1.h"
#include "TGraphErrors.h"

// ─────────────────────────────────────────────────────────────────────────────
// CalibrationPoint - Single energy vs channel calibration point
// ─────────────────────────────────────────────────────────────────────────────
struct CalibrationPoint {
    double energy = 0.0;        // keV (reference value from source or standard)
    double energyErr = 0.0;     // uncertainty in energy
    double channel = 0.0;       // channel position (from fitted peak)
    double channelErr = 0.0;    // uncertainty in channel position
    
    std::string sourceHistogram;  // which source measurement this came from
    std::string isotope;          // isotope ID
    bool   enabled = true;        // include in fit if true
};

// ─────────────────────────────────────────────────────────────────────────────
// LinearCalibrationModel - E = a + b*ch linear energy calibration
// ─────────────────────────────────────────────────────────────────────────────
class LinearCalibrationModel {
public:
    LinearCalibrationModel() = default;
    
    // Fit linear calibration: E = a + b*ch from calibration points
    // Returns true on success, sets paramA and paramB
    bool FitLinear(const std::vector<CalibrationPoint>& points,
                   double& paramA, double& paramB,
                   double& paramAErr, double& paramBErr,
                   double& chi2ndf);
    
    // Apply calibration: convert channel to energy
    double ChannelToEnergy(double channel) const {
        if (paramB == 0.0) return paramA;  // prevent divide by zero
        return paramA + paramB * channel;
    }
    
    // Reverse: convert energy to channel
    double EnergyToChannel(double energy) const {
        if (std::abs(paramB) < 1e-10) return 0.0;
        return (energy - paramA) / paramB;
    }
    
    // Store calibration to file
    bool SaveToFile(const std::string& filepath) const;
    
    // Load calibration from file
    bool LoadFromFile(const std::string& filepath);
    
    // Get fitted function
    TF1* GetFitFunction(const std::string& name = "calib") const;
    
    // Get residuals (energy_measured - energy_predicted)
    std::vector<double> GetResiduals(const std::vector<CalibrationPoint>& points) const;
    
    // Properties
    double GetParamA() const { return paramA; }
    double GetParamB() const { return paramB; }
    double GetParamAErr() const { return paramAErr; }
    double GetParamBErr() const { return paramBErr; }
    double GetChi2NDF() const { return chi2ndf; }
    
    void SetParams(double a, double b, double aErr = 0.0, double bErr = 0.0) {
        paramA = a; paramB = b; paramAErr = aErr; paramBErr = bErr;
    }
    
    bool IsValid() const { 
        return std::abs(paramB) > 1e-10;  // b must be non-zero
    }
    
private:
    double paramA = 0.0;      // offset
    double paramB = 1.0;      // slope (keV per channel)
    double paramAErr = 0.0;   // error on A
    double paramBErr = 0.0;   // error on B
    double chi2ndf = -1.0;    // fit quality metric
};

// ─────────────────────────────────────────────────────────────────────────────
// CalibrationBuilder - Aggregator for building and managing calibrations
// ─────────────────────────────────────────────────────────────────────────────
class CalibrationBuilder {
public:
    CalibrationBuilder() = default;
    
    // Add a calibration point
    void AddPoint(const CalibrationPoint& point) {
        points_.push_back(point);
    }
    
    // Add multiple points at once
    void AddPoints(const std::vector<CalibrationPoint>& points) {
        points_.insert(points_.end(), points.begin(), points.end());
    }
    
    // Remove a point by index
    void RemovePoint(size_t index) {
        if (index < points_.size()) {
            points_.erase(points_.begin() + index);
        }
    }
    
    // Clear all points
    void Clear() { 
        points_.clear();
        model_ = LinearCalibrationModel();
    }
    
    // Fit the calibration model
    bool Fit();
    
    // Get the fitted model
    const LinearCalibrationModel& GetModel() const { return model_; }
    LinearCalibrationModel& GetModel() { return model_; }
    
    // Access points
    const std::vector<CalibrationPoint>& GetPoints() const { return points_; }
    std::vector<CalibrationPoint>& GetPoints() { return points_; }
    
    // Get point count
    size_t GetPointCount() const { return points_.size(); }
    size_t GetEnabledPointCount() const;
    
    // Toggle point enabled state
    void SetPointEnabled(size_t index, bool enabled) {
        if (index < points_.size()) {
            points_[index].enabled = enabled;
        }
    }
    
    // Get TGraph for visualization
    TGraphErrors* GetGraph(const std::string& name = "calib_graph") const;
    
    // Save builder state (points and fitted model)
    bool SaveToFile(const std::string& filepath) const;
    
    // Load builder state
    bool LoadFromFile(const std::string& filepath);
    
private:
    std::vector<CalibrationPoint> points_;
    LinearCalibrationModel model_;
};

#endif // CALIBRATIONMODEL_H
