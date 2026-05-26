#include "CalibrationModel.h"
#include "TAxis.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <numeric>
#include <cmath>

// ─────────────────────────────────────────────────────────────────────────────
// LinearCalibrationModel Implementation
// ─────────────────────────────────────────────────────────────────────────────

bool LinearCalibrationModel::FitLinear(const std::vector<CalibrationPoint>& points,
                                       double& paramA, double& paramB,
                                       double& paramAErr, double& paramBErr,
                                       double& chi2ndf) {
    if (points.size() < 2) return false;
    
    // Simple least-squares fit: E = a + b*ch
    // Collect enabled points only
    std::vector<CalibrationPoint> enabledPoints;
    for (const auto& p : points) {
        if (p.enabled) enabledPoints.push_back(p);
    }
    
    if (enabledPoints.size() < 2) return false;
    
    int n = enabledPoints.size();
    double sumX = 0.0, sumY = 0.0, sumXY = 0.0, sumX2 = 0.0, sumY2 = 0.0;
    double sumWX = 0.0, sumWY = 0.0, sumWXY = 0.0, sumWX2 = 0.0, sumW = 0.0;
    
    // Weighted least squares (using 1/err² as weights)
    for (const auto& p : enabledPoints) {
        double w = 1.0;  // default weight
        if (p.energyErr > 0) w = 1.0 / (p.energyErr * p.energyErr);
        
        double x = p.channel;
        double y = p.energy;
        
        sumX += x;
        sumY += y;
        sumXY += x * y;
        sumX2 += x * x;
        sumY2 += y * y;
        
        sumW += w;
        sumWX += w * x;
        sumWY += w * y;
        sumWXY += w * x * y;
        sumWX2 += w * x * x;
    }
    
    // Use weighted fitting
    double denom = sumW * sumWX2 - sumWX * sumWX;
    if (std::abs(denom) < 1e-10) return false;
    
    paramB = (sumW * sumWXY - sumWX * sumWY) / denom;
    paramA = (sumWY - paramB * sumWX) / sumW;
    
    // Compute chi²
    double chi2 = 0.0;
    for (const auto& p : enabledPoints) {
        double w = 1.0;
        if (p.energyErr > 0) w = 1.0 / (p.energyErr * p.energyErr);
        
        double predicted = paramA + paramB * p.channel;
        double residual = p.energy - predicted;
        chi2 += w * residual * residual;
    }
    
    int dof = n - 2;  // 2 parameters (a, b)
    chi2ndf = (dof > 0) ? chi2 / dof : chi2;
    
    // Parameter errors (simplified: from covariance matrix)
    double s2 = chi2 / dof;
    paramAErr = std::sqrt(s2 * sumWX2 / denom);
    paramBErr = std::sqrt(s2 * sumW / denom);
    
    this->paramA = paramA;
    this->paramB = paramB;
    this->paramAErr = paramAErr;
    this->paramBErr = paramBErr;
    this->chi2ndf = chi2ndf;
    
    return true;
}

bool LinearCalibrationModel::SaveToFile(const std::string& filepath) const {
    std::ofstream outfile(filepath);
    if (!outfile.is_open()) return false;
    
    outfile << paramA << " " << paramB << " " << paramAErr << " " << paramBErr << " " << chi2ndf << "\n";
    outfile.close();
    return true;
}

bool LinearCalibrationModel::LoadFromFile(const std::string& filepath) {
    std::ifstream infile(filepath);
    if (!infile.is_open()) return false;
    
    infile >> paramA >> paramB >> paramAErr >> paramBErr >> chi2ndf;
    infile.close();
    return true;
}

TF1* LinearCalibrationModel::GetFitFunction(const std::string& name) const {
    auto* f1 = new TF1(name.c_str(), "[0] + [1]*x", 0, 4096);
    f1->SetParameter(0, paramA);
    f1->SetParameter(1, paramB);
    return f1;
}

std::vector<double> LinearCalibrationModel::GetResiduals(
    const std::vector<CalibrationPoint>& points) const {
    std::vector<double> residuals;
    for (const auto& p : points) {
        if (!p.enabled) continue;
        double predicted = paramA + paramB * p.channel;
        residuals.push_back(p.energy - predicted);
    }
    return residuals;
}

// ─────────────────────────────────────────────────────────────────────────────
// CalibrationBuilder Implementation
// ─────────────────────────────────────────────────────────────────────────────

size_t CalibrationBuilder::GetEnabledPointCount() const {
    return std::count_if(points_.begin(), points_.end(),
                        [](const CalibrationPoint& p) { return p.enabled; });
}

bool CalibrationBuilder::Fit() {
    if (GetEnabledPointCount() < 2) return false;
    
    double a = 0, b = 1, aErr = 0, bErr = 0, chi2ndf = 0;
    if (model_.FitLinear(points_, a, b, aErr, bErr, chi2ndf)) {
        model_.SetParams(a, b, aErr, bErr);
        return true;
    }
    return false;
}

TGraphErrors* CalibrationBuilder::GetGraph(const std::string& name) const {
    auto* graph = new TGraphErrors();
    graph->SetName(name.c_str());
    graph->SetTitle("Energy vs Channel Calibration");
    graph->GetXaxis()->SetTitle("Channel");
    graph->GetYaxis()->SetTitle("Energy (keV)");
    
    int pointIndex = 0;
    for (size_t i = 0; i < points_.size(); i++) {
        if (points_[i].enabled) {
            graph->SetPoint(pointIndex, points_[i].channel, points_[i].energy);
            graph->SetPointError(pointIndex, points_[i].channelErr, points_[i].energyErr);
            pointIndex++;
        }
    }
    
    return graph;
}

bool CalibrationBuilder::SaveToFile(const std::string& filepath) const {
    std::ofstream outfile(filepath);
    if (!outfile.is_open()) return false;
    
    // Save model first
    outfile << "# MODEL\n";
    outfile << model_.GetParamA() << " " << model_.GetParamB() << " ";
    outfile << model_.GetParamAErr() << " " << model_.GetParamBErr() << " ";
    outfile << model_.GetChi2NDF() << "\n";
    
    // Save points
    outfile << "# POINTS\n";
    for (const auto& p : points_) {
        outfile << p.energy << " " << p.energyErr << " ";
        outfile << p.channel << " " << p.channelErr << " ";
        outfile << p.sourceHistogram << " " << p.isotope << " ";
        outfile << (p.enabled ? "1" : "0") << "\n";
    }
    
    outfile.close();
    return true;
}

bool CalibrationBuilder::LoadFromFile(const std::string& filepath) {
    std::ifstream infile(filepath);
    if (!infile.is_open()) return false;
    
    Clear();
    
    std::string line;
    bool readingModel = false;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') {
            if (line.find("MODEL") != std::string::npos) readingModel = true;
            else if (line.find("POINTS") != std::string::npos) readingModel = false;
            continue;
        }
        
        if (readingModel) {
            // Parse model line
            double a, b, aErr, bErr, chi2ndf;
            std::istringstream iss(line);
            iss >> a >> b >> aErr >> bErr >> chi2ndf;
            model_.SetParams(a, b, aErr, bErr);
        } else {
            // Parse point line
            CalibrationPoint p;
            std::istringstream iss(line);
            int enabled;
            iss >> p.energy >> p.energyErr >> p.channel >> p.channelErr 
                >> p.sourceHistogram >> p.isotope >> enabled;
            p.enabled = (enabled != 0);
            AddPoint(p);
        }
    }
    
    infile.close();
    return true;
}
