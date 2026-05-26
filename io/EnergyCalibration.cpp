#include "EnergyCalibration.h"
#include "TH1.h"
#include "TH1F.h"
#include <algorithm>
#include <fstream>
#include <sstream>
#include <iostream>

EnergyCalibration::~EnergyCalibration() {
    // TF1 objects will be deleted by ROOT when the histogram is deleted
    // Just clear the map
    calibrations_.clear();
}

void EnergyCalibration::AddCalibration(const std::string& name,
                                      const LinearCalibrationModel& model) {
    calibrations_[name] = model;
}

const LinearCalibrationModel* EnergyCalibration::GetCalibration(
    const std::string& name) const {
    auto it = calibrations_.find(name);
    if (it != calibrations_.end()) {
        return &it->second;
    }
    return nullptr;
}

LinearCalibrationModel* EnergyCalibration::GetCalibration(const std::string& name) {
    auto it = calibrations_.find(name);
    if (it != calibrations_.end()) {
        return &it->second;
    }
    return nullptr;
}

std::vector<std::string> EnergyCalibration::GetCalibrationNames() const {
    std::vector<std::string> names;
    for (const auto& p : calibrations_) {
        names.push_back(p.first);
    }
    return names;
}

TH1* EnergyCalibration::ApplyCalibrationToHistogram(
    TH1* hist, const std::string& calibName, const std::string& newHistName) {
    
    if (!hist) return nullptr;
    
    const auto* calib = GetCalibration(calibName);
    if (!calib || !calib->IsValid()) return nullptr;
    
    // Create energy bin edges from channel bin edges
    auto energyBinEdges = ComputeEnergyBinEdges(hist, *calib);
    if (energyBinEdges.empty()) return nullptr;
    
    // Create new histogram with energy bins
    std::string outName = newHistName.empty() ? 
        std::string(hist->GetName()) + "_calib" : newHistName;
    
    TH1* newHist = new TH1F(outName.c_str(), hist->GetTitle(),
                           energyBinEdges.size() - 1,
                           &energyBinEdges[0]);
    newHist->GetXaxis()->SetTitle("Energy (keV)");
    
    // Copy bin contents and errors
    int nBins = hist->GetNbinsX();
    for (int i = 1; i <= nBins; i++) {
        newHist->SetBinContent(i, hist->GetBinContent(i));
        newHist->SetBinError(i, hist->GetBinError(i));
    }
    
    return newHist;
}

bool EnergyCalibration::ApplyCalibrationInPlace(TH1*& hist,
                                               const std::string& calibName) {
    if (!hist) return false;
    
    TH1* calibHist = ApplyCalibrationToHistogram(hist, calibName);
    if (!calibHist) return false;
    
    // Replace original
    delete hist;
    hist = calibHist;
    return true;
}

std::vector<double> EnergyCalibration::ChannelsToEnergies(
    const std::vector<double>& channels,
    const std::string& calibName) const {
    
    std::vector<double> energies;
    const auto* calib = GetCalibration(calibName);
    if (!calib || !calib->IsValid()) return energies;
    
    for (double ch : channels) {
        energies.push_back(calib->ChannelToEnergy(ch));
    }
    return energies;
}

bool EnergyCalibration::LoadFromFile(const std::string& dirPath) {
    // Load all calibration files from directory
    // Format: calib_<name>.txt containing A, B, Aerr, Berr, chi2ndf
    
    // For now, implement a simple file format
    std::ifstream infile(dirPath + "/calibrations.txt");
    if (!infile.is_open()) return false;
    
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;
        
        std::istringstream iss(line);
        std::string name;
        double a, b, aErr, bErr, chi2ndf;
        
        iss >> name >> a >> b >> aErr >> bErr >> chi2ndf;
        
        LinearCalibrationModel model;
        model.SetParams(a, b, aErr, bErr);
        AddCalibration(name, model);
    }
    
    infile.close();
    return true;
}

bool EnergyCalibration::SaveToFile(const std::string& dirPath) const {
    // Save all calibrations to directory
    std::ofstream outfile(dirPath + "/calibrations.txt");
    if (!outfile.is_open()) return false;
    
    outfile << "# name A B AErr BErr Chi2NDF\n";
    for (const auto& p : calibrations_) {
        outfile << p.first << " ";
        outfile << p.second.GetParamA() << " ";
        outfile << p.second.GetParamB() << " ";
        outfile << p.second.GetParamAErr() << " ";
        outfile << p.second.GetParamBErr() << " ";
        outfile << p.second.GetChi2NDF() << "\n";
    }
    
    outfile.close();
    return true;
}

std::vector<double> EnergyCalibration::ComputeEnergyBinEdges(
    const TH1* hist, const LinearCalibrationModel& calib) const {
    
    std::vector<double> energyEdges;
    int nBins = hist->GetNbinsX();
    
    // Get the bin edges from the original histogram
    for (int i = 0; i <= nBins; i++) {
        double channelEdge = hist->GetBinLowEdge(i);
        if (i == nBins) {
            // For last edge, use upper edge of last bin
            channelEdge = hist->GetBinLowEdge(nBins) + hist->GetBinWidth(nBins);
        }
        
        double energyEdge = calib.ChannelToEnergy(channelEdge);
        energyEdges.push_back(energyEdge);
    }
    
    // Sort in ascending order
    std::sort(energyEdges.begin(), energyEdges.end());
    
    return energyEdges;
}
