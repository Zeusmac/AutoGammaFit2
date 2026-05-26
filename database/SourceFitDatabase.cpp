#include "SourceFitDatabase.h"
#include <sstream>
#include <algorithm>

bool SourceFitDatabase::Load(const std::string& cachePath) {
    std::ifstream infile(cachePath);
    if (!infile.is_open()) {
        Debug::Log(Debug::CACHE, std::string("SourceFitDatabase: No cache found at ") + cachePath);
        return false;
    }

    entries_.clear();
    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty() || line[0] == '#') continue;

        SourcePeakEntry entry;
        if (ParseLine(line, entry)) {
            entries_.push_back(entry);
        }
    }

    infile.close();
    Debug::Log(Debug::CACHE, std::string("SourceFitDatabase: Loaded ") +
               std::to_string(entries_.size()) + " source peak entries");
    return true;
}

bool SourceFitDatabase::Save(const std::string& cachePath) const {
    std::ofstream outfile(cachePath);
    if (!outfile.is_open()) {
        Debug::Log(Debug::CACHE, std::string("SourceFitDatabase: Cannot write to ") + cachePath);
        return false;
    }
    
    outfile << "# Source Fit Database\n";
    outfile << "# key | energy | channel | channelErr | sigma | sigmaErr | area | areaErr ";
    outfile << "| chi2ndf | isotope | sourceHistogram | fitDate | fitMethod\n";
    
    for (const auto& entry : entries_) {
        outfile << SerializeEntry(entry) << "\n";
    }
    
    outfile.close();
    Debug::Log(Debug::CACHE, std::string("SourceFitDatabase: Saved ") +
               std::to_string(entries_.size()) + " source peak entries");
    return true;
}

void SourceFitDatabase::AddEntry(const SourcePeakEntry& entry) {
    entries_.push_back(entry);
}

const SourcePeakEntry* SourceFitDatabase::GetEntry(const std::string& sourceHist, 
                                                   double energy, double toleranceKeV) const {
    for (const auto& entry : entries_) {
        if (entry.sourceHistogram == sourceHist && 
            std::abs(entry.energy - energy) < toleranceKeV) {
            return &entry;
        }
    }
    return nullptr;
}

std::vector<const SourcePeakEntry*> SourceFitDatabase::GetEntriesForSource(
    const std::string& sourceHist) const {
    std::vector<const SourcePeakEntry*> result;
    for (const auto& entry : entries_) {
        if (entry.sourceHistogram == sourceHist) {
            result.push_back(&entry);
        }
    }
    return result;
}

std::vector<const SourcePeakEntry*> SourceFitDatabase::GetEntriesForIsotope(
    const std::string& isotope) const {
    std::vector<const SourcePeakEntry*> result;
    for (const auto& entry : entries_) {
        if (entry.isotope == isotope) {
            result.push_back(&entry);
        }
    }
    return result;
}

std::vector<std::string> SourceFitDatabase::GetUniqueSourceHistograms() const {
    std::vector<std::string> result;
    for (const auto& entry : entries_) {
        if (std::find(result.begin(), result.end(), entry.sourceHistogram) == result.end()) {
            result.push_back(entry.sourceHistogram);
        }
    }
    return result;
}

std::vector<std::string> SourceFitDatabase::GetUniqueIsotopes() const {
    std::vector<std::string> result;
    for (const auto& entry : entries_) {
        if (!entry.isotope.empty() && 
            std::find(result.begin(), result.end(), entry.isotope) == result.end()) {
            result.push_back(entry.isotope);
        }
    }
    return result;
}

bool SourceFitDatabase::ParseLine(const std::string& line, SourcePeakEntry& entry) {
    std::istringstream iss(line);
    std::string token;
    
    // Parse: key | energy | channel | channelErr | sigma | sigmaErr | area | areaErr 
    //        | chi2ndf | isotope | sourceHistogram | fitDate | fitMethod
    try {
        std::getline(iss, entry.key, '|');
        iss >> entry.energy;
        iss.ignore();
        iss >> entry.centroidChannel;
        iss.ignore();
        iss >> entry.centroidChannelErr;
        iss.ignore();
        iss >> entry.sigma;
        iss.ignore();
        iss >> entry.sigmaErr;
        iss.ignore();
        iss >> entry.area;
        iss.ignore();
        iss >> entry.areaErr;
        iss.ignore();
        iss >> entry.chi2ndf;
        iss.ignore();
        std::getline(iss, entry.isotope, '|');
        std::getline(iss, entry.sourceHistogram, '|');
        std::getline(iss, entry.fitDate, '|');
        std::getline(iss, entry.fitMethod, '|');
        return true;
    } catch (...) {
        return false;
    }
}

std::string SourceFitDatabase::SerializeEntry(const SourcePeakEntry& entry) const {
    std::ostringstream oss;
    oss << entry.key << " | ";
    oss << entry.energy << " | ";
    oss << entry.centroidChannel << " | ";
    oss << entry.centroidChannelErr << " | ";
    oss << entry.sigma << " | ";
    oss << entry.sigmaErr << " | ";
    oss << entry.area << " | ";
    oss << entry.areaErr << " | ";
    oss << entry.chi2ndf << " | ";
    oss << entry.isotope << " | ";
    oss << entry.sourceHistogram << " | ";
    oss << entry.fitDate << " | ";
    oss << entry.fitMethod;
    return oss.str();
}
