#include "GammaDB.h"
#include "Debug.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <map>
#include <string>

using namespace std;

bool GammaDB::Load(const string& filename) {
    ifstream in(filename);
    if (!in.is_open()) {
        cerr << "ERROR: Cannot open DB file: " << filename << "\n";
        return false;
    }

    db.clear();

    string line;
    int skipped = 0;
    while (getline(in, line)) {

        if (line.empty()) continue;
        if (line[0] == '#') continue;

        string iso;
        double energy;
        double intensity = 100.0;

        stringstream ss(line);

        if (!(ss >> iso >> energy)) {
            ++skipped;
            Debug::Log(Debug::GAMMADB, "Skipping bad line: " + line);
            continue;
        }

        bool hasInt = static_cast<bool>(ss >> intensity);

        db.push_back({iso, energy, intensity, hasInt});
        Debug::Log(Debug::GAMMADB,
            "Loaded: " + iso +
            "  E=" + std::to_string(energy) +
            "  I=" + std::to_string(intensity));
    }

    cout << "[GammaDB] Loaded " << db.size() << " lines"
         << (skipped ? "  (" + std::to_string(skipped) + " skipped)" : "")
         << " from " << filename << "\n";
    return true;
}

std::vector<GammaMatch> GammaDB::Match(double energy, double fwhmKeV) const {

    if (fwhmKeV <= 0) fwhmKeV = 5.0;

    double sigma  = fwhmKeV / 2.3548200450309493;
    double window = 3.0 * fwhmKeV;

    std::vector<GammaMatch> matches;

    for (const auto& l : db) {
        double dE = std::fabs(l.energy - energy);
        if (dE > window) continue;

        double score = l.intensity * std::exp(-0.5 * (dE / sigma) * (dE / sigma));
        matches.push_back({l.isotope, l.energy, dE, l.intensity, score});
    }

    std::sort(matches.begin(), matches.end(),
              [](const GammaMatch& a, const GammaMatch& b) {
                  return a.score > b.score;
              });

    Debug::Log(Debug::GAMMADB,
        "Match E=" + std::to_string(energy) +
        "  FWHM=" + std::to_string(fwhmKeV) +
        "  window=+-" + std::to_string(window) +
        "  found=" + std::to_string(matches.size()) +
        (matches.empty() ? "" : "  best=" + matches[0].isotope));

    return matches;
}

std::vector<IsotopeCandidate> GammaDB::IdentifyIsotopes(
    const std::vector<double>& observedEnergies,
    double fwhmKeV) const
{
    if (fwhmKeV <= 0) fwhmKeV = 5.0;

    double sigma  = fwhmKeV / 2.3548200450309493;
    double window = 3.0 * fwhmKeV;

    struct Accum {
        int    nMatched = 0;
        double totalScore = 0.0;
        std::vector<double> matchedEnergies;
    };
    std::map<std::string, Accum> acc;

    for (const auto& l : db) {
        for (double obs : observedEnergies) {
            double dE = std::fabs(l.energy - obs);
            if (dE > window) continue;

            double score = l.intensity * std::exp(-0.5 * (dE / sigma) * (dE / sigma));
            auto& a = acc[l.isotope];
            a.nMatched++;
            a.totalScore += score;
            a.matchedEnergies.push_back(obs);
        }
    }

    std::vector<IsotopeCandidate> candidates;
    candidates.reserve(acc.size());

    for (auto& kv : acc) {
        IsotopeCandidate c;
        c.isotope         = kv.first;
        c.nMatched        = kv.second.nMatched;
        c.totalScore      = kv.second.totalScore;
        c.matchedEnergies = kv.second.matchedEnergies;
        candidates.push_back(std::move(c));
    }

    std::sort(candidates.begin(), candidates.end(),
              [](const IsotopeCandidate& a, const IsotopeCandidate& b) {
                  if (a.nMatched != b.nMatched) return a.nMatched > b.nMatched;
                  return a.totalScore > b.totalScore;
              });

    Debug::Log(Debug::GAMMADB,
        "IdentifyIsotopes: " + std::to_string(observedEnergies.size()) +
        " observed peaks  -> " + std::to_string(candidates.size()) +
        " candidates" +
        (candidates.empty() ? "" :
            "  #1=" + candidates[0].isotope +
            " lines=" + std::to_string(candidates[0].nMatched)));

    return candidates;
}
