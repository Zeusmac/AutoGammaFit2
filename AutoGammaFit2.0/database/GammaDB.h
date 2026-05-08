#ifndef GAMMADB_H
#define GAMMADB_H

#include <vector>
#include <string>
#include <map>

struct GammaLine {
    std::string isotope;
    double energy;
    double intensity = 100.0;
};

// Single-peak match, scored by intensity × Gaussian proximity
struct GammaMatch {
    std::string isotope;
    double energy;
    double deltaE;
    double intensity;
    double score;   // higher is better
};

// Multi-line isotope candidate, scored across all observed peaks
struct IsotopeCandidate {
    std::string isotope;
    int    nMatched;    // lines matched
    double totalScore;  // sum of per-line scores
    std::vector<double> matchedEnergies;  // observed energies that matched
};

class GammaDB {
public:
    bool Load(const std::string& filename);

    // Single-peak match; fwhmKeV sets the Gaussian search window (use res.FWHM(E))
    std::vector<GammaMatch> Match(double energy, double fwhmKeV) const;

    // Multi-line scorer: rank isotopes by how many of their lines appear in
    // the list of observed fitted energies. fwhmKeV is the typical resolution.
    std::vector<IsotopeCandidate> IdentifyIsotopes(
        const std::vector<double>& observedEnergies,
        double fwhmKeV) const;

    std::vector<GammaLine> db;
};

#endif
