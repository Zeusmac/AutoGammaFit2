#ifndef NUCLEARDATA_H
#define NUCLEARDATA_H

#include <string>
#include <vector>
#include <map>
#include <cmath>

// ─────────────────────────────────────────────────────────────────────────────
// Nuclear level / gamma / beta structures from IAEA LiveChart
// ─────────────────────────────────────────────────────────────────────────────

struct NucLevel {
    double energy   = 0.0;   // keV
    std::string jpi;         // spin-parity
    double halflife_s = 0.0; // seconds (0=stable or unknown)
    std::string hl_str;      // human-readable half-life
};

struct NucGamma {
    double start_level = 0.0; // level energy (keV)
    double energy      = 0.0; // gamma energy (keV)
    double intensity   = 0.0; // relative intensity
    std::string multipolarity;
};

struct NucBetaBranch {
    double endpoint   = 0.0;  // keV
    double intensity  = 0.0;  // %
    std::string daughter_jpi;
};

struct NucIsotope {
    int    A = 0;
    int    Z = 0;
    int    N = 0;
    std::string symbol;   // e.g. "Kr"
    std::string jpi;      // ground state spin-parity
    std::string hl_str;   // human-readable half-life
    double halflife_s = 0.0;  // seconds; 0 = stable or unknown
    double qBetaMinus = 0.0;  // keV
    std::string decayMode1;   // primary decay mode: "B-", "B+", "EC", "A", "B-N", etc.
    double decay1frac = 0.0;  // fraction (0-100%) of primary decay mode

    std::vector<NucLevel>      levels;
    std::vector<NucGamma>      gammas;
    std::vector<NucBetaBranch> betaBranches;

    std::string fetchDate;
    bool _valid = false;

    bool valid() const { return _valid; }
};

// NUBASE entry with beta-delayed neutron emission info
struct NubaseEntry {
    int A = 0;
    int Z = 0;
    double mass_excess = 0.0;  // keV
    std::string jpi;
    double halflife_s = 0.0;
    std::string hl_str;
    double brBetaN  = 0.0;  // branching ratio for beta-n  (%)
    double brBeta2N = 0.0;  // branching ratio for beta-2n (%)
};

// AME table entry
struct AMEEntry {
    int    A = 0;
    int    Z = 0;
    double mass_excess = 0.0;   // keV/c^2
    double binding_per_A = 0.0; // keV/nucleon
    double qBetaMinus = 0.0;    // keV
    double qAlpha = 0.0;        // keV
};

// ─────────────────────────────────────────────────────────────────────────────
// Element symbol <-> Z lookup tables
// ─────────────────────────────────────────────────────────────────────────────

inline const char* kSymbolTable[] = {
    "n",  "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",
    "Ne", "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",
    "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu",
    "Zn", "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",
    "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr",
    "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm",
    "Yb", "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au",
    "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac",
    "Th", "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",
    "Fm", "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt",
    "Ds", "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
};
static constexpr int kSymbolTableSize = 119;

inline int NucSymbolToZ(const std::string& sym) {
    if (sym.empty()) return -1;
    for (int z = 0; z < kSymbolTableSize; z++) {
        if (sym == kSymbolTable[z]) return z;
    }
    return -1;
}

inline std::string NucZToSymbol(int Z) {
    if (Z < 0 || Z >= kSymbolTableSize) return "?";
    return kSymbolTable[Z];
}

#endif // NUCLEARDATA_H
