#ifndef LEVELSCHEME_H
#define LEVELSCHEME_H

#include "NuclearData.h"
#include <string>
#include <vector>
#include <cmath>

// ─────────────────────────────────────────────────────────────────────────────
// Enumerations
// ─────────────────────────────────────────────────────────────────────────────

enum class LSBetaType {
    Unknown,
    Fermi,           // ΔJ=0, Δπ=0 (allowed Fermi / superallowed)
    GamowTeller,     // ΔJ=0,1 (not 0→0), Δπ=0
    Mixed,           // mixed Fermi + GT
    FirstForbidden,  // Δπ=1: ΔJ=0,1 or unique (ΔJ=2, Δπ=0 "unique 1st")
    Forbidden        // higher-order forbidden
};
std::string     LSBetaTypeStr(LSBetaType t);
LSBetaType      LSBetaTypeFromStr(const std::string& s);

enum class LSMultipol {
    Unknown,
    E1, E2, E3, E4,
    M1, M2, M3,
    M1E2,   // mixed M1+E2
    E1M2    // mixed E1+M2
};
std::string     LSMultipolStr(LSMultipol m);
LSMultipol      LSMultipolFromStr(const std::string& s);

// ─────────────────────────────────────────────────────────────────────────────
// LSLevel — one nuclear excited state (id=0 is ground state convention)
// ─────────────────────────────────────────────────────────────────────────────
struct LSLevel {
    int    id          = 0;
    double energy      = 0.0;      // excitation energy (keV)
    std::string jpi;               // "0+", "2+", "3/2-", etc.
    double halflife_s  = -1.0;     // isomeric HL; -1 = prompt

    // β-feeding from parent decay TO this level
    double betaFeedPct = 0.0;      // % of all β decays feeding this level
    double betaFeedErr = 0.0;
    LSBetaType betaType = LSBetaType::Unknown;
    double betaEndptKeV = 0.0;     // β⁻ endpoint to this level (keV)

    // Derived quantities (filled by ComputeAllLogFt)
    double logFt    = -1.0;        // -1 = not computed
    double B_GT     = -1.0;
    double B_Fermi  = -1.0;

    // Filled by CascadeBalance
    double population = 0.0;       // total in-flow: β-feed + cascade-in
};

// ─────────────────────────────────────────────────────────────────────────────
// LSTransition — one gamma transition between two levels
// ─────────────────────────────────────────────────────────────────────────────
struct LSTransition {
    int    id           = 0;
    int    from_id      = -1;     // higher-energy LSLevel id
    int    to_id        = -1;     // lower-energy  LSLevel id (0 = ground state)
    double energy       = 0.0;   // measured gamma energy (keV)
    double intensity    = 0.0;   // relative or area-based intensity
    double intensityErr = 0.0;
    LSMultipol  multipol = LSMultipol::Unknown;
    std::string multipolStr;      // user-readable: "E2", "M1+E2", etc.
    double mixingRatio  = 0.0;   // δ for mixed transitions

    // Link to a fitted peak in a cache file
    std::string linkedCache;
    double      linkedEnergy = 0.0;  // 0 = not linked
};

// ─────────────────────────────────────────────────────────────────────────────
// LevelSchemeData — full user-editable level scheme for one isotope
// ─────────────────────────────────────────────────────────────────────────────
struct LevelSchemeData {
    std::string isoID;         // daughter nucleus: e.g. "44S"
    std::string parentID;      // β⁻ parent: e.g. "44P"
    std::string parentJpi;     // parent ground-state Jπ
    double parentHL_s   = -1.0;
    double parentQbeta  = 0.0; // keV
    int    parentZ      = 0;   // proton number of parent (daughter Z = parentZ+1)

    std::vector<LSLevel>      levels;
    std::vector<LSTransition> transitions;

    int _nextLevelId = 1;  // 0 reserved for ground state
    int _nextTransId = 0;

    // ── Mutators ──────────────────────────────────────────────────────────────
    int  AddLevel(double energy, const std::string& jpi);  // returns id
    int  AddTransition(int from_id, int to_id, double energy);
    void RemoveLevel(int id);
    void RemoveTransition(int id);
    void SortLevels();      // sort levels by energy ascending
    void Clear();

    // ── Accessors ─────────────────────────────────────────────────────────────
    LSLevel*       FindLevel(int id);
    const LSLevel* FindLevel(int id) const;
    LSTransition*  FindTransition(int id);
    const LSLevel* FindLevelByEnergy(double e, double tol = 0.5) const;
    std::vector<const LSTransition*> TransitionsFrom(int level_id) const;
    std::vector<const LSTransition*> TransitionsTo(int level_id)   const;

    bool empty() const { return levels.empty(); }

    // ── Branching Ratios ──────────────────────────────────────────────────────
    // For a given level: fraction of intensity going to each lower level.
    struct BranchResult {
        int    to_id;
        int    trans_id;
        double absIntensity;   // measured intensity for this branch
        double branchingRatio; // fraction (0-1) of total outflow
        double relativeBR;     // ratio to strongest branch (strongest = 1.0)
    };
    std::vector<BranchResult> BranchingRatios(int from_id) const;

    // ── Cascade Balance (infrastructure for complex solver) ───────────────────
    // For each level: check that β-feed + incoming-gamma ≈ outgoing-gamma.
    // This is the intensity-balance equation; residual > 0 indicates missing
    // feedings or incorrectly linked peaks.
    struct BalanceRow {
        int    level_id;
        double betaIn;     // β-feeding intensity entering this level
        double gammaIn;    // sum of intensities of γ transitions feeding this level
        double gammaOut;   // sum of intensities of γ transitions leaving this level
        double residual;   // |betaIn + gammaIn - gammaOut| / max(betaIn+gammaIn, 1e-10)
    };
    std::vector<BalanceRow> CascadeBalance() const;

    // ── Log ft & Transition Strength ─────────────────────────────────────────
    // Phase-space integral approximation (Feenberg-Trigg, ±0.3 in log ft):
    //   W0 = (Q_β - E_level) / 511.0 + 1.0   [endpoint in m_e c² units]
    //   log10(f) ≈ 4.0·log10(W0-1) + C(Z)
    //   C(Z) ≈ 0.78 + 0.012·Z^(2/3)   [empirical Coulomb correction]
    //   log ft = log10(f) + log10(T½_parent_s)
    // Accurate to ~0.3 for non-unique allowed transitions, Z < 50.
    static double ApproxLogF(double Q_beta_keV, double E_level_keV, int Z_daughter);

    // Fill logFt, B_GT, B_Fermi for every level with betaFeedPct > 0.
    // Requires parentHL_s, parentQbeta, parentZ to be set.
    void ComputeAllLogFt();

    // B(GT) = K / (ft · g_A²),  K=6146.5 s,  g_A=1.2754 (free nucleon)
    static double LogFtToB_GT(double logft)
    { return 6146.5 / (std::pow(10.0, logft) * 1.2754 * 1.2754); }

    // B(F) = K / ft  (pure Fermi, g_V=1, isospin-zero limit)
    static double LogFtToB_F(double logft)
    { return 6146.5 / std::pow(10.0, logft); }

    // Parse a Jπ string into 2J (integer) and parity (+1 or -1).
    // Returns false if the string cannot be parsed.
    static bool ParseJpi(const std::string& s, int& twoJ, int& parity);

    // Classify a β transition from Jπ strings.
    static LSBetaType ClassifyBeta(const std::string& parentJpi,
                                   const std::string& daughterJpi);

    // ── Seed from NNDC data ───────────────────────────────────────────────────
    // Populate levels and transitions from a fetched NucIsotope.
    // Optionally also provide the parent to seed β-feeding info.
    void SeedFromNNDC(const NucIsotope& iso,
                      const NucIsotope* parent  = nullptr,
                      int               parentZ = 0);

    // ── Persistence ───────────────────────────────────────────────────────────
    bool Save(const std::string& path) const;
    bool Load(const std::string& path);
};

#endif // LEVELSCHEME_H
