#include "LevelScheme.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <cstdio>

// ─────────────────────────────────────────────────────────────────────────────
// Enum helpers
// ─────────────────────────────────────────────────────────────────────────────

std::string LSBetaTypeStr(LSBetaType t) {
    switch (t) {
        case LSBetaType::Fermi:          return "Fermi";
        case LSBetaType::GamowTeller:    return "GT";
        case LSBetaType::Mixed:          return "Mixed";
        case LSBetaType::FirstForbidden: return "1F";
        case LSBetaType::Forbidden:      return "Forbidden";
        default:                         return "Unknown";
    }
}

LSBetaType LSBetaTypeFromStr(const std::string& s) {
    if (s == "Fermi")    return LSBetaType::Fermi;
    if (s == "GT")       return LSBetaType::GamowTeller;
    if (s == "Mixed")    return LSBetaType::Mixed;
    if (s == "1F")       return LSBetaType::FirstForbidden;
    if (s == "Forbidden")return LSBetaType::Forbidden;
    return LSBetaType::Unknown;
}

std::string LSMultipolStr(LSMultipol m) {
    switch (m) {
        case LSMultipol::E1:   return "E1";
        case LSMultipol::E2:   return "E2";
        case LSMultipol::E3:   return "E3";
        case LSMultipol::E4:   return "E4";
        case LSMultipol::M1:   return "M1";
        case LSMultipol::M2:   return "M2";
        case LSMultipol::M3:   return "M3";
        case LSMultipol::M1E2: return "M1+E2";
        case LSMultipol::E1M2: return "E1+M2";
        default:               return "Unknown";
    }
}

LSMultipol LSMultipolFromStr(const std::string& s) {
    if (s == "E1")     return LSMultipol::E1;
    if (s == "E2")     return LSMultipol::E2;
    if (s == "E3")     return LSMultipol::E3;
    if (s == "E4")     return LSMultipol::E4;
    if (s == "M1")     return LSMultipol::M1;
    if (s == "M2")     return LSMultipol::M2;
    if (s == "M3")     return LSMultipol::M3;
    if (s == "M1+E2" || s == "M1E2") return LSMultipol::M1E2;
    if (s == "E1+M2" || s == "E1M2") return LSMultipol::E1M2;
    return LSMultipol::Unknown;
}

// ─────────────────────────────────────────────────────────────────────────────
// Mutators
// ─────────────────────────────────────────────────────────────────────────────

int LevelSchemeData::AddLevel(double energy, const std::string& jpi)
{
    int id = (energy == 0.0 && levels.empty()) ? 0 : _nextLevelId++;
    if (id == 0) {
        // Ground state: check if already present
        for (auto& lv : levels) { if (lv.id == 0) return 0; }
    }
    LSLevel lv;
    lv.id     = id;
    lv.energy = energy;
    lv.jpi    = jpi;
    levels.push_back(lv);
    SortLevels();
    return id;
}

int LevelSchemeData::AddTransition(int from_id, int to_id, double energy)
{
    LSTransition tr;
    tr.id      = _nextTransId++;
    tr.from_id = from_id;
    tr.to_id   = to_id;
    tr.energy  = energy;
    transitions.push_back(tr);
    return tr.id;
}

void LevelSchemeData::RemoveLevel(int id)
{
    levels.erase(std::remove_if(levels.begin(), levels.end(),
        [id](const LSLevel& lv){ return lv.id == id; }), levels.end());
    // Also remove orphaned transitions
    transitions.erase(std::remove_if(transitions.begin(), transitions.end(),
        [id](const LSTransition& tr){
            return tr.from_id == id || tr.to_id == id; }),
        transitions.end());
}

void LevelSchemeData::RemoveTransition(int id)
{
    transitions.erase(std::remove_if(transitions.begin(), transitions.end(),
        [id](const LSTransition& tr){ return tr.id == id; }), transitions.end());
}

void LevelSchemeData::SortLevels()
{
    std::sort(levels.begin(), levels.end(),
        [](const LSLevel& a, const LSLevel& b){ return a.energy < b.energy; });
}

void LevelSchemeData::Clear()
{
    levels.clear();
    transitions.clear();
    _nextLevelId = 1;
    _nextTransId = 0;
}

// ─────────────────────────────────────────────────────────────────────────────
// Accessors
// ─────────────────────────────────────────────────────────────────────────────

LSLevel* LevelSchemeData::FindLevel(int id) {
    for (auto& lv : levels) if (lv.id == id) return &lv;
    return nullptr;
}
const LSLevel* LevelSchemeData::FindLevel(int id) const {
    for (const auto& lv : levels) if (lv.id == id) return &lv;
    return nullptr;
}
LSTransition* LevelSchemeData::FindTransition(int id) {
    for (auto& tr : transitions) if (tr.id == id) return &tr;
    return nullptr;
}
const LSLevel* LevelSchemeData::FindLevelByEnergy(double e, double tol) const {
    const LSLevel* best = nullptr;
    double bestD = tol;
    for (const auto& lv : levels) {
        double d = std::abs(lv.energy - e);
        if (d < bestD) { bestD = d; best = &lv; }
    }
    return best;
}
std::vector<const LSTransition*> LevelSchemeData::TransitionsFrom(int level_id) const {
    std::vector<const LSTransition*> out;
    for (const auto& tr : transitions)
        if (tr.from_id == level_id) out.push_back(&tr);
    return out;
}
std::vector<const LSTransition*> LevelSchemeData::TransitionsTo(int level_id) const {
    std::vector<const LSTransition*> out;
    for (const auto& tr : transitions)
        if (tr.to_id == level_id) out.push_back(&tr);
    return out;
}

// ─────────────────────────────────────────────────────────────────────────────
// Branching Ratios
// ─────────────────────────────────────────────────────────────────────────────

std::vector<LevelSchemeData::BranchResult>
LevelSchemeData::BranchingRatios(int from_id) const
{
    std::vector<BranchResult> out;
    auto trs = TransitionsFrom(from_id);
    if (trs.empty()) return out;

    double total = 0.0;
    for (const auto* tr : trs) total += tr->intensity;
    if (total <= 0.0) {
        // Cannot compute ratios without intensities — return structure only
        for (const auto* tr : trs)
            out.push_back({tr->to_id, tr->id, 0.0, 0.0, 0.0});
        return out;
    }

    double strongest = 0.0;
    for (const auto* tr : trs)
        if (tr->intensity > strongest) strongest = tr->intensity;

    for (const auto* tr : trs) {
        BranchResult r;
        r.to_id          = tr->to_id;
        r.trans_id       = tr->id;
        r.absIntensity   = tr->intensity;
        r.branchingRatio = tr->intensity / total;
        r.relativeBR     = (strongest > 0.0) ? tr->intensity / strongest : 0.0;
        out.push_back(r);
    }
    // Sort by intensity descending
    std::sort(out.begin(), out.end(),
        [](const BranchResult& a, const BranchResult& b){
            return a.absIntensity > b.absIntensity; });
    return out;
}

// ─────────────────────────────────────────────────────────────────────────────
// Cascade Balance
// ─────────────────────────────────────────────────────────────────────────────

std::vector<LevelSchemeData::BalanceRow> LevelSchemeData::CascadeBalance() const
{
    std::vector<BalanceRow> out;
    for (const auto& lv : levels) {
        BalanceRow row;
        row.level_id = lv.id;
        row.betaIn   = lv.betaFeedPct;

        row.gammaIn = 0.0;
        for (const auto& tr : transitions)
            if (tr.to_id == lv.id) row.gammaIn += tr.intensity;

        row.gammaOut = 0.0;
        for (const auto& tr : transitions)
            if (tr.from_id == lv.id) row.gammaOut += tr.intensity;

        double inflow = row.betaIn + row.gammaIn;
        row.residual  = (inflow > 1e-10)
            ? std::abs(inflow - row.gammaOut) / inflow
            : (row.gammaOut > 0 ? 1.0 : 0.0);
        out.push_back(row);
    }
    return out;
}

// ─────────────────────────────────────────────────────────────────────────────
// Log ft / Transition Strength
// ─────────────────────────────────────────────────────────────────────────────

double LevelSchemeData::ApproxLogF(double Q_beta_keV, double E_level_keV, int Z_daughter)
{
    static constexpr double me_keV = 511.0;
    double E0 = Q_beta_keV - E_level_keV;   // endpoint kinetic energy (keV)
    if (E0 <= 0.0) return -99.0;            // energetically forbidden
    double W0 = E0 / me_keV + 1.0;          // endpoint in units of m_e c²
    if (W0 <= 1.0) return -99.0;

    // Feenberg-Trigg approximation with Coulomb correction
    // C(Z) ≈ 0.78 + 0.012·Z^(2/3) — empirical, valid for Z < 60
    double Cz   = 0.78 + 0.012 * std::pow(static_cast<double>(Z_daughter), 2.0/3.0);
    double logf  = 4.0 * std::log10(W0 - 1.0) + Cz;
    return logf;
}

void LevelSchemeData::ComputeAllLogFt()
{
    if (parentHL_s <= 0.0 || parentQbeta <= 0.0) return;
    int Zd = parentZ + 1;  // daughter Z for β⁻

    double logT = std::log10(parentHL_s);

    for (auto& lv : levels) {
        if (lv.betaFeedPct <= 0.0) { lv.logFt = -1.0; lv.B_GT = -1.0; lv.B_Fermi = -1.0; continue; }

        double logf = ApproxLogF(parentQbeta, lv.energy, Zd);
        if (logf <= -90.0) { lv.logFt = -1.0; lv.B_GT = -1.0; lv.B_Fermi = -1.0; continue; }

        // Partial half-life: T½_partial = T½ / (betaFeedPct/100)
        double feedFrac = lv.betaFeedPct / 100.0;
        if (feedFrac <= 0.0) { lv.logFt = -1.0; continue; }
        double T_partial  = parentHL_s / feedFrac;
        double logft      = logf + std::log10(T_partial);
        lv.logFt   = logft;
        lv.B_GT    = LogFtToB_GT(logft);
        lv.B_Fermi = LogFtToB_F(logft);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Jπ parsing and β-transition classification
// ─────────────────────────────────────────────────────────────────────────────

bool LevelSchemeData::ParseJpi(const std::string& s, int& twoJ, int& parity)
{
    if (s.empty()) return false;
    twoJ   = 0;
    parity = 0;

    std::string spin = s;
    if (spin.back() == '+') { parity = +1; spin.pop_back(); }
    else if (spin.back() == '-') { parity = -1; spin.pop_back(); }
    else return false;  // no parity given

    // spin string is now e.g. "0", "2", "3/2", "5/2"
    auto slash = spin.find('/');
    if (slash == std::string::npos) {
        // Integer spin
        try { twoJ = std::stoi(spin) * 2; return true; }
        catch (...) { return false; }
    } else {
        // Half-integer spin: numerator / 2
        try {
            int num = std::stoi(spin.substr(0, slash));
            int den = std::stoi(spin.substr(slash + 1));
            if (den != 2) return false;
            twoJ = num;
            return true;
        } catch (...) { return false; }
    }
}

LSBetaType LevelSchemeData::ClassifyBeta(const std::string& parentJpi,
                                          const std::string& daughterJpi)
{
    int twoJp, pip, twoJd, pid;
    if (!ParseJpi(parentJpi, twoJp, pip) || !ParseJpi(daughterJpi, twoJd, pid))
        return LSBetaType::Unknown;

    int deltaJ  = std::abs(twoJp - twoJd);  // in units of ½ (so 2 = ΔJ=1)
    int deltaPi = pip * pid;                  // +1 = same parity, -1 = opposite

    // Unique first forbidden: ΔJ = 2, Δπ = -1
    if (deltaPi == -1 && deltaJ == 2) return LSBetaType::FirstForbidden;
    // Regular first forbidden: ΔJ = 0 or 2, Δπ = -1
    if (deltaPi == -1 && (deltaJ == 0 || deltaJ == 2)) return LSBetaType::FirstForbidden;
    // Fermi: ΔJ=0, Δπ=+1
    if (deltaPi == +1 && deltaJ == 0) {
        // 0→0 is pure Fermi (if not 0→0 it could be mixed)
        return (twoJp == 0 && twoJd == 0) ? LSBetaType::Fermi : LSBetaType::Mixed;
    }
    // GT: ΔJ=0 or 2, Δπ=+1 (but not 0→0)
    if (deltaPi == +1 && (deltaJ == 0 || deltaJ == 2))
        return LSBetaType::GamowTeller;

    return LSBetaType::Forbidden;
}

// ─────────────────────────────────────────────────────────────────────────────
// Seed from NNDC NucIsotope data
// ─────────────────────────────────────────────────────────────────────────────

void LevelSchemeData::SeedFromNNDC(const NucIsotope& iso,
                                    const NucIsotope* parent,
                                    int pZ)
{
    Clear();
    isoID  = std::to_string(iso.A) + iso.symbol;

    if (parent) {
        parentID    = std::to_string(parent->A) + parent->symbol;
        parentJpi   = parent->jpi;
        parentHL_s  = parent->halflife_s;
        parentQbeta = parent->qBetaMinus;
        parentZ     = pZ;
    }

    // Ground state (always id=0)
    {
        LSLevel gs;
        gs.id     = 0;
        gs.energy = 0.0;
        gs.jpi    = iso.jpi;
        levels.push_back(gs);
    }

    // Excited levels
    for (const auto& nl : iso.levels) {
        if (nl.energy <= 0.0) continue;
        if (FindLevelByEnergy(nl.energy, 0.1)) continue;  // skip duplicates
        LSLevel lv;
        lv.id        = _nextLevelId++;
        lv.energy    = nl.energy;
        lv.jpi       = nl.jpi;
        lv.halflife_s = nl.halflife_s;
        levels.push_back(lv);
    }
    SortLevels();

    // Gamma transitions
    for (const auto& ng : iso.gammas) {
        if (ng.energy <= 0.0) continue;
        const LSLevel* fromLv = FindLevelByEnergy(ng.start_level, 1.0);
        if (!fromLv) continue;
        double eTo = ng.start_level - ng.energy;
        const LSLevel* toLv = FindLevelByEnergy(std::max(eTo, 0.0), 1.0);
        if (!toLv) toLv = FindLevel(0);  // fallback to ground state
        if (!toLv) continue;

        LSTransition tr;
        tr.id           = _nextTransId++;
        tr.from_id      = fromLv->id;
        tr.to_id        = toLv->id;
        tr.energy       = ng.energy;
        tr.intensity    = ng.intensity;
        tr.intensityErr = 0.0;
        tr.multipol     = LSMultipolFromStr(ng.multipolarity);
        tr.multipolStr  = ng.multipolarity;
        transitions.push_back(tr);
    }

    // Beta branches → seed β-feeding on levels
    if (parent && parent->qBetaMinus > 0.0) {
        for (const auto& bb : parent->betaBranches) {
            if (bb.intensity <= 0.0) continue;
            // Level energy ≈ Q_β - endpoint
            double lvE = parent->qBetaMinus - bb.endpoint;
            const LSLevel* lv = FindLevelByEnergy(std::max(lvE, 0.0), 100.0);
            if (!lv) lv = FindLevel(0);
            if (!lv) continue;
            LSLevel* lvMut = FindLevel(lv->id);
            if (!lvMut) continue;
            lvMut->betaFeedPct  += bb.intensity;
            lvMut->betaEndptKeV  = bb.endpoint;
            lvMut->betaType = ClassifyBeta(parent->jpi, lv->jpi);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Save / Load  (.lsdat — human-readable line-based format)
// ─────────────────────────────────────────────────────────────────────────────

bool LevelSchemeData::Save(const std::string& path) const
{
    std::ofstream out(path);
    if (!out.is_open()) return false;

    out << "# AutoGammaFit LevelSchemeData v1\n";
    out << "ISO "       << isoID       << "\n";
    out << "PARENT "    << parentID    << "\n";
    out << "PARENTJPI " << parentJpi   << "\n";
    out << "PARENTHL "  << parentHL_s  << "\n";
    out << "PARENTQB "  << parentQbeta << "\n";
    out << "PARENTZ "   << parentZ     << "\n";

    // LEVEL  id  energy  jpi  halflife_s  betaFeedPct  betaFeedErr
    //        betaTypeStr  betaEndptKeV  logFt  B_GT  B_Fermi
    for (const auto& lv : levels) {
        out << "LEVEL "
            << lv.id          << " "
            << lv.energy      << " "
            << (lv.jpi.empty() ? "?" : lv.jpi) << " "
            << lv.halflife_s  << " "
            << lv.betaFeedPct << " "
            << lv.betaFeedErr << " "
            << LSBetaTypeStr(lv.betaType) << " "
            << lv.betaEndptKeV << " "
            << lv.logFt       << " "
            << lv.B_GT        << " "
            << lv.B_Fermi     << "\n";
    }

    // TRANS  id  from_id  to_id  energy  intensity  intensityErr
    //        multipolStr  mixingRatio  linkedCache  linkedEnergy
    for (const auto& tr : transitions) {
        std::string lc = tr.linkedCache.empty() ? "(none)" : tr.linkedCache;
        out << "TRANS "
            << tr.id           << " "
            << tr.from_id      << " "
            << tr.to_id        << " "
            << tr.energy       << " "
            << tr.intensity    << " "
            << tr.intensityErr << " "
            << (tr.multipolStr.empty() ? "?" : tr.multipolStr) << " "
            << tr.mixingRatio  << " "
            << lc              << " "
            << tr.linkedEnergy << "\n";
    }
    return true;
}

bool LevelSchemeData::Load(const std::string& path)
{
    std::ifstream in(path);
    if (!in.is_open()) return false;

    Clear();
    int maxLevelId = -1;
    int maxTransId = -1;

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string tag;
        ss >> tag;

        if      (tag == "ISO")       { ss >> isoID; }
        else if (tag == "PARENT")    { ss >> parentID; }
        else if (tag == "PARENTJPI") { ss >> parentJpi; }
        else if (tag == "PARENTHL")  { ss >> parentHL_s; }
        else if (tag == "PARENTQB")  { ss >> parentQbeta; }
        else if (tag == "PARENTZ")   { ss >> parentZ; }
        else if (tag == "LEVEL") {
            LSLevel lv;
            std::string btStr;
            ss >> lv.id >> lv.energy >> lv.jpi >> lv.halflife_s
               >> lv.betaFeedPct >> lv.betaFeedErr >> btStr
               >> lv.betaEndptKeV >> lv.logFt >> lv.B_GT >> lv.B_Fermi;
            lv.betaType = LSBetaTypeFromStr(btStr);
            levels.push_back(lv);
            if (lv.id > maxLevelId) maxLevelId = lv.id;
        }
        else if (tag == "TRANS") {
            LSTransition tr;
            std::string lc;
            ss >> tr.id >> tr.from_id >> tr.to_id >> tr.energy
               >> tr.intensity >> tr.intensityErr >> tr.multipolStr
               >> tr.mixingRatio >> lc >> tr.linkedEnergy;
            if (lc == "(none)") lc.clear();
            tr.linkedCache = lc;
            tr.multipol    = LSMultipolFromStr(tr.multipolStr);
            transitions.push_back(tr);
            if (tr.id > maxTransId) maxTransId = tr.id;
        }
    }

    _nextLevelId = maxLevelId + 1;
    _nextTransId = maxTransId + 1;
    SortLevels();
    return !levels.empty();
}
