#ifndef FITDATABASE_H
#define FITDATABASE_H

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <limits>

#include "TF1.h"
#include "TH1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "Debug.h"

// Per-bin residual metrics computed over the fit window.
struct ResidualMetrics {
    double rms     = std::numeric_limits<double>::max();
    double maxPull = std::numeric_limits<double>::max();
    int    count   = 0;  // bins with err > 0 used in computation
};

struct FitEntry {
    std::string        key;
    std::vector<double> params;
    std::vector<double> paramErrors;   // MINUIT/MINOS errors, same length as params (0 if unknown)
    double chi2ndf     = std::numeric_limits<double>::max();
    double residualRMS = std::numeric_limits<double>::max();
    double maxPull     = std::numeric_limits<double>::max();
    // ROOT fit method flags actually used: "L"=log-likelihood, "M"=IMPROVE,
    // "E"=MINOS, combinations thereof.  Empty string = default chi2/least-squares.
    std::string fitMethod;
    // Fit window [xlo, xhi] in keV — stored so OverlayFitPeaks redraws the
    // exact same range that was used during fitting.  Both are 0 on old entries.
    double xlo = 0.0;
    double xhi = 0.0;
    // Isotope annotation (optional).  label = e.g. "Co-60"; classification =
    // "Parent", "Daughter", "Granddaughter", "Beta-n", "Beta-2n",
    // "Background", or "Custom:<user text>".  Both empty on unannotated entries.
    std::string label;
    std::string classification;

    // Per-Gaussian labels for wide multi-Gaussian entries (peak spread > 4 keV).
    // Empty vector → use .label / .classification for all Gaussians (legacy /
    // tight-cluster entries).  When non-empty, size == number of Gaussians.
    std::vector<std::string> peakLabels;
    std::vector<std::string> peakClassifications;

    // Return the effective label for Gaussian index gi (-1 = whole entry).
    const std::string& PeakLabel(int gi) const {
        if (gi >= 0 && gi < (int)peakLabels.size() && !peakLabels[gi].empty())
            return peakLabels[gi];
        return label;
    }
    const std::string& PeakClass(int gi) const {
        if (gi >= 0 && gi < (int)peakClassifications.size() && !peakClassifications[gi].empty())
            return peakClassifications[gi];
        return classification;
    }
};

class FitDatabase {
public:
    // When false, SeedFromDB always returns false so every fit starts from
    // scratch — but StoreIfBetter still runs, so the cache is still updated.
    bool useCachedSeeds = true;

    // When true, AdaptiveFitter skips MIGRAD entirely and returns the cached
    // TF1 parameters directly.  Groups with no cache entry are skipped.
    // StoreIfBetter is NOT called (cache is read-only in this mode).
    bool cacheOnly = false;

    // Background subtraction settings used when these fits were produced.
    // Stored in the cache file so the GUI can reproduce the matching
    // histogram view automatically.
    bool        bgSubtracted = false;
    int         bgIterations = 14;

    // ROOT file this cache was created from (absolute path).
    // Empty on old cache files that predate this field.
    std::string rootFile;

    static std::string MakeKey(const std::vector<double>& peaks) {
        std::ostringstream oss;
        for (size_t i = 0; i < peaks.size(); i++) {
            if (i) oss << "_";
            oss << std::fixed << std::setprecision(1) << peaks[i];
        }
        return oss.str();
    }

    // Compute pull = (data − fit) / sigma per bin over [xlo, xhi].
    // Returns {rms-of-pulls, max-|pull|}.
    static ResidualMetrics ComputeResiduals(TH1* h, TF1* f,
                                            double xlo, double xhi) {
        double sumSq = 0.0;
        double maxP  = 0.0;
        int    count = 0;
        int binLo = h->FindBin(xlo);
        int binHi = h->FindBin(xhi);
        for (int b = binLo; b <= binHi; b++) {
            double err = h->GetBinError(b);
            if (err <= 0.0) continue;
            double x    = h->GetBinCenter(b);
            double pull = (h->GetBinContent(b) - f->Eval(x)) / err;
            sumSq += pull * pull;
            if (std::abs(pull) > maxP) maxP = std::abs(pull);
            ++count;
        }
        ResidualMetrics m;
        m.rms     = (count > 0) ? std::sqrt(sumSq / count) : 0.0;
        m.maxPull = maxP;
        m.count   = count;
        return m;
    }

    // chi²/ndf = sumSq / (count - nparams), falling back to sumSq/count when
    // there are too few bins.  rms² = sumSq/count, so we un-multiply here.
    static double Chi2Ndf(const ResidualMetrics& rm, int nparams) {
        if (rm.rms >= 1.0e6 || rm.count <= 0)
            return std::numeric_limits<double>::max();
        double sumSq = rm.rms * rm.rms * rm.count;
        int    dof   = rm.count - nparams;
        // When DOF <= 0 (too few bins for the number of parameters — common for
        // closely-spaced double peaks in a narrow window) fall back to chi2/count
        // rather than returning a sentinel.  The value is still meaningful as a
        // pull-RMS squared; callers can check rm.count vs nparams if they need
        // to flag the underdetermined case.
        return (dof > 0) ? sumSq / dof : rm.rms * rm.rms;
    }

    // Like SeedFromDB but always seeds regardless of useCachedSeeds.
    bool ForceSeedFromDB(const std::string& key, TF1* f) {
        auto it = entries_.find(key);
        if (it == entries_.end()) return false;
        const auto& params = it->second.params;
        if (static_cast<int>(params.size()) != f->GetNpar()) return false;
        for (int i = 0; i < f->GetNpar(); i++)
            f->SetParameter(i, params[i]);
        Debug::Log(Debug::DB, "[FitDatabase] Cache-only seed: " + key);
        return true;
    }

    bool SeedFromDB(const std::string& key, TF1* f) {
        if (!useCachedSeeds) {
            Debug::Log(Debug::DB, "[FitDatabase] Seeding disabled — fitting from scratch for " + key);
            return false;
        }
        auto it = entries_.find(key);
        if (it == entries_.end()) {
            Debug::Log(Debug::DB, "[FitDatabase] No cached seed for " + key);
            return false;
        }
        const auto& params = it->second.params;
        if (static_cast<int>(params.size()) != f->GetNpar()) {
            Debug::Log(Debug::DB, "[FitDatabase] Param count mismatch for " + key + " — ignoring cache");
            return false;
        }
        Debug::Log(Debug::DB, "[FitDatabase] Seeding from cache: " + key +
                   "  cached chi2/ndf=" + std::to_string(it->second.chi2ndf) +
                   "  maxPull=" + std::to_string(it->second.maxPull));
        for (int i = 0; i < f->GetNpar(); i++)
            f->SetParameter(i, params[i]);
        return true;
    }

    static FitEntry MakeEntry(const std::string& key,
                               TFitResultPtr& r,
                               TF1* f,
                               int /*nUsed*/,
                               const ResidualMetrics& rm = ResidualMetrics{})
    {
        FitEntry e;
        e.key         = key;
        // Use Pearson chi2/ndf with proper degrees of freedom when available.
        // r->Chi2() from a log-likelihood fit returns -2*NLL, not chi2.
        if (rm.rms < 1.0e6)
            e.chi2ndf = Chi2Ndf(rm, f->GetNpar());
        else if (r.Get() && r->Ndf() > 0)
            e.chi2ndf = r->Chi2() / r->Ndf();
        else
            e.chi2ndf = std::numeric_limits<double>::max();
        e.residualRMS = rm.rms;
        e.maxPull     = rm.maxPull;
        int n = f->GetNpar();
        e.params.resize(n);
        e.paramErrors.resize(n);
        for (int i = 0; i < n; i++) {
            e.params[i]      = f->GetParameter(i);
            e.paramErrors[i] = f->GetParError(i);
        }
        return e;
    }

    // Composite score: chi2/ndf + penalty * maxPull^2.
    // A large single-bin pull (missed peak) penalises an otherwise-OK chi2.
    // Falls back to chi2/ndf alone when residual metrics are unavailable.
    static double CompositeScore(const FitEntry& e) {
        if (e.chi2ndf >= std::numeric_limits<double>::max())
            return std::numeric_limits<double>::max();
        double mp = (e.maxPull < 1.0e6) ? e.maxPull : 0.0;
        return e.chi2ndf + 0.08 * mp * mp;
    }

    // Always store, overriding any existing entry with the same key.
    void ForceStore(const std::string& key, const FitEntry& candidate) {
        Debug::Log(Debug::DB, "[FitDatabase] Force-store (manual override): " + key);
        entries_[key] = candidate;
    }

    // Remove all entries whose seed peaks (parsed from their key) are within
    // tolKeV of ANY energy in peaks.  Used so a manual grouped fit can evict
    // the individual AutoFit entries it supersedes.
    void RemoveOverlapping(const std::vector<double>& peaks, double tolKeV) {
        auto it = entries_.begin();
        while (it != entries_.end()) {
            std::vector<double> seedEs;
            std::istringstream ss(it->first);
            std::string tok;
            while (std::getline(ss, tok, '_')) {
                try { seedEs.push_back(std::stod(tok)); } catch (...) {}
            }
            bool overlap = false;
            for (double s : seedEs)
                for (double m : peaks)
                    if (std::abs(s - m) < tolKeV) { overlap = true; break; }
            if (overlap) {
                Debug::Log(Debug::DB, "[FitDatabase] Removing overlapping entry: " + it->first);
                it = entries_.erase(it);
            } else {
                ++it;
            }
        }
    }

    void StoreIfBetter(const std::string& key, const FitEntry& candidate) {
        auto it = entries_.find(key);
        if (it == entries_.end()) {
            Debug::Log(Debug::DB, "[FitDatabase] New entry stored: " + key +
                       "  chi2/ndf=" + std::to_string(candidate.chi2ndf) +
                       "  maxPull=" + std::to_string(candidate.maxPull));
            entries_[key] = candidate;
        } else if (CompositeScore(candidate) < CompositeScore(it->second)) {
            Debug::Log(Debug::DB, "[FitDatabase] Improved entry: " + key +
                       "  score " + std::to_string(CompositeScore(it->second)) +
                       " -> " + std::to_string(CompositeScore(candidate)));
            FitEntry updated = candidate;
            if (updated.label.empty())          updated.label          = it->second.label;
            if (updated.classification.empty()) updated.classification = it->second.classification;
            entries_[key] = std::move(updated);
        } else {
            Debug::Log(Debug::DB, "[FitDatabase] Kept existing entry: " + key +
                       "  cached score=" + std::to_string(CompositeScore(it->second)) +
                       "  new score=" + std::to_string(CompositeScore(candidate)));
        }
    }

    // Load persisted fit cache from file. Silently succeeds if file doesn't
    // exist yet (first run). Format (backward-compatible):
    //   # META bg_subtracted 0|1 bg_iterations N
    //   key npar chi2ndf p0 p1 ... pN-1 [residualRMS maxPull]
    bool Load(const std::string& filename) {
        std::ifstream in(filename);
        if (!in.is_open()) {
            Debug::Log(Debug::DB, "[FitDatabase] No cache file found at " + filename + " — starting fresh");
            return false;
        }
        int count = 0;
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty()) continue;
            // META header: "# META key value ..."
            if (line.size() >= 7 && line.substr(0, 7) == "# META ") {
                std::istringstream ms(line.substr(7));
                std::string tok;
                while (ms >> tok) {
                    if (tok == "bg_subtracted") {
                        int v = 0; ms >> v; bgSubtracted = (v != 0);
                    } else if (tok == "bg_iterations") {
                        ms >> bgIterations;
                    } else if (tok == "root_file") {
                        // rest of line is the path (may contain spaces)
                        std::string rest;
                        std::getline(ms, rest);
                        if (!rest.empty() && rest[0] == ' ') rest = rest.substr(1);
                        rootFile = rest;
                        break;  // root_file is always last — nothing follows
                    }
                }
                continue;
            }
            if (line[0] == '#') continue;
            std::istringstream ss(line);
            FitEntry e;
            int npar = 0;
            if (!(ss >> e.key >> npar >> e.chi2ndf)) continue;
            e.params.resize(npar);
            bool ok = true;
            for (int i = 0; i < npar; i++) {
                if (!(ss >> e.params[i])) { ok = false; break; }
            }
            if (!ok) continue;
            // Optional trailing fields — silently ignored on old cache files
            ss >> e.residualRMS;
            ss >> e.maxPull;
            {
                std::string meth;
                if (ss >> meth && meth != "-") e.fitMethod = meth;
            }
            ss >> e.xlo;
            ss >> e.xhi;
            // Optional label and classification — empty on old cache files
            {
                std::string lbl, cls;
                if (ss >> lbl && lbl != "-") e.label = lbl;
                if (ss >> cls && cls != "-") e.classification = cls;
            }
            // Optional parameter errors — absent on old cache files
            {
                int ne = 0;
                if (ss >> ne && ne == npar) {
                    e.paramErrors.resize(ne);
                    for (int i = 0; i < ne; i++) {
                        if (!(ss >> e.paramErrors[i])) { e.paramErrors.clear(); break; }
                    }
                }
            }
            // Optional per-Gaussian labels / classifications
            {
                int npl = 0;
                if (ss >> npl && npl > 0) {
                    e.peakLabels.resize(npl);
                    for (int i = 0; i < npl; i++) {
                        std::string t; if (ss >> t && t != "-") e.peakLabels[i] = t;
                    }
                }
                int npc = 0;
                if (ss >> npc && npc > 0) {
                    e.peakClassifications.resize(npc);
                    for (int i = 0; i < npc; i++) {
                        std::string t; if (ss >> t && t != "-") e.peakClassifications[i] = t;
                    }
                }
            }
            entries_[e.key] = std::move(e);
            ++count;
        }
        Debug::Log(Debug::DB, "[FitDatabase] Loaded " + std::to_string(count) +
                   " cached fits from " + filename +
                   "  bg_subtracted=" + std::to_string(bgSubtracted) +
                   "  bg_iterations=" + std::to_string(bgIterations));
        return true;
    }

    // Write all entries atomically. Format:
    //   # META bg_subtracted 0|1 bg_iterations N
    //   key npar chi2ndf p0 p1 ... pN-1 residualRMS maxPull
    void Save(const std::string& filename) const {
        std::string tmp = filename + ".tmp";
        std::ofstream out(tmp);
        if (!out.is_open()) {
            std::cerr << "[FitDatabase] Cannot write cache to " << tmp << "\n";
            return;
        }
        out << "# META bg_subtracted " << (bgSubtracted ? 1 : 0)
            << " bg_iterations " << bgIterations;
        if (!rootFile.empty())
            out << " root_file " << rootFile;
        out << "\n";
        out << std::fixed << std::setprecision(8);
        for (const auto& kv : entries_) {
            const auto& e = kv.second;
            out << e.key << " " << e.params.size() << " " << e.chi2ndf;
            for (double p : e.params) out << " " << p;
            out << " " << e.residualRMS << " " << e.maxPull;
            out << " " << (e.fitMethod.empty() ? "-" : e.fitMethod);
            out << " " << e.xlo << " " << e.xhi;
            out << " " << (e.label.empty()          ? "-" : e.label);
            out << " " << (e.classification.empty() ? "-" : e.classification);
            // Parameter errors (ne followed by ne values)
            int ne = (int)e.paramErrors.size();
            out << " " << ne;
            for (int i = 0; i < ne; i++) out << " " << e.paramErrors[i];
            // Per-Gaussian labels / classifications (npl l0 l1 ... npc c0 c1 ...)
            int npl = (int)e.peakLabels.size();
            out << " " << npl;
            for (int i = 0; i < npl; i++)
                out << " " << (e.peakLabels[i].empty() ? "-" : e.peakLabels[i]);
            int npc = (int)e.peakClassifications.size();
            out << " " << npc;
            for (int i = 0; i < npc; i++)
                out << " " << (e.peakClassifications[i].empty() ? "-" : e.peakClassifications[i]);
            out << "\n";
        }
        out.close();
        std::rename(tmp.c_str(), filename.c_str());
        Debug::Log(Debug::DB, "[FitDatabase] Saved " + std::to_string(entries_.size()) +
                   " cached fits to " + filename);
    }

    const std::map<std::string, FitEntry>& GetEntries() const { return entries_; }

    // Remove a single entry by key. Returns true if the key existed.
    bool Remove(const std::string& key) {
        auto it = entries_.find(key);
        if (it == entries_.end()) return false;
        entries_.erase(it);
        return true;
    }

private:
    std::map<std::string, FitEntry> entries_;
};

#endif // FITDATABASE_H
