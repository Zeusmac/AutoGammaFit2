#ifndef FITDATABASE_H
#define FITDATABASE_H

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <sstream>
#include <iomanip>
#include <limits>

#include "TF1.h"
#include "TFitResult.h"
#include "TFitResultPtr.h"
#include "Debug.h"

struct FitEntry {
    std::string key;
    std::vector<double> params;
    double chi2ndf = std::numeric_limits<double>::max();
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

    static std::string MakeKey(const std::vector<double>& peaks) {
        std::ostringstream oss;
        for (size_t i = 0; i < peaks.size(); i++) {
            if (i) oss << "_";
            oss << std::fixed << std::setprecision(1) << peaks[i];
        }
        return oss.str();
    }

    // Like SeedFromDB but always seeds regardless of useCachedSeeds.
    // Used by cache-only mode where we never run MIGRAD at all.
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
                   "  cached chi2/ndf=" + std::to_string(it->second.chi2ndf));
        for (int i = 0; i < f->GetNpar(); i++)
            f->SetParameter(i, params[i]);
        return true;
    }

    static FitEntry MakeEntry(const std::string& key,
                               TFitResultPtr& r,
                               TF1* f,
                               int /*nUsed*/)
    {
        FitEntry e;
        e.key     = key;
        e.chi2ndf = (r.Get() && r->Ndf() > 0)
                    ? r->Chi2() / r->Ndf()
                    : std::numeric_limits<double>::max();
        int n = f->GetNpar();
        e.params.resize(n);
        for (int i = 0; i < n; i++)
            e.params[i] = f->GetParameter(i);
        return e;
    }

    void StoreIfBetter(const std::string& key, const FitEntry& candidate) {
        auto it = entries_.find(key);
        if (it == entries_.end()) {
            Debug::Log(Debug::DB, "[FitDatabase] New entry stored: " + key +
                       "  chi2/ndf=" + std::to_string(candidate.chi2ndf));
            entries_[key] = candidate;
        } else if (candidate.chi2ndf < it->second.chi2ndf) {
            Debug::Log(Debug::DB, "[FitDatabase] Improved entry: " + key +
                       "  " + std::to_string(it->second.chi2ndf) +
                       " -> " + std::to_string(candidate.chi2ndf));
            entries_[key] = candidate;
        } else {
            Debug::Log(Debug::DB, "[FitDatabase] Kept existing entry: " + key +
                       "  cached=" + std::to_string(it->second.chi2ndf) +
                       "  new=" + std::to_string(candidate.chi2ndf));
        }
    }

    // Load persisted fit cache from file. Silently succeeds if file doesn't
    // exist yet (first run).
    bool Load(const std::string& filename) {
        std::ifstream in(filename);
        if (!in.is_open()) {
            Debug::Log(Debug::DB, "[FitDatabase] No cache file found at " + filename + " — starting fresh");
            return false;
        }
        int count = 0;
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            FitEntry e;
            int npar = 0;
            if (!(ss >> e.key >> npar >> e.chi2ndf)) continue;
            e.params.resize(npar);
            for (int i = 0; i < npar; i++) {
                if (!(ss >> e.params[i])) { npar = 0; break; }
            }
            if (npar > 0) { entries_[e.key] = std::move(e); ++count; }
        }
        Debug::Log(Debug::DB, "[FitDatabase] Loaded " + std::to_string(count) +
                   " cached fits from " + filename);
        return true;
    }

    // Write all entries to file atomically via a temp file so a crash
    // mid-write doesn't corrupt the cache.
    void Save(const std::string& filename) const {
        std::string tmp = filename + ".tmp";
        std::ofstream out(tmp);
        if (!out.is_open()) {
            std::cerr << "[FitDatabase] Cannot write cache to " << tmp << "\n";
            return;
        }
        out << std::fixed << std::setprecision(8);
        for (const auto& kv : entries_) {
            const auto& e = kv.second;
            out << e.key << " " << e.params.size() << " " << e.chi2ndf;
            for (double p : e.params) out << " " << p;
            out << "\n";
        }
        out.close();
        std::rename(tmp.c_str(), filename.c_str());
        Debug::Log(Debug::DB, "[FitDatabase] Saved " + std::to_string(entries_.size()) +
                   " cached fits to " + filename);
    }

private:
    std::map<std::string, FitEntry> entries_;
};

#endif // FITDATABASE_H
