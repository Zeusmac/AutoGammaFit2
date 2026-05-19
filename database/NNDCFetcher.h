#ifndef NNDCFETCHER_H
#define NNDCFETCHER_H

#include "NuclearData.h"

#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <map>
#include <algorithm>
#include <functional>
#include <cstdlib>
#include <sys/stat.h>

// ─────────────────────────────────────────────────────────────────────────────
// NNDCFetcher — fetch nuclear data from IAEA LiveChart REST API
//
// API base: https://nds.iaea.org/relnsd/v1/data
// Uses wget (curl NOT available)
// Offline mode: LoadFromFile() reads locally-saved CSV exports
// ─────────────────────────────────────────────────────────────────────────────

class NNDCFetcher {
public:

    // ─── Public API ──────────────────────────────────────────────────────────

    // Fetch all available data for a given isotope. Returns true if at least
    // ground-state data was successfully retrieved (from network or cache).
    static bool Fetch(int A, const std::string& symbol,
                      NucIsotope& iso,
                      const std::string& cacheDir)
    {
        iso.A = A; iso.symbol = symbol;
        iso.Z = NucSymbolToZ(symbol);
        iso.N = A - iso.Z;
        iso._valid = false;

        EnsureDir(cacheDir);
        std::string cachePath = cacheDir + "/" + std::to_string(A) + symbol + ".nucdat";

        // Try cache first
        if (LoadCache(cachePath, iso)) return true;

        // Fetch from IAEA
        std::string base = "https://nds.iaea.org/relnsd/v1/data?fields=";
        std::string nuc  = "&nuclides=" + std::to_string(A) + symbol;

        // Ground states (half-life, jpi, decay modes)
        {
            std::string url = base + "ground_states" + nuc;
            std::string csv = WgetToString(url);
            if (!csv.empty()) ParseGroundStates(csv, iso);
        }
        // Gammas
        {
            std::string url = base + "gammas" + nuc;
            std::string csv = WgetToString(url);
            if (!csv.empty()) ParseGammas(csv, iso);
        }
        // Levels
        {
            std::string url = base + "levels" + nuc;
            std::string csv = WgetToString(url);
            if (!csv.empty()) ParseLevels(csv, iso);
        }
        // Beta branches
        {
            std::string url = base + "beta_minus" + nuc;
            std::string csv = WgetToString(url);
            if (!csv.empty()) ParseBeta(csv, iso);
        }

        if (!iso.jpi.empty() || !iso.gammas.empty() || iso.halflife_s >= 0) {
            iso._valid = true;
            SaveCache(cachePath, iso);
        }
        return iso._valid;
    }

    // Load a locally-saved IAEA CSV file for a single isotope.
    // Detects the file type from the CSV header row and calls the appropriate
    // parse method. Saves a new cache file on success.
    static bool LoadFromFile(const std::string& path,
                             int A, const std::string& symbol,
                             NucIsotope& iso,
                             const std::string& cacheDir)
    {
        std::ifstream f(path);
        if (!f.is_open()) return false;
        std::string content((std::istreambuf_iterator<char>(f)),
                             std::istreambuf_iterator<char>());
        if (content.empty() || content.find(',') == std::string::npos) return false;

        iso.A = A; iso.symbol = symbol;
        iso.Z = NucSymbolToZ(symbol); iso.N = A - iso.Z;
        iso._valid = false;

        // Detect file type from header row
        auto nl = content.find('\n');
        std::string firstLine = (nl != std::string::npos)
                                ? content.substr(0, nl) : content;
        bool parsed = false;
        if (firstLine.find("half_life") != std::string::npos ||
            firstLine.find("jp")        != std::string::npos)
        { ParseGroundStates(content, iso); parsed = true; }
        if (firstLine.find("start_level_energy") != std::string::npos ||
            firstLine.find("gamma_") != std::string::npos ||
            firstLine.find("energy") != std::string::npos)
        { ParseGammas(content, iso); parsed = true; }
        if (firstLine.find("level_energy") != std::string::npos)
        { ParseLevels(content, iso); parsed = true; }
        if (firstLine.find("endpoint_energy") != std::string::npos)
        { ParseBeta(content, iso); parsed = true; }

        if (!parsed) return false;
        iso._valid = true;

        // Save cache
        EnsureDir(cacheDir);
        std::string cachePath = cacheDir + "/" + std::to_string(A) + symbol + ".nucdat";
        SaveCache(cachePath, iso);
        return true;
    }

    // Load from cache file (returns false if not found or invalid)
    static bool LoadCache(const std::string& path, NucIsotope& iso) {
        std::ifstream f(path);
        if (!f.is_open()) return false;

        std::string line;
        bool valid = false;
        iso.levels.clear();
        iso.gammas.clear();
        iso.betaBranches.clear();

        while (std::getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            std::istringstream ss(line);
            std::string tag;
            ss >> tag;
            if (tag == "A")          { ss >> iso.A; }
            else if (tag == "Z")     { ss >> iso.Z; }
            else if (tag == "N")     { ss >> iso.N; }
            else if (tag == "SYM")   { ss >> iso.symbol; }
            else if (tag == "JPI")   { std::getline(ss, iso.jpi); TrimInPlace(iso.jpi); }
            else if (tag == "HL_S")  { ss >> iso.halflife_s; }
            else if (tag == "HL_STR"){ std::getline(ss, iso.hl_str); TrimInPlace(iso.hl_str); }
            else if (tag == "QB-")   { ss >> iso.qBetaMinus; }
            else if (tag == "DM1")   { std::getline(ss, iso.decayMode1); TrimInPlace(iso.decayMode1); }
            else if (tag == "DM1F")  { ss >> iso.decay1frac; }
            else if (tag == "DATE")  { ss >> iso.fetchDate; }
            else if (tag == "VALID") { int v = 0; ss >> v; valid = (v != 0); }
            else if (tag == "LEVEL") {
                NucLevel lv;
                ss >> lv.energy >> lv.halflife_s;
                std::getline(ss, lv.jpi); TrimInPlace(lv.jpi);
                iso.levels.push_back(lv);
            }
            else if (tag == "GAMMA") {
                NucGamma gm;
                ss >> gm.start_level >> gm.energy >> gm.intensity;
                iso.gammas.push_back(gm);
            }
            else if (tag == "BETA") {
                NucBetaBranch bb;
                ss >> bb.endpoint >> bb.intensity;
                iso.betaBranches.push_back(bb);
            }
        }
        iso._valid = valid;
        return valid;
    }

    // Save to cache file
    static bool SaveCache(const std::string& path, const NucIsotope& iso) {
        std::ofstream f(path);
        if (!f.is_open()) return false;
        f << "# NucIsotope cache — do not edit manually\n";
        f << "A "      << iso.A          << "\n";
        f << "Z "      << iso.Z          << "\n";
        f << "N "      << iso.N          << "\n";
        f << "SYM "    << iso.symbol     << "\n";
        f << "JPI "    << iso.jpi        << "\n";
        f << "HL_S "   << iso.halflife_s << "\n";
        f << "HL_STR " << iso.hl_str     << "\n";
        f << "QB- "    << iso.qBetaMinus << "\n";
        f << "DM1 "    << iso.decayMode1 << "\n";
        f << "DM1F "   << iso.decay1frac << "\n";
        f << "DATE "   << iso.fetchDate  << "\n";
        f << "VALID "  << (iso._valid ? 1 : 0) << "\n";
        for (const auto& lv : iso.levels)
            f << "LEVEL " << lv.energy << " " << lv.halflife_s << " " << lv.jpi << "\n";
        for (const auto& gm : iso.gammas)
            f << "GAMMA " << gm.start_level << " " << gm.energy << " " << gm.intensity << "\n";
        for (const auto& bb : iso.betaBranches)
            f << "BETA " << bb.endpoint << " " << bb.intensity << "\n";
        return true;
    }

    // ─── Parse methods (also called by LoadFromFile) ─────────────────────────

    static void ParseGroundStates(const std::string& csv, NucIsotope& iso) {
        auto rows = ParseCSV(csv);
        if (rows.size() < 2) return;
        const auto& hdr = rows[0];
        auto get = [&](const std::string& field) -> std::string {
            for (size_t i = 0; i < hdr.size(); i++) {
                if (hdr[i] == field && rows.size() > 1 && i < rows[1].size())
                    return rows[1][i];
            }
            return "";
        };

        iso.jpi        = get("jp");
        iso.hl_str     = get("half_life");
        iso.decayMode1 = get("decay_1");

        // Parse decay_1_%
        std::string d1pct = get("decay_1_%");
        if (!d1pct.empty()) {
            std::string d1clean = d1pct;
            d1clean.erase(
                std::remove_if(d1clean.begin(), d1clean.end(),
                    [](char c){ return c=='<'||c=='>'||c=='~'||c==' '; }),
                d1clean.end());
            try { iso.decay1frac = std::stod(d1clean); } catch (...) {}
        }

        // Parse halflife in seconds
        std::string hl_s_str = get("half_life_sec");
        if (hl_s_str.empty()) hl_s_str = get("halflife_sec");
        if (!hl_s_str.empty()) {
            try { iso.halflife_s = std::stod(hl_s_str); } catch (...) {}
        } else if (iso.hl_str == "STABLE" || iso.hl_str == "stable") {
            iso.halflife_s = 0.0;
        }

        // Q-value beta minus
        std::string qbm = get("q_beta_minus");
        if (!qbm.empty()) { try { iso.qBetaMinus = std::stod(qbm); } catch (...) {} }
    }

    static void ParseGammas(const std::string& csv, NucIsotope& iso) {
        auto rows = ParseCSV(csv);
        if (rows.size() < 2) return;
        const auto& hdr = rows[0];

        // Find column indices
        int iStartLv = ColIdx(hdr, "start_level_energy");
        int iEnergy  = ColIdx(hdr, "energy");
        if (iEnergy < 0) iEnergy = ColIdx(hdr, "gamma_energy");
        int iIntens  = ColIdx(hdr, "intensity");
        if (iIntens < 0) iIntens = ColIdx(hdr, "gamma_intensity");

        for (size_t r = 1; r < rows.size(); r++) {
            const auto& row = rows[r];
            if (row.empty()) continue;
            NucGamma gm;
            if (iStartLv >= 0 && (size_t)iStartLv < row.size())
                try { gm.start_level = std::stod(row[iStartLv]); } catch (...) {}
            if (iEnergy >= 0 && (size_t)iEnergy < row.size())
                try { gm.energy = std::stod(row[iEnergy]); } catch (...) {}
            if (iIntens >= 0 && (size_t)iIntens < row.size())
                try { gm.intensity = std::stod(row[iIntens]); } catch (...) {}
            if (gm.energy > 0) iso.gammas.push_back(gm);
        }
    }

    static void ParseLevels(const std::string& csv, NucIsotope& iso) {
        auto rows = ParseCSV(csv);
        if (rows.size() < 2) return;
        const auto& hdr = rows[0];

        int iEnergy  = ColIdx(hdr, "level_energy");
        if (iEnergy < 0) iEnergy = ColIdx(hdr, "energy");
        int iJpi     = ColIdx(hdr, "jp");
        int iHL      = ColIdx(hdr, "half_life_sec");
        if (iHL < 0)  iHL = ColIdx(hdr, "halflife_sec");

        for (size_t r = 1; r < rows.size(); r++) {
            const auto& row = rows[r];
            if (row.empty()) continue;
            NucLevel lv;
            if (iEnergy >= 0 && (size_t)iEnergy < row.size())
                try { lv.energy = std::stod(row[iEnergy]); } catch (...) {}
            if (iJpi >= 0 && (size_t)iJpi < row.size())
                lv.jpi = row[iJpi];
            if (iHL >= 0 && (size_t)iHL < row.size())
                try { lv.halflife_s = std::stod(row[iHL]); } catch (...) {}
            iso.levels.push_back(lv);
        }
    }

    static void ParseBeta(const std::string& csv, NucIsotope& iso) {
        auto rows = ParseCSV(csv);
        if (rows.size() < 2) return;
        const auto& hdr = rows[0];

        int iEndpt   = ColIdx(hdr, "endpoint_energy");
        int iIntens  = ColIdx(hdr, "intensity");
        if (iIntens < 0) iIntens = ColIdx(hdr, "beta_intensity");

        for (size_t r = 1; r < rows.size(); r++) {
            const auto& row = rows[r];
            if (row.empty()) continue;
            NucBetaBranch bb;
            if (iEndpt >= 0 && (size_t)iEndpt < row.size())
                try { bb.endpoint = std::stod(row[iEndpt]); } catch (...) {}
            if (iIntens >= 0 && (size_t)iIntens < row.size())
                try { bb.intensity = std::stod(row[iIntens]); } catch (...) {}
            if (bb.endpoint > 0 || bb.intensity > 0) iso.betaBranches.push_back(bb);
        }
    }

    // ─── Helpers ─────────────────────────────────────────────────────────────

    static void EnsureDir(const std::string& dir) {
        mkdir(dir.c_str(), 0755);
    }

private:

    static std::string WgetToString(const std::string& url) {
        // Write to a temp file, then read it
        std::string tmpfile = "/tmp/nndc_fetch_tmp_" + std::to_string(getpid()) + ".csv";
        std::string cmd = "wget -q -O " + tmpfile + " \"" + url + "\" 2>/dev/null";
        int ret = system(cmd.c_str());
        if (ret != 0) { remove(tmpfile.c_str()); return ""; }
        std::ifstream f(tmpfile);
        std::string content((std::istreambuf_iterator<char>(f)),
                             std::istreambuf_iterator<char>());
        remove(tmpfile.c_str());
        return content;
    }

    // Parse CSV into a 2D vector of strings (first row = headers)
    static std::vector<std::vector<std::string>> ParseCSV(const std::string& csv) {
        std::vector<std::vector<std::string>> result;
        std::istringstream ss(csv);
        std::string line;
        while (std::getline(ss, line)) {
            // Strip trailing \r
            if (!line.empty() && line.back() == '\r') line.pop_back();
            if (line.empty()) continue;
            std::vector<std::string> row;
            std::istringstream ls(line);
            std::string cell;
            while (std::getline(ls, cell, ',')) {
                // Trim whitespace
                size_t s = cell.find_first_not_of(" \t");
                size_t e = cell.find_last_not_of(" \t");
                row.push_back((s == std::string::npos) ? "" : cell.substr(s, e - s + 1));
            }
            result.push_back(row);
        }
        return result;
    }

    static int ColIdx(const std::vector<std::string>& hdr, const std::string& name) {
        for (int i = 0; i < (int)hdr.size(); i++)
            if (hdr[i] == name) return i;
        return -1;
    }

    static void TrimInPlace(std::string& s) {
        size_t first = s.find_first_not_of(" \t\r\n");
        if (first == std::string::npos) { s.clear(); return; }
        size_t last  = s.find_last_not_of(" \t\r\n");
        s = s.substr(first, last - first + 1);
    }
};

#endif // NNDCFETCHER_H
