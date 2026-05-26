#ifndef UISETTINGS_H
#define UISETTINGS_H

#include <string>
#include <unordered_map>
#include <fstream>
#include <sstream>

// ─────────────────────────────────────────────────────────────────────────────
// UISettings  -  flat key/value store persisted to gamma_gui.conf
//
// Format: one "key=value" per line; blank lines and '#' comments ignored.
// All values stored as strings; typed helpers below convert on access.
//
// Usage:
//   UISettings s; s.Load("gamma_gui.conf");
//   s.Set("last_file", "/data/run.root");
//   int n = s.GetInt("bg_iters", 5);
//   s.Save("gamma_gui.conf");
// ─────────────────────────────────────────────────────────────────────────────
class UISettings {
public:
    void Load(const std::string& path) {
        path_ = path;
        std::ifstream f(path);
        if (!f) return;
        std::string line;
        while (std::getline(f, line)) {
            if (line.empty() || line[0] == '#') continue;
            auto eq = line.find('=');
            if (eq == std::string::npos) continue;
            map_[line.substr(0, eq)] = line.substr(eq + 1);
        }
    }

    bool Save(const std::string& path = "") const {
        const std::string& dest = path.empty() ? path_ : path;
        if (dest.empty()) return false;
        std::string tmp = dest + ".tmp";
        std::ofstream f(tmp);
        if (!f) return false;
        for (const auto& [k, v] : map_)
            f << k << "=" << v << "\n";
        f.close();
        return std::rename(tmp.c_str(), dest.c_str()) == 0;
    }

    void Set(const std::string& key, const std::string& val) { map_[key] = val; }
    void Set(const std::string& key, int val)    { map_[key] = std::to_string(val); }
    void Set(const std::string& key, double val) {
        std::ostringstream ss; ss << val; map_[key] = ss.str();
    }
    void Set(const std::string& key, bool val)   { map_[key] = val ? "1" : "0"; }

    std::string Get(const std::string& key, const std::string& def = "") const {
        auto it = map_.find(key);
        return it != map_.end() ? it->second : def;
    }
    int    GetInt(const std::string& key, int def = 0) const {
        auto it = map_.find(key);
        if (it == map_.end()) return def;
        try { return std::stoi(it->second); } catch (...) { return def; }
    }
    double GetDouble(const std::string& key, double def = 0.0) const {
        auto it = map_.find(key);
        if (it == map_.end()) return def;
        try { return std::stod(it->second); } catch (...) { return def; }
    }
    bool   GetBool(const std::string& key, bool def = false) const {
        auto it = map_.find(key);
        if (it == map_.end()) return def;
        return it->second == "1" || it->second == "true";
    }

    bool Has(const std::string& key) const { return map_.count(key) > 0; }

private:
    std::string path_;
    std::unordered_map<std::string, std::string> map_;
};

#endif // UISETTINGS_H
