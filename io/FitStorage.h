#ifndef FITSTORAGE_H
#define FITSTORAGE_H

#include <map>
#include <vector>
#include <string>

#include "PeakFitResult.h"

class FitStorage {
private:
    std::map<std::string, std::vector<PeakFitResult>> data;

public:
    void Save(const std::string& name, const PeakFitResult& p) {
        data[name].push_back(p);
    }

    const std::vector<PeakFitResult>* Load(const std::string& name) const {
        if (data.count(name)) return &data.at(name);
        return nullptr;
    }
};

#endif