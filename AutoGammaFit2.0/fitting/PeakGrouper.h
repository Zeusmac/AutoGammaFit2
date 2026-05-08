#ifndef PEAKGROUPER_H
#define PEAKGROUPER_H

#include <vector>
#include <cmath>
#include <iostream>
#include <string>
#include "Debug.h"

struct PeakGroup {
    std::vector<double> energies;
};

class PeakGrouper {
public:
    template<typename ResModel>
    static std::vector<PeakGroup>
    Group(const std::vector<double>& peaks, const ResModel& res) {

        std::vector<PeakGroup> groups;
        if (peaks.empty()) return groups;

        PeakGroup current;
        current.energies.push_back(peaks[0]);

        for (size_t i = 1; i < peaks.size(); i++) {

            double prev = current.energies.back();
            double E = peaks[i];

            double sigma = res.Sigma(prev);
            double gap   = std::fabs(E - prev);
            if (gap < 4.0 * sigma) {
                Debug::Log(Debug::GROUPER,
                    "Merge " + std::to_string(prev) + " + " + std::to_string(E) +
                    "  gap=" + std::to_string(gap) +
                    "  4sigma=" + std::to_string(4.0 * sigma));
                current.energies.push_back(E);
            } else {
                Debug::Log(Debug::GROUPER,
                    "Split at " + std::to_string(E) +
                    "  gap=" + std::to_string(gap) +
                    "  4sigma=" + std::to_string(4.0 * sigma) +
                    "  -> new group");
                groups.push_back(current);
                current = PeakGroup();
                current.energies.push_back(E);
            }
        }

        groups.push_back(current);
        Debug::Log(Debug::GROUPER,
            "Total groups=" + std::to_string(groups.size()) +
            "  from " + std::to_string(peaks.size()) + " peaks");
        return groups;
    }
};

#endif