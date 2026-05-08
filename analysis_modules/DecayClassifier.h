#ifndef DECAYCLASSIFIER_H
#define DECAYCLASSIFIER_H

#include <vector>

enum DecayType {
    PARENT,
    DAUGHTER,
    BACKGROUND
};

class DecayClassifier {
public:
    DecayType Classify(const std::vector<double>& t,
                       const std::vector<double>& c) {

        if (t.size() < 3) return BACKGROUND;

        if (c.front() > c.back()) return PARENT;
        if (c[1] > c.front()) return DAUGHTER;

        return BACKGROUND;
    }
};

#endif