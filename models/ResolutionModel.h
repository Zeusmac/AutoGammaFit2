#ifndef RESOLUTIONMODEL_H
#define RESOLUTIONMODEL_H

#include <cmath>
#include <string>
#include "Debug.h"

class ResolutionModel {
public:
    double a = 0.5;
    double b = 0.02;
    double c = 0.0;  // kept for cache compat; not used in FWHM (Fano model has no E² term)

    // optional smoothing factor
    double alpha = 0.2;

    // Fano model: FWHM² = η² + F·ε·E  →  FWHM(E) = sqrt(a + b·E)
    // a = electronic noise², b = Fano × charge-creation energy (~0.002 for HPGe)
    double FWHM(double E) const {
        double fwhm_sq = a + b * E;
        return (fwhm_sq > 0.0) ? std::sqrt(fwhm_sq) : 0.0;
    }

    double Sigma(double E) const {
        return FWHM(E) / 2.3548200450309493;
    }

    // -----------------------------
    // AUTO UPDATE FROM FITS (Quadrature Based)
    // -----------------------------
    void UpdateFromQuadratureFit(double a_new, double b_new) {
        Debug::Log(Debug::RESMODEL,
            "Update: old a=" + std::to_string(a) + "  b=" + std::to_string(b));
        Debug::Log(Debug::RESMODEL,
            "        raw a=" + std::to_string(a_new) + "  b=" + std::to_string(b_new));

        a = (1.0 - alpha) * a + alpha * a_new;
        b = (1.0 - alpha) * b + alpha * b_new;
        c = 0.0;

        Debug::Log(Debug::RESMODEL,
            "     smooth a=" + std::to_string(a) +
            "  b=" + std::to_string(b) +
            "  alpha=" + std::to_string(alpha));
    }
};

#endif