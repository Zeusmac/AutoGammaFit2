#ifndef RESOLUTIONMODEL_H
#define RESOLUTIONMODEL_H

#include <cmath>
#include <string>
#include "Debug.h"

class ResolutionModel {
public:
    double a = 0.5;
    double b = 0.02;
    double c = 1e-8;

    // optional smoothing factor
    double alpha = 0.2;

    // Model: FWHM^2 = a + b*E + c*E^2  (standard quadrature detector model)
    double FWHM(double E) const {
        double fwhm_sq = a + b * E + c * E * E;
        return (fwhm_sq > 0) ? std::sqrt(fwhm_sq) : 0.0;
    }

    double Sigma(double E) const {
        return FWHM(E) / 2.3548200450309493;
    }

    // -----------------------------
    // AUTO UPDATE FROM FITS (Quadrature Based)
    // -----------------------------
    void UpdateFromQuadratureFit(double a_new, double b_new, double c_new) {
        Debug::Log(Debug::RESMODEL,
            "Update: old a=" + std::to_string(a) +
            "  b=" + std::to_string(b) +
            "  c=" + std::to_string(c));
        Debug::Log(Debug::RESMODEL,
            "        raw a=" + std::to_string(a_new) +
            "  b=" + std::to_string(b_new) +
            "  c=" + std::to_string(c_new));

        a = (1 - alpha) * a + alpha * a_new;
        b = (1 - alpha) * b + alpha * b_new;
        c = (1 - alpha) * c + alpha * c_new;

        Debug::Log(Debug::RESMODEL,
            "     smooth a=" + std::to_string(a) +
            "  b=" + std::to_string(b) +
            "  c=" + std::to_string(c) +
            "  alpha=" + std::to_string(alpha));
    }
};

#endif