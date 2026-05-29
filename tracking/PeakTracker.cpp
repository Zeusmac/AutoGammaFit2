#include "TF1.h"
#include "TGraphErrors.h"
#include "PeakTracker.h"
#include "ResolutionModel.h"
#include "RootFileManager.h"
#include "Debug.h"

void PeakTracker::Add(double energy,
                      double time,
                      double counts_val,
                      double error,
                      double sigma,
                      double sigma_err,
                      int    fit_status,
                      double fit_edm)
{
    if (!std::isfinite(energy) ||
        !std::isfinite(time) ||
        !std::isfinite(counts_val))
        return;

    double key = GetKey(energy);

    PeakPoint p;
    p.energy     = energy;
    p.time       = time;
    p.counts     = counts_val;
    p.error      = error;
    p.sigma      = sigma;
    p.sigma_err  = sigma_err;
    p.fit_status = fit_status;
    p.fit_edm    = fit_edm;

    data[key].push_back(p);

    Debug::Log(Debug::TRACKER,
        "Add E=" + std::to_string(energy) +
        "  sigma=" + std::to_string(sigma) +
        "  sigma_err=" + std::to_string(sigma_err) +
        "  status=" + std::to_string(fit_status) +
        "  edm=" + std::to_string(fit_edm));
}

double PeakTracker::GetKey(double energy, double tol_override)
{
    double tol = (tol_override > 0) ? tol_override : energy_tol;

    if (!std::isfinite(energy)) return -1;

    for (auto &kv : data) {
        if (std::fabs(kv.first - energy) < tol)
            return kv.first;
    }

    return energy;
}


void PeakTracker::FitResolutionModel(ResolutionModel& model, TFile* fout, const std::string& tag)
{
    std::vector<double> E, S, dS;

    for (const auto& kv : data) {

        double energy = kv.first;
        const auto& vec = kv.second;

        double wsum = 0.0;
        double mean = 0.0;

        for (const auto& p : vec) {

            if (p.sigma <= 0 || p.sigma_err <= 0) continue;
            if (2.355 * p.sigma > 15.0) {
                Debug::Log(Debug::TRACKER,
                    "FitRes E=" + std::to_string(energy) +
                    "  skip: FWHM=" + std::to_string(2.355*p.sigma) + " > 15 keV");
                continue;
            }
            if (p.sigma_err < 0.001 * p.sigma) {
                Debug::Log(Debug::TRACKER,
                    "FitRes E=" + std::to_string(energy) +
                    "  skip: boundary-convergence artifact");
                continue;
            }
            if (p.fit_edm > 0.01) {
                Debug::Log(Debug::TRACKER,
                    "FitRes E=" + std::to_string(energy) +
                    "  skip: EDM=" + std::to_string(p.fit_edm) + " > 0.01");
                continue;
            }

            double w = 1.0 / (p.sigma_err * p.sigma_err);
            mean += w * p.sigma;
            wsum += w;
        }

        if (wsum == 0) continue;

        mean /= wsum;
        double err = std::sqrt(1.0 / wsum);

        constexpr double kSigToFWHM = 2.3548200450309493;
        double fwhm     = mean * kSigToFWHM;
        double fwhm_err = err  * kSigToFWHM;

        if (fwhm > 15.0) {
            Debug::Log(Debug::TRACKER,
                "FitRes E=" + std::to_string(energy) +
                "  skip post-avg: FWHM=" + std::to_string(fwhm) + " > 15 keV");
            continue;
        }
        if (fwhm_err / fwhm > 0.4) {
            Debug::Log(Debug::TRACKER,
                "FitRes E=" + std::to_string(energy) +
                "  skip post-avg: rel_err=" + std::to_string(fwhm_err/fwhm) + " > 0.4");
            continue;
        }

        Debug::Log(Debug::TRACKER,
            "FitRes using E=" + std::to_string(energy) +
            "  FWHM=" + std::to_string(fwhm) +
            "  err=" + std::to_string(fwhm_err));

        E.push_back(energy);
        S.push_back(fwhm);
        dS.push_back(fwhm_err);
    }

    Debug::Log(Debug::TRACKER,
        "FitResolutionModel tag=" + tag +
        "  nPoints=" + std::to_string(E.size()));

    if (E.size() < 3) {
        Debug::Log(Debug::TRACKER, "FitRes: too few points (<3), skipping fit");
        return;
    }

    // ---- FWHM graph ----
    auto g = new TGraphErrors(E.size(), &E[0], &S[0], nullptr, &dS[0]);
    g->SetName(("FWHMvsEnergy_" + tag).c_str());

    // ---- Fit with physical lower bounds so the curve can't invert ----
    auto f = new TF1(("FWHM_fit" + tag).c_str(), "sqrt([0] + [1]*x)", 0, 3000);
    f->SetParameters(1.0, 0.002);
    f->SetParLimits(0, 1e-4, 1000.0);
    f->SetParLimits(1, 0.0,  10.0);

    g->Fit(f, "R Q B");

    double a = f->GetParameter(0);
    double b = f->GetParameter(1);

    Debug::Log(Debug::TRACKER,
        "FitRes result: a=" + std::to_string(a) + "  b=" + std::to_string(b));

    // ---- Save to ROOT file ----
    RootFileManager::SaveGraph(fout, g);
    RootFileManager::SaveObject(fout, f, "FWHM_fit_par");

    model.UpdateFromQuadratureFit(a, b);
}
