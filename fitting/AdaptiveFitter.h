#ifndef ADAPTIVE_FITTER_H
#define ADAPTIVE_FITTER_H
 
// ---------------------------------------------------------------------------
// AdaptiveFitter.h  (v3)
//
// Improvements over v2:
//   - DB-seeded initial parameters (uses stored fit as starting point)
//   - StoreIfBetter integration (auto-replaces DB entry when fit improves)
//   - Smarter amplitude seed: uses histogram integral in ±2σ window
//   - Iterative sigma loosening on failure (3-stage retry)
//   - Background seed from linear interpolation across fit window
//   - Peak-height normalisation guard (prevents absurdly tall seeds)
//   - Fit range now ±6σ per edge peak (was 5×FWHM, inconsistent at low E)
//   - Per-peak convergence check: flags peaks whose fitted mean drifted >5 keV
// ---------------------------------------------------------------------------
 
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <cmath>
 
#include "TF1.h"
#include "TH1.h"
#include "TFitResult.h"
#include "TVirtualFitter.h"
#include "Debug.h"
 
#include "FitDatabase.h"

// Defined at file scope so AdaptiveFitter::FitGroup can use Config{} as a
// default argument — a nested struct with in-class initializers cannot be
// used as a default argument inside its own enclosing class (C++ §11.3.6).
// Peaks known to be physically broader than detector resolution alone
// (e.g. 511 keV annihilation peak with Doppler broadening). When the seed
// energy falls within ±tolerance keV of the listed energy, AdaptiveFitter
// tries a shared-mean double-Gaussian in addition to the normal model and
// picks whichever wins the AIC comparison.
struct DoubleGaussOverride {
    double energy;     // keV
    double tolerance;  // ±keV
    double broadFrac;  // sigma_broad = broadFrac × sigma_narrow
};

struct AdaptiveFitterConfig {
    double kPenalty       = 0.05;
    double sigmaLoFrac    = 0.3;
    double sigmaHiFrac    = 3.0;
    double ampLoFrac      = 0.05;
    double ampHiFrac      = 8.0;
    double meanWindowKeV  = 8.0;
    double driftWarnKeV   = 5.0;
    int    maxPeaksPerGrp = 3;
    bool   useLogLikelihood = true;   // "L" option; uncheck for chi2/least-squares
    bool   useImprove       = false;  // "M" option: IMPROVE after MIGRAD

    // Default: treat the 511 keV annihilation peak as double-Gaussian.
    // broadFrac=3 means the broad component is 3× wider than the narrow one.
    std::vector<DoubleGaussOverride> doubleGaussPeaks = {
        {511.0, 5.0, 3.0}
    };
};

// ---------------------------------------------------------------------------
class AdaptiveFitter {
public:
 
    static constexpr double kFWHMtoSigma = 1.0 / 2.3548200450309493;
 
    // -----------------------------------------------------------------------
    // Sigma_ — converts whatever ResModel::FWHM(E) returns into sigma.
    //   All internal code must go through this; never call res.FWHM directly.
    // -----------------------------------------------------------------------
    template<typename ResModel>
    static double Sigma_(const ResModel& res, double E) {
        return res.FWHM(E) * kFWHMtoSigma;
    }
 
    using Config = AdaptiveFitterConfig;
 
    // -----------------------------------------------------------------------
    // Build — composite model: n×gaus + pol1
    // -----------------------------------------------------------------------
    static TF1* Build(int nPeaks, double xmin, double xmax) {
        std::string formula;
        for (int i = 0; i < nPeaks; i++)
            formula += "gaus(" + std::to_string(3*i) + ")+";
        formula += "pol1(" + std::to_string(3*nPeaks) + ")";
        return new TF1("model", formula.c_str(), xmin, xmax);
    }

    // -----------------------------------------------------------------------
    // BuildDoubleGauss — shared-mean double Gaussian + pol1 (7 parameters)
    //   [0] A_narrow   [1] mean (shared)   [2] sigma_narrow
    //   [3] A_broad                         [4] sigma_broad
    //   [5] bg0        [6] bg1
    // Both Gaussians share [1] so they cannot drift to different energies.
    // -----------------------------------------------------------------------
    static TF1* BuildDoubleGauss(double xmin, double xmax) {
        const char* formula =
            "[0]*exp(-0.5*((x-[1])/[2])^2)"
            "+[3]*exp(-0.5*((x-[1])/[4])^2)"
            "+[5]+[6]*x";
        return new TF1("model_dg", formula, xmin, xmax);
    }
 
    // -----------------------------------------------------------------------
    // FitGroup — main entry point.
    //
    // Template ResModel must provide:
    //   double FWHM(double E)  const
    //   double Sigma(double E) const
    //
    // db  — optional pointer to FitDatabase; if non-null, seeds are loaded
    //       from the DB and the result is stored back if it wins.
    // -----------------------------------------------------------------------
    template<typename ResModel>
    static TF1* FitGroup(TH1*                       h,
                         const std::vector<double>&  peaks,
                         const ResModel&             res,
                         TFitResultPtr&              bestFit,
                         int&                        nUsed,
                         FitDatabase*                db     = nullptr,
                         const Config&               cfg    = Config{})
    {
        const int maxPeaks = std::min((int)peaks.size(), cfg.maxPeaksPerGrp);
        const std::string key = FitDatabase::MakeKey(peaks);

        const double sigL  = Sigma_(res, peaks.front());
        const double sigR  = Sigma_(res, peaks.back());
        const double xmin  = peaks.front() - 4.0 * sigL;
        const double xmax  = peaks.back()  + 4.0 * sigR;

        Debug::Log(Debug::FITTER, "FitGroup key=" + key +
                   "  xmin=" + std::to_string(xmin) +
                   "  xmax=" + std::to_string(xmax) +
                   "  maxPeaks=" + std::to_string(maxPeaks));

        // ---- Cache-only mode: skip MIGRAD, return stored parameters directly ----
        if (db && db->cacheOnly) {
            Debug::Log(Debug::FITTER, "Cache-only mode for key=" + key);
            return ReconstructFromCache_(db, key, xmin, xmax, nUsed);
        }

        TF1*   best     = nullptr;
        double bestChi2 = std::numeric_limits<double>::max();

        // Pre-compute background estimate once for the whole window
        double bgL = SafeBinContent_(h, xmin);
        double bgR = SafeBinContent_(h, xmax);

        for (int n = 1; n <= maxPeaks; n++) {
 
            TF1* f = Build(n, xmin, xmax);
            f->SetName(("model_" + std::to_string(n) + "pk").c_str());
 
            // ---- try to seed from DB ----
            bool seededFromDB = false;
            if (db) seededFromDB = db->SeedFromDB(key, f);

            if (!seededFromDB) {
                Debug::Log(Debug::FITTER, std::to_string(n) + "-peak: fresh seed (no DB entry)");
                InitParams_(f, n, peaks, res, h, bgL, bgR, xmin, xmax, cfg);
            } else {
                Debug::Log(Debug::FITTER, std::to_string(n) + "-peak: seeded from DB");
                EnforceLimits_(f, n, peaks, res, cfg);
            }

            // ---- fit: 3-stage retry with progressively looser sigma ----
            TFitResultPtr r = TryFit_(h, f, false, cfg.useLogLikelihood);
            Debug::Log(Debug::FITTER, std::to_string(n) + "-peak stage1"
                       "  status=" + std::to_string(r.Get() ? r->Status() : -1) +
                       "  edm=" + std::to_string(r.Get() ? r->Edm() : -1.0));

            if (!r.Get() || r->Status() != 0) {
                Debug::Log(Debug::FITTER, std::to_string(n) + "-peak stage2: relaxing sigma/amp");
                for (int i = 0; i < n; i++) {
                    double sig = Sigma_(res, peaks[i]);
                    double A   = std::max(std::abs(f->GetParameter(3*i)), 1.0);
                    SetParSafe_(f, 3*i+2, f->GetParameter(3*i+2), 0.1*sig, 8.0*sig);
                    SetParSafe_(f, 3*i+0, A, 0.01*A, 50.0*A);
                }
                r = TryFit_(h, f, /*improve=*/true, cfg.useLogLikelihood);
                Debug::Log(Debug::FITTER, std::to_string(n) + "-peak stage2"
                           "  status=" + std::to_string(r.Get() ? r->Status() : -1));
            }

            if (!r.Get() || r->Status() != 0) {
                Debug::Log(Debug::FITTER, std::to_string(n) + "-peak stage3: freeing sigma");
                for (int i = 0; i < n; i++) {
                    SetParSafe_(f, 3*i+2, f->GetParameter(3*i+2), 0.01, 50.0);
                }
                r = TryFit_(h, f, /*improve=*/true, cfg.useLogLikelihood);
                Debug::Log(Debug::FITTER, std::to_string(n) + "-peak stage3"
                           "  status=" + std::to_string(r.Get() ? r->Status() : -1));
            }

            Debug::LogFitFailure(n,
                r.Get() ? r->Status() : -1,
                r.Get() && r->Ndf() > 0 ? r->Chi2()/r->Ndf() : -1.0);

            if (!r.Get() || r->Ndf() <= 0) { delete f; continue; }

            CheckDrift_(f, n, peaks, cfg.driftWarnKeV);

            // AIC-style selection — Pearson chi2/ndf with proper DOF
            ResidualMetrics rmN = FitDatabase::ComputeResiduals(h, f, xmin, xmax);
            double chi2ndf   = (rmN.rms < 1.0e6) ? FitDatabase::Chi2Ndf(rmN, 3*n + 2)
                                                  : r->Chi2() / r->Ndf();
            double penalised = chi2ndf + cfg.kPenalty * (n - 1);

            Debug::Log(Debug::FITTER, std::to_string(n) + "-peak"
                       "  chi2/ndf=" + std::to_string(chi2ndf) +
                       "  penalised=" + std::to_string(penalised) +
                       "  bestSoFar=" + std::to_string(bestChi2));

            if (penalised < bestChi2) {
                if (best && best != f) delete best;
                bestChi2 = penalised;
                best     = f;
                bestFit  = r;
                nUsed    = n;
                Debug::Log(Debug::FITTER, "  -> new best: " + std::to_string(n) + "-peak model");
            } else {
                delete f;
            }
        }
 
        // ---- double-Gaussian attempt for known broad peaks (e.g. 511 keV) ----
        for (const auto& peak : peaks) {
            const DoubleGaussOverride* ovr = MatchDoubleGauss_(peak, cfg);
            if (!ovr) continue;

            TF1* dg = BuildDoubleGauss(xmin, xmax);
            dg->SetName("model_dg");
            InitDoubleGaussParams_(dg, peak, res, h, bgL, bgR,
                                   xmin, xmax, cfg, *ovr);

            TFitResultPtr r = TryFit_(h, dg, false, cfg.useLogLikelihood);
            if (!r.Get() || r->Ndf() <= 0) { delete dg; break; }

            CheckDrift_DG_(dg, peak, cfg.driftWarnKeV);

            // AIC penalty: 2 extra free parameters vs single-Gaussian+pol1
            ResidualMetrics rmDG = FitDatabase::ComputeResiduals(h, dg, xmin, xmax);
            double chi2ndf   = (rmDG.rms < 1.0e6) ? FitDatabase::Chi2Ndf(rmDG, 7)
                                                   : r->Chi2() / r->Ndf();
            double penalised = chi2ndf + cfg.kPenalty * 2.0;

            Debug::Log(Debug::FITTER, "[FITTER] Double-Gauss 511: chi2/ndf=" +
                       std::to_string(chi2ndf));

            if (penalised < bestChi2) {
                if (best) delete best;
                bestChi2 = penalised;
                best     = dg;
                bestFit  = r;
                nUsed    = 1;
            } else {
                delete dg;
            }
            break;  // only one double-Gauss peak per group
        }

        // ---- Residual-guided peak discovery (post-fit) ----------------------
        // Scan (data − best_fit)/σ for bins that exceed 3σ. If a significant
        // cluster is found, add the missed peak, expand the model, and refit.
        // The AIC penalty prevents spurious expansions.
        if (best && bestFit.Get() && bestFit->Ndf() > 0) {
            best = ResidualExpand_(h, best, res, bestFit, bestChi2, nUsed,
                                   xmin, xmax, cfg);
        }

        // ---- Compute residual metrics for composite cache comparison --------
        ResidualMetrics rm;
        if (best)
            rm = FitDatabase::ComputeResiduals(h, best, xmin, xmax);

        // ---- warn if any parameter hit its bound (errors will be zero) ----
        if (best) {
            int npar = best->GetNpar();
            for (int i = 0; i < npar; i++) {
                double val = best->GetParameter(i);
                double lo, hi;
                best->GetParLimits(i, lo, hi);
                if (hi <= lo) continue;
                double range = hi - lo;
                if (std::abs(val - lo) / range < 0.002)
                    Debug::Log(Debug::FITTER,
                        Form("[BoundHit] %s  par[%d]=%s  hit LOWER bound %.5g",
                             key.c_str(), i, best->GetParName(i), lo));
                else if (std::abs(val - hi) / range < 0.002)
                    Debug::Log(Debug::FITTER,
                        Form("[BoundHit] %s  par[%d]=%s  hit UPPER bound %.5g",
                             key.c_str(), i, best->GetParName(i), hi));
            }
        }

        // ---- store in DB if better ----
        if (best && db) {
            FitEntry candidate = FitDatabase::MakeEntry(key, bestFit, best, nUsed, rm);
            db->StoreIfBetter(key, candidate);
        }

        return best;
    }

private:
 
    // -----------------------------------------------------------------------
    // InitParams_ — physics-aware parameter initialisation
    // -----------------------------------------------------------------------
    template<typename ResModel>
    static void InitParams_(TF1* f, int n,
                            const std::vector<double>& peaks,
                            const ResModel& res,
                            TH1* h,
                            double bgL, double bgR,
                            double xmin, double xmax,
                            const Config& cfg)
    {
        double xrange = xmax - xmin;
 
        for (int i = 0; i < n; i++) {
            double E     = peaks[i];
            double sigma = Sigma_(res, E);

            // Amplitude seed: direct bin height minus interpolated background.
            // ROOT's gaus amplitude parameter IS the peak height, so reading the
            // bin at E is the most direct estimate and can't be thrown off by
            // background integration errors the way the integral method can.
            double bgAtE = LinearInterp_(bgL, bgR, xmin, xmax, E);
            double A     = std::max(h->GetBinContent(h->FindBin(E)) - bgAtE, 1.0);

            Debug::Log(Debug::FITTER, "[FITTER] Peak " + std::to_string(i) +
                       " E=" + std::to_string(E) +
                       " A_seed=" + std::to_string(A) +
                       " sigma_seed=" + std::to_string(sigma));

            SetParSafe_(f, 3*i+0, A,     cfg.ampLoFrac  * A,     cfg.ampHiFrac  * A);
            SetParSafe_(f, 3*i+1, E,     E - cfg.meanWindowKeV,  E + cfg.meanWindowKeV);
            SetParSafe_(f, 3*i+2, sigma, cfg.sigmaLoFrac * sigma, cfg.sigmaHiFrac * sigma);
        }

        // Background: seed from edge interpolation; allow negative to handle
        // over-subtracted residuals from TSpectrum background removal.
        double bgMag  = std::max({std::abs(bgL), std::abs(bgR), 10.0});
        double slope  = (xrange > 0) ? (bgR - bgL) / xrange : 0.0;
        SetParSafe_(f, 3*n,   bgL,   -bgMag * 5,  bgMag * 10);
        SetParSafe_(f, 3*n+1, slope, -1e3,         1e3);
    }
 
    // -----------------------------------------------------------------------
    // EnforceLimits_ — apply limits to a DB-seeded TF1 (don't change values)
    // -----------------------------------------------------------------------
    template<typename ResModel>
    static void EnforceLimits_(TF1* f, int n,
                               const std::vector<double>& peaks,
                               const ResModel& res,
                               const Config& cfg)
    {
        for (int i = 0; i < n; i++) {
            double E     = peaks[i];
            double sigma = Sigma_(res, E);
            double A     = std::max(std::abs(f->GetParameter(3*i)), 1.0);
            double mean  = f->GetParameter(3*i+1);
            double sig   = f->GetParameter(3*i+2);

            // Clamp cached values into their new limits via SetParSafe_
            SetParSafe_(f, 3*i+0, A,    cfg.ampLoFrac   * A,     cfg.ampHiFrac   * A);
            SetParSafe_(f, 3*i+1, mean, E - cfg.meanWindowKeV,   E + cfg.meanWindowKeV);
            SetParSafe_(f, 3*i+2, sig,  cfg.sigmaLoFrac * sigma,  cfg.sigmaHiFrac * sigma);
        }
        double bg0    = f->GetParameter(3*n);
        double bg1    = f->GetParameter(3*n+1);
        double bgMag2 = std::max(std::abs(bg0), 10.0);
        SetParSafe_(f, 3*n,   bg0, -bgMag2 * 5,  bgMag2 * 10);
        SetParSafe_(f, 3*n+1, bg1, -1e3,          1e3);
    }

    // -----------------------------------------------------------------------
    // TryFit_ — single fit attempt
    //   "R"  use TF1 range
    //   "S"  return TFitResultPtr
    //   "Q"  quiet
    //   "L"  log-likelihood (better for low counts)
    //   "B"  respect parameter limits (critical!)
    //   "M"  IMPROVE: search for better minimum after MIGRAD — only on
    //        retry stages; using it on every attempt causes mnimpr to report
    //        status 4 even when MIGRAD itself converged cleanly.
    // -----------------------------------------------------------------------
    static TFitResultPtr TryFit_(TH1* h, TF1* f, bool improve = false,
                                  bool logLik = true) {
        std::string opts = "R S Q B";
        if (logLik)  opts += " L";
        if (improve) opts += " M";
        return h->Fit(f, opts.c_str());
    }
 
    // -----------------------------------------------------------------------
    // ReconstructFromCache_ — cache-only mode.
    // Tries each standard n-peak model (n=1..3) then double-Gauss.
    // Returns the first TF1 whose parameter count matches the cache entry,
    // with parameters loaded from the cache.  Returns nullptr if no match.
    // nUsed is set to the number of Gaussian peaks in the returned model.
    // bestFit is intentionally left empty — no MIGRAD was run.
    // -----------------------------------------------------------------------
    static TF1* ReconstructFromCache_(FitDatabase*       db,
                                      const std::string& key,
                                      double             xmin,
                                      double             xmax,
                                      int&               nUsed)
    {
        // Standard models: npar = 3*n + 2
        for (int n = 1; n <= 3; n++) {
            TF1* f = Build(n, xmin, xmax);
            f->SetName(("model_cached_" + std::to_string(n) + "pk").c_str());
            if (db->ForceSeedFromDB(key, f)) {
                nUsed = n;
                return f;
            }
            delete f;
        }
        // Double-Gaussian model: 7 parameters
        TF1* dg = BuildDoubleGauss(xmin, xmax);
        dg->SetName("model_cached_dg");
        if (db->ForceSeedFromDB(key, dg)) {
            nUsed = 1;
            return dg;
        }
        delete dg;

        Debug::Log(Debug::FITTER, "[FITTER] Cache-only: no entry for " + key + " — skipping");
        return nullptr;
    }

    // -----------------------------------------------------------------------
    // SetParSafe_ — set value AND limits together, clamping value into [lo,hi].
    // Eliminates the ROOT Info "lower/upper bounds outside current parameter
    // value" warning, which fires whenever SetParLimits is called with the
    // current parameter value already outside the requested range.
    // -----------------------------------------------------------------------
    static void SetParSafe_(TF1* f, int par, double val, double lo, double hi) {
        if (lo > hi) std::swap(lo, hi);
        double v = std::max(lo, std::min(hi, val));
        f->SetParameter(par, v);
        f->SetParLimits(par, lo, hi);
    }

    // -----------------------------------------------------------------------
    // Helpers
    // -----------------------------------------------------------------------
    static double SafeBinContent_(TH1* h, double x) {
        int bin = h->FindBin(x);
        bin = std::max(1, std::min(bin, h->GetNbinsX()));
        return h->GetBinContent(bin);  // allow negative: bg-subtracted spectra can have negative residuals
    }
 
    static double LinearInterp_(double yL, double yR,
                                 double xL, double xR, double x) {
        if (xR <= xL) return yL;
        return yL + (yR - yL) * (x - xL) / (xR - xL);
    }
 
    static void CheckDrift_(TF1* f, int n,
                             const std::vector<double>& peaks,
                             double warnKeV)
    {
        for (int i = 0; i < n; i++) {
            double fitted  = f->GetParameter(3*i+1);
            double nominal = peaks[i];
            if (std::abs(fitted - nominal) > warnKeV)
                std::cerr << "[AdaptiveFitter] Peak " << i
                          << " drifted: nominal=" << nominal
                          << " fitted=" << fitted << " keV\n";
        }
    }

    // -----------------------------------------------------------------------
    // Double-Gaussian helpers
    // -----------------------------------------------------------------------

    static const DoubleGaussOverride* MatchDoubleGauss_(
            double E, const Config& cfg)
    {
        for (const auto& ovr : cfg.doubleGaussPeaks)
            if (std::abs(E - ovr.energy) <= ovr.tolerance)
                return &ovr;
        return nullptr;
    }

    template<typename ResModel>
    static void InitDoubleGaussParams_(TF1* f,
                                       double E,
                                       const ResModel& res,
                                       TH1* h,
                                       double bgL, double bgR,
                                       double xmin, double xmax,
                                       const Config& cfg,
                                       const DoubleGaussOverride& ovr)
    {
        double sigN  = Sigma_(res, E);                // narrow sigma
        double sigB  = sigN * ovr.broadFrac;          // broad sigma
        double bgAtE = LinearInterp_(bgL, bgR, xmin, xmax, E);
        double total = std::max(h->GetBinContent(h->FindBin(E)) - bgAtE, 1.0);

        // Split amplitude: narrow core gets ~70%, broad skirt ~30%
        double AN = 0.7 * total;
        double AB = 0.3 * total;

        double bgMag = std::max({std::abs(bgL), std::abs(bgR), 10.0});
        double slope = (xmax > xmin) ? (bgR - bgL) / (xmax - xmin) : 0.0;

        SetParSafe_(f, 0, AN,    cfg.ampLoFrac * AN,    cfg.ampHiFrac * AN);
        SetParSafe_(f, 1, E,     E - cfg.meanWindowKeV, E + cfg.meanWindowKeV);
        SetParSafe_(f, 2, sigN,  cfg.sigmaLoFrac * sigN, cfg.sigmaHiFrac * sigN);
        SetParSafe_(f, 3, AB,    cfg.ampLoFrac * AB,    cfg.ampHiFrac * AB);
        // broad sigma must stay wider than narrow; upper cap = 10× narrow
        SetParSafe_(f, 4, sigB,  sigN, 10.0 * sigN);
        SetParSafe_(f, 5, bgL,   -bgMag * 5, bgMag * 10);
        SetParSafe_(f, 6, slope, -1e3, 1e3);
    }

    static void CheckDrift_DG_(TF1* f, double nominal, double warnKeV) {
        double fitted = f->GetParameter(1);
        if (std::abs(fitted - nominal) > warnKeV)
            std::cerr << "[AdaptiveFitter] Double-Gauss mean drifted: nominal="
                      << nominal << " fitted=" << fitted << " keV\n";
    }

    // -----------------------------------------------------------------------
    // ScanResidualPeaks_ — find clusters of bins where (data−fit)/σ > threshold.
    // Returns the energy of the highest-pull bin in each cluster.
    // -----------------------------------------------------------------------
    static std::vector<double> ScanResidualPeaks_(TH1* h, TF1* f,
                                                   double xmin, double xmax,
                                                   double threshold)
    {
        std::vector<double> peaks;
        int binLo = h->FindBin(xmin);
        int binHi = h->FindBin(xmax);

        bool   inCluster = false;
        double bestPull  = 0.0;
        double bestX     = -1.0;

        for (int b = binLo; b <= binHi; b++) {
            double err = h->GetBinError(b);
            if (err <= 0.0) continue;
            double x    = h->GetBinCenter(b);
            double pull = (h->GetBinContent(b) - f->Eval(x)) / err;

            if (pull > threshold) {
                if (!inCluster || pull > bestPull) {
                    bestPull = pull;
                    bestX    = x;
                }
                inCluster = true;
            } else if (inCluster) {
                peaks.push_back(bestX);
                inCluster = false;
                bestPull  = 0.0;
                bestX     = -1.0;
            }
        }
        if (inCluster && bestX > 0.0) peaks.push_back(bestX);
        return peaks;
    }

    // -----------------------------------------------------------------------
    // ResidualExpand_ — iteratively add missed peaks found in residuals.
    // Up to maxPasses rounds; each round costs one MIGRAD call.
    // AIC penalty guards against over-expansion.
    // -----------------------------------------------------------------------
    template<typename ResModel>
    static TF1* ResidualExpand_(TH1*           h,
                                 TF1*           current,
                                 const ResModel& res,
                                 TFitResultPtr&  bestFit,
                                 double&         bestScore,
                                 int&            nUsed,
                                 double          xmin,
                                 double          xmax,
                                 const Config&   cfg,
                                 int             maxPasses = 2)
    {
        for (int pass = 0; pass < maxPasses; pass++) {
            std::vector<double> newPeaks =
                ScanResidualPeaks_(h, current, xmin, xmax, 3.0);
            if (newPeaks.empty()) break;

            // Collect currently fitted peak means
            std::vector<double> expanded;
            for (int i = 0; i < nUsed; i++)
                expanded.push_back(current->GetParameter(3*i + 1));
            for (double E : newPeaks)
                expanded.push_back(E);
            std::sort(expanded.begin(), expanded.end());

            // Merge peaks closer than 2 keV (avoid duplicating existing peaks)
            auto last = std::unique(expanded.begin(), expanded.end(),
                [](double a, double b){ return std::abs(a - b) < 2.0; });
            expanded.erase(last, expanded.end());

            int newN = (int)expanded.size();
            if (newN == nUsed) break;                       // nothing changed
            if (newN > cfg.maxPeaksPerGrp + maxPasses) break; // cap total

            TF1* f = Build(newN, xmin, xmax);
            f->SetName(("model_res_" + std::to_string(newN) + "pk").c_str());

            double bgL = SafeBinContent_(h, xmin);
            double bgR = SafeBinContent_(h, xmax);
            InitParams_(f, newN, expanded, res, h, bgL, bgR, xmin, xmax, cfg);

            TFitResultPtr r = TryFit_(h, f, false, cfg.useLogLikelihood);
            if (!r.Get() || r->Ndf() <= 0) { delete f; break; }

            ResidualMetrics rmEx = FitDatabase::ComputeResiduals(h, f, xmin, xmax);
            double chi2ndf   = (rmEx.rms < 1.0e6) ? FitDatabase::Chi2Ndf(rmEx, 3*newN + 2)
                                                   : r->Chi2() / r->Ndf();
            double penalised = chi2ndf + cfg.kPenalty * (newN - 1);

            Debug::Log(Debug::FITTER,
                "[ResidualExpand] pass=" + std::to_string(pass) +
                "  newN=" + std::to_string(newN) +
                "  chi2/ndf=" + std::to_string(chi2ndf) +
                "  penalised=" + std::to_string(penalised) +
                "  prev=" + std::to_string(bestScore));

            if (penalised < bestScore) {
                delete current;
                current   = f;
                bestFit   = r;
                bestScore = penalised;
                nUsed     = newN;
                Debug::Log(Debug::FITTER,
                    "[ResidualExpand] -> accepted " + std::to_string(newN) + "-peak model");
            } else {
                delete f;
                break;
            }
        }
        return current;
    }
};
 
#endif // ADAPTIVE_FITTER_H
