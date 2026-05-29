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
    // Enhanced peak-shape models — each adds parameters per Gaussian and competes
    // with the plain Gaussian via chi²/ndf.  Both default on for better areas.
    bool   useStepBg   = true;   // Compton-step background under each peak
    bool   useExpoTail = true;   // HYPERMET low-energy exponential tail

    // ML sigma constraint — filled by PeakFitter when useMLSigma is set in BgOptions.
    // One entry per peak in the group (keV); NaN or <=0 means fall back to default bounds.
    std::vector<double> mlSigmaHints;
    double              mlSigmaBandFrac = 0.08;  // ±8% band around ML prediction

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
    // Build — n×gaus  [+ Compton step]  [+ HYPERMET tail]  + pol1
    //
    // Parameter layout (matches GammaFitGUI_shared.h BuildNGaussFormulaEx):
    //   [3i], [3i+1], [3i+2]      A_i, E_i, sig_i        (i = 0..n-1)
    //   [3n],  [3n+1]             bg0, bg1
    //   [3n+2+i]                  step_i                  (comptonStep only)
    //   [3n+2+nStep+i]            T_i  (tail amplitude)   (expoTail only)
    //   [3n+2+nStep+n+i]          beta_i (tail slope)     (expoTail only)
    // -----------------------------------------------------------------------
    static TF1* Build(int nPeaks, double xmin, double xmax,
                      bool comptonStep = false, bool expoTail = false) {
        std::string formula = BuildFormula_(nPeaks, comptonStep, expoTail);
        return new TF1("model", formula.c_str(), xmin, xmax);
    }

    static int NTotalPars(int n, bool step, bool tail) {
        return 3*n + 2 + (step ? n : 0) + (tail ? 2*n : 0);
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
        bool   bestHasStep = false;
        bool   bestHasTail = false;

        // Pre-compute background estimate once for the whole window
        double bgL = SafeBinContent_(h, xmin);
        double bgR = SafeBinContent_(h, xmax);

        // Variants to try: (useStep, useTail).  The plain model is always tried;
        // enhanced models compete via chi²/ndf (extra params reduce ndf naturally).
        std::vector<std::pair<bool,bool>> variants;
        variants.push_back({false, false});
        if (cfg.useStepBg)   variants.push_back({true,  false});
        if (cfg.useExpoTail) variants.push_back({false, true });
        if (cfg.useStepBg && cfg.useExpoTail) variants.push_back({true, true});

        for (int n = 1; n <= maxPeaks; n++) {
            for (auto [useStep, useTail] : variants) {
                TF1* f = Build(n, xmin, xmax, useStep, useTail);
                f->SetName(Form("model_%dpk%s%s", n,
                                useStep ? "_step" : "",
                                useTail ? "_tail" : ""));

                // ---- try to seed from DB ----
                bool seededFromDB = false;
                if (db) seededFromDB = db->SeedFromDB(key, f);

                if (!seededFromDB) {
                    Debug::Log(Debug::FITTER, std::to_string(n) + "-peak"
                               + (useStep ? "+step" : "")
                               + (useTail ? "+tail" : "")
                               + ": fresh seed");
                    InitParams_(f, n, peaks, res, h, bgL, bgR, xmin, xmax, cfg);
                    InitStepTailParams_(f, n, useStep, useTail);
                    RelaxBGForStepTail_(f, n, useStep, useTail);
                } else {
                    Debug::Log(Debug::FITTER, std::to_string(n) + "-peak"
                               + (useStep ? "+step" : "")
                               + (useTail ? "+tail" : "")
                               + ": seeded from DB");
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

                // AIC-style selection — Pearson chi2/ndf; extra parameters in
                // step/tail models reduce ndf, naturally penalising overfit.
                int nModelPars = NTotalPars(n, useStep, useTail);
                ResidualMetrics rmN = FitDatabase::ComputeResiduals(h, f, xmin, xmax);
                double chi2ndf   = (rmN.rms < 1.0e6) ? FitDatabase::Chi2Ndf(rmN, nModelPars)
                                                      : r->Chi2() / r->Ndf();

                double areaRatio    = ComputeAreaRatio_(f, h, 3*n, n, xmin, xmax);
                double ratioPenalty = 0.0;
                if (areaRatio > 0.0) {
                    if      (areaRatio < 0.3 || areaRatio > 3.0) ratioPenalty = 10.0;
                    else if (areaRatio < 0.5 || areaRatio > 2.0) ratioPenalty =  3.0;
                    else if (areaRatio < 0.7 || areaRatio > 1.5) ratioPenalty =  1.0;
                }
                double penalised = chi2ndf + cfg.kPenalty * (n - 1) + ratioPenalty;

                Debug::Log(Debug::FITTER, std::to_string(n) + "-peak"
                           + (useStep ? "+step" : "") + (useTail ? "+tail" : "") +
                           "  chi2/ndf=" + std::to_string(chi2ndf) +
                           "  area_ratio=" + (areaRatio > 0 ? std::to_string(areaRatio) : "n/a") +
                           "  penalised=" + std::to_string(penalised) +
                           "  bestSoFar=" + std::to_string(bestChi2));

                if (penalised < bestChi2) {
                    if (best && best != f) delete best;
                    bestChi2    = penalised;
                    best        = f;
                    bestFit     = r;
                    nUsed       = n;
                    bestHasStep = useStep;
                    bestHasTail = useTail;
                    Debug::Log(Debug::FITTER, "  -> new best: " + std::to_string(n) + "-peak"
                               + (useStep ? "+step" : "") + (useTail ? "+tail" : ""));
                } else {
                    delete f;
                }
            }  // variants
        }  // n
 
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
        if (best && bestFit.Get() && bestFit->Ndf() > 0) {
            best = ResidualExpand_(h, best, res, bestFit, bestChi2, nUsed,
                                   xmin, xmax, cfg, bestHasStep, bestHasTail);
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

        // ---- ensure background never goes negative over the fit window ----
        // Skip for step/tail models: they have positive contributions on the left
        // that keep the TOTAL model non-negative even when bg0 < 0.  Applying the
        // clamp post-fit changes MIGRAD's parameters, making the drawn TF1 not
        // match the fit (the visual "wrong background" the user sees).
        if (best && !bestHasStep && !bestHasTail) {
            int npar = best->GetNpar();
            bool isDG = (npar == 7);
            int bgIdx = isDG ? 5 : 3 * nUsed;
            ClampNonNegativeBG_(best, bgIdx, xmin, xmax);
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

            SetParSafe_(f, 3*i+0, A, cfg.ampLoFrac * A, cfg.ampHiFrac * A);
            SetParSafe_(f, 3*i+1, E, E - cfg.meanWindowKeV, E + cfg.meanWindowKeV);
            {
                double sigSeed, sigLo, sigHi;
                if (i < (int)cfg.mlSigmaHints.size() &&
                    !std::isnan(cfg.mlSigmaHints[i]) && cfg.mlSigmaHints[i] > 0.0) {
                    sigSeed = cfg.mlSigmaHints[i];
                    sigLo   = sigSeed * (1.0 - cfg.mlSigmaBandFrac);
                    sigHi   = sigSeed * (1.0 + cfg.mlSigmaBandFrac);
                } else {
                    sigSeed = sigma;
                    sigLo   = cfg.sigmaLoFrac * sigma;
                    sigHi   = cfg.sigmaHiFrac * sigma;
                }
                SetParSafe_(f, 3*i+2, sigSeed, sigLo, sigHi);
            }
        }

        double bgMag  = std::max({std::abs(bgL), std::abs(bgR), 10.0});
        double slope  = (xrange > 0) ? (bgR - bgL) / xrange : 0.0;
        SetParSafe_(f, 3*n,   bgL,   0.0,   bgMag * 10);
        SetParSafe_(f, 3*n+1, slope, -1e3,  1e3);
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
            SetParSafe_(f, 3*i+0, A,    cfg.ampLoFrac * A,     cfg.ampHiFrac * A);
            SetParSafe_(f, 3*i+1, mean, E - cfg.meanWindowKeV, E + cfg.meanWindowKeV);
            {
                double sigLo, sigHi;
                if (i < (int)cfg.mlSigmaHints.size() &&
                    !std::isnan(cfg.mlSigmaHints[i]) && cfg.mlSigmaHints[i] > 0.0) {
                    sigLo = cfg.mlSigmaHints[i] * (1.0 - cfg.mlSigmaBandFrac);
                    sigHi = cfg.mlSigmaHints[i] * (1.0 + cfg.mlSigmaBandFrac);
                } else {
                    sigLo = cfg.sigmaLoFrac * sigma;
                    sigHi = cfg.sigmaHiFrac * sigma;
                }
                SetParSafe_(f, 3*i+2, sig, sigLo, sigHi);
            }
        }
        double bg0    = f->GetParameter(3*n);
        double bg1    = f->GetParameter(3*n+1);
        double bgMag2 = std::max(std::abs(bg0), 10.0);
        SetParSafe_(f, 3*n,   bg0, 0.0,   bgMag2 * 10);
        SetParSafe_(f, 3*n+1, bg1, -1e3,  1e3);
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
        // Try all model variants in order of increasing complexity.
        struct Variant { bool step; bool tail; const char* tag; };
        static constexpr Variant kVariants[] = {
            {false, false, ""},
            {true,  false, "_step"},
            {false, true,  "_tail"},
            {true,  true,  "_step_tail"},
        };
        for (const auto& v : kVariants) {
            for (int n = 1; n <= 3; n++) {
                TF1* f = Build(n, xmin, xmax, v.step, v.tail);
                f->SetName(Form("model_cached_%dpk%s", n, v.tag));
                if (db->ForceSeedFromDB(key, f)) {
                    nUsed = n;
                    return f;
                }
                delete f;
            }
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
    // ClampNonNegativeBG_ — raise bg0 so the linear background stays >= 0
    // over the entire fit window.  The slope is unchanged; only the constant
    // term shifts up by the amount needed to bring the minimum to zero.
    // Called on the winning TF1 after MIGRAD, before storing in the cache.
    // npar_bg: index of bg0 in the parameter list (= 3*n for n-Gaussian model,
    //          or 5 for the double-Gaussian model).
    // -----------------------------------------------------------------------
    static void ClampNonNegativeBG_(TF1* f, int npar_bg, double xmin, double xmax) {
        double bg0   = f->GetParameter(npar_bg);
        double bg1   = f->GetParameter(npar_bg + 1);
        double bgMin = std::min(bg0 + bg1 * xmin, bg0 + bg1 * xmax);
        if (bgMin < 0.0) {
            f->SetParameter(npar_bg, bg0 - bgMin);  // bgMin is negative, so this raises bg0
        }
    }

    // -----------------------------------------------------------------------
    // ComputeAreaRatio_ — integral-of-fit / observed-counts-above-background.
    // Returns -1 if the denominator is too small to be reliable.
    // bgIdx: parameter index of bg0 (3*n for standard model, 5 for DG).
    // nGauss: number of Gaussian peaks in the model.
    // -----------------------------------------------------------------------
    static double ComputeAreaRatio_(TF1* f, TH1* h,
                                    int bgIdx, int nGauss,
                                    double xmin, double xmax)
    {
        double bg0 = f->GetParameter(bgIdx);
        double bg1 = f->GetParameter(bgIdx + 1);

        // Observed counts above the fitted linear background
        double obsAboveBg = 0.0;
        int b1 = h->FindBin(xmin), b2 = h->FindBin(xmax);
        for (int b = b1; b <= b2; b++)
            obsAboveBg += h->GetBinContent(b) - (bg0 + bg1 * h->GetBinCenter(b));

        if (obsAboveBg < 1.0) return -1.0;

        // Total fitted Gaussian area (counts)
        double fittedArea = 0.0;
        for (int i = 0; i < nGauss; i++) {
            double A   = f->GetParameter(3*i);
            double E   = f->GetParameter(3*i + 1);
            double sig = f->GetParameter(3*i + 2);
            double bw  = h->GetBinWidth(h->FindBin(E));
            if (bw <= 0.0) bw = 1.0;
            fittedArea += A * sig * std::sqrt(2.0 * M_PI) / bw;
        }

        return fittedArea / obsAboveBg;
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
    // BuildFormula_ — private copy of BuildNGaussFormulaEx (from shared.h)
    // kept here so fitting code has no GUI header dependency.
    // -----------------------------------------------------------------------
    static std::string BuildFormula_(int n, bool step, bool tail) {
        auto idx = [](int v) { return "[" + std::to_string(v) + "]"; };
        std::string f;
        for (int i = 0; i < n; i++) {
            if (i) f += "+";
            int p = 3*i;
            f += idx(p) + "*exp(-0.5*((x-" + idx(p+1) + ")/" + idx(p+2) + ")^2)";
        }
        f += "+" + idx(3*n) + "+" + idx(3*n+1) + "*x";
        if (step) {
            for (int i = 0; i < n; i++) {
                int sIdx = 3*n + 2 + i;
                f += "+" + idx(sIdx) + "*TMath::Erfc((x-"
                   + idx(3*i+1) + ")/(" + idx(3*i+2) + "*1.41421356))";
            }
        }
        if (tail) {
            int nStep = step ? n : 0;
            for (int i = 0; i < n; i++) {
                int tIdx = 3*n + 2 + nStep + i;
                int bIdx = 3*n + 2 + nStep + n + i;
                // HYPERMET: T/2 * exp((x-E)/β + σ²/(2β²)) * erfc((x-E)/(σ√2) + σ/(β√2))
                f += "+" + idx(tIdx) + "/2*exp((x-" + idx(3*i+1) + ")/" + idx(bIdx)
                   + "+" + idx(3*i+2) + "*" + idx(3*i+2)
                   + "/(2*" + idx(bIdx) + "*" + idx(bIdx) + "))"
                   + "*TMath::Erfc((x-" + idx(3*i+1) + ")/(" + idx(3*i+2) + "*1.41421356)"
                   + "+" + idx(3*i+2) + "/(" + idx(bIdx) + "*1.41421356))";
            }
        }
        return f;
    }

    // -----------------------------------------------------------------------
    // InitStepTailParams_ — seed and bound step/tail parameters after the
    // standard Gaussian + background parameters have been initialised.
    // -----------------------------------------------------------------------
    static void InitStepTailParams_(TF1* f, int n, bool step, bool tail) {
        if (!step && !tail) return;
        int nStep = step ? n : 0;
        for (int i = 0; i < n; i++) {
            double A   = std::max(std::abs(f->GetParameter(3*i)), 1.0);
            double sig = std::max(f->GetParameter(3*i + 2), 0.01);
            if (step) {
                int sIdx = 3*n + 2 + i;
                // Step height ~5% of peak amplitude; bounded [0, half peak height].
                SetParSafe_(f, sIdx, 0.05 * A, 0.0, 0.5 * A);
            }
            if (tail) {
                int tIdx = 3*n + 2 + nStep + i;
                int bIdx = 3*n + 2 + nStep + n + i;
                // Tail amplitude ~10% of peak; tail slope ≈ sigma (physical range: 0.1–5σ).
                SetParSafe_(f, tIdx, 0.1 * A,  0.0,         A);
                SetParSafe_(f, bIdx, sig,       0.1 * sig,   5.0 * sig);
            }
        }
    }

    // -----------------------------------------------------------------------
    // RelaxBGForStepTail_ — called immediately after InitStepTailParams_.
    //
    // Each Compton step adds ≈ 2H at the left edge of the fit window (erfc→2)
    // and each HYPERMET tail adds ≈ T near the peak.  If bg0 is constrained
    // ≥ 0 (as InitParams_ sets it), MIGRAD cannot compensate for this positive
    // left-side contribution, forcing the background slope wildly negative or
    // suppressing the step/tail to zero.
    //
    // The fix: lower bg0's seed and lower bound by the expected step+tail
    // contribution, giving MIGRAD the freedom to find the true minimum.
    // ClampNonNegativeBG_ is NOT applied to step/tail models because the
    // post-hoc parameter change would make the drawn TF1 inconsistent with
    // the MIGRAD result (which is what the user sees as "wrong background").
    // -----------------------------------------------------------------------
    static void RelaxBGForStepTail_(TF1* f, int n, bool step, bool tail) {
        if (!step && !tail) return;
        int nStep = step ? n : 0;

        double stepContrib = 0.0;   // erfc → 2 at far left
        double tailContrib  = 0.0;  // conservative: ≤ T near peak
        for (int i = 0; i < n; i++) {
            if (step) stepContrib += 2.0 * f->GetParameter(3*n + 2 + i);
            if (tail) tailContrib += f->GetParameter(3*n + 2 + nStep + i);
        }
        double totalContrib = stepContrib + tailContrib;
        if (totalContrib <= 0.0) return;

        double lo, hi;
        f->GetParLimits(3*n, lo, hi);
        double bg0    = f->GetParameter(3*n);
        double newLo  = -totalContrib * 2.0;          // 2× safety margin
        double newBg0 = bg0 - totalContrib;            // seed corrected
        SetParSafe_(f, 3*n, std::max(newBg0, newLo), newLo, hi);
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
                                 bool            hasStep  = false,
                                 bool            hasTail  = false,
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

            TF1* f = Build(newN, xmin, xmax, hasStep, hasTail);
            f->SetName(Form("model_res_%dpk%s%s", newN,
                            hasStep ? "_step" : "", hasTail ? "_tail" : ""));

            double bgL = SafeBinContent_(h, xmin);
            double bgR = SafeBinContent_(h, xmax);
            InitParams_(f, newN, expanded, res, h, bgL, bgR, xmin, xmax, cfg);
            InitStepTailParams_(f, newN, hasStep, hasTail);
            RelaxBGForStepTail_(f, newN, hasStep, hasTail);

            TFitResultPtr r = TryFit_(h, f, false, cfg.useLogLikelihood);
            if (!r.Get() || r->Ndf() <= 0) { delete f; break; }

            int nModelPars = NTotalPars(newN, hasStep, hasTail);
            ResidualMetrics rmEx = FitDatabase::ComputeResiduals(h, f, xmin, xmax);
            double chi2ndf   = (rmEx.rms < 1.0e6) ? FitDatabase::Chi2Ndf(rmEx, nModelPars)
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
