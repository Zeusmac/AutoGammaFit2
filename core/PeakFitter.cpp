#include "PeakFitter.h"
#include "PeakGrouper.h"
#include "AdaptiveFitter.h"
#include "PeakTracker.h"

#include "Debug.h"
#include "TSpectrum.h"
#include "TROOT.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFitResult.h"
#include "RootFileManager.h"

#ifdef HAS_ONNX
#include "MLInference.h"
#include <memory>

// File-scoped model instances loaded lazily on first useML call.
static std::unique_ptr<OrtMLModel> s_bgModel;
static std::unique_ptr<OrtMLModel> s_peakModel;
static std::unique_ptr<OrtMLModel> s_centroidModel;
static std::unique_ptr<OrtMLModel> s_sigmaModel;
static std::string                  s_bgModelPath;
static std::string                  s_peakModelPath;
static std::string                  s_centroidModelPath;
static std::string                  s_sigmaModelPath;

void PeakFitter::ReloadMLModels()
{
    s_bgModel.reset();
    s_peakModel.reset();
    s_centroidModel.reset();
    s_sigmaModel.reset();
    s_bgModelPath.clear();
    s_peakModelPath.clear();
    s_centroidModelPath.clear();
    s_sigmaModelPath.clear();
}

TH1* PeakFitter::GetMLBackground(TH1* raw, const std::string& modelPath)
{
    if (!raw) return nullptr;
    if (!s_bgModel || s_bgModelPath != modelPath) {
        try {
            s_bgModel     = std::make_unique<OrtMLModel>(modelPath);
            s_bgModelPath = modelPath;
        } catch (const std::exception& ex) {
            Debug::Log(Debug::PEAKFITTER,
                std::string("GetMLBackground: model load failed: ") + ex.what());
            return nullptr;
        }
    }
    const int nBins = raw->GetNbinsX();
    std::vector<float> rawV(nBins, 0.0f);
    float rawMax = 0.0f;
    for (int b = 1; b <= nBins; ++b) {
        float v  = static_cast<float>(std::max(0.0, raw->GetBinContent(b)));
        rawV[b-1] = v;
        rawMax    = std::max(rawMax, v);
    }
    if (rawMax == 0.0f) rawMax = 1.0f;
    for (float& v : rawV) v /= rawMax;

    TH1* bg = static_cast<TH1*>(raw->Clone(Form("%s_mlbg_gui", raw->GetName())));
    bg->SetDirectory(nullptr);
    bg->Reset();
    constexpr int BG_WIN  = 101;
    constexpr int BG_HALF = BG_WIN / 2;
    std::vector<float> win(BG_WIN, 0.0f);
    for (int b = 0; b < nBins; ++b) {
        for (int k = 0; k < BG_WIN; ++k) {
            int idx = b - BG_HALF + k;
            win[k]  = (idx >= 0 && idx < nBins) ? rawV[idx] : 0.0f;
        }
        float bgEst = s_bgModel->RunScalar(win) * rawMax;
        bg->SetBinContent(b + 1, std::max(0.0f, bgEst));
    }
    return bg;
}
TH1* PeakFitter::GetMLPeakProbabilities(TH1* bgSub, const std::string& modelPath)
{
    if (!bgSub) return nullptr;
    if (!s_peakModel || s_peakModelPath != modelPath) {
        try {
            s_peakModel     = std::make_unique<OrtMLModel>(modelPath);
            s_peakModelPath = modelPath;
        } catch (const std::exception& ex) {
            Debug::Log(Debug::PEAKFITTER,
                std::string("GetMLPeakProbabilities: model load failed: ") + ex.what());
            return nullptr;
        }
    }
    const int nBins = bgSub->GetNbinsX();
    std::vector<float> sub(nBins, 0.0f);
    float subMax = 0.0f;
    for (int b = 1; b <= nBins; ++b) {
        float v  = static_cast<float>(std::max(0.0, bgSub->GetBinContent(b)));
        sub[b-1] = v;
        subMax   = std::max(subMax, v);
    }
    if (subMax == 0.0f) subMax = 1.0f;
    for (float& v : sub) v /= subMax;

    TH1* prob = static_cast<TH1*>(bgSub->Clone(
                    Form("%s_mlpkprob_gui", bgSub->GetName())));
    prob->SetDirectory(nullptr);
    prob->Reset();
    constexpr int PK_WIN  = 51;
    constexpr int PK_HALF = PK_WIN / 2;
    std::vector<float> win(PK_WIN, 0.0f);
    for (int b = 0; b < nBins; ++b) {
        for (int k = 0; k < PK_WIN; ++k) {
            int idx = b - PK_HALF + k;
            win[k]  = (idx >= 0 && idx < nBins) ? sub[idx] : 0.0f;
        }
        prob->SetBinContent(b + 1, s_peakModel->RunScalar(win));
    }
    return prob;
}
double PeakFitter::GetMLCentroidOffset(TH1* h, double E, double sig,
                                        const std::string& modelPath)
{
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    if (!s_centroidModel || s_centroidModelPath != modelPath) {
        try {
            s_centroidModel     = std::make_unique<OrtMLModel>(modelPath);
            s_centroidModelPath = modelPath;
        } catch (const std::exception& ex) {
            Debug::Log(Debug::PEAKFITTER,
                std::string("GetMLCentroidOffset: model load failed: ") + ex.what());
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    const int nBins = h->GetNbinsX();
    int peakBin = h->FindBin(E);

    constexpr int CEN_WIN  = 51;
    constexpr int CEN_HALF = CEN_WIN / 2;

    // Extract normalised 51-bin window centred on the peak bin
    std::vector<float> win(CEN_WIN, 0.0f);
    float winMax = 0.0f;
    for (int k = 0; k < CEN_WIN; ++k) {
        int idx = peakBin - CEN_HALF + k;
        float v = (idx >= 1 && idx <= nBins)
                  ? static_cast<float>(std::max(0.0, h->GetBinContent(idx)))
                  : 0.0f;
        win[k]  = v;
        winMax  = std::max(winMax, v);
    }
    if (winMax == 0.0f) return std::numeric_limits<double>::quiet_NaN();
    for (float& v : win) v /= winMax;

    // Model outputs a value in [-0.5, 0.5] bins; convert to keV
    float offsetBins = s_centroidModel->RunScalar(win);
    double bw = h->GetBinWidth(peakBin);
    return static_cast<double>(offsetBins) * bw;
}

double PeakFitter::GetMLSigmaWidth(TH1* h, double E, const std::string& modelPath)
{
    if (!h) return std::numeric_limits<double>::quiet_NaN();
    if (!s_sigmaModel || s_sigmaModelPath != modelPath) {
        try {
            s_sigmaModel     = std::make_unique<OrtMLModel>(modelPath);
            s_sigmaModelPath = modelPath;
        } catch (const std::exception& ex) {
            Debug::Log(Debug::PEAKFITTER,
                std::string("GetMLSigmaWidth: model load failed: ") + ex.what());
            return std::numeric_limits<double>::quiet_NaN();
        }
    }
    const int nBins  = h->GetNbinsX();
    int peakBin = h->FindBin(E);

    constexpr int CEN_WIN  = 51;
    constexpr int CEN_HALF = CEN_WIN / 2;

    std::vector<float> win(CEN_WIN, 0.0f);
    float winMax = 0.0f;
    for (int k = 0; k < CEN_WIN; ++k) {
        int idx = peakBin - CEN_HALF + k;
        float v = (idx >= 1 && idx <= nBins)
                  ? static_cast<float>(std::max(0.0, h->GetBinContent(idx)))
                  : 0.0f;
        win[k]  = v;
        winMax  = std::max(winMax, v);
    }
    if (winMax == 0.0f) return std::numeric_limits<double>::quiet_NaN();
    for (float& v : win) v /= winMax;

    float sigmaBins = s_sigmaModel->RunScalar(win);   // sigma in bins (Softplus → positive)
    double bw = h->GetBinWidth(peakBin);
    return static_cast<double>(sigmaBins) * bw;        // sigma in keV
}

#else
void PeakFitter::ReloadMLModels() {}   // no-op when ONNX not compiled in
TH1* PeakFitter::GetMLBackground(TH1*, const std::string&)               { return nullptr; }
TH1* PeakFitter::GetMLPeakProbabilities(TH1*, const std::string&)         { return nullptr; }
double PeakFitter::GetMLCentroidOffset(TH1*, double, double, const std::string&)
    { return std::numeric_limits<double>::quiet_NaN(); }
double PeakFitter::GetMLSigmaWidth(TH1*, double, const std::string&)
    { return std::numeric_limits<double>::quiet_NaN(); }
#endif

#include <fstream>
#include <algorithm>
#include <limits>
#include <set>
#include "TMatrixDSym.h"

using namespace std;

PeakFitter::PeakFitter(GammaDB& db,
                       PeakTracker* tracker,
                       ResolutionModel& res,
                       FitStorage& storage,
                       FitDatabase* fitdb)
    : db(db), tracker(tracker), res(res), storage(storage), fitdb(fitdb)
{}

void PeakFitter::FitHistogram(TH1* h,
                              TFile* fout,
                              bool enableTracking,
                              TCanvas* extCanvas,
                              BgOptions bg,
                              const std::vector<double>& forcedSeeds)
{
    if (!h) return;
    std::string hname = h->GetName();
    Debug::Log(Debug::PEAKFITTER, "=== FitHistogram: " + hname + " ===");

    // -------------------------------------------------
    // Canvas + text
    // -------------------------------------------------
    TCanvas* c;
    if (extCanvas) {
        c = extCanvas;
        c->SetName(Form("c_%s", hname.c_str()));
        c->SetTitle(hname.c_str());
        c->Clear();
    } else {
        // No external canvas — create one silently in batch mode so no
        // X11 window appears.  Batch state is restored immediately after.
        Bool_t wasBatch = gROOT->IsBatch();
        gROOT->SetBatch(kTRUE);
        c = new TCanvas(Form("c_%s", hname.c_str()),
                        hname.c_str(),
                        1200, 800);
        gROOT->SetBatch(wasBatch);
    }
    c->cd();

    TLatex latex;
    latex.SetTextSize(0.022);

    std::ofstream out("../Gamma_fits/" +
                      hname +
                      "_fit.txt");

    // -------------------------------------------------
    // Working histogram
    // -------------------------------------------------
    TH1* h_work =
        (TH1*)h->Clone(Form("%s_work", hname.c_str()));

    h_work->Sumw2();

    // -------------------------------------------------
    // Background subtraction + peak finding
    // -------------------------------------------------
    std::vector<double> peaks;

    // TSpectrum path — original behaviour, also used as fallback when HAS_ONNX
    // is not compiled in or model files are missing.
    auto runTSpectrum = [&]() {
        TSpectrum spectrum(1000);
        if (bg.subtractBg) {
            TH1* bkg = spectrum.Background(h_work, bg.iterations);
            h_work->Add(bkg, -1);
        }
        if (forcedSeeds.empty()) {
            int nPeaks = spectrum.Search(h_work, bg.tspecSigma, "", bg.tspecThresh);
            Debug::Log(Debug::PEAKFITTER,
                "TSpectrum found " + std::to_string(nPeaks) + " peaks in " + hname +
                "  sigma=" + std::to_string(bg.tspecSigma) +
                "  thresh=" + std::to_string(bg.tspecThresh));
            peaks.assign(spectrum.GetPositionX(), spectrum.GetPositionX() + nPeaks);
        } else {
            peaks = forcedSeeds;
            Debug::Log(Debug::PEAKFITTER,
                "Source mode: using " + std::to_string(peaks.size()) +
                " forced seed energies in " + hname);
        }
    };

    TH1* bkgML = nullptr;  // ML background histogram; kept alive for net area (Method 4)
    bool mlUsed = false;
    if (bg.useML) {
#ifdef HAS_ONNX
        // ── Lazy-load BG model ────────────────────────────────────────────
        bool modelsOk = true;
        if (!s_bgModel || s_bgModelPath != bg.mlBgModelPath) {
            try {
                s_bgModel     = std::make_unique<OrtMLModel>(bg.mlBgModelPath);
                s_bgModelPath = bg.mlBgModelPath;
            } catch (const std::exception& ex) {
                Debug::Log(Debug::PEAKFITTER,
                    std::string("ML BG model load failed (") +
                    bg.mlBgModelPath + "): " + ex.what());
                modelsOk = false;
            }
        }
        // ── Lazy-load peak model ──────────────────────────────────────────
        if (modelsOk && (!s_peakModel || s_peakModelPath != bg.mlPeakModelPath)) {
            try {
                s_peakModel     = std::make_unique<OrtMLModel>(bg.mlPeakModelPath);
                s_peakModelPath = bg.mlPeakModelPath;
            } catch (const std::exception& ex) {
                Debug::Log(Debug::PEAKFITTER,
                    std::string("ML peak model load failed (") +
                    bg.mlPeakModelPath + "): " + ex.what());
                modelsOk = false;
            }
        }

        if (modelsOk) {
            const int nBins = h_work->GetNbinsX();

            // ── 1. Normalize raw spectrum to [0,1] ────────────────────────
            std::vector<float> raw(nBins, 0.0f);
            float rawMax = 0.0f;
            for (int b = 1; b <= nBins; ++b) {
                float v = static_cast<float>(std::max(0.0, h_work->GetBinContent(b)));
                raw[b - 1] = v;
                rawMax     = std::max(rawMax, v);
            }
            if (rawMax == 0.0f) rawMax = 1.0f;
            for (float& v : raw) v /= rawMax;

            // ── 2. Slide BG MLP — 101-bin window, zero-pad at edges ───────
            constexpr int BG_WIN  = 101;
            constexpr int BG_HALF = BG_WIN / 2;
            bkgML = static_cast<TH1*>(h_work->Clone(
                        Form("%s_bkg_ml", hname.c_str())));
            bkgML->SetDirectory(nullptr);
            bkgML->Reset();
            std::vector<float> win(BG_WIN, 0.0f);
            for (int b = 0; b < nBins; ++b) {
                for (int k = 0; k < BG_WIN; ++k) {
                    int idx = b - BG_HALF + k;
                    win[k]  = (idx >= 0 && idx < nBins) ? raw[idx] : 0.0f;
                }
                float bgEst = s_bgModel->RunScalar(win) * rawMax;
                bkgML->SetBinContent(b + 1, std::max(0.0f, bgEst));
            }

            // ── 3. Subtract background, floor negatives ───────────────────
            h_work->Add(bkgML, -1.0);
            // bkgML kept alive — deleted at function end for net area (Method 4)
            for (int b = 1; b <= nBins; ++b)
                if (h_work->GetBinContent(b) < 0.0) h_work->SetBinContent(b, 0.0);

            if (!forcedSeeds.empty()) {
                peaks  = forcedSeeds;
                mlUsed = true;
                Debug::Log(Debug::PEAKFITTER,
                    "ML BG done; using " + std::to_string(peaks.size()) +
                    " forced seed energies in " + hname);
            } else {
                // ── 4. Normalize BG-subtracted spectrum ───────────────────
                std::vector<float> sub(nBins, 0.0f);
                float subMax = 0.0f;
                for (int b = 1; b <= nBins; ++b) {
                    float v  = static_cast<float>(std::max(0.0, h_work->GetBinContent(b)));
                    sub[b-1] = v;
                    subMax   = std::max(subMax, v);
                }
                if (subMax == 0.0f) subMax = 1.0f;
                for (float& v : sub) v /= subMax;

                // ── 5. Slide peak MLP — 51-bin window → probability array ─
                constexpr int PK_WIN  = 51;
                constexpr int PK_HALF = PK_WIN / 2;
                std::vector<float> prob(nBins, 0.0f);
                std::vector<float> pkw(PK_WIN, 0.0f);
                for (int b = 0; b < nBins; ++b) {
                    for (int k = 0; k < PK_WIN; ++k) {
                        int idx = b - PK_HALF + k;
                        pkw[k]  = (idx >= 0 && idx < nBins) ? sub[idx] : 0.0f;
                    }
                    prob[b] = s_peakModel->RunScalar(pkw);
                }

                // ── 6. Non-maximum suppression ─────────────────────────────
                const float thresh = bg.mlPeakThresh;
                for (int b = 0; b < nBins; ++b) {
                    if (prob[b] < thresh) continue;
                    bool isMax = true;
                    for (int k = -PK_HALF; k <= PK_HALF && isMax; ++k) {
                        int idx = b + k;
                        if (idx >= 0 && idx < nBins && idx != b && prob[idx] > prob[b])
                            isMax = false;
                    }
                    if (isMax)
                        peaks.push_back(h_work->GetBinCenter(b + 1));
                }
                mlUsed = true;
                Debug::Log(Debug::PEAKFITTER,
                    "ML found " + std::to_string(peaks.size()) +
                    " peaks in " + hname);
            }
        }
#else
        Debug::Log(Debug::PEAKFITTER,
            "ML mode requested but HAS_ONNX not compiled in; falling back to TSpectrum");
#endif
    }

    if (!mlUsed) runTSpectrum();

    h_work->SetTitle(hname.c_str());
    h_work->GetXaxis()->SetTitle("Energy (keV)");
    h_work->GetYaxis()->SetTitle("Counts / 1ms per bin");
    h_work->SetMinimum(0);
    h_work->SetLineColor(kBlack);
    h_work->SetMarkerSize(0);
    h_work->Draw("hist");
    h_work->Draw("E1 same");
    c->Modified();
    c->Update();

    Debug::DumpHistogram(fout, h_work);

    std::sort(peaks.begin(), peaks.end());

    // -------------------------------------------------
    // Peak grouping
    // -------------------------------------------------
    auto groups = PeakGrouper::Group(peaks, res);

    out << "=====================================\n";
    out << "Histogram: " << hname << "\n";
    out << "Found peaks: " << peaks.size() << "\n";
    out << "Fit groups : " << groups.size() << "\n";
    out << "=====================================\n\n";

    int peakIndex = 0;

    std::vector<double> fittedEnergies;  // collect for IdentifyIsotopes at end
    double fwhmSum = 0.0;
    int    fwhmCount = 0;
    // Tracks DB lines already assigned to a peak so each line is claimed once.
    std::set<std::string> claimedLines;

    // =================================================
    // LOOP OVER GROUPS
    // =================================================
    for (size_t g = 0; g < groups.size(); g++) {

        
        auto& group = groups[g];

        Debug::LogGroup(g, group.energies);

        out << "\n=====================================\n";
        out << "GROUP " << g
            << " (" << group.energies.size()
            << " peaks)\n";

        for (double e : group.energies)
            out << "  Seed Peak: " << e << " keV\n";

        out << "=====================================\n";

        // -------------------------------------------------
        // S/N pre-screen: skip groups with insufficient signal-to-noise.
        // Computed on the bg-subtracted working histogram as:
        //   SN = Σ(peak window ±2σ) / √(scaled sideband |counts|)
        // References: gnuScope GetSumFF; Helmer & McCullagh 1979
        // -------------------------------------------------
        if (bg.snMinRatio > 0.0) {
            double groupSN = 0.0;
            for (double E : group.energies) {
                double sigma  = AdaptiveFitter::Sigma_(res, E);
                int bPlo = h_work->FindBin(E - 2.0*sigma);
                int bPhi = h_work->FindBin(E + 2.0*sigma);
                int bSlo1 = h_work->FindBin(E - 7.0*sigma);
                int bSlo2 = h_work->FindBin(E - 4.0*sigma);
                int bShi1 = h_work->FindBin(E + 4.0*sigma);
                int bShi2 = h_work->FindBin(E + 7.0*sigma);

                double sumPeak = 0.0;
                for (int b = bPlo; b <= bPhi; b++)
                    sumPeak += std::max(h_work->GetBinContent(b), 0.0);

                double sumSB = 0.0; int nSB = 0;
                for (int b = bSlo1; b <= bSlo2; b++) { sumSB += std::abs(h_work->GetBinContent(b)); nSB++; }
                for (int b = bShi1; b <= bShi2; b++) { sumSB += std::abs(h_work->GetBinContent(b)); nSB++; }

                int nPeak = bPhi - bPlo + 1;
                double noise = (nSB > 0 && nPeak > 0)
                    ? std::sqrt(sumSB / nSB * nPeak)
                    : std::sqrt(std::max(sumPeak, 1.0));
                double sn = (noise > 0.0) ? sumPeak / noise : 0.0;
                groupSN = std::max(groupSN, sn);
            }
            if (groupSN < bg.snMinRatio) {
                out << "  [S/N=" << groupSN << " < " << bg.snMinRatio << "] group skipped\n";
                Debug::Log(Debug::PEAKFITTER, "Group " + std::to_string(g)
                    + " skipped: S/N=" + std::to_string(groupSN)
                    + " < " + std::to_string(bg.snMinRatio));
                continue;
            }
        }

        // -------------------------------------------------
        // Adaptive fit
        // -------------------------------------------------
        TFitResultPtr fitResult;

        int nUsed = 1;

        AdaptiveFitterConfig fitCfg;
        fitCfg.useLogLikelihood = bg.useLogLikelihood;
        fitCfg.useImprove       = bg.useImprove;

        if (bg.useMLSigma) {
            for (const auto& E : group.energies) {
                double sigML = GetMLSigmaWidth(h_work, E, bg.mlSigmaModelPath);
                fitCfg.mlSigmaHints.push_back(sigML);
            }
        }

        TF1* model =
            AdaptiveFitter::FitGroup(
                h_work,
                group.energies,
                res,
                fitResult,
                nUsed,
                fitdb,
                fitCfg);

        if (!model) {
            Debug::Log(Debug::PEAKFITTER,
                "Group " + std::to_string(g) + ": group fit failed — retrying each peak individually");
            // Fallback: fit each peak alone using TSpectrum sigma × bin-width as
            // the resolution estimate, so every TSpectrum-found peak gets stored.
            if (fitdb) {
                for (double E : group.energies) {
                    double bw     = h_work->GetBinWidth(h_work->FindBin(E));
                    ConstantResModel cres;
                    cres.fwhm = bg.tspecSigma * bw * 2.3548200450309493;
                    TFitResultPtr fr;
                    int nu = 1;
                    TF1* fm = AdaptiveFitter::FitGroup(h_work, {E}, cres, fr, nu, fitdb);
                    Debug::Log(Debug::PEAKFITTER,
                        "  Fallback peak E=" + std::to_string(E) +
                        "  fwhm_est=" + std::to_string(cres.fwhm) +
                        (fm ? "  stored" : "  failed"));
                    delete fm;
                }
            }
            continue;
        }

        Debug::DrawFitComponents(fout, h_work, model,
                         Form("fit_group_%ld", g));

        // -------------------------------------------------
        // Fit quality — Pearson chi2/ndf from pull RMS (count-independent).
        // r->Chi2() from a log-likelihood fit returns -2*NLL, not chi2.
        // -------------------------------------------------
        const double sigL = res.FWHM(group.energies.front()) / 2.3548;
        const double sigR = res.FWHM(group.energies.back())  / 2.3548;
        const double xlo  = group.energies.front() - 4.0 * sigL;
        const double xhi  = group.energies.back()  + 4.0 * sigR;
        ResidualMetrics rm = FitDatabase::ComputeResiduals(h_work, model, xlo, xhi);
        double chi2ndf = (rm.rms < 1.0e6) ? FitDatabase::Chi2Ndf(rm, model->GetNpar()) : -1.0;

        Debug::Log(Debug::PEAKFITTER,
            "Group " + std::to_string(g) +
            "  nUsed=" + std::to_string(nUsed) +
            "  chi2/ndf=" + std::to_string(chi2ndf) +
            "  status=" + std::to_string(fitResult.Get() ? fitResult->Status() : -1) +
            "  edm=" + std::to_string(fitResult.Get() ? fitResult->Edm() : -1.0));

        // -------------------------------------------------
        // Background parameters
        // DG model (7 params) has bg at [5],[6]; all other models (step, tail,
        // standard) have bg at [3*nUsed], [3*nUsed+1].
        // -------------------------------------------------
        bool isDG    = (model->GetNpar() == 7);
        int bgOffset = isDG ? 5 : 3 * nUsed;

        double bg0 = model->GetParameter(bgOffset);
        double bg1 = model->GetParameter(bgOffset + 1);

        // =================================================
        // LOOP OVER INDIVIDUAL PEAKS
        // =================================================
        for (int i = 0; i < nUsed; i++) {
            // -------------------------------------------------
            // Parameter indices
            // -------------------------------------------------
            int p0 = 3 * i;

            double A =
                model->GetParameter(p0);

            double E =
                model->GetParameter(p0 + 1);

            double sig =
                model->GetParameter(p0 + 2);

            double Aerr =
                model->GetParError(p0);

            double Eerr =
                model->GetParError(p0 + 1);

            double sigerr =
                model->GetParError(p0 + 2);

            Debug::LogPeak(g, i, E, sig);
            Debug::Log(Debug::PEAKFITTER,
                "  Peak " + std::to_string(i) +
                "  A=" + std::to_string(A) +
                " ± " + std::to_string(Aerr) +
                "  E=" + std::to_string(E) +
                " ± " + std::to_string(Eerr) +
                "  sig=" + std::to_string(sig) +
                " ± " + std::to_string(sigerr));

            // -------------------------------------------------
            // Derived quantities
            // -------------------------------------------------
            double FWHM = 2.355 * sig;

            // Gaussian amplitude A is in counts/bin; integral over energy axis
            // is A*sigma*sqrt(2pi) keV, so divide by bin width (keV/bin) to get counts.
            double bw = h_work->GetBinWidth(h_work->FindBin(E));
            if (bw <= 0.0) bw = 1.0;

            double counts = A * sig * std::sqrt(2.0 * TMath::Pi()) / bw;

            // Full error propagation including A-sigma covariance from MIGRAD.
            // Cov(A,sigma) is typically negative (amplitude and width anticorrelate),
            // which slightly reduces the area uncertainty vs. ignoring the term.
            double counts_err = 0.0;
            if (A > 0 && sig > 0) {
                double relA   = Aerr   / A;
                double relSig = sigerr / sig;
                double cov_Asig = 0.0;
                if (fitResult.Get() && fitResult->IsValid())
                    cov_Asig = fitResult->GetCovarianceMatrix()(p0, p0 + 2);
                double varRel = relA * relA + relSig * relSig
                                + 2.0 * cov_Asig / (A * sig);
                if (varRel > 0.0)
                    counts_err = counts * std::sqrt(varRel);
            }

            // -------------------------------------------------
            // SNR estimate
            // -------------------------------------------------
            double bkgCounts =
                std::abs((bg0 + bg1 * E) * (5.0 * sig)) / bw;

            double SNR =
                (bkgCounts > 0)
                ? counts / std::sqrt(bkgCounts)
                : 0;

            // -------------------------------------------------
            // Net peak area (model-independent sum method)
            // When ML BG available use raw h - bkgML; else h_work - polynomial.
            // Statistical uncertainty from raw histogram h (full Poisson counts).
            // -------------------------------------------------
            double netArea    = 0.0;
            double netAreaErr = 0.0;
            {
                double xLo = E - 3.0 * sig;
                double xHi = E + 3.0 * sig;
                int bLo = h_work->FindBin(xLo + 0.5 * bw);
                int bHi = h_work->FindBin(xHi - 0.5 * bw);
                bLo = std::max(bLo, 1);
                bHi = std::min(bHi, h_work->GetNbinsX());
                double rawVar = 0.0;
                for (int b = bLo; b <= bHi; b++) {
                    double x = h_work->GetBinCenter(b);
                    if (bkgML)
                        netArea += h->GetBinContent(b) - bkgML->GetBinContent(b);
                    else
                        netArea += h_work->GetBinContent(b) - (bg0 + bg1 * x);
                    rawVar += std::max(h->GetBinContent(b), 0.0);
                }
                netAreaErr = (rawVar > 0.0) ? std::sqrt(rawVar) : 0.0;
            }

            // Flag if Gaussian integral and net area disagree by >5 %.
            // Possible causes: non-Gaussian tail, bad background, unresolved doublet.
            double areaDiscrepancy = (counts > 0.0)
                ? std::abs(netArea - counts) / counts : 0.0;
            bool   areaWarning = (areaDiscrepancy > 0.05);

            // -------------------------------------------------
            // Total area via TF1::Integral (Method 1)
            // Integrates the full model (Gaussian + step + tail) over ±5σ and
            // subtracts the polynomial BG contribution.  For singlets exact;
            // for multi-peak groups a small tail from neighbours may be included.
            // -------------------------------------------------
            double totalCounts = 0.0, totalCounts_err = 0.0;
            {
                double xlo5 = E - 5.0 * sig, xhi5 = E + 5.0 * sig;
                double bgI = (bg0 * (xhi5 - xlo5)
                              + 0.5 * bg1 * (xhi5*xhi5 - xlo5*xlo5)) / bw;
                totalCounts = model->Integral(xlo5, xhi5) / bw - bgI;
                if (fitResult.Get() && fitResult->IsValid()) {
                    TMatrixDSym cov = fitResult->GetCovarianceMatrix();
                    totalCounts_err = std::abs(
                        model->IntegralError(xlo5, xhi5,
                                             nullptr,
                                             cov.GetMatrixArray())) / bw;
                }
            }

            // -------------------------------------------------
            // COM centroid cross-check (Method 2)
            // Barycenter over ±2σ window on BG-subtracted h_work.
            // -------------------------------------------------
            double xCOM    = E;
            bool   comValid = false;
            if (SNR >= 3.0) {
                double sumW = 0.0, sumXW = 0.0;
                int bLo2 = std::max(h_work->FindBin(E - 2.0*sig), 1);
                int bHi2 = std::min(h_work->FindBin(E + 2.0*sig), h_work->GetNbinsX());
                for (int b = bLo2; b <= bHi2; b++) {
                    double net = h_work->GetBinContent(b);
                    if (net > 0.0) {
                        sumXW += h_work->GetBinCenter(b) * net;
                        sumW  += net;
                    }
                }
                if (sumW > 0.0) { xCOM = sumXW / sumW; comValid = true; }
            }
            bool comShift = comValid && (std::fabs(xCOM - E) > 2.0);

            // -------------------------------------------------
            // Data FWHM via half-maximum crossing (Method 3, singlets only)
            // -------------------------------------------------
            double dataFWHM      = 0.0;
            bool   dataFWHMValid = false;
            if (nUsed == 1 && SNR >= 3.0) {
                int    peakBin = h_work->FindBin(E);
                double peakH   = h_work->GetBinContent(peakBin);
                double halfH   = peakH * 0.5;
                if (halfH > 0.0) {
                    double xLeft = h_work->GetBinCenter(std::max(peakBin - 1, 1));
                    for (int b = peakBin - 1; b >= 1; b--) {
                        if (h_work->GetBinContent(b) <= halfH) {
                            double y1 = h_work->GetBinContent(b);
                            double y2 = h_work->GetBinContent(b + 1);
                            double x1 = h_work->GetBinCenter(b);
                            double x2 = h_work->GetBinCenter(b + 1);
                            double d  = y2 - y1;
                            xLeft = (std::fabs(d) > 0.0)
                                ? x1 + (halfH - y1) / d * (x2 - x1)
                                : (x1 + x2) * 0.5;
                            break;
                        }
                    }
                    double xRight = h_work->GetBinCenter(
                                        std::min(peakBin + 1, h_work->GetNbinsX()));
                    for (int b = peakBin + 1; b <= h_work->GetNbinsX(); b++) {
                        if (h_work->GetBinContent(b) <= halfH) {
                            double y1 = h_work->GetBinContent(b - 1);
                            double y2 = h_work->GetBinContent(b);
                            double x1 = h_work->GetBinCenter(b - 1);
                            double x2 = h_work->GetBinCenter(b);
                            double d  = y2 - y1;
                            xRight = (std::fabs(d) > 0.0)
                                ? x1 + (halfH - y1) / d * (x2 - x1)
                                : (x1 + x2) * 0.5;
                            break;
                        }
                    }
                    dataFWHM      = xRight - xLeft;
                    dataFWHMValid = (dataFWHM > 0.1 * sig);
                }
            }
            bool dataFWHMMismatch = dataFWHMValid
                && (FWHM > 0.0) && (std::fabs(dataFWHM - FWHM) / FWHM > 0.10);

            // -------------------------------------------------
            // ML centroid refinement (Method 5, singlets only)
            // -------------------------------------------------
            double mlCentroid    = E;
            bool   mlCentValid   = false;
#ifdef HAS_ONNX
            if (nUsed == 1 && bg.useML && SNR >= 3.0) {
                double offset = GetMLCentroidOffset(h_work, E, sig,
                                                    bg.mlCentroidModelPath);
                if (std::isfinite(offset)) {
                    mlCentroid = E + offset;
                    mlCentValid = true;
                }
            }
#endif

            // -------------------------------------------------
            // Matching — each DB line claimed by at most one peak
            // -------------------------------------------------
            double fwhmAtE = res.FWHM(E);
            auto matches = db.Match(E, fwhmAtE);
            // Remove DB lines already claimed by an earlier peak
            matches.erase(
                std::remove_if(matches.begin(), matches.end(),
                    [&](const GammaMatch& m) {
                        return claimedLines.count(
                            m.isotope + Form("_%.1f", m.energy));
                    }),
                matches.end());
            // Claim the best remaining line
            if (!matches.empty())
                claimedLines.insert(matches[0].isotope + Form("_%.1f", matches[0].energy));

            fittedEnergies.push_back(E);
            fwhmSum += fwhmAtE;
            fwhmCount++;

            Debug::Log(Debug::PEAKFITTER,
                "  counts=" + std::to_string(counts) +
                " ± " + std::to_string(counts_err) +
                "  net=" + std::to_string(netArea) +
                " ± " + std::to_string(netAreaErr) +
                (areaWarning ? "  [AREA_MISMATCH]" : "") +
                "  SNR=" + std::to_string(SNR) +
                "  FWHM=" + std::to_string(2.355*sig) +
                "  nMatches=" + std::to_string(matches.size()) +
                (matches.empty() ? "" : "  best=" + matches[0].isotope));

            // -------------------------------------------------
            // Label
            // -------------------------------------------------
            std::string label =
                Form("%.2f", E);

            if (!matches.empty()) {

                label += " ";
                label += matches[0].isotope;
            }

            // =================================================
            // OUTPUT
            // =================================================
            out << "\n-------------------------------------\n";
            out << "Peak " << peakIndex++ << "\n";
            out << "-------------------------------------\n";

            out << "Fit Group Size : "
                << nUsed << "\n";

            out << "Fit Status     : "
                << fitResult->Status() << "\n";

            out << "EDM            : "
                << fitResult->Edm() << "\n";

            out << "Chi2/NDF       : "
                << chi2ndf << "\n";

            out << "Amplitude      : "
                << A << " ± "
                << Aerr << "\n";

            out << "Energy         : "
                << E << " ± "
                << Eerr << " keV\n";

            if (comValid)
                out << "COM Centroid   : "
                    << xCOM << " keV"
                    << (comShift
                        ? Form("  [SHIFT %.2f keV]", xCOM - E)
                        : "")
                    << "\n";

            if (mlCentValid)
                out << "ML Centroid    : "
                    << mlCentroid << " keV\n";

            out << "Sigma          : "
                << sig << " ± "
                << sigerr << "\n";

            out << "FWHM           : "
                << FWHM << "\n";

            if (dataFWHMValid)
                out << "Data FWHM      : "
                    << dataFWHM << " keV"
                    << (dataFWHMMismatch
                        ? Form("  [MISMATCH %.1f%%]",
                               std::fabs(dataFWHM - FWHM) / FWHM * 100.0)
                        : "")
                    << "\n";

            out << "Peak Counts    : "
                << counts << " ± "
                << counts_err << "\n";

            out << "Total Counts   : "
                << totalCounts << " ± "
                << totalCounts_err
                << (nUsed > 1 ? "  [MULTI-PEAK: approx]" : "")
                << "\n";

            out << "Net Area       : "
                << netArea << " ± "
                << netAreaErr
                << (bkgML ? "  [ML-BG]" : "")
                << (areaWarning
                    ? Form("  [MISMATCH %.1f%%]", areaDiscrepancy * 100.0)
                    : "")
                << "\n";

            out << "SNR            : "
                << SNR << "\n";

            out << "BG Constant    : "
                << bg0 << "\n";

            out << "BG Slope       : "
                << bg1 << "\n";

            // Step and tail parameters (present only in enhanced models).
            if (!isDG) {
                int npar   = model->GetNpar();
                int nStep  = (npar - 3*nUsed - 2);       // total extra params
                // step only: nStep == nUsed; tail only: nStep == 2*nUsed; both: 3*nUsed
                bool hasStep = (npar == 4*nUsed + 2 || npar == 6*nUsed + 2);
                bool hasTail = (npar == 5*nUsed + 2 || npar == 6*nUsed + 2);
                (void)nStep;
                if (hasStep) {
                    int sIdx = 3*nUsed + 2 + i;
                    out << "Step Height    : " << model->GetParameter(sIdx) << "\n";
                }
                if (hasTail) {
                    int nS   = hasStep ? nUsed : 0;
                    int tIdx = 3*nUsed + 2 + nS + i;
                    int bIdx = 3*nUsed + 2 + nS + nUsed + i;
                    out << "Tail Amplitude : " << model->GetParameter(tIdx) << "\n";
                    out << "Tail Slope     : " << model->GetParameter(bIdx) << "\n";
                }
            }

            // -------------------------------------------------
            // Matches
            // -------------------------------------------------
            out << "Possible Matches:\n";

            if (matches.empty()) {

                out << "   None\n";

            } else {

                for (auto& m : matches) {

                    out << "   "
                        << m.isotope
                        << " | "
                        << m.energy
                        << " keV"
                        << " | dE = "
                        << m.deltaE
                        << " keV"
                        << " | I = " << m.intensity
                        << " | score = " << m.score
                        << "\n";
                }

                out << "Best Match: "
                    << matches[0].isotope
                    << " ("
                    << matches[0].energy
                    << " keV)\n";
            }

            // -------------------------------------------------
            // Save fit
            // -------------------------------------------------
            storage.Save(
                hname,
                {
                    E,
                    sig,
                    counts,
                    chi2ndf,
                    nUsed
                });

            // -------------------------------------------------
            // Tracking
            // -------------------------------------------------
            if (enableTracking && tracker) {

                tracker->Add(
                    E,
                    0,
                    counts,
                    counts_err,
                    sig,
                    sigerr,
                    fitResult->Status(),
                    fitResult->Edm());
            }

            // -------------------------------------------------
            // Label
            // -------------------------------------------------
            double y = model->Eval(E);

            latex.DrawLatex(
                E,
                y * 1.08,
                label.c_str());
        }

        model->Draw("same");

        // Flush the GUI event queue so the canvas repaints after each group
        if (extCanvas) {
            c->Modified();
            c->Update();
            gSystem->ProcessEvents();
        }
    }
    // -------------------------------------------------
    // Multi-line isotope identification for full spectrum
    // -------------------------------------------------
    if (!fittedEnergies.empty()) {

        double avgFWHM = (fwhmCount > 0) ? fwhmSum / fwhmCount : 5.0;
        auto candidates = db.IdentifyIsotopes(fittedEnergies, avgFWHM);

        out << "\n=====================================\n";
        out << "ISOTOPE IDENTIFICATION (multi-line)\n";
        out << "=====================================\n";

        int rank = 1;
        for (auto& c : candidates) {
            if (rank > 10) break;
            out << rank++ << ". "
                << c.isotope
                << "  lines=" << c.nMatched
                << "  score=" << c.totalScore
                << "  matched keV:";
            for (double e : c.matchedEnergies)
                out << " " << e;
            out << "\n";
        }
    }

    tracker->FitResolutionModel(res, fout, hname);
    c->Modified();
    c->Update();

 
    // -------------------------------------------------
    // Save canvas to peakFits directory
    // -------------------------------------------------
    RootFileManager::SaveCanvas(fout, c);
    RootFileManager::SaveHistogram(fout, h);
    RootFileManager::SaveBgSub(fout, h_work);
 
    fout->cd();



    if (bkgML) { delete bkgML; bkgML = nullptr; }
    out.close();
}