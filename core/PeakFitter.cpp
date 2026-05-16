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


#include <fstream>
#include <algorithm>
#include <set>

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
    TSpectrum spectrum(1000);

    // Background subtraction is independent of peak search
    if (bg.subtractBg) {
        TH1* bkg = spectrum.Background(h_work, bg.iterations);
        h_work->Add(bkg, -1);
    }

    // Peak positions: either from TSpectrum or caller-supplied forced seeds
    std::vector<double> peaks;
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
        // -------------------------------------------------
        int bgOffset = 3 * nUsed;

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

            double counts =
                A * sig *
                std::sqrt(2 * TMath::Pi());

            double counts_err = 0;

            if (A > 0 && sig > 0) {

                counts_err =
                    counts *
                    std::sqrt(
                        std::pow(Aerr / A, 2) +
                        std::pow(sigerr / sig, 2));
            }

            // -------------------------------------------------
            // SNR estimate
            // -------------------------------------------------
            double bkgCounts =
                std::abs(
                    bg0 * (5 * sig) +
                    bg1 * E * (5 * sig));

            double SNR =
                (bkgCounts > 0)
                ? counts / std::sqrt(bkgCounts)
                : 0;

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

            out << "Sigma          : "
                << sig << " ± "
                << sigerr << "\n";

            out << "FWHM           : "
                << FWHM << "\n";

            out << "Peak Counts    : "
                << counts << " ± "
                << counts_err << "\n";

            out << "SNR            : "
                << SNR << "\n";

            out << "BG Constant    : "
                << bg0 << "\n";

            out << "BG Slope       : "
                << bg1 << "\n";

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



    out.close();   
   
}