#include "PeakFitter.h"
#include "PeakGrouper.h"
#include "AdaptiveFitter.h"
#include "PeakTracker.h"

#include "Debug.h"
#include "TSpectrum.h"
#include "TLatex.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TFitResult.h"
#include "RootFileManager.h"


#include <fstream>
#include <algorithm>

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
                              TCanvas* extCanvas)
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
        c = new TCanvas(Form("c_%s", hname.c_str()),
                        hname.c_str(),
                        1200, 800);
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
    // Peak search
    // -------------------------------------------------
    TSpectrum spectrum(1000);

    int nPeaks =
        spectrum.Search(h_work, 2, "", 0.02);

    Debug::Log(Debug::PEAKFITTER,
        "TSpectrum found " + std::to_string(nPeaks) + " peaks in " + hname);

    TH1* bkg = spectrum.Background(h_work, 14);
    h_work->Add(bkg, -1);

    Debug::DumpHistogram(fout, h_work);

    std::vector<double> peaks(
        spectrum.GetPositionX(),
        spectrum.GetPositionX() + nPeaks);

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
        // Adaptive fit
        // -------------------------------------------------
        TFitResultPtr fitResult;

        int nUsed = 1;

        TF1* model =
            AdaptiveFitter::FitGroup(
                h_work,
                group.energies,
                res,
                fitResult,
                nUsed,
                fitdb);

        if (!model) {
            Debug::Log(Debug::PEAKFITTER,
                "Group " + std::to_string(g) + ": no model returned — skipping");
            continue;
        }

        Debug::DrawFitComponents(fout, h_work, model,
                         Form("fit_group_%ld", g));

        // -------------------------------------------------
        // Fit quality
        // -------------------------------------------------
        double chi2ndf = -1;

        if (fitResult.Get() && fitResult->Ndf() > 0) {
            chi2ndf = fitResult->Chi2() / fitResult->Ndf();
        }

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
            // Matching
            // -------------------------------------------------
            double fwhmAtE = res.FWHM(E);
            auto matches = db.Match(E, fwhmAtE);

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