#ifndef ADAPTIVE_FITTER_H
#define ADAPTIVE_FITTER_H

#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "TH1.h"
#include "TF1.h"
#include "TSpectrum.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "TMath.h"
#include "Debug.h"

class AdaptiveFitter {

public:

    // =========================================================
    // Build N-peak Gaussian + linear background model
    // =========================================================
    static TF1* Build(int nPeaks,
                      double xmin,
                      double xmax)
    {
        std::string formula;

        for (int i = 0; i < nPeaks; i++) {

            formula +=
                "gaus(" +
                std::to_string(3*i) +
                ")+";
        }

        formula +=
            "pol1(" +
            std::to_string(3*nPeaks) +
            ")";

        TF1* f =
            new TF1(
                Form("adaptiveFit_%d", nPeaks),
                formula.c_str(),
                xmin,
                xmax);

        return f;
    }

    // =========================================================
    // Adaptive multi-peak fitting
    // =========================================================
    template<typename ResModel>
    static TF1* FitGroup(
        TH1* h,
        const std::vector<double>& peaks,
        const ResModel& res,
        TFitResultPtr& bestFit,
        int& nUsed)
    {
        if (!h || peaks.empty())
            return nullptr;

        // -----------------------------------------------------
        // Estimate expected detector resolution
        // -----------------------------------------------------
        double meanE = 0;

        for (auto p : peaks)
            meanE += p;

        meanE /= peaks.size();

        double expectedSigma =
            res.Sigma(meanE);

        double expectedFWHM =
            2.355 * expectedSigma;

        // -----------------------------------------------------
        // Estimate group width
        // -----------------------------------------------------

        // -----------------------------------------------------
        // Decide allowed complexity
        // -----------------------------------------------------
        int maxPeaks = peaks.size();
        //if (maxPeaks > 1){
        double groupWidth = peaks.back() - peaks.front();
        

        // if (groupWidth >
        //     0.8 * expectedFWHM)
        // {cd ,
        //     maxPeaks = 2;
        // }

        // if (groupWidth >
        //     1.8 * expectedFWHM)
        // {
        //     maxPeaks = 3;
        // }

        maxPeaks =
            std::min(
                maxPeaks,
                (int)peaks.size());

        // -----------------------------------------------------
        // Debug
        // -----------------------------------------------------
        Debug::LogResolutionDecision(groupWidth,
                                    expectedFWHM,
                                    maxPeaks);

        // -----------------------------------------------------
        // Best-fit tracking
        // -----------------------------------------------------
        TF1* best = nullptr;

        double bestChi2 = 1e99;

        nUsed = 1;

        double previousChi2 = 1e99;

        // =====================================================
        // Try 1,2,3 peak models
        // =====================================================
        for (int n = 1;
             n <= maxPeaks;
             n++)
        {
            double xmin =
                peaks.front()
                - 10.0 * expectedFWHM;

            double xmax =
                peaks.back()
                + 10.0 * expectedFWHM;

            TF1* f =
                Build(
                    n,
                    xmin,
                    xmax);

            // -------------------------------------------------
            // Initialize peaks
            // -------------------------------------------------
            for (int i = 0; i < n; i++) {

                double E =
                    peaks[i];

                double sigma =
                    fabs(res.Sigma(E));

                double A =
                    h->GetBinContent(
                        h->FindBin(E));

                // ---------------------------------------------
                // Amplitude
                // ---------------------------------------------
                f->SetParameter(
                    3*i,
                    A);

                f->SetParLimits(
                    3*i,
                    0,
                    50*A);

                // ---------------------------------------------
                // Centroid
                // ---------------------------------------------
                f->SetParameter(
                    3*i+1,
                    E);

                f->SetParLimits(
                    3*i+1,
                    E - 10*sigma,
                    E + 10*sigma);

                // ---------------------------------------------
                // Sigma
                // ---------------------------------------------
                f->SetParameter(
                    3*i+2,
                    sigma);

                f->SetParLimits(
                    3*i+2,
                    0.3*sigma,
                    1.5*sigma);
            }

            // -------------------------------------------------
            // Background initialization
            // -------------------------------------------------
            int bgIndex =
                3*n;

            f->SetParameter(
                bgIndex,
                0);

            f->SetParameter(
                bgIndex+1,
                0);

            // -------------------------------------------------
            // Perform fit
            // -------------------------------------------------
            TFitResultPtr r =
                h->Fit(
                    f,
                    "R S");

            Debug::LogFitFailure(n,r ? r->Status() : -1);
        
            // -------------------------------------------------
            // Evaluate fit quality
            // -------------------------------------------------
            double chi2 = 1e99;

            if (r)
            {
                chi2 = r->Chi2();

                if (r->Ndf() > 0) chi2 /= r->Ndf();
            }   

            double improvement =
                (previousChi2 - chi2)
                / previousChi2;

            // -------------------------------------------------
            // Measure broadening
            // -------------------------------------------------
            double avgSigma = 0;

            for (int i = 0; i < n; i++)
                avgSigma +=
                    f->GetParameter(3*i+2);

            avgSigma /= n;

            double measuredFWHM =
                2.355 * avgSigma;

            double resolutionRatio =
                measuredFWHM
                / expectedFWHM;

            // -------------------------------------------------
            // Debug output
            // -------------------------------------------------
            Debug::LogFitAttempt(n,
                            chi2,
                            improvement,
                            resolutionRatio,
                            r ? r->Status() : -1);

            // -------------------------------------------------
            // Decide whether to accept model
            // -------------------------------------------------

            bool accept = false;

            // Always accept first fit
            if (n == 1)
                accept = true;

            // Require significant improvement
            else if (improvement > 0.15)
                accept = true;

            // Force escalation if broadened
            if (resolutionRatio > 1.4)
                accept = true;

            // -------------------------------------------------
            // Save best model
            // -------------------------------------------------
            if (accept &&
                chi2 < bestChi2)
            {
                bestChi2 = chi2;

                best = f;

                bestFit = r;

                nUsed = n;
            }

            previousChi2 = chi2;
        }

        // =====================================================
        // Residual peak check (optional future upgrade)
        // =====================================================
        // Could add:
        //
        // 1. build residual histogram
        // 2. run TSpectrum on residual
        // 3. escalate fit automatically
        //
        // =====================================================

        return best;
    }
};

#endif