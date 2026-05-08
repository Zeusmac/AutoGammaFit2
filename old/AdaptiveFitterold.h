#ifndef ADAPTIVEFITTER_H
#define ADAPTIVEFITTER_H

#include "TF1.h"
#include "TH1.h"
#include "TFitResultPtr.h"
#include "TFitResult.h"
#include "Debug.h"
#include <iostream>
#include <vector>

class AdaptiveFitter {
public:

    static TF1* Build(int nPeaks, double xmin, double xmax) {

        std::string formula;

        for (int i = 0; i < nPeaks; i++)
            formula += "gaus(" + std::to_string(3*i) + ")+";

        formula += "pol1(" + std::to_string(3*nPeaks) + ")";

        return new TF1("model", formula.c_str(), xmin, xmax);
    }

    template<typename ResModel>
    static TF1* FitGroup(TH1* h,
                         const std::vector<double>& peaks,
                         const ResModel& res,
                         TFitResultPtr& bestFit,
                         int& nUsed) {

        int maxPeaks = std::min((int)peaks.size(), 3);
        std::cout << maxPeaks << std::endl;
        TF1* best = nullptr;
        double bestChi2 = 1e99;
        double fwhm = res.FWHM(peaks[0]);        

        TF1* f = Build(maxPeaks, peaks.front()-5, peaks.back()+5);
            
            // initialize params
        for (int i = 0; i < maxPeaks; i++) {
            double E = peaks[i];
            double sigma = res.Sigma(E);
            double A = h->GetBinContent(h->FindBin(E));

            f->SetParameter(3*i+0, A);
            f->SetParameter(3*i+1, E);
            f->SetParameter(3*i+2, sigma);
            
            f->SetParLimits(3*i+2,0,10);
            f->SetParLimits(3*i+1,E - 8,E + 8);
            f->SetParLimits(3*i+0, 0.5*A, 10*A);
        }
        f->SetParameter(3*maxPeaks, 0);
        f->SetParLimits(3*maxPeaks, 0, 1e7);
        f->SetParLimits(3*maxPeaks+1, -1e4, 0);
        TFitResultPtr r = h->Fit(f, "R S");
        Debug::LogFitFailure(maxPeaks,r ? r->Status() : -1);
        double chi2 = (r && r->Ndf() > 0) ? r->Chi2() / r->Ndf(): 0.0;
        bestFit = r;
        nUsed = maxPeaks;
        best = f;
            // if (chi2 < bestChi2) {
            //     bestChi2 = chi2;
            //     best = f;
            //     bestFit = r;
            //     nUsed = n;
            // }
        

    return best;
    }

};

#endif