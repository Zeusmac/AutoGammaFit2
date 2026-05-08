#ifndef PEAKFITRESULT_H
#define PEAKFITRESULT_H

struct PeakFitResult {
    double energy;
    double sigma;
    double area;
    double chi2;
    int nPeaks;
};

#endif