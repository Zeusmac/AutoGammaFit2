#include "GammaFitGUI.h"
#include "TApplication.h"
#include "TGClient.h"
#include <iostream>

int main(int argc, char** argv)
{
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            std::cout <<
"Usage:\n"
"  gamma_gui [input.root]\n"
"\n"
"  Launches the interactive AutoGammaFit 2.0 GUI.\n"
"\n"
"  AutoFit tab\n"
"    Open ROOT File    -  browse for a .root file containing TH1 spectra\n"
"    Histogram list    -  click to preview; single-click selects for fitting\n"
"    Cache-Only        -  skip MIGRAD; draw stored fit parameters directly\n"
"    Use Cached Seeds  -  warm-start MIGRAD from the best previous parameters\n"
"    Debug sections    -  enable per-module verbose output (stdout)\n"
"    Run AutoFit (selected)  -  fit the highlighted histogram\n"
"    Run AutoFit (ALL)       -  fit every histogram in sequence\n"
"\n"
"  Manual Fit tab\n"
"    Load to Canvas    -  display the chosen histogram\n"
"    Click on spectrum  -  places peak marker and seeds parameters\n"
"    Seed from Model   -  fills sigma from resolution model, amp from bin height\n"
"    Preview           -  draw a Gaussian at the current parameters (no fitting)\n"
"    Run Fit           -  run MIGRAD on the region Energy +/- RangexSigma\n"
"    Accept & Save     -  write the fit parameters to the histogram cache file\n"
"    Reject            -  discard the result and restore the histogram\n"
"\n"
"  Output\n"
"    AllGammaFits.root          -  ROOT file with fit canvases and histograms\n"
"    fit_caches/<name>.dat      -  per-histogram parameter caches\n"
"    ../Gamma_fits/<name>.txt   -  plain-text per-peak results\n"
"\n";
            return 0;
        }
    }

    TApplication app("GammaFitApp", &argc, argv);
    new GammaFitGUI(gClient->GetRoot(), 1400, 920);
    app.Run();
    return 0;
}
