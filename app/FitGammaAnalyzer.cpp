// app/FitGammaAnalyzer.cpp

#include "PeakFitter.h"
#include "GammaDB.h"
#include "TKey.h"
#include "TFile.h"
#include "TH1.h"
#include "FitStorage.h"
#include "PeakTracker.h"
#include "Debug.h"
#include <iostream>
#include <string>
#include <sys/stat.h>

static const char* kCacheDir = "fit_caches";

void RunAnalysis(const char* rootfile, TFile* fout, bool cacheOnly) {

    // Create cache directory if it doesn't exist
    mkdir(kCacheDir, 0755);

    GammaDB db;
    ResolutionModel res;
    FitStorage storage;
    PeakTracker tracker;

    if (!db.Load("../Isotope_energys.txt")) return;

    PeakFitter fitter(db, &tracker, res, storage, nullptr);

    TFile *file = new TFile(rootfile);
    if (!file || file->IsZombie()) return;

    TIter next(file->GetListOfKeys());
    TKey *key;
    for (int iter = 0; iter < 1; iter++) {
        while ((key = (TKey*)next())) {
            TObject *obj = key->ReadObj();

            fout->cd();

            if (!obj->InheritsFrom("TH1")) continue;

            TH1* h = (TH1*)obj;
            std::string hname = h->GetName();
            std::string cacheFile = std::string(kCacheDir) + "/fit_cache_" + hname + ".dat";

            FitDatabase fitdb;
            fitdb.Load(cacheFile);
            fitdb.cacheOnly = cacheOnly;

            fitter.SetFitDatabase(&fitdb);
            fitter.FitHistogram(h, fout, true);

            // In cache-only mode the cache is read-only — don't overwrite it
            if (!cacheOnly) fitdb.Save(cacheFile);
        }
    }

    fout->Close();

    std::cout << "Analysis complete\n";
}


static void PrintHelp(const char* prog) {
    std::cout <<
"Usage:\n"
"  " << prog << " <input.root> [options]\n"
"\n"
"Arguments:\n"
"  <input.root>          ROOT file containing TH1 gamma spectra to fit.\n"
"\n"
"Options:\n"
"  -h, --help            Show this help message and exit.\n"
"\n"
"  -c, --cache-only      Skip MIGRAD entirely. Draw fit results directly from\n"
"                        the per-histogram cache files (fit_caches/).  Groups\n"
"                        with no cached entry are silently skipped. The cache\n"
"                        files are NOT updated in this mode.\n"
"\n"
"  --debug=<sections>    Enable debug output for one or more comma-separated\n"
"                        sections.  Available sections:\n"
"\n"
"                          FITTER     AdaptiveFitter: seed source (DB vs fresh),\n"
"                                     each retry stage (stage1/2/3), AIC model\n"
"                                     comparison, parameter bounds.\n"
"\n"
"                          GROUPER    PeakGrouper: every merge/split decision,\n"
"                                     gap vs 4-sigma threshold, total groups.\n"
"\n"
"                          TRACKER    PeakTracker: why each point is filtered\n"
"                                     (FWHM cut, EDM cut, boundary artifact),\n"
"                                     weighted mean per energy, FWHM fit result,\n"
"                                     resolution model update.\n"
"\n"
"                          DB         FitDatabase: every cache load, save, seed\n"
"                                     attempt, and store-if-better decision.\n"
"\n"
"                          GAMMADB    GammaDB: isotope line loading, per-peak\n"
"                                     match results with score, IdentifyIsotopes\n"
"                                     top candidates.\n"
"\n"
"                          RESMODEL   ResolutionModel: old parameters, raw fit\n"
"                                     values, and exponentially-smoothed result\n"
"                                     for each histogram update.\n"
"\n"
"                          FILEIO     RootFileManager: every object written to\n"
"                                     AllGammaFits.root (canvas, histogram,\n"
"                                     graph, debug fits).\n"
"\n"
"                          PEAKFITTER PeakFitter main loop: TSpectrum peak count,\n"
"                                     per-group fit quality (chi2, status, EDM),\n"
"                                     per-peak E/sigma/counts/SNR/best match.\n"
"\n"
"                        Examples:\n"
"                          --debug=FITTER\n"
"                          --debug=FITTER,DB,TRACKER\n"
"\n"
"  --debug-all           Enable ALL debug sections at once.\n"
"\n"
"Output:\n"
"  AllGammaFits.root     ROOT file with fit canvases, raw histograms,\n"
"                        background-subtracted histograms, and FWHM vs Energy\n"
"                        graphs for each input histogram.\n"
"\n"
"  fit_caches/           Directory of per-histogram cache files\n"
"                        (fit_cache_<name>.dat).  Each file stores the best\n"
"                        converged fit parameters seen so far.  Warm-starts\n"
"                        subsequent runs and feeds --cache-only mode.\n"
"\n"
"  ../Gamma_fits/        Text files (<name>_fit.txt) with per-peak results:\n"
"                        energy, sigma, FWHM, counts, SNR, chi2/NDF, EDM,\n"
"                        isotope matches, and multi-line isotope identification.\n"
"\n"
"Examples:\n"
"  " << prog << " run42.root\n"
"  " << prog << " run42.root --cache-only\n"
"  " << prog << " run42.root --debug=FITTER,DB\n"
"  " << prog << " run42.root --cache-only --debug=PEAKFITTER\n"
"  " << prog << " run42.root --debug-all\n"
"\n";
}

int main(int argc, char** argv)
{
    // Check for --help before anything else
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") {
            PrintHelp(argv[0]);
            return 0;
        }
    }

    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <input.root> [options]\n"
                  << "       " << argv[0] << " --help  for full option list\n";
        return 1;
    }

    bool cacheOnly = false;
    for (int i = 2; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--cache-only" || arg == "-c") {
            cacheOnly = true;
        } else if (arg == "--debug-all") {
            Debug::SetAll(true);
        } else if (arg.rfind("--debug=", 0) == 0) {
            Debug::EnableFromString(arg.substr(8));
        } else {
            std::cerr << "Unknown option: " << arg
                      << "  (run with --help for usage)\n";
        }
    }

    Debug::PrintConfig();

    if (cacheOnly)
        std::cout << "Cache-only mode: using stored fits, no MIGRAD.\n";

    TFile* fout = new TFile("AllGammaFits.root", "RECREATE");

    RunAnalysis(argv[1], fout, cacheOnly);

    std::cout << "Analysis complete\n";

    fout->Close();
    return 0;
}
