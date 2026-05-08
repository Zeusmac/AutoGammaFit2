#ifndef ROOTFILEMANAGER_H
#define ROOTFILEMANAGER_H

class TFile;
class TDirectory;
class TObject;
class TCanvas;
class TH1;
class TGraph;
class TGraphErrors;

#include <string>

namespace RootFileManager {

    // --------------------------------------------------------
    // Generic directory creation
    // --------------------------------------------------------
    TDirectory* GetOrCreate(TFile* fout,
                            const std::string& dirName);

    // --------------------------------------------------------
    // Save generic ROOT object
    // --------------------------------------------------------
    void SaveObject(TFile* fout,
                    TObject* obj,
                    const std::string& dirName);

    // --------------------------------------------------------
    // Save canvas
    // --------------------------------------------------------
    void SaveCanvas(TFile* fout,
                    TCanvas* canvas,
                    const std::string& dirName = "peakFits");

    // --------------------------------------------------------
    // Save histogram
    // --------------------------------------------------------
    void SaveHistogram(TFile* fout,
                       TH1* hist,
                       const std::string& dirName = "rawHistograms");

    // --------------------------------------------------------
    // Save graph
    // --------------------------------------------------------
    void SaveGraph(TFile* fout,
                   TGraphErrors* graph,
                   const std::string& dirName = "FWHMvsEnergy");
    // --------------------------------------------------------
    // Background sub
    // --------------------------------------------------------
    void SaveBgSub(TFile* fout,
                   TH1* Hist,
                   const std::string& dirName = "Bg_subhist");

}

#endif