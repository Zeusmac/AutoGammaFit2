#include "RootFileManager.h"
#include "Debug.h"

#include "TFile.h"
#include "TDirectory.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TGraphErrors.h"

using namespace std;

// ============================================================
// Get or create directory
// ============================================================

TDirectory* RootFileManager::GetOrCreate(
    TFile* fout,
    const std::string& dirName)
{
    if (!fout)
        return nullptr;

    TDirectory* dir =
        fout->GetDirectory(dirName.c_str());

    if (!dir)
        dir = fout->mkdir(dirName.c_str());

    return dir;
}

// ============================================================
// Generic object saver
// ============================================================

void RootFileManager::SaveObject(
    TFile* fout,
    TObject* obj,
    const std::string& dirName)
{
    if (!fout || !obj)
        return;

    TDirectory* dir =
        GetOrCreate(fout, dirName);

    if (!dir)
        return;

    dir->cd();
    dir->WriteObject(obj, obj->GetName());
    fout->cd();

    Debug::Log(Debug::FILEIO,
        "Saved " + std::string(obj->ClassName()) +
        " \"" + std::string(obj->GetName()) +
        "\"  -> dir=" + dirName);
}

// ============================================================
// Canvas saver
// ============================================================

void RootFileManager::SaveCanvas(
    TFile* fout,
    TCanvas* canvas,
    const std::string& dirName)
{
    SaveObject(fout,
               canvas,
               dirName);
}

// ============================================================
// Histogram saver
// ============================================================

void RootFileManager::SaveHistogram(
    TFile* fout,
    TH1* hist,
    const std::string& dirName)
{
    SaveObject(fout,
               hist,
               dirName);
}

// ============================================================
// Graph saver
// ============================================================

void RootFileManager::SaveGraph(
    TFile* fout,
    TGraphErrors* graph,
    const std::string& dirName)
{
    SaveObject(fout,
               graph,
               dirName);
}
// ============================================================
// Bg_subHist saver
// ============================================================

void RootFileManager::SaveBgSub(
    TFile* fout,
    TH1* Hist,
    const std::string& dirName)
{
    SaveObject(fout,
               Hist,
               dirName);
}