#ifndef GAMMAFITGUI_H
#define GAMMAFITGUI_H

// ROOT GUI
#include "TGFrame.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"
#include "TGListBox.h"
#include "TGComboBox.h"
#include "TGTab.h"
#include "TGTextView.h"
#include "TGStatusBar.h"
#include "TRootEmbeddedCanvas.h"

// Project
#include "GammaDB.h"
#include "ResolutionModel.h"
#include "FitStorage.h"
#include "PeakTracker.h"
#include "FitDatabase.h"
#include "Debug.h"

class TH1;
class TFile;
class TCanvas;
class TF1;

// ─────────────────────────────────────────────────────────────────────────────
// GammaFitGUI
//
// Main interactive analysis window.
//
// Layout
//   Left panel  (310 px)  — TGTab with "AutoFit" and "Manual Fit" control tabs
//   Right panel (rest)    — shared spectrum canvas + log strip at the bottom
//   Status bar  (bottom)  — current file / action
//
// Extending the GUI
//   • AutoFit options  → BuildAutoFitTab()   (add checkboxes / number entries)
//   • Manual analysis  → BuildManualFitTab() (add a new TGGroupFrame at the bottom)
//   • New tab          → ctrlTab->AddTab() in the constructor + a BuildXxxTab()
// ─────────────────────────────────────────────────────────────────────────────
class GammaFitGUI : public TGMainFrame {
public:
    explicit GammaFitGUI(const TGWindow* p,
                         UInt_t w = 1400, UInt_t h = 920);
    ~GammaFitGUI() override;

    // ── File & histogram ──────────────────────────────────────────────────────
    void OnOpenFile();
    void OnOpenIsotopeDB();
    void OnHistogramSelected(Int_t id);   // listbox in AutoFit tab

    // ── AutoFit tab ───────────────────────────────────────────────────────────
    void OnRunSelected();
    void OnRunAll();
    void OnDebugAllOn();
    void OnDebugAllOff();

    // ── Manual Fit tab ────────────────────────────────────────────────────────
    void OnLoadManual();
    void OnCanvasEvent(Int_t event, Int_t px, Int_t py, TObject* obj);
    void OnSeedParams();
    void OnPreview();
    void OnManualFit();
    void OnAcceptFit();
    void OnRejectFit();

    // ── Utilities ─────────────────────────────────────────────────────────────
    void AppendLog(const std::string& msg);
    void SetStatus(const std::string& msg);

    ClassDefOverride(GammaFitGUI, 0)

private:
    // ── Analysis objects ──────────────────────────────────────────────────────
    GammaDB       db_;
    ResolutionModel res_;
    FitStorage    storage_;
    PeakTracker   tracker_;
    bool          dbLoaded_ = false;

    // ── Run-time state ────────────────────────────────────────────────────────
    TFile*       inputFile_   = nullptr;
    std::string  inputPath_;
    std::vector<std::string> histNames_;
    std::string  currentHist_;
    TH1*         rawHist_    = nullptr;
    TF1*         manualTF1_  = nullptr;
    double       clickEnergy_ = 500.0;

    // ── Shared widgets ────────────────────────────────────────────────────────
    TGStatusBar*         statusBar_ = nullptr;
    TGTextView*          logView_   = nullptr;
    TRootEmbeddedCanvas* canvas_    = nullptr;

    // ── AutoFit tab widgets ───────────────────────────────────────────────────
    std::string    isotopePath_;
    TGLabel*       fileLbl_      = nullptr;
    TGLabel*       isotopeLbl_   = nullptr;
    TGListBox*     histList_     = nullptr;
    TGCheckButton* cacheOnlyChk_ = nullptr;
    TGCheckButton* useSeedsChk_  = nullptr;
    TGCheckButton* debugChk_[8]  = {};

    // ── Manual Fit tab widgets ────────────────────────────────────────────────
    TGComboBox*    manualCombo_  = nullptr;
    TGNumberEntry* mEnergy_      = nullptr;
    TGNumberEntry* mSigma_       = nullptr;
    TGNumberEntry* mAmp_         = nullptr;
    TGNumberEntry* mBg0_         = nullptr;
    TGNumberEntry* mBg1_         = nullptr;
    TGNumberEntry* mRange_       = nullptr;
    TGLabel*       mResultLbl_   = nullptr;

    // ── Build helpers ─────────────────────────────────────────────────────────
    void BuildAutoFitTab  (TGCompositeFrame* parent);
    void BuildManualFitTab(TGCompositeFrame* parent);
    void PopulateHistWidgets();
    void SyncDebugToggles();
    void DrawOnCanvas(TH1* h, TF1* fit = nullptr);
    void RunFitOnHistogram(const std::string& hname);
    std::string CacheFileFor(const std::string& hname) const;
};

#endif
