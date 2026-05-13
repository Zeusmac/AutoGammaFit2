#include "GammaFitGUI.h"
#include "GammaFitGUI_shared.h"

ClassImp(GammaFitGUI)

#include "TGFileDialog.h"
#include "TGMsgBox.h"
#include "TGTextEntry.h"
#include "TMath.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TPad.h"
#include "TH1.h"
#include "TH1D.h"
#include "TFile.h"
#include "TKey.h"
#include "TF1.h"
#include "TLine.h"
#include "TLatex.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"

#include "PeakFitter.h"
#include "TSpectrum.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <numeric>
#include <cstdio>
#include <set>
#include <sys/stat.h>

// ─────────────────────────────────────────────────────────────────────────────
// Constructor
// ─────────────────────────────────────────────────────────────────────────────
GammaFitGUI::GammaFitGUI(const TGWindow* p, UInt_t w, UInt_t h)
    : TGMainFrame(p, w, h)
{
    SetWindowName("AutoGammaFit 2.0  —  Interactive Fitting GUI");
    SetCleanup(kDeepCleanup);

    // ── Main horizontal split: controls (left) | canvas+log (right) ──────────
    TGHorizontalFrame* main = new TGHorizontalFrame(this, w, h - 25);
    AddFrame(main, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    // ── LEFT control panel ────────────────────────────────────────────────────
    TGVerticalFrame* left = new TGVerticalFrame(main, 312, h - 25);
    main->AddFrame(left, new TGLayoutHints(kLHintsLeft | kLHintsExpandY, 2, 0, 2, 2));

    TGTab* ctrlTab = new TGTab(left, 310, h - 28);
    left->AddFrame(ctrlTab, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    BuildAutoFitTab   (ctrlTab->AddTab("AutoFit"));
    BuildSourceTab    (ctrlTab->AddTab("Source"));
    BuildManualFitTab (ctrlTab->AddTab("Manual Fit"));
    BuildFWHMTab      (ctrlTab->AddTab("FWHM"));
    BuildDecayTab     (ctrlTab->AddTab("Decay"));
    BuildFitResultsTab(ctrlTab->AddTab("Fit Results"));
    BuildIsotopesTab  (ctrlTab->AddTab("Isotopes"));

    // ── RIGHT: canvas + log ───────────────────────────────────────────────────
    TGVerticalFrame* right = new TGVerticalFrame(main, w - 318, h - 25);
    main->AddFrame(right, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 2, 2, 2));

    canvas_ = new TRootEmbeddedCanvas("GUICanvas", right, w - 322, h - 230);
    right->AddFrame(canvas_, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));

    // Connect mouse clicks for manual fit peak placement
    canvas_->GetCanvas()->Connect(
        "ProcessedEvent(Int_t,Int_t,Int_t,TObject*)",
        "GammaFitGUI", this,
        "OnCanvasEvent(Int_t,Int_t,Int_t,TObject*)");

    // Button row between canvas and log
    {
        TGHorizontalFrame* btnRow = new TGHorizontalFrame(right);
        right->AddFrame(btnRow, new TGLayoutHints(kLHintsRight, 2, 2, 0, 0));

        TGTextButton* clearCacheBtn = new TGTextButton(btnRow, "  Clear Cache  ");
        btnRow->AddFrame(clearCacheBtn, new TGLayoutHints(kLHintsRight, 2, 2, 1, 1));
        clearCacheBtn->Connect("Clicked()", "GammaFitGUI", this, "OnClearHistCache()");
        clearCacheBtn->SetToolTipText("Delete ALL cache entries for the currently displayed histogram");

        TGTextButton* savePlotBtn = new TGTextButton(btnRow, "  Save Plot  ");
        btnRow->AddFrame(savePlotBtn, new TGLayoutHints(kLHintsRight, 2, 2, 1, 1));
        savePlotBtn->Connect("Clicked()", "GammaFitGUI", this, "OnSavePlot()");
        savePlotBtn->SetToolTipText(
            "Save the current canvas to a file.\n"
            "Supported formats: .pdf  .png  .root");
    }

    // Log strip below canvas
    TGGroupFrame* logGrp = new TGGroupFrame(right, "Log");
    right->AddFrame(logGrp, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
    logView_ = new TGTextView(logGrp, w - 330, 155);
    logGrp->AddFrame(logView_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    // ── Status bar ────────────────────────────────────────────────────────────
    statusBar_ = new TGStatusBar(this, w, 22, kHorizontalFrame);
    Int_t parts[] = {55, 30, 15};
    statusBar_->SetParts(parts, 3);
    AddFrame(statusBar_, new TGLayoutHints(kLHintsBottom | kLHintsExpandX));

    // ── Finalise ──────────────────────────────────────────────────────────────
    MapSubwindows();
    Resize(GetDefaultSize());
    MapWindow();

    // Load isotope database
    isotopePath_ = kIsotopeDBDefault;
    if (db_.Load(isotopePath_)) {
        dbLoaded_ = true;
        AppendLog("Isotope DB loaded: " + isotopePath_ +
                  "  (" + std::to_string(db_.db.size()) + " lines)");
        isotopeLbl_->SetText(isotopePath_.c_str());
        PopulateIsoDbList();
    } else {
        AppendLog("WARNING: Could not load isotope DB from " +
                  isotopePath_ + " — matches will be empty");
        isotopeLbl_->SetText("(not loaded)");
    }

    AppendLog("AutoGammaFit 2.0 GUI ready.  Open a ROOT file to begin.");
    SetStatus("Ready");
}

GammaFitGUI::~GammaFitGUI()
{
    if (rawHistOwned_) { delete rawHist_; rawHist_ = nullptr; }
    if (inputFile_)   { inputFile_->Close();   delete inputFile_; }
    if (srcRootFile_) { srcRootFile_->Close(); delete srcRootFile_; }
    delete bgTF1_;
    delete viewHist_;
    delete fwhmTF1_;
    for (TF1* c : fitComponents_) delete c;
    fitComponents_.clear();
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildAutoFitTab
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::BuildAutoFitTab(TGCompositeFrame* p)
{
    // ── File ──────────────────────────────────────────────────────────────────
    TGGroupFrame* fg = new TGGroupFrame(p, "File");
    p->AddFrame(fg, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

    TGTextButton* openBtn = new TGTextButton(fg, "Open ROOT File...");
    fg->AddFrame(openBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    openBtn->Connect("Clicked()", "GammaFitGUI", this, "OnOpenFile()");
    openBtn->SetToolTipText("Browse for a ROOT file containing TH1 gamma spectra");

    fileLbl_ = new TGLabel(fg, "No file loaded");
    fg->AddFrame(fileLbl_, new TGLayoutHints(kLHintsLeft, 4, 4, 0, 4));

    TGTextButton* isoBtn = new TGTextButton(fg, "Open Isotope DB...");
    fg->AddFrame(isoBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 2));
    isoBtn->Connect("Clicked()", "GammaFitGUI", this, "OnOpenIsotopeDB()");
    isoBtn->SetToolTipText("Browse for the isotope energy database text file");

    isotopeLbl_ = new TGLabel(fg, "(not loaded)");
    fg->AddFrame(isotopeLbl_, new TGLayoutHints(kLHintsLeft, 4, 4, 0, 4));

    // ── Histogram list ────────────────────────────────────────────────────────
    TGGroupFrame* hg = new TGGroupFrame(p, "Histograms");
    p->AddFrame(hg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    histList_ = new TGListBox(hg, 100);
    histList_->Resize(285, 130);
    hg->AddFrame(histList_, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
    histList_->Connect("Selected(Int_t)", "GammaFitGUI", this, "OnHistogramSelected(Int_t)");
    histList_->SetMultipleSelections(kFALSE);

    // ── Classification row ────────────────────────────────────────────────────
    {
        TGHorizontalFrame* cr = new TGHorizontalFrame(hg);
        hg->AddFrame(cr, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        cr->AddFrame(new TGLabel(cr, "Class:"),
                     new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        histClassCombo_ = new TGComboBox(cr, 900);
        histClassCombo_->AddEntry("Gamma Spectrum", 1);
        histClassCombo_->AddEntry("Decay Curve",    2);
        histClassCombo_->AddEntry("2D Histogram",   3);
        histClassCombo_->Select(1, kFALSE);
        histClassCombo_->Resize(140, 22);
        cr->AddFrame(histClassCombo_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        TGTextButton* setClassBtn = new TGTextButton(cr, "Set");
        cr->AddFrame(setClassBtn, new TGLayoutHints(kLHintsLeft));
        setClassBtn->Connect("Clicked()", "GammaFitGUI", this, "OnHistClassSet()");
        setClassBtn->SetToolTipText("Assign the chosen classification to the selected histogram");
    }

    // ── Fit options ───────────────────────────────────────────────────────────
    TGGroupFrame* og = new TGGroupFrame(p, "Fit Options");
    p->AddFrame(og, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    TGTextButton* loadCacheSel = new TGTextButton(og, "Load Cache  (selected)");
    og->AddFrame(loadCacheSel, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 1));
    loadCacheSel->Connect("Clicked()", "GammaFitGUI", this, "OnLoadCacheSelected()");
    loadCacheSel->SetToolTipText("Display cached fits on the selected histogram without re-running AutoFit");

    TGTextButton* loadCacheAll = new TGTextButton(og, "Load Cache  (ALL histograms)");
    og->AddFrame(loadCacheAll, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 4));
    loadCacheAll->Connect("Clicked()", "GammaFitGUI", this, "OnLoadCacheAll()");
    loadCacheAll->SetToolTipText("Add all histograms that have a cache file to the Fit Results list");

    useSeedsChk_ = new TGCheckButton(og, "Use Cached Seeds");
    useSeedsChk_->SetState(kButtonDown);
    og->AddFrame(useSeedsChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 1, 2));
    useSeedsChk_->SetToolTipText(
        "Warm-start MIGRAD from the best previously converged parameters.\n"
        "Uncheck to always fit from scratch.");

    autoLogLikChk_ = new TGCheckButton(og, "Log-likelihood  (L)");
    autoLogLikChk_->SetState(kButtonDown);
    og->AddFrame(autoLogLikChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 1));
    autoLogLikChk_->SetToolTipText(
        "Use Poisson log-likelihood (better for low-count bins).\n"
        "Uncheck to use chi2 / least-squares (more reliable chi2/ndf for high-count spectra).");

    autoImprovChk_ = new TGCheckButton(og, "IMPROVE  (M)");
    og->AddFrame(autoImprovChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 4));
    autoImprovChk_->SetToolTipText("Run IMPROVE after MIGRAD to search for a better minimum");

    // ── Background subtraction ────────────────────────────────────────────────
    TGGroupFrame* bgGrpAuto = new TGGroupFrame(p, "Background Subtraction");
    p->AddFrame(bgGrpAuto, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    bgSubtractChk_ = new TGCheckButton(bgGrpAuto, "Enable TSpectrum background removal");
    bgSubtractChk_->SetState(kButtonDown);
    bgGrpAuto->AddFrame(bgSubtractChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
    bgSubtractChk_->SetToolTipText(
        "Subtract the estimated background from the spectrum before fitting.\n"
        "Disable if the spectrum is already background-free.");

    {
        TGHorizontalFrame* iterRow = new TGHorizontalFrame(bgGrpAuto);
        bgGrpAuto->AddFrame(iterRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
        TGLabel* iterLbl = new TGLabel(iterRow, "Iterations:");
        iterRow->AddFrame(iterLbl, new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
        bgIterEntry_ = new TGNumberEntry(iterRow, 14, 4, -1,
            TGNumberFormat::kNESInteger,
            TGNumberFormat::kNEAPositive,
            TGNumberFormat::kNELLimitMinMax, 2, 100);
        iterRow->AddFrame(bgIterEntry_, new TGLayoutHints(kLHintsLeft));
        bgIterEntry_->GetNumberEntry()->SetToolTipText("Number of TSpectrum background estimation iterations (default 14)");
    }

    {
        TGHorizontalFrame* btnRow = new TGHorizontalFrame(bgGrpAuto);
        bgGrpAuto->AddFrame(btnRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));

        TGTextButton* bgSel = new TGTextButton(btnRow, " Preview (selected) ");
        btnRow->AddFrame(bgSel, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        bgSel->Connect("Clicked()", "GammaFitGUI", this, "OnApplyBgSelected()");
        bgSel->SetToolTipText("Draw the selected histogram with background subtracted on the canvas (preview only).");

        TGTextButton* bgAll = new TGTextButton(btnRow, " Apply to All (Gamma) ");
        btnRow->AddFrame(bgAll, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        bgAll->Connect("Clicked()", "GammaFitGUI", this, "OnApplyBgAll()");
        bgAll->SetToolTipText(
            "Save background-subtracted histograms classified as 'Gamma Spectrum' to AllGammaFits.root.");

        TGTextButton* bgReset = new TGTextButton(btnRow, " Reset ");
        btnRow->AddFrame(bgReset, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
        bgReset->Connect("Clicked()", "GammaFitGUI", this, "OnResetBgSub()");
        bgReset->SetToolTipText("Reset to defaults: 14 iterations, subtraction enabled.");
    }

    // ── Peak Finding ──────────────────────────────────────────────────────────
    TGGroupFrame* pfg = new TGGroupFrame(p, "Peak Finding");
    p->AddFrame(pfg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    {
        TGHorizontalFrame* r1 = new TGHorizontalFrame(pfg);
        pfg->AddFrame(r1, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 1));
        r1->AddFrame(new TGLabel(r1, "Sigma (bins):"),
                     new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
        tspecSigmaEntry_ = new TGNumberEntry(r1, 2.0, 5, -1,
            TGNumberFormat::kNESRealOne,
            TGNumberFormat::kNEAPositive,
            TGNumberFormat::kNELLimitMinMax, 0.5, 20.0);
        r1->AddFrame(tspecSigmaEntry_, new TGLayoutHints(kLHintsLeft));
        tspecSigmaEntry_->GetNumberEntry()->SetToolTipText(
            "TSpectrum peak width in bins. Lower = finds narrower peaks.");
    }
    {
        TGHorizontalFrame* r2 = new TGHorizontalFrame(pfg);
        pfg->AddFrame(r2, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 4));
        r2->AddFrame(new TGLabel(r2, "Threshold:"),
                     new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
        tspecThreshEntry_ = new TGNumberEntry(r2, 0.02, 6, -1,
            TGNumberFormat::kNESRealThree,
            TGNumberFormat::kNEAPositive,
            TGNumberFormat::kNELLimitMinMax, 0.001, 0.99);
        r2->AddFrame(tspecThreshEntry_, new TGLayoutHints(kLHintsLeft));
        tspecThreshEntry_->GetNumberEntry()->SetToolTipText(
            "Minimum peak height as a fraction of the tallest peak (0.001–0.99).\n"
            "Lower = finds weaker peaks. Default 0.02 = peaks ≥ 2% of max.");
    }

    // ── Run ───────────────────────────────────────────────────────────────────
    TGGroupFrame* rg = new TGGroupFrame(p, "Run");
    p->AddFrame(rg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

    TGTextButton* runSel = new TGTextButton(rg, "Run AutoFit  (selected)");
    rg->AddFrame(runSel, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    runSel->Connect("Clicked()", "GammaFitGUI", this, "OnRunSelected()");
    runSel->SetToolTipText("Fit the histogram currently selected in the list");

    TGTextButton* runAll = new TGTextButton(rg, "Run AutoFit  (ALL Gamma histograms)");
    rg->AddFrame(runAll, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
    runAll->Connect("Clicked()", "GammaFitGUI", this, "OnRunAll()");
    runAll->SetToolTipText("Fit all histograms classified as 'Gamma Spectrum'");

    // ── Custom Projection ──────────────────────────────────────────────────────
    TGGroupFrame* cpg = new TGGroupFrame(p, "Custom Projection");
    p->AddFrame(cpg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

    {
        TGHorizontalFrame* r1 = new TGHorizontalFrame(cpg);
        cpg->AddFrame(r1, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 1));
        r1->AddFrame(new TGLabel(r1, "TH2:"),
                     new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        custProjTh2Combo_ = new TGComboBox(r1, 910);
        custProjTh2Combo_->Resize(200, 22);
        r1->AddFrame(custProjTh2Combo_, new TGLayoutHints(kLHintsExpandX));
    }
    {
        TGHorizontalFrame* r2 = new TGHorizontalFrame(cpg);
        cpg->AddFrame(r2, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
        r2->AddFrame(new TGLabel(r2, "Project:"),
                     new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        custProjAxisCombo_ = new TGComboBox(r2, 911);
        custProjAxisCombo_->AddEntry("X  (cut on Y axis)", 1);
        custProjAxisCombo_->AddEntry("Y  (cut on X axis)", 2);
        custProjAxisCombo_->Select(1, kFALSE);
        custProjAxisCombo_->Resize(170, 22);
        r2->AddFrame(custProjAxisCombo_, new TGLayoutHints(kLHintsLeft));
    }
    {
        TGHorizontalFrame* r3 = new TGHorizontalFrame(cpg);
        cpg->AddFrame(r3, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
        r3->AddFrame(new TGLabel(r3, "Range:"),
                     new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        custProjLo_ = new TGNumberEntry(r3, 0.0, 7, -1, TGNumberFormat::kNESRealThree);
        custProjLo_->SetWidth(75);
        r3->AddFrame(custProjLo_, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
        r3->AddFrame(new TGLabel(r3, "to"), new TGLayoutHints(kLHintsCenterY, 2, 2, 0, 0));
        custProjHi_ = new TGNumberEntry(r3, 100.0, 7, -1, TGNumberFormat::kNESRealThree);
        custProjHi_->SetWidth(75);
        r3->AddFrame(custProjHi_, new TGLayoutHints(kLHintsLeft, 2, 0, 0, 0));
    }
    {
        TGHorizontalFrame* r4 = new TGHorizontalFrame(cpg);
        cpg->AddFrame(r4, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 2));
        r4->AddFrame(new TGLabel(r4, "Name:"),
                     new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        custProjName_ = new TGTextEntry(r4, "myProj");
        custProjName_->SetWidth(140);
        r4->AddFrame(custProjName_, new TGLayoutHints(kLHintsExpandX, 0, 4, 0, 0));
        TGTextButton* addBtn = new TGTextButton(r4, "Add");
        r4->AddFrame(addBtn, new TGLayoutHints(kLHintsRight));
        addBtn->Connect("Clicked()", "GammaFitGUI", this, "OnAddCustomProjection()");
        addBtn->SetToolTipText("Create this projection and add it to the histogram list");
    }

    // ── Debug toggles ─────────────────────────────────────────────────────────
    TGGroupFrame* dg = new TGGroupFrame(p, "Debug Sections");
    p->AddFrame(dg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

    static const char* kSecNames[8] = {
        "FITTER", "GROUPER", "TRACKER", "DB",
        "GAMMADB", "RESMODEL", "FILEIO", "PEAKFITTER"
    };
    static const char* kSecTips[8] = {
        "AdaptiveFitter: seeds, retry stages, AIC",
        "PeakGrouper: merge/split decisions",
        "PeakTracker: filtering, FWHM fit, model update",
        "FitDatabase: cache load/save/seed/store",
        "GammaDB: loading, match results, isotope scoring",
        "ResolutionModel: parameter updates",
        "RootFileManager: every ROOT file write",
        "PeakFitter main loop: per-peak results"
    };

    TGHorizontalFrame* dcol = nullptr;
    for (int i = 0; i < 8; ++i) {
        if (i % 2 == 0) {
            dcol = new TGHorizontalFrame(dg);
            dg->AddFrame(dcol, new TGLayoutHints(kLHintsExpandX));
        }
        debugChk_[i] = new TGCheckButton(dcol, kSecNames[i], 200 + i);
        dcol->AddFrame(debugChk_[i], new TGLayoutHints(kLHintsLeft, 2, 6, 1, 1));
        debugChk_[i]->SetToolTipText(kSecTips[i]);
    }

    TGHorizontalFrame* dbRow = new TGHorizontalFrame(dg);
    dg->AddFrame(dbRow, new TGLayoutHints(kLHintsLeft, 2, 2, 4, 2));

    TGTextButton* allOn = new TGTextButton(dbRow, " All On ");
    dbRow->AddFrame(allOn, new TGLayoutHints(kLHintsLeft, 2, 6, 0, 0));
    allOn->Connect("Clicked()", "GammaFitGUI", this, "OnDebugAllOn()");

    TGTextButton* allOff = new TGTextButton(dbRow, " All Off ");
    dbRow->AddFrame(allOff, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
    allOff->Connect("Clicked()", "GammaFitGUI", this, "OnDebugAllOff()");
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildSourceTab
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::BuildSourceTab(TGCompositeFrame* p)
{
    // ── Spectrum (separate ROOT file) ──────────────────────────────────────────
    TGGroupFrame* hg = new TGGroupFrame(p, "Spectrum");
    p->AddFrame(hg, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

    // Row 1: open ROOT file
    {
        TGHorizontalFrame* row = new TGHorizontalFrame(hg);
        hg->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        TGTextButton* openBtn = new TGTextButton(row, " Open ROOT File ");
        row->AddFrame(openBtn, new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        openBtn->Connect("Clicked()", "GammaFitGUI", this, "OnOpenSourceRootFile()");
        openBtn->SetToolTipText("Open a ROOT file containing source spectra");

        srcRootFileLbl_ = new TGLabel(row, "(no file)");
        row->AddFrame(srcRootFileLbl_, new TGLayoutHints(kLHintsLeft | kLHintsCenterY));
    }

    // Row 2: histogram combo + Load Cache / AutoFit buttons
    {
        TGHorizontalFrame* row = new TGHorizontalFrame(hg);
        hg->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));

        srcHistCombo_ = new TGComboBox(row, 900);
        srcHistCombo_->Resize(130, 22);
        row->AddFrame(srcHistCombo_, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 0, 4, 0, 0));

        TGTextButton* cacheBtn = new TGTextButton(row, " Load Cache ");
        row->AddFrame(cacheBtn, new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        cacheBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLoadSourceCache()");
        cacheBtn->SetToolTipText("Load cached fit results for the selected histogram");

        TGTextButton* fitBtn = new TGTextButton(row, " AutoFit ");
        row->AddFrame(fitBtn, new TGLayoutHints(kLHintsCenterY, 0, 0, 0, 0));
        fitBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRunSourceAutoFit()");
        fitBtn->SetToolTipText(
            "Run AutoFit on the selected source histogram.\n"
            "Peaks are seeded from the loaded source description file,\n"
            "not discovered by TSpectrum.");
    }

    // Row 3: background subtraction options
    {
        TGHorizontalFrame* row = new TGHorizontalFrame(hg);
        hg->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));

        srcBgSubChk_ = new TGCheckButton(row, "BG subtract");
        srcBgSubChk_->SetState(kButtonDown);
        row->AddFrame(srcBgSubChk_, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 8, 0, 0));
        srcBgSubChk_->SetToolTipText("Apply TSpectrum background subtraction before fitting");

        TGLabel* itLbl = new TGLabel(row, "Iterations:");
        row->AddFrame(itLbl, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 4, 0, 0));
        srcBgIterEntry_ = new TGNumberEntry(row, 14, 4, -1,
            TGNumberFormat::kNESInteger,
            TGNumberFormat::kNEAPositive,
            TGNumberFormat::kNELLimitMinMax, 1, 100);
        row->AddFrame(srcBgIterEntry_, new TGLayoutHints(kLHintsLeft));
    }

    // ── Source file ────────────────────────────────────────────────────────────
    TGGroupFrame* sg = new TGGroupFrame(p, "Source Description");
    p->AddFrame(sg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    {
        TGHorizontalFrame* row = new TGHorizontalFrame(sg);
        sg->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        TGTextButton* loadBtn = new TGTextButton(row, " Load File ");
        row->AddFrame(loadBtn, new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        loadBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLoadSourceFile()");
        loadBtn->SetToolTipText(
            "Load a .txt source description file.\n"
            "Format (one keyword per line):\n"
            "  isotope   Co-60\n"
            "  activity  37000    # Bq at caldate\n"
            "  halflife  1925.5   # days\n"
            "  caldate   2020-01-15\n"
            "  measdate  2025-03-10\n"
            "  livetime  3600.0   # seconds\n"
            "  # energy_keV  branching_ratio\n"
            "  1173.228  0.9985\n"
            "  1332.492  0.9998");

        srcFileLbl_ = new TGLabel(row, "(no file loaded)");
        row->AddFrame(srcFileLbl_, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 0, 0, 0));
    }

    srcInfoLbl_ = new TGLabel(sg, "  —");
    sg->AddFrame(srcInfoLbl_, new TGLayoutHints(kLHintsLeft, 4, 4, 0, 2));

    {
        TGHorizontalFrame* row = new TGHorizontalFrame(sg);
        sg->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        TGLabel* lbl = new TGLabel(row, "Activity unit:");
        row->AddFrame(lbl, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 4, 0, 0));
        srcActivityUnit_ = new TGComboBox(row, 910);
        srcActivityUnit_->AddEntry("Bq",  1);
        srcActivityUnit_->AddEntry("µCi (x 37000 Bq)", 2);
        srcActivityUnit_->Select(1, kFALSE);
        srcActivityUnit_->Resize(160, 22);
        row->AddFrame(srcActivityUnit_, new TGLayoutHints(kLHintsLeft));
        srcActivityUnit_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                  "OnActivityUnitChanged(Int_t)");
    }

    {
        TGHorizontalFrame* row = new TGHorizontalFrame(sg);
        sg->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));
        TGLabel* lbl = new TGLabel(row, "Live time (s):");
        row->AddFrame(lbl, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 4, 0, 0));
        srcLiveTime_ = new TGNumberEntry(row, 1.0, 10, -1,
            TGNumberFormat::kNESReal,
            TGNumberFormat::kNEAPositive,
            TGNumberFormat::kNELLimitMinMax, 1e-6, 1e9);
        row->AddFrame(srcLiveTime_, new TGLayoutHints(kLHintsExpandX));
    }

    // ── Peak assignments ────────────────────────────────────────────────────────
    TGGroupFrame* ag = new TGGroupFrame(p, "Source Lines & Assignments");
    p->AddFrame(ag, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    srcLineList_ = new TGListBox(ag, 920);
    srcLineList_->Resize(290, 110);
    ag->AddFrame(srcLineList_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    {
        TGHorizontalFrame* row = new TGHorizontalFrame(ag);
        ag->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));

        TGLabel* lbl = new TGLabel(row, "Fitted E (keV):");
        row->AddFrame(lbl, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 4, 0, 0));
        srcManualE_ = new TGNumberEntry(row, 0.0, 8, -1,
            TGNumberFormat::kNESReal,
            TGNumberFormat::kNEAPositive,
            TGNumberFormat::kNELLimitMinMax, 0.0, 1e6);
        row->AddFrame(srcManualE_, new TGLayoutHints(kLHintsExpandX));
    }

    {
        TGHorizontalFrame* btnRow = new TGHorizontalFrame(ag);
        ag->AddFrame(btnRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));

        TGTextButton* autoBtn = new TGTextButton(btnRow, " Auto Identify ");
        btnRow->AddFrame(autoBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        autoBtn->Connect("Clicked()", "GammaFitGUI", this, "OnAutoIdentify()");
        autoBtn->SetToolTipText("Match each source line to the nearest cached fitted peak (within 3 FWHM)");

        TGTextButton* manBtn = new TGTextButton(btnRow, " Manual Assign ");
        btnRow->AddFrame(manBtn, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
        manBtn->Connect("Clicked()", "GammaFitGUI", this, "OnManualAssign()");
        manBtn->SetToolTipText(
            "Assign the energy typed above to the selected source line.\n"
            "Finds the closest cached peak within 5 FWHM of that energy.");
    }

    // ── Plots ──────────────────────────────────────────────────────────────────
    TGGroupFrame* pg = new TGGroupFrame(p, "Plots");
    p->AddFrame(pg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

    TGTextButton* calibBtn = new TGTextButton(pg, "Energy Calibration (fitted − reference)");
    pg->AddFrame(calibBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    calibBtn->Connect("Clicked()", "GammaFitGUI", this, "OnShowEnergyCalib()");
    calibBtn->SetToolTipText("Plot (fitted E − reference E) vs reference E for all assigned lines");

    TGTextButton* effBtn = new TGTextButton(pg, "Efficiency  vs  Energy");
    pg->AddFrame(effBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
    effBtn->Connect("Clicked()", "GammaFitGUI", this, "OnShowEfficiency()");
    effBtn->SetToolTipText(
        "Plot absolute full-energy peak efficiency = counts / (A_meas * BR * liveTime).\n"
        "Source activity is decay-corrected from caldate to measdate.");

    TGTextButton* fwhmBtn = new TGTextButton(pg, "FWHM  vs  Energy");
    pg->AddFrame(fwhmBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
    fwhmBtn->Connect("Clicked()", "GammaFitGUI", this, "OnShowSourceFWHM()");
    fwhmBtn->SetToolTipText("Load FWHM data from cache and draw FWHM vs Energy on the shared canvas");

}

// ─────────────────────────────────────────────────────────────────────────────
// BuildManualFitTab
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::BuildManualFitTab(TGCompositeFrame* p)
{
    // Wrap everything in a scrollable canvas so content isn't cut off
    TGCanvas* sc = new TGCanvas(p, 308, 860, kSunkenFrame);
    p->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    TGCompositeFrame* cf = new TGCompositeFrame(sc->GetViewPort(), 295, 10, kVerticalFrame);
    sc->SetContainer(cf);
    p = cf;   // all widgets added to the scrollable container from here on

    // ── Histogram selector ────────────────────────────────────────────────────
    TGGroupFrame* hg = new TGGroupFrame(p, "Histogram");
    p->AddFrame(hg, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

    {
        TGHorizontalFrame* comboRow = new TGHorizontalFrame(hg);
        hg->AddFrame(comboRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        manualCombo_ = new TGComboBox(comboRow, 300);
        manualCombo_->Resize(215, 22);
        comboRow->AddFrame(manualCombo_,
                           new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 0, 4, 0, 0));

        TGTextButton* loadBtn = new TGTextButton(comboRow, " Load ");
        comboRow->AddFrame(loadBtn, new TGLayoutHints(kLHintsCenterY, 0, 0, 0, 0));
        loadBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLoadManual()");
        loadBtn->SetToolTipText("Draw the selected histogram in the chosen view mode.");
    }

    // View-mode selector
    {
        TGHorizontalFrame* vrow = new TGHorizontalFrame(hg);
        hg->AddFrame(vrow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
        vrow->AddFrame(new TGLabel(vrow, "View:"),
                       new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        histViewCombo_ = new TGComboBox(vrow, 700);
        histViewCombo_->AddEntry("Raw histogram",              1);
        histViewCombo_->AddEntry("Background subtracted",      2);
        histViewCombo_->AddEntry("BG sub + AutoFit peaks",     3);
        histViewCombo_->Select(1, kFALSE);
        histViewCombo_->Resize(200, 22);
        vrow->AddFrame(histViewCombo_, new TGLayoutHints(kLHintsExpandX));
        histViewCombo_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                "OnHistViewChanged(Int_t)");
    }

    TGLabel* inst = new TGLabel(p, "  ◀ Click on spectrum to place peak markers");
    p->AddFrame(inst, new TGLayoutHints(kLHintsLeft, 4, 4, 0, 2));

    // ── Peak parameters ───────────────────────────────────────────────────────
    TGGroupFrame* pg = new TGGroupFrame(p, "Fit Parameters");
    p->AddFrame(pg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    // Helper lambda: one labelled number-entry row
    auto addRow = [&](TGGroupFrame* grp, const char* lbl,
                      TGNumberEntry*& entry, double val,
                      double lo, double hi)
    {
        TGHorizontalFrame* row = new TGHorizontalFrame(grp);
        grp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 0));
        TGLabel* l = new TGLabel(row, lbl);
        row->AddFrame(l, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 6, 0, 0));
        l->SetWidth(95);
        entry = new TGNumberEntry(row, val, 9, -1,
            TGNumberFormat::kNESRealThree,
            TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMinMax, lo, hi);
        row->AddFrame(entry, new TGLayoutHints(kLHintsRight | kLHintsExpandX));
    };

    addRow(pg, "Last E (keV)", mEnergy_, 500.0,  0.0,  10000.0);
    addRow(pg, "Sigma  (keV)", mSigma_,    1.5,  0.05,    50.0);
    addRow(pg, "Amplitude",    mAmp_,    500.0,  0.0,    1.0e9);
    addRow(pg, "BG const",     mBg0_,      0.0, -1.0e6,  1.0e6);
    addRow(pg, "BG slope",     mBg1_,      0.0, -1.0e4,  1.0e4);
    addRow(pg, "Range (×σ)",   mRange_,    4.0,  1.0,     20.0);

    TGTextButton* seedBtn = new TGTextButton(pg, "Seed from Resolution Model");
    pg->AddFrame(seedBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 6, 2));
    seedBtn->Connect("Clicked()", "GammaFitGUI", this, "OnSeedParams()");
    seedBtn->SetToolTipText("Fill Sigma from resolution model at the last clicked energy.\nFill Amplitude from the histogram bin height.");

    // ── Background region ─────────────────────────────────────────────────────
    TGGroupFrame* bgGrp = new TGGroupFrame(p, "Background Region");
    p->AddFrame(bgGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    {
        TGHorizontalFrame* bgRow = new TGHorizontalFrame(bgGrp);
        bgGrp->AddFrame(bgRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        bgRow->AddFrame(new TGLabel(bgRow, "Lo:"),
                        new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        mBgLo_ = new TGNumberEntry(bgRow, 0.0, 7, -1,
            TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMinMax, 0.0, 10000.0);
        bgRow->AddFrame(mBgLo_, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));

        bgRow->AddFrame(new TGLabel(bgRow, "Hi:"),
                        new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        mBgHi_ = new TGNumberEntry(bgRow, 0.0, 7, -1,
            TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMinMax, 0.0, 10000.0);
        bgRow->AddFrame(mBgHi_, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
    }
    {
        TGHorizontalFrame* bgBtnRow = new TGHorizontalFrame(bgGrp);
        bgGrp->AddFrame(bgBtnRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));

        TGTextButton* setBgBtn = new TGTextButton(bgBtnRow, "Set from Canvas");
        bgBtnRow->AddFrame(setBgBtn, new TGLayoutHints(kLHintsLeft, 2, 4, 0, 0));
        setBgBtn->Connect("Clicked()", "GammaFitGUI", this, "OnSetBgFromCanvas()");
        setBgBtn->SetToolTipText("Click twice on the spectrum to define the background region [Lo, Hi]");

        TGTextButton* fitBgBtn = new TGTextButton(bgBtnRow, "Fit Background");
        bgBtnRow->AddFrame(fitBgBtn, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));
        fitBgBtn->Connect("Clicked()", "GammaFitGUI", this, "OnFitBackground()");
        fitBgBtn->SetToolTipText("Fit a linear background in the selected region and seed BG parameters");

        TGTextButton* clearBgBtn = new TGTextButton(bgBtnRow, "Clear");
        bgBtnRow->AddFrame(clearBgBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
        clearBgBtn->Connect("Clicked()", "GammaFitGUI", this, "OnClearBackground()");
        clearBgBtn->SetToolTipText("Clear the background fit and reset Lo/Hi to zero");
    }

    // ── Fit range (2-click canvas selection) ─────────────────────────────────
    TGGroupFrame* rangeGrp = new TGGroupFrame(p, "Fit Range");
    p->AddFrame(rangeGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    {
        TGHorizontalFrame* rangeRow = new TGHorizontalFrame(rangeGrp);
        rangeGrp->AddFrame(rangeRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        rangeRow->AddFrame(new TGLabel(rangeRow, "Lo:"),
                           new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        mFitLo_ = new TGNumberEntry(rangeRow, 0.0, 7, -1,
            TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMinMax, 0.0, 10000.0);
        rangeRow->AddFrame(mFitLo_, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));

        rangeRow->AddFrame(new TGLabel(rangeRow, "Hi:"),
                           new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        mFitHi_ = new TGNumberEntry(rangeRow, 0.0, 7, -1,
            TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMinMax, 0.0, 10000.0);
        rangeRow->AddFrame(mFitHi_, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
    }
    {
        TGHorizontalFrame* rangeBtnRow = new TGHorizontalFrame(rangeGrp);
        rangeGrp->AddFrame(rangeBtnRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));

        TGTextButton* setRangeBtn = new TGTextButton(rangeBtnRow, "Set from Canvas");
        rangeBtnRow->AddFrame(setRangeBtn, new TGLayoutHints(kLHintsLeft, 2, 4, 0, 0));
        setRangeBtn->Connect("Clicked()", "GammaFitGUI", this, "OnSetRangeFromCanvas()");
        setRangeBtn->SetToolTipText("Click twice on the spectrum to define the fit range [Lo, Hi]");

        TGTextButton* clearRangeBtn = new TGTextButton(rangeBtnRow, "Auto (clear)");
        rangeBtnRow->AddFrame(clearRangeBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
        clearRangeBtn->Connect("Clicked()", "GammaFitGUI", this, "OnClearFitRange()");
        clearRangeBtn->SetToolTipText("Reset to automatic range: peaks ± Range×σ");
    }

    // ── Peak navigation ───────────────────────────────────────────────────────
    TGGroupFrame* navGrp = new TGGroupFrame(p, "Peak Navigation");
    p->AddFrame(navGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    {
        TGHorizontalFrame* navBtnRow = new TGHorizontalFrame(navGrp);
        navGrp->AddFrame(navBtnRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        TGTextButton* prevBtn = new TGTextButton(navBtnRow, "◀ Prev");
        navBtnRow->AddFrame(prevBtn, new TGLayoutHints(kLHintsLeft, 2, 4, 0, 0));
        prevBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPrevPeak()");
        prevBtn->SetToolTipText("Zoom to the previous cached peak");

        TGTextButton* nextBtn = new TGTextButton(navBtnRow, "Next ▶");
        navBtnRow->AddFrame(nextBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
        nextBtn->Connect("Clicked()", "GammaFitGUI", this, "OnNextPeak()");
        nextBtn->SetToolTipText("Zoom to the next cached peak");

        TGTextButton* zoomInBtn = new TGTextButton(navBtnRow, " + ");
        navBtnRow->AddFrame(zoomInBtn, new TGLayoutHints(kLHintsLeft, 8, 2, 0, 0));
        zoomInBtn->Connect("Clicked()", "GammaFitGUI", this, "OnZoomIn()");
        zoomInBtn->SetToolTipText("Narrow the x-axis view around the current peak");

        TGTextButton* zoomOutBtn = new TGTextButton(navBtnRow, " − ");
        navBtnRow->AddFrame(zoomOutBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
        zoomOutBtn->Connect("Clicked()", "GammaFitGUI", this, "OnZoomOut()");
        zoomOutBtn->SetToolTipText("Widen the x-axis view around the current peak");
    }
    peakNavLbl_ = new TGLabel(navGrp, "No peaks loaded");
    navGrp->AddFrame(peakNavLbl_, new TGLayoutHints(kLHintsLeft, 4, 4, 0, 2));

    TGTextButton* delEntryBtn = new TGTextButton(navGrp, "Remove Current Entry from Cache");
    navGrp->AddFrame(delEntryBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
    delEntryBtn->Connect("Clicked()", "GammaFitGUI", this, "OnDeleteCacheEntry()");
    delEntryBtn->SetToolTipText("Delete the currently displayed cache entry from the fit cache file");

    // ── Peaks ─────────────────────────────────────────────────────────────────
    TGGroupFrame* pkGrp = new TGGroupFrame(p, "Peaks");
    p->AddFrame(pkGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    TGHorizontalFrame* pkRow = new TGHorizontalFrame(pkGrp);
    pkGrp->AddFrame(pkRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    addPeakChk_ = new TGCheckButton(pkRow, "Accumulate peaks");
    addPeakChk_->SetState(kButtonDown);
    pkRow->AddFrame(addPeakChk_, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 8, 0, 0));
    addPeakChk_->SetToolTipText("When checked, each canvas click adds to the peak list.\nWhen unchecked, each click replaces the previous peak.");

    TGTextButton* clrBtn = new TGTextButton(pkRow, "Clear All");
    pkRow->AddFrame(clrBtn, new TGLayoutHints(kLHintsRight, 0, 2, 0, 0));
    clrBtn->Connect("Clicked()", "GammaFitGUI", this, "OnClearPeaks()");

    TGTextButton* rmBtn = new TGTextButton(pkRow, "Remove");
    pkRow->AddFrame(rmBtn, new TGLayoutHints(kLHintsRight, 0, 4, 0, 0));
    rmBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRemovePeak()");
    rmBtn->SetToolTipText("Remove the selected peak from the list");

    peakListBox_ = new TGListBox(pkGrp, 400);
    peakListBox_->Resize(285, 62);
    pkGrp->AddFrame(peakListBox_, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));

    // ── Peak Statistics ───────────────────────────────────────────────────────
    TGGroupFrame* statsGrp = new TGGroupFrame(p, "Peak Statistics");
    p->AddFrame(statsGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    peakStatsView_ = new TGTextView(statsGrp, 285, 110);
    statsGrp->AddFrame(peakStatsView_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    // ── Residuals ─────────────────────────────────────────────────────────────
    TGGroupFrame* resGrp = new TGGroupFrame(p, "Residuals");
    p->AddFrame(resGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    residualChk_ = new TGCheckButton(resGrp, "Show residuals  (data-fit)/sigma");
    resGrp->AddFrame(residualChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
    residualChk_->Connect("Clicked()", "GammaFitGUI", this, "OnToggleResiduals()");
    residualChk_->SetToolTipText("Split the canvas and show (data-fit)/sigma below the spectrum");

    {
        TGHorizontalFrame* resRow = new TGHorizontalFrame(resGrp);
        resGrp->AddFrame(resRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
        resRow->AddFrame(new TGLabel(resRow, "Fit:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        residualCombo_ = new TGComboBox(resRow, 600);
        residualCombo_->Resize(210, 22);
        resRow->AddFrame(residualCombo_, new TGLayoutHints(kLHintsExpandX));
        residualCombo_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                "OnSelectResidualFit(Int_t)");
    }

    // ── Fit options ───────────────────────────────────────────────────────────
    TGGroupFrame* fitOptGrp = new TGGroupFrame(p, "Fit Options  (default: chi2 / least squares)");
    p->AddFrame(fitOptGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    mFitLogLikChk_ = new TGCheckButton(fitOptGrp, "Log-likelihood  (L)");
    fitOptGrp->AddFrame(mFitLogLikChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 0));
    mFitLogLikChk_->SetToolTipText("Use Poisson log-likelihood instead of chi2 (better for low-count bins)");

    mFitImprovChk_ = new TGCheckButton(fitOptGrp, "IMPROVE  (M)");
    fitOptGrp->AddFrame(mFitImprovChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));
    mFitImprovChk_->SetToolTipText("Run IMPROVE after MIGRAD to search for a better minimum");

    mFitMinosChk_ = new TGCheckButton(fitOptGrp, "MINOS errors  (E)");
    fitOptGrp->AddFrame(mFitMinosChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 4));
    mFitMinosChk_->SetToolTipText("Compute asymmetric MINOS errors (slower, more accurate near parameter boundaries)");

    // ── Fit actions ───────────────────────────────────────────────────────────
    TGGroupFrame* ag = new TGGroupFrame(p, "Actions");
    p->AddFrame(ag, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    TGHorizontalFrame* btnRow1 = new TGHorizontalFrame(ag);
    ag->AddFrame(btnRow1, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    TGTextButton* prevBtn = new TGTextButton(btnRow1, "  Preview  ");
    btnRow1->AddFrame(prevBtn, new TGLayoutHints(kLHintsLeft, 2, 6, 0, 0));
    prevBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPreview()");
    prevBtn->SetToolTipText("Draw a Gaussian at the current parameters without running MIGRAD");

    TGTextButton* fitBtn = new TGTextButton(btnRow1, "  Run Fit  ");
    btnRow1->AddFrame(fitBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
    fitBtn->Connect("Clicked()", "GammaFitGUI", this, "OnManualFit()");
    fitBtn->SetToolTipText("Run MIGRAD on the fit window defined by Energy ± Range×Sigma");

    mResultLbl_ = new TGLabel(ag, "No fit yet");
    ag->AddFrame(mResultLbl_, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

    mFitParamsView_ = new TGTextView(ag, 290, 90);
    ag->AddFrame(mFitParamsView_, new TGLayoutHints(kLHintsExpandX, 4, 4, 0, 4));

    // ── Peak label + classification (saved with cache entry) ──────────────────
    {
        TGHorizontalFrame* lblRow = new TGHorizontalFrame(ag);
        ag->AddFrame(lblRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 0));
        lblRow->AddFrame(new TGLabel(lblRow, "Label:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        mPeakLabelCombo_ = new TGComboBox(lblRow, 921);
        mPeakLabelCombo_->AddEntry("(none)", 1);
        mPeakLabelCombo_->Select(1, kFALSE);
        mPeakLabelCombo_->Resize(200, 22);
        lblRow->AddFrame(mPeakLabelCombo_, new TGLayoutHints(kLHintsExpandX));
    }
    {
        TGHorizontalFrame* clsRow = new TGHorizontalFrame(ag);
        ag->AddFrame(clsRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 0));
        clsRow->AddFrame(new TGLabel(clsRow, "Class:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        mPeakClass_ = new TGComboBox(clsRow, 920);
        mPeakClass_->AddEntry("(none)",               1);
        mPeakClass_->AddEntry("Parent",               2);
        mPeakClass_->AddEntry("Daughter",             3);
        mPeakClass_->AddEntry("Granddaughter",        4);
        mPeakClass_->AddEntry("Beta-n Daughter",      5);
        mPeakClass_->AddEntry("Beta-2n Daughter",     6);
        mPeakClass_->AddEntry("Beta-n Granddaughter", 7);
        mPeakClass_->AddEntry("Beta-2n Granddaughter",8);
        mPeakClass_->AddEntry("Background",           9);
        mPeakClass_->AddEntry("Custom",               10);
        mPeakClass_->AddEntry("X-ray",                11);
        mPeakClass_->Select(1, kFALSE);
        mPeakClass_->Resize(150, 22);
        clsRow->AddFrame(mPeakClass_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        mPeakCustom_ = new TGTextEntry(clsRow, "");
        mPeakCustom_->SetToolTipText("Custom sub-classification name (used when Class = Custom)");
        clsRow->AddFrame(mPeakCustom_, new TGLayoutHints(kLHintsExpandX));
    }

    TGHorizontalFrame* btnRow2 = new TGHorizontalFrame(ag);
    ag->AddFrame(btnRow2, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 4));

    TGTextButton* accBtn = new TGTextButton(btnRow2, "Accept & Save to Cache");
    btnRow2->AddFrame(accBtn, new TGLayoutHints(kLHintsExpandX, 2, 6, 0, 0));
    accBtn->Connect("Clicked()", "GammaFitGUI", this, "OnAcceptFit()");
    accBtn->SetToolTipText("Write the current fit parameters to this histogram's cache file");

    TGTextButton* rejBtn = new TGTextButton(btnRow2, "Reject");
    btnRow2->AddFrame(rejBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
    rejBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRejectFit()");
    rejBtn->SetToolTipText("Discard the current manual fit and redraw the histogram");

    TGTextButton* loadCacheBtn = new TGTextButton(ag, "Load Cache onto Histogram");
    ag->AddFrame(loadCacheBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 2));
    loadCacheBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLoadCache()");
    loadCacheBtn->SetToolTipText("Reload the cache from disk and overlay all stored fits on the current histogram");

    TGTextButton* scanBtn = new TGTextButton(ag, "Parameter Scan Plot");
    ag->AddFrame(scanBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
    scanBtn->Connect("Clicked()", "GammaFitGUI", this, "OnParameterScan()");
    scanBtn->SetToolTipText("Open a popup canvas: chi2/ndf vs each fit parameter scanned around the best-fit value");
}

// ─────────────────────────────────────────────────────────────────────────────
// File management
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnOpenFile()
{
    static const char* kTypes[] = {
        "ROOT files", "*.root",
        "All files",  "*",
        nullptr, nullptr
    };
    TGFileInfo fi;
    fi.fFileTypes = kTypes;
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if (!fi.fFilename) return;

    if (inputFile_) { inputFile_->Close(); delete inputFile_; inputFile_ = nullptr; }

    inputPath_ = fi.fFilename;
    inputFile_ = TFile::Open(inputPath_.c_str(), "READ");
    if (!inputFile_ || inputFile_->IsZombie()) {
        AppendLog("ERROR: Cannot open " + inputPath_);
        inputFile_ = nullptr;
        return;
    }

    histNames_.clear();
    th2Names_.clear();
    projParent_.clear();
    histClass_.clear();
    customProjDefs_.clear();
    TIter next(inputFile_->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        if (!obj) continue;
        std::string name = obj->GetName();
        if (obj->InheritsFrom("TH2")) {
            histNames_.push_back(name);
            th2Names_.insert(name);
            // Virtual projection entries — created on demand when selected
            histNames_.push_back(name + "_px");
            histNames_.push_back(name + "_py");
            projParent_[name + "_px"] = name;
            projParent_[name + "_py"] = name;
        } else if (obj->InheritsFrom("TH1")) {
            histNames_.push_back(name);
        }
    }

    LoadMetadata();
    PopulateHistWidgets();
    PopulateCustProjTh2Combo();

    // Pre-populate Fit Results list with any histogram that already has a cache
    fittedHists_.clear();
    fitResultsList_->RemoveAll();
    for (const auto& hname : histNames_) {
        struct stat st;
        if (stat(CacheFileFor(hname).c_str(), &st) == 0) {
            fittedHists_.push_back(hname);
            fitResultsList_->AddEntry(hname.c_str(), (Int_t)fittedHists_.size());
        }
    }
    fitResultsList_->MapSubwindows(); fitResultsList_->Layout();

    std::string display = inputPath_;
    if (display.size() > 38)
        display = "..." + display.substr(display.size() - 35);
    fileLbl_->SetText(display.c_str());

    AppendLog("Opened: " + inputPath_ +
              "  (" + std::to_string(histNames_.size()) + " histograms)");
    SetStatus("File: " + inputPath_);
}

void GammaFitGUI::OnOpenIsotopeDB()
{
    static const char* kTypes[] = {
        "Text files", "*.txt",
        "All files",  "*",
        nullptr, nullptr
    };
    TGFileInfo fi;
    fi.fFileTypes = kTypes;
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if (!fi.fFilename) return;

    isotopePath_ = fi.fFilename;
    dbLoaded_ = db_.Load(isotopePath_);
    if (dbLoaded_) {
        AppendLog("Isotope DB loaded: " + isotopePath_ +
                  "  (" + std::to_string(db_.db.size()) + " lines)");
        if (isoDbSearch_) isoDbSearch_->SetText("");
        PopulateIsoDbList();
    } else {
        AppendLog("WARNING: Could not load isotope DB from " +
                  isotopePath_ + " — matches will be empty");
    }

    std::string display = isotopePath_;
    if (display.size() > 38)
        display = "..." + display.substr(display.size() - 35);
    isotopeLbl_->SetText(display.c_str());
}

void GammaFitGUI::PopulateHistWidgets()
{
    histList_->RemoveAll();
    manualCombo_->RemoveAll();
    fwhmCombo_->RemoveAll();
    for (size_t i = 0; i < histNames_.size(); ++i) {
        const std::string& name = histNames_[i];
        std::string display = name;
        bool isTH2   = th2Names_.count(name) > 0;
        bool isProj  = projParent_.find(name) != projParent_.end();
        bool isCProj = customProjDefs_.find(name) != customProjDefs_.end();
        if (isTH2) {
            display = "[2D] " + name;
        } else if (isCProj) {
            display = "[Custom] " + name;
        } else if (isProj) {
            auto& pp = projParent_.at(name);
            bool isX = name.size() >= 3 && name.substr(name.size()-3) == "_px";
            display = (isX ? "[ProjX] " : "[ProjY] ") + pp;
        } else {
            std::string cls = ClassOf(name);
            if (cls == "Decay") display = "[Decay] " + name;
        }
        histList_->AddEntry(display.c_str(), (Int_t)i + 1);
        if (!isTH2) {
            manualCombo_->AddEntry(display.c_str(), (Int_t)i + 1);
            fwhmCombo_->AddEntry(display.c_str(), (Int_t)i + 1);
        }
    }
    histList_->MapSubwindows();    histList_->Layout();
    manualCombo_->MapSubwindows(); manualCombo_->Layout();
    fwhmCombo_->MapSubwindows();   fwhmCombo_->Layout();
    PopulateDecayTh2Combo();
}

void GammaFitGUI::OnHistogramSelected(Int_t id)
{
    if (!inputFile_ || id < 1 || (size_t)id > histNames_.size()) return;
    currentHist_ = histNames_[id - 1];

    if (rawHistOwned_) { delete rawHist_; rawHist_ = nullptr; rawHistOwned_ = false; }
    rawHist_ = LoadHistFromFile(currentHist_, rawHistOwned_);

    // Sync classification combo
    if (histClassCombo_) {
        std::string cls = ClassOf(currentHist_);
        if      (cls == "Decay") histClassCombo_->Select(2, kFALSE);
        else if (cls == "2D")    histClassCombo_->Select(3, kFALSE);
        else                     histClassCombo_->Select(1, kFALSE);
    }

    DrawOnCanvas(rawHist_);
    SetStatus("Selected: " + currentHist_);
}

void GammaFitGUI::OnLoadManual()
{
    Int_t id = manualCombo_->GetSelected();
    if (id < 1) { AppendLog("Select a histogram from the dropdown first."); return; }
    OnHistogramSelected(id);

    manualPeaks_.clear();
    peakListBox_->RemoveAll();
    peakListBox_->MapSubwindows(); peakListBox_->Layout();
    delete manualTF1_; manualTF1_ = nullptr;
    delete bgTF1_;     bgTF1_     = nullptr;
    for (TF1* fc : fitComponents_) delete fc;
    fitComponents_.clear();
    bgClickMode_    = false;
    bgClickCount_   = 0;
    rangeClickMode_  = false;
    rangeClickCount_ = 0;
    mFitLo_->SetNumber(0.0);
    mFitHi_->SetNumber(0.0);
    mResultLbl_->SetText("No fit yet");
    viewXmin_ = 0.0;
    viewXmax_ = 0.0;
    peakNavIdx_ = 0;

    // Apply the selected view mode (raw / bg-sub / bg-sub + peaks)
    OnHistViewChanged(histViewCombo_->GetSelected());
    AppendLog("Loaded: " + currentHist_);

    PopulateNavAndResidual();
}

void GammaFitGUI::OnLoadCache()
{
    if (!rawHist_ || currentHist_.empty()) {
        AppendLog("Load a histogram first.");
        return;
    }

    // Rebuild bg-sub if the current view mode calls for it
    int viewMode = histViewCombo_ ? histViewCombo_->GetSelected() : 1;
    delete viewHist_; viewHist_ = nullptr;
    if (viewMode == 2 || viewMode == 3) {
        int iters = bgIterEntry_ ? (int)bgIterEntry_->GetNumber() : 14;
        viewHist_ = MakeBgSubHist(rawHist_, true, iters);
    }
    TH1* h = viewHist_ ? viewHist_ : rawHist_;

    TCanvas* c = canvas_->GetCanvas();
    c->Clear();
    c->cd();
    ApplyHistStyle(h, currentHist_.c_str());
    if (!currentHist_.empty() && ClassOf(currentHist_) == "Decay")
        h->GetXaxis()->SetTitle("Time (ms)");
    h->SetLineColor(kBlack);
    h->SetMarkerSize(0);
    if (viewXmin_ < viewXmax_)
        h->GetXaxis()->SetRangeUser(viewXmin_, viewXmax_);
    else
        h->GetXaxis()->UnZoom();
    SetYMaxFromVisible(h);
    h->Draw("hist");
    h->Draw("E1 same");

    // Always overlay all cached fits regardless of view mode
    OverlayFitPeaks(currentHist_, c);

    // Keep any current manual TF1 visible on top, with components
    if (manualTF1_) {
        manualTF1_->SetLineColor(kRed);
        manualTF1_->SetLineWidth(2);
        DrawFitComponents(c, manualTF1_);  // bg (green dashed) + per-Gaussian (blue dashed)
        manualTF1_->Draw("same");          // total fit on top in red
        DrawPeakLabels(manualTF1_);
    }

    c->Modified(); c->Update();

    // Refresh nav + residual list so newly accepted manual fits appear
    PopulateNavAndResidual();

    // If residuals are on, replace the plain draw with a split residual canvas
    if (residualsOn_) {
        TF1* f = nullptr;
        if (manualTF1_) {
            f = manualTF1_;
        } else {
            Int_t sel = residualCombo_->GetSelected();
            if (sel >= 1 && sel <= (Int_t)peakNavKeys_.size())
                f = BuildFromCacheKey(peakNavKeys_[sel - 1]);
            else if (!peakNavKeys_.empty())
                f = BuildFromCacheKey(peakNavKeys_[0]);
        }
        if (f) {
            DrawWithResiduals(h, f, f->GetXmin(), f->GetXmax());
            if (f != manualTF1_) delete f;
        }
    }

    // Ensure currentHist_ appears in the Fit Results list
    auto it = std::find(fittedHists_.begin(), fittedHists_.end(), currentHist_);
    if (it == fittedHists_.end()) {
        fittedHists_.push_back(currentHist_);
        fitResultsList_->AddEntry(currentHist_.c_str(), (Int_t)fittedHists_.size());
        fitResultsList_->MapSubwindows(); fitResultsList_->Layout();
    }

    AppendLog("Cache loaded: " + currentHist_);
    SetStatus("Cache loaded: " + currentHist_);
}

// ─────────────────────────────────────────────────────────────────────────────
// AutoFit
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::SyncDebugToggles()
{
    for (int i = 0; i < 8; ++i)
        Debug::Set(static_cast<Debug::Section>(i), debugChk_[i]->IsOn());
    Debug::PrintConfig();
}

void GammaFitGUI::OnDebugAllOn()
{
    for (int i = 0; i < 8; ++i) debugChk_[i]->SetState(kButtonDown);
    Debug::SetAll(true);
    AppendLog("Debug: ALL sections enabled");
}

void GammaFitGUI::OnDebugAllOff()
{
    for (int i = 0; i < 8; ++i) debugChk_[i]->SetState(kButtonUp);
    Debug::SetAll(false);
    AppendLog("Debug: all sections disabled");
}

std::string GammaFitGUI::CacheFileFor(const std::string& hname) const
{
    return std::string(kCacheDir) + "/fit_cache_" + hname + ".dat";
}

TH1* GammaFitGUI::LoadProjection(TFile* f,
                                  const std::string& projName,
                                  const std::map<std::string, std::string>& projParent,
                                  const TGTextEntry* xLblEntry,
                                  const TGTextEntry* yLblEntry) const
{
    auto it = projParent.find(projName);
    if (it == projParent.end() || !f) return nullptr;

    TH2* h2 = (TH2*)f->Get(it->second.c_str());
    if (!h2) return nullptr;

    bool isX = projName.size() >= 3 && projName.substr(projName.size()-3) == "_px";
    TH1* proj = isX ? h2->ProjectionX(projName.c_str())
                    : h2->ProjectionY(projName.c_str());
    if (!proj) return nullptr;

    proj->SetDirectory(nullptr);

    // Apply axis labels: the projected axis keeps its TH2 axis label
    std::string xLbl = isX
        ? (xLblEntry ? std::string(xLblEntry->GetText()) : "X")
        : (yLblEntry ? std::string(yLblEntry->GetText()) : "Y");
    if (!xLbl.empty()) proj->GetXaxis()->SetTitle(xLbl.c_str());
    proj->GetYaxis()->SetTitle("Counts");

    return proj;
}

TH1* GammaFitGUI::LoadHistFromFile(const std::string& hname, bool& owned) const
{
    owned = false;
    if (!inputFile_) return nullptr;

    // Custom projection with user-defined cut range
    auto cit = customProjDefs_.find(hname);
    if (cit != customProjDefs_.end()) {
        const CustomProjDef& def = cit->second;
        TH2* h2 = (TH2*)inputFile_->Get(def.th2Name.c_str());
        if (!h2) return nullptr;
        TH1* proj = nullptr;
        if (def.projX) {
            int b1 = h2->GetYaxis()->FindBin(def.lo);
            int b2 = h2->GetYaxis()->FindBin(def.hi);
            proj = h2->ProjectionX(hname.c_str(), b1, b2);
        } else {
            int b1 = h2->GetXaxis()->FindBin(def.lo);
            int b2 = h2->GetXaxis()->FindBin(def.hi);
            proj = h2->ProjectionY(hname.c_str(), b1, b2);
        }
        if (proj) { proj->SetDirectory(nullptr); owned = true; }
        return proj;
    }

    // Auto projection (full axis)
    if (projParent_.find(hname) != projParent_.end()) {
        TH1* h = LoadProjection(inputFile_, hname, projParent_, nullptr, nullptr);
        owned = (h != nullptr);
        return h;
    }

    return (TH1*)inputFile_->Get(hname.c_str());
}

void GammaFitGUI::OnApplyTh2Labels()
{
    if (!rawHist_) return;
    if (rawHist_->InheritsFrom("TH2")) {
        if (th2XLabelEntry_) rawHist_->GetXaxis()->SetTitle(th2XLabelEntry_->GetText());
        if (th2YLabelEntry_) rawHist_->GetYaxis()->SetTitle(th2YLabelEntry_->GetText());
        DrawOnCanvas(rawHist_);
        return;
    }
    // Projection — recreate with updated label
    if (projParent_.count(currentHist_)) {
        if (rawHistOwned_) { delete rawHist_; rawHist_ = nullptr; rawHistOwned_ = false; }
        rawHist_ = LoadProjection(inputFile_, currentHist_, projParent_,
                                  th2XLabelEntry_, th2YLabelEntry_);
        rawHistOwned_ = (rawHist_ != nullptr);
    } else {
        if (th2XLabelEntry_ && rawHist_) rawHist_->GetXaxis()->SetTitle(th2XLabelEntry_->GetText());
        if (th2YLabelEntry_ && rawHist_) rawHist_->GetYaxis()->SetTitle(th2YLabelEntry_->GetText());
    }
    DrawOnCanvas(rawHist_);
}

void GammaFitGUI::OnApplySrcTh2Labels()
{
    if (!rawHist_) return;
    if (rawHist_->InheritsFrom("TH2")) {
        if (srcTh2XLabelEntry_) rawHist_->GetXaxis()->SetTitle(srcTh2XLabelEntry_->GetText());
        if (srcTh2YLabelEntry_) rawHist_->GetYaxis()->SetTitle(srcTh2YLabelEntry_->GetText());
        DrawOnCanvas(rawHist_);
        return;
    }
    if (srcProjParent_.count(currentHist_)) {
        if (rawHistOwned_) { delete rawHist_; rawHist_ = nullptr; rawHistOwned_ = false; }
        rawHist_ = LoadProjection(srcRootFile_, currentHist_, srcProjParent_,
                                  srcTh2XLabelEntry_, srcTh2YLabelEntry_);
        rawHistOwned_ = (rawHist_ != nullptr);
    }
    DrawOnCanvas(rawHist_);
}

void GammaFitGUI::RunFitOnHistogram(const std::string& hname,
                                    TFile* overrideFile,
                                    const std::vector<double>& forcedSeeds)
{
    TFile*       srcFile = overrideFile ? overrideFile : inputFile_;
    const std::string& srcPath = overrideFile ? srcRootPath_ : inputPath_;

    if (!srcFile) { AppendLog("No ROOT file loaded."); return; }

    currentHist_ = hname;
    if (rawHistOwned_) { delete rawHist_; rawHist_ = nullptr; rawHistOwned_ = false; }
    if (!overrideFile) {
        rawHist_ = LoadHistFromFile(hname, rawHistOwned_);
    } else {
        auto sit = srcProjParent_.find(hname);
        if (sit != srcProjParent_.end()) {
            rawHist_ = LoadProjection(overrideFile, hname, srcProjParent_, nullptr, nullptr);
            rawHistOwned_ = (rawHist_ != nullptr);
        } else {
            rawHist_ = (TH1*)overrideFile->Get(hname.c_str());
        }
    }

    SyncDebugToggles();
    mkdir(kCacheDir, 0755);

    bool useSeeds = useSeedsChk_->IsOn();

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(hname));

    if (!fitdb.rootFile.empty() && fitdb.rootFile != srcPath)
        AppendLog("WARNING: cache for " + hname + " was built from a different file");

    // Apply cached resolution model so grouping + seeding use up-to-date σ(E)
    {
        auto it = fitdb.GetEntries().find(kResolutionKey);
        if (it != fitdb.GetEntries().end() && it->second.params.size() == 3) {
            res_.a = it->second.params[0];
            res_.b = it->second.params[1];
            res_.c = it->second.params[2];
            AppendLog("Resolution model from cache: a=" + Fmt(res_.a, 4) +
                      "  b=" + Fmt(res_.b, 4) +
                      "  c=" + Fmt(res_.c, 8));
        }
    }

    fitdb.useCachedSeeds = useSeeds;

    // Re-create tracker and storage per run so residuals don't accumulate
    PeakTracker runTracker;
    FitStorage  runStorage;
    PeakFitter  fitter(db_, &runTracker, res_, runStorage, &fitdb);

    bool hOwned = false;
    TH1* h = nullptr;
    if (!overrideFile) {
        h = LoadHistFromFile(hname, hOwned);
    } else {
        auto sit = srcProjParent_.find(hname);
        if (sit != srcProjParent_.end()) {
            h = LoadProjection(overrideFile, hname, srcProjParent_, nullptr, nullptr);
            hOwned = (h != nullptr);
        } else {
            h = (TH1*)overrideFile->Get(hname.c_str());
        }
    }
    if (!h) { AppendLog("ERROR: histogram '" + hname + "' not found"); return; }

    const std::string outPath = "AllGammaFits.root";
    TFile* fout = TFile::Open(outPath.c_str(), "UPDATE");
    if (!fout || fout->IsZombie()) {
        delete fout;
        fout = TFile::Open(outPath.c_str(), "RECREATE");
    }

    PeakFitter::BgOptions bgOpts;
    if (overrideFile) {
        // Source tab: use source-specific BG controls
        bgOpts.subtractBg = srcBgSubChk_ && srcBgSubChk_->IsOn();
        bgOpts.iterations = srcBgIterEntry_ ? (int)srcBgIterEntry_->GetNumber() : 14;
        AppendLog("Running Source AutoFit: " + hname +
                  (forcedSeeds.empty() ? "" :
                   "  [" + std::to_string(forcedSeeds.size()) + " seeds]") +
                  (bgOpts.subtractBg ? "  [bg sub]" : ""));
    } else {
        bgOpts.subtractBg       = bgSubtractChk_->IsOn();
        bgOpts.iterations       = (int)bgIterEntry_->GetNumber();
        bgOpts.tspecSigma       = tspecSigmaEntry_  ? tspecSigmaEntry_->GetNumber()  : 2.0;
        bgOpts.tspecThresh      = tspecThreshEntry_ ? tspecThreshEntry_->GetNumber() : 0.02;
        bgOpts.useLogLikelihood = autoLogLikChk_ ? autoLogLikChk_->IsOn() : true;
        bgOpts.useImprove       = autoImprovChk_ ? autoImprovChk_->IsOn() : false;
        AppendLog("Running AutoFit: " + hname +
                  (useSeeds ? "  [seeds on]" : "  [fresh start]") +
                  "  sigma=" + Fmt(bgOpts.tspecSigma, 1) +
                  "  thresh=" + Fmt(bgOpts.tspecThresh, 3));
    }
    SetStatus("Fitting: " + hname + " ...");

    fitter.FitHistogram(h, fout, true, nullptr, bgOpts, forcedSeeds);
    if (hOwned) { delete h; h = nullptr; }

    fitdb.bgSubtracted = bgOpts.subtractBg;
    fitdb.bgIterations = bgOpts.iterations;
    fitdb.rootFile     = srcPath;
    fitdb.Save(CacheFileFor(hname));
    fout->Close();
    delete fout;

    PopulateNavAndResidual();

    bool alreadyListed = false;
    for (const auto& n : fittedHists_)
        if (n == hname) { alreadyListed = true; break; }
    if (!alreadyListed) {
        fittedHists_.push_back(hname);
        fitResultsList_->AddEntry(hname.c_str(), (Int_t)fittedHists_.size());
        fitResultsList_->MapSubwindows(); fitResultsList_->Layout();
    }
    ShowFitResult(hname);

    if (!overrideFile) {
        for (size_t i = 0; i < histNames_.size(); i++) {
            if (histNames_[i] == hname) {
                manualCombo_->Select((Int_t)i + 1, kFALSE);
                break;
            }
        }
    }

    AppendLog("Done: " + hname);
    SetStatus("Finished: " + hname);
}

void GammaFitGUI::OnRunSelected()
{
    Int_t id = histList_->GetSelected();
    if (id < 1 || (size_t)id > histNames_.size()) {
        AppendLog("Select a histogram from the list first.");
        return;
    }
    RunFitOnHistogram(histNames_[id - 1]);
}

void GammaFitGUI::OnRunAll()
{
    if (histNames_.empty()) {
        AppendLog("No histograms loaded — open a ROOT file first.");
        return;
    }
    int nGamma = 0;
    for (const auto& name : histNames_)
        if (ClassOf(name) == "Gamma") nGamma++;
    AppendLog("=== AutoFit ALL Gamma: " + std::to_string(nGamma) +
              " / " + std::to_string(histNames_.size()) + " histograms ===");
    for (const auto& name : histNames_) {
        if (ClassOf(name) != "Gamma") continue;
        RunFitOnHistogram(name);
    }
    AppendLog("=== All Gamma histograms complete ===");
}

void GammaFitGUI::OnLoadCacheSelected()
{
    if (!inputFile_) { AppendLog("No ROOT file loaded."); return; }
    Int_t id = histList_->GetSelected();
    if (id < 1 || (size_t)id > histNames_.size()) {
        AppendLog("Select a histogram from the list first.");
        return;
    }
    const std::string& hname = histNames_[id - 1];
    currentHist_ = hname;
    if (rawHistOwned_) { delete rawHist_; rawHist_ = nullptr; rawHistOwned_ = false; }
    rawHist_ = LoadHistFromFile(hname, rawHistOwned_);
    if (!rawHist_) { AppendLog("ERROR: histogram '" + hname + "' not found"); return; }

    // Read bg settings from the cache and apply them to the view controls so
    // that OnLoadCache() displays the same histogram that was used during fitting.
    {
        FitDatabase fitdb;
        fitdb.Load(CacheFileFor(hname));
        if (fitdb.bgSubtracted) {
            // Mode 2 = background subtracted; fits are always overlaid by OnLoadCache
            histViewCombo_->Select(2, kFALSE);
            bgSubtractChk_->SetState(kButtonDown);
            bgIterEntry_->SetNumber(fitdb.bgIterations);
            AppendLog("Cache bg settings: subtracted=" +
                      std::to_string(fitdb.bgSubtracted) +
                      "  iterations=" + std::to_string(fitdb.bgIterations));
        } else {
            histViewCombo_->Select(1, kFALSE);
            bgSubtractChk_->SetState(kButtonUp);
        }
    }

    // Sync Manual Fit combo so the histogram is selected there too
    for (size_t i = 0; i < histNames_.size(); i++) {
        if (histNames_[i] == hname) {
            manualCombo_->Select((Int_t)i + 1, kFALSE);
            break;
        }
    }
    OnLoadCache();
}

void GammaFitGUI::OnLoadCacheAll()
{
    if (!inputFile_) { AppendLog("No ROOT file loaded."); return; }
    if (histNames_.empty()) { AppendLog("No histograms loaded."); return; }

    AppendLog("=== Load Cache ALL: " + std::to_string(histNames_.size()) + " histograms ===");
    int loaded = 0;
    for (const auto& hname : histNames_) {
        std::string cacheFile = CacheFileFor(hname);
        std::ifstream test(cacheFile);
        if (!test.is_open()) continue;   // no cache for this histogram
        test.close();

        if (th2Names_.count(hname)) continue;   // TH2 parents are never fitted directly

        // Register in Fit Results list (avoid duplicates)
        auto it = std::find(fittedHists_.begin(), fittedHists_.end(), hname);
        if (it == fittedHists_.end()) {
            fittedHists_.push_back(hname);
            fitResultsList_->AddEntry(hname.c_str(), (Int_t)fittedHists_.size());
        }
        ++loaded;
    }
    fitResultsList_->MapSubwindows();
    fitResultsList_->Layout();

    // Show the last histogram that was processed (if any)
    if (loaded > 0) {
        // point currentHist_ to the last selected / first available
        if (currentHist_.empty()) {
            currentHist_ = fittedHists_.front();
            if (rawHistOwned_) { delete rawHist_; rawHistOwned_ = false; }
            rawHist_ = LoadHistFromFile(currentHist_, rawHistOwned_);
        }
        AppendLog("Cache loaded for " + std::to_string(loaded) + " histograms.");
        SetStatus("Cache loaded: " + std::to_string(loaded) + " histograms.");
    } else {
        AppendLog("No cache files found for any loaded histogram.");
        SetStatus("No cache files found.");
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Manual Fit
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnCanvasEvent(Int_t event, Int_t px, Int_t py, TObject* /*obj*/)
{
    if (event != kButton1Down) return;

    // ── FWHM point toggle mode ────────────────────────────────────────────────
    if (fwhmRemoveModeChk_ && fwhmRemoveModeChk_->IsOn() && !fwhmAllX_.empty()) {
        TCanvas* c = canvas_->GetCanvas();
        c->cd();
        double xClick = c->PadtoX(c->AbsPixeltoX(px));
        double yClick = c->PadtoY(c->AbsPixeltoY(py));
        double xRange = c->GetUxmax() - c->GetUxmin();
        double yRange = c->GetUymax() - c->GetUymin();
        if (xRange <= 0.0 || yRange <= 0.0) return;

        bool showSigma = fwhmShowSigmaChk_ && fwhmShowSigmaChk_->IsOn();
        bool showRes   = fwhmShowResChk_   && fwhmShowResChk_->IsOn();

        int    nearIdx  = -1;
        double nearDist = 0.06;  // 6% of axis range threshold
        for (size_t i = 0; i < fwhmAllX_.size(); i++) {
            double fwhm = fwhmAllY_[i];
            double E    = fwhmAllX_[i];
            double y = showRes   ? 100.0 * fwhm / E
                     : showSigma ? fwhm / 2.3548
                     : fwhm;
            double dx = (E - xClick) / xRange;
            double dy = (y  - yClick) / yRange;
            double dist = std::sqrt(dx*dx + dy*dy);
            if (dist < nearDist) { nearDist = dist; nearIdx = (int)i; }
        }
        if (nearIdx >= 0) {
            fwhmExcluded_[nearIdx] = !fwhmExcluded_[nearIdx];
            AppendLog(std::string(fwhmExcluded_[nearIdx] ? "Excluded" : "Restored") +
                      " FWHM point: E=" + Fmt(fwhmAllX_[nearIdx], 1) +
                      " keV  FWHM=" + Fmt(fwhmAllY_[nearIdx], 3) + " keV");
            RedrawFWHM();

            // Auto-persist exclusion list to cache so it survives re-load
            if (!fwhmHistName_.empty()) {
                FitDatabase edb;
                edb.Load(CacheFileFor(fwhmHistName_));
                FitEntry eExcl;
                eExcl.key = kExcludedFwhmKey;
                for (size_t ei = 0; ei < fwhmExcluded_.size(); ei++)
                    if (fwhmExcluded_[ei]) eExcl.params.push_back(fwhmAllX_[ei]);
                edb.ForceStore(kExcludedFwhmKey, eExcl);
                mkdir(kCacheDir, 0755);
                edb.Save(CacheFileFor(fwhmHistName_));
            }
        }
        return;
    }

    if (!rawHist_) { AppendLog("Load a histogram first (Manual Fit tab → Load to Canvas)."); return; }

    TCanvas* c = canvas_->GetCanvas();

    // When residuals are on the canvas is split into padTop (spectrum) and
    // padBot (pulls).  All coordinate reads and marker draws must target
    // padTop; otherwise coordinates are wrong and lines are invisible.
    TPad* targetPad = (TPad*)c->FindObject("padTop");
    if (!targetPad) targetPad = (TPad*)c;
    targetPad->cd();

    double energy = targetPad->PadtoX(targetPad->AbsPixeltoX(px));
    double xmin   = rawHist_->GetXaxis()->GetXmin();
    double xmax   = rawHist_->GetXaxis()->GetXmax();
    if (energy < xmin || energy > xmax) return;

    double ylo = targetPad->GetUymin();
    double yhi = targetPad->GetUymax() * 0.95;

    // ── Fit range click mode (magenta) ────────────────────────────────────────
    if (rangeClickMode_) {
        TLine* line = new TLine(energy, ylo, energy, yhi);
        line->SetLineColor(kMagenta + 1);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");
        c->Modified(); c->Update();

        if (rangeClickCount_ == 0) {
            mFitLo_->SetNumber(energy);
            rangeClickCount_++;
            AppendLog("Fit range low: " + Fmt(energy, 2) + " keV — now click high edge");
        } else {
            double lo = mFitLo_->GetNumber();
            if (energy < lo) {
                mFitHi_->SetNumber(lo);
                mFitLo_->SetNumber(energy);
            } else {
                mFitHi_->SetNumber(energy);
            }
            rangeClickMode_  = false;
            rangeClickCount_ = 0;
            AppendLog("Fit range set: [" + Fmt(mFitLo_->GetNumber(), 2) +
                      ", " + Fmt(mFitHi_->GetNumber(), 2) + "] keV");
        }
        return;
    }

    // ── Background region click mode ──────────────────────────────────────────
    if (bgClickMode_) {
        TLine* line = new TLine(energy, ylo, energy, yhi);
        line->SetLineColor(kGreen + 2);
        line->SetLineStyle(2);
        line->SetLineWidth(2);
        line->Draw("same");
        c->Modified(); c->Update();

        if (bgClickCount_ == 0) {
            mBgLo_->SetNumber(energy);
            bgClickCount_++;
            AppendLog("BG low edge: " + Fmt(energy, 2) + " keV — now click high edge");
        } else {
            double lo = mBgLo_->GetNumber();
            if (energy < lo) {
                mBgHi_->SetNumber(lo);
                mBgLo_->SetNumber(energy);
            } else {
                mBgHi_->SetNumber(energy);
            }
            bgClickMode_  = false;
            bgClickCount_ = 0;
            AppendLog("BG region set: [" + Fmt(mBgLo_->GetNumber(), 2) +
                      ", " + Fmt(mBgHi_->GetNumber(), 2) + "] keV");
        }
        return;
    }

    // ── Peak placement ────────────────────────────────────────────────────────
    clickEnergy_ = energy;
    mEnergy_->SetNumber(energy);

    if (!addPeakChk_->IsOn()) {
        manualPeaks_.clear();
        peakListBox_->RemoveAll();
    }

    manualPeaks_.push_back(energy);
    peakListBox_->AddEntry(Form("%.2f keV", energy), (Int_t)manualPeaks_.size());
    peakListBox_->MapSubwindows(); peakListBox_->Layout();

    TLine* line = new TLine(energy, ylo, energy, yhi);
    line->SetLineColor(kBlue + 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");
    c->Modified(); c->Update();

    OnSeedParams();

    if (addPeakChk_->IsOn())
        AppendLog("Peak added: " + Fmt(energy, 2) + " keV  (total: " +
                  std::to_string(manualPeaks_.size()) + ")  — add more or click Run Fit");
    else
        AppendLog("Peak: " + Fmt(energy, 2) + " keV  — adjust parameters then click Run Fit");
}

void GammaFitGUI::OnSeedParams()
{
    double E = mEnergy_->GetNumber();
    if (E <= 0) return;

    mSigma_->SetNumber(res_.Sigma(E));

    if (rawHist_) {
        double amp = rawHist_->GetBinContent(rawHist_->FindBin(E));
        mAmp_->SetNumber(std::max(amp, 1.0));
    }
    mBg0_->SetNumber(0.0);
    mBg1_->SetNumber(0.0);
}

void GammaFitGUI::OnPreview()
{
    if (!rawHist_) { AppendLog("Load a histogram first."); return; }
    if (manualPeaks_.empty()) { AppendLog("Click on at least one peak first."); return; }

    int    n    = (int)manualPeaks_.size();
    double nSig = mRange_->GetNumber();
    double bg0  = mBg0_->GetNumber();
    double bg1  = mBg1_->GetNumber();

    double fitLo = mFitLo_->GetNumber();
    double fitHi = mFitHi_->GetNumber();
    double xmin, xmax;
    if (fitLo > 0 && fitHi > fitLo) {
        xmin = fitLo;
        xmax = fitHi;
    } else {
        xmin = manualPeaks_.front() - nSig * res_.Sigma(manualPeaks_.front());
        xmax = manualPeaks_.back()  + nSig * res_.Sigma(manualPeaks_.back());
    }

    delete manualTF1_;
    for (TF1* fc : fitComponents_) delete fc;
    fitComponents_.clear();
    manualTF1_ = new TF1("manual_preview", BuildNGaussFormula(n).c_str(), xmin, xmax);

    for (int i = 0; i < n; i++) {
        double E   = manualPeaks_[i];
        double sig = res_.Sigma(E);
        double A   = std::max(rawHist_->GetBinContent(rawHist_->FindBin(E)), 1.0);
        manualTF1_->SetParameter(3*i,   A);
        manualTF1_->SetParameter(3*i+1, E);
        manualTF1_->SetParameter(3*i+2, sig);
    }
    manualTF1_->SetParameter(3*n,   bg0);
    manualTF1_->SetParameter(3*n+1, bg1);
    manualTF1_->SetLineColor(kOrange + 1);
    manualTF1_->SetLineStyle(2);
    manualTF1_->SetLineWidth(2);

    TCanvas* c = canvas_->GetCanvas();
    c->cd();
    manualTF1_->Draw("same");
    c->Modified(); c->Update();

    AppendLog("Preview: " + std::to_string(n) + " peak(s)  window ["
              + Fmt(xmin,1) + ", " + Fmt(xmax,1) + "] keV");
}

void GammaFitGUI::OnManualFit()
{
    if (!rawHist_) { AppendLog("Load a histogram first."); return; }
    if (manualPeaks_.empty()) { AppendLog("Click on at least one peak first."); return; }

    int    n    = (int)manualPeaks_.size();
    double nSig = mRange_->GetNumber();
    double bg0  = mBg0_->GetNumber();
    double bg1  = mBg1_->GetNumber();

    double fitLo = mFitLo_->GetNumber();
    double fitHi = mFitHi_->GetNumber();
    double xmin, xmax;
    if (fitLo > 0 && fitHi > fitLo) {
        xmin = fitLo;
        xmax = fitHi;
    } else {
        xmin = manualPeaks_.front() - nSig * res_.Sigma(manualPeaks_.front());
        xmax = manualPeaks_.back()  + nSig * res_.Sigma(manualPeaks_.back());
    }

    // Determine which histogram to fit: bg-subtracted when view mode 2 or 3.
    int viewMode = histViewCombo_ ? histViewCombo_->GetSelected() : 1;
    if (viewMode == 2 || viewMode == 3) {
        if (!viewHist_) {
            int iters = bgIterEntry_ ? (int)bgIterEntry_->GetNumber() : 14;
            viewHist_ = MakeBgSubHist(rawHist_, true, iters);
        }
    } else {
        delete viewHist_; viewHist_ = nullptr;
    }
    TH1* fitHist = viewHist_ ? viewHist_ : rawHist_;

    delete manualTF1_;
    for (TF1* fc : fitComponents_) delete fc;
    fitComponents_.clear();
    manualTF1_ = new TF1("manual_fit", BuildNGaussFormula(n).c_str(), xmin, xmax);

    for (int i = 0; i < n; i++) {
        double E        = manualPeaks_[i];
        double sigModel = res_.Sigma(E);   // resolution-model prediction
        double sig      = sigModel;
        double A        = std::max(fitHist->GetBinContent(fitHist->FindBin(E)), 1.0);
        manualTF1_->SetParName(3*i,   Form("A_%d",   i+1));
        manualTF1_->SetParName(3*i+1, Form("E_%d",   i+1));
        manualTF1_->SetParName(3*i+2, Form("sig_%d", i+1));
        manualTF1_->SetParameter(3*i,     A);
        manualTF1_->SetParLimits(3*i,     std::max(A * 0.01, 1.0), A * 20.0);
        manualTF1_->SetParameter(3*i+1,   E);
        manualTF1_->SetParLimits(3*i+1,   E - 8.0, E + 8.0);
        manualTF1_->SetParameter(3*i+2,   sig);
        // For multi-peak fits keep sigmas tightly bound to the resolution model
        // so neighbouring Gaussians don't absorb each other.  Single-peak fits
        // allow a wider range in case the model is slightly off.
        if (n > 1)
            manualTF1_->SetParLimits(3*i+2, sigModel * 0.5, sigModel * 2.0);
        else
            manualTF1_->SetParLimits(3*i+2, sigModel * 0.2, sigModel * 4.0);
    }
    manualTF1_->SetParName(3*n,   "bg0");
    manualTF1_->SetParName(3*n+1, "bg1");
    manualTF1_->SetParameter(3*n,   bg0);
    manualTF1_->SetParameter(3*n+1, bg1);

    // Build fit option string from checkboxes
    std::string fitOpts = "R S Q B 0";
    if (mFitLogLikChk_ && mFitLogLikChk_->IsOn()) fitOpts += " L";
    if (mFitImprovChk_ && mFitImprovChk_->IsOn()) fitOpts += " M";
    if (mFitMinosChk_  && mFitMinosChk_ ->IsOn()) fitOpts += " E";

    AppendLog("Running manual fit: " + std::to_string(n) + " peak(s)"
              + "  window [" + Fmt(xmin,1) + ", " + Fmt(xmax,1) + "] keV"
              + "  opts=" + fitOpts
              + (viewHist_ ? "  [bg-sub]" : "  [raw]"));

    // "0" prevents ROOT from drawing/updating the canvas internally during the fit.
    // We redraw explicitly afterward so the zoom and view mode are preserved.
    TFitResultPtr r = fitHist->Fit(manualTF1_, fitOpts.c_str());

    // Update shared parameter displays
    mBg0_->SetNumber(manualTF1_->GetParameter(3*n));
    mBg1_->SetNumber(manualTF1_->GetParameter(3*n+1));
    if (n == 1) {
        mEnergy_->SetNumber(manualTF1_->GetParameter(1));
        mSigma_ ->SetNumber(manualTF1_->GetParameter(2));
        mAmp_   ->SetNumber(manualTF1_->GetParameter(0));
    }

    // r->Ndf() is unreliable for log-likelihood fits (returns 0 when bg-sub
    // histograms have negative bins that the Poisson fitter skips).
    // Compute chi²/ndf directly from pull residuals — valid for any fit method.
    double chi2ndf = -1.0;
    {
        double chi2 = 0.0;
        int    cnt  = 0;
        int b1 = fitHist->FindBin(xmin), b2 = fitHist->FindBin(xmax);
        for (int b = b1; b <= b2; b++) {
            double err = fitHist->GetBinError(b);
            if (err <= 0.0) continue;
            double pull = (fitHist->GetBinContent(b)
                           - manualTF1_->Eval(fitHist->GetBinCenter(b))) / err;
            chi2 += pull * pull;
            ++cnt;
        }
        int npar = manualTF1_->GetNpar();
        if (cnt > npar)
            chi2ndf = chi2 / (cnt - npar);
        else if (cnt > 0)
            chi2ndf = chi2 / cnt;   // DOF <= 0: fall back to chi2/count (no DOF correction)
    }
    lastManualChi2ndf_ = chi2ndf;   // carry to OnAcceptFit
    int status = r.Get() ? r->Status() : -1;

    std::string result = std::to_string(n) + " peak(s):";
    for (int i = 0; i < n; i++) {
        double fE   = manualTF1_->GetParameter(3*i+1);
        double fSig = manualTF1_->GetParameter(3*i+2);
        result += "  [" + Fmt(fE, 2) + " keV  FWHM=" + Fmt(2.355*fSig, 2) + "]";
    }
    result += "  chi2/ndf=" + Fmt(chi2ndf, 2) + "  status=" + std::to_string(status);

    mResultLbl_->SetText(result.c_str());

    // Fill parameter uncertainty table
    if (mFitParamsView_) {
        mFitParamsView_->Clear();
        for (int i = 0; i < n; i++) {
            std::string hdr = (n > 1) ? "Peak " + std::to_string(i+1) + ":" : "Peak:";
            mFitParamsView_->AddLine(hdr.c_str());
            double A    = manualTF1_->GetParameter(3*i);
            double E    = manualTF1_->GetParameter(3*i+1);
            double sig  = manualTF1_->GetParameter(3*i+2);
            double Aerr = manualTF1_->GetParError(3*i);
            double Eerr = manualTF1_->GetParError(3*i+1);
            double serr = manualTF1_->GetParError(3*i+2);
            mFitParamsView_->AddLine(("  A   = " + Fmt(A,3) + " +/- " + Fmt(Aerr,3)).c_str());
            mFitParamsView_->AddLine(("  E   = " + Fmt(E,4) + " +/- " + Fmt(Eerr,4) + " keV").c_str());
            mFitParamsView_->AddLine(("  sig = " + Fmt(sig,4) + " +/- " + Fmt(serr,4) + " keV").c_str());
            mFitParamsView_->AddLine(("  FWHM= " + Fmt(2.355*sig,4) + " +/- " + Fmt(2.355*serr,4) + " keV").c_str());
        }
        double bg0    = manualTF1_->GetParameter(3*n);
        double bg1    = manualTF1_->GetParameter(3*n+1);
        double bg0err = manualTF1_->GetParError(3*n);
        double bg1err = manualTF1_->GetParError(3*n+1);
        mFitParamsView_->AddLine(("  bg0 = " + Fmt(bg0,3) + " +/- " + Fmt(bg0err,3)).c_str());
        mFitParamsView_->AddLine(("  bg1 = " + Fmt(bg1,5) + " +/- " + Fmt(bg1err,5)).c_str());
        mFitParamsView_->AddLine(("chi2/ndf = " + Fmt(chi2ndf,3) + "  status=" + std::to_string(status)).c_str());
        mFitParamsView_->MapSubwindows(); mFitParamsView_->Layout();
    }

    AppendLog("Manual fit: " + result);

    // Populate label combo with DB matches for each fitted peak
    if (mPeakLabelCombo_ && dbLoaded_) {
        mPeakLabelCombo_->RemoveAll();
        mPeakLabelCombo_->AddEntry("(none)", 1);
        int cid = 2;
        std::set<std::string> seen;
        for (int i = 0; i < n; i++) {
            double Efit = manualTF1_->GetParameter(3*i+1);
            double fwhm = res_.FWHM(Efit);
            auto matches = db_.Match(Efit, fwhm);
            for (const auto& m : matches) {
                if (seen.count(m.isotope)) continue;
                seen.insert(m.isotope);
                mPeakLabelCombo_->AddEntry(m.isotope.c_str(), cid++);
            }
        }
        mPeakLabelCombo_->Select(1, kFALSE);
        mPeakLabelCombo_->MapSubwindows(); mPeakLabelCombo_->Layout();
        // Auto-select X-ray class if all peaks are < 100 keV
        if (mPeakClass_) {
            bool allXray = true;
            for (int i = 0; i < n; i++)
                if (manualTF1_->GetParameter(3*i+1) >= 100.0) { allXray = false; break; }
            if (allXray) mPeakClass_->Select(11, kFALSE);
        }
    }

    manualTF1_->SetLineColor(kRed);
    manualTF1_->SetLineWidth(2);

    // Update peak statistics display
    {
        TH1* dispHist = viewHist_ ? viewHist_ : rawHist_;
        UpdatePeakStats(manualTF1_, dispHist,
                        manualTF1_->GetXmin(), manualTF1_->GetXmax());
    }

    // Add current fit to residual combo (entry N+1)
    PopulateNavAndResidual();

    if (residualsOn_) {
        TH1* dispHist = viewHist_ ? viewHist_ : rawHist_;
        DrawWithResiduals(dispHist, manualTF1_,
                          manualTF1_->GetXmin(), manualTF1_->GetXmax());
    } else {
        // Redraw via the view-mode handler so zoom and histogram version are correct.
        // OnHistViewChanged draws manualTF1_ automatically when it is non-null.
        OnHistViewChanged(histViewCombo_->GetSelected());
    }

    SetStatus("Manual fit done  chi2/ndf=" + Fmt(chi2ndf, 2));
}

// ─────────────────────────────────────────────────────────────────────────────
// OnParameterScan  — popup chi2/ndf vs each fit parameter
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnParameterScan()
{
    if (!manualTF1_) { AppendLog("Run a manual fit first."); return; }
    if (!rawHist_)   { AppendLog("No histogram loaded."); return; }

    TH1*   fitHist = viewHist_ ? viewHist_ : rawHist_;
    int    npar    = manualTF1_->GetNpar();
    double xmin    = manualTF1_->GetXmin();
    double xmax    = manualTF1_->GetXmax();

    // Reference chi2/ndf at best-fit (pull-based, matches OnManualFit)
    double chi2best = 0.0;
    int    refCnt   = 0;
    {
        int b1 = fitHist->FindBin(xmin), b2 = fitHist->FindBin(xmax);
        for (int b = b1; b <= b2; b++) {
            double err = fitHist->GetBinError(b);
            if (err <= 0.0) continue;
            double pull = (fitHist->GetBinContent(b)
                           - manualTF1_->Eval(fitHist->GetBinCenter(b))) / err;
            chi2best += pull * pull;
            ++refCnt;
        }
        if (refCnt > npar) chi2best /= (refCnt - npar);
    }

    // Save best-fit values and errors
    std::vector<double> bestVals(npar), bestErrs(npar);
    for (int p = 0; p < npar; p++) {
        bestVals[p] = manualTF1_->GetParameter(p);
        bestErrs[p] = manualTF1_->GetParError(p);
    }

    // Grid layout: up to 3 columns
    int ncols = std::min(npar, 3);
    int nrows = (npar + ncols - 1) / ncols;
    int cw    = std::min(ncols * 380, 1200);
    int ch    = std::max(nrows * 300, 300);

    TCanvas* scanC = new TCanvas("param_scan",
        Form("Parameter Scan — %s", currentHist_.c_str()), cw, ch);
    scanC->Divide(ncols, nrows);

    static const int kSteps = 41;

    for (int p = 0; p < npar; p++) {
        double val0  = bestVals[p];
        double err0  = bestErrs[p];
        double delta = (err0 > 0.0) ? 3.0 * err0 : std::max(std::abs(val0) * 0.3, 1e-3);
        double lo    = val0 - delta;
        double hi    = val0 + delta;

        std::vector<double> xs(kSteps), ys(kSteps);
        for (int s = 0; s < kSteps; s++) {
            double v = lo + (hi - lo) * s / (kSteps - 1);
            xs[s] = v;

            // Hold all parameters at best-fit, fix parameter p to v
            for (int q = 0; q < npar; q++)
                manualTF1_->SetParameter(q, bestVals[q]);
            manualTF1_->FixParameter(p, v);

            double chi2 = 0.0;
            int    cnt2 = 0;
            int b1 = fitHist->FindBin(xmin), b2 = fitHist->FindBin(xmax);
            for (int b = b1; b <= b2; b++) {
                double e = fitHist->GetBinError(b);
                if (e <= 0.0) continue;
                double pull = (fitHist->GetBinContent(b)
                               - manualTF1_->Eval(fitHist->GetBinCenter(b))) / e;
                chi2 += pull * pull;
                ++cnt2;
            }
            ys[s] = (cnt2 > npar) ? chi2 / (cnt2 - npar) : chi2;
        }

        // Restore parameter to free state at best-fit value
        manualTF1_->ReleaseParameter(p);
        manualTF1_->SetParameter(p, val0);

        scanC->cd(p + 1);
        gPad->SetTopMargin(0.14);
        gPad->SetBottomMargin(0.16);
        gPad->SetLeftMargin(0.18);
        gPad->SetRightMargin(0.05);

        TGraph* gr = new TGraph(kSteps, xs.data(), ys.data());
        const char* pname = manualTF1_->GetParName(p);
        gr->SetTitle(Form("%s;%s;#chi^{2}/ndf", pname, pname));
        gr->SetLineColor(kBlue + 1);
        gr->SetLineWidth(2);
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.5);
        gr->Draw("AL");

        // Ensure y range includes chi2best+1 and the scan range
        double ymin = gr->GetHistogram()->GetMinimum();
        double ymax = std::max(gr->GetHistogram()->GetMaximum(),
                               chi2best + 1.5);
        gr->GetHistogram()->SetMinimum(ymin);
        gr->GetHistogram()->SetMaximum(ymax);
        gr->GetHistogram()->GetXaxis()->SetTitleSize(0.06);
        gr->GetHistogram()->GetXaxis()->SetLabelSize(0.05);
        gr->GetHistogram()->GetYaxis()->SetTitleSize(0.06);
        gr->GetHistogram()->GetYaxis()->SetLabelSize(0.05);
        gr->GetHistogram()->GetYaxis()->SetTitleOffset(1.3);
        gr->Draw("AL");

        // Vertical red dashed line at best-fit value
        TLine* vline = new TLine(val0, ymin, val0, ymax);
        vline->SetLineColor(kRed);
        vline->SetLineStyle(2);
        vline->SetLineWidth(2);
        vline->Draw();

        // Horizontal green dotted line at chi2best + 1 (1-sigma boundary)
        double oneSig = chi2best + 1.0;
        TLine* hline = new TLine(lo, oneSig, hi, oneSig);
        hline->SetLineColor(kGreen + 2);
        hline->SetLineStyle(3);
        hline->SetLineWidth(2);
        hline->Draw();
    }

    scanC->Modified();
    scanC->Update();
    AppendLog("Parameter scan: " + std::to_string(npar) +
              " params  chi2/ndf_best=" + Fmt(chi2best, 3));
}

void GammaFitGUI::OnAcceptFit()
{
    if (!manualTF1_) {
        AppendLog("No fit to save — run a manual fit first.");
        return;
    }
    if (currentHist_.empty()) {
        AppendLog("No histogram selected — load a histogram first.");
        return;
    }

    // Key is built from user-clicked manualPeaks_ so it matches the positions
    // TSpectrum would have found.  Fitted means may differ slightly but the
    // clicked seeds are what the rest of the cache is indexed by.
    int npar = manualTF1_->GetNpar();
    bool isDG = (npar == 7 && (int)manualPeaks_.size() <= 1);
    int nPeaks = isDG ? 1 : (int)manualPeaks_.size();
    if (!isDG && nPeaks == 0 && (npar - 2) % 3 == 0)
        nPeaks = (npar - 2) / 3;

    // Sort peaks ascending so MakeKey is consistent with AutoFit ordering
    std::vector<double> seedsForKey = manualPeaks_;
    std::sort(seedsForKey.begin(), seedsForKey.end());
    std::string key       = FitDatabase::MakeKey(seedsForKey);
    std::string cacheFile = CacheFileFor(currentHist_);

    FitDatabase fitdb;
    fitdb.Load(cacheFile);
    if (fitdb.GetEntries().empty()) {
        fitdb.bgSubtracted = false;
        fitdb.bgIterations = (int)bgIterEntry_->GetNumber();
    }

    // Remove all existing cache entries that overlap with any manual peak.
    // Tolerance = 8 keV so that both individual entries for a combined fit
    // and entries whose TSpectrum seed drifted slightly are all evicted.
    double tolKeV = 8.0;
    fitdb.RemoveOverlapping(seedsForKey, tolKeV);

    // Build the new entry with honest residual metrics so future seeding and
    // the OverlayFitPeaks quality check both see realistic numbers.
    FitEntry e;
    e.key = key;
    e.params.resize(npar);
    e.paramErrors.resize(npar);
    for (int i = 0; i < npar; ++i) {
        e.params[i]      = manualTF1_->GetParameter(i);
        e.paramErrors[i] = manualTF1_->GetParError(i);
    }

    double xlo = manualTF1_->GetXmin();
    double xhi = manualTF1_->GetXmax();
    TH1* fitHist = viewHist_ ? viewHist_ : rawHist_;
    auto rm = FitDatabase::ComputeResiduals(fitHist, manualTF1_, xlo, xhi);
    e.residualRMS = rm.rms;
    e.maxPull     = rm.maxPull;
    // Use the chi2/ndf that OnManualFit computed (chi2 / (N_bins - N_par)),
    // which matches what the user sees in mResultLbl_ and the residual plot.
    e.chi2ndf     = (lastManualChi2ndf_ > 0.0) ? lastManualChi2ndf_ : rm.rms * rm.rms;
    e.xlo         = xlo;
    e.xhi         = xhi;

    // Record which fit method was used
    {
        std::string meth;
        if (mFitLogLikChk_ && mFitLogLikChk_->IsOn()) meth += "L";
        if (mFitImprovChk_ && mFitImprovChk_->IsOn()) meth += "M";
        if (mFitMinosChk_  && mFitMinosChk_ ->IsOn()) meth += "E";
        e.fitMethod = meth;  // empty = default chi2 / least squares
    }

    // Record isotope label from combo
    if (mPeakLabelCombo_) {
        TGLBEntry* le = mPeakLabelCombo_->GetSelectedEntry();
        std::string lbl = le ? le->GetTitle() : "";
        if (lbl != "(none)" && !lbl.empty()) e.label = lbl;
    }
    // Record classification; auto-assign X-ray for unclassified peaks < 100 keV
    if (mPeakClass_) {
        int sel = mPeakClass_->GetSelected();
        std::string cust = mPeakCustom_ ? mPeakCustom_->GetText() : "";
        e.classification = ClassToString(sel, cust);
    }
    if (e.classification.empty()) {
        // Check if the primary fitted peak is below 100 keV
        int nPeaks2 = ((int)e.params.size() - 2) / 3;
        for (int pi = 0; pi < nPeaks2; pi++) {
            double Efit = e.params[3 * pi + 1];
            if (Efit > 0 && Efit < 100.0) { e.classification = "X-ray"; break; }
        }
    }

    // Force-store: manual fits always win regardless of score comparison.
    fitdb.ForceStore(key, e);
    fitdb.rootFile = inputPath_;
    mkdir(kCacheDir, 0755);
    fitdb.Save(cacheFile);

    // Append to the Gamma_fits text report
    {
        std::ofstream ftxt("../Gamma_fits/" + currentHist_ + "_fit.txt", std::ios::app);
        if (ftxt.is_open()) {
            ftxt << "\n=====================================\n";
            ftxt << "MANUAL FIT (saved to cache)\n";
            ftxt << "key: " << key << "\n";
            ftxt << "chi2/ndf (approx): " << e.chi2ndf << "\n";
            ftxt << "=====================================\n";
            for (int i = 0; i < nPeaks; i++) {
                double fE   = manualTF1_->GetParameter(isDG ? 1 : 3*i + 1);
                double fSig = manualTF1_->GetParameter(isDG ? 2 : 3*i + 2);
                double fA   = manualTF1_->GetParameter(isDG ? 0 : 3*i);
                ftxt << "Peak " << i
                     << ":  E="      << fE
                     << " keV  sigma=" << fSig
                     << "  FWHM="    << 2.355 * fSig
                     << "  A="       << fA << "\n";
            }
        }
    }

    AppendLog("Cache saved (manual override):"
              "  hist=" + currentHist_ +
              "  key="  + key +
              "  chi2/ndf≈" + Fmt(e.chi2ndf, 2));
    SetStatus("Saved: " + currentHist_ + "  key=" + key);
}

void GammaFitGUI::OnRejectFit()
{
    delete manualTF1_; manualTF1_ = nullptr;
    delete bgTF1_;     bgTF1_     = nullptr;
    for (TF1* fc : fitComponents_) delete fc;
    fitComponents_.clear();
    mResultLbl_->SetText("Fit rejected.");
    OnHistViewChanged(histViewCombo_->GetSelected());
    AppendLog("Manual fit rejected — canvas cleared.");
}

// ─────────────────────────────────────────────────────────────────────────────
// Multi-peak / background helpers
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnClearPeaks()
{
    manualPeaks_.clear();
    peakListBox_->RemoveAll();
    peakListBox_->MapSubwindows(); peakListBox_->Layout();
    AppendLog("Peak list cleared.");
}

void GammaFitGUI::OnRemovePeak()
{
    Int_t id = peakListBox_->GetSelected();
    if (id < 1 || (size_t)id > manualPeaks_.size()) {
        AppendLog("Select a peak in the list first.");
        return;
    }
    double removed = manualPeaks_[id - 1];
    manualPeaks_.erase(manualPeaks_.begin() + (id - 1));

    // Rebuild listbox with fresh sequential IDs
    peakListBox_->RemoveAll();
    for (size_t i = 0; i < manualPeaks_.size(); i++)
        peakListBox_->AddEntry(Form("%.2f keV", manualPeaks_[i]), (Int_t)i + 1);
    peakListBox_->MapSubwindows(); peakListBox_->Layout();

    AppendLog("Removed peak at " + Fmt(removed, 2) + " keV  (" +
              std::to_string(manualPeaks_.size()) + " remaining)");
}

void GammaFitGUI::OnSetBgFromCanvas()
{
    bgClickMode_  = true;
    bgClickCount_ = 0;
    AppendLog("Click once on the canvas for the background low edge, once for high edge.");
}

void GammaFitGUI::OnFitBackground()
{
    if (!rawHist_) { AppendLog("Load a histogram first."); return; }
    double lo = mBgLo_->GetNumber();
    double hi = mBgHi_->GetNumber();
    if (lo >= hi) { AppendLog("Set a valid background region (Lo < Hi) first."); return; }

    // Match the histogram used by OnManualFit: bg-sub when view mode 2 or 3
    int viewMode = histViewCombo_ ? histViewCombo_->GetSelected() : 1;
    if ((viewMode == 2 || viewMode == 3) && !viewHist_) {
        int iters = bgIterEntry_ ? (int)bgIterEntry_->GetNumber() : 14;
        viewHist_ = MakeBgSubHist(rawHist_, true, iters);
    } else if (viewMode == 1) {
        delete viewHist_; viewHist_ = nullptr;
    }
    TH1* fitHist = viewHist_ ? viewHist_ : rawHist_;

    delete bgTF1_;
    bgTF1_ = new TF1("bg_manual", "[0]+[1]*x", lo, hi);
    fitHist->Fit(bgTF1_, "R S Q B 0");

    bgTF1_->SetLineColor(kGreen + 2);
    bgTF1_->SetLineStyle(2);
    bgTF1_->SetLineWidth(2);

    mBg0_->SetNumber(bgTF1_->GetParameter(0));
    mBg1_->SetNumber(bgTF1_->GetParameter(1));

    // Redraw via view-mode handler to preserve zoom and histogram version,
    // then overlay the background fit line on top.
    OnHistViewChanged(histViewCombo_->GetSelected());
    TCanvas* c = canvas_->GetCanvas();
    c->cd();
    bgTF1_->Draw("same");
    c->Modified(); c->Update();

    AppendLog("Background fit in [" + Fmt(lo, 1) + ", " + Fmt(hi, 1) + "] keV:"
              + "  bg0=" + Fmt(bgTF1_->GetParameter(0), 2)
              + "  bg1=" + Fmt(bgTF1_->GetParameter(1), 5));
}

void GammaFitGUI::OnClearBackground()
{
    // Cancel any in-progress canvas click session
    bgClickMode_  = false;
    bgClickCount_ = 0;

    // Delete the background TF1 and redraw without it
    delete bgTF1_; bgTF1_ = nullptr;

    // Reset Lo/Hi fields
    if (mBgLo_) mBgLo_->SetNumber(0.0);
    if (mBgHi_) mBgHi_->SetNumber(0.0);

    // Redraw the histogram so the green background line disappears
    if (rawHist_)
        OnHistViewChanged(histViewCombo_ ? histViewCombo_->GetSelected() : 1);

    AppendLog("Background region cleared.");
}

// ─────────────────────────────────────────────────────────────────────────────
// Fit Results tab
// ─────────────────────────────────────────────────────────────────────────────
TH1* GammaFitGUI::MakeBgSubHist(TH1* raw, bool doSubtract, int iterations)
{
    TH1* h = (TH1*)raw->Clone(Form("%s_bgsub_gui", raw->GetName()));
    h->SetDirectory(0);
    h->Sumw2();
    if (doSubtract) {
        TSpectrum sp;
        TH1* bkg = sp.Background(h, iterations);
        h->Add(bkg, -1);
    }
    return h;
}

void GammaFitGUI::OverlayFitPeaks(const std::string& hname, TCanvas* c)
{
    FitDatabase fitdb;
    if (!fitdb.Load(CacheFileFor(hname))) return;

    if (!fitdb.rootFile.empty() && !inputPath_.empty() && fitdb.rootFile != inputPath_)
        AppendLog("WARNING: cache for " + hname + " was built from " +
                  fitdb.rootFile + " — current file is " + inputPath_);

    // Build candidates sorted by composite score (best first).
    // Greedy draw: skip any entry whose seeds are already claimed within ±3σ.
    struct Candidate {
        std::vector<double> peaks;  // seed energies parsed from cache key
        const FitEntry*     entry;
    };
    std::vector<Candidate> candidates;
    for (const auto& [key, entry] : fitdb.GetEntries()) {
        std::vector<double> peakE;
        std::istringstream ss(key);
        std::string tok;
        while (std::getline(ss, tok, '_')) {
            try { peakE.push_back(std::stod(tok)); } catch (...) {}
        }
        if (!peakE.empty()) candidates.push_back({peakE, &entry});
    }
    std::sort(candidates.begin(), candidates.end(),
        [](const Candidate& a, const Candidate& b){
            return FitDatabase::CompositeScore(*a.entry) <
                   FitDatabase::CompositeScore(*b.entry);
        });

    std::vector<double> claimed;
    static const int kFitColors[] = { kRed, kBlue+1, kGreen+2, kMagenta+1, kOrange+1 };
    int colorIdx = 0;

    TH1* dispH = viewHist_ ? viewHist_ : rawHist_;

    TLatex lbl;
    lbl.SetTextSize(0.022);
    lbl.SetTextAlign(21);

    // Global set of already-labeled energies — prevents two labels within 1 FWHM
    // of each other regardless of which cache group they come from.
    std::vector<double> labeledEs;

    auto labelTooClose = [&](double E) {
        double sep = res_.FWHM(E);
        for (double de : labeledEs)
            if (std::abs(E - de) < sep) return true;
        return false;
    };

    for (const auto& cand : candidates) {
        // Skip if any seed is within 3σ of an already-claimed seed
        bool overlap = false;
        for (double E : cand.peaks) {
            double tol = 3.0 * res_.Sigma(E);
            for (double drawn : claimed)
                if (std::abs(E - drawn) < tol) { overlap = true; break; }
            if (overlap) break;
        }
        if (overlap) continue;

        // Claim seeds now (regardless of fit quality) so later entries don't
        // produce duplicate labels for the same physical peak.
        for (double E : cand.peaks) claimed.push_back(E);

        double xlo, xhi;
        if (cand.entry->xlo < cand.entry->xhi) {
            xlo = cand.entry->xlo;
            xhi = cand.entry->xhi;
        } else {
            xlo = cand.peaks.front() - 5.0 * res_.Sigma(cand.peaks.front());
            xhi = cand.peaks.back()  + 5.0 * res_.Sigma(cand.peaks.back());
        }

        TF1* f = RebuildFromEntry(*cand.entry, xlo, xhi);
        if (f) {
            // Transfer ownership to the pad so the TF1 survives until the
            // canvas is cleared or written to file (needed for csave->Write()).
            f->SetBit(kCanDelete);
            f->SetLineColor(kFitColors[colorIdx++ % 5]);
            f->SetLineWidth(2);
            f->Draw("same");
        }

        // Collect label positions from Gaussian components (not seeds).
        // Iterating components avoids the seed→Gaussian mapping that would
        // assign the same Gaussian mean to multiple seeds in the same group
        // (e.g. double-Gaussian with a shared mean parameter).
        std::vector<double> compEs;
        if (f) {
            int npar = f->GetNpar();
            if (npar == 7) {
                compEs.push_back(f->GetParameter(1));
            } else if (npar >= 5 && (npar - 2) % 3 == 0) {
                int n = (npar - 2) / 3;
                for (int i = 0; i < n; i++)
                    compEs.push_back(f->GetParameter(3*i + 1));
            }
        } else {
            compEs = cand.peaks;   // no fit — fall back to seed energies
        }

        for (double E : compEs) {
            double yLabel = f ? f->Eval(E)
                              : (dispH ? dispH->GetBinContent(dispH->FindBin(E)) : 0.0);
            if (yLabel <= 0.0) continue;
            if (labelTooClose(E)) continue;   // another label is already within 1 FWHM
            labeledEs.push_back(E);

            std::string label = Form("%.1f", E);
            if (dbLoaded_) {
                auto matches = db_.Match(E, res_.FWHM(E));
                if (!matches.empty())
                    label += " (" + matches[0].isotope + ")";
            }
            lbl.DrawLatex(E, yLabel * 1.08, label.c_str());
        }
        // f is now owned by the pad (kCanDelete) — do NOT delete here
    }
    c->Modified(); c->Update();
}

// ─────────────────────────────────────────────────────────────────────────────
// Manual Fit — histogram view mode
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnHistViewChanged(Int_t id)
{
    if (!rawHist_) return;

    delete viewHist_; viewHist_ = nullptr;

    TCanvas* c = canvas_->GetCanvas();
    c->Clear();
    c->cd();

    TH1* h = rawHist_;

    if (id == 2 || id == 3) {
        int iters = bgIterEntry_ ? (int)bgIterEntry_->GetNumber() : 14;
        viewHist_ = MakeBgSubHist(rawHist_, true, iters);
        h = viewHist_;
    }

    ApplyHistStyle(h, currentHist_.empty() ? nullptr : currentHist_.c_str());
    if (!currentHist_.empty() && ClassOf(currentHist_) == "Decay")
        h->GetXaxis()->SetTitle("Time (ms)");
    h->SetLineColor(kBlack);
    h->SetMarkerSize(0);
    if (viewXmin_ < viewXmax_)
        h->GetXaxis()->SetRangeUser(viewXmin_, viewXmax_);
    else
        h->GetXaxis()->UnZoom();
    SetYMaxFromVisible(h);
    h->Draw("hist");
    h->Draw("E1 same");

    if (id == 3 && !currentHist_.empty())
        OverlayFitPeaks(currentHist_, c);

    if (manualTF1_) {
        manualTF1_->SetLineColor(kRed);
        manualTF1_->SetLineWidth(2);
        manualTF1_->Draw("same");
        DrawPeakLabels(manualTF1_);
    }

    c->Modified(); c->Update();
}

void GammaFitGUI::OnSetRangeFromCanvas()
{
    rangeClickMode_  = true;
    rangeClickCount_ = 0;
    AppendLog("Click once on the canvas for the fit low edge, once for the high edge. (magenta lines)");
}

void GammaFitGUI::OnClearFitRange()
{
    rangeClickMode_  = false;
    rangeClickCount_ = 0;
    mFitLo_->SetNumber(0.0);
    mFitHi_->SetNumber(0.0);
    AppendLog("Fit range cleared — will use peaks ± Range×σ.");
}

// ─────────────────────────────────────────────────────────────────────────────
// Peak navigation
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::PopulateNavAndResidual()
{
    peakNavKeys_.clear();
    FitDatabase fitdb;
    if (fitdb.Load(CacheFileFor(currentHist_))) {
        for (const auto& kv : fitdb.GetEntries()) {
            // Skip internal keys like __RESOLUTION__
            if (kv.first.size() >= 2 && kv.first[0] == '_' && kv.first[1] == '_') continue;
            peakNavKeys_.push_back(kv.first);
        }
        // Sort by first peak energy
        std::sort(peakNavKeys_.begin(), peakNavKeys_.end(),
            [](const std::string& a, const std::string& b){
                try { return std::stod(a) < std::stod(b); } catch (...) { return a < b; }
            });
    }

    // Update nav label
    if (peakNavKeys_.empty())
        peakNavLbl_->SetText("No cached peaks");
    else
        peakNavLbl_->SetText(Form("Peak 1/%d: %s keV",
                                   (int)peakNavKeys_.size(),
                                   peakNavKeys_[0].c_str()));

    // Rebuild residual combo, preserving the current selection
    Int_t prevSel = residualCombo_->GetSelected();
    Int_t maxId   = (Int_t)peakNavKeys_.size() + (manualTF1_ ? 1 : 0);

    residualCombo_->RemoveAll();
    for (size_t i = 0; i < peakNavKeys_.size(); i++)
        residualCombo_->AddEntry(("Cache: " + peakNavKeys_[i] + " keV").c_str(),
                                  (Int_t)i + 1);
    if (manualTF1_)
        residualCombo_->AddEntry("Current Manual Fit",
                                  (Int_t)peakNavKeys_.size() + 1);
    residualCombo_->MapSubwindows(); residualCombo_->Layout();

    // Restore previous selection if it still exists; fall back to 1
    if (prevSel > 0 && prevSel <= maxId)
        residualCombo_->Select(prevSel, kFALSE);
    else if (!peakNavKeys_.empty())
        residualCombo_->Select(1, kFALSE);
}

TF1* GammaFitGUI::RebuildFromEntry(const FitEntry& entry, double xlo, double xhi)
{
    int  npar = (int)entry.params.size();
    bool isDG = (npar == 7);
    TF1* f    = nullptr;

    // Use a unique name each time so ROOT's global function registry doesn't
    // delete a previous TF1 that the canvas is still holding a pointer to.
    static int sCtr = 0;
    std::string fname = "res_tf1_" + std::to_string(sCtr++);

    if (isDG) {
        f = new TF1(fname.c_str(),
            "[0]*exp(-0.5*((x-[1])/[2])^2)"
            "+[3]*exp(-0.5*((x-[1])/[4])^2)+[5]+[6]*x",
            xlo, xhi);
    } else if (npar >= 5 && (npar - 2) % 3 == 0) {
        int n = (npar - 2) / 3;
        f = new TF1(fname.c_str(), BuildNGaussFormula(n).c_str(), xlo, xhi);
    }
    if (!f) return nullptr;
    for (int i = 0; i < npar; i++)
        f->SetParameter(i, entry.params[i]);
    return f;
}

TF1* GammaFitGUI::BuildFromCacheKey(const std::string& key)
{
    FitDatabase fitdb;
    if (!fitdb.Load(CacheFileFor(currentHist_))) return nullptr;
    const auto& entries = fitdb.GetEntries();
    auto it = entries.find(key);
    if (it == entries.end()) return nullptr;

    // Parse peak energies from key to determine zoom window
    std::vector<double> peakE;
    std::istringstream ss(key);
    std::string tok;
    while (std::getline(ss, tok, '_')) {
        try { peakE.push_back(std::stod(tok)); } catch (...) {}
    }
    if (peakE.empty()) return nullptr;

    double xlo, xhi;
    if (it->second.xlo < it->second.xhi) {
        xlo = it->second.xlo;
        xhi = it->second.xhi;
    } else {
        xlo = peakE.front() - 6.0 * res_.Sigma(peakE.front());
        xhi = peakE.back()  + 6.0 * res_.Sigma(peakE.back());
        if (rawHist_) {
            xlo = std::max(xlo, rawHist_->GetXaxis()->GetXmin());
            xhi = std::min(xhi, rawHist_->GetXaxis()->GetXmax());
        }
    }
    return RebuildFromEntry(it->second, xlo, xhi);
}

void GammaFitGUI::NavigateToPeak(int idx)
{
    if (peakNavKeys_.empty()) {
        AppendLog("No cached peaks to navigate — run AutoFit or accept a manual fit first.");
        return;
    }
    peakNavIdx_ = std::max(0, std::min(idx, (int)peakNavKeys_.size() - 1));
    const std::string& key = peakNavKeys_[peakNavIdx_];

    // Load cache to get bg settings so we display the matching histogram version
    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(currentHist_));

    // Parse first peak energy for zoom
    double E = 0.0;
    try { E = std::stod(key.substr(0, key.find('_'))); } catch (...) {}

    if (E > 0.0 && rawHist_) {
        double sig = res_.Sigma(E);
        viewXmin_ = std::max(E - 8.0 * sig, rawHist_->GetXaxis()->GetXmin());
        viewXmax_ = std::min(E + 8.0 * sig, rawHist_->GetXaxis()->GetXmax());
    }

    // Show cached chi²/ndf in navigation label
    {
        const auto& entries = fitdb.GetEntries();
        auto it = entries.find(key);
        if (it != entries.end() && it->second.chi2ndf < 1.0e6) {
            peakNavLbl_->SetText(Form("Peak %d/%d: %s keV  chi2/ndf=%.2f",
                peakNavIdx_ + 1, (int)peakNavKeys_.size(),
                key.c_str(), it->second.chi2ndf));
        } else {
            peakNavLbl_->SetText(Form("Peak %d/%d: %s keV",
                peakNavIdx_ + 1, (int)peakNavKeys_.size(), key.c_str()));
        }
    }

    // Sync residual combo
    residualCombo_->Select(peakNavIdx_ + 1, kFALSE);

    // Reconstruct TF1 for this key
    TF1* f = BuildFromCacheKey(key);
    if (f) {
        f->SetLineColor(kRed);
        f->SetLineWidth(2);
    }

    // Respect the view combo — not the cache's bg flag — so the Manual Fit tab
    // always shows what the user selected.  Keep the result in viewHist_ so
    // the canvas never holds a dangling pointer when ROOT calls gPad->Update().
    int viewMode = histViewCombo_->GetSelected();
    delete viewHist_;
    if (viewMode == 2 || viewMode == 3) {
        int iters = bgIterEntry_ ? (int)bgIterEntry_->GetNumber() : 14;
        viewHist_ = MakeBgSubHist(rawHist_, true, iters);
    } else {
        viewHist_ = nullptr;
    }
    TH1* dispHist = viewHist_ ? viewHist_ : rawHist_;

    if (residualsOn_ && f) {
        DrawWithResiduals(dispHist, f, f->GetXmin(), f->GetXmax());
    } else {
        TCanvas* c = canvas_->GetCanvas();
        c->Clear();
        c->cd();
        ApplyHistStyle(dispHist, currentHist_.c_str());
        if (viewXmin_ < viewXmax_)
            dispHist->GetXaxis()->SetRangeUser(viewXmin_, viewXmax_);
        SetYMaxFromVisible(dispHist);
        dispHist->SetLineColor(kBlack);
        dispHist->SetMarkerSize(0);
        dispHist->Draw("hist");
        dispHist->Draw("E1 same");
        if (f) { f->Draw("same"); DrawPeakLabels(f); }
        c->Modified(); c->Update();
    }

    AppendLog("Peak " + std::to_string(peakNavIdx_+1) + "/" +
              std::to_string(peakNavKeys_.size()) + ": " + key + " keV");

    if (f) UpdatePeakStats(f, dispHist, f->GetXmin(), f->GetXmax());
    delete f;
}

void GammaFitGUI::OnPrevPeak()
{
    NavigateToPeak(peakNavIdx_ - 1);
}

void GammaFitGUI::OnNextPeak()
{
    NavigateToPeak(peakNavIdx_ + 1);
}

void GammaFitGUI::OnZoomIn()
{
    if (!rawHist_) return;

    // If no zoom is set yet, initialise from full histogram range first
    if (viewXmin_ >= viewXmax_) {
        viewXmin_ = rawHist_->GetXaxis()->GetXmin();
        viewXmax_ = rawHist_->GetXaxis()->GetXmax();
    }
    double centre = 0.5 * (viewXmin_ + viewXmax_);
    double half   = 0.5 * (viewXmax_ - viewXmin_) * 0.60; // shrink by 40%
    double minHalf = 2.0 * res_.Sigma(centre);             // floor: 2σ
    if (half < minHalf) half = minHalf;
    viewXmin_ = std::max(centre - half, rawHist_->GetXaxis()->GetXmin());
    viewXmax_ = std::min(centre + half, rawHist_->GetXaxis()->GetXmax());
    OnHistViewChanged(histViewCombo_->GetSelected());
}

void GammaFitGUI::OnZoomOut()
{
    if (!rawHist_) return;

    if (viewXmin_ >= viewXmax_) {
        // Already at full range — nothing to widen
        return;
    }
    double centre = 0.5 * (viewXmin_ + viewXmax_);
    double half   = 0.5 * (viewXmax_ - viewXmin_) * 1.70; // expand by 70%
    double xlo = centre - half;
    double xhi = centre + half;
    // If we've expanded to cover the full range, unzoom completely
    if (xlo <= rawHist_->GetXaxis()->GetXmin() &&
        xhi >= rawHist_->GetXaxis()->GetXmax()) {
        viewXmin_ = 0.0;
        viewXmax_ = 0.0;
    } else {
        viewXmin_ = std::max(xlo, rawHist_->GetXaxis()->GetXmin());
        viewXmax_ = std::min(xhi, rawHist_->GetXaxis()->GetXmax());
    }
    OnHistViewChanged(histViewCombo_->GetSelected());
}

void GammaFitGUI::OnDeleteCacheEntry()
{
    if (currentHist_.empty()) { AppendLog("Load a histogram first."); return; }
    if (peakNavKeys_.empty() || peakNavIdx_ < 0 ||
        peakNavIdx_ >= (int)peakNavKeys_.size()) {
        AppendLog("No cache entry selected — use Prev/Next to navigate to one.");
        return;
    }

    std::string key       = peakNavKeys_[peakNavIdx_];
    std::string cacheFile = CacheFileFor(currentHist_);

    FitDatabase fitdb;
    if (!fitdb.Load(cacheFile)) { AppendLog("No cache file found."); return; }

    if (!fitdb.Remove(key)) {
        AppendLog("Entry \"" + key + "\" not found in cache.");
        return;
    }
    fitdb.rootFile = inputPath_;
    fitdb.Save(cacheFile);

    AppendLog("Removed cache entry: " + key);
    SetStatus("Removed: " + key);

    // Refresh display — navigate to the previous entry if possible
    peakNavIdx_ = std::max(0, peakNavIdx_ - 1);
    PopulateNavAndResidual();
    OnLoadCache();
}

void GammaFitGUI::OnClearHistCache()
{
    if (currentHist_.empty()) { AppendLog("No histogram selected."); return; }

    std::string cacheFile = CacheFileFor(currentHist_);
    struct stat st;
    if (stat(cacheFile.c_str(), &st) != 0) {
        AppendLog("[ClearCache] No cache file for " + currentHist_);
        return;
    }

    Int_t ret = 0;
    new TGMsgBox(gClient->GetRoot(), this,
                 "Confirm Clear Cache",
                 ("Delete ALL cache entries for:\n" + currentHist_).c_str(),
                 kMBIconQuestion, kMBYes | kMBNo, &ret);
    if (ret != kMBYes) return;

    std::remove(cacheFile.c_str());

    // Remove from Fit Results list if present
    auto it = std::find(fittedHists_.begin(), fittedHists_.end(), currentHist_);
    if (it != fittedHists_.end()) {
        fittedHists_.erase(it);
        fitResultsList_->RemoveAll();
        for (size_t i = 0; i < fittedHists_.size(); i++)
            fitResultsList_->AddEntry(fittedHists_[i].c_str(), (Int_t)i + 1);
        fitResultsList_->MapSubwindows(); fitResultsList_->Layout();
    }

    peakNavKeys_.clear();
    peakNavIdx_ = 0;
    PopulateNavAndResidual();
    RedrawCurrent();

    AppendLog("[ClearCache] Cleared cache for " + currentHist_);
    SetStatus("Cache cleared: " + currentHist_);
}

// ─────────────────────────────────────────────────────────────────────────────
// Residuals
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::DrawWithResiduals(TH1* h, TF1* fit, double xlo, double xhi)
{
    TCanvas* c = canvas_->GetCanvas();
    c->Clear();

    // Top pad: spectrum + fit (70%)
    TPad* padTop = new TPad("padTop", "", 0.0, 0.30, 1.0, 1.0);
    padTop->SetBottomMargin(0.14);
    padTop->SetLeftMargin(0.12);
    padTop->Draw();
    padTop->cd();

    if (h) {
        if (viewXmin_ < viewXmax_)
            h->GetXaxis()->SetRangeUser(viewXmin_, viewXmax_);
        ApplyHistStyle(h, currentHist_.empty() ? nullptr : currentHist_.c_str());
        h->GetXaxis()->SetLabelSize(0.055);
        h->GetXaxis()->SetTitleSize(0.055);
        h->GetXaxis()->SetTitleOffset(0.85);
        h->GetYaxis()->SetLabelSize(0.05);
        h->GetYaxis()->SetTitleSize(0.05);
        h->SetLineColor(kBlack);
        h->SetMarkerSize(0);
        SetYMaxFromVisible(h);
        h->Draw("hist");
        h->Draw("E1 same");
        if (fit) {
            fit->SetLineColor(kRed);
            fit->SetLineWidth(2);
            fit->Draw("same");
            DrawPeakLabels(fit);
        }
    }

    // Bottom pad: (data-fit)/sigma (30%)
    c->cd();
    TPad* padBot = new TPad("padBot", "", 0.0, 0.0, 1.0, 0.30);
    padBot->SetTopMargin(0.02);
    padBot->SetBottomMargin(0.30);
    padBot->SetLeftMargin(0.12);
    padBot->Draw();
    padBot->cd();

    if (h && fit) {
        double lo = xlo;
        double hi = xhi;
        int binLo = h->FindBin(lo);
        int binHi = h->FindBin(hi);
        int nBins = binHi - binLo + 1;

        if (nBins > 0) {
            TH1D* hres = new TH1D("hres_display",
                                   ";Energy (keV);(data#minusfit)/#sigma",
                                   nBins,
                                   h->GetBinLowEdge(binLo),
                                   h->GetBinLowEdge(binHi + 1));
            for (int b = binLo; b <= binHi; b++) {
                double err  = h->GetBinError(b);
                double pull = (err > 0.0)
                    ? (h->GetBinContent(b) - fit->Eval(h->GetBinCenter(b))) / err
                    : 0.0;
                hres->SetBinContent(b - binLo + 1, pull);
            }
            hres->SetLineColor(kBlack);
            hres->SetFillColorAlpha(kAzure - 9, 0.6);
            hres->GetXaxis()->SetTitle("Energy (keV)");
            hres->GetYaxis()->SetTitle("(data - fit) / #sigma");
            hres->GetYaxis()->SetTitleSize(0.12);
            hres->GetYaxis()->SetTitleOffset(0.45);
            hres->GetYaxis()->SetLabelSize(0.10);
            hres->GetXaxis()->SetTitleSize(0.12);
            hres->GetXaxis()->SetTitleOffset(0.9);
            hres->GetXaxis()->SetLabelSize(0.10);
            hres->SetStats(0);
            hres->Draw("hist");

            // Reference lines: 0, ±1σ (blue dashed), ±3σ (orange dotted)
            double refLo = hres->GetXaxis()->GetXmin();
            double refHi = hres->GetXaxis()->GetXmax();
            static const double kLevels[] = { 0.0, 1.0, -1.0, 3.0, -3.0 };
            for (double lev : kLevels) {
                TLine* l = new TLine(refLo, lev, refHi, lev);
                if (lev == 0.0) {
                    l->SetLineColor(kRed);   l->SetLineWidth(2);
                } else if (std::abs(lev) < 2.0) {
                    l->SetLineColor(kBlue);  l->SetLineStyle(2);
                } else {
                    l->SetLineColor(kOrange+1); l->SetLineStyle(3);
                }
                l->Draw("same");
            }

            // Goodness-of-fit label in top-right corner of residual pad
            int    npar    = fit->GetNpar();
            int    ndf     = nBins - npar;
            double sumSq   = 0.0;
            double maxPull = 0.0;
            for (int b = 1; b <= nBins; b++) {
                double p = hres->GetBinContent(b);
                sumSq += p * p;
                if (std::abs(p) > maxPull) maxPull = std::abs(p);
            }
            double chi2ndf = (ndf > 0) ? sumSq / ndf
                           : (nBins > 0 ? sumSq / nBins : -1.0);
            const char* chi2tag = (ndf > 0) ? "#chi^{2}/ndf" : "#chi^{2}/N (low DOF)";

            TLatex gofLbl;
            gofLbl.SetNDC();
            gofLbl.SetTextSize(0.14);
            gofLbl.SetTextAlign(33); // top-right
            gofLbl.DrawLatex(0.90, 0.95,
                Form("%s = %.2f    max|pull| = %.1f", chi2tag, chi2ndf, maxPull));
        }
    }

    c->cd();
    c->Modified();
    c->Update();
}

void GammaFitGUI::RedrawCurrent()
{
    if (!rawHist_) return;
    TF1* fit = manualTF1_;  // prefer manual fit if it exists

    if (residualsOn_) {
        // Pick the TF1 from residual combo selection if no manual fit
        if (!fit && !peakNavKeys_.empty()) {
            Int_t sel = residualCombo_->GetSelected();
            if (sel >= 1 && sel <= (Int_t)peakNavKeys_.size())
                fit = BuildFromCacheKey(peakNavKeys_[sel - 1]);
        }
        if (fit) {
            DrawWithResiduals(rawHist_, fit, fit->GetXmin(), fit->GetXmax());
            if (fit != manualTF1_) delete fit;
            return;
        }
    }

    // Normal single-pad draw
    TCanvas* c = canvas_->GetCanvas();
    c->Clear();
    c->cd();
    if (viewXmin_ < viewXmax_)
        rawHist_->GetXaxis()->SetRangeUser(viewXmin_, viewXmax_);
    else
        rawHist_->GetXaxis()->UnZoom();
    rawHist_->SetLineColor(kBlack);
    rawHist_->SetMarkerSize(0);
    rawHist_->Draw("hist");
    rawHist_->Draw("E1 same");
    if (fit) {
        fit->SetLineColor(kRed);
        fit->SetLineWidth(2);
        fit->Draw("same");
    }
    c->Modified(); c->Update();
}

void GammaFitGUI::OnToggleResiduals()
{
    residualsOn_ = residualChk_->IsOn();
    if (residualsOn_) {
        AppendLog("Residuals ON — select a fit in the combo box to inspect.");
        Int_t sel = residualCombo_->GetSelected();
        if (sel >= 1)
            OnSelectResidualFit(sel);
        else
            OnHistViewChanged(histViewCombo_->GetSelected()); // show correct hist with no fit overlay
    } else {
        AppendLog("Residuals OFF.");
        OnHistViewChanged(histViewCombo_->GetSelected());
    }
}

void GammaFitGUI::OnSelectResidualFit(Int_t id)
{
    if (!rawHist_) return;

    TF1* f = nullptr;

    if (id >= 1 && id <= (Int_t)peakNavKeys_.size()) {
        f = BuildFromCacheKey(peakNavKeys_[id - 1]);

        // Zoom to this peak
        double E = 0.0;
        try { E = std::stod(peakNavKeys_[id-1].substr(0,
                            peakNavKeys_[id-1].find('_'))); } catch (...) {}
        if (E > 0.0) {
            double sig = res_.Sigma(E);
            viewXmin_ = std::max(E - 8.0*sig, rawHist_->GetXaxis()->GetXmin());
            viewXmax_ = std::min(E + 8.0*sig, rawHist_->GetXaxis()->GetXmax());
        }

        // Keep the display histogram alive in viewHist_ (view-combo governs).
        int viewMode = histViewCombo_->GetSelected();
        delete viewHist_;
        if (viewMode == 2 || viewMode == 3) {
            int iters = bgIterEntry_ ? (int)bgIterEntry_->GetNumber() : 14;
            viewHist_ = MakeBgSubHist(rawHist_, true, iters);
        } else {
            viewHist_ = nullptr;
        }
    } else if (manualTF1_) {
        f = manualTF1_;
    }

    if (!f) return;

    TH1* dispHist = viewHist_ ? viewHist_ : rawHist_;

    if (residualsOn_)
        DrawWithResiduals(dispHist, f, f->GetXmin(), f->GetXmax());

    if (f != manualTF1_) delete f;
}

// ─────────────────────────────────────────────────────────────────────────────
// Utilities
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::DrawFitComponents(TCanvas* /*c*/, TF1* f)
{
    // Delete components from the previous draw cycle
    for (TF1* c : fitComponents_) delete c;
    fitComponents_.clear();

    if (!f) return;
    int npar = f->GetNpar();
    if (npar < 5 || (npar - 2) % 3 != 0) return;  // not a standard N-Gaussian + linear BG
    int    n   = (npar - 2) / 3;
    double xlo = f->GetXmin();
    double xhi = f->GetXmax();
    double bg0 = f->GetParameter(3*n);
    double bg1 = f->GetParameter(3*n + 1);

    // Background: green dashed — always shown
    TF1* bgComp = new TF1(Form("fitcomp_bg_%d", (int)(size_t)f & 0xFFFF),
                          "[0]+[1]*x", xlo, xhi);
    bgComp->SetParameter(0, bg0);
    bgComp->SetParameter(1, bg1);
    bgComp->SetLineColor(kGreen + 2);
    bgComp->SetLineStyle(2);
    bgComp->SetLineWidth(2);
    fitComponents_.push_back(bgComp);
    bgComp->Draw("same");

    // Individual Gaussians + BG: blue dashed (only when n > 1)
    for (int i = 0; i < n && n > 1; i++) {
        double A   = f->GetParameter(3*i);
        double E   = f->GetParameter(3*i + 1);
        double sig = f->GetParameter(3*i + 2);
        TF1* gComp = new TF1(Form("fitcomp_g%d_%d", i, (int)(size_t)f & 0xFFFF),
                             "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x", xlo, xhi);
        gComp->SetParameter(0, A);
        gComp->SetParameter(1, E);
        gComp->SetParameter(2, sig);
        gComp->SetParameter(3, bg0);
        gComp->SetParameter(4, bg1);
        gComp->SetLineColor(kBlue + 1);
        gComp->SetLineStyle(2);
        gComp->SetLineWidth(1);
        fitComponents_.push_back(gComp);
        gComp->Draw("same");
    }
}

void GammaFitGUI::DrawPeakLabels(TF1* f)
{
    if (!f) return;
    int npar = f->GetNpar();

    std::vector<double> peakEs;
    bool isDG = (npar == 7);
    if (isDG) {
        peakEs.push_back(f->GetParameter(1));
    } else if (npar >= 5 && (npar - 2) % 3 == 0) {
        int n = (npar - 2) / 3;
        for (int i = 0; i < n; i++)
            peakEs.push_back(f->GetParameter(3*i + 1));
    } else {
        return;
    }

    // Sort by amplitude descending so the tallest peak wins when two labels
    // would land within 1 FWHM of each other.
    std::sort(peakEs.begin(), peakEs.end(), [&](double a, double b){
        return f->Eval(a) > f->Eval(b);
    });

    TLatex lbl;
    lbl.SetTextSize(0.022);
    lbl.SetTextAlign(21);

    std::vector<double> drawn;
    for (double E : peakEs) {
        double sep = res_.FWHM(E);
        bool tooClose = false;
        for (double de : drawn)
            if (std::abs(E - de) < sep) { tooClose = true; break; }
        if (tooClose) continue;
        drawn.push_back(E);

        double yTop = f->Eval(E);
        std::string label = Form("%.1f", E);
        if (dbLoaded_) {
            auto matches = db_.Match(E, res_.FWHM(E));
            if (!matches.empty())
                label += " (" + matches[0].isotope + ")";
        }
        lbl.DrawLatex(E, yTop * 1.08, label.c_str());
    }
}

void GammaFitGUI::UpdatePeakStats(TF1* f, TH1* h, double xlo, double xhi)
{
    if (!peakStatsView_) return;
    peakStatsView_->Clear();

    if (!f || !h) {
        peakStatsView_->MapSubwindows(); peakStatsView_->Layout();
        return;
    }

    int npar = f->GetNpar();
    int n    = (npar - 2) / 3;
    if (n < 1 || npar != 3*n + 2) {
        peakStatsView_->AddLine("Invalid fit model");
        peakStatsView_->MapSubwindows(); peakStatsView_->Layout();
        return;
    }

    // Recompute chi2/ndf and p-value from histogram residuals over the window
    double chi2 = 0.0;
    int    cnt  = 0;
    int b1 = h->FindBin(xlo), b2 = h->FindBin(xhi);
    for (int b = b1; b <= b2; b++) {
        double err = h->GetBinError(b);
        if (err <= 0.0) continue;
        double pull = (h->GetBinContent(b) - f->Eval(h->GetBinCenter(b))) / err;
        chi2 += pull * pull;
        ++cnt;
    }
    int    ndf     = cnt - npar;
    double chi2ndf = (ndf > 0) ? chi2 / ndf : (cnt > 0 ? chi2 / cnt : -1.0);
    double pval    = (ndf > 0) ? TMath::Prob(chi2, ndf) : -1.0;

    // Total raw counts in the fit window
    double totalCounts = 0.0;
    for (int b = b1; b <= b2; b++)
        totalCounts += h->GetBinContent(b);

    double bg0 = f->GetParameter(3*n);
    double bg1 = f->GetParameter(3*n+1);

    for (int i = 0; i < n; i++) {
        double A   = f->GetParameter(3*i);
        double E   = f->GetParameter(3*i+1);
        double sig = f->GetParameter(3*i+2);

        double fwhm   = 2.3548 * sig;
        double counts = A * sig * std::sqrt(2.0 * TMath::Pi());
        double bgUnder = std::abs((bg0 + bg1 * E) * 5.0 * sig);
        double snr = (bgUnder > 0) ? counts / std::sqrt(bgUnder) : 0.0;
        double frac = (totalCounts > 0) ? counts / totalCounts : 0.0;

        if (n > 1)
            peakStatsView_->AddLine(Form("--- Peak %d (%.2f keV) ---", i+1, E));
        peakStatsView_->AddLine(Form("  FWHM    = %.3f keV", fwhm));
        peakStatsView_->AddLine(Form("  Counts  = %.1f", counts));
        peakStatsView_->AddLine(Form("  SNR     = %.2f", snr));
        peakStatsView_->AddLine(Form("  Pk/Tot  = %.4f", frac));
    }

    if (ndf > 0) {
        peakStatsView_->AddLine(Form("chi2/ndf = %.3f  (%d dof)", chi2ndf, ndf));
        peakStatsView_->AddLine(Form("p-value  = %s",
            pval < 1e-4 ? Form("%.3e", pval) : Form("%.4f", pval)));
    }

    peakStatsView_->MapSubwindows();
    peakStatsView_->Layout();
}

void GammaFitGUI::DrawOnCanvas(TH1* h, TF1* fit)
{
    if (!h) return;
    TCanvas* c = canvas_->GetCanvas();
    c->cd();
    c->Clear();

    // TH2 — colour palette view (no fitting)
    if (h->InheritsFrom("TH2")) {
        h->SetTitle(h->GetName());
        if (th2XLabelEntry_ && std::string(th2XLabelEntry_->GetText()).size())
            h->GetXaxis()->SetTitle(th2XLabelEntry_->GetText());
        if (th2YLabelEntry_ && std::string(th2YLabelEntry_->GetText()).size())
            h->GetYaxis()->SetTitle(th2YLabelEntry_->GetText());
        h->Draw("COLZ");
        c->Modified(); c->Update();
        return;
    }

    ApplyHistStyle(h, h->GetName());
    if (!currentHist_.empty() && ClassOf(currentHist_) == "Decay")
        h->GetXaxis()->SetTitle("Time (ms)");
    h->SetLineColor(kBlack);
    h->SetMarkerSize(0);
    h->Draw("hist");
    h->Draw("E1 same");
    if (fit) {
        fit->SetLineColor(kRed);
        fit->SetLineWidth(2);
        fit->Draw("same");
    }
    c->Modified();
    c->Update();
}

void GammaFitGUI::OnSaveFWHMCanvas()
{
    if (fwhmAllX_.empty()) {
        AppendLog("No FWHM data loaded — load a histogram first.");
        return;
    }

    static const char* kTypes[] = {
        "ROOT files", "*.root",
        "All files",  "*",
        nullptr,      nullptr
    };
    TGFileInfo fi;
    fi.fFileTypes = kTypes;
    fi.fIniDir    = StrDup(".");
    new TGFileDialog(gClient->GetRoot(), this, kFDSave, &fi);
    if (!fi.fFilename) return;

    std::string fname = fi.fFilename;
    if (fname.size() < 5 ||
        fname.substr(fname.size() - 5) != ".root")
        fname += ".root";

    TFile* fout = TFile::Open(fname.c_str(), "UPDATE");
    if (!fout || fout->IsZombie()) {
        delete fout;
        fout = TFile::Open(fname.c_str(), "RECREATE");
    }
    if (!fout || fout->IsZombie()) {
        AppendLog("ERROR: cannot create " + fname);
        delete fout;
        return;
    }

    bool showSigma = fwhmShowSigmaChk_ && fwhmShowSigmaChk_->IsOn();
    bool showStat  = fwhmStatLineChk_  && fwhmStatLineChk_->IsOn();
    bool showRes   = fwhmShowResChk_   && fwhmShowResChk_->IsOn();

    Bool_t wasBatch = gROOT->IsBatch();
    gROOT->SetBatch(kTRUE);
    TCanvas* csave = new TCanvas(
        Form("c_fwhm_%s", fwhmHistName_.c_str()),
        ("FWHM: " + fwhmHistName_).c_str(), 800, 600);
    gROOT->SetBatch(wasBatch);
    csave->cd();
    DrawFWHMToCanvas(csave, showSigma, showStat, showRes);
    csave->Modified(); csave->Update();

    fout->cd();
    csave->Write();
    fout->Close(); delete fout;
    delete csave;

    AppendLog("FWHM canvas saved: " + fname);
    SetStatus("Saved: " + fname);
}

void GammaFitGUI::OnSavePlot()
{
    TCanvas* c = canvas_->GetCanvas();
    if (!c) { AppendLog("No canvas to save."); return; }

    static const char* kTypes[] = {
        "PDF files",  "*.pdf",
        "PNG files",  "*.png",
        "ROOT files", "*.root",
        "All files",  "*",
        nullptr,      nullptr
    };
    TGFileInfo fi;
    fi.fFileTypes = kTypes;
    fi.fIniDir    = StrDup(".");
    new TGFileDialog(gClient->GetRoot(), this, kFDSave, &fi);
    if (!fi.fFilename) return;

    std::string fname = fi.fFilename;
    // If no extension, default to PDF
    if (fname.find('.') == std::string::npos) fname += ".pdf";

    if (fname.size() >= 5 &&
        fname.substr(fname.size() - 5) == ".root")
    {
        TFile* fout = TFile::Open(fname.c_str(), "RECREATE");
        if (!fout || fout->IsZombie()) {
            AppendLog("ERROR: cannot create " + fname);
            delete fout;
            return;
        }
        c->Write();
        fout->Close(); delete fout;
    } else {
        c->Print(fname.c_str());
    }
    AppendLog("Canvas saved: " + fname);
    SetStatus("Saved: " + fname);
}

// ─────────────────────────────────────────────────────────────────────────────
// Canvas Annotations
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnApplyCanvasAnnotations()
{
    TCanvas* c = canvas_->GetCanvas();
    if (!c) { AppendLog("No canvas to annotate."); return; }

    c->cd();

    // Find the first TH1 drawn in the canvas to update title / axis labels
    TH1* h = nullptr;
    TIter nextObj(c->GetListOfPrimitives());
    TObject* obj;
    while ((obj = nextObj())) {
        if (obj->InheritsFrom(TH1::Class())) { h = (TH1*)obj; break; }
        // Also check inside TPad children
        if (obj->InheritsFrom(TPad::Class())) {
            TIter nextPad(((TPad*)obj)->GetListOfPrimitives());
            TObject* po;
            while ((po = nextPad())) {
                if (po->InheritsFrom(TH1::Class())) { h = (TH1*)po; break; }
            }
            if (h) break;
        }
    }

    if (h) {
        if (frTitleEntry_  && std::string(frTitleEntry_ ->GetText()).size())
            h->SetTitle(frTitleEntry_->GetText());
        if (frXLabelEntry_ && std::string(frXLabelEntry_->GetText()).size())
            h->GetXaxis()->SetTitle(frXLabelEntry_->GetText());
        if (frYLabelEntry_ && std::string(frYLabelEntry_->GetText()).size())
            h->GetYaxis()->SetTitle(frYLabelEntry_->GetText());
    }

    // Free-form TLatex annotation
    if (frAnnotText_) {
        std::string txt = frAnnotText_->GetText();
        if (!txt.empty()) {
            double ax   = frAnnotX_    ? frAnnotX_   ->GetNumber() : 0.15;
            double ay   = frAnnotY_    ? frAnnotY_   ->GetNumber() : 0.85;
            double asiz = frAnnotSize_ ? frAnnotSize_->GetNumber() : 0.040;
            TLatex* ltx = new TLatex();
            ltx->SetNDC();
            ltx->SetTextSize(asiz);
            ltx->DrawLatex(ax, ay, txt.c_str());
        }
    }

    c->Modified();
    c->Update();
    AppendLog("Canvas annotations applied.");
}

// ─────────────────────────────────────────────────────────────────────────────
// AutoFit — BG-apply helpers
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnApplyBgSelected()
{
    Int_t id = histList_->GetSelected();
    if (!inputFile_ || id < 1 || (size_t)id > histNames_.size()) {
        AppendLog("Select a histogram in the AutoFit list first.");
        return;
    }
    const std::string& hname = histNames_[id - 1];
    bool bgRawOwned = false;
    TH1* raw = LoadHistFromFile(hname, bgRawOwned);
    if (!raw) { AppendLog("Histogram not found: " + hname); return; }

    int iters = (int)bgIterEntry_->GetNumber();
    TH1* bgsub = MakeBgSubHist(raw, true, iters);
    if (bgRawOwned) { delete raw; raw = nullptr; }
    if (!bgsub) { AppendLog("Background subtraction failed for " + hname); return; }

    ApplyHistStyle(bgsub, hname.c_str());
    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();
    bgsub->SetBit(kCanDelete);
    bgsub->SetLineColor(kBlack);
    bgsub->SetMarkerSize(0);
    SetYMaxFromVisible(bgsub);
    bgsub->Draw("hist");
    bgsub->Draw("E1 same");
    c->Modified(); c->Update();

    AppendLog("BG subtracted (preview): " + hname + "  iters=" + std::to_string(iters));
    SetStatus("BG sub preview: " + hname);
}

void GammaFitGUI::OnApplyBgAll()
{
    if (!inputFile_) { AppendLog("No ROOT file loaded."); return; }
    if (histNames_.empty()) { AppendLog("No histograms loaded."); return; }

    int iters = (int)bgIterEntry_->GetNumber();
    const std::string outPath = "AllGammaFits.root";
    TFile* fout = TFile::Open(outPath.c_str(), "UPDATE");
    if (!fout || fout->IsZombie()) {
        delete fout;
        fout = TFile::Open(outPath.c_str(), "RECREATE");
    }

    int nDone = 0;
    for (const auto& hname : histNames_) {
        if (ClassOf(hname) != "Gamma") continue;
        bool batchOwned = false;
        TH1* raw = LoadHistFromFile(hname, batchOwned);
        if (!raw) continue;
        TH1* bgsub = MakeBgSubHist(raw, true, iters);
        if (batchOwned) { delete raw; raw = nullptr; }
        if (!bgsub) continue;
        fout->cd();
        bgsub->Write(Form("%s_bgsub", hname.c_str()), TObject::kOverwrite);
        delete bgsub;
        nDone++;
    }
    fout->Close(); delete fout;
    AppendLog("BG subtraction saved for " + std::to_string(nDone) + " histograms"
              "  iters=" + std::to_string(iters) + "  -> " + outPath);
    SetStatus("BG subtraction done: " + std::to_string(nDone) + " histograms");
}

void GammaFitGUI::OnResetBgSub()
{
    if (bgSubtractChk_) bgSubtractChk_->SetState(kButtonDown);
    if (bgIterEntry_)   bgIterEntry_->SetNumber(14);

    Int_t id = histList_->GetSelected();
    if (!inputFile_ || id < 1 || (size_t)id > histNames_.size()) {
        AppendLog("Background subtraction reset to defaults (14 iterations, enabled).");
        return;
    }
    const std::string& hname = histNames_[id - 1];
    bool owned = false;
    TH1* raw = LoadHistFromFile(hname, owned);
    if (!raw) { AppendLog("Histogram not found: " + hname); return; }

    ApplyHistStyle(raw, hname.c_str());
    if (ClassOf(hname) == "Decay") raw->GetXaxis()->SetTitle("Time (ms)");
    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();
    raw->SetLineColor(kBlack);
    raw->SetMarkerSize(0);
    SetYMaxFromVisible(raw);
    raw->Draw("hist");
    raw->Draw("E1 same");
    c->Modified(); c->Update();

    if (owned) delete raw;
    AppendLog("BG sub preview cleared — showing raw histogram: " + hname);
    SetStatus("Raw: " + hname);
}

// Source tab — helpers and slot implementations
// ─────────────────────────────────────────────────────────────────────────────

static int DateToJDN(const std::string& s)
{
    int y = 0, m = 0, d = 0;
    std::sscanf(s.c_str(), "%d-%d-%d", &y, &m, &d);
    int a  = (14 - m) / 12;
    int yr = y + 4800 - a;
    int mo = m + 12 * a - 3;
    return d + (153 * mo + 2) / 5 + 365 * yr + yr / 4 - yr / 100 + yr / 400 - 32045;
}

double GammaFitGUI::ComputeDecayedActivity() const
{
    if (srcActivity_ <= 0.0) return 0.0;
    if (srcHalflifeDays_ <= 0.0 || srcCalDate_.empty() || srcMeasDate_.empty())
        return srcActivity_;
    double dt = (double)(DateToJDN(srcMeasDate_) - DateToJDN(srcCalDate_));
    return srcActivity_ * std::exp(-std::log(2.0) * dt / srcHalflifeDays_);
}

void GammaFitGUI::ExtractPeaksFromCache(const std::string& hname)
{
    srcPeakEs_.clear();
    srcPeakCounts_.clear();
    srcPeakCountsErr_.clear();

    FitDatabase fitdb;
    if (!fitdb.Load(CacheFileFor(hname))) return;

    for (const auto& [key, entry] : fitdb.GetEntries()) {
        if (key.size() >= 2 && key[0] == '_' && key[1] == '_') continue;
        int npar = (int)entry.params.size();
        if (npar < 5 || (npar - 2) % 3 != 0) continue;
        int n = (npar - 2) / 3;
        for (int i = 0; i < n; i++) {
            double A   = entry.params[3 * i];
            double E   = entry.params[3 * i + 1];
            double sig = entry.params[3 * i + 2];
            if (A <= 0.0 || E <= 0.0 || sig <= 0.0) continue;
            double counts = A * sig * std::sqrt(2.0 * TMath::Pi());
            srcPeakEs_.push_back(E);
            srcPeakCounts_.push_back(counts);
            srcPeakCountsErr_.push_back(std::sqrt(std::max(counts, 1.0)));
        }
    }

    // Sort by energy ascending
    size_t N = srcPeakEs_.size();
    std::vector<size_t> idx(N);
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
        return srcPeakEs_[a] < srcPeakEs_[b];
    });
    std::vector<double> se, sc, serr;
    se.reserve(N); sc.reserve(N); serr.reserve(N);
    for (size_t i : idx) {
        se.push_back(srcPeakEs_[i]);
        sc.push_back(srcPeakCounts_[i]);
        serr.push_back(srcPeakCountsErr_[i]);
    }
    srcPeakEs_       = std::move(se);
    srcPeakCounts_   = std::move(sc);
    srcPeakCountsErr_= std::move(serr);

    AppendLog("Extracted " + std::to_string(srcPeakEs_.size()) +
              " fitted peaks from cache for " + hname);
}

void GammaFitGUI::PopulateSourceList()
{
    srcLineList_->RemoveAll();

    double A_meas   = ComputeDecayedActivity();
    double liveTime = srcLiveTime_->GetNumber();

    for (size_t i = 0; i < srcLines_.size(); i++) {
        const SourceLine& sl = srcLines_[i];
        std::ostringstream oss;
        oss << std::fixed << std::setprecision(3)
            << sl.energy << " keV | I=" << sl.intensity;

        if (sl.assigned >= 0 && sl.assigned < (int)srcPeakEs_.size()) {
            double fitE   = srcPeakEs_[sl.assigned];
            double counts = srcPeakCounts_[sl.assigned];
            oss << " | fit: " << fitE << " keV";
            if (A_meas > 0.0 && liveTime > 0.0 && sl.intensity > 0.0) {
                double eff = counts / (A_meas * sl.intensity * liveTime);
                oss << " | ε=" << std::scientific << std::setprecision(2) << eff;
            }
        } else {
            oss << " | (unassigned)";
        }
        srcLineList_->AddEntry(oss.str().c_str(), (Int_t)i + 1);
    }
    srcLineList_->MapSubwindows(); srcLineList_->Layout();
}

void GammaFitGUI::UpdateSourceInfoLabel()
{
    bool isuCi  = (srcActivityUnit_->GetSelected() == 2);
    double displayA = isuCi ? srcActivityRaw_ : srcActivity_;
    std::string unitStr = isuCi ? " µCi" : " Bq";

    std::ostringstream info;
    info << (srcIsotope_.empty() ? "?" : srcIsotope_);
    if (srcActivity_ > 0.0)
        info << std::fixed << std::setprecision(4) << "  A=" << displayA << unitStr;
    if (!srcCalDate_.empty())  info << "  cal=" << srcCalDate_;
    if (!srcMeasDate_.empty()) info << "  meas=" << srcMeasDate_;
    if (srcHalflifeDays_ > 0.0 && srcActivity_ > 0.0) {
        double Am = ComputeDecayedActivity();  // always Bq
        double displayAm = isuCi ? Am / 37000.0 : Am;
        info << std::scientific << std::setprecision(2)
             << "  A_meas=" << displayAm << unitStr;
    }
    srcInfoLbl_->SetText(info.str().c_str());
}

void GammaFitGUI::OnActivityUnitChanged(Int_t /*id*/)
{
    bool isuCi = (srcActivityUnit_->GetSelected() == 2);
    // Recompute the Bq value from the raw magnitude and the new unit
    srcActivity_ = srcActivityRaw_ * (isuCi ? 37000.0 : 1.0);
    UpdateSourceInfoLabel();
    // Refresh the list so efficiency numbers update
    PopulateSourceList();
}

void GammaFitGUI::OnLoadSourceFile()
{
    static const char* kSrcTypes[] = {
        "Text files", "*.txt",
        "All files",  "*",
        nullptr,      nullptr
    };
    TGFileInfo fi;
    fi.fFileTypes = kSrcTypes;
    fi.fIniDir    = StrDup(".");
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if (!fi.fFilename) return;

    std::ifstream in(fi.fFilename);
    if (!in.is_open()) {
        AppendLog("Cannot open source file: " + std::string(fi.fFilename));
        return;
    }

    // Reset state
    srcIsotope_.clear();
    srcActivityRaw_   = 0.0;
    srcActivity_      = 0.0;
    srcHalflifeDays_  = 0.0;
    srcCalDate_.clear();
    srcMeasDate_.clear();
    srcLines_.clear();

    std::string line;
    while (std::getline(in, line)) {
        // Strip inline comments
        auto cpos = line.find('#');
        if (cpos != std::string::npos) line = line.substr(0, cpos);
        size_t start = line.find_first_not_of(" \t\r\n");
        if (start == std::string::npos) continue;
        line = line.substr(start);

        std::istringstream ss(line);
        std::string tok;
        if (!(ss >> tok)) continue;

        if      (tok == "isotope")  { ss >> srcIsotope_; }
        else if (tok == "activity") {
            double val;
            if (ss >> val) {
                srcActivityRaw_ = val;
                // Check for optional unit token in the file
                std::string unit;
                if (ss >> unit) {
                    // Normalise: lowercase, drop special chars
                    std::string u; for (char ch : unit) u += std::tolower((unsigned char)ch);
                    if (u == "uci" || u == "muci" || u == "µci" ||
                        u == "microcurie" || u == "microcuries") {
                        srcActivityUnit_->Select(2, kFALSE);
                    } else {
                        srcActivityUnit_->Select(1, kFALSE);  // Bq or unknown → Bq
                    }
                }
                // Convert to Bq using the (possibly just-set) combo
                bool isuCi = (srcActivityUnit_->GetSelected() == 2);
                srcActivity_ = srcActivityRaw_ * (isuCi ? 37000.0 : 1.0);
            }
        }
        else if (tok == "halflife") { ss >> srcHalflifeDays_; }
        else if (tok == "caldate")  { ss >> srcCalDate_; }
        else if (tok == "measdate") { ss >> srcMeasDate_; }
        else if (tok == "livetime") {
            double lt;
            if (ss >> lt) srcLiveTime_->SetNumber(lt);
        } else {
            // Try two-number data line: energy_keV  branching_ratio
            try {
                double e = std::stod(tok);
                double br;
                if (ss >> br && e > 0.0 && br > 0.0 && br <= 1.0)
                    srcLines_.push_back({e, br, -1});
            } catch (...) {}
        }
    }
    in.close();

    // Update file label (basename only)
    std::string path = fi.fFilename;
    size_t slash = path.rfind('/');
    if (slash != std::string::npos) path = path.substr(slash + 1);
    srcFileLbl_->SetText(path.c_str());

    UpdateSourceInfoLabel();
    AppendLog("Source file loaded: " + std::string(fi.fFilename) +
              "  isotope=" + srcIsotope_ +
              "  lines=" + std::to_string(srcLines_.size()));
    PopulateSourceList();
}

void GammaFitGUI::OnOpenSourceRootFile()
{
    static const char* kRootTypes[] = {
        "ROOT files", "*.root", "All files", "*", nullptr, nullptr
    };
    TGFileInfo fi;
    fi.fFileTypes = kRootTypes;
    fi.fIniDir    = StrDup(".");
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if (!fi.fFilename) return;

    if (srcRootFile_) { srcRootFile_->Close(); delete srcRootFile_; srcRootFile_ = nullptr; }
    srcRootPath_ = fi.fFilename;
    srcRootFile_ = TFile::Open(srcRootPath_.c_str(), "READ");
    if (!srcRootFile_ || srcRootFile_->IsZombie()) {
        AppendLog("ERROR: cannot open " + srcRootPath_);
        srcRootFile_ = nullptr; return;
    }

    srcHistNames_.clear();
    srcTh2Names_.clear();
    srcProjParent_.clear();
    TIter next(srcRootFile_->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        if (!obj) continue;
        std::string name = obj->GetName();
        if (obj->InheritsFrom("TH2")) {
            srcHistNames_.push_back(name);
            srcTh2Names_.insert(name);
            srcHistNames_.push_back(name + "_px");
            srcHistNames_.push_back(name + "_py");
            srcProjParent_[name + "_px"] = name;
            srcProjParent_[name + "_py"] = name;
        } else if (obj->InheritsFrom("TH1")) {
            srcHistNames_.push_back(name);
        }
    }
    PopulateSrcHistCombo();

    std::string display = srcRootPath_;
    size_t slash = display.rfind('/');
    if (slash != std::string::npos) display = display.substr(slash + 1);
    srcRootFileLbl_->SetText(display.c_str());
    AppendLog("Source ROOT file: " + srcRootPath_ +
              "  (" + std::to_string(srcHistNames_.size()) + " histograms)");
}

void GammaFitGUI::PopulateSrcHistCombo()
{
    srcHistCombo_->RemoveAll();
    for (size_t i = 0; i < srcHistNames_.size(); i++) {
        const std::string& name = srcHistNames_[i];
        std::string display = name;
        if (srcTh2Names_.count(name))
            display = "[2D] " + name;
        else {
            auto it = srcProjParent_.find(name);
            if (it != srcProjParent_.end()) {
                bool isX = name.size() >= 3 && name.substr(name.size()-3) == "_px";
                display = (isX ? "[ProjX] " : "[ProjY] ") + it->second;
            }
        }
        // Only add non-TH2 entries to the source histogram combo for fitting
        if (!srcTh2Names_.count(name))
            srcHistCombo_->AddEntry(display.c_str(), (Int_t)i + 1);
    }
    srcHistCombo_->MapSubwindows(); srcHistCombo_->Layout();
}

std::string GammaFitGUI::SourceAnalysisFileFor(const std::string& hname) const
{
    return std::string(kCacheDir) + "/fit_cache_" + hname + "_source.root";
}

void GammaFitGUI::SaveSourceAnalysis()
{
    if (srcHist_.empty()) return;
    mkdir(kCacheDir, 0755);
    std::string fname = SourceAnalysisFileFor(srcHist_);
    TFile* fout = TFile::Open(fname.c_str(), "RECREATE");
    if (!fout || fout->IsZombie()) {
        AppendLog("ERROR: cannot write source analysis to " + fname);
        delete fout; return;
    }
    fout->cd();

    // Energy calibration graph
    {
        std::vector<double> refEs, residuals;
        for (const auto& sl : srcLines_) {
            if (sl.assigned < 0 || sl.assigned >= (int)srcPeakEs_.size()) continue;
            refEs.push_back(sl.energy);
            residuals.push_back(srcPeakEs_[sl.assigned] - sl.energy);
        }
        if (!refEs.empty()) {
            TGraph* gr = new TGraph((int)refEs.size(), refEs.data(), residuals.data());
            gr->SetName("energy_calib");
            gr->SetTitle("Energy Calibration;Reference Energy (keV);Fitted #minus Reference (keV)");
            gr->Write(); delete gr;
        }
    }

    // Efficiency graph
    {
        double A_meas = ComputeDecayedActivity();
        double liveTime = srcLiveTime_->GetNumber();
        bool absolute = (A_meas > 0.0 && liveTime > 0.0);
        std::vector<double> efEs, eff, effErr, exErr;
        for (const auto& sl : srcLines_) {
            if (sl.assigned < 0 || sl.assigned >= (int)srcPeakEs_.size()) continue;
            if (sl.intensity <= 0.0) continue;
            double counts = srcPeakCounts_[sl.assigned];
            double countsErr = srcPeakCountsErr_[sl.assigned];
            double denom = absolute ? (A_meas * sl.intensity * liveTime) : sl.intensity;
            if (denom <= 0.0) continue;
            efEs.push_back(srcPeakEs_[sl.assigned]);
            eff.push_back(counts / denom);
            effErr.push_back(countsErr / denom);
            exErr.push_back(0.0);
        }
        if (!efEs.empty()) {
            TGraphErrors* gr = new TGraphErrors((int)efEs.size(), efEs.data(), eff.data(),
                                                exErr.data(), effErr.data());
            gr->SetName("efficiency");
            std::string yTitle = absolute ? "Absolute efficiency" : "Relative efficiency";
            gr->SetTitle(Form("Efficiency;Energy (keV);%s", yTitle.c_str()));
            gr->Write(); delete gr;
        }
    }

    // FWHM data and model
    if (!fwhmAllX_.empty() && fwhmHistName_ == srcHist_) {
        std::vector<double> x, y, ex, ey;
        for (size_t i = 0; i < fwhmAllX_.size(); i++) {
            if (fwhmExcluded_[i]) continue;
            x.push_back(fwhmAllX_[i]); y.push_back(fwhmAllY_[i]);
            ex.push_back(0.0);         ey.push_back(0.0);
        }
        if (!x.empty()) {
            TGraphErrors* gr = new TGraphErrors((int)x.size(), x.data(), y.data(),
                                                ex.data(), ey.data());
            gr->SetName("fwhm_data");
            gr->SetTitle("FWHM vs Energy;Energy (keV);FWHM (keV)");
            gr->Write(); delete gr;
        }
        if (fwhmTF1_) fwhmTF1_->Write("fwhm_model");
    }

    fout->Close(); delete fout;
    AppendLog("Source analysis saved: " + fname);
}

void GammaFitGUI::OnRunSourceAutoFit()
{
    if (!srcRootFile_) {
        AppendLog("Open a source ROOT file first.");
        return;
    }
    Int_t id = srcHistCombo_->GetSelected();
    if (id < 1 || (size_t)id > srcHistNames_.size()) {
        AppendLog("Select a source histogram first.");
        return;
    }
    if (srcLines_.empty()) {
        AppendLog("Load a source description file first — its energies are used as fit seeds.");
        return;
    }
    srcHist_ = srcHistNames_[id - 1];
    for (auto& sl : srcLines_) sl.assigned = -1;

    std::vector<double> seeds;
    for (const auto& sl : srcLines_) seeds.push_back(sl.energy);

    RunFitOnHistogram(srcHist_, srcRootFile_, seeds);

    ExtractPeaksFromCache(srcHist_);
    PopulateSourceList();
}

void GammaFitGUI::OnLoadSourceCache()
{
    if (!srcRootFile_) {
        AppendLog("Open a source ROOT file first.");
        return;
    }
    Int_t id = srcHistCombo_->GetSelected();
    if (id < 1 || (size_t)id > srcHistNames_.size()) {
        AppendLog("Select a source histogram first.");
        return;
    }
    srcHist_     = srcHistNames_[id - 1];
    currentHist_ = srcHist_;

    if (rawHistOwned_) { delete rawHist_; rawHist_ = nullptr; rawHistOwned_ = false; }

    auto srcProjIt = srcProjParent_.find(srcHist_);
    if (srcProjIt != srcProjParent_.end()) {
        rawHist_ = LoadProjection(srcRootFile_, srcHist_, srcProjParent_,
                                  srcTh2XLabelEntry_, srcTh2YLabelEntry_);
        rawHistOwned_ = (rawHist_ != nullptr);
    } else {
        rawHist_ = (TH1*)srcRootFile_->Get(srcHist_.c_str());
        rawHistOwned_ = false;
    }
    if (!rawHist_) {
        AppendLog("Histogram not found in source file: " + srcHist_);
        return;
    }

    OnLoadCache();  // draws cached fits on canvas (uses currentHist_ / rawHist_)

    for (auto& sl : srcLines_) sl.assigned = -1;
    ExtractPeaksFromCache(srcHist_);
    PopulateSourceList();

    // Report whether a saved source analysis file exists
    std::string af = SourceAnalysisFileFor(srcHist_);
    struct stat st;
    if (stat(af.c_str(), &st) == 0)
        AppendLog("Source analysis cache exists: " + af);

    AppendLog("Source cache loaded: " + srcHist_ +
              "  " + std::to_string(srcPeakEs_.size()) + " peaks");
}

void GammaFitGUI::OnAutoIdentify()
{
    if (srcLines_.empty()) {
        AppendLog("Load a source description file first.");
        return;
    }
    if (srcPeakEs_.empty()) {
        AppendLog("Load cache or run AutoFit first to get fitted peaks.");
        return;
    }

    int nAssigned = 0;
    for (auto& sl : srcLines_) {
        sl.assigned = -1;
        double tol     = 3.0 * res_.FWHM(sl.energy);
        double bestDiff = tol;
        int    bestIdx  = -1;
        for (size_t j = 0; j < srcPeakEs_.size(); j++) {
            double diff = std::abs(srcPeakEs_[j] - sl.energy);
            if (diff < bestDiff) { bestDiff = diff; bestIdx = (int)j; }
        }
        if (bestIdx >= 0) { sl.assigned = bestIdx; nAssigned++; }
    }
    PopulateSourceList();
    AppendLog("Auto-identified " + std::to_string(nAssigned) +
              " / " + std::to_string(srcLines_.size()) + " source lines");
}

void GammaFitGUI::OnManualAssign()
{
    if (srcLines_.empty()) {
        AppendLog("Load a source description file first.");
        return;
    }
    if (srcPeakEs_.empty()) {
        AppendLog("Load cache or run AutoFit first.");
        return;
    }

    Int_t selId = srcLineList_->GetSelected();
    if (selId < 1 || (size_t)selId > srcLines_.size()) {
        AppendLog("Select a source line from the list first.");
        return;
    }

    double targetE = srcManualE_->GetNumber();
    if (targetE <= 0.0) {
        AppendLog("Enter a fitted energy > 0 keV.");
        return;
    }

    double tol     = 5.0 * res_.FWHM(targetE);
    int    bestIdx = -1;
    double best    = tol;
    for (size_t j = 0; j < srcPeakEs_.size(); j++) {
        double diff = std::abs(srcPeakEs_[j] - targetE);
        if (diff < best) { best = diff; bestIdx = (int)j; }
    }

    if (bestIdx < 0) {
        AppendLog("No fitted peak within 5 FWHM of " + Fmt(targetE) + " keV");
        return;
    }

    srcLines_[selId - 1].assigned = bestIdx;
    PopulateSourceList();
    AppendLog("Assigned source line " + Fmt(srcLines_[selId - 1].energy) +
              " keV → fitted peak " + Fmt(srcPeakEs_[bestIdx]) + " keV");
}

void GammaFitGUI::OnShowEnergyCalib()
{
    std::vector<double> refEs, residuals;
    for (const auto& sl : srcLines_) {
        if (sl.assigned < 0 || sl.assigned >= (int)srcPeakEs_.size()) continue;
        refEs.push_back(sl.energy);
        residuals.push_back(srcPeakEs_[sl.assigned] - sl.energy);
    }

    if (refEs.empty()) {
        AppendLog("No assigned peaks — run Auto Identify or Manual Assign first.");
        return;
    }

    TCanvas* c = canvas_->GetCanvas();
    c->Clear();
    c->cd();
    c->SetLogx(0); c->SetLogy(0);

    int n = (int)refEs.size();
    TGraph* gr = new TGraph(n, refEs.data(), residuals.data());
    gr->SetBit(kCanDelete);
    std::string iso = srcIsotope_.empty() ? "Source" : srcIsotope_;
    gr->SetTitle(Form("Energy Calibration — %s;Reference Energy (keV);Fitted #minus Reference (keV)",
                      iso.c_str()));
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.2);
    gr->SetMarkerColor(kBlue + 1);
    gr->Draw("AP");

    double xlo = *std::min_element(refEs.begin(), refEs.end()) * 0.95;
    double xhi = *std::max_element(refEs.begin(), refEs.end()) * 1.05;
    TLine* zeroLine = new TLine(xlo, 0.0, xhi, 0.0);
    zeroLine->SetBit(kCanDelete);
    zeroLine->SetLineColor(kRed);
    zeroLine->SetLineStyle(2);
    zeroLine->Draw("same");

    c->Modified(); c->Update();
    AppendLog("Energy calibration plot: " + std::to_string(n) + " assigned lines");
    SaveSourceAnalysis();
}

void GammaFitGUI::OnShowEfficiency()
{
    std::vector<double> efEs, eff, effErr;
    double A_meas   = ComputeDecayedActivity();
    double liveTime = srcLiveTime_->GetNumber();
    bool   absolute = (A_meas > 0.0 && liveTime > 0.0);

    for (const auto& sl : srcLines_) {
        if (sl.assigned < 0 || sl.assigned >= (int)srcPeakEs_.size()) continue;
        if (sl.intensity <= 0.0) continue;
        double counts    = srcPeakCounts_[sl.assigned];
        double countsErr = srcPeakCountsErr_[sl.assigned];
        double denom = absolute ? (A_meas * sl.intensity * liveTime) : sl.intensity;
        if (denom <= 0.0) continue;
        efEs.push_back(srcPeakEs_[sl.assigned]);
        eff.push_back(counts / denom);
        effErr.push_back(countsErr / denom);
    }

    if (efEs.empty()) {
        AppendLog("No assigned peaks — run Auto Identify or Manual Assign first.");
        return;
    }

    TCanvas* c = canvas_->GetCanvas();
    c->Clear();
    c->cd();
    c->SetLogx(1); c->SetLogy(1);

    int n = (int)efEs.size();
    std::vector<double> exErr(n, 0.0);
    TGraphErrors* gr = new TGraphErrors(n, efEs.data(), eff.data(),
                                         exErr.data(), effErr.data());
    gr->SetBit(kCanDelete);
    std::string iso    = srcIsotope_.empty() ? "Source" : srcIsotope_;
    std::string yTitle = absolute ? "Absolute efficiency  #varepsilon"
                                  : "Relative efficiency (counts/BR)";
    gr->SetTitle(Form("Detection Efficiency — %s;Energy (keV);%s",
                      iso.c_str(), yTitle.c_str()));
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.2);
    gr->SetMarkerColor(kBlue + 1);
    gr->SetLineColor(kBlue + 1);
    gr->Draw("AP");

    if (!absolute)
        AppendLog("Note: no source activity set — showing relative efficiency");

    c->Modified(); c->Update();

    std::string msg = "Efficiency plot: " + std::to_string(n) + " points";
    if (absolute)
        msg += Form("  A_meas=%.2e Bq  t=%.1f s", A_meas, liveTime);
    AppendLog(msg);
    SaveSourceAnalysis();
}

void GammaFitGUI::OnShowSourceFWHM()
{
    Int_t id = srcHistCombo_->GetSelected();
    if (id < 1 || (size_t)id > srcHistNames_.size()) {
        AppendLog("Select a source histogram first.");
        return;
    }
    // Sync srcHist_ so SaveSourceAnalysis knows which file to write
    if ((size_t)id <= srcHistNames_.size())
        srcHist_ = srcHistNames_[id - 1];
    // Sync the FWHM tab combo to the same histogram and call existing logic
    fwhmCombo_->Select(id, kFALSE);
    OnLoadFWHM();
    SaveSourceAnalysis();
}

void GammaFitGUI::AppendLog(const std::string& msg)
{
    logView_->AddLine(msg.c_str());
    logView_->ShowBottom();
    std::cout << "[GUI] " << msg << "\n";
}

void GammaFitGUI::SetStatus(const std::string& msg)
{
    statusBar_->SetText(msg.c_str(), 0);
    if (!inputPath_.empty())
        statusBar_->SetText(inputPath_.c_str(), 1);
}

// ─────────────────────────────────────────────────────────────────────────────
// Histogram classification helpers
// ─────────────────────────────────────────────────────────────────────────────

std::string GammaFitGUI::ClassOf(const std::string& name) const
{
    if (th2Names_.count(name)) return "2D";
    auto it = histClass_.find(name);
    return (it != histClass_.end()) ? it->second : "Gamma";
}

std::string GammaFitGUI::MetadataFileFor() const
{
    if (inputPath_.empty()) return "";
    // Sanitise path: replace non-alphanumeric with '_'
    std::string tag = inputPath_;
    for (char& c : tag)
        if (!std::isalnum((unsigned char)c) && c != '.') c = '_';
    return std::string(kCacheDir) + "/metadata_" + tag + ".txt";
}

void GammaFitGUI::SaveMetadata() const
{
    std::string path = MetadataFileFor();
    if (path.empty()) return;
    mkdir(kCacheDir, 0755);
    std::ofstream out(path);
    if (!out.is_open()) return;
    out << "# classifications\n";
    for (const auto& kv : histClass_)
        out << "class\t" << kv.first << "\t" << kv.second << "\n";
    out << "# custom projections\n";
    for (const auto& kv : customProjDefs_) {
        const auto& d = kv.second;
        out << "cproj\t" << kv.first << "\t" << d.th2Name
            << "\t" << (d.projX ? 1 : 0)
            << "\t" << std::fixed << std::setprecision(6) << d.lo
            << "\t" << d.hi << "\n";
    }
}

void GammaFitGUI::LoadMetadata()
{
    std::string path = MetadataFileFor();
    if (path.empty()) return;
    std::ifstream in(path);
    if (!in.is_open()) return;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string tag;
        ss >> tag;
        if (tag == "class") {
            std::string name, cls;
            ss >> name >> cls;
            if (!name.empty() && !cls.empty()) histClass_[name] = cls;
        } else if (tag == "cproj") {
            std::string name, th2name;
            int pxFlag = 1;
            double lo = 0, hi = 0;
            ss >> name >> th2name >> pxFlag >> lo >> hi;
            if (!name.empty() && !th2name.empty()) {
                CustomProjDef def;
                def.th2Name = th2name;
                def.projX   = (pxFlag != 0);
                def.lo      = lo;
                def.hi      = hi;
                customProjDefs_[name] = def;
                // Re-add to histNames_ if not already present
                if (std::find(histNames_.begin(), histNames_.end(), name) == histNames_.end())
                    histNames_.push_back(name);
            }
        }
    }
}

void GammaFitGUI::OnHistClassSet()
{
    Int_t id = histList_->GetSelected();
    if (id < 1 || (size_t)id > histNames_.size()) {
        AppendLog("Select a histogram in the list first."); return;
    }
    const std::string& name = histNames_[id - 1];
    if (th2Names_.count(name)) {
        AppendLog("[Class] 2D histograms are classified automatically."); return;
    }
    int sel = histClassCombo_ ? histClassCombo_->GetSelected() : 1;
    std::string cls = (sel == 2) ? "Decay" : (sel == 3) ? "2D" : "Gamma";
    histClass_[name] = cls;
    SaveMetadata();
    PopulateHistWidgets();
    AppendLog("[Class] " + name + " → " + cls);
    SetStatus("Classified: " + name + " = " + cls);
}

// ─────────────────────────────────────────────────────────────────────────────
// Custom projection helpers and slot
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::PopulateCustProjTh2Combo()
{
    if (!custProjTh2Combo_) return;
    custProjTh2Combo_->RemoveAll();
    int idx = 1;
    for (const auto& name : th2Names_)
        custProjTh2Combo_->AddEntry(name.c_str(), idx++);
    custProjTh2Combo_->MapSubwindows();
    custProjTh2Combo_->Layout();
}

void GammaFitGUI::OnAddCustomProjection()
{
    if (!inputFile_) { AppendLog("[CustomProj] No ROOT file loaded."); return; }

    // Gather inputs
    TGLBEntry* th2Sel = custProjTh2Combo_ ? custProjTh2Combo_->GetSelectedEntry() : nullptr;
    if (!th2Sel) { AppendLog("[CustomProj] Select a TH2 first."); return; }
    std::string th2name = th2Sel->GetTitle();

    int axisId = custProjAxisCombo_ ? custProjAxisCombo_->GetSelected() : 1;
    bool projX = (axisId == 1);  // Project onto X → cut on Y
    double lo  = custProjLo_ ? custProjLo_->GetNumber() : 0.0;
    double hi  = custProjHi_ ? custProjHi_->GetNumber() : 0.0;
    if (lo >= hi) { AppendLog("[CustomProj] Range lo must be < hi."); return; }

    std::string name = custProjName_ ? std::string(custProjName_->GetText()) : "";
    if (name.empty()) { AppendLog("[CustomProj] Enter a name."); return; }

    // Reject duplicates
    if (std::find(histNames_.begin(), histNames_.end(), name) != histNames_.end()) {
        AppendLog("[CustomProj] Name '" + name + "' already in use."); return;
    }

    // Validate by creating the projection
    TH2* h2 = (TH2*)inputFile_->Get(th2name.c_str());
    if (!h2) { AppendLog("[CustomProj] TH2 '" + th2name + "' not found."); return; }

    TH1* test = nullptr;
    if (projX) {
        int b1 = h2->GetYaxis()->FindBin(lo);
        int b2 = h2->GetYaxis()->FindBin(hi);
        test = h2->ProjectionX(name.c_str(), b1, b2);
    } else {
        int b1 = h2->GetXaxis()->FindBin(lo);
        int b2 = h2->GetXaxis()->FindBin(hi);
        test = h2->ProjectionY(name.c_str(), b1, b2);
    }
    if (!test || test->GetEntries() == 0) {
        delete test;
        AppendLog("[CustomProj] Projection is empty — check range and TH2."); return;
    }
    delete test;

    // Store definition
    CustomProjDef def;
    def.th2Name = th2name;
    def.projX   = projX;
    def.lo      = lo;
    def.hi      = hi;
    customProjDefs_[name] = def;
    histNames_.push_back(name);

    SaveMetadata();
    PopulateHistWidgets();
    PopulateDecayTh2Combo();

    AppendLog(Form("[CustomProj] Added '%s' from %s  [%.3g, %.3g]  Proj%s",
                   name.c_str(), th2name.c_str(), lo, hi, projX ? "X" : "Y"));
    SetStatus("Added custom projection: " + name);
}
