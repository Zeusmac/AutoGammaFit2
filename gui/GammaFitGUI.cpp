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
#include "TPaveText.h"
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
#include <dirent.h>

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
        clearCacheBtn->SetToolTipText("Delete the live cache for the current histogram (archive is preserved)");

        TGTextButton* restoreArchBtn = new TGTextButton(btnRow, "  Restore Archive  ");
        btnRow->AddFrame(restoreArchBtn, new TGLayoutHints(kLHintsRight, 2, 2, 1, 1));
        restoreArchBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRestoreArchivedCache()");
        restoreArchBtn->SetToolTipText("Restore a previously archived cache for the current histogram");

        TGTextButton* archiveCacheBtn = new TGTextButton(btnRow, "  Archive Cache  ");
        btnRow->AddFrame(archiveCacheBtn, new TGLayoutHints(kLHintsRight, 2, 2, 1, 1));
        archiveCacheBtn->Connect("Clicked()", "GammaFitGUI", this, "OnArchiveHistCache()");
        archiveCacheBtn->SetToolTipText("Copy the current cache to fit_caches/archive/ with a timestamp");

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
void GammaFitGUI::BuildAutoFitTab(TGCompositeFrame* tab)
{
    // ── Scrollable wrapper ─────────────────────────────────────────────────────
    TGCanvas* scroller = new TGCanvas(tab, 10, 600, kSunkenFrame);
    tab->AddFrame(scroller, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 0, 0, 0, 0));
    TGCompositeFrame* p = new TGCompositeFrame(scroller->GetViewPort(), 1, 1, kVerticalFrame);
    scroller->SetContainer(p);

    // ── File ──────────────────────────────────────────────────────────────────
    TGGroupFrame* fg = new TGGroupFrame(p, "File");
    p->AddFrame(fg, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

    TGTextButton* openBtn = new TGTextButton(fg, "Open ROOT File...");
    fg->AddFrame(openBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    openBtn->Connect("Clicked()", "GammaFitGUI", this, "OnOpenFile()");
    openBtn->SetToolTipText("Browse for a ROOT file containing TH1 gamma spectra");

    // Recent files row
    {
        TGHorizontalFrame* rr = new TGHorizontalFrame(fg);
        fg->AddFrame(rr, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
        recentCombo_ = new TGComboBox(rr, 950);
        recentCombo_->AddEntry("(no recent files)", 1);
        recentCombo_->Select(1, kFALSE);
        recentCombo_->Resize(185, 22);
        rr->AddFrame(recentCombo_, new TGLayoutHints(kLHintsExpandX, 0, 4, 0, 0));
        TGTextButton* recentBtn = new TGTextButton(rr, "Open");
        rr->AddFrame(recentBtn, new TGLayoutHints(kLHintsRight));
        recentBtn->Connect("Clicked()", "GammaFitGUI", this, "OnOpenRecent()");
        recentBtn->SetToolTipText("Open the selected recent ROOT file");
        LoadRecentFiles();
        RefreshRecentCombo();
    }

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
        cr->AddFrame(new TGLabel(cr, "Histogram Type:"),
                     new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        histClassCombo_ = new TGComboBox(cr, 900);
        histClassCombo_->AddEntry("Gamma Spectrum", 1);
        histClassCombo_->AddEntry("Decay Curve",    2);
        histClassCombo_->AddEntry("2D Histogram",   3);
        histClassCombo_->AddEntry("Background",     4);
        histClassCombo_->Select(1, kFALSE);
        histClassCombo_->Resize(140, 22);
        cr->AddFrame(histClassCombo_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        TGTextButton* setClassBtn = new TGTextButton(cr, "Set");
        cr->AddFrame(setClassBtn, new TGLayoutHints(kLHintsLeft));
        setClassBtn->Connect("Clicked()", "GammaFitGUI", this, "OnHistClassSet()");
        setClassBtn->SetToolTipText("Assign the chosen classification to the selected histogram");
    }

    // Delete histogram row
    {
        TGHorizontalFrame* dr = new TGHorizontalFrame(hg);
        hg->AddFrame(dr, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
        TGTextButton* delBtn = new TGTextButton(dr, "Delete Selected from File");
        delBtn->SetForegroundColor(0xCC0000);
        dr->AddFrame(delBtn, new TGLayoutHints(kLHintsExpandX));
        delBtn->Connect("Clicked()", "GammaFitGUI", this, "OnDeleteHistogram()");
        delBtn->SetToolTipText(
            "Permanently removes the selected histogram from the ROOT file.\n"
            "Virtual histograms (projections, bg-subtracted) are removed from\n"
            "the session only — they are not stored in the file.");
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
    rg->AddFrame(runAll, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
    runAll->Connect("Clicked()", "GammaFitGUI", this, "OnRunAll()");
    runAll->SetToolTipText("Fit all histograms classified as 'Gamma Spectrum'");

    TGTextButton* clearCacheAutoFit = new TGTextButton(rg, "Clear Cache  (selected)");
    rg->AddFrame(clearCacheAutoFit, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    clearCacheAutoFit->Connect("Clicked()", "GammaFitGUI", this, "OnClearHistCache()");
    clearCacheAutoFit->SetToolTipText("Delete the live cache for the selected histogram (archive is preserved)");

    TGTextButton* archiveCacheAutoFit = new TGTextButton(rg, "Archive Cache  (selected)");
    rg->AddFrame(archiveCacheAutoFit, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
    archiveCacheAutoFit->Connect("Clicked()", "GammaFitGUI", this, "OnArchiveHistCache()");
    archiveCacheAutoFit->SetToolTipText("Copy the selected histogram's cache to fit_caches/archive/ with a timestamp");

    TGTextButton* restoreArchAutoFit = new TGTextButton(rg, "Restore Archived Cache...");
    rg->AddFrame(restoreArchAutoFit, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
    restoreArchAutoFit->Connect("Clicked()", "GammaFitGUI", this, "OnRestoreArchivedCache()");
    restoreArchAutoFit->SetToolTipText("Restore a previously archived cache for the selected histogram");

    TGTextButton* transferBtn = new TGTextButton(rg, "Transfer Cache From...");
    rg->AddFrame(transferBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));
    transferBtn->Connect("Clicked()", "GammaFitGUI", this, "OnTransferCache()");
    transferBtn->SetToolTipText("Copy fit entries from another histogram's cache into the selected histogram's cache");

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

    // ── Background Histogram Subtraction ─────────────────────────────────────
    TGGroupFrame* bsg = new TGGroupFrame(p, "Background Subtraction (Histogram)");
    p->AddFrame(bsg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

    // Row: source histogram
    {
        TGHorizontalFrame* r1 = new TGHorizontalFrame(bsg);
        bsg->AddFrame(r1, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 1));
        r1->AddFrame(new TGLabel(r1, "Source:"),
                     new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        bgSubSrcCombo_ = new TGComboBox(r1, 930);
        bgSubSrcCombo_->Resize(200, 22);
        r1->AddFrame(bgSubSrcCombo_, new TGLayoutHints(kLHintsExpandX));
    }
    // Row: background histogram (only Background-typed)
    {
        TGHorizontalFrame* r2 = new TGHorizontalFrame(bsg);
        bsg->AddFrame(r2, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
        r2->AddFrame(new TGLabel(r2, "Background:"),
                     new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        bgSubBgCombo_ = new TGComboBox(r2, 931);
        bgSubBgCombo_->Resize(200, 22);
        r2->AddFrame(bgSubBgCombo_, new TGLayoutHints(kLHintsExpandX));
    }
    // Row: scale factor + result name
    {
        TGHorizontalFrame* r3 = new TGHorizontalFrame(bsg);
        bsg->AddFrame(r3, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
        r3->AddFrame(new TGLabel(r3, "Scale:"),
                     new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        bgSubScaleEntry_ = new TGNumberEntry(r3, 1.0, 6, -1,
            TGNumberFormat::kNESRealThree,
            TGNumberFormat::kNEAPositive,
            TGNumberFormat::kNELLimitMin, 0.0);
        bgSubScaleEntry_->SetWidth(70);
        r3->AddFrame(bgSubScaleEntry_, new TGLayoutHints(kLHintsLeft, 0, 8, 0, 0));
        r3->AddFrame(new TGLabel(r3, "Name:"),
                     new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        bgSubNameEntry_ = new TGTextEntry(r3, "subtracted");
        bgSubNameEntry_->SetWidth(100);
        r3->AddFrame(bgSubNameEntry_, new TGLayoutHints(kLHintsExpandX));
    }
    {
        TGTextButton* subBtn = new TGTextButton(bsg, "Create Subtracted Histogram");
        bsg->AddFrame(subBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));
        subBtn->Connect("Clicked()", "GammaFitGUI", this, "OnSubtractHistogram()");
        subBtn->SetToolTipText(
            "Creates: result = source - scale * background\n"
            "Sumw2() is applied so error bars are correctly propagated.\n"
            "Tag source histograms as 'Background' type to populate the background list.");
    }

    // ── Rebin Histogram ───────────────────────────────────────────────────────
    TGGroupFrame* rbg = new TGGroupFrame(p, "Rebin Histogram");
    p->AddFrame(rbg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

    {
        TGHorizontalFrame* r1 = new TGHorizontalFrame(rbg);
        rbg->AddFrame(r1, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        r1->AddFrame(new TGLabel(r1, "Rebin Factor:"),
                     new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
        rebinEntry_ = new TGNumberEntry(r1, 1, 4, -1,
            TGNumberFormat::kNESInteger,
            TGNumberFormat::kNEAPositive,
            TGNumberFormat::kNELLimitMinMax, 1, 1000);
        rebinEntry_->SetWidth(60);
        r1->AddFrame(rebinEntry_, new TGLayoutHints(kLHintsLeft));
        rebinEntry_->GetNumberEntry()->SetToolTipText(
            "Merge N bins into one. 1 = no rebin. Stored per histogram in metadata.");
    }
    {
        TGHorizontalFrame* r2 = new TGHorizontalFrame(rbg);
        rbg->AddFrame(r2, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
        TGTextButton* prevBtn = new TGTextButton(r2, " Preview ");
        r2->AddFrame(prevBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        prevBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRebinPreview()");
        prevBtn->SetToolTipText("Show rebinned histogram on canvas without storing the factor");
        TGTextButton* applyBtn = new TGTextButton(r2, " Apply ");
        r2->AddFrame(applyBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        applyBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRebinApply()");
        applyBtn->SetToolTipText("Store this rebin factor for the selected histogram (saved to cache metadata)");
        TGTextButton* resetBtn = new TGTextButton(r2, " Reset ");
        r2->AddFrame(resetBtn, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
        resetBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRebinReset()");
        resetBtn->SetToolTipText("Remove stored rebin factor and redisplay at original binning");
    }

    // ── S/N Pre-screen ────────────────────────────────────────────────────────
    TGGroupFrame* sng = new TGGroupFrame(p, "S/N Pre-screen (AutoFit)");
    p->AddFrame(sng, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

    {
        TGHorizontalFrame* r = new TGHorizontalFrame(sng);
        sng->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));
        r->AddFrame(new TGLabel(r, "Min S/N ratio:"),
                    new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
        snThreshEntry_ = new TGNumberEntry(r, 0.0, 6, -1,
            TGNumberFormat::kNESRealOne,
            TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMin, 0.0);
        r->AddFrame(snThreshEntry_, new TGLayoutHints(kLHintsLeft));
        snThreshEntry_->GetNumberEntry()->SetToolTipText(
            "Minimum S/N for a peak group to be fitted (0 = all groups fitted).\n"
            "S/N = signal in peak window / sqrt(scaled sideband noise).\n"
            "Set to 3-5 to suppress weak/spurious peaks.");
    }

    // ── Efficiency Correction ─────────────────────────────────────────────────
    TGGroupFrame* effg = new TGGroupFrame(p, "Efficiency Correction");
    p->AddFrame(effg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

    {
        TGLabel* effLbl = new TGLabel(effg, "ln(ε) = a − b·ln(E) + c·ln(E)² − d/E²");
        effg->AddFrame(effLbl, new TGLayoutHints(kLHintsLeft, 4, 4, 2, 2));
    }
    {
        TGHorizontalFrame* ra = new TGHorizontalFrame(effg);
        effg->AddFrame(ra, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 1));
        ra->AddFrame(new TGLabel(ra, "a:"), new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        effA_ = new TGNumberEntry(ra, 0.0, 8, -1,
            TGNumberFormat::kNESRealThree, TGNumberFormat::kNEAAnyNumber);
        effA_->SetWidth(80); ra->AddFrame(effA_, new TGLayoutHints(kLHintsLeft));
        ra->AddFrame(new TGLabel(ra, "  b:"), new TGLayoutHints(kLHintsCenterY, 4, 4, 0, 0));
        effB_ = new TGNumberEntry(ra, 0.0, 8, -1,
            TGNumberFormat::kNESRealThree, TGNumberFormat::kNEAAnyNumber);
        effB_->SetWidth(80); ra->AddFrame(effB_, new TGLayoutHints(kLHintsLeft));
    }
    {
        TGHorizontalFrame* rb = new TGHorizontalFrame(effg);
        effg->AddFrame(rb, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
        rb->AddFrame(new TGLabel(rb, "c:"), new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        effC_ = new TGNumberEntry(rb, 0.0, 8, -1,
            TGNumberFormat::kNESRealThree, TGNumberFormat::kNEAAnyNumber);
        effC_->SetWidth(80); rb->AddFrame(effC_, new TGLayoutHints(kLHintsLeft));
        rb->AddFrame(new TGLabel(rb, "  d:"), new TGLayoutHints(kLHintsCenterY, 4, 4, 0, 0));
        effD_ = new TGNumberEntry(rb, 0.0, 8, -1,
            TGNumberFormat::kNESRealThree, TGNumberFormat::kNEAAnyNumber);
        effD_->SetWidth(80); rb->AddFrame(effD_, new TGLayoutHints(kLHintsLeft));
    }
    {
        TGTextButton* applyEffBtn = new TGTextButton(effg, "Apply (update peak stats)");
        effg->AddFrame(applyEffBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));
        applyEffBtn->Connect("Clicked()", "GammaFitGUI", this, "OnApplyEfficiency()");
        applyEffBtn->SetToolTipText(
            "Store efficiency parameters. Corrected area = raw area / eps(E) will\n"
            "appear in peak statistics after the next manual fit.");
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

    TGTextButton* calibBtn = new TGTextButton(pg, "Energy Calibration (fitted - reference)");
    pg->AddFrame(calibBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    calibBtn->Connect("Clicked()", "GammaFitGUI", this, "OnShowEnergyCalib()");
    calibBtn->SetToolTipText("Plot (fitted E - reference E) vs reference E for all assigned lines");

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

    // ── Channel → keV Calibration ─────────────────────────────────────────────
    TGGroupFrame* calg = new TGGroupFrame(p, "Channel → keV Calibration");
    p->AddFrame(calg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

    {
        TGLabel* calLbl = new TGLabel(calg, "E(ch) = a + b·ch + c·ch²");
        calg->AddFrame(calLbl, new TGLayoutHints(kLHintsLeft, 4, 4, 2, 2));
    }
    {
        TGHorizontalFrame* ra = new TGHorizontalFrame(calg);
        calg->AddFrame(ra, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
        ra->AddFrame(new TGLabel(ra, "a (offset):"),
                     new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        calibA_ = new TGNumberEntry(ra, 0.0, 9, -1,
            TGNumberFormat::kNESRealThree, TGNumberFormat::kNEAAnyNumber);
        calibA_->SetWidth(90); ra->AddFrame(calibA_, new TGLayoutHints(kLHintsLeft));
    }
    {
        TGHorizontalFrame* rb = new TGHorizontalFrame(calg);
        calg->AddFrame(rb, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
        rb->AddFrame(new TGLabel(rb, "b (keV/ch): "),
                     new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        calibB_ = new TGNumberEntry(rb, 1.0, 9, -1,
            TGNumberFormat::kNESRealThree, TGNumberFormat::kNEAAnyNumber);
        calibB_->SetWidth(90); rb->AddFrame(calibB_, new TGLayoutHints(kLHintsLeft));
    }
    {
        TGHorizontalFrame* rc = new TGHorizontalFrame(calg);
        calg->AddFrame(rc, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
        rc->AddFrame(new TGLabel(rc, "c (keV/ch²):"),
                     new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        calibC_ = new TGNumberEntry(rc, 0.0, 9, -1,
            TGNumberFormat::kNESRealThree, TGNumberFormat::kNEAAnyNumber);
        calibC_->SetWidth(90); rc->AddFrame(calibC_, new TGLayoutHints(kLHintsLeft));
    }
    {
        TGTextButton* applyCalibBtn = new TGTextButton(calg, "Apply to Selected Histogram");
        calg->AddFrame(applyCalibBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));
        applyCalibBtn->Connect("Clicked()", "GammaFitGUI", this, "OnApplyCalibration()");
        applyCalibBtn->SetToolTipText(
            "Rescale the x-axis of the selected histogram using E = a + b·ch + c·ch².\n"
            "Creates a clone named '<hist>_cal' with the calibrated axis.\n"
            "Coefficients can be read from the Energy Calibration residual plot.");
    }

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

    TGLabel* inst = new TGLabel(p, "  <-- Click on spectrum to place peak markers");
    p->AddFrame(inst, new TGLayoutHints(kLHintsLeft, 4, 4, 0, 2));

    // ── Peak parameters (popup) ───────────────────────────────────────────────
    TGTextButton* fitParamBtn = new TGTextButton(p, "Fit Parameters...");
    p->AddFrame(fitParamBtn, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));
    fitParamBtn->Connect("Clicked()", "GammaFitGUI", this, "OnShowFitParams()");
    fitParamBtn->SetToolTipText("Edit fit seed values and parameter bounds");

    // Pre-create the dialog now so mEnergy_, mSigma_, etc. are valid immediately
    // (OnShowFitParams builds the dialog without mapping when fitParamDlg_ is null)
    OnShowFitParams();

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
    mBgFlatChk_ = new TGCheckButton(bgGrp, "Flat background  (constant, BG slope = 0)");
    bgGrp->AddFrame(mBgFlatChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
    mBgFlatChk_->SetToolTipText("Fix the linear BG slope to zero — use a flat (constant) background");

    // BG anchor regions — two off-peak regions used to seed bg0 and bg1
    bgGrp->AddFrame(new TGLabel(bgGrp, "Anchor regions (seed BG from two off-peak areas):"),
                    new TGLayoutHints(kLHintsLeft, 2, 2, 4, 1));
    {
        TGHorizontalFrame* ar1 = new TGHorizontalFrame(bgGrp);
        bgGrp->AddFrame(ar1, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
        ar1->AddFrame(new TGLabel(ar1, "Left  Lo:"),
                      new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        mBgAnch1Lo_ = new TGNumberEntry(ar1, 0.0, 7, -1,
            TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMinMax, 0.0, 10000.0);
        ar1->AddFrame(mBgAnch1Lo_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        ar1->AddFrame(new TGLabel(ar1, "Hi:"),
                      new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        mBgAnch1Hi_ = new TGNumberEntry(ar1, 0.0, 7, -1,
            TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMinMax, 0.0, 10000.0);
        ar1->AddFrame(mBgAnch1Hi_, new TGLayoutHints(kLHintsLeft));
    }
    {
        TGHorizontalFrame* ar2 = new TGHorizontalFrame(bgGrp);
        bgGrp->AddFrame(ar2, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
        ar2->AddFrame(new TGLabel(ar2, "Right Lo:"),
                      new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        mBgAnch2Lo_ = new TGNumberEntry(ar2, 0.0, 7, -1,
            TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMinMax, 0.0, 10000.0);
        ar2->AddFrame(mBgAnch2Lo_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        ar2->AddFrame(new TGLabel(ar2, "Hi:"),
                      new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        mBgAnch2Hi_ = new TGNumberEntry(ar2, 0.0, 7, -1,
            TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMinMax, 0.0, 10000.0);
        ar2->AddFrame(mBgAnch2Hi_, new TGLayoutHints(kLHintsLeft));
    }
    {
        TGTextButton* anchBtn = new TGTextButton(bgGrp, "Apply Anchors → seed bg0, bg1");
        bgGrp->AddFrame(anchBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));
        anchBtn->Connect("Clicked()", "GammaFitGUI", this, "OnApplyBgAnchors()");
        anchBtn->SetToolTipText(
            "Compute mean counts in each anchor region and fit a line through the\n"
            "two region centres. Seeds bg0 and bg1 for the next manual fit.\n"
            "Anchor regions should be off-peak (no gamma lines).");
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
        clearRangeBtn->SetToolTipText("Reset to automatic range: peaks +/- Range*sigma");
    }

    // ── Peak navigation ───────────────────────────────────────────────────────
    TGGroupFrame* navGrp = new TGGroupFrame(p, "Peak Navigation");
    p->AddFrame(navGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    {
        TGHorizontalFrame* navBtnRow = new TGHorizontalFrame(navGrp);
        navGrp->AddFrame(navBtnRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        TGTextButton* prevBtn = new TGTextButton(navBtnRow, "<< Prev");
        navBtnRow->AddFrame(prevBtn, new TGLayoutHints(kLHintsLeft, 2, 4, 0, 0));
        prevBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPrevPeak()");
        prevBtn->SetToolTipText("Zoom to the previous cached peak");

        TGTextButton* nextBtn = new TGTextButton(navBtnRow, "Next >>");
        navBtnRow->AddFrame(nextBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
        nextBtn->Connect("Clicked()", "GammaFitGUI", this, "OnNextPeak()");
        nextBtn->SetToolTipText("Zoom to the next cached peak");

        TGTextButton* zoomInBtn = new TGTextButton(navBtnRow, " + ");
        navBtnRow->AddFrame(zoomInBtn, new TGLayoutHints(kLHintsLeft, 8, 2, 0, 0));
        zoomInBtn->Connect("Clicked()", "GammaFitGUI", this, "OnZoomIn()");
        zoomInBtn->SetToolTipText("Narrow the x-axis view around the current peak");

        TGTextButton* zoomOutBtn = new TGTextButton(navBtnRow, " - ");
        navBtnRow->AddFrame(zoomOutBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
        zoomOutBtn->Connect("Clicked()", "GammaFitGUI", this, "OnZoomOut()");
        zoomOutBtn->SetToolTipText("Widen the x-axis view around the current peak");
    }
    peakNavLbl_ = new TGLabel(navGrp, "No peaks loaded");
    navGrp->AddFrame(peakNavLbl_, new TGLayoutHints(kLHintsLeft, 4, 4, 0, 2));

    // X range row — manual zoom to a specific energy range
    {
        TGHorizontalFrame* xrRow = new TGHorizontalFrame(navGrp);
        navGrp->AddFrame(xrRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
        xrRow->AddFrame(new TGLabel(xrRow, "X range:"),
                        new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        navXMinEntry_ = new TGNumberEntry(xrRow, 0.0, 7, -1,
            TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber);
        xrRow->AddFrame(navXMinEntry_, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
        xrRow->AddFrame(new TGLabel(xrRow, "–"),
                        new TGLayoutHints(kLHintsCenterY, 2, 2, 0, 0));
        navXMaxEntry_ = new TGNumberEntry(xrRow, 0.0, 7, -1,
            TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber);
        xrRow->AddFrame(navXMaxEntry_, new TGLayoutHints(kLHintsLeft, 2, 4, 0, 0));
        TGTextButton* goBtn = new TGTextButton(xrRow, " Go ");
        goBtn->Connect("Clicked()", "GammaFitGUI", this, "OnNavXRangeGo()");
        goBtn->SetToolTipText("Zoom to the specified x range");
        xrRow->AddFrame(goBtn, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
    }

    // Fit selector dropdown — selecting a fit zooms and draws exactly like Prev/Next
    {
        TGHorizontalFrame* resRow = new TGHorizontalFrame(navGrp);
        navGrp->AddFrame(resRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        resRow->AddFrame(new TGLabel(resRow, "Select fit:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        residualCombo_ = new TGComboBox(resRow, 600);
        residualCombo_->Resize(195, 22);
        resRow->AddFrame(residualCombo_, new TGLayoutHints(kLHintsExpandX));
        residualCombo_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                "OnSelectResidualFit(Int_t)");
    }

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
    addPeakChk_->Connect("Toggled(Bool_t)", "GammaFitGUI", this, "OnClearPeaks()");

    TGTextButton* clrBtn = new TGTextButton(pkRow, "Clear Clicks");
    pkRow->AddFrame(clrBtn, new TGLayoutHints(kLHintsRight, 0, 2, 0, 0));
    clrBtn->SetToolTipText("Clear peak list and remove all click markers from the canvas");
    clrBtn->Connect("Clicked()", "GammaFitGUI", this, "OnClearPeaks()");

    TGTextButton* rmBtn = new TGTextButton(pkRow, "Remove");
    pkRow->AddFrame(rmBtn, new TGLayoutHints(kLHintsRight, 0, 4, 0, 0));
    rmBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRemovePeak()");
    rmBtn->SetToolTipText("Remove the selected peak from the list");

    peakListBox_ = new TGListBox(pkGrp, 400);
    peakListBox_->Resize(285, 62);
    pkGrp->AddFrame(peakListBox_, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));

    // ── Fit & Statistics (merged: actions + fit selector + residuals + options) ─
    TGGroupFrame* statsGrp = new TGGroupFrame(p, "Fit & Statistics");
    p->AddFrame(statsGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    // Preview / Run Fit buttons
    {
        TGHorizontalFrame* fitBtnRow = new TGHorizontalFrame(statsGrp);
        statsGrp->AddFrame(fitBtnRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 2));
        TGTextButton* prevBtn2 = new TGTextButton(fitBtnRow, "  Preview  ");
        fitBtnRow->AddFrame(prevBtn2, new TGLayoutHints(kLHintsLeft, 2, 6, 0, 0));
        prevBtn2->Connect("Clicked()", "GammaFitGUI", this, "OnPreview()");
        prevBtn2->SetToolTipText("Draw a Gaussian at the current parameters without running MIGRAD");
        TGTextButton* fitBtn2 = new TGTextButton(fitBtnRow, "  Run Fit  ");
        fitBtnRow->AddFrame(fitBtn2, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
        fitBtn2->Connect("Clicked()", "GammaFitGUI", this, "OnManualFit()");
        fitBtn2->SetToolTipText("Run MIGRAD on the fit window defined by Energy +/- Range*Sigma");
    }

    // Show residuals checkbox
    residualChk_ = new TGCheckButton(statsGrp, "Show residuals  (data-fit)/sigma");
    statsGrp->AddFrame(residualChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
    residualChk_->Connect("Clicked()", "GammaFitGUI", this, "OnToggleResiduals()");
    residualChk_->SetToolTipText("Split the canvas and show (data-fit)/sigma below the spectrum");

    // Fit statistics text view
    peakStatsView_ = new TGTextView(statsGrp, 285, 155);
    statsGrp->AddFrame(peakStatsView_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    // Fit options checkboxes
    mFitLogLikChk_ = new TGCheckButton(statsGrp, "Log-likelihood  (L)");
    statsGrp->AddFrame(mFitLogLikChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 0));
    mFitLogLikChk_->SetToolTipText("Use Poisson log-likelihood instead of chi2 (better for low-count bins)");

    mFitImprovChk_ = new TGCheckButton(statsGrp, "IMPROVE  (M)");
    statsGrp->AddFrame(mFitImprovChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));
    mFitImprovChk_->SetToolTipText("Run IMPROVE after MIGRAD to search for a better minimum");

    mFitMinosChk_ = new TGCheckButton(statsGrp, "MINOS errors  (E)");
    statsGrp->AddFrame(mFitMinosChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));
    mFitMinosChk_->SetToolTipText("Compute asymmetric MINOS errors (slower, more accurate near parameter boundaries)");

    mShowCompChk_ = new TGCheckButton(statsGrp, "Show fit components  (BG + Gaussians)");
    statsGrp->AddFrame(mShowCompChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 2));
    mShowCompChk_->SetToolTipText("Overlay the background (green dashed) and individual Gaussian components (blue dashed) when fit is drawn");
    mShowCompChk_->Connect("Toggled(Bool_t)", "GammaFitGUI", this, "OnShowCompToggled()");

    showIsoLabelsChk_ = new TGCheckButton(statsGrp, "Show isotope matches on peaks");
    showIsoLabelsChk_->SetState(kButtonDown);
    statsGrp->AddFrame(showIsoLabelsChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 4));
    showIsoLabelsChk_->SetToolTipText(
        "Show the matched isotope name above each peak energy label.\n"
        "Uncheck to display energy values only.");
    showIsoLabelsChk_->Connect("Clicked()", "GammaFitGUI", this, "OnToggleIsoLabels()");

    mBgQuadChk_ = new TGCheckButton(statsGrp, "Quadratic background  (B0 + B1·x + B2·x²)");
    statsGrp->AddFrame(mBgQuadChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));
    mBgQuadChk_->SetToolTipText(
        "Add a quadratic (parabolic) background term B2·x².\n"
        "Useful when the Compton continuum has visible curvature under the peak.");

    mComptonStepChk_ = new TGCheckButton(statsGrp, "Compton step  (Erfc term per peak)");
    statsGrp->AddFrame(mComptonStepChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));
    mComptonStepChk_->SetToolTipText(
        "Add an Erfc step function for each Gaussian: step_i·erfc((x-E_i)/(sig_i·√2)).\n"
        "Models photons scattered into the continuum below the full-energy peak.\n"
        "Ref: Helmer & McCullagh, Nucl. Instrum. Methods 168 (1979) 593.");

    mTieWidthsChk_ = new TGCheckButton(statsGrp, "Tie widths to resolution model");
    statsGrp->AddFrame(mTieWidthsChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 4));
    mTieWidthsChk_->SetToolTipText(
        "Fix each Gaussian sigma to the FWHM-model prediction σ(E).\n"
        "Enforces physically consistent detector resolution for doublet fits.\n"
        "Requires a valid FWHM model — load from FWHM tab first.");

    // Peak label + label-pick mode
    {
        TGHorizontalFrame* lblRow = new TGHorizontalFrame(statsGrp);
        statsGrp->AddFrame(lblRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
        lblRow->AddFrame(new TGLabel(lblRow, "Label:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        mPeakLabelCombo_ = new TGComboBox(lblRow, 921);
        mPeakLabelCombo_->AddEntry("(none)", 1);
        mPeakLabelCombo_->Select(1, kFALSE);
        mPeakLabelCombo_->Resize(200, 22);
        lblRow->AddFrame(mPeakLabelCombo_, new TGLayoutHints(kLHintsExpandX));
    }
    {
        TGHorizontalFrame* pickRow = new TGHorizontalFrame(statsGrp);
        statsGrp->AddFrame(pickRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
        labelPickChk_ = new TGCheckButton(pickRow, "Label Pick Mode");
        labelPickChk_->SetToolTipText(
            "When checked, clicking on the canvas selects the nearest cached peak.\n"
            "Change its label in the combo above then click Apply Label.");
        pickRow->AddFrame(labelPickChk_, new TGLayoutHints(kLHintsCenterY, 0, 8, 0, 0));
        applyLabelBtn_ = new TGTextButton(pickRow, "Apply Label");
        applyLabelBtn_->Connect("Clicked()", "GammaFitGUI", this, "OnApplyPickedLabel()");
        applyLabelBtn_->SetToolTipText("Write the selected label to the picked peak's cache entry");
        pickRow->AddFrame(applyLabelBtn_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        labelPickInfo_ = new TGLabel(pickRow, "(no peak picked)");
        pickRow->AddFrame(labelPickInfo_, new TGLayoutHints(kLHintsCenterY | kLHintsLeft));
    }

    // Accept / Reject / Load Cache / Parameter Scan
    {
        TGHorizontalFrame* accRow = new TGHorizontalFrame(statsGrp);
        statsGrp->AddFrame(accRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        TGTextButton* accBtn = new TGTextButton(accRow, "Accept & Save");
        accRow->AddFrame(accBtn, new TGLayoutHints(kLHintsExpandX, 2, 4, 0, 0));
        accBtn->Connect("Clicked()", "GammaFitGUI", this, "OnAcceptFit()");
        accBtn->SetToolTipText("Write the current fit parameters to this histogram's cache file");
        TGTextButton* rejBtn = new TGTextButton(accRow, "Reject");
        accRow->AddFrame(rejBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
        rejBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRejectFit()");
        rejBtn->SetToolTipText("Discard the current manual fit and redraw the histogram");
    }

    {
        TGTextButton* loadCacheBtn = new TGTextButton(statsGrp, "Load Cache onto Histogram");
        statsGrp->AddFrame(loadCacheBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
        loadCacheBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLoadCache()");
        loadCacheBtn->SetToolTipText("Reload the cache from disk and overlay all stored fits on the current histogram");

        TGTextButton* scanBtn = new TGTextButton(statsGrp, "Parameter Scan Plot");
        statsGrp->AddFrame(scanBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
        scanBtn->Connect("Clicked()", "GammaFitGUI", this, "OnParameterScan()");
        scanBtn->SetToolTipText("Open a popup canvas: chi2/ndf vs each fit parameter scanned around the best-fit value");
    }
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

    AddToRecentFiles(inputPath_);

    histNames_.clear();
    th2Names_.clear();
    projParent_.clear();
    histClass_.clear();
    customProjDefs_.clear();
    bgSubtractDefs_.clear();
    rebinFactors_.clear();
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
    PopulateBgSubCombos();

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
        bool isBgSub = bgSubtractDefs_.count(name) > 0;
        if (isTH2) {
            display = "[2D] " + name;
        } else if (isCProj) {
            display = "[Custom] " + name;
        } else if (isProj) {
            auto& pp = projParent_.at(name);
            bool isX = name.size() >= 3 && name.substr(name.size()-3) == "_px";
            display = (isX ? "[ProjX] " : "[ProjY] ") + pp;
        } else if (isBgSub) {
            display = "[BgSub] " + name;
        } else {
            std::string cls = ClassOf(name);
            if      (cls == "Decay")      display = "[Decay] " + name;
            else if (cls == "Background") display = "[BG] " + name;
        }
        // Prepend rebin label when factor > 1
        {
            auto rbit = rebinFactors_.find(name);
            if (rbit != rebinFactors_.end() && rbit->second > 1)
                display = "[Rb:" + std::to_string(rbit->second) + "] " + display;
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
    schematicDrawn_ = false;  // new histogram — don't auto-redraw schematic

    if (rawHistOwned_) { delete rawHist_; rawHist_ = nullptr; rawHistOwned_ = false; }
    rawHist_ = LoadHistFromFile(currentHist_, rawHistOwned_);

    // Sync classification combo
    if (histClassCombo_) {
        std::string cls = ClassOf(currentHist_);
        if      (cls == "Decay")      histClassCombo_->Select(2, kFALSE);
        else if (cls == "2D")         histClassCombo_->Select(3, kFALSE);
        else if (cls == "Background") histClassCombo_->Select(4, kFALSE);
        else                          histClassCombo_->Select(1, kFALSE);
    }

    // Sync rebin entry to stored factor for this histogram
    if (rebinEntry_) {
        auto rbit = rebinFactors_.find(currentHist_);
        rebinEntry_->SetNumber(rbit != rebinFactors_.end() ? rbit->second : 1);
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
    if (mResultLbl_) mResultLbl_->SetText("No fit yet");
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
        bool showComp = mShowCompChk_ && mShowCompChk_->IsDown();
        manualTF1_->SetLineColor(kRed);
        manualTF1_->SetLineWidth(2);
        DrawFitComponents(c, manualTF1_);
        if (!showComp)
            manualTF1_->Draw("same");
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

std::string GammaFitGUI::CacheDirFor() const
{
    if (inputPath_.empty()) return std::string(kCacheDir);
    std::string fname = inputPath_;
    auto slash = fname.find_last_of("/\\");
    if (slash != std::string::npos)
        fname = fname.substr(slash + 1);
    for (char& c : fname)
        if (!std::isalnum((unsigned char)c) && c != '.') c = '_';
    return std::string(kCacheDir) + "/" + fname;
}

std::string GammaFitGUI::ArchiveDirFor() const
{
    return CacheDirFor() + "/archive";
}

void GammaFitGUI::EnsureCacheDir() const
{
    ::mkdir(kCacheDir, 0755);
    ::mkdir(CacheDirFor().c_str(), 0755);
}

std::string GammaFitGUI::CacheFileFor(const std::string& hname) const
{
    return CacheDirFor() + "/fit_cache_" + hname + ".dat";
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

    TH1* result = nullptr;

    // Background-subtracted virtual histogram
    auto bsit = bgSubtractDefs_.find(hname);
    if (bsit != bgSubtractDefs_.end()) {
        const BgSubtractDef& bsd = bsit->second;
        bool srcOwned = false, bgOwned = false;
        TH1* src = LoadHistFromFile(bsd.srcName, srcOwned);
        TH1* bg  = LoadHistFromFile(bsd.bgName,  bgOwned);
        if (src && bg) {
            result = (TH1*)src->Clone(hname.c_str());
            result->SetDirectory(nullptr);
            result->Sumw2();
            result->Add(bg, -bsd.scale);
            owned = true;
        }
        if (srcOwned) delete src;
        if (bgOwned)  delete bg;
    }
    // Custom projection with user-defined cut range
    else {
        auto cit = customProjDefs_.find(hname);
        if (cit != customProjDefs_.end()) {
            const CustomProjDef& def = cit->second;
            TH2* h2 = (TH2*)inputFile_->Get(def.th2Name.c_str());
            if (h2) {
                if (def.projX) {
                    int b1 = h2->GetYaxis()->FindBin(def.lo);
                    int b2 = h2->GetYaxis()->FindBin(def.hi);
                    result = h2->ProjectionX(hname.c_str(), b1, b2);
                } else {
                    int b1 = h2->GetXaxis()->FindBin(def.lo);
                    int b2 = h2->GetXaxis()->FindBin(def.hi);
                    result = h2->ProjectionY(hname.c_str(), b1, b2);
                }
                if (result) { result->SetDirectory(nullptr); owned = true; }
            }
        }
        // Auto projection (full axis)
        else if (projParent_.find(hname) != projParent_.end()) {
            result = LoadProjection(inputFile_, hname, projParent_, nullptr, nullptr);
            owned = (result != nullptr);
        }
        // Plain histogram from file
        else {
            result = (TH1*)inputFile_->Get(hname.c_str());
        }
    }

    // Apply stored rebin factor (1D only)
    if (result && !result->InheritsFrom("TH2")) {
        auto rbit = rebinFactors_.find(hname);
        if (rbit != rebinFactors_.end() && rbit->second > 1) {
            int n = rbit->second;
            if (!owned) {
                result = (TH1*)result->Clone((hname + "_rb").c_str());
                result->SetDirectory(nullptr);
                owned = true;
            }
            // Sumw2 before Rebin ensures Poisson errors add in quadrature
            result->Sumw2();
            result->Rebin(n);
        }
    }

    return result;
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
    EnsureCacheDir();

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
        bgOpts.snMinRatio       = snThreshEntry_  ? snThreshEntry_->GetNumber() : 0.0;
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
                EnsureCacheDir();
                edb.Save(CacheFileFor(fwhmHistName_));
            }
        }
        return;
    }

    if (!rawHist_) { AppendLog("Load a histogram first (Manual Fit tab -> Load to Canvas)."); return; }

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

    // ── Label pick mode — click to select nearest cached Gaussian ────────────
    if (labelPickChk_ && labelPickChk_->IsOn()) {
        FitDatabase fitdb;
        if (!fitdb.Load(CacheFileFor(currentHist_))) {
            AppendLog("No cache for this histogram — run a fit first.");
            return;
        }
        double bestDist = 1e9;
        std::string bestKey;
        int bestGi = -1;
        for (const auto& kv : fitdb.GetEntries()) {
            if (kv.first.empty() || kv.first[0] == '_') continue;
            const FitEntry& fe = kv.second;
            int npar = (int)fe.params.size();
            if (npar < 5 || (npar-2)%3 != 0) continue;
            int n = (npar-2)/3;
            for (int gi = 0; gi < n; gi++) {
                double E = fe.params[3*gi+1];
                double dist = std::fabs(E - energy);
                if (dist < bestDist) { bestDist = dist; bestKey = kv.first; bestGi = (n>1) ? gi : -1; }
            }
        }
        if (bestKey.empty() || bestDist > 20.0) {
            AppendLog(Form("No cached peak within 20 keV of %.1f — click closer to a peak.", energy));
            return;
        }
        labelPickKey_      = bestKey;
        labelPickGaussIdx_ = bestGi;
        const FitEntry& picked = fitdb.GetEntries().at(bestKey);
        std::string curLabel = (bestGi >= 0) ? picked.PeakLabel(bestGi) : picked.label;
        double pickedE = picked.params[3*(bestGi>=0?bestGi:0)+1];

        // Populate combo and pre-select current label
        if (mPeakLabelCombo_ && dbLoaded_) {
            mPeakLabelCombo_->RemoveAll();
            mPeakLabelCombo_->AddEntry("(none)", 1);
            int cid = 2;
            if (!curLabel.empty()) mPeakLabelCombo_->AddEntry(curLabel.c_str(), cid++);
            for (const auto& gl : db_.db) {
                if (gl.isotope == curLabel) continue;
                mPeakLabelCombo_->AddEntry(gl.isotope.c_str(), cid++);
            }
            mPeakLabelCombo_->Select(curLabel.empty() ? 1 : 2, kFALSE);
            mPeakLabelCombo_->MapSubwindows(); mPeakLabelCombo_->Layout();
        }
        std::string info = Form("%.1f keV", pickedE);
        if (!curLabel.empty()) info += " [" + curLabel + "]";
        if (labelPickInfo_) labelPickInfo_->SetText(info.c_str());
        AppendLog("Label pick: " + info + (bestGi>=0 ? Form("  Gaussian %d", bestGi+1) : ""));
        return;
    }

    // ── Peak placement ────────────────────────────────────────────────────────
    clickEnergy_ = energy;
    mEnergy_->SetNumber(energy);

    if (!addPeakChk_->IsOn()) {
        manualPeaks_.clear();
        peakListBox_->RemoveAll();
        peakListBox_->MapSubwindows(); peakListBox_->Layout();
        if (rawHist_) {
            RedrawView();
            ylo = targetPad->GetUymin();
            yhi = targetPad->GetUymax() * 0.95;
        }
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

    bool previewQuad  = mBgQuadChk_      && mBgQuadChk_->IsDown();
    bool previewStep  = mComptonStepChk_ && mComptonStepChk_->IsDown();

    delete manualTF1_;
    for (TF1* fc : fitComponents_) delete fc;
    fitComponents_.clear();
    manualTF1_ = new TF1("manual_preview",
                         BuildNGaussFormulaEx(n, previewQuad, previewStep).c_str(),
                         xmin, xmax);

    for (int i = 0; i < n; i++) {
        double E   = manualPeaks_[i];
        double sig = res_.Sigma(E);
        double A   = std::max(rawHist_->GetBinContent(rawHist_->FindBin(E)), 1.0);
        manualTF1_->SetParameter(3*i,   A);
        manualTF1_->SetParameter(3*i+1, E);
        manualTF1_->SetParameter(3*i+2, sig);
    }
    manualTF1_->SetParameter(3*n,   bg0);
    if (mBgFlatChk_ && mBgFlatChk_->IsDown())
        manualTF1_->SetParameter(3*n+1, 0.0);
    else
        manualTF1_->SetParameter(3*n+1, bg1);
    if (previewQuad)
        manualTF1_->SetParameter(3*n+2, 0.0);
    if (previewStep) {
        int nbg = NBgPars(previewQuad);
        for (int i = 0; i < n; i++) {
            double A = std::max(rawHist_->GetBinContent(rawHist_->FindBin(manualPeaks_[i])), 1.0);
            manualTF1_->SetParameter(3*n + nbg + i, 0.05 * A);
        }
    }
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

    bool useQuadBG    = mBgQuadChk_      && mBgQuadChk_->IsDown();
    bool useCompStep  = mComptonStepChk_ && mComptonStepChk_->IsDown();
    bool useTieWidths = mTieWidthsChk_   && mTieWidthsChk_->IsDown();

    delete manualTF1_;
    for (TF1* fc : fitComponents_) delete fc;
    fitComponents_.clear();
    manualTF1_ = new TF1("manual_fit",
                         BuildNGaussFormulaEx(n, useQuadBG, useCompStep).c_str(),
                         xmin, xmax);

    double ampLo  = mAmpLoFrac_  ? mAmpLoFrac_ ->GetNumber() : 0.01;
    double ampHi  = mAmpHiFrac_  ? mAmpHiFrac_ ->GetNumber() : 20.0;
    double sigLo  = mSigLoFrac_  ? mSigLoFrac_ ->GetNumber() : 0.2;
    double sigHi  = mSigHiFrac_  ? mSigHiFrac_ ->GetNumber() : 4.0;
    double eWin   = mEnergyWin_  ? mEnergyWin_ ->GetNumber() : 8.0;

    for (int i = 0; i < n; i++) {
        double E        = manualPeaks_[i];
        double sigModel = res_.Sigma(E);
        double sig      = sigModel;
        double A        = std::max(fitHist->GetBinContent(fitHist->FindBin(E)), 1.0);
        manualTF1_->SetParName(3*i,   Form("A_%d",   i+1));
        manualTF1_->SetParName(3*i+1, Form("E_%d",   i+1));
        manualTF1_->SetParName(3*i+2, Form("sig_%d", i+1));

        manualTF1_->SetParameter(3*i,     A);
        manualTF1_->SetParLimits(3*i,     std::max(A * ampLo, 1.0), A * ampHi);
        manualTF1_->SetParameter(3*i+1,   E);
        manualTF1_->SetParLimits(3*i+1,   E - eWin, E + eWin);
        if (useTieWidths) {
            manualTF1_->SetParameter(3*i+2, sigModel);
            manualTF1_->FixParameter(3*i+2, sigModel);
        } else {
            manualTF1_->SetParameter(3*i+2,   sig);
            manualTF1_->SetParLimits(3*i+2,   sigModel * sigLo, sigModel * sigHi);
        }
    }
    manualTF1_->SetParName(3*n,   "bg0");
    manualTF1_->SetParName(3*n+1, "bg1");
    manualTF1_->SetParameter(3*n,   bg0);
    manualTF1_->SetParameter(3*n+1, bg1);
    if (mBgFlatChk_ && mBgFlatChk_->IsDown())
        manualTF1_->FixParameter(3*n+1, 0.0);
    if (useQuadBG) {
        manualTF1_->SetParName(3*n+2, "bg2");
        manualTF1_->SetParameter(3*n+2, 0.0);
        double bg2Lim = (xmax > xmin)
            ? std::max(std::abs(bg0) / ((xmax-xmin)*(xmax-xmin)), 1e-8)
            : 1e-6;
        manualTF1_->SetParLimits(3*n+2, -bg2Lim*10.0, bg2Lim*10.0);
    }
    if (useCompStep) {
        int nbg = NBgPars(useQuadBG);
        for (int i = 0; i < n; i++) {
            int sIdx = 3*n + nbg + i;
            double A = manualTF1_->GetParameter(3*i);
            manualTF1_->SetParName(sIdx, Form("step_%d", i+1));
            manualTF1_->SetParameter(sIdx, 0.05 * A);
            manualTF1_->SetParLimits(sIdx, 0.0, A);
        }
    }

    // Snapshot old stats before fitting so UpdatePeakStats can show old vs new
    peakStatsOld_ = peakStatsCurrent_;

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
    lastManualEdm_ = (r.Get() && r->IsValid()) ? r->Edm() : -1.0;

    // Clamp background so it never goes negative over the fit window
    {
        double bg0   = manualTF1_->GetParameter(3*n);
        double bg1   = manualTF1_->GetParameter(3*n+1);
        double bgMin = std::min(bg0 + bg1 * xmin, bg0 + bg1 * xmax);
        if (bgMin < 0.0)
            manualTF1_->SetParameter(3*n, bg0 - bgMin);
    }

    // Check whether any parameter hit its bound — zero errors are a symptom
    {
        int npar = manualTF1_->GetNpar();
        for (int i = 0; i < npar; i++) {
            double val = manualTF1_->GetParameter(i);
            double lo, hi;
            manualTF1_->GetParLimits(i, lo, hi);
            if (hi <= lo) continue;  // fixed or no limits
            double range = hi - lo;
            const char* pname = manualTF1_->GetParName(i);
            if (std::abs(val - lo) / range < 0.002)
                AppendLog(Form("  WARNING: %s hit LOWER bound  (val=%.5g  bound=%.5g)",
                               pname, val, lo));
            else if (std::abs(val - hi) / range < 0.002)
                AppendLog(Form("  WARNING: %s hit UPPER bound  (val=%.5g  bound=%.5g)",
                               pname, val, hi));
        }
    }

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

    if (mResultLbl_) mResultLbl_->SetText(result.c_str());

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

    // Build canvas title with EDM if available
    std::string scanTitle = Form("Parameter Scan  #chi^{2}/ndf = %.3f", chi2best);
    if (lastManualEdm_ >= 0.0)
        scanTitle += Form("   EDM = %.3g", lastManualEdm_);
    scanTitle += Form("  [%s]", currentHist_.c_str());

    TCanvas* scanC = new TCanvas("param_scan", scanTitle.c_str(), cw, ch);
    scanC->Divide(ncols, nrows);

    // Helper: convert parameter name to ROOT TLatex notation
    auto parLatex = [](const char* nm) -> std::string {
        std::string s = nm;
        if (s == "bg0") return "b_{0}";
        if (s == "bg1") return "b_{1}";
        if (s.size() > 4 && s.substr(0,4) == "sig_") return "#sigma_{" + s.substr(4) + "}";
        if (s.size() > 2 && s[0] == 'A' && s[1] == '_') return "A_{" + s.substr(2) + "}";
        if (s.size() > 2 && s[0] == 'E' && s[1] == '_') return "E_{" + s.substr(2) + "} (keV)";
        return s;
    };

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
        std::string pLatex = parLatex(pname);
        gr->SetTitle(Form(";%s;#chi^{2}/ndf", pLatex.c_str()));
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

        // Best-fit value label in top-left corner of each pad
        TLatex* lbl = new TLatex(0.21, 0.88,
            Form("%s = %.4g #pm %.3g", pLatex.c_str(), val0, err0));
        lbl->SetNDC();
        lbl->SetTextSize(0.055);
        lbl->SetTextColor(kRed);
        lbl->Draw();

        // EDM annotation in the top-right of the first pad only
        if (p == 0 && lastManualEdm_ >= 0.0) {
            TLatex* edmLbl = new TLatex(0.96, 0.88,
                Form("EDM = %.3g", lastManualEdm_));
            edmLbl->SetNDC();
            edmLbl->SetTextSize(0.055);
            edmLbl->SetTextAlign(31);
            edmLbl->SetTextColor(kGray + 2);
            edmLbl->Draw();
        }
    }

    scanC->Modified();
    scanC->Update();
    std::string edmStr = (lastManualEdm_ >= 0.0) ? Form("  EDM=%.3g", lastManualEdm_) : "";
    AppendLog("Parameter scan: " + std::to_string(npar) +
              " params  chi2/ndf=" + Fmt(chi2best, 3) + edmStr);
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
    if (fitHist) {
        int midBin = fitHist->FindBin(0.5 * (xlo + xhi));
        double bw  = fitHist->GetBinWidth(midBin);
        if (bw > 0.0) e.binWidth = bw;
    }

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
    // Auto-populate classification from labelClassMap_ (set via Isotopes tab)
    if (!e.label.empty() && labelClassMap_.count(e.label))
        e.classification = labelClassMap_.at(e.label);
    if (e.classification.empty()) {
        // X-ray fallback for low-energy unclassified peaks
        int nPeaks2 = ((int)e.params.size() - 2) / 3;
        for (int pi = 0; pi < nPeaks2; pi++) {
            double Efit = e.params[3 * pi + 1];
            if (Efit > 0 && Efit < 100.0) { e.classification = "X-ray"; break; }
        }
    }

    // Force-store: manual fits always win regardless of score comparison.
    fitdb.ForceStore(key, e);
    fitdb.rootFile = inputPath_;
    EnsureCacheDir();
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
              "  chi2/ndf~" + Fmt(e.chi2ndf, 2));
    SetStatus("Saved: " + currentHist_ + "  key=" + key);
}

void GammaFitGUI::OnRejectFit()
{
    delete manualTF1_; manualTF1_ = nullptr;
    delete bgTF1_;     bgTF1_     = nullptr;
    for (TF1* fc : fitComponents_) delete fc;
    fitComponents_.clear();
    if (mResultLbl_) mResultLbl_->SetText("Fit rejected.");
    OnHistViewChanged(histViewCombo_->GetSelected());
    AppendLog("Manual fit rejected — canvas cleared.");
}

void GammaFitGUI::OnShowCompToggled()
{
    // Immediate canvas refresh so checking/unchecking takes effect right away
    if (!rawHist_) return;
    OnHistViewChanged(histViewCombo_ ? histViewCombo_->GetSelected() : 1);
}

// ─────────────────────────────────────────────────────────────────────────────
// Multi-peak / background helpers
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnClearPeaks()
{
    manualPeaks_.clear();
    peakListBox_->RemoveAll();
    peakListBox_->MapSubwindows(); peakListBox_->Layout();
    if (rawHist_) RedrawView();
    AppendLog("Peak list and click markers cleared.");
}

// Lightweight redraw — reuses existing viewHist_ rather than rerunning MakeBgSubHist.
// Use this instead of OnHistViewChanged when only click markers need clearing.
void GammaFitGUI::RedrawView()
{
    if (!rawHist_) return;
    TCanvas* c = canvas_->GetCanvas();
    c->Clear();
    c->cd();

    TH1* h = viewHist_ ? viewHist_ : rawHist_;
    ApplyHistStyle(h, currentHist_.empty() ? nullptr : currentHist_.c_str());
    if (!currentHist_.empty() && ClassOf(currentHist_) == "Decay")
        h->GetXaxis()->SetTitle("Time (ms)");
    h->SetLineColor(kBlack);
    h->SetMarkerSize(0);
    h->SetStats(0);
    if (viewXmin_ < viewXmax_)
        h->GetXaxis()->SetRangeUser(viewXmin_, viewXmax_);
    else
        h->GetXaxis()->UnZoom();
    SetYMaxFromVisible(h);
    h->Draw("hist");
    h->Draw("E1 same");

    Int_t mode = histViewCombo_ ? histViewCombo_->GetSelected() : 1;
    if (mode == 3 && !currentHist_.empty())
        OverlayFitPeaks(currentHist_, c);

    if (manualTF1_) {
        bool showComp = mShowCompChk_ && mShowCompChk_->IsDown();
        manualTF1_->SetLineColor(kRed);
        manualTF1_->SetLineWidth(2);
        DrawFitComponents(c, manualTF1_);
        if (!showComp)
            manualTF1_->Draw("same");
        DrawPeakLabels(manualTF1_);
    }

    // Integral statistics overlay
    {
        TPaveText* pt = new TPaveText(0.55, 0.75, 0.98, 0.98, "NDC");
        pt->SetBorderSize(1);
        pt->SetFillColor(kWhite);
        pt->SetFillStyle(1001);
        pt->SetTextAlign(12);
        pt->SetTextSize(0.030);
        pt->SetBit(kCanDelete);
        pt->AddText(Form("Integral = %.0f", h->Integral()));
        pt->AddText(Form("Entries  = %.0f", h->GetEntries()));
        pt->AddText(Form("Mean     = %.2f keV", h->GetMean()));
        pt->AddText(Form("RMS      = %.2f keV", h->GetRMS()));
        auto rbit = rebinFactors_.find(currentHist_);
        if (rbit != rebinFactors_.end() && rbit->second > 1)
            pt->AddText(Form("Rebin x%d", rbit->second));
        pt->Draw();
    }

    c->Modified(); c->Update();
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

// ─────────────────────────────────────────────────────────────────────────────
// Staggered peak labels — shared by OverlayFitPeaks and DrawPeakLabels
// ─────────────────────────────────────────────────────────────────────────────
namespace {
struct PkLblData {
    double E;
    double peakY;       // data y at peak (fit or histogram)
    std::string iso;    // isotope name, may be empty
    std::string ekeV;   // energy string e.g. "1234.5"
};

void DrawStaggedLabels(const std::vector<PkLblData>& in, TVirtualPad* pad)
{
    if (in.empty() || !pad) return;

    std::vector<const PkLblData*> pts;
    for (auto& p : in) pts.push_back(&p);
    std::sort(pts.begin(), pts.end(), [](auto* a, auto* b){ return a->E < b->E; });

    double xlo   = pad->GetUxmin(), xhi = pad->GetUxmax();
    double xSpan = xhi - xlo;
    if (xSpan <= 0) return;

    double ylo = pad->GetUymin(), yhi = pad->GetUymax();
    bool   logY = (pad->GetLogy() != 0);

    // Peak y → NDC fraction (0 = bottom, 1 = top of pad)
    auto toNDC = [&](double dataY) -> double {
        if (logY) {
            if (dataY <= 0) return 0.0;
            double lmin = (ylo > 0) ? std::log10(ylo) : -1.0;
            double lmax = std::log10(yhi);
            return (std::log10(dataY) - lmin) / (lmax - lmin);
        }
        return (yhi > ylo) ? (dataY - ylo) / (yhi - ylo) : 0.0;
    };

    // NDC fraction → data y (for drawing TLatex / TLine)
    auto fromNDC = [&](double f) -> double {
        f = std::max(0.0, std::min(f, 0.999));
        if (logY) {
            double lmin = (ylo > 0) ? std::log10(ylo) : -1.0;
            double lmax = std::log10(yhi);
            return std::pow(10.0, lmin + f * (lmax - lmin));
        }
        return ylo + f * (yhi - ylo);
    };

    // NDC constants
    // kRowH: vertical space (NDC) reserved for one text row (textSize + leading).
    // kMinGap: clearance from peak top to bottom of energy label.
    // kXConf: NDC x-fraction within which two labels can conflict vertically.
    const double kRowH   = 0.032;
    const double kMinGap = 0.018;
    const double kXConf  = 0.10;
    const double kMaxY   = 0.96;

    // Greedy y-assignment: process left-to-right, tracking placed label extents.
    struct Placed { double xNDC, yBot, yTop; };
    std::vector<Placed>  placed;
    std::vector<double>  yBotVec(pts.size(), -1.0);

    for (int i = 0; i < (int)pts.size(); i++) {
        double E = pts[i]->E;
        if (E < xlo || E > xhi) continue;

        double xNDC  = (E - xlo) / xSpan;
        double pNDC  = toNDC(pts[i]->peakY);
        bool   hasIso = !pts[i]->iso.empty();
        double lblH  = hasIso ? 2.0 * kRowH : kRowH;

        // Start the label just above the peak
        double yBot = pNDC + kMinGap;

        // Push up until no overlap with any nearby placed label
        bool changed = true;
        while (changed) {
            changed = false;
            for (auto& pl : placed) {
                if (std::abs(pl.xNDC - xNDC) > kXConf) continue;
                double myTop = yBot + lblH;
                if (myTop > pl.yBot && yBot < pl.yTop) {
                    yBot    = pl.yTop + kMinGap * 0.4;
                    changed = true;
                }
            }
        }

        // Clamp to pad top
        if (yBot + lblH > kMaxY) yBot = kMaxY - lblH;
        yBot = std::max(yBot, kMinGap);

        yBotVec[i] = yBot;
        placed.push_back({xNDC, yBot, yBot + lblH});
    }

    TLatex lbl;
    lbl.SetTextSize(0.022);
    lbl.SetTextAlign(21);  // centre-bottom

    for (int i = 0; i < (int)pts.size(); i++) {
        double E = pts[i]->E;
        if (E < xlo || E > xhi || yBotVec[i] < 0) continue;

        double yBot   = yBotVec[i];
        double peakTp = pts[i]->peakY;

        // Dotted stem from just above peak top to just below energy label
        if (peakTp > 0) {
            double stemBot = fromNDC(toNDC(peakTp) + 0.004);
            double stemTop = fromNDC(yBot - 0.004);
            if (stemBot < stemTop) {
                TLine* stem = new TLine(E, stemBot, E, stemTop);
                stem->SetLineColor(kGray + 1);
                stem->SetLineWidth(1);
                stem->SetLineStyle(3);
                stem->SetBit(kCanDelete);
                stem->Draw("same");
            }
        }

        // Energy label (bottom row)
        lbl.SetTextColor(kBlack);
        lbl.DrawLatex(E, fromNDC(yBot), pts[i]->ekeV.c_str());

        // Isotope label (row above)
        if (!pts[i]->iso.empty()) {
            lbl.SetTextColor(kAzure + 2);
            lbl.DrawLatex(E, fromNDC(yBot + kRowH), pts[i]->iso.c_str());
        }
    }
}
} // namespace

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

    // Global set of already-labeled energies — prevents two labels within 1 FWHM
    // of each other regardless of which cache group they come from.
    std::vector<double> labeledEs;
    std::vector<PkLblData> pkLabels;

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

            // Draw background component (always visible)
            int npar = f->GetNpar();
            if (npar >= 5 && (npar - 2) % 3 == 0) {
                int np = (npar - 2) / 3;
                double bg0 = f->GetParameter(3*np);
                double bg1 = f->GetParameter(3*np + 1);
                TF1* bgf = new TF1(Form("ofp_bg_%p", (void*)f), "[0]+[1]*x", xlo, xhi);
                bgf->SetParameter(0, bg0);
                bgf->SetParameter(1, bg1);
                bgf->SetLineColor(kGreen + 2);
                bgf->SetLineStyle(2);
                bgf->SetLineWidth(2);
                bgf->SetBit(kCanDelete);
                bgf->Draw("same");

                // Individual Gaussian components when option is on and n > 1
                if (np > 1 && mShowCompChk_ && mShowCompChk_->IsDown()) {
                    for (int gi = 0; gi < np; gi++) {
                        double A   = f->GetParameter(3*gi);
                        double E   = f->GetParameter(3*gi + 1);
                        double sig = f->GetParameter(3*gi + 2);
                        TF1* gf = new TF1(Form("ofp_g%d_%p", gi, (void*)f),
                                          "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x", xlo, xhi);
                        gf->SetParameter(0, A);
                        gf->SetParameter(1, E);
                        gf->SetParameter(2, sig);
                        gf->SetParameter(3, bg0);
                        gf->SetParameter(4, bg1);
                        gf->SetLineColor(kBlue + 1);
                        gf->SetLineStyle(2);
                        gf->SetLineWidth(1);
                        gf->SetBit(kCanDelete);
                        gf->Draw("same");
                    }
                }
            }
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

        for (int gi = 0; gi < (int)compEs.size(); gi++) {
            double E = compEs[gi];
            double yLabel = f ? f->Eval(E)
                              : (dispH ? dispH->GetBinContent(dispH->FindBin(E)) : 0.0);
            if (yLabel <= 0.0) continue;
            if (labelTooClose(E)) continue;
            labeledEs.push_back(E);

            // Centroid uncertainty from fit parameter error (0 if no fit)
            double eErr = 0.0;
            if (f) {
                int npar = f->GetNpar();
                if (npar == 7)
                    eErr = f->GetParError(1);
                else if (npar >= 5 && (npar - 2) % 3 == 0)
                    eErr = f->GetParError(3*gi + 1);
            }

            // Prefer per-Gaussian label, then entry-level label, then DB match
            std::string isoName;
            if (showIsoLabels_) {
                if (gi < (int)cand.entry->peakLabels.size() && !cand.entry->peakLabels[gi].empty())
                    isoName = cand.entry->peakLabels[gi];
                else if (!cand.entry->label.empty())
                    isoName = cand.entry->label;
                else if (dbLoaded_) {
                    auto matches = db_.Match(E, res_.FWHM(E));
                    if (!matches.empty()) isoName = matches[0].isotope;
                }
            }

            pkLabels.push_back({E, yLabel, isoName, NNDCFormat(E, eErr)});
        }
        // f is now owned by the pad (kCanDelete) — do NOT delete here
    }
    DrawStaggedLabels(pkLabels, c->cd());
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
        bool showComp = mShowCompChk_ && mShowCompChk_->IsDown();
        manualTF1_->SetLineColor(kRed);
        manualTF1_->SetLineWidth(2);
        DrawFitComponents(c, manualTF1_);  // BG always + Gaussians if mShowCompChk_ on
        if (!showComp)
            manualTF1_->Draw("same");      // total fit in red only when components hidden
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
    AppendLog("Fit range cleared — will use peaks +/- Range*sigma.");
}

// ─────────────────────────────────────────────────────────────────────────────
// Fit Parameters popup
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnFitParamClose()
{
    if (fitParamDlg_) fitParamDlg_->UnmapWindow();
}

void GammaFitGUI::OnTransferCache()
{
    // Destination = histogram selected in the AutoFit list
    Int_t selId = histList_ ? histList_->GetSelected() : -1;
    if (selId < 1 || (size_t)selId > histNames_.size()) {
        AppendLog("Select a destination histogram in the Histograms list first.");
        return;
    }
    const std::string destHist = histNames_[selId - 1];

    // Each entry: display label + full path to the .dat file
    struct SrcEntry { std::string label; std::string path; };
    std::vector<SrcEntry> sources;

    // Helper: scan a directory for fit_cache_*.dat files
    auto scanLive = [&](const std::string& scanDir, const std::string& tag) {
        DIR* d = opendir(scanDir.c_str());
        if (!d) return;
        struct dirent* ent;
        const std::string pfx = "fit_cache_", sfx = ".dat";
        while ((ent = readdir(d)) != nullptr) {
            std::string fn = ent->d_name;
            if (fn.size() <= pfx.size() + sfx.size()) continue;
            if (fn.substr(0, pfx.size()) != pfx) continue;
            if (fn.substr(fn.size() - sfx.size()) != sfx) continue;
            std::string hname = fn.substr(pfx.size(), fn.size() - pfx.size() - sfx.size());
            if (hname == destHist) continue;
            std::string lbl = tag.empty() ? hname : hname + "  [" + tag + "]";
            sources.push_back({lbl, scanDir + "/" + fn});
        }
        closedir(d);
    };

    // Helper: scan archive dir for fit_cache_<hname>_<timestamp>.dat
    auto scanArchive = [&](const std::string& archDir) {
        DIR* d = opendir(archDir.c_str());
        if (!d) return;
        struct dirent* ent;
        const std::string pfx = "fit_cache_", sfx = ".dat";
        // Timestamp suffix is always _YYYYMMDD_HHMMSS = 16 chars
        const int tsLen = 16;
        while ((ent = readdir(d)) != nullptr) {
            std::string fn = ent->d_name;
            if ((int)fn.size() <= (int)(pfx.size() + sfx.size()) + tsLen) continue;
            if (fn.substr(0, pfx.size()) != pfx) continue;
            if (fn.substr(fn.size() - sfx.size()) != sfx) continue;
            std::string stem = fn.substr(pfx.size(), fn.size() - pfx.size() - sfx.size());
            if ((int)stem.size() <= tsLen) continue;
            std::string ts    = stem.substr(stem.size() - tsLen); // _YYYYMMDD_HHMMSS
            std::string hname = stem.substr(0, stem.size() - tsLen);
            sources.push_back({hname + "  [archived " + ts.substr(1) + "]",
                                archDir + "/" + fn});
        }
        closedir(d);
    };

    scanLive(CacheDirFor(), "");                   // current per-file live caches
    scanLive(std::string(kCacheDir), "legacy");    // old flat caches from before subdir scheme
    scanArchive(ArchiveDirFor());                  // archived caches for this file

    // De-duplicate by path
    std::sort(sources.begin(), sources.end(),
              [](const SrcEntry& a, const SrcEntry& b){ return a.label < b.label; });
    sources.erase(std::unique(sources.begin(), sources.end(),
              [](const SrcEntry& a, const SrcEntry& b){ return a.path == b.path; }),
              sources.end());

    if (sources.empty()) {
        AppendLog("[Transfer] No source caches found (live, legacy, or archived).");
        return;
    }

    // Build picker dialog
    TGTransientFrame* dlg = new TGTransientFrame(gClient->GetRoot(), this, 460, 10, kVerticalFrame);
    dlg->SetWindowName("Transfer Cache From");
    dlg->SetCleanup(kDeepCleanup);

    dlg->AddFrame(new TGLabel(dlg, Form("Destination: %s", destHist.c_str())),
                  new TGLayoutHints(kLHintsLeft, 8, 8, 8, 2));
    dlg->AddFrame(new TGLabel(dlg, "Source  (live, legacy, or archived):"),
                  new TGLayoutHints(kLHintsLeft, 8, 8, 0, 4));

    TGComboBox* srcCombo = new TGComboBox(dlg, 500);
    srcCombo->Resize(430, 22);
    dlg->AddFrame(srcCombo, new TGLayoutHints(kLHintsExpandX, 8, 8, 0, 4));
    for (int i = 0; i < (int)sources.size(); i++)
        srcCombo->AddEntry(sources[i].label.c_str(), i + 1);
    srcCombo->Select(1, kFALSE);

    TGCheckButton* mergeChk = new TGCheckButton(dlg, "Merge  (keep existing entries in destination)");
    dlg->AddFrame(mergeChk, new TGLayoutHints(kLHintsLeft, 8, 8, 4, 4));
    mergeChk->SetState(kButtonDown);

    TGHorizontalFrame* btnRow = new TGHorizontalFrame(dlg);
    dlg->AddFrame(btnRow, new TGLayoutHints(kLHintsCenterX, 8, 8, 4, 8));
    TGTextButton* okBtn  = new TGTextButton(btnRow, "  Transfer  ");
    TGTextButton* canBtn = new TGTextButton(btnRow, "  Cancel  ");
    btnRow->AddFrame(okBtn,  new TGLayoutHints(kLHintsLeft, 0, 8, 0, 0));
    btnRow->AddFrame(canBtn, new TGLayoutHints(kLHintsLeft));
    // OK sets a user bit on the dialog so we can detect it after WaitForUnmap
    okBtn->Connect("Clicked()", "TObject", dlg, "SetBit(UInt_t=16384)");
    okBtn->Connect("Clicked()", "TGFrame", dlg, "UnmapWindow()");
    canBtn->Connect("Clicked()", "TGFrame", dlg, "UnmapWindow()");

    dlg->MapSubwindows();
    dlg->Resize(dlg->GetDefaultSize());
    dlg->CenterOnParent();
    dlg->MapWindow();

    gClient->WaitForUnmap(dlg);

    Int_t srcSelId = srcCombo->GetSelected();
    bool  merge    = mergeChk->IsOn();
    bool  confirmed = dlg->TestBit(BIT(14));  // BIT(14) = 16384, set only when OK clicked
    dlg->DeleteWindow();

    if (!confirmed || srcSelId < 1 || srcSelId > (int)sources.size()) return;

    const std::string srcPath = sources[srcSelId - 1].path;
    const std::string srcLabel = sources[srcSelId - 1].label;

    // Load source cache from whatever path was selected
    FitDatabase srcDb;
    if (!srcDb.Load(srcPath)) {
        AppendLog("[Transfer] Could not load cache: " + srcPath);
        return;
    }

    // Load destination cache
    FitDatabase dstDb;
    dstDb.Load(CacheFileFor(destHist));  // ok if missing

    if (!merge) {
        // Full replace: wipe destination (keep only internal _ keys), then copy everything
        std::vector<std::string> toRemove;
        for (const auto& kv : dstDb.GetEntries())
            if (kv.first.empty() || kv.first[0] != '_')
                toRemove.push_back(kv.first);
        for (const auto& k : toRemove)
            dstDb.Remove(k);
        // Copy bg metadata from source so OnLoadCacheSelected picks it up correctly
        dstDb.bgSubtracted = srcDb.bgSubtracted;
        dstDb.bgIterations = srcDb.bgIterations;
    }

    int nCopied = 0, nSkipped = 0;
    for (const auto& kv : srcDb.GetEntries()) {
        if (!kv.first.empty() && kv.first[0] == '_') continue;  // skip internal keys
        if (merge && dstDb.GetEntries().count(kv.first)) {
            ++nSkipped;
            continue;
        }
        dstDb.ForceStore(kv.first, kv.second);
        ++nCopied;
    }

    EnsureCacheDir();
    dstDb.Save(CacheFileFor(destHist));
    AppendLog(Form("[Transfer] '%s' -> '%s': copied %d, skipped %d (merge=%s)",
                   srcLabel.c_str(), destHist.c_str(), nCopied, nSkipped,
                   merge ? "on" : "off"));

    // If the destination is currently displayed, re-apply bg settings then refresh
    if (destHist == currentHist_) {
        if (dstDb.bgSubtracted) {
            histViewCombo_->Select(2, kFALSE);
            bgSubtractChk_->SetState(kButtonDown);
            bgIterEntry_->SetNumber(dstDb.bgIterations);
        } else {
            histViewCombo_->Select(1, kFALSE);
            bgSubtractChk_->SetState(kButtonUp);
        }
        OnLoadCache();
        PopulateNavAndResidual();
    }
}

void GammaFitGUI::OnSeedBoundsFromModel()
{
    if (!mEnergyWin_ || !mSigLoFrac_ || !mSigHiFrac_ || !mAmpLoFrac_ || !mAmpHiFrac_) return;
    double E    = mEnergy_ ? mEnergy_->GetNumber() : 500.0;
    double fwhm = res_.FWHM(E);
    if (fwhm <= 0.0) { AppendLog("Load FWHM model first (FWHM tab)."); return; }
    // Energy window = 3 × FWHM — wide enough to catch drift, tight enough to avoid neighbours
    mEnergyWin_->SetNumber(3.0 * fwhm);
    // Reset sigma fracs to defaults
    mSigLoFrac_->SetNumber(0.2);
    mSigHiFrac_->SetNumber(4.0);
    // Reset amp fracs to defaults
    mAmpLoFrac_->SetNumber(0.01);
    mAmpHiFrac_->SetNumber(20.0);
    AppendLog(Form("Bounds seeded at E=%.1f keV:  E-win=+/-%.3f keV  sig=[0.2,4.0]x model  amp=[0.01,20]x seed",
                   E, 3.0 * fwhm));
}

void GammaFitGUI::OnShowFitParams()
{
    if (fitParamDlg_) {
        if (fitParamDlg_->IsMapped()) { fitParamDlg_->RaiseWindow(); return; }
        fitParamDlg_->MapWindow();
        return;
    }

    fitParamDlg_ = new TGMainFrame(gClient->GetRoot(), 360, 10, kVerticalFrame);
    fitParamDlg_->SetWindowName("Fit Parameters");
    fitParamDlg_->SetCleanup(kNoCleanup);
    fitParamDlg_->Connect("CloseWindow()", "GammaFitGUI", this, "OnFitParamClose()");

    // Helper: add a row with label + value entry + optional bounds entries
    auto addRow = [&](TGCompositeFrame* grp, const char* lbl,
                      TGNumberEntry*& val, double v, double lo, double hi)
    {
        TGHorizontalFrame* row = new TGHorizontalFrame(grp);
        grp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 0));
        TGLabel* l = new TGLabel(row, lbl);
        l->SetWidth(90);
        row->AddFrame(l, new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        val = new TGNumberEntry(row, v, 9, -1,
            TGNumberFormat::kNESRealThree,
            TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMinMax, lo, hi);
        row->AddFrame(val, new TGLayoutHints(kLHintsExpandX));
    };

    auto addFrac = [&](TGCompositeFrame* grp, const char* lbl,
                       TGNumberEntry*& lo, double loV,
                       TGNumberEntry*& hi, double hiV,
                       const char* suffix)
    {
        TGHorizontalFrame* row = new TGHorizontalFrame(grp);
        grp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 0));
        TGLabel* l = new TGLabel(row, lbl);
        l->SetWidth(90);
        row->AddFrame(l, new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        lo = new TGNumberEntry(row, loV, 6, -1,
            TGNumberFormat::kNESRealThree, TGNumberFormat::kNEAPositive,
            TGNumberFormat::kNELLimitMinMax, 0.001, 100.0);
        row->AddFrame(lo, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
        row->AddFrame(new TGLabel(row, "to"),
                      new TGLayoutHints(kLHintsCenterY, 2, 2, 0, 0));
        hi = new TGNumberEntry(row, hiV, 6, -1,
            TGNumberFormat::kNESRealThree, TGNumberFormat::kNEAPositive,
            TGNumberFormat::kNELLimitMinMax, 0.001, 100.0);
        row->AddFrame(hi, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 0));
        row->AddFrame(new TGLabel(row, suffix),
                      new TGLayoutHints(kLHintsCenterY, 2, 0, 0, 0));
    };

    // ── Seed values ──────────────────────────────────────────────────────────
    TGGroupFrame* seedGrp = new TGGroupFrame(fitParamDlg_, "Seed Values");
    fitParamDlg_->AddFrame(seedGrp, new TGLayoutHints(kLHintsExpandX, 6, 6, 6, 4));
    addRow(seedGrp, "Last E (keV)", mEnergy_, 500.0,  0.0,  10000.0);
    addRow(seedGrp, "Sigma  (keV)", mSigma_,    1.5,  0.05,    50.0);
    addRow(seedGrp, "Amplitude",    mAmp_,    500.0,  0.0,    1.0e9);
    addRow(seedGrp, "BG const",     mBg0_,      0.0, -1.0e6,  1.0e6);
    addRow(seedGrp, "BG slope",     mBg1_,      0.0, -1.0e4,  1.0e4);
    addRow(seedGrp, "Range (xSig)", mRange_,    4.0,  1.0,     20.0);

    TGTextButton* seedBtn = new TGTextButton(seedGrp, "Seed from Resolution Model");
    seedGrp->AddFrame(seedBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 6, 4));
    seedBtn->Connect("Clicked()", "GammaFitGUI", this, "OnSeedParams()");
    seedBtn->SetToolTipText("Fill Sigma from resolution model at the last clicked energy.\nFill Amplitude from the histogram bin height.");

    // ── Fit bounds ────────────────────────────────────────────────────────────
    TGGroupFrame* bndGrp = new TGGroupFrame(fitParamDlg_, "Parameter Bounds");
    fitParamDlg_->AddFrame(bndGrp, new TGLayoutHints(kLHintsExpandX, 6, 6, 4, 4));

    // Energy: symmetric window
    {
        TGHorizontalFrame* row = new TGHorizontalFrame(bndGrp);
        bndGrp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 0));
        TGLabel* l = new TGLabel(row, "E bounds");
        l->SetWidth(90);
        row->AddFrame(l, new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        row->AddFrame(new TGLabel(row, "+/-"),
                      new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        mEnergyWin_ = new TGNumberEntry(row, 8.0, 6, -1,
            TGNumberFormat::kNESRealOne, TGNumberFormat::kNEAPositive,
            TGNumberFormat::kNELLimitMinMax, 0.1, 500.0);
        row->AddFrame(mEnergyWin_, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
        row->AddFrame(new TGLabel(row, "keV"), new TGLayoutHints(kLHintsCenterY));
    }

    addFrac(bndGrp, "Sigma bounds", mSigLoFrac_, 0.2, mSigHiFrac_, 4.0, "x model");
    addFrac(bndGrp, "Amp bounds",   mAmpLoFrac_, 0.01, mAmpHiFrac_, 20.0, "x seed");

    TGTextButton* seedBndBtn = new TGTextButton(bndGrp, "Seed Bounds from Model");
    bndGrp->AddFrame(seedBndBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 6, 4));
    seedBndBtn->Connect("Clicked()", "GammaFitGUI", this, "OnSeedBoundsFromModel()");
    seedBndBtn->SetToolTipText(
        "Set energy window = 3 x FWHM at current energy.\n"
        "Reset sigma and amplitude fractions to defaults.");

    // ── Close button ──────────────────────────────────────────────────────────
    TGTextButton* closeBtn = new TGTextButton(fitParamDlg_, "  Close  ");
    fitParamDlg_->AddFrame(closeBtn, new TGLayoutHints(kLHintsCenterX, 6, 6, 4, 8));
    closeBtn->Connect("Clicked()", "GammaFitGUI", this, "OnFitParamClose()");

    fitParamDlg_->MapSubwindows();
    fitParamDlg_->Resize(fitParamDlg_->GetDefaultSize());
    // Don't map here — caller decides whether to show
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

    // Reconstruct TF1 for this key and restore parameter errors from cache
    TF1* f = BuildFromCacheKey(key);
    if (f) {
        f->SetLineColor(kRed);
        f->SetLineWidth(2);
        const auto& entries = fitdb.GetEntries();
        auto eit = entries.find(key);
        if (eit != entries.end()) {
            const auto& errs = eit->second.paramErrors;
            if ((int)errs.size() == f->GetNpar())
                for (int i = 0; i < f->GetNpar(); i++)
                    f->SetParError(i, errs[i]);
        }
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
        if (f) {
            f->Draw("same");
            DrawPeakLabels(f);
            // Draw background component (pad-owned so it survives canvas updates)
            int npar = f->GetNpar();
            if (npar >= 5 && (npar - 2) % 3 == 0) {
                int np = (npar - 2) / 3;
                double bg0 = f->GetParameter(3*np);
                double bg1 = f->GetParameter(3*np + 1);
                TF1* bgf = new TF1(Form("nav_bg_%p", (void*)f), "[0]+[1]*x",
                                   f->GetXmin(), f->GetXmax());
                bgf->SetParameter(0, bg0);
                bgf->SetParameter(1, bg1);
                bgf->SetLineColor(kGreen + 2);
                bgf->SetLineStyle(2);
                bgf->SetLineWidth(2);
                bgf->SetBit(kCanDelete);
                bgf->Draw("same");
            }
        }
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

void GammaFitGUI::OnApplyPickedLabel()
{
    if (labelPickKey_.empty()) {
        AppendLog("No peak picked — enable Label Pick Mode and click a peak first.");
        return;
    }
    if (!mPeakLabelCombo_) return;

    TGLBEntry* le = mPeakLabelCombo_->GetSelectedEntry();
    std::string newLabel = le ? le->GetTitle() : "";
    if (newLabel == "(none)") newLabel = "";

    std::string cacheFile = CacheFileFor(currentHist_);
    FitDatabase fitdb;
    fitdb.Load(cacheFile);
    if (!fitdb.GetEntries().count(labelPickKey_)) {
        AppendLog("Picked entry no longer in cache — refresh and re-pick.");
        return;
    }
    FitEntry e = fitdb.GetEntries().at(labelPickKey_);  // mutable copy

    if (labelPickGaussIdx_ >= 0) {
        int n = ((int)e.params.size() - 2) / 3;
        if ((int)e.peakLabels.size() < n) e.peakLabels.resize(n);
        e.peakLabels[labelPickGaussIdx_] = newLabel;
        if (!newLabel.empty() && labelClassMap_.count(newLabel)) {
            if ((int)e.peakClassifications.size() < n) e.peakClassifications.resize(n);
            e.peakClassifications[labelPickGaussIdx_] = labelClassMap_.at(newLabel);
        }
    } else {
        e.label = newLabel;
        if (!newLabel.empty() && labelClassMap_.count(newLabel))
            e.classification = labelClassMap_.at(newLabel);
    }

    fitdb.ForceStore(labelPickKey_, e);
    EnsureCacheDir();
    fitdb.Save(cacheFile);

    // Refresh label info and redraw
    std::string info = Form("%.1f keV", e.params[3*(labelPickGaussIdx_>=0?labelPickGaussIdx_:0)+1]);
    if (!newLabel.empty()) info += " [" + newLabel + "]";
    if (labelPickInfo_) labelPickInfo_->SetText(info.c_str());

    OnHistViewChanged(histViewCombo_ ? histViewCombo_->GetSelected() : 1);
    if (schematicDrawn_) DrawDecaySchematic(canvas_->GetCanvas());
    AppendLog("Label applied: " + (newLabel.empty() ? "(cleared)" : newLabel) +
              "  →  peak at " + info);
}

void GammaFitGUI::OnNavXRangeGo()
{
    if (!rawHist_ || !navXMinEntry_ || !navXMaxEntry_) return;
    double lo = navXMinEntry_->GetNumber();
    double hi = navXMaxEntry_->GetNumber();
    if (lo == 0.0 && hi == 0.0) {
        viewXmin_ = 0.0;
        viewXmax_ = 0.0;
        AppendLog("X range reset to full spectrum.");
    } else {
        if (lo >= hi) { AppendLog("X range: low must be less than high."); return; }
        viewXmin_ = std::max(lo, rawHist_->GetXaxis()->GetXmin());
        viewXmax_ = std::min(hi, rawHist_->GetXaxis()->GetXmax());
        AppendLog(Form("X range set: [%.1f, %.1f] keV", viewXmin_, viewXmax_));
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

void GammaFitGUI::OnArchiveHistCache()
{
    if (currentHist_.empty()) { AppendLog("No histogram selected."); return; }

    std::string cacheFile = CacheFileFor(currentHist_);
    struct stat st;
    if (stat(cacheFile.c_str(), &st) != 0) {
        AppendLog("[Archive] No cache file for " + currentHist_);
        return;
    }

    EnsureCacheDir();
    mkdir(ArchiveDirFor().c_str(), 0755);

    // Timestamp: YYYYMMDD_HHMMSS
    time_t now = time(nullptr);
    struct tm* tm_info = localtime(&now);
    char ts[20];
    strftime(ts, sizeof(ts), "%Y%m%d_%H%M%S", tm_info);

    std::string destFile = ArchiveDirFor() + "/fit_cache_" + currentHist_ + "_" + ts + ".dat";

    if (gSystem->CopyFile(cacheFile.c_str(), destFile.c_str(), kTRUE) != 0) {
        AppendLog("[Archive] ERROR: could not write " + destFile);
        return;
    }

    AppendLog("[Archive] Saved cache to " + destFile);
    SetStatus("Archived: " + currentHist_ + " (" + ts + ")");
}

void GammaFitGUI::OnRestoreArchivedCache()
{
    if (currentHist_.empty()) { AppendLog("No histogram selected."); return; }

    // Scan archive directory for files matching this histogram
    std::vector<std::string> archives;
    DIR* dir = opendir(ArchiveDirFor().c_str());
    if (dir) {
        struct dirent* ent;
        const std::string prefix = "fit_cache_" + currentHist_ + "_";
        const std::string suffix = ".dat";
        while ((ent = readdir(dir)) != nullptr) {
            std::string fname = ent->d_name;
            if (fname.size() <= prefix.size() + suffix.size()) continue;
            if (fname.substr(0, prefix.size()) != prefix) continue;
            if (fname.substr(fname.size() - suffix.size()) != suffix) continue;
            archives.push_back(fname);
        }
        closedir(dir);
    }
    std::sort(archives.rbegin(), archives.rend()); // newest first

    if (archives.empty()) {
        AppendLog("[Restore] No archived caches found for " + currentHist_);
        return;
    }

    // Build picker dialog
    TGTransientFrame* dlg = new TGTransientFrame(gClient->GetRoot(), this, 460, 10, kVerticalFrame);
    dlg->SetWindowName("Restore Archived Cache");
    dlg->SetCleanup(kDeepCleanup);

    dlg->AddFrame(new TGLabel(dlg, Form("Histogram: %s", currentHist_.c_str())),
                  new TGLayoutHints(kLHintsLeft, 8, 8, 8, 2));
    dlg->AddFrame(new TGLabel(dlg, "Select archive to restore (newest first):"),
                  new TGLayoutHints(kLHintsLeft, 8, 8, 0, 4));

    TGComboBox* archCombo = new TGComboBox(dlg, 600);
    archCombo->Resize(430, 22);
    dlg->AddFrame(archCombo, new TGLayoutHints(kLHintsExpandX, 8, 8, 0, 4));
    for (int i = 0; i < (int)archives.size(); i++) {
        // Strip prefix and suffix to show just the timestamp
        const std::string prefix = "fit_cache_" + currentHist_ + "_";
        std::string ts = archives[i].substr(prefix.size(),
                                            archives[i].size() - prefix.size() - 4);
        archCombo->AddEntry(ts.c_str(), i + 1);
    }
    archCombo->Select(1, kFALSE);

    dlg->AddFrame(new TGLabel(dlg, "This will overwrite the current live cache."),
                  new TGLayoutHints(kLHintsLeft, 8, 8, 0, 4));

    TGHorizontalFrame* btnRow = new TGHorizontalFrame(dlg);
    dlg->AddFrame(btnRow, new TGLayoutHints(kLHintsCenterX, 8, 8, 4, 8));
    TGTextButton* okBtn  = new TGTextButton(btnRow, "  Restore  ");
    TGTextButton* canBtn = new TGTextButton(btnRow, "  Cancel  ");
    btnRow->AddFrame(okBtn,  new TGLayoutHints(kLHintsLeft, 0, 8, 0, 0));
    btnRow->AddFrame(canBtn, new TGLayoutHints(kLHintsLeft));
    okBtn->Connect("Clicked()", "TObject", dlg, "SetBit(UInt_t=16384)");
    okBtn->Connect("Clicked()", "TGFrame", dlg, "UnmapWindow()");
    canBtn->Connect("Clicked()", "TGFrame", dlg, "UnmapWindow()");

    dlg->MapSubwindows();
    dlg->Resize(dlg->GetDefaultSize());
    dlg->CenterOnParent();
    dlg->MapWindow();

    gClient->WaitForUnmap(dlg);

    Int_t selId    = archCombo->GetSelected();
    bool confirmed = dlg->TestBit(BIT(14));
    dlg->DeleteWindow();

    if (!confirmed || selId < 1 || (size_t)selId > archives.size()) return;

    std::string srcFile  = ArchiveDirFor() + "/" + archives[selId - 1];
    std::string destFile = CacheFileFor(currentHist_);

    EnsureCacheDir();
    if (gSystem->CopyFile(srcFile.c_str(), destFile.c_str(), kTRUE) != 0) {
        AppendLog("[Restore] ERROR: could not restore " + srcFile);
        return;
    }

    AppendLog("[Restore] Restored " + archives[selId - 1] + " → " + destFile);

    // Reload the restored cache onto the canvas
    OnLoadCache();
    SetStatus("Cache restored: " + currentHist_);
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
    if (id >= 1 && id <= (Int_t)peakNavKeys_.size())
        NavigateToPeak(id - 1);
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
    double xlo = f->GetXmin();
    double xhi = f->GetXmax();
    double bg0, bg1;
    int n = 0;
    bool isDG = (npar == 7);

    if (isDG) {
        // Double-Gaussian: [0] A_narrow [1] mean [2] sig_narrow [3] A_broad [4] sig_broad [5] bg0 [6] bg1
        bg0 = f->GetParameter(5);
        bg1 = f->GetParameter(6);
    } else if (npar >= 5 && (npar - 2) % 3 == 0) {
        n   = (npar - 2) / 3;
        bg0 = f->GetParameter(3*n);
        bg1 = f->GetParameter(3*n + 1);
    } else {
        return;  // unrecognised model
    }

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

    // Individual Gaussians + BG: blue dashed (only when n > 1 and option is on; not for DG model)
    for (int i = 0; !isDG && i < n && n > 1 && (mShowCompChk_ && mShowCompChk_->IsDown()); i++) {
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
    if (!f || !gPad) return;
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

    std::vector<PkLblData> pkLabels;
    std::vector<double> drawn;
    for (int i = 0; i < (int)peakEs.size(); i++) {
        double E = peakEs[i];
        double sep = res_.FWHM(E);
        bool tooClose = false;
        for (double de : drawn)
            if (std::abs(E - de) < sep) { tooClose = true; break; }
        if (tooClose) continue;
        drawn.push_back(E);

        // Centroid uncertainty from fit parameter
        int parIdx = isDG ? 1 : 3*i + 1;
        double eErr = f->GetParError(parIdx);

        double yTop = f->Eval(E);
        std::string isoName;
        if (showIsoLabels_ && dbLoaded_) {
            auto matches = db_.Match(E, res_.FWHM(E));
            if (!matches.empty()) isoName = matches[0].isotope;
        }
        pkLabels.push_back({E, yTop, isoName, NNDCFormat(E, eErr)});
    }
    DrawStaggedLabels(pkLabels, gPad);
}

void GammaFitGUI::UpdatePeakStats(TF1* f, TH1* h, double xlo, double xhi)
{
    if (!peakStatsView_) return;

    auto buildLines = [&](TF1* tf, TH1* th, double lo, double hi) -> std::vector<std::string> {
        std::vector<std::string> lines;
        if (!tf || !th) return lines;
        int npar = tf->GetNpar();
        // Detect model layout: standard (3n+2), quadBG (3n+3), comptonStep (4n+2), both (4n+3)
        int n = -1; bool hasQuadBG = false;
        for (int nc = 1; nc <= 12; nc++) {
            if      (npar == 3*nc + 2) { n = nc; break; }
            else if (npar == 3*nc + 3) { n = nc; hasQuadBG = true; break; }
            else if (npar == 4*nc + 2) { n = nc; break; }               // compton step only
            else if (npar == 4*nc + 3) { n = nc; hasQuadBG = true; break; } // compton + quad
        }
        if (n < 1) { lines.push_back("Invalid fit model"); return lines; }
        int bgBase = 3*n;

        // Chi2/ndf and p-value from histogram residuals
        double chi2 = 0.0; int cnt = 0;
        int b1 = th->FindBin(lo), b2 = th->FindBin(hi);
        for (int b = b1; b <= b2; b++) {
            double err = th->GetBinError(b);
            if (err <= 0.0) continue;
            double pull = (th->GetBinContent(b) - tf->Eval(th->GetBinCenter(b))) / err;
            chi2 += pull * pull; ++cnt;
        }
        int    ndf     = cnt - npar;
        double chi2ndf = (ndf > 0) ? chi2 / ndf : (cnt > 0 ? chi2 / cnt : -1.0);
        double pval    = (ndf > 0) ? TMath::Prob(chi2, ndf) : -1.0;

        // Total counts in fit window (for local P/T)
        double windowCounts = 0.0;
        for (int b = b1; b <= b2; b++) windowCounts += th->GetBinContent(b);

        // Total counts in the full raw spectrum (for global P/T)
        double spectrumCounts = 0.0;
        if (rawHist_) {
            for (int b = 1; b <= rawHist_->GetNbinsX(); b++)
                spectrumCounts += rawHist_->GetBinContent(b);
        }

        double bg0    = tf->GetParameter(bgBase);
        double bg1    = tf->GetParameter(bgBase+1);
        double bg0err = tf->GetParError(bgBase);
        double bg1err = tf->GetParError(bgBase+1);
        double bg2    = hasQuadBG ? tf->GetParameter(bgBase+2) : 0.0;
        double bg2err = hasQuadBG ? tf->GetParError(bgBase+2)  : 0.0;

        // ── Area ratio: fitted Gaussian integral vs observed counts above BG ──
        {
            double obsAboveBg = 0.0;
            for (int b = b1; b <= b2; b++) {
                double x = th->GetBinCenter(b);
                obsAboveBg += th->GetBinContent(b) - (bg0 + bg1*x + bg2*x*x);
            }

            double totalFitArea = 0.0;
            for (int i = 0; i < n; i++) {
                double A   = tf->GetParameter(3*i);
                double E   = tf->GetParameter(3*i+1);
                double sig = tf->GetParameter(3*i+2);
                double bw  = th->GetBinWidth(th->FindBin(E));
                if (bw <= 0.0) bw = 1.0;
                totalFitArea += A * sig * std::sqrt(2.0 * TMath::Pi()) / bw;
            }

            if (obsAboveBg > 1.0) {
                double ratio = totalFitArea / obsAboveBg;
                const char* grade =
                    (ratio >= 0.90 && ratio <= 1.10) ? "GOOD" :
                    (ratio >= 0.70 && ratio <= 1.30) ? "fair" : "*** POOR ***";
                lines.push_back(Form("  Area ratio (fit/obs) = %.3f  %s", ratio, grade));
                lines.push_back(Form("    fit area = %.0f   obs above BG = %.0f",
                                     totalFitArea, obsAboveBg));
            } else {
                lines.push_back("  Area ratio (fit/obs) = n/a  (obs counts too low)");
            }
            lines.push_back("  ----------------------------------------");
        }

        for (int i = 0; i < n; i++) {
            double A    = tf->GetParameter(3*i);
            double E    = tf->GetParameter(3*i+1);
            double sig  = tf->GetParameter(3*i+2);
            double Aerr = tf->GetParError(3*i);
            double Eerr = tf->GetParError(3*i+1);
            double serr = tf->GetParError(3*i+2);

            double fwhm    = 2.3548 * sig;
            double fwhmerr = 2.3548 * serr;
            // bin-width correction: ROOT fits counts/bin so area needs /binWidth
            double binWidth = th->GetBinWidth(th->FindBin(E));
            if (binWidth <= 0.0) binWidth = 1.0;
            double peakArea = A * sig * std::sqrt(2.0 * TMath::Pi()) / binWidth;
            double areaErr  = (A > 0 && sig > 0)
                ? peakArea * std::sqrt(std::pow(Aerr/A, 2) + std::pow(serr/sig, 2))
                : 0.0;
            // Background counts under ±2.5 sigma of the peak (uses full model)
            double bgAtE  = bg0 + bg1*E + bg2*E*E;
            double bgUnder = std::abs(bgAtE / binWidth * 5.0 * sig);
            // SNR = peak area / sqrt(background counts in peak region)
            double snr      = (bgUnder > 0) ? peakArea / std::sqrt(bgUnder) : 0.0;
            double ptLocal  = (windowCounts  > 0) ? peakArea / windowCounts  : 0.0;
            double ptGlobal = (spectrumCounts > 0) ? peakArea / spectrumCounts : 0.0;

            if (n > 1)
                lines.push_back(Form("--- Peak %d (%.2f keV) ---", i+1, E));
            lines.push_back(Form("  Amplitude     = %.3g +/- %.3g", A,    Aerr));
            lines.push_back(Form("  E (keV)       = %.4f +/- %.4f", E,    Eerr));
            lines.push_back(Form("  sigma (keV)   = %.4f +/- %.4f", sig,  serr));
            lines.push_back(Form("  FWHM (keV)    = %.4f +/- %.4f", fwhm, fwhmerr));
            lines.push_back(Form("  Peak area     = %.1f +/- %.1f counts", peakArea, areaErr));
            if (spectrumCounts > 0 && peakArea > spectrumCounts)
                lines.push_back("  *** WARNING: area > spectrum total — fit likely bad ***");
            else if (peakArea < 0)
                lines.push_back("  *** WARNING: negative area — check amplitude/sigma ***");
            lines.push_back(Form("  SNR (Np/sqrtNbg) = %.2f",        snr));
            lines.push_back(Form("  P/T (window)  = %.4f  [%.0f / %.0f]", ptLocal,  peakArea, windowCounts));
            lines.push_back(Form("  P/T (all)     = %.4f  [%.0f / %.0f]", ptGlobal, peakArea, spectrumCounts));

            // Efficiency-corrected area (if parameters are loaded)
            if (effB_ || effC_ || effD_) {
                double ea = effA_ ? effA_->GetNumber() : 0.0;
                double eb = effB_ ? effB_->GetNumber() : 0.0;
                double ec = effC_ ? effC_->GetNumber() : 0.0;
                double ed = effD_ ? effD_->GetNumber() : 0.0;
                if (eb != 0.0 || ec != 0.0 || ed != 0.0) {
                    double lnE   = (E > 0) ? std::log(E) : 0.0;
                    double lnEps = ea - eb*lnE + ec*lnE*lnE - (E > 0 ? ed/(E*E) : 0.0);
                    double eps   = std::exp(lnEps);
                    if (eps > 0 && std::isfinite(eps)) {
                        double corrArea = peakArea / eps;
                        double corrErr  = areaErr  / eps;
                        lines.push_back(Form("  Eff-corr area = %.1f +/- %.1f  (eps=%.4e)",
                                             corrArea, corrErr, eps));
                    }
                }
            }
        }

        if (hasQuadBG) {
            lines.push_back(Form("  bg0           = %.4g +/- %.4g", bg0, bg0err));
            lines.push_back(Form("  bg1           = %.5g +/- %.5g", bg1, bg1err));
            lines.push_back(Form("  bg2           = %.5g +/- %.5g  (quadratic)", bg2, bg2err));
        } else if (bg1 != 0.0) {
            lines.push_back(Form("  bg0           = %.4g +/- %.4g", bg0, bg0err));
            lines.push_back(Form("  bg1           = %.5g +/- %.5g", bg1, bg1err));
        } else {
            lines.push_back(Form("  bg0           = %.4g +/- %.4g  (flat)", bg0, bg0err));
        }

        if (ndf > 0) {
            lines.push_back(Form("  chi2/ndf      = %.3f  (%d DOF)", chi2ndf, ndf));
            lines.push_back(Form("  p-value      = %s",
                pval < 1e-4 ? Form("%.3e", pval) : Form("%.4f", pval)));
        }
        return lines;
    };

    peakStatsCurrent_ = buildLines(f, h, xlo, xhi);

    peakStatsView_->Clear();

    // If we have old stats (from before RunFit), show them above a separator
    if (!peakStatsOld_.empty()) {
        for (const auto& line : peakStatsOld_)
            peakStatsView_->AddLine(line.c_str());
        peakStatsView_->AddLine("  ----------------------------------------");
        peakStatsOld_.clear();
    }

    for (const auto& line : peakStatsCurrent_)
        peakStatsView_->AddLine(line.c_str());

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

        // Downsample for display: COLZ draws every bin as a rectangle.
        // A 2000×2000 TH2 = 4M rectangles, which freezes the GUI.
        // Rebin to at most kMaxAxis bins per dimension before painting.
        TH2* h2 = (TH2*)h;
        int nX = h2->GetNbinsX(), nY = h2->GetNbinsY();
        const int kMaxAxis = 500;
        int rx = std::max(1, nX / kMaxAxis);
        int ry = std::max(1, nY / kMaxAxis);

        if (rx > 1 || ry > 1) {
            // Downsample for display speed — original data is untouched.
            TH2* disp = (TH2*)h2->Clone("_th2_disp_tmp");
            disp->SetDirectory(nullptr);
            disp->Rebin2D(rx, ry);
            disp->SetBit(kCanDelete);
            // Title makes the downsampling impossible to miss
            disp->SetTitle(Form("%s  [DISPLAY ONLY: %dx%d rebin — original %dx%d bins]",
                                h->GetName(), rx, ry, nX, nY));
            disp->Draw("COLZ");
            SetStatus(Form("WARNING: %s displayed at %dx%d rebin (original %d x %d bins)",
                           h->GetName(), rx, ry, nX, nY));
        } else {
            h->Draw("COLZ");
        }

        c->Modified(); c->Update();
        return;
    }

    ApplyHistStyle(h, h->GetName());
    if (!currentHist_.empty() && ClassOf(currentHist_) == "Decay")
        h->GetXaxis()->SetTitle("Time (ms)");
    h->SetLineColor(kBlack);
    h->SetMarkerSize(0);
    h->SetStats(0);
    SetYMaxFromVisible(h);
    h->Draw("hist");
    h->Draw("E1 same");
    if (fit) {
        fit->SetLineColor(kRed);
        fit->SetLineWidth(2);
        fit->Draw("same");
    }

    // Integral statistics overlay (top-right corner)
    {
        TPaveText* pt = new TPaveText(0.55, 0.75, 0.98, 0.98, "NDC");
        pt->SetBorderSize(1);
        pt->SetFillColor(kWhite);
        pt->SetFillStyle(1001);
        pt->SetTextAlign(12);
        pt->SetTextSize(0.030);
        pt->SetBit(kCanDelete);
        pt->AddText(Form("Integral = %.0f", h->Integral()));
        pt->AddText(Form("Entries  = %.0f", h->GetEntries()));
        pt->AddText(Form("Mean     = %.2f keV", h->GetMean()));
        pt->AddText(Form("RMS      = %.2f keV", h->GetRMS()));
        auto rbit = rebinFactors_.find(currentHist_);
        if (rbit != rebinFactors_.end() && rbit->second > 1)
            pt->AddText(Form("Rebin x%d", rbit->second));
        pt->Draw();
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

    // Discard any bg-subtracted view and redraw the raw histogram that is
    // already stored in rawHist_.  Loading a fresh copy and then deleting it
    // caused a dangling pointer inside the canvas, producing a blank display.
    delete viewHist_; viewHist_ = nullptr;

    if (rawHist_) {
        RedrawView();
        AppendLog("BG sub reset — showing raw histogram: " + currentHist_);
        SetStatus("Raw: " + currentHist_);
    } else {
        AppendLog("Background subtraction reset to defaults (14 iterations, enabled).");
    }
}

void GammaFitGUI::OnDeleteHistogram()
{
    Int_t id = histList_->GetSelected();
    if (id < 1 || (size_t)id > histNames_.size()) {
        AppendLog("[Delete] Select a histogram first."); return;
    }
    const std::string name = histNames_[id - 1];

    // Virtual histograms (projections, bg-subtracted) live only in the session.
    bool isVirtual = customProjDefs_.count(name) || projParent_.count(name)
                     || bgSubtractDefs_.count(name);

    // Confirm — this is destructive for file-backed histograms.
    Int_t ret = kMBNo;
    std::string msg = isVirtual
        ? "Remove '" + name + "' from this session?\n(It is not stored in the file.)"
        : "Permanently delete '" + name + "' from the ROOT file?\n\nThis cannot be undone.";
    new TGMsgBox(gClient->GetRoot(), this,
                 isVirtual ? "Remove Virtual Histogram" : "Delete from File",
                 msg.c_str(),
                 isVirtual ? kMBIconQuestion : kMBIconExclamation,
                 kMBYes | kMBNo, &ret);
    if (ret != kMBYes) return;

    // ── Delete from file (real histograms only) ──────────────────────────────
    if (!isVirtual) {
        if (!inputFile_) { AppendLog("[Delete] No file open."); return; }
        std::string path = inputFile_->GetName();

        // Close READ handle, reopen for UPDATE
        inputFile_->Close(); delete inputFile_; inputFile_ = nullptr;
        TFile* uf = TFile::Open(path.c_str(), "UPDATE");
        if (!uf || uf->IsZombie()) {
            AppendLog("[Delete] Cannot open file for writing: " + path);
            // Reopen read-only so the session remains usable
            inputFile_ = TFile::Open(path.c_str(), "READ");
            return;
        }

        // Delete all cycles of the key
        uf->Delete((name + ";*").c_str());
        uf->Purge();   // reclaim disk space
        uf->Close(); delete uf;

        // Reopen in READ mode
        inputFile_ = TFile::Open(path.c_str(), "READ");
        if (!inputFile_ || inputFile_->IsZombie()) {
            AppendLog("[Delete] Re-open failed after deletion: " + path);
            inputFile_ = nullptr;
        }

        // Also remove auto-generated projection siblings (name_px / name_py)
        // if this was a TH2
        if (th2Names_.count(name)) {
            for (const std::string& sibling : {name + "_px", name + "_py"}) {
                histNames_.erase(std::remove(histNames_.begin(), histNames_.end(), sibling), histNames_.end());
                projParent_.erase(sibling);
                histClass_.erase(sibling);
            }
            th2Names_.erase(name);
        }

        AppendLog("[Delete] Removed '" + name + "' from " + path);
    } else {
        // Remove derived definitions
        customProjDefs_.erase(name);
        projParent_.erase(name);
        bgSubtractDefs_.erase(name);
        AppendLog("[Delete] Removed virtual histogram '" + name + "' from session.");
    }

    // ── Clean up session state ───────────────────────────────────────────────
    histNames_.erase(std::remove(histNames_.begin(), histNames_.end(), name), histNames_.end());
    histClass_.erase(name);
    rebinFactors_.erase(name);

    // Also remove any virtual histograms that depended on this one
    std::vector<std::string> toRemove;
    for (const auto& kv : bgSubtractDefs_)
        if (kv.second.srcName == name || kv.second.bgName == name)
            toRemove.push_back(kv.first);
    for (const auto& kv : customProjDefs_)
        if (kv.second.th2Name == name)
            toRemove.push_back(kv.first);
    for (const auto& dep : toRemove) {
        histNames_.erase(std::remove(histNames_.begin(), histNames_.end(), dep), histNames_.end());
        bgSubtractDefs_.erase(dep);
        customProjDefs_.erase(dep);
        histClass_.erase(dep);
        rebinFactors_.erase(dep);
        AppendLog("[Delete] Also removed dependent histogram '" + dep + "'.");
    }

    // Clear canvas if the deleted histogram was on display
    if (currentHist_ == name) {
        if (rawHistOwned_) { delete rawHist_; rawHistOwned_ = false; }
        rawHist_ = nullptr;
        delete viewHist_; viewHist_ = nullptr;
        currentHist_.clear();
        canvas_->GetCanvas()->Clear();
        canvas_->GetCanvas()->Modified();
        canvas_->GetCanvas()->Update();
    }

    SaveMetadata();
    PopulateHistWidgets();
    PopulateBgSubCombos();
    PopulateCustProjTh2Combo();
    SetStatus("Deleted: " + name);
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
    return CacheDirFor() + "/fit_cache_" + hname + "_source.root";
}

void GammaFitGUI::SaveSourceAnalysis()
{
    if (srcHist_.empty()) return;
    EnsureCacheDir();
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
              " keV -> fitted peak " + Fmt(srcPeakEs_[bestIdx]) + " keV");
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
    if (bgSubtractDefs_.count(name)) {
        // Background-subtracted histograms inherit their source's class unless
        // the user has explicitly overridden it
        auto it = histClass_.find(name);
        if (it != histClass_.end()) return it->second;
        // Default to the source histogram's class
        const std::string& src = bgSubtractDefs_.at(name).srcName;
        auto sit = histClass_.find(src);
        return (sit != histClass_.end()) ? sit->second : "Gamma";
    }
    auto it = histClass_.find(name);
    return (it != histClass_.end()) ? it->second : "Gamma";
}

std::string GammaFitGUI::MetadataFileFor() const
{
    if (inputPath_.empty()) return "";
    return CacheDirFor() + "/metadata.txt";
}

void GammaFitGUI::SaveMetadata() const
{
    std::string path = MetadataFileFor();
    if (path.empty()) return;
    EnsureCacheDir();
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
    out << "# background subtractions\n";
    for (const auto& kv : bgSubtractDefs_) {
        const auto& d = kv.second;
        out << "bgsub\t" << kv.first << "\t" << d.srcName
            << "\t" << d.bgName
            << "\t" << std::fixed << std::setprecision(6) << d.scale << "\n";
    }
    out << "# rebin factors\n";
    for (const auto& kv : rebinFactors_)
        if (kv.second > 1)
            out << "rebin\t" << kv.first << "\t" << kv.second << "\n";
}

void GammaFitGUI::LoadMetadata()
{
    std::string path = MetadataFileFor();
    if (path.empty()) return;
    std::ifstream in(path);
    if (!in.is_open()) {
        AppendLog("[Metadata] No metadata file found — no virtual histograms restored.");
        AppendLog("[Metadata]   (looked for: " + path + ")");
        return;
    }
    int nClasses = 0, nProj = 0, nBgSub = 0, nRebin = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string tag;
        ss >> tag;
        if (tag == "class") {
            std::string name, cls;
            ss >> name >> cls;
            if (!name.empty() && !cls.empty()) { histClass_[name] = cls; nClasses++; }
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
                if (std::find(histNames_.begin(), histNames_.end(), name) == histNames_.end())
                    histNames_.push_back(name);
                nProj++;
            }
        } else if (tag == "bgsub") {
            std::string name, srcName, bgName;
            double scale = 1.0;
            ss >> name >> srcName >> bgName >> scale;
            if (!name.empty() && !srcName.empty() && !bgName.empty()) {
                BgSubtractDef def;
                def.srcName = srcName;
                def.bgName  = bgName;
                def.scale   = scale;
                bgSubtractDefs_[name] = def;
                if (std::find(histNames_.begin(), histNames_.end(), name) == histNames_.end())
                    histNames_.push_back(name);
                nBgSub++;
            }
        } else if (tag == "rebin") {
            std::string name; int factor = 1;
            ss >> name >> factor;
            if (!name.empty() && factor > 1) {
                rebinFactors_[name] = factor;
                nRebin++;
            }
        }
    }
    AppendLog(Form("[Metadata] Loaded: %d types, %d custom projections, %d bg-subtracted, %d rebin",
                   nClasses, nProj, nBgSub, nRebin));
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
    std::string cls;
    if      (sel == 2) cls = "Decay";
    else if (sel == 3) cls = "2D";
    else if (sel == 4) cls = "Background";
    else               cls = "Gamma";
    histClass_[name] = cls;
    SaveMetadata();
    PopulateHistWidgets();
    PopulateBgSubCombos();
    AppendLog("[Type] " + name + " -> " + cls);
    SetStatus("Histogram type: " + name + " = " + cls);
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

// ─────────────────────────────────────────────────────────────────────────────
// Background subtraction helpers
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::PopulateBgSubCombos()
{
    if (!bgSubSrcCombo_ || !bgSubBgCombo_) return;
    bgSubSrcCombo_->RemoveAll();
    bgSubBgCombo_->RemoveAll();
    int srcIdx = 1, bgIdx = 1;
    for (const auto& name : histNames_) {
        if (th2Names_.count(name)) continue;  // skip 2D
        bgSubSrcCombo_->AddEntry(name.c_str(), srcIdx++);
        if (ClassOf(name) == "Background")
            bgSubBgCombo_->AddEntry(name.c_str(), bgIdx++);
    }
    bgSubSrcCombo_->MapSubwindows(); bgSubSrcCombo_->Layout();
    bgSubBgCombo_->MapSubwindows();  bgSubBgCombo_->Layout();
}

void GammaFitGUI::OnSubtractHistogram()
{
    if (!inputFile_) { AppendLog("[BgSub] No ROOT file loaded."); return; }

    // --- gather inputs ---
    int srcSel = bgSubSrcCombo_ ? bgSubSrcCombo_->GetSelected() : -1;
    int bgSel  = bgSubBgCombo_  ? bgSubBgCombo_->GetSelected()  : -1;
    if (srcSel < 1 || bgSel < 1) {
        AppendLog("[BgSub] Select source and background histograms."); return;
    }

    // Map combo IDs back to names (combos were filled in order, skipping 2D)
    auto nameFromCombo = [&](TGComboBox* cb, int sel) -> std::string {
        TGLBEntry* entry = cb->GetListBox()->GetEntry(sel);
        return entry ? std::string(entry->GetTitle()) : "";
    };
    std::string srcName = nameFromCombo(bgSubSrcCombo_, srcSel);
    std::string bgName  = nameFromCombo(bgSubBgCombo_,  bgSel);
    if (srcName.empty() || bgName.empty()) {
        AppendLog("[BgSub] Could not resolve histogram names."); return;
    }
    if (srcName == bgName) {
        AppendLog("[BgSub] Source and background cannot be the same histogram."); return;
    }

    double scale = bgSubScaleEntry_ ? bgSubScaleEntry_->GetNumber() : 1.0;
    if (scale < 0.0) scale = 1.0;

    std::string outName = bgSubNameEntry_ ? std::string(bgSubNameEntry_->GetText()) : "subtracted";
    if (outName.empty()) outName = srcName + "_bgsub";

    // Check for duplicate name
    if (std::find(histNames_.begin(), histNames_.end(), outName) != histNames_.end()) {
        AppendLog("[BgSub] Name '" + outName + "' already exists. Choose a different name.");
        return;
    }

    // --- load and subtract ---
    bool srcOwned = false, bgOwned = false;
    TH1* src = LoadHistFromFile(srcName, srcOwned);
    TH1* bg  = LoadHistFromFile(bgName,  bgOwned);
    if (!src) { AppendLog("[BgSub] Failed to load source: " + srcName); return; }
    if (!bg)  {
        if (srcOwned) delete src;
        AppendLog("[BgSub] Failed to load background: " + bgName); return;
    }

    // Verify compatibility before subtracting
    if (src->GetNbinsX() != bg->GetNbinsX() ||
        std::abs(src->GetXaxis()->GetXmin() - bg->GetXaxis()->GetXmin()) > 1e-6 ||
        std::abs(src->GetXaxis()->GetXmax() - bg->GetXaxis()->GetXmax()) > 1e-6) {
        AppendLog(Form("[BgSub] ERROR: Incompatible histograms — cannot subtract."));
        AppendLog(Form("  Source:     %d bins  [%.3g, %.3g]",
                       src->GetNbinsX(),
                       src->GetXaxis()->GetXmin(), src->GetXaxis()->GetXmax()));
        AppendLog(Form("  Background: %d bins  [%.3g, %.3g]",
                       bg->GetNbinsX(),
                       bg->GetXaxis()->GetXmin(), bg->GetXaxis()->GetXmax()));
        if (srcOwned) delete src;
        if (bgOwned)  delete bg;
        return;
    }

    double srcIntegral = src->Integral();
    double bgIntegral  = bg->Integral();

    TH1* result = (TH1*)src->Clone(outName.c_str());
    result->SetDirectory(nullptr);
    result->Sumw2();
    result->Add(bg, -scale);

    double resIntegral = result->Integral();
    int negBins = 0;
    for (int b = 1; b <= result->GetNbinsX(); b++)
        if (result->GetBinContent(b) < 0) negBins++;

    AppendLog(Form("[BgSub] Source integral:      %.0f", srcIntegral));
    AppendLog(Form("[BgSub] Background integral:  %.0f  (scale=%.4g → subtract %.0f)",
                   bgIntegral, scale, scale * bgIntegral));
    AppendLog(Form("[BgSub] Result integral:      %.0f  (expected %.0f)",
                   resIntegral, srcIntegral - scale * bgIntegral));
    if (negBins > 0)
        AppendLog(Form("[BgSub] Note: %d bins went negative (background > source in those bins)",
                       negBins));

    if (srcOwned) delete src;
    if (bgOwned)  delete bg;

    // --- store definition ---
    BgSubtractDef def;
    def.srcName = srcName;
    def.bgName  = bgName;
    def.scale   = scale;
    bgSubtractDefs_[outName] = def;
    histNames_.push_back(outName);

    // Inherit source classification
    auto cit = histClass_.find(srcName);
    if (cit != histClass_.end())
        histClass_[outName] = cit->second;

    SaveMetadata();
    PopulateHistWidgets();
    PopulateBgSubCombos();

    AppendLog(Form("[BgSub] Created '%s' = '%s' - %.4g * '%s'",
                   outName.c_str(), srcName.c_str(), scale, bgName.c_str()));
    SetStatus("Background subtracted: " + outName);

    // Display the result
    delete result;  // we'll reload via LoadHistFromFile for consistency
    bool owned = false;
    TH1* display = LoadHistFromFile(outName, owned);
    if (display) {
        if (rawHistOwned_) { delete rawHist_; rawHist_ = nullptr; rawHistOwned_ = false; }
        rawHist_     = display;
        rawHistOwned_ = owned;
        currentHist_  = outName;
        DrawOnCanvas(rawHist_);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Recent files
// ─────────────────────────────────────────────────────────────────────────────
// Rebin slots
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnRebinPreview()
{
    if (!rawHist_ || currentHist_.empty()) {
        AppendLog("[Rebin] Load a histogram first."); return;
    }
    if (rawHist_->InheritsFrom("TH2")) {
        AppendLog("[Rebin] Cannot rebin a 2D histogram."); return;
    }
    int n = rebinEntry_ ? (int)rebinEntry_->GetNumber() : 1;
    if (n <= 1) { DrawOnCanvas(rawHist_); return; }

    TH1* preview = (TH1*)rawHist_->Clone("_rebin_preview_");
    preview->SetDirectory(nullptr);
    preview->SetBit(kCanDelete);
    // Sumw2 ensures Poisson errors add in quadrature when merging bins
    preview->Sumw2();
    preview->Rebin(n);
    preview->SetTitle(Form("%s  [PREVIEW: rebin x%d]", rawHist_->GetName(), n));
    DrawOnCanvas(preview);
    AppendLog(Form("[Rebin] Preview rebin x%d for %s — use Apply to store.",
                   n, currentHist_.c_str()));
    SetStatus(Form("PREVIEW rebin x%d: %s", n, currentHist_.c_str()));
}

void GammaFitGUI::OnRebinApply()
{
    if (currentHist_.empty() || !rawHist_) {
        AppendLog("[Rebin] Load a histogram first."); return;
    }
    if (rawHist_->InheritsFrom("TH2")) {
        AppendLog("[Rebin] Cannot rebin a 2D histogram."); return;
    }
    int n = rebinEntry_ ? (int)rebinEntry_->GetNumber() : 1;
    if (n < 1) n = 1;

    if (n == 1) {
        rebinFactors_.erase(currentHist_);
        AppendLog("[Rebin] Rebin factor cleared for " + currentHist_);
    } else {
        rebinFactors_[currentHist_] = n;
        AppendLog(Form("[Rebin] Stored rebin x%d for %s", n, currentHist_.c_str()));
    }

    // Reload histogram with the new factor applied
    if (rawHistOwned_) { delete rawHist_; rawHist_ = nullptr; rawHistOwned_ = false; }
    rawHist_ = LoadHistFromFile(currentHist_, rawHistOwned_);
    delete viewHist_; viewHist_ = nullptr;

    SaveMetadata();
    PopulateHistWidgets();
    DrawOnCanvas(rawHist_);
    SetStatus(n > 1 ? Form("Rebin x%d applied: %s", n, currentHist_.c_str())
                    : ("Rebin cleared: " + currentHist_));
}

void GammaFitGUI::OnRebinReset()
{
    if (currentHist_.empty()) return;
    rebinFactors_.erase(currentHist_);
    if (rebinEntry_) rebinEntry_->SetNumber(1);

    if (rawHistOwned_) { delete rawHist_; rawHist_ = nullptr; rawHistOwned_ = false; }
    rawHist_ = LoadHistFromFile(currentHist_, rawHistOwned_);
    delete viewHist_; viewHist_ = nullptr;

    SaveMetadata();
    PopulateHistWidgets();
    DrawOnCanvas(rawHist_);
    AppendLog("[Rebin] Reset to original binning: " + currentHist_);
    SetStatus("Rebin reset: " + currentHist_);
}

void GammaFitGUI::OnToggleIsoLabels()
{
    showIsoLabels_ = showIsoLabelsChk_ && showIsoLabelsChk_->IsOn();
    if (rawHist_) RedrawView();
}

// ─────────────────────────────────────────────────────────────────────────────
// OnApplyBgAnchors — seed bg0, bg1 from two off-peak anchor regions
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnApplyBgAnchors()
{
    if (!rawHist_) { AppendLog("Load a histogram first."); return; }
    if (!mBgAnch1Lo_ || !mBgAnch1Hi_ || !mBgAnch2Lo_ || !mBgAnch2Hi_) return;

    double lo1 = mBgAnch1Lo_->GetNumber(), hi1 = mBgAnch1Hi_->GetNumber();
    double lo2 = mBgAnch2Lo_->GetNumber(), hi2 = mBgAnch2Hi_->GetNumber();

    if (hi1 <= lo1 || hi2 <= lo2) {
        AppendLog("Anchor regions are invalid (Hi must be > Lo for both regions).");
        return;
    }
    if (std::abs((lo1+hi1)*0.5 - (lo2+hi2)*0.5) < 1.0) {
        AppendLog("Anchor regions overlap — use two distinct off-peak regions.");
        return;
    }

    TH1* h = viewHist_ ? viewHist_ : rawHist_;

    auto regionMean = [&](double lo, double hi) -> std::pair<double,double> {
        int b1 = h->FindBin(lo + h->GetBinWidth(1)*0.01);
        int b2 = h->FindBin(hi - h->GetBinWidth(1)*0.01);
        if (b1 > b2) std::swap(b1, b2);
        double sum = 0.0; int nb = 0;
        for (int b = b1; b <= b2; b++) { sum += h->GetBinContent(b); nb++; }
        return { (lo + hi) * 0.5, (nb > 0) ? sum / nb : 0.0 };
    };

    auto [x1, y1] = regionMean(lo1, hi1);
    auto [x2, y2] = regionMean(lo2, hi2);

    double bg1val = (std::abs(x2 - x1) > 1e-6) ? (y2 - y1) / (x2 - x1) : 0.0;
    double bg0val = y1 - bg1val * x1;

    mBg0_->SetNumber(bg0val);
    mBg1_->SetNumber(bg1val);

    AppendLog(Form("BG anchors applied: bg0=%.2f  bg1=%.5f  (regions [%.0f,%.0f] and [%.0f,%.0f])",
                   bg0val, bg1val, lo1, hi1, lo2, hi2));
}

// ─────────────────────────────────────────────────────────────────────────────
// OnApplyCalibration — apply E = a + b·ch + c·ch² to selected histogram
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnApplyCalibration()
{
    if (!rawHist_) { AppendLog("Load a histogram first."); return; }
    if (!calibA_ || !calibB_ || !calibC_) return;

    double aParam = calibA_->GetNumber();
    double bParam = calibB_->GetNumber();
    double cParam = calibC_->GetNumber();

    if (std::abs(bParam) < 1e-10 && std::abs(cParam) < 1e-20) {
        AppendLog("Calibration requires b ≠ 0 (slope in keV/channel).");
        return;
    }

    TH1* src  = rawHist_;
    int  nbin = src->GetNbinsX();
    std::string newName = currentHist_ + "_cal";

    // Build non-uniform bin edges mapped through the calibration polynomial
    std::vector<double> edges(nbin + 1);
    for (int i = 0; i <= nbin; i++) {
        double ch = (i < nbin)
            ? src->GetXaxis()->GetBinLowEdge(i + 1)
            : src->GetXaxis()->GetBinUpEdge(nbin);
        edges[i] = aParam + bParam * ch + cParam * ch * ch;
    }

    if (edges.front() > edges.back())
        AppendLog("WARNING: negative slope — x-axis will be reversed.");

    TH1D* cal = new TH1D(newName.c_str(),
                         (std::string(src->GetTitle()) + " (keV)").c_str(),
                         nbin, edges.data());
    cal->SetDirectory(nullptr);
    cal->GetXaxis()->SetTitle("Energy (keV)");
    cal->GetYaxis()->SetTitle(src->GetYaxis()->GetTitle());
    cal->Sumw2();

    for (int i = 1; i <= nbin; i++) {
        cal->SetBinContent(i, src->GetBinContent(i));
        cal->SetBinError(i, src->GetBinError(i));
    }

    histNames_.push_back(newName);
    histClass_[newName] = histClass_.count(currentHist_) ? histClass_[currentHist_] : "Gamma Spectrum";

    if (rawHistOwned_) { delete rawHist_; rawHistOwned_ = false; }
    rawHist_     = cal;
    rawHistOwned_ = true;
    currentHist_ = newName;
    PopulateHistWidgets();
    DrawOnCanvas(rawHist_);

    AppendLog(Form("Calibration applied: E = %.4g + %.6g·ch + %.4g·ch²", aParam, bParam, cParam));
    AppendLog("  Created '" + newName + "'  range ["
              + Fmt(cal->GetXaxis()->GetXmin(), 1) + ", "
              + Fmt(cal->GetXaxis()->GetXmax(), 1) + "] keV");
}

// ─────────────────────────────────────────────────────────────────────────────
// OnApplyEfficiency — store efficiency model parameters for peak stats display
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnApplyEfficiency()
{
    if (!effA_ || !effB_ || !effC_ || !effD_) return;

    double a = effA_->GetNumber(), b = effB_->GetNumber();
    double c = effC_->GetNumber(), d = effD_->GetNumber();

    if (b == 0.0 && c == 0.0 && d == 0.0) {
        AppendLog("Efficiency: all non-offset parameters are zero — model has no E-dependence.");
        return;
    }
    AppendLog(Form("Efficiency model: ln(ε) = %.4g - %.4g·ln(E) + %.4g·ln(E)² - %.4g/E²",
                   a, b, c, d));
    AppendLog("Eff-corrected areas will appear in peak statistics after the next fit.");

    // Refresh stats display if a fit is already shown
    if (manualTF1_ && rawHist_) {
        TH1* dispHist = viewHist_ ? viewHist_ : rawHist_;
        UpdatePeakStats(manualTF1_, dispHist,
                        manualTF1_->GetXmin(), manualTF1_->GetXmax());
    }
}

// ─────────────────────────────────────────────────────────────────────────────

static std::string RecentFilePath()
{
    const char* home = gSystem->HomeDirectory();
    return std::string(home ? home : ".") + "/.autogammafit_recent.txt";
}

static constexpr int kMaxRecent = 10;

void GammaFitGUI::LoadRecentFiles()
{
    recentFiles_.clear();
    std::ifstream in(RecentFilePath());
    std::string line;
    while (std::getline(in, line)) {
        // Strip \r and trailing spaces so \r\n files don't corrupt the path
        while (!line.empty() && (line.back() == '\r' || line.back() == ' '))
            line.pop_back();
        if (!line.empty())
            recentFiles_.push_back(line);
    }
}

void GammaFitGUI::SaveRecentFiles() const
{
    std::ofstream out(RecentFilePath());
    for (const auto& p : recentFiles_)
        out << p << "\n";
}

void GammaFitGUI::AddToRecentFiles(const std::string& path)
{
    if (path.empty()) return;
    // Remove duplicate if already present
    recentFiles_.erase(
        std::remove(recentFiles_.begin(), recentFiles_.end(), path),
        recentFiles_.end());
    recentFiles_.insert(recentFiles_.begin(), path);
    if ((int)recentFiles_.size() > kMaxRecent)
        recentFiles_.resize(kMaxRecent);
    SaveRecentFiles();
    RefreshRecentCombo();
}

void GammaFitGUI::RefreshRecentCombo()
{
    if (!recentCombo_) return;
    recentCombo_->RemoveAll();
    if (recentFiles_.empty()) {
        recentCombo_->AddEntry("(no recent files)", 1);
        recentCombo_->Select(1, kFALSE);
    } else {
        for (int i = 0; i < (int)recentFiles_.size(); ++i) {
            // Show only the filename portion, keep full path as title via id
            std::string display = recentFiles_[i];
            auto slash = display.rfind('/');
            if (slash != std::string::npos) display = display.substr(slash + 1);
            recentCombo_->AddEntry(display.c_str(), i + 1);
        }
        recentCombo_->Select(1, kFALSE);
    }
    recentCombo_->MapSubwindows();
    recentCombo_->Layout();
}

void GammaFitGUI::OnOpenRecent()
{
    if (recentFiles_.empty()) { AppendLog("No recent files."); return; }
    int sel = recentCombo_ ? recentCombo_->GetSelected() : -1;
    if (sel < 1 || sel > (int)recentFiles_.size()) {
        AppendLog("Select a recent file from the dropdown first."); return;
    }
    std::string path = recentFiles_[sel - 1];

    if (inputFile_) { inputFile_->Close(); delete inputFile_; inputFile_ = nullptr; }
    inputPath_ = path;
    inputFile_ = TFile::Open(inputPath_.c_str(), "READ");
    if (!inputFile_ || inputFile_->IsZombie()) {
        AppendLog("ERROR: Cannot open " + inputPath_);
        inputFile_ = nullptr;
        return;
    }

    AddToRecentFiles(inputPath_);

    histNames_.clear(); th2Names_.clear(); projParent_.clear();
    histClass_.clear(); customProjDefs_.clear(); bgSubtractDefs_.clear();
    rebinFactors_.clear();

    TIter next(inputFile_->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        if (!obj) continue;
        std::string name = obj->GetName();
        if (obj->InheritsFrom("TH2")) {
            histNames_.push_back(name);
            th2Names_.insert(name);
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
    PopulateBgSubCombos();

    fittedHists_.clear();
    fitResultsList_->RemoveAll();
    for (const auto& hname : histNames_) {
        struct stat st;
        if (stat(CacheFileFor(hname).c_str(), &st) == 0)
            fittedHists_.push_back(hname);
    }
    for (size_t i = 0; i < fittedHists_.size(); ++i)
        fitResultsList_->AddEntry(fittedHists_[i].c_str(), (Int_t)i + 1);
    fitResultsList_->MapSubwindows(); fitResultsList_->Layout();

    std::string display = inputPath_;
    if (display.size() > 45) display = "..." + display.substr(display.size() - 42);
    fileLbl_->SetText(display.c_str());
    AppendLog("Opened (recent): " + inputPath_);
    SetStatus("Loaded: " + inputPath_);
}
