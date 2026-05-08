#include "GammaFitGUI.h"

#include "TGFileDialog.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TFile.h"
#include "TKey.h"
#include "TF1.h"
#include "TLine.h"
#include "TROOT.h"
#include "TFitResult.h"

#include "PeakFitter.h"

#include <iostream>
#include <sstream>
#include <iomanip>
#include <limits>
#include <sys/stat.h>

static const char* kCacheDir           = "fit_caches";
static const char* kIsotopeDBDefault   = "../Isotope_energys.txt";

static std::string Fmt(double v, int n = 3) {
    std::ostringstream ss;
    ss << std::fixed << std::setprecision(n) << v;
    return ss.str();
}

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

    BuildAutoFitTab  (ctrlTab->AddTab("AutoFit"));
    BuildManualFitTab(ctrlTab->AddTab("Manual Fit"));
    // ── Add more tabs here: ctrlTab->AddTab("Analysis") → BuildAnalysisTab()

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
        AppendLog("Isotope DB loaded: " + isotopePath_);
        isotopeLbl_->SetText(isotopePath_.c_str());
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
    if (inputFile_) { inputFile_->Close(); delete inputFile_; }
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
    p->AddFrame(hg, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 4, 4, 2, 2));

    histList_ = new TGListBox(hg, 100);
    histList_->Resize(285, 130);
    hg->AddFrame(histList_, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
    histList_->Connect("Selected(Int_t)", "GammaFitGUI", this, "OnHistogramSelected(Int_t)");
    histList_->SetMultipleSelections(kFALSE);

    // ── Fit options ───────────────────────────────────────────────────────────
    TGGroupFrame* og = new TGGroupFrame(p, "Fit Options");
    p->AddFrame(og, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    cacheOnlyChk_ = new TGCheckButton(og, "Cache-Only  (skip MIGRAD)");
    og->AddFrame(cacheOnlyChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 1));
    cacheOnlyChk_->SetToolTipText(
        "Use previously saved fit parameters directly.\n"
        "Histograms with no cache entry are skipped.");

    useSeedsChk_ = new TGCheckButton(og, "Use Cached Seeds");
    useSeedsChk_->SetState(kButtonDown);
    og->AddFrame(useSeedsChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 1, 4));
    useSeedsChk_->SetToolTipText(
        "Warm-start MIGRAD from the best previously converged parameters.\n"
        "Uncheck to always fit from scratch.");

    // ── Debug toggles ─────────────────────────────────────────────────────────
    TGGroupFrame* dg = new TGGroupFrame(p, "Debug Sections");
    p->AddFrame(dg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

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

    // ── Run ───────────────────────────────────────────────────────────────────
    TGGroupFrame* rg = new TGGroupFrame(p, "Run");
    p->AddFrame(rg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

    TGTextButton* runSel = new TGTextButton(rg, "Run AutoFit  (selected)");
    rg->AddFrame(runSel, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    runSel->Connect("Clicked()", "GammaFitGUI", this, "OnRunSelected()");
    runSel->SetToolTipText("Fit the histogram currently selected in the list");

    TGTextButton* runAll = new TGTextButton(rg, "Run AutoFit  (ALL histograms)");
    rg->AddFrame(runAll, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
    runAll->Connect("Clicked()", "GammaFitGUI", this, "OnRunAll()");
    runAll->SetToolTipText("Fit every histogram in the file in sequence");
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildManualFitTab
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::BuildManualFitTab(TGCompositeFrame* p)
{
    // ── Histogram selector ────────────────────────────────────────────────────
    TGGroupFrame* hg = new TGGroupFrame(p, "Histogram");
    p->AddFrame(hg, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

    manualCombo_ = new TGComboBox(hg, 300);
    manualCombo_->Resize(285, 22);
    hg->AddFrame(manualCombo_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    TGTextButton* loadBtn = new TGTextButton(hg, "Load to Canvas");
    hg->AddFrame(loadBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));
    loadBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLoadManual()");
    loadBtn->SetToolTipText("Draw the selected histogram. Then click on a peak to place the fit marker.");

    TGLabel* inst = new TGLabel(p, "  ◀ Click on spectrum to place peak marker");
    p->AddFrame(inst, new TGLayoutHints(kLHintsLeft, 4, 4, 0, 2));

    // ── Peak parameters ───────────────────────────────────────────────────────
    TGGroupFrame* pg = new TGGroupFrame(p, "Peak Parameters");
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

    addRow(pg, "Energy (keV)", mEnergy_, 500.0,  0.0,  10000.0);
    addRow(pg, "Sigma  (keV)", mSigma_,    1.5,  0.05,    50.0);
    addRow(pg, "Amplitude",    mAmp_,    500.0,  0.0,    1.0e9);
    addRow(pg, "BG const",     mBg0_,      0.0, -1.0e6,  1.0e6);
    addRow(pg, "BG slope",     mBg1_,      0.0, -1.0e4,  1.0e4);
    addRow(pg, "Range (×σ)",   mRange_,    4.0,  1.0,     20.0);

    TGTextButton* seedBtn = new TGTextButton(pg, "Seed from Resolution Model");
    pg->AddFrame(seedBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 6, 2));
    seedBtn->Connect("Clicked()", "GammaFitGUI", this, "OnSeedParams()");
    seedBtn->SetToolTipText(
        "Fill Sigma from the current resolution model at the energy above.\n"
        "Fill Amplitude from the histogram bin height.");

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

    TGHorizontalFrame* btnRow2 = new TGHorizontalFrame(ag);
    ag->AddFrame(btnRow2, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));

    TGTextButton* accBtn = new TGTextButton(btnRow2, "Accept & Save to Cache");
    btnRow2->AddFrame(accBtn, new TGLayoutHints(kLHintsExpandX, 2, 6, 0, 0));
    accBtn->Connect("Clicked()", "GammaFitGUI", this, "OnAcceptFit()");
    accBtn->SetToolTipText("Write the current fit parameters to this histogram's cache file");

    TGTextButton* rejBtn = new TGTextButton(btnRow2, "Reject");
    btnRow2->AddFrame(rejBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
    rejBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRejectFit()");
    rejBtn->SetToolTipText("Discard the current manual fit and redraw the histogram");

    // ── Extension point ───────────────────────────────────────────────────────
    // To add a new analysis tool: create a TGGroupFrame here and connect buttons
    // to new slot methods declared in GammaFitGUI.h.
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
    TIter next(inputFile_->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        TObject* obj = key->ReadObj();
        if (obj && obj->InheritsFrom("TH1"))
            histNames_.push_back(obj->GetName());
    }

    PopulateHistWidgets();

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
        AppendLog("Isotope DB loaded: " + isotopePath_);
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
    for (size_t i = 0; i < histNames_.size(); ++i) {
        histList_->AddEntry(histNames_[i].c_str(), (Int_t)i + 1);
        manualCombo_->AddEntry(histNames_[i].c_str(), (Int_t)i + 1);
    }
    histList_->MapSubwindows();   histList_->Layout();
    manualCombo_->MapSubwindows(); manualCombo_->Layout();
}

void GammaFitGUI::OnHistogramSelected(Int_t id)
{
    if (!inputFile_ || id < 1 || (size_t)id > histNames_.size()) return;
    currentHist_ = histNames_[id - 1];

    rawHist_ = (TH1*)inputFile_->Get(currentHist_.c_str());
    DrawOnCanvas(rawHist_);
    SetStatus("Selected: " + currentHist_);
}

void GammaFitGUI::OnLoadManual()
{
    Int_t id = manualCombo_->GetSelected();
    if (id < 1) { AppendLog("Select a histogram from the dropdown first."); return; }
    OnHistogramSelected(id);
    AppendLog("Loaded to canvas: " + currentHist_ +
              "  — click on a peak to place the fit marker");
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

void GammaFitGUI::RunFitOnHistogram(const std::string& hname)
{
    if (!inputFile_) { AppendLog("No ROOT file loaded."); return; }

    SyncDebugToggles();
    mkdir(kCacheDir, 0755);

    bool cacheOnly = cacheOnlyChk_->IsOn();
    bool useSeeds  = useSeedsChk_->IsOn();

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(hname));
    fitdb.cacheOnly      = cacheOnly;
    fitdb.useCachedSeeds = useSeeds;

    // Re-create tracker and storage per run so residuals don't accumulate
    PeakTracker runTracker;
    FitStorage  runStorage;
    PeakFitter  fitter(db_, &runTracker, res_, runStorage, &fitdb);

    TH1* h = (TH1*)inputFile_->Get(hname.c_str());
    if (!h) { AppendLog("ERROR: histogram '" + hname + "' not found"); return; }

    const std::string outPath = "AllGammaFits.root";
    TFile* fout = TFile::Open(outPath.c_str(), "UPDATE");
    if (!fout || fout->IsZombie()) {
        delete fout;
        fout = TFile::Open(outPath.c_str(), "RECREATE");
    }

    AppendLog("Running AutoFit: " + hname +
              (cacheOnly ? "  [cache-only]" : "") +
              (useSeeds  ? "  [seeds on]"   : "  [fresh start]"));
    SetStatus("Fitting: " + hname + " ...");

    // Pass the embedded canvas so the fit is drawn live in the GUI
    fitter.FitHistogram(h, fout, true, canvas_->GetCanvas());

    canvas_->GetCanvas()->Modified();
    canvas_->GetCanvas()->Update();

    if (!cacheOnly) fitdb.Save(CacheFileFor(hname));
    fout->Close();
    delete fout;

    // Update the resolution model from this run's tracker
    res_ = runTracker.GetData().empty() ? res_ : res_;  // model updated inside FitHistogram via tracker

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
    AppendLog("=== AutoFit ALL: " +
              std::to_string(histNames_.size()) + " histograms ===");
    for (const auto& name : histNames_)
        RunFitOnHistogram(name);
    AppendLog("=== All histograms complete ===");
}

// ─────────────────────────────────────────────────────────────────────────────
// Manual Fit
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnCanvasEvent(Int_t event, Int_t px, Int_t /*py*/, TObject* /*obj*/)
{
    if (event != kButton1Down) return;
    if (!rawHist_) { AppendLog("Load a histogram first (Manual Fit tab → Load to Canvas)."); return; }

    TCanvas* c = canvas_->GetCanvas();
    c->cd();

    // Convert pixel x to the histogram's energy axis
    double energy = c->PadtoX(c->AbsPixeltoX(px));

    double xmin = rawHist_->GetXaxis()->GetXmin();
    double xmax = rawHist_->GetXaxis()->GetXmax();
    if (energy < xmin || energy > xmax) return;

    clickEnergy_ = energy;
    mEnergy_->SetNumber(energy);

    // Draw a dashed vertical marker at the clicked energy
    double ylo = c->GetUymin();
    double yhi = c->GetUymax() * 0.95;
    TLine* line = new TLine(energy, ylo, energy, yhi);
    line->SetLineColor(kBlue + 1);
    line->SetLineStyle(2);
    line->SetLineWidth(2);
    line->Draw("same");
    c->Modified();
    c->Update();

    OnSeedParams();

    AppendLog("Peak marker: " + Fmt(energy, 2) + " keV  — adjust parameters then click Run Fit");
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

    double E    = mEnergy_->GetNumber();
    double sig  = mSigma_->GetNumber();
    double A    = mAmp_->GetNumber();
    double bg0  = mBg0_->GetNumber();
    double bg1  = mBg1_->GetNumber();
    double nSig = mRange_->GetNumber();

    delete manualTF1_;
    manualTF1_ = new TF1("manual_preview",
        "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x",
        E - nSig * sig, E + nSig * sig);
    manualTF1_->SetParameters(A, E, sig, bg0, bg1);
    manualTF1_->SetLineColor(kOrange + 1);
    manualTF1_->SetLineStyle(2);
    manualTF1_->SetLineWidth(2);

    TCanvas* c = canvas_->GetCanvas();
    c->cd();
    manualTF1_->Draw("same");
    c->Modified();
    c->Update();

    AppendLog("Preview: E=" + Fmt(E,2) + "  sig=" + Fmt(sig) + "  A=" + Fmt(A,1));
}

void GammaFitGUI::OnManualFit()
{
    if (!rawHist_) { AppendLog("Load a histogram first."); return; }

    double E    = mEnergy_->GetNumber();
    double sig  = mSigma_->GetNumber();
    double A    = mAmp_->GetNumber();
    double bg0  = mBg0_->GetNumber();
    double bg1  = mBg1_->GetNumber();
    double nSig = mRange_->GetNumber();
    double xmin = E - nSig * sig;
    double xmax = E + nSig * sig;

    delete manualTF1_;
    manualTF1_ = new TF1("manual_fit",
        "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]+[4]*x", xmin, xmax);
    manualTF1_->SetParameters(A, E, sig, bg0, bg1);
    manualTF1_->SetParLimits(0, std::max(A * 0.01, 1.0), A * 20.0);
    manualTF1_->SetParLimits(1, E - 8.0,  E + 8.0);
    manualTF1_->SetParLimits(2, sig * 0.2, sig * 4.0);

    AppendLog("Running manual fit at E=" + Fmt(E,2) +
              "  window [" + Fmt(xmin,1) + ", " + Fmt(xmax,1) + "] keV");

    TFitResultPtr r = rawHist_->Fit(manualTF1_, "R S Q L B");

    // Update entry boxes with fitted values
    mEnergy_->SetNumber(manualTF1_->GetParameter(1));
    mSigma_ ->SetNumber(manualTF1_->GetParameter(2));
    mAmp_   ->SetNumber(manualTF1_->GetParameter(0));
    mBg0_   ->SetNumber(manualTF1_->GetParameter(3));
    mBg1_   ->SetNumber(manualTF1_->GetParameter(4));

    double chi2ndf = (r.Get() && r->Ndf() > 0) ? r->Chi2() / r->Ndf() : -1.0;
    int    status  = r.Get() ? r->Status() : -1;

    std::string result =
        "E="    + Fmt(manualTF1_->GetParameter(1), 2) +
        " keV  sig=" + Fmt(manualTF1_->GetParameter(2)) +
        "  FWHM=" + Fmt(2.355 * manualTF1_->GetParameter(2), 2) +
        "  chi2/ndf=" + Fmt(chi2ndf, 2) +
        "  status=" + std::to_string(status);

    mResultLbl_->SetText(result.c_str());
    AppendLog("Manual fit result: " + result);

    manualTF1_->SetLineColor(kRed);
    manualTF1_->SetLineWidth(2);

    TCanvas* c = canvas_->GetCanvas();
    c->cd();
    manualTF1_->Draw("same");
    c->Modified();
    c->Update();

    SetStatus("Manual fit: " + result);
}

void GammaFitGUI::OnAcceptFit()
{
    if (!manualTF1_) {
        AppendLog("No fit to save — run a manual fit first.");
        return;
    }
    if (currentHist_.empty()) {
        AppendLog("No histogram selected.");
        return;
    }

    double E = manualTF1_->GetParameter(1);
    std::string key = FitDatabase::MakeKey({E});
    std::string cacheFile = CacheFileFor(currentHist_);

    FitDatabase fitdb;
    fitdb.Load(cacheFile);

    FitEntry e;
    e.key = key;
    int n = manualTF1_->GetNpar();
    e.params.resize(n);
    for (int i = 0; i < n; ++i)
        e.params[i] = manualTF1_->GetParameter(i);
    // Use max chi2ndf so the auto-fitter will replace this if it does better
    e.chi2ndf = std::numeric_limits<double>::max();

    fitdb.StoreIfBetter(key, e);
    mkdir(kCacheDir, 0755);
    fitdb.Save(cacheFile);

    AppendLog("Saved to cache: key=" + key + "  file=" + cacheFile);
    SetStatus("Manual fit saved to cache.");
}

void GammaFitGUI::OnRejectFit()
{
    delete manualTF1_;
    manualTF1_ = nullptr;
    mResultLbl_->SetText("Fit rejected.");
    DrawOnCanvas(rawHist_);
    AppendLog("Manual fit rejected — canvas cleared.");
}

// ─────────────────────────────────────────────────────────────────────────────
// Utilities
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::DrawOnCanvas(TH1* h, TF1* fit)
{
    if (!h) return;
    TCanvas* c = canvas_->GetCanvas();
    c->cd();
    c->Clear();
    h->SetLineColor(kBlack);
    h->Draw("hist");
    if (fit) {
        fit->SetLineColor(kRed);
        fit->SetLineWidth(2);
        fit->Draw("same");
    }
    c->Modified();
    c->Update();
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
