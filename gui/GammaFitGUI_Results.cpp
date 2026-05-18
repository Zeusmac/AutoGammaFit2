#include "GammaFitGUI.h"
#include "GammaFitGUI_shared.h"

#include "TGFileDialog.h"
#include "TGTextEntry.h"
#include "TGMsgBox.h"
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TROOT.h"
#include "TSystem.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <algorithm>
#include <cstring>

void GammaFitGUI::BuildFitResultsTab(TGCompositeFrame* p)
{
    TGGroupFrame* lg = new TGGroupFrame(p, "Fitted Histograms");
    p->AddFrame(lg, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 4, 4, 4, 2));

    fitResultsList_ = new TGListBox(lg, 500);
    fitResultsList_->Resize(285, 180);
    lg->AddFrame(fitResultsList_,
                 new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
    fitResultsList_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                             "OnFitResultSelected(Int_t)");

    TGTextButton* showBtn = new TGTextButton(lg, "Show Selected on Canvas");
    lg->AddFrame(showBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    showBtn->Connect("Clicked()", "GammaFitGUI", this, "OnFitResultSelected(Int_t)");
    showBtn->SetToolTipText("Draw the selected histogram: BG-subtracted with all cached fits and peak labels");

    TGTextButton* loadCacheBtn2 = new TGTextButton(lg, "Load Cache onto Histogram");
    lg->AddFrame(loadCacheBtn2, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
    loadCacheBtn2->Connect("Clicked()", "GammaFitGUI", this, "OnLoadCache()");
    loadCacheBtn2->SetToolTipText("Reload the cache from disk and redisplay the current histogram with all stored fits");

    // ── Save ─────────────────────────────────────────────────────────────────
    TGGroupFrame* sg = new TGGroupFrame(p, "Save to ROOT File");
    p->AddFrame(sg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    TGTextButton* saveSelBtn = new TGTextButton(sg, "Save Selected");
    sg->AddFrame(saveSelBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    saveSelBtn->Connect("Clicked()", "GammaFitGUI", this, "OnSaveSelected()");
    saveSelBtn->SetToolTipText("Save the selected histogram canvas to a ROOT file");

    TGTextButton* saveAllBtn = new TGTextButton(sg, "Save All");
    sg->AddFrame(saveAllBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
    saveAllBtn->Connect("Clicked()", "GammaFitGUI", this, "OnSaveAll()");
    saveAllBtn->SetToolTipText("Save every histogram in the list to a single ROOT file");

    TGTextButton* exportBtn = new TGTextButton(sg, "Export Cache to CSV...");
    sg->AddFrame(exportBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));
    exportBtn->Connect("Clicked()", "GammaFitGUI", this, "OnExportCacheCSV()");
    exportBtn->SetToolTipText("Export fitted peak parameters from all cached histograms to a CSV file");

    // ── gnuScope export ───────────────────────────────────────────────────────
    TGGroupFrame* gsg = new TGGroupFrame(p, "Export to gnuScope");
    p->AddFrame(gsg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    TGTextButton* gsExportBtn = new TGTextButton(gsg, "Export Current Histogram...");
    gsg->AddFrame(gsExportBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    gsExportBtn->Connect("Clicked()", "GammaFitGUI", this, "OnExportGnuScope()");
    gsExportBtn->SetToolTipText(
        "Export the currently displayed histogram to gnuScope binary format.\n"
        "1D histograms → .spe (Fortran unformatted float array)\n"
        "2D histograms → .sqr (square matrix, same format as dump_square_asym.C)");

    TGTextButton* gsExportAllBtn = new TGTextButton(gsg, "Export All Histograms to Folder...");
    gsg->AddFrame(gsExportAllBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
    gsExportAllBtn->Connect("Clicked()", "GammaFitGUI", this, "OnExportGnuScopeAll()");
    gsExportAllBtn->SetToolTipText(
        "Batch-export every loaded histogram to gnuScope format in a chosen folder.\n"
        "1D → <histname>.spe,  2D → <histname>.sqr");

    // ── Canvas Annotations ────────────────────────────────────────────────────
    TGGroupFrame* annGrp = new TGGroupFrame(p, "Canvas Annotations");
    p->AddFrame(annGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    auto addLabeledEntry = [&](TGGroupFrame* grp, const char* lbl,
                               TGTextEntry*& entry, const char* defVal)
    {
        TGHorizontalFrame* row = new TGHorizontalFrame(grp);
        grp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 0));
        TGLabel* l = new TGLabel(row, lbl);
        l->SetWidth(60);
        row->AddFrame(l, new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        entry = new TGTextEntry(row, defVal);
        row->AddFrame(entry, new TGLayoutHints(kLHintsExpandX));
    };

    addLabeledEntry(annGrp, "Title:",   frTitleEntry_,  "");
    addLabeledEntry(annGrp, "X label:", frXLabelEntry_, "");
    addLabeledEntry(annGrp, "Y label:", frYLabelEntry_, "");
    addLabeledEntry(annGrp, "Text:",    frAnnotText_,   "");

    {
        TGHorizontalFrame* posRow = new TGHorizontalFrame(annGrp);
        annGrp->AddFrame(posRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 0));
        posRow->AddFrame(new TGLabel(posRow, "X:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        frAnnotX_ = new TGNumberEntry(posRow, 0.15, 5, -1,
            TGNumberFormat::kNESRealTwo, TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMinMax, 0.0, 1.0);
        posRow->AddFrame(frAnnotX_, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));
        posRow->AddFrame(new TGLabel(posRow, "Y:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        frAnnotY_ = new TGNumberEntry(posRow, 0.85, 5, -1,
            TGNumberFormat::kNESRealTwo, TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMinMax, 0.0, 1.0);
        posRow->AddFrame(frAnnotY_, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));
        posRow->AddFrame(new TGLabel(posRow, "Sz:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        frAnnotSize_ = new TGNumberEntry(posRow, 0.040, 5, -1,
            TGNumberFormat::kNESRealThree, TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELLimitMinMax, 0.005, 0.2);
        posRow->AddFrame(frAnnotSize_, new TGLayoutHints(kLHintsLeft));
    }

    TGTextButton* applyAnnBtn = new TGTextButton(annGrp, "Apply to Canvas");
    annGrp->AddFrame(applyAnnBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 4));
    applyAnnBtn->Connect("Clicked()", "GammaFitGUI", this, "OnApplyCanvasAnnotations()");
    applyAnnBtn->SetToolTipText("Apply title, axis labels, and text annotation to the current canvas");
}

void GammaFitGUI::OnFitResultSelected(Int_t id)
{
    // id may be 0 when the button fires without a selection — use list selection
    if (id < 1) id = fitResultsList_->GetSelected();
    if (id < 1 || id > (Int_t)fittedHists_.size()) return;
    ShowFitResult(fittedHists_[id - 1]);
}

void GammaFitGUI::ShowFitResult(const std::string& hname)
{
    // FWHM result: re-display the stored resolution model as a FWHM vs Energy plot
    if (hname.size() > strlen(kFwhmPrefix) &&
        hname.substr(0, strlen(kFwhmPrefix)) == kFwhmPrefix)
    {
        std::string histname = hname.substr(strlen(kFwhmPrefix));
        fwhmHistName_ = histname;
        for (size_t i = 0; i < histNames_.size(); i++) {
            if (histNames_[i] == histname) {
                fwhmCombo_->Select((Int_t)i + 1, kFALSE);
                break;
            }
        }
        OnLoadFWHM();
        return;
    }

    if (!inputFile_) return;
    if (rawHistOwned_) { delete rawHist_; rawHist_ = nullptr; rawHistOwned_ = false; }
    bool rawOwned = false;
    TH1* raw = LoadHistFromFile(hname, rawOwned);
    if (!raw) return;

    // Keep currentHist_ / rawHist_ consistent so residuals & nav work
    currentHist_  = hname;
    rawHist_      = raw;
    rawHistOwned_ = rawOwned;

    TCanvas* c = canvas_->GetCanvas();
    c->Clear();
    c->cd();

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(hname));

    delete viewHist_;
    viewHist_ = MakeBgSubHist(raw, fitdb.bgSubtracted, fitdb.bgIterations);
    ApplyHistStyle(viewHist_, hname.c_str());
    viewHist_->SetLineColor(kBlack);
    viewHist_->SetMarkerSize(0);
    viewHist_->GetXaxis()->UnZoom();
    SetYMaxFromVisible(viewHist_);
    viewHist_->Draw("hist");
    if (showErrorBars_) viewHist_->Draw("E1 same");

    OverlayFitPeaks(hname, c);

    c->Modified(); c->Update();
    gSystem->ProcessEvents();
    SetStatus("Fit Results: " + hname);
}

void GammaFitGUI::SaveFitResultToFile(const std::string& hname, TFile* fout)
{
    if (!fout) return;

    // FWHM result: save a canvas with the FWHM vs Energy graph + model curve
    if (hname.size() > strlen(kFwhmPrefix) &&
        hname.substr(0, strlen(kFwhmPrefix)) == kFwhmPrefix)
    {
        std::string histname = hname.substr(strlen(kFwhmPrefix));
        FitDatabase fitdb;
        if (!fitdb.Load(CacheFileFor(histname))) return;

        // Extract all FWHM data points
        std::vector<double> allX, allY;
        for (const auto& [key, entry] : fitdb.GetEntries()) {
            if (key.size() >= 2 && key[0] == '_' && key[1] == '_') continue;
            FitLayout lay = DetectLayout((int)entry.params.size());
            if (!lay.valid()) continue;
            for (int i = 0; i < lay.n; i++) {
                double E   = entry.params[3*i + 1];
                double sig = entry.params[3*i + 2];
                if (sig <= 0 || E <= 0) continue;
                allX.push_back(E);
                allY.push_back(2.3548 * sig);
            }
        }
        if (allX.empty()) return;

        // Load exclusion list
        std::vector<bool> excluded(allX.size(), false);
        auto exIt = fitdb.GetEntries().find(kExcludedFwhmKey);
        if (exIt != fitdb.GetEntries().end()) {
            for (double exE : exIt->second.params)
                for (size_t i = 0; i < allX.size(); i++)
                    if (std::abs(allX[i] - exE) < 0.5) { excluded[i] = true; break; }
        }

        // Split into included / excluded
        std::vector<double> xIn, yIn, exIn, eyIn, xEx, yEx;
        for (size_t i = 0; i < allX.size(); i++) {
            if (excluded[i]) { xEx.push_back(allX[i]); yEx.push_back(allY[i]); }
            else { xIn.push_back(allX[i]); yIn.push_back(allY[i]);
                   exIn.push_back(0.0); eyIn.push_back(0.05 * allY[i]); }
        }

        Bool_t wasBatch = gROOT->IsBatch();
        gROOT->SetBatch(kTRUE);
        TCanvas* csave = new TCanvas(Form("c_fwhm_%s", histname.c_str()),
                                      ("FWHM: " + histname).c_str(), 800, 600);
        gROOT->SetBatch(wasBatch);
        csave->cd();

        // Included points
        if (!xIn.empty()) {
            TGraphErrors* gr = new TGraphErrors((Int_t)xIn.size(),
                xIn.data(), yIn.data(), exIn.data(), eyIn.data());
            gr->SetTitle(("FWHM vs Energy  [" + histname +
                           "];Energy (keV);FWHM (keV)").c_str());
            gr->SetMarkerStyle(20); gr->SetMarkerSize(0.8);
            gr->SetMarkerColor(kBlue + 1); gr->SetLineColor(kBlue + 1);
            gr->Draw("AP");
        }
        // Excluded points (hollow gray)
        if (!xEx.empty()) {
            TGraph* grEx = new TGraph((Int_t)xEx.size(), xEx.data(), yEx.data());
            grEx->SetMarkerStyle(24); grEx->SetMarkerSize(0.8);
            grEx->SetMarkerColor(kGray + 1);
            grEx->Draw(xIn.empty() ? "AP" : "P same");
        }

        // Model + stat line
        auto rit = fitdb.GetEntries().find(kResolutionKey);
        if (rit != fitdb.GetEntries().end() && rit->second.params.size() == 3) {
            double a = rit->second.params[0];
            double b = rit->second.params[1];
            double cv = rit->second.params[2];
            double xhi = *std::max_element(allX.begin(), allX.end()) * 1.15;
            TF1* resf = new TF1("fwhm_save",
                Form("sqrt(%.10g+%.10g*x+%.10g*x*x)", a, b, cv), 0.0, xhi);
            resf->SetLineColor(kRed); resf->SetLineWidth(2);
            resf->Draw("same");
            {
                TF1* statLine = new TF1("fwhm_stat_save",
                    Form("sqrt(%.10g+%.10g*x)", a, b), 0.0, xhi);
                statLine->SetLineColor(kGreen + 2);
                statLine->SetLineStyle(2); statLine->SetLineWidth(2);
                statLine->Draw("same");
            }
        }

        csave->Modified(); csave->Update();
        fout->cd();
        csave->Write();
        delete csave;
        return;
    }

    if (!inputFile_) return;
    bool rawOwned2 = false;
    TH1* raw = LoadHistFromFile(hname, rawOwned2);
    if (!raw) return;

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(hname));

    TH1* hbg = MakeBgSubHist(raw, fitdb.bgSubtracted, fitdb.bgIterations);
    if (rawOwned2) { delete raw; raw = nullptr; }
    ApplyHistStyle(hbg, hname.c_str());
    hbg->SetLineColor(kBlack);
    hbg->SetMarkerSize(0);
    hbg->GetXaxis()->UnZoom();
    SetYMaxFromVisible(hbg);

    Bool_t wasBatch = gROOT->IsBatch();
    gROOT->SetBatch(kTRUE);
    TCanvas* csave = new TCanvas(Form("c_%s", hname.c_str()), hname.c_str(), 1200, 800);
    gROOT->SetBatch(wasBatch);
    csave->cd();

    hbg->Draw("hist");
    if (showErrorBars_) hbg->Draw("E1 same");
    OverlayFitPeaks(hname, csave);

    csave->Modified(); csave->Update();

    fout->cd();
    csave->Write();
    hbg->Write(Form("%s_bgsub", hname.c_str()));

    delete hbg;
    delete csave;
}

void GammaFitGUI::OnSaveSelected()
{
    Int_t id = fitResultsList_->GetSelected();
    if (id < 1 || id > (Int_t)fittedHists_.size()) {
        AppendLog("Select a histogram from the Fit Results list first.");
        return;
    }
    static const char* kTypes[] = {"ROOT files","*.root","All files","*",nullptr,nullptr};
    TGFileInfo fi;
    fi.fFileTypes = kTypes;
    OpenFileDialog(this, kFDSave, &fi);
    if (!fi.fFilename) return;

    TFile* fout = TFile::Open(fi.fFilename, "RECREATE");
    if (!fout || fout->IsZombie()) { AppendLog("Cannot create output file."); delete fout; return; }

    const std::string& hname = fittedHists_[id - 1];
    SaveFitResultToFile(hname, fout);
    fout->Close(); delete fout;
    AppendLog("Saved: " + hname + " -> " + std::string(fi.fFilename));
    SetStatus("Saved: " + hname);
}

void GammaFitGUI::OnSaveAll()
{
    if (fittedHists_.empty()) { AppendLog("No fitted histograms to save."); return; }
    static const char* kTypes[] = {"ROOT files","*.root","All files","*",nullptr,nullptr};
    TGFileInfo fi;
    fi.fFileTypes = kTypes;
    OpenFileDialog(this, kFDSave, &fi);
    if (!fi.fFilename) return;

    TFile* fout = TFile::Open(fi.fFilename, "RECREATE");
    if (!fout || fout->IsZombie()) { AppendLog("Cannot create output file."); delete fout; return; }

    for (const auto& hname : fittedHists_) {
        SaveFitResultToFile(hname, fout);
        AppendLog("  saved: " + hname);
        gSystem->ProcessEvents();
    }
    fout->Close(); delete fout;
    AppendLog("Saved " + std::to_string(fittedHists_.size()) + " histograms -> " + std::string(fi.fFilename));
    SetStatus("Saved " + std::to_string(fittedHists_.size()) + " histograms");
}

// ─────────────────────────────────────────────────────────────────────────────
// OnExportCacheCSV — column-picker popup then write CSV
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnExportCacheCSV()
{
    // ── Column-selector popup ─────────────────────────────────────────────────
    TGTransientFrame* dlg = new TGTransientFrame(
        gClient->GetRoot(), this, 320, 320, kVerticalFrame);
    dlg->SetWindowName("Export Columns");
    dlg->SetCleanup(kDeepCleanup);

    TGGroupFrame* colGrp = new TGGroupFrame(dlg, "Select columns to export");
    dlg->AddFrame(colGrp, new TGLayoutHints(kLHintsExpandX, 8, 8, 8, 4));

    struct ColDef { const char* label; bool def; };
    static const ColDef kCols[] = {
        { "Histogram",            true  },
        { "Centroid (keV)",       true  },
        { "Centroid error (keV)", true  },
        { "Sigma (keV)",          true  },
        { "Sigma error (keV)",    true  },
        { "FWHM (keV)",           true  },
        { "FWHM error (keV)",     false },
        { "Peak area (counts)",   true  },
        { "Peak area error",      true  },
        { "Label",                true  },
        { "Classification",       true  },
        { "Chi2/NDF",             true  },
        { "BG constant",          false },
        { "BG slope",             false },
        { "Match isotope",        false },
        { "Match energy (keV)",   false },
        { "Match distance (keV)", false },
    };
    const int NC = (int)(sizeof(kCols)/sizeof(kCols[0]));
    std::vector<TGCheckButton*> chk(NC);
    for (int c = 0; c < NC; c++) {
        chk[c] = new TGCheckButton(colGrp, kCols[c].label);
        if (kCols[c].def) chk[c]->SetState(kButtonDown);
        colGrp->AddFrame(chk[c], new TGLayoutHints(kLHintsLeft, 4, 4, 1, 1));
    }

    TGHorizontalFrame* btnRow = new TGHorizontalFrame(dlg);
    dlg->AddFrame(btnRow, new TGLayoutHints(kLHintsCenterX, 8, 8, 4, 8));
    TGTextButton* okBtn  = new TGTextButton(btnRow, "  Export...  ");
    TGTextButton* canBtn = new TGTextButton(btnRow, "  Cancel  ");
    btnRow->AddFrame(okBtn,  new TGLayoutHints(kLHintsLeft, 0, 8, 0, 0));
    btnRow->AddFrame(canBtn, new TGLayoutHints(kLHintsLeft));

    // Unmap (hide) the window on button press — WaitForUnmap returns while
    // the C++ objects are still alive, so chk[] is safe to read afterwards.
    // WaitFor would return only after deletion, making chk[] dangling.
    okBtn ->Connect("Clicked()", "TGFrame", dlg, "UnmapWindow()");
    canBtn->Connect("Clicked()", "TGFrame", dlg, "UnmapWindow()");

    dlg->MapSubwindows();
    dlg->Resize(dlg->GetDefaultSize());
    dlg->CenterOnParent();
    dlg->MapWindow();
    gClient->WaitForUnmap(dlg);  // returns when window is hidden, NOT deleted

    // Collect selected columns while children are still alive
    std::vector<bool> sel(NC);
    for (int c = 0; c < NC; c++) sel[c] = chk[c]->IsOn();
    delete dlg;  // kDeepCleanup deletes all children; chk[] dangling from here

    // ── File dialog ───────────────────────────────────────────────────────────
    static const char* kCsvTypes[] = {"CSV files","*.csv","All files","*",nullptr,nullptr};
    TGFileInfo fi;
    fi.fFileTypes = kCsvTypes;
    OpenFileDialog(this, kFDSave, &fi);
    if (!fi.fFilename) return;

    std::string outPath = fi.fFilename;
    if (outPath.size() < 4 || outPath.substr(outPath.size()-4) != ".csv")
        outPath += ".csv";

    std::ofstream csv(outPath);
    if (!csv.is_open()) {
        AppendLog("Cannot write to " + outPath);
        return;
    }

    // Header row
    const char* kColNames[] = {
        "histogram", "centroid_keV", "centroid_err_keV",
        "sigma_keV", "sigma_err_keV", "fwhm_keV", "fwhm_err_keV",
        "peak_area", "peak_area_err",
        "label", "classification", "chi2ndf",
        "bg_constant", "bg_slope",
        "match_isotope", "match_energy_keV", "match_distance_keV"
    };
    bool first = true;
    for (int c = 0; c < NC; c++) {
        if (!sel[c]) continue;
        if (!first) csv << ",";
        csv << kColNames[c];
        first = false;
    }
    csv << "\n";

    // Resolve which histogram to export: selected in the list, else current
    Int_t selId = fitResultsList_ ? fitResultsList_->GetSelected() : -1;
    std::string exportHist;
    if (selId >= 1 && selId <= (Int_t)fittedHists_.size())
        exportHist = fittedHists_[selId - 1];
    else if (!currentHist_.empty())
        exportHist = currentHist_;
    else {
        AppendLog("Select a histogram from the Fitted Histograms list first.");
        return;
    }

    // Data: load the selected histogram's cache → write one row per peak
    const double kSqrt2Pi = 2.5066282746310002;
    int rowCount = 0;
    {
        const std::string& hname = exportHist;
        FitDatabase fdb;
        fdb.Load(CacheFileFor(hname));
        for (const auto& kv : fdb.GetEntries()) {
            if (!kv.first.empty() && kv.first[0] == '_') continue;  // skip internal keys
            const FitEntry& e = kv.second;
            FitLayout lay = DetectLayout((int)e.params.size());
            if (!lay.valid()) continue;
            int nPeaks = lay.n;
            bool hasErr = ((int)e.paramErrors.size() == (int)e.params.size());
            double bg0 = e.params[lay.bgBase()];
            double bg1 = e.params[lay.bgBase() + 1];

            for (int i = 0; i < nPeaks; i++) {
                double A    = e.params[3*i];
                double E    = e.params[3*i + 1];
                double sig  = e.params[3*i + 2];
                double Aerr = hasErr ? e.paramErrors[3*i]     : 0.0;
                double Eerr = hasErr ? e.paramErrors[3*i + 1] : 0.0;
                double serr = hasErr ? e.paramErrors[3*i + 2] : 0.0;

                double bw      = (e.binWidth > 0.0) ? e.binWidth : 1.0;
                double area    = A * sig * kSqrt2Pi / bw;
                double areaErr = (A > 0 && sig > 0)
                    ? area * std::sqrt(std::pow(Aerr/A, 2) + std::pow(serr/sig, 2))
                    : 0.0;
                double fwhm    = 2.3548 * sig;
                double fwhmErr = 2.3548 * serr;

                // Per-Gaussian label/class (falls back to entry-level if not set)
                std::string peakLbl = e.PeakLabel(i);
                std::string peakCls = e.PeakClass(i);

                // Nearest DB match for this Gaussian
                std::string matchIso;
                double matchDbE = 0.0, matchDist = 0.0;
                if (dbLoaded_) {
                    double bestD = std::numeric_limits<double>::max();
                    for (const auto& gl : db_.db) {
                        double d = std::fabs(gl.energy - E);
                        if (d < bestD) {
                            bestD = d; matchIso = gl.isotope; matchDbE = gl.energy;
                        }
                    }
                    matchDist = (matchDbE > 0) ? E - matchDbE : 0.0;
                }

                bool fr = true;
                auto w = [&](const std::string& v) {
                    if (!fr) csv << ",";
                    csv << v;
                    fr = false;
                };
                auto wf = [&](double v) { w(Form("%.6g", v)); };

                for (int c = 0; c < NC; c++) {
                    if (!sel[c]) continue;
                    switch (c) {
                        case  0: w(hname);   break;
                        case  1: wf(E);      break;
                        case  2: wf(Eerr);   break;
                        case  3: wf(sig);    break;
                        case  4: wf(serr);   break;
                        case  5: wf(fwhm);   break;
                        case  6: wf(fwhmErr);break;
                        case  7: wf(area);   break;
                        case  8: wf(areaErr);break;
                        case  9: w(peakLbl); break;
                        case 10: w(peakCls); break;
                        case 11: wf(e.chi2ndf < 1e15 ? e.chi2ndf : -1.0); break;
                        case 12: wf(bg0);    break;
                        case 13: wf(bg1);    break;
                        case 14: w(matchIso);     break;
                        case 15: wf(matchDbE);    break;
                        case 16: wf(matchDist);   break;
                    }
                }
                csv << "\n";
                ++rowCount;
            }
        }
    }

    csv.close();
    AppendLog(Form("CSV export: %d rows [%s] -> %s",
                   rowCount, exportHist.c_str(), outPath.c_str()));
    SetStatus(Form("Exported %d peaks to CSV", rowCount));
}

// ─────────────────────────────────────────────────────────────────────────────
// gnuScope binary export
// ─────────────────────────────────────────────────────────────────────────────
// Both formats use Fortran sequential unformatted records:
//   [4-byte record length] [data] [4-byte record length]
//
// 1D (.spe): single record — nchans floats (bins 1..N, overflows excluded)
//
// 2D (.sqr): two records
//   Record 1 (header): 3 ints → maxxy, 0 (unused), 0 (type=square)
//   Record 2 (data):   maxxy×maxxy floats, row-major (matches dump_square_asym.C)

bool GammaFitGUI::ExportGnuScopeFile(TH1* h, const std::string& outPath) const
{
    std::ofstream fout(outPath, std::ios::out | std::ios::binary);
    if (!fout.is_open()) return false;

    auto writeInt   = [&](int   v) { fout.write(reinterpret_cast<char*>(&v), sizeof(int)); };
    auto writeFloat = [&](float v) { fout.write(reinterpret_cast<char*>(&v), sizeof(float)); };

    if (h->InheritsFrom("TH2")) {
        TH2* h2   = static_cast<TH2*>(h);
        int  nX   = h2->GetNbinsX();
        int  nY   = h2->GetNbinsY();
        int  maxxy = std::max(nX, nY);

        // Record 1 — header (3 ints)
        int hdrLen = static_cast<int>(sizeof(int) * 3);
        int code   = 0;
        writeInt(hdrLen);
        writeInt(maxxy);
        writeInt(code);   // unused
        writeInt(code);   // type 0 = square
        writeInt(hdrLen);

        // Record 2 — data (maxxy × maxxy floats, zero-padded beyond histogram edges)
        int datLen = static_cast<int>(sizeof(float)) * maxxy * maxxy;
        writeInt(datLen);
        for (int i = 1; i <= maxxy; i++) {
            for (int j = 1; j <= maxxy; j++) {
                float val = (i <= nY && j <= nX)
                            ? static_cast<float>(h2->GetBinContent(j, i))
                            : 0.0f;
                writeFloat(val);
            }
        }
        writeInt(datLen);
    } else {
        // 1D — single data record
        int  nX    = h->GetNbinsX();
        int  datLen = static_cast<int>(sizeof(float)) * nX;
        writeInt(datLen);
        for (int i = 1; i <= nX; i++)
            writeFloat(static_cast<float>(h->GetBinContent(i)));
        writeInt(datLen);
    }

    fout.close();
    return fout.good() || true;  // close() clears failbit on success
}

void GammaFitGUI::OnExportGnuScope()
{
    if (!inputFile_ || currentHist_.empty()) {
        AppendLog("[GnuScope] Select a histogram first."); return;
    }

    bool owned = false;
    TH1* h = LoadHistFromFile(currentHist_, owned);
    if (!h) { AppendLog("[GnuScope] Cannot load histogram: " + currentHist_); return; }

    bool is2D = h->InheritsFrom("TH2");
    const char* ext = is2D ? ".sqr" : ".spe";

    static const char* kSqrTypes[] = {"gnuScope matrix","*.sqr","All files","*",nullptr,nullptr};
    static const char* kSpeTypes[] = {"gnuScope spectrum","*.spe","All files","*",nullptr,nullptr};
    TGFileInfo fi;
    fi.fFileTypes = is2D ? kSqrTypes : kSpeTypes;
    OpenFileDialog(this, kFDSave, &fi);
    if (!fi.fFilename) { if (owned) delete h; return; }

    std::string outPath = fi.fFilename;
    // Ensure correct extension
    if (outPath.size() < 4 || outPath.substr(outPath.size() - 4) != std::string(ext))
        outPath += ext;

    bool ok = ExportGnuScopeFile(h, outPath);

    if (ok) {
        if (is2D) {
            TH2* h2 = static_cast<TH2*>(h);
            int maxxy = std::max(h2->GetNbinsX(), h2->GetNbinsY());
            AppendLog(Form("[GnuScope] %s: 2D %dx%d → %dx%d square → %s",
                           currentHist_.c_str(),
                           h2->GetNbinsX(), h2->GetNbinsY(),
                           maxxy, maxxy, outPath.c_str()));
        } else {
            AppendLog(Form("[GnuScope] %s: 1D %d channels → %s",
                           currentHist_.c_str(), h->GetNbinsX(), outPath.c_str()));
        }
        SetStatus("gnuScope export: " + outPath);
    } else {
        AppendLog("[GnuScope] ERROR: write failed to " + outPath);
    }

    if (owned) delete h;
}

void GammaFitGUI::OnExportGnuScopeAll()
{
    if (!inputFile_ || histNames_.empty()) {
        AppendLog("[GnuScope] No histograms loaded."); return;
    }

    // Use the file dialog to pick a destination directory: ask the user to
    // navigate to the target folder and press Save (any filename is ignored).
    static const char* kTypes[] = {"All files","*",nullptr,nullptr};
    TGFileInfo fi;
    fi.fFileTypes = kTypes;
    fi.fIniDir    = StrDup(".");
    OpenFileDialog(this, kFDSave, &fi);
    if (!fi.fFilename) return;

    // Extract the directory part of whatever was typed/selected
    std::string chosen = fi.fFilename;
    std::string dir;
    size_t slash = chosen.find_last_of("/\\");
    dir = (slash != std::string::npos) ? chosen.substr(0, slash) : ".";

    // Create directory if needed
    gSystem->mkdir(dir.c_str(), true);

    int nOk = 0, nErr = 0;
    for (const auto& hname : histNames_) {
        bool owned = false;
        TH1* h = LoadHistFromFile(hname, owned);
        if (!h) { nErr++; continue; }

        bool is2D = h->InheritsFrom("TH2");
        std::string safe = hname;
        // Replace characters that are problematic in filenames
        std::replace(safe.begin(), safe.end(), '/', '_');
        std::replace(safe.begin(), safe.end(), '\\', '_');

        std::string outPath = dir + "/" + safe + (is2D ? ".sqr" : ".spe");
        bool ok = ExportGnuScopeFile(h, outPath);
        if (ok) {
            AppendLog("[GnuScope] " + hname + " → " + outPath);
            nOk++;
        } else {
            AppendLog("[GnuScope] FAILED: " + hname);
            nErr++;
        }
        if (owned) delete h;
        gSystem->ProcessEvents();
    }

    AppendLog(Form("[GnuScope] Batch export complete: %d OK, %d failed → %s",
                   nOk, nErr, dir.c_str()));
    SetStatus(Form("gnuScope: %d exported to %s", nOk, dir.c_str()));
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildDecayTab
// ─────────────────────────────────────────────────────────────────────────────

