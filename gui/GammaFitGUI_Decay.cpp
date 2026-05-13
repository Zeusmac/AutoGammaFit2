#include "GammaFitGUI.h"
#include "GammaFitGUI_shared.h"

#include "TGTextEntry.h"
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFitResult.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TFile.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <sys/stat.h>

void GammaFitGUI::BuildDecayTab(TGCompositeFrame* p)
{
    TGCanvas* sc = new TGCanvas(p, 308, 860, kSunkenFrame);
    p->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    TGCompositeFrame* cf = new TGCompositeFrame(sc->GetViewPort(), 295, 10, kVerticalFrame);
    sc->SetContainer(cf);
    p = cf;

    // ── TH2 selection ─────────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "2D Histogram");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

        decayTh2Combo_ = new TGComboBox(grp, 800);
        decayTh2Combo_->Resize(280, 22);
        grp->AddFrame(decayTh2Combo_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        decayTh2Combo_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                "OnDecayTh2Changed(Int_t)");

        TGHorizontalFrame* axisRow = new TGHorizontalFrame(grp);
        grp->AddFrame(axisRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        axisRow->AddFrame(new TGLabel(axisRow, "Gamma axis:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
        decayGammaAxisCombo_ = new TGComboBox(axisRow, 801);
        decayGammaAxisCombo_->AddEntry("X axis", 1);
        decayGammaAxisCombo_->AddEntry("Y axis", 2);
        decayGammaAxisCombo_->Select(1, kFALSE);
        decayGammaAxisCombo_->Resize(100, 22);
        axisRow->AddFrame(decayGammaAxisCombo_, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
        decayGammaAxisCombo_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                      "OnDecayTh2Changed(Int_t)");
    }

    // ── Peak selection ────────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Fitted Peaks");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        TGHorizontalFrame* row = new TGHorizontalFrame(grp);
        grp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        row->AddFrame(new TGLabel(row, "σ window:"),
                      new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        decaySigRangeEntry_ = new TGNumberEntry(row, 1.0, 5, -1,
                                                TGNumberFormat::kNESRealTwo,
                                                TGNumberFormat::kNEAPositive);
        decaySigRangeEntry_->SetWidth(60);
        row->AddFrame(decaySigRangeEntry_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        row->AddFrame(new TGLabel(row, "σ"), new TGLayoutHints(kLHintsCenterY, 0, 8, 0, 0));

        TGTextButton* refreshBtn = new TGTextButton(row, "Refresh");
        row->AddFrame(refreshBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        refreshBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRefreshDecayPeaks()");

        TGTextButton* previewBtn = new TGTextButton(row, "Preview");
        row->AddFrame(previewBtn, new TGLayoutHints(kLHintsLeft));
        previewBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPreviewDecay()");
        previewBtn->SetToolTipText("Show decay projection for selected peak without fitting");

        // Rebin row
        TGHorizontalFrame* rebinRow = new TGHorizontalFrame(grp);
        grp->AddFrame(rebinRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
        rebinRow->AddFrame(new TGLabel(rebinRow, "Rebin:"),
                           new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        decayRebinEntry_ = new TGNumberEntry(rebinRow, 1, 4, -1,
                                             TGNumberFormat::kNESInteger,
                                             TGNumberFormat::kNEAPositive,
                                             TGNumberFormat::kNELLimitMinMax, 1, 1024);
        decayRebinEntry_->SetWidth(55);
        rebinRow->AddFrame(decayRebinEntry_, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));
        rebinRow->AddFrame(new TGLabel(rebinRow, "(1 = original)"),
                           new TGLayoutHints(kLHintsCenterY));

        decayPeakList_ = new TGListBox(grp, 802);
        decayPeakList_->Resize(280, 100);
        grp->AddFrame(decayPeakList_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        decayPeakList_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                "OnDecayPeakSelected(Int_t)");

        // Label / class for selected peak
        TGHorizontalFrame* lblRow = new TGHorizontalFrame(grp);
        grp->AddFrame(lblRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 0));
        lblRow->AddFrame(new TGLabel(lblRow, "Label:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        decayLabelEntry_ = new TGTextEntry(lblRow, "");
        decayLabelEntry_->SetWidth(120);
        lblRow->AddFrame(decayLabelEntry_, new TGLayoutHints(kLHintsLeft, 0, 8, 0, 0));

        lblRow->AddFrame(new TGLabel(lblRow, "Class:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        decayLabelClassCombo_ = new TGComboBox(lblRow, 805);
        decayLabelClassCombo_->AddEntry("(none)",               1);
        decayLabelClassCombo_->AddEntry("Parent",               2);
        decayLabelClassCombo_->AddEntry("Daughter",             3);
        decayLabelClassCombo_->AddEntry("Granddaughter",        4);
        decayLabelClassCombo_->AddEntry("Beta-n Daughter",      5);
        decayLabelClassCombo_->AddEntry("Beta-2n Daughter",     6);
        decayLabelClassCombo_->AddEntry("Beta-n Granddaughter", 7);
        decayLabelClassCombo_->AddEntry("Beta-2n Granddaughter",8);
        decayLabelClassCombo_->AddEntry("Background",           9);
        decayLabelClassCombo_->AddEntry("X-ray",               10);
        decayLabelClassCombo_->Select(1, kFALSE);
        decayLabelClassCombo_->Resize(130, 22);
        lblRow->AddFrame(decayLabelClassCombo_, new TGLayoutHints(kLHintsLeft));

        TGTextButton* applyLblBtn = new TGTextButton(grp, "Apply Label to Peak");
        grp->AddFrame(applyLblBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        applyLblBtn->Connect("Clicked()", "GammaFitGUI", this, "OnDecayApplyLabel()");
        applyLblBtn->SetToolTipText("Save label and class to the gamma-projection cache entry for this peak");

        TGTextButton* loadDecCacheBtn = new TGTextButton(grp, "Load Decay Cache");
        grp->AddFrame(loadDecCacheBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));
        loadDecCacheBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLoadDecayCache()");
        loadDecCacheBtn->SetToolTipText("Reload saved decay fit results for the current TH2");
    }

    // ── Model ─────────────────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Decay Model");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        TGHorizontalFrame* modelRow = new TGHorizontalFrame(grp);
        grp->AddFrame(modelRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        modelRow->AddFrame(new TGLabel(modelRow, "Model:"),
                           new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
        decayModelCombo_ = new TGComboBox(modelRow, 803);
        decayModelCombo_->AddEntry("Parent + BG",         1);
        decayModelCombo_->AddEntry("Daughter + BG",       2);
        decayModelCombo_->AddEntry("Granddaughter + BG",  3);
        decayModelCombo_->AddEntry("Background only",     4);
        decayModelCombo_->Select(1, kFALSE);
        decayModelCombo_->Resize(160, 22);
        modelRow->AddFrame(decayModelCombo_, new TGLayoutHints(kLHintsLeft));
        decayModelCombo_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                  "OnDecayModelChanged(Int_t)");

        // Equation label — updated by OnDecayModelChanged
        decayEquationLabel_ = new TGLabel(grp,
            "f(t) = A*exp(-ln2*t/T_P) + BG");
        decayEquationLabel_->SetTextJustify(kTextLeft);
        grp->AddFrame(decayEquationLabel_,
                      new TGLayoutHints(kLHintsExpandX, 6, 2, 0, 4));

        // Half-life seeds
        auto AddHalfLifeRow = [&](TGCompositeFrame* par, const char* lbl,
                                  TGNumberEntry*& entry, TGCheckButton*& fix) {
            TGHorizontalFrame* row = new TGHorizontalFrame(par);
            par->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
            row->AddFrame(new TGLabel(row, lbl),
                          new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            entry = new TGNumberEntry(row, 100.0, 10, -1,
                                      TGNumberFormat::kNESRealFour,
                                      TGNumberFormat::kNEAPositive);
            entry->SetWidth(90);
            row->AddFrame(entry, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            fix = new TGCheckButton(row, "Fix");
            row->AddFrame(fix, new TGLayoutHints(kLHintsCenterY));
        };

        AddHalfLifeRow(grp, "T½ Parent:  ", decayThalfP_, decayFixP_);
        AddHalfLifeRow(grp, "T½ Daughter:", decayThalfD_, decayFixD_);
        AddHalfLifeRow(grp, "T½ GDaughter:", decayThalfG_, decayFixG_);

        // Fit range
        {
            TGHorizontalFrame* rr = new TGHorizontalFrame(grp);
            grp->AddFrame(rr, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 0));
            rr->AddFrame(new TGLabel(rr, "Fit from:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            decayFitLo_ = new TGNumberEntry(rr, 0.0, 8, -1,
                                            TGNumberFormat::kNESRealFour,
                                            TGNumberFormat::kNEAAnyNumber);
            decayFitLo_->SetWidth(80);
            rr->AddFrame(decayFitLo_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            rr->AddFrame(new TGLabel(rr, "to:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            decayFitHi_ = new TGNumberEntry(rr, 0.0, 8, -1,
                                            TGNumberFormat::kNESRealFour,
                                            TGNumberFormat::kNEAAnyNumber);
            decayFitHi_->SetWidth(80);
            rr->AddFrame(decayFitHi_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            decayFitFullRange_ = new TGCheckButton(rr, "Full");
            decayFitFullRange_->SetState(kButtonDown);
            decayFitFullRange_->SetToolTipText("Use full time-axis range");
            rr->AddFrame(decayFitFullRange_, new TGLayoutHints(kLHintsCenterY));
        }

        TGTextButton* fitBtn = new TGTextButton(grp, "  Fit Decay  ");
        grp->AddFrame(fitBtn, new TGLayoutHints(kLHintsCenterX, 0, 0, 6, 2));
        fitBtn->Connect("Clicked()", "GammaFitGUI", this, "OnFitDecay()");
    }

    // ── Peak Counts vs Time ───────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Peak Counts vs Time");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        TGHorizontalFrame* row = new TGHorizontalFrame(grp);
        grp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        row->AddFrame(new TGLabel(row, "Slice width:"),
                      new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        decaySliceWidthEntry_ = new TGNumberEntry(row, 10.0, 8, -1,
                                                  TGNumberFormat::kNESRealFour,
                                                  TGNumberFormat::kNEAPositive);
        decaySliceWidthEntry_->SetWidth(80);
        row->AddFrame(decaySliceWidthEntry_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        row->AddFrame(new TGLabel(row, "ms"),
                      new TGLayoutHints(kLHintsCenterY, 0, 0, 0, 0));

        TGTextButton* plotBtn = new TGTextButton(grp, "Plot Peak Counts vs Time");
        grp->AddFrame(plotBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        plotBtn->Connect("Clicked()", "GammaFitGUI", this, "OnMakePeakCountVsTime()");
        plotBtn->SetToolTipText(
            "For each time slice (x>0), project the gamma axis and fit the selected peak; plot counts vs time");
    }

    // ── Results ───────────────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Fit Results");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 4, 4, 2, 4));
        decayResultView_ = new TGTextView(grp, 280, 140);
        grp->AddFrame(decayResultView_, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// PopulateDecayTh2Combo
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::PopulateDecayTh2Combo()
{
    decayTh2Combo_->RemoveAll();
    int idx = 1;
    for (const auto& name : th2Names_) {
        decayTh2Combo_->AddEntry(name.c_str(), idx++);
    }
    decayTh2Combo_->MapSubwindows();
    decayTh2Combo_->Layout();
}

// ─────────────────────────────────────────────────────────────────────────────
// Decay tab slots
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnDecayTh2Changed(Int_t /*id*/)
{
    decayPeakEs_.clear();
    decayPeakSigs_.clear();
    decayPeakKeys_.clear();
    decayFitStore_.clear();
    decayGammaProjName_.clear();
    decayPeakList_->RemoveAll();
    decayPeakList_->MapSubwindows();
    decayPeakList_->Layout();
    if (decayLabelEntry_) decayLabelEntry_->SetText("");
    if (decayLabelClassCombo_) decayLabelClassCombo_->Select(1, kFALSE);
}

void GammaFitGUI::OnDecayModelChanged(Int_t id)
{
    if (!decayEquationLabel_) return;
    const char* eq = "";
    switch (id) {
        case 1:
            eq = "f(t) = A*exp(-ln2*t/T_P) + BG";
            break;
        case 2:
            eq = "f(t) = A*lD/(lD-lP)*(exp(-lP*t)-exp(-lD*t))+BG  [l=ln2/T]";
            break;
        case 3:
            eq = "f(t) = A*lP*lD*lG*Sum[exp(-li*t)/((lj-li)(lk-li))]+BG";
            break;
        case 4:
            eq = "f(t) = BG  (constant background)";
            break;
        default:
            eq = "";
    }
    decayEquationLabel_->SetText(eq);
    decayEquationLabel_->Layout();
}

void GammaFitGUI::OnLoadDecayCache()
{
    if (decayTh2Name_.empty()) {
        AppendLog("[Decay] Run 'Refresh' first to set the active TH2");
        return;
    }
    LoadDecayFitCache();
    // Re-select current peak to restore its values from the freshly loaded cache
    if (decayPeakList_) {
        Int_t sel = decayPeakList_->GetSelected();
        if (sel >= 1) OnDecayPeakSelected(sel);
    }
    SetStatus("Decay cache loaded: " + decayTh2Name_);
}

static int ClassToDecayComboIdx(const std::string& cls)
{
    if (cls == "Parent")               return 2;
    if (cls == "Daughter")             return 3;
    if (cls == "Granddaughter")        return 4;
    if (cls == "Beta-n Daughter"  || cls == "Beta-n")  return 5;
    if (cls == "Beta-2n Daughter" || cls == "Beta-2n") return 6;
    if (cls == "Beta-n Granddaughter")  return 7;
    if (cls == "Beta-2n Granddaughter") return 8;
    if (cls == "Background")           return 9;
    if (cls == "X-ray")                return 10;
    return 1;
}

void GammaFitGUI::OnDecayPeakSelected(Int_t id)
{
    if (id < 1 || (size_t)id > decayPeakKeys_.size()) return;
    if (decayGammaProjName_.empty()) return;

    // Populate label/class from gamma projection cache
    FitDatabase fdb;
    fdb.Load(CacheFileFor(decayGammaProjName_));
    const auto& entries = fdb.GetEntries();
    auto it = entries.find(decayPeakKeys_[id - 1]);
    if (it != entries.end()) {
        if (decayLabelEntry_)
            decayLabelEntry_->SetText(it->second.label.c_str());
        if (decayLabelClassCombo_)
            decayLabelClassCombo_->Select(
                ClassToDecayComboIdx(it->second.classification), kFALSE);
    }

    // Restore stored decay fit half-lives, label and class for this peak
    double E = decayPeakEs_[id - 1];
    auto fit = decayFitStore_.find(E);
    if (fit != decayFitStore_.end()) {
        const DecayFitResult& r = fit->second;
        if (decayModelCombo_) decayModelCombo_->Select(r.model, kFALSE);
        // Half-life parameter indices depend on model
        // model 1: [A, T_P, BG]  → T_P at index 1
        // model 2: [A, T_P, T_D, BG] → T_P=1, T_D=2
        // model 3: [A, T_P, T_D, T_G, BG] → T_P=1, T_D=2, T_G=3
        if (r.params.size() >= 2 && decayThalfP_)
            decayThalfP_->SetNumber(r.params[1]);
        if (r.params.size() >= 3 && r.model >= 2 && decayThalfD_)
            decayThalfD_->SetNumber(r.params[2]);
        if (r.params.size() >= 4 && r.model >= 3 && decayThalfG_)
            decayThalfG_->SetNumber(r.params[3]);
        // Restore label/class from decay cache (overrides gamma proj cache)
        if (!r.label.empty() && decayLabelEntry_)
            decayLabelEntry_->SetText(r.label.c_str());
        if (!r.classification.empty() && decayLabelClassCombo_)
            decayLabelClassCombo_->Select(ClassToDecayComboIdx(r.classification), kFALSE);
        // Restore sigma window
        if (r.Nsig > 0 && decaySigRangeEntry_)
            decaySigRangeEntry_->SetNumber(r.Nsig);
    }
}

void GammaFitGUI::OnRefreshDecayPeaks()
{
    decayPeakEs_.clear();
    decayPeakSigs_.clear();
    decayPeakKeys_.clear();
    decayPeakList_->RemoveAll();

    TGLBEntry* selE = decayTh2Combo_->GetSelectedEntry();
    if (!selE) { AppendLog("[Decay] No TH2 selected"); return; }
    decayTh2Name_ = selE->GetTitle();

    int axisId = decayGammaAxisCombo_->GetSelected();
    decayGammaProjName_ = decayTh2Name_ + (axisId == 1 ? "_px" : "_py");

    FitDatabase fdb;
    fdb.Load(CacheFileFor(decayGammaProjName_));

    int listIdx = 1;
    for (const auto& kv : fdb.GetEntries()) {
        const FitEntry& e = kv.second;
        if (e.params.size() < 5) continue;
        int nPeaks = ((int)e.params.size() - 2) / 3;
        for (int i = 0; i < nPeaks; i++) {
            double Ep  = e.params[3 * i + 1];
            double sig = std::abs(e.params[3 * i + 2]);
            if (sig <= 0 || sig > 100) continue;
            decayPeakEs_.push_back(Ep);
            decayPeakSigs_.push_back(sig);
            decayPeakKeys_.push_back(kv.first);
            std::string lbl = Form("%.2f keV  (σ=%.3f)", Ep, sig);
            if (!e.label.empty()) lbl = e.label + "  " + lbl;
            if (!e.classification.empty()) lbl += "  [" + e.classification + "]";
            decayPeakList_->AddEntry(lbl.c_str(), listIdx++);
        }
    }
    decayPeakList_->MapSubwindows();
    decayPeakList_->Layout();
    LoadDecayFitCache();
    AppendLog(Form("[Decay] %d peaks loaded from cache of %s",
                   (int)decayPeakEs_.size(), decayGammaProjName_.c_str()));
}

void GammaFitGUI::OnFitDecay()
{
    if (!inputFile_) { AppendLog("[Decay] No ROOT file loaded"); return; }
    if (decayTh2Name_.empty()) { AppendLog("[Decay] Run 'Refresh' first"); return; }

    Int_t peakSel = decayPeakList_->GetSelected();
    if (peakSel < 1 || (size_t)peakSel > decayPeakEs_.size()) {
        AppendLog("[Decay] Select a peak first"); return;
    }
    double E    = decayPeakEs_[peakSel - 1];
    double sig  = decayPeakSigs_[peakSel - 1];
    double Nsig = decaySigRangeEntry_->GetNumber();

    TH2* h2 = (TH2*)inputFile_->Get(decayTh2Name_.c_str());
    if (!h2) { AppendLog("[Decay] TH2 not found in file"); return; }

    int axisId = decayGammaAxisCombo_->GetSelected();
    TH1* hDecay = nullptr;
    double eMin = E - Nsig * sig;
    double eMax = E + Nsig * sig;

    if (axisId == 1) {
        // Gamma on X → project onto Y axis
        int b1 = h2->GetXaxis()->FindBin(eMin);
        int b2 = h2->GetXaxis()->FindBin(eMax);
        hDecay = h2->ProjectionY(Form("hDecay_%.1f", E), b1, b2);
    } else {
        // Gamma on Y → project onto X axis
        int b1 = h2->GetYaxis()->FindBin(eMin);
        int b2 = h2->GetYaxis()->FindBin(eMax);
        hDecay = h2->ProjectionX(Form("hDecay_%.1f", E), b1, b2);
    }

    if (!hDecay || hDecay->GetEntries() == 0) {
        AppendLog("[Decay] Projection is empty — check TH2 and peak selection");
        return;
    }
    hDecay->SetDirectory(nullptr);

    // Rebin
    int rebin = decayRebinEntry_ ? (int)decayRebinEntry_->GetNumber() : 1;
    if (rebin > 1) hDecay->Rebin(rebin);

    // Restrict to t > 0
    double tAxisMin = std::max(0.0, hDecay->GetXaxis()->GetXmin());
    hDecay->GetXaxis()->SetRangeUser(tAxisMin, hDecay->GetXaxis()->GetXmax());

    TCanvas* c = canvas_->GetCanvas();
    c->Clear();
    c->cd();
    hDecay->GetXaxis()->SetTitle("Time (ms)");
    hDecay->GetYaxis()->SetTitle("Counts");
    hDecay->SetLineColor(kBlack);
    hDecay->SetMarkerSize(0);
    hDecay->Draw("hist E");
    c->Modified(); c->Update();
    gSystem->ProcessEvents();

    // ── Build fit function ────────────────────────────────────────────────────
    int modelId = decayModelCombo_->GetSelected();
    bool fullRange = !decayFitFullRange_ || decayFitFullRange_->IsOn();
    double xlo = fullRange ? tAxisMin : std::max(0.0, decayFitLo_->GetNumber());
    double xhi = fullRange ? hDecay->GetXaxis()->GetXmax() : decayFitHi_->GetNumber();
    if (xlo >= xhi) { xlo = tAxisMin; xhi = hDecay->GetXaxis()->GetXmax(); }

    TF1* fDecay = nullptr;

    if (modelId == 4) {
        fDecay = new TF1("fDecayBG", "[0]", xlo, xhi);
        fDecay->SetParName(0, "BG");
        fDecay->SetParameter(0, std::max(0.0, hDecay->GetMinimum()));

    } else if (modelId == 1) {
        fDecay = new TF1("fDecayP",
                         "[0]*exp(-0.6931471805599453*x/[1]) + [2]",
                         xlo, xhi);
        fDecay->SetParName(0, "A");
        fDecay->SetParName(1, "T_{1/2}");
        fDecay->SetParName(2, "BG");
        double T = decayThalfP_->GetNumber();
        fDecay->SetParameters(std::max(1.0, hDecay->GetMaximum() - hDecay->GetMinimum()),
                              T > 0 ? T : (xhi - xlo) / 2.0,
                              std::max(0.0, hDecay->GetMinimum()));
        fDecay->SetParLimits(1, 1e-6, 1e12);
        if (decayFixP_->IsOn()) fDecay->FixParameter(1, T);

    } else if (modelId == 2) {
        auto lambdaD = [](double* x, double* p) -> double {
            double t  = x[0];
            double A  = p[0];
            double T1 = p[1], T2 = p[2];
            double BG = p[3];
            if (T1 <= 0 || T2 <= 0) return BG;
            const double ln2 = 0.6931471805599453;
            double lP = ln2 / T1, lD = ln2 / T2;
            double denom = lD - lP;
            if (std::abs(denom) < 1e-10 * std::max(lP, lD))
                return A * lP * t * TMath::Exp(-lP * t) + BG;
            return A * lD / denom * (TMath::Exp(-lP * t) - TMath::Exp(-lD * t)) + BG;
        };
        fDecay = new TF1("fDecayD", lambdaD, xlo, xhi, 4);
        fDecay->SetParName(0, "A");
        fDecay->SetParName(1, "T_{1/2}^{P}");
        fDecay->SetParName(2, "T_{1/2}^{D}");
        fDecay->SetParName(3, "BG");
        double TP = decayThalfP_->GetNumber();
        double TD = decayThalfD_->GetNumber();
        fDecay->SetParameters(hDecay->GetMaximum(),
                              TP > 0 ? TP : (xhi - xlo) / 4.0,
                              TD > 0 ? TD : (xhi - xlo) / 2.0,
                              std::max(0.0, hDecay->GetMinimum()));
        fDecay->SetParLimits(1, 1e-6, 1e12);
        fDecay->SetParLimits(2, 1e-6, 1e12);
        if (decayFixP_->IsOn()) fDecay->FixParameter(1, TP);
        if (decayFixD_->IsOn()) fDecay->FixParameter(2, TD);

    } else {
        // Granddaughter
        auto lambdaG = [](double* x, double* p) -> double {
            double t  = x[0];
            double A  = p[0];
            double T1 = p[1], T2 = p[2], T3 = p[3];
            double BG = p[4];
            if (T1 <= 0 || T2 <= 0 || T3 <= 0) return BG;
            const double ln2 = 0.6931471805599453;
            double lP = ln2/T1, lD = ln2/T2, lG = ln2/T3;
            double eps = 1e-10 * std::max({lP, lD, lG});
            if (std::abs(lD-lP) < eps || std::abs(lG-lP) < eps || std::abs(lG-lD) < eps)
                return A * TMath::Exp(-lP * t) + BG;
            double t1 = TMath::Exp(-lP*t) / ((lD-lP) * (lG-lP));
            double t2 = TMath::Exp(-lD*t) / ((lP-lD) * (lG-lD));
            double t3 = TMath::Exp(-lG*t) / ((lP-lG) * (lD-lG));
            return A * lP * lD * lG * (t1 + t2 + t3) + BG;
        };
        fDecay = new TF1("fDecayG", lambdaG, xlo, xhi, 5);
        fDecay->SetParName(0, "A");
        fDecay->SetParName(1, "T_{1/2}^{P}");
        fDecay->SetParName(2, "T_{1/2}^{D}");
        fDecay->SetParName(3, "T_{1/2}^{G}");
        fDecay->SetParName(4, "BG");
        double TP = decayThalfP_->GetNumber();
        double TD = decayThalfD_->GetNumber();
        double TG = decayThalfG_->GetNumber();
        fDecay->SetParameters(hDecay->GetMaximum(),
                              TP > 0 ? TP : (xhi-xlo)/8.0,
                              TD > 0 ? TD : (xhi-xlo)/4.0,
                              TG > 0 ? TG : (xhi-xlo)/2.0,
                              std::max(0.0, hDecay->GetMinimum()));
        fDecay->SetParLimits(1, 1e-6, 1e12);
        fDecay->SetParLimits(2, 1e-6, 1e12);
        fDecay->SetParLimits(3, 1e-6, 1e12);
        if (decayFixP_->IsOn()) fDecay->FixParameter(1, TP);
        if (decayFixD_->IsOn()) fDecay->FixParameter(2, TD);
        if (decayFixG_->IsOn()) fDecay->FixParameter(3, TG);
    }

    fDecay->SetLineColor(kRed);
    fDecay->SetLineWidth(2);

    TFitResultPtr fitRes = hDecay->Fit(fDecay, "SR");

    hDecay->Draw("hist E");
    fDecay->Draw("same");
    c->Modified(); c->Update();

    // ── Report ────────────────────────────────────────────────────────────────
    decayResultView_->Clear();
    decayResultView_->AddLine(Form("Peak: %.2f keV  window: ±%.1fσ  [%.2f, %.2f]",
                                   E, Nsig, eMin, eMax));
    double chi2ndf = (fitRes.Get() && fitRes->Ndf() > 0)
                     ? fitRes->Chi2() / fitRes->Ndf() : -1.0;
    decayResultView_->AddLine(Form("Chi2/NDF: %.4g  status: %d",
                                   chi2ndf, fitRes.Get() ? fitRes->Status() : -1));
    for (int i = 0; i < fDecay->GetNpar(); i++) {
        decayResultView_->AddLine(Form("  %-16s = %11.5g  ±  %.4g",
            fDecay->GetParName(i), fDecay->GetParameter(i), fDecay->GetParError(i)));
    }
    // Store fit result in memory + save to disk
    {
        DecayFitResult r;
        r.peakE          = E;
        r.model          = modelId;
        r.chi2ndf        = chi2ndf;
        r.status         = fitRes.Get() ? fitRes->Status() : -1;
        r.histName       = decayTh2Name_;
        r.eMin           = eMin;
        r.eMax           = eMax;
        r.Nsig           = Nsig;
        // Preserve existing label/class if already set
        auto existing = decayFitStore_.find(E);
        if (existing != decayFitStore_.end()) {
            r.label          = existing->second.label;
            r.classification = existing->second.classification;
        }
        // Overwrite with whatever is shown in the decay label/class widgets
        if (decayLabelEntry_) {
            std::string lbl = decayLabelEntry_->GetText();
            if (!lbl.empty()) r.label = lbl;
        }
        if (decayLabelClassCombo_) {
            static const char* kClsNames[] = {
                "", "Parent", "Daughter", "Granddaughter",
                "Beta-n Daughter", "Beta-2n Daughter",
                "Beta-n Granddaughter", "Beta-2n Granddaughter",
                "Background", "X-ray"
            };
            int ci = decayLabelClassCombo_->GetSelected() - 1;
            if (ci >= 1 && ci <= 9) r.classification = kClsNames[ci];
        }
        int np = fDecay->GetNpar();
        r.params.resize(np); r.errors.resize(np);
        for (int i = 0; i < np; i++) {
            r.params[i] = fDecay->GetParameter(i);
            r.errors[i] = fDecay->GetParError(i);
        }
        decayFitStore_[E] = r;
        SaveDecayFitCache();
    }

    // Auto-classify and save to gamma projection cache
    if (fitRes.Get() && fitRes->Status() == 0) {
        static const char* kModelCls[] = { "", "Parent", "Daughter", "Granddaughter", "Background" };
        std::string suggestedCls = (modelId >= 1 && modelId <= 4) ? kModelCls[modelId] : "";

        if (!suggestedCls.empty()) {
            // Update class combo
            if (decayLabelClassCombo_)
                decayLabelClassCombo_->Select(ClassToDecayComboIdx(suggestedCls), kFALSE);

            // Save to gamma projection cache automatically
            if (!decayGammaProjName_.empty() && peakSel >= 1 &&
                (size_t)peakSel <= decayPeakKeys_.size()) {
                FitDatabase fdb;
                fdb.Load(CacheFileFor(decayGammaProjName_));
                const auto& ents = fdb.GetEntries();
                auto eit = ents.find(decayPeakKeys_[peakSel - 1]);
                if (eit != ents.end()) {
                    FitEntry fe = eit->second;
                    if (fe.classification.empty())
                        fe.classification = suggestedCls;
                    fdb.ForceStore(decayPeakKeys_[peakSel - 1], fe);
                    mkdir(kCacheDir, 0755);
                    fdb.Save(CacheFileFor(decayGammaProjName_));
                }
            }
            decayResultView_->AddLine(
                Form("  → Auto-classified as: %s", suggestedCls.c_str()));
        }
    }
    decayResultView_->ShowBottom();

    AppendLog(Form("[Decay] Fit done: peak=%.2f keV  chi2/ndf=%.3f", E, chi2ndf));
    SetStatus(Form("Decay fit: %.2f keV  chi2/ndf=%.3f", E, chi2ndf));
}

// ─────────────────────────────────────────────────────────────────────────────
// Decay fit cache helpers
// ─────────────────────────────────────────────────────────────────────────────

std::string GammaFitGUI::DecayCacheFileFor() const
{
    if (decayTh2Name_.empty()) return "";
    return std::string(kCacheDir) + "/decay_fits_" + decayTh2Name_ + ".dat";
}

void GammaFitGUI::SaveDecayFitCache()
{
    std::string path = DecayCacheFileFor();
    if (path.empty()) return;
    mkdir(kCacheDir, 0755);
    std::ofstream out(path);
    if (!out.is_open()) return;
    out << std::fixed << std::setprecision(8);
    for (const auto& kv : decayFitStore_) {
        const DecayFitResult& r = kv.second;
        int np = (int)r.params.size();
        // format: peakE model npar [pi ei]... chi2ndf status eMin eMax Nsig histName label classification
        out << r.peakE << " " << r.model << " " << np;
        for (int i = 0; i < np; i++) out << " " << r.params[i] << " " << r.errors[i];
        out << " " << r.chi2ndf << " " << r.status;
        out << " " << r.eMin << " " << r.eMax << " " << r.Nsig;
        // string fields: replace spaces with underscores, use "-" for empty
        auto esc = [](const std::string& s) -> std::string {
            if (s.empty()) return "-";
            std::string t = s;
            for (char& c : t) if (c == ' ') c = '_';
            return t;
        };
        out << " " << esc(r.histName)
            << " " << esc(r.label)
            << " " << esc(r.classification)
            << "\n";
    }
}

void GammaFitGUI::LoadDecayFitCache()
{
    decayFitStore_.clear();
    std::string path = DecayCacheFileFor();
    if (path.empty()) return;
    std::ifstream in(path);
    if (!in.is_open()) return;
    auto unescape = [](const std::string& s) -> std::string {
        if (s == "-") return "";
        std::string t = s;
        for (char& c : t) if (c == '_') c = ' ';
        return t;
    };
    std::string line;
    int count = 0;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        DecayFitResult r;
        int np = 0;
        if (!(ss >> r.peakE >> r.model >> np)) continue;
        r.params.resize(np); r.errors.resize(np);
        bool ok = true;
        for (int i = 0; i < np; i++) {
            if (!(ss >> r.params[i] >> r.errors[i])) { ok = false; break; }
        }
        if (!ok) continue;
        ss >> r.chi2ndf >> r.status;
        // Optional new fields (backward compatible with old cache files)
        ss >> r.eMin >> r.eMax >> r.Nsig;
        std::string hn, lb, cl;
        if (ss >> hn) r.histName       = unescape(hn);
        if (ss >> lb) r.label          = unescape(lb);
        if (ss >> cl) r.classification = unescape(cl);
        decayFitStore_[r.peakE] = r;
        ++count;
    }
    AppendLog(Form("[Decay] Loaded %d cached decay fit(s) from %s", count, path.c_str()));
}

// ─────────────────────────────────────────────────────────────────────────────
// OnPreviewDecay — show decay projection without fitting
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnPreviewDecay()
{
    if (!inputFile_) { AppendLog("[Decay] No ROOT file loaded"); return; }
    if (decayTh2Name_.empty()) { AppendLog("[Decay] Run 'Refresh' first"); return; }

    Int_t peakSel = decayPeakList_->GetSelected();
    if (peakSel < 1 || (size_t)peakSel > decayPeakEs_.size()) {
        AppendLog("[Decay] Select a peak first"); return;
    }
    double E    = decayPeakEs_[peakSel - 1];
    double sig  = decayPeakSigs_[peakSel - 1];
    double Nsig = decaySigRangeEntry_->GetNumber();
    double eMin = E - Nsig * sig;
    double eMax = E + Nsig * sig;

    TH2* h2 = (TH2*)inputFile_->Get(decayTh2Name_.c_str());
    if (!h2) { AppendLog("[Decay] TH2 not found"); return; }

    int axisId = decayGammaAxisCombo_->GetSelected();
    TH1* hDecay = nullptr;
    if (axisId == 1) {
        int b1 = h2->GetXaxis()->FindBin(eMin);
        int b2 = h2->GetXaxis()->FindBin(eMax);
        hDecay = h2->ProjectionY(Form("hPrev_%.1f", E), b1, b2);
    } else {
        int b1 = h2->GetYaxis()->FindBin(eMin);
        int b2 = h2->GetYaxis()->FindBin(eMax);
        hDecay = h2->ProjectionX(Form("hPrev_%.1f", E), b1, b2);
    }
    if (!hDecay) { AppendLog("[Decay] Projection failed"); return; }
    hDecay->SetDirectory(nullptr);

    // Rebin before display
    int rebin = decayRebinEntry_ ? (int)decayRebinEntry_->GetNumber() : 1;
    if (rebin > 1) hDecay->Rebin(rebin);

    // Restrict to t > 0
    double tLo = std::max(0.0, hDecay->GetXaxis()->GetXmin());
    hDecay->GetXaxis()->SetRangeUser(tLo, hDecay->GetXaxis()->GetXmax());

    hDecay->GetXaxis()->SetTitle("Time (ms)");
    hDecay->GetYaxis()->SetTitle("Counts");
    hDecay->SetLineColor(kBlack);
    hDecay->SetMarkerSize(0);

    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();
    hDecay->Draw("hist E");
    c->Modified(); c->Update();
    gSystem->ProcessEvents();

    SetStatus(Form("Preview: %.2f keV  ±%.1fσ  [%.2f, %.2f]", E, Nsig, eMin, eMax));
}

// ─────────────────────────────────────────────────────────────────────────────
// OnDecayApplyLabel — save label/class from Decay tab fields to gamma cache
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnDecayApplyLabel()
{
    Int_t peakSel = decayPeakList_->GetSelected();
    if (peakSel < 1 || (size_t)peakSel > decayPeakKeys_.size()) {
        AppendLog("[Decay] Select a peak first"); return;
    }
    if (decayGammaProjName_.empty()) { AppendLog("[Decay] Refresh first"); return; }

    FitDatabase fdb;
    fdb.Load(CacheFileFor(decayGammaProjName_));
    const auto& entries = fdb.GetEntries();
    const std::string& key = decayPeakKeys_[peakSel - 1];
    auto it = entries.find(key);
    if (it == entries.end()) { AppendLog("[Decay] Cache entry not found"); return; }

    FitEntry e = it->second;
    if (decayLabelEntry_) {
        std::string lbl = decayLabelEntry_->GetText();
        if (!lbl.empty()) e.label = lbl;
    }
    if (decayLabelClassCombo_) {
        int sel = decayLabelClassCombo_->GetSelected();
        // Map this combo's indices to classification strings
        static const char* kClsNames[] = {
            "", "Parent", "Daughter", "Granddaughter",
            "Beta-n Daughter", "Beta-2n Daughter",
            "Beta-n Granddaughter", "Beta-2n Granddaughter",
            "Background", "X-ray"
        };
        if (sel >= 2 && sel <= 10)
            e.classification = kClsNames[sel - 1];
    }
    fdb.ForceStore(key, e);
    mkdir(kCacheDir, 0755);
    fdb.Save(CacheFileFor(decayGammaProjName_));
    AppendLog("[Decay] Saved label=" + e.label + "  class=" + e.classification +
              "  → " + key);
    OnRefreshDecayPeaks();
}

// ─────────────────────────────────────────────────────────────────────────────
// OnMakePeakCountVsTime
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnMakePeakCountVsTime()
{
    if (!inputFile_) { AppendLog("[Decay] No ROOT file loaded"); return; }
    if (decayTh2Name_.empty()) { AppendLog("[Decay] Run 'Refresh' first"); return; }

    Int_t peakSel = decayPeakList_->GetSelected();
    if (peakSel < 1 || (size_t)peakSel > decayPeakEs_.size()) {
        AppendLog("[Decay] Select a peak first"); return;
    }

    double E      = decayPeakEs_[peakSel - 1];
    double sig    = decayPeakSigs_[peakSel - 1];
    double Nsig   = decaySigRangeEntry_->GetNumber();
    double sliceW = decaySliceWidthEntry_->GetNumber();
    if (sliceW <= 0) { AppendLog("[Decay] Slice width must be > 0"); return; }

    TH2* h2 = (TH2*)inputFile_->Get(decayTh2Name_.c_str());
    if (!h2) { AppendLog("[Decay] TH2 not found"); return; }

    // Identify time axis
    int axisId = decayGammaAxisCombo_->GetSelected();
    TAxis* timeAxis = (axisId == 1) ? h2->GetYaxis() : h2->GetXaxis();

    double tMin = std::max(0.0, timeAxis->GetXmin());
    double tMax = timeAxis->GetXmax();
    if (tMax <= tMin) { AppendLog("[Decay] Time axis has no positive range"); return; }

    // Gamma energy window
    double eMin = E - Nsig * sig;
    double eMax = E + Nsig * sig;
    const double kSqrt2Pi = 2.5066282746310002;

    // Open output ROOT file for slice histograms (recreate each run)
    mkdir(kCacheDir, 0755);
    std::string sliceRootPath = std::string(kCacheDir) + "/decay_slices_" +
                                decayTh2Name_ + Form("_%.2fkeV.root", E);
    TFile* sliceFile = TFile::Open(sliceRootPath.c_str(), "RECREATE");

    std::vector<double> tCenter, counts, countsErr;

    for (double t = tMin; t < tMax; t += sliceW) {
        double t1 = t, t2 = std::min(t + sliceW, tMax);
        int bt1 = timeAxis->FindBin(t1);
        int bt2 = timeAxis->FindBin(t2 - 1e-9 * sliceW);

        TH1* hSlice = nullptr;
        if (axisId == 1) {
            hSlice = h2->ProjectionX(Form("slice_t%.3f", t), bt1, bt2);
        } else {
            hSlice = h2->ProjectionY(Form("slice_t%.3f", t), bt1, bt2);
        }
        if (!hSlice || hSlice->GetEntries() == 0) { delete hSlice; continue; }
        hSlice->SetDirectory(nullptr);
        hSlice->SetTitle(Form("t=[%.2f,%.2f] ms  E=%.2f keV", t1, t2, E));

        double peakH = hSlice->GetBinContent(hSlice->FindBin(E));
        if (peakH <= 0) { delete hSlice; continue; }

        TF1 fPeak("fPk", "gaus(0)+pol1(3)", eMin, eMax);
        fPeak.SetParameter(0, peakH);
        fPeak.SetParameter(1, E);
        fPeak.SetParameter(2, sig);
        fPeak.SetParLimits(0, 0.001 * peakH, 100.0 * peakH);
        fPeak.SetParLimits(1, E - 5*sig, E + 5*sig);
        fPeak.SetParLimits(2, 0.3*sig, 5.0*sig);

        // Use log-likelihood for low-statistics slices (Poisson regime)
        double windowCounts = hSlice->Integral(hSlice->FindBin(eMin),
                                               hSlice->FindBin(eMax));
        const char* fitOpts = (windowCounts < 20.0) ? "R S Q B L" : "R S Q B";
        TFitResultPtr fr = hSlice->Fit(&fPeak, fitOpts);

        // Save slice histogram (TF1 is attached to hSlice by Fit())
        if (sliceFile && sliceFile->IsOpen()) {
            sliceFile->cd();
            hSlice->Write();
        }
        delete hSlice;

        if (!fr.Get() || fr->Status() != 0) continue;
        double A    = fPeak.GetParameter(0);
        double sFit = fPeak.GetParameter(2);
        if (A <= 0 || sFit <= 0) continue;

        double cnt    = A * sFit * kSqrt2Pi;
        double cntErr = (cnt > 0) ? std::sqrt(cnt) : 1.0;

        tCenter.push_back(0.5 * (t1 + t2));
        counts.push_back(cnt);
        countsErr.push_back(cntErr);
    }

    if (tCenter.empty()) {
        if (sliceFile) { sliceFile->Close(); delete sliceFile; }
        AppendLog("[Decay] No successful slice fits — try wider window or larger slice width");
        return;
    }

    // Plot
    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();

    int N = (int)tCenter.size();
    TGraphErrors* gr = new TGraphErrors(N);
    gr->SetName(Form("peakCounts_%.2fkeV", E));
    gr->SetTitle(Form("Peak counts vs time  (%.2f keV  ±%.1f#sigma);Time (ms);Peak counts", E, Nsig));
    for (int i = 0; i < N; i++) {
        gr->SetPoint(i, tCenter[i], counts[i]);
        gr->SetPointError(i, 0.5*sliceW, countsErr[i]);
    }
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(0.8);
    gr->SetLineColor(kBlue+1);
    gr->SetMarkerColor(kBlue+1);
    gr->Draw("AP");
    c->Modified(); c->Update();
    gSystem->ProcessEvents();

    if (sliceFile && sliceFile->IsOpen()) {
        sliceFile->cd();
        gr->Write("peakCountsVsTime");
        sliceFile->Close();
        delete sliceFile;
        sliceFile = nullptr;
    }

    AppendLog(Form("[Decay] Peak counts plot: %.2f keV  %d points  slice=%.2f ms  saved to %s",
                   E, N, sliceW, sliceRootPath.c_str()));
    SetStatus(Form("Peak counts vs time: %.2f keV  %d slices", E, N));
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildIsotopesTab
// ─────────────────────────────────────────────────────────────────────────────

