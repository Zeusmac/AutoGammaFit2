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
#include "TLine.h"
#include "TLegend.h"

#include "TKey.h"
#include "TStyle.h"

#include <fstream>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <sys/stat.h>

// ─────────────────────────────────────────────────────────────────────────────
// File-scope decay model helpers
// ─────────────────────────────────────────────────────────────────────────────

// Map model ID -> signal type (1=Parent, 2=Daughter, 3=Granddaughter, 4=BGonly)
static int DecaySignalType(int modelId) {
    if (modelId == 4) return 4;
    if (modelId == 1) return 1;
    if (modelId == 2 || modelId == 5 || modelId == 6) return 2;
    if (modelId == 3 || modelId == 7 || modelId == 8) return 3;
    return 1;
}
// Number of signal parameters (not counting BG)
static int DecaySigNpar(int sigType) {
    switch (sigType) {
        case 1: return 2;  // A, T_P
        case 2: return 3;  // A, T_P, T_D
        case 3: return 4;  // A, T_P, T_D, T_G
        default: return 0; // BG only
    }
}
// Number of BG parameters
static int DecayBGNpar(int bgType) {
    if (bgType == 2) return 3;  // BG_flat, A_bg, T_bg
    if (bgType == 3) return 2;  // A_bg, T_bg
    return 1;                   // flat: BG
}

// Forward declaration (defined just before OnFitDecay)
static void ApplyExtraCuts(TH2* h2, int axisId, TH1* hBase,
                            const std::vector<std::pair<double,double>>& cuts);

// Factory: build a TF1 for any (modelId, bgType) combination.
// Parameter layout: [signal params...] [BG params...]
//   signal params: sig==1 -> A, T_P
//                  sig==2 -> A, T_P, T_D
//                  sig==3 -> A, T_P, T_D, T_G
//                  sig==4 -> (none)
//   BG params:     bgType==1 -> BG
//                  bgType==2 -> BG_flat, A_bg, T_bg
//                  bgType==3 -> A_bg, T_bg
static TF1* BuildDecayTF1(const char* name, int modelId, int bgType,
                           double xlo, double xhi)
{
    const int  sig  = DecaySignalType(modelId);
    const int  nSig = DecaySigNpar(sig);
    const int  nBG  = DecayBGNpar(bgType);
    const int  npar = nSig + nBG;
    const double ln2 = 0.6931471805599453;

    TF1* f = new TF1(name,
        [sig, nSig, bgType, ln2](double* x, double* p) -> double {
            double t   = x[0];
            double val = 0.0;
            // Signal
            if (sig == 1 && p[1] > 0) {
                val = p[0] * TMath::Exp(-ln2 * t / p[1]);
            } else if (sig == 2 && p[1] > 0 && p[2] > 0) {
                double lP = ln2/p[1], lD = ln2/p[2];
                double den = lD - lP;
                val = (std::abs(den) < 1e-10 * std::max(lP,lD))
                    ? p[0] * lP * t * TMath::Exp(-lP*t)
                    : p[0] * lD / den * (TMath::Exp(-lP*t) - TMath::Exp(-lD*t));
            } else if (sig == 3 && p[1] > 0 && p[2] > 0 && p[3] > 0) {
                double lP = ln2/p[1], lD = ln2/p[2], lG = ln2/p[3];
                double eps = 1e-10 * std::max({lP, lD, lG});
                if (std::abs(lD-lP)<eps || std::abs(lG-lP)<eps || std::abs(lG-lD)<eps)
                    val = p[0] * TMath::Exp(-lP*t);
                else {
                    double t1 = TMath::Exp(-lP*t)/((lD-lP)*(lG-lP));
                    double t2 = TMath::Exp(-lD*t)/((lP-lD)*(lG-lD));
                    double t3 = TMath::Exp(-lG*t)/((lP-lG)*(lD-lG));
                    val = p[0] * lP * lD * lG * (t1 + t2 + t3);
                }
            }
            // Background
            if (bgType == 1) {
                val += p[nSig];
            } else if (bgType == 2) {
                double Tbg = p[nSig+2];
                val += p[nSig] + (Tbg > 0 ? p[nSig+1] * TMath::Exp(-ln2*t/Tbg) : 0.0);
            } else {
                double Tbg = p[nSig+1];
                val += (Tbg > 0 ? p[nSig] * TMath::Exp(-ln2*t/Tbg) : 0.0);
            }
            return val;
        },
        xlo, xhi, npar);

    // Par names
    static const char* kSigNames[][4] = {
        {"A", "T_{1/2}^{P}", "", ""},        // sig==1
        {"A", "T_{1/2}^{P}", "T_{1/2}^{D}", ""},    // sig==2
        {"A", "T_{1/2}^{P}", "T_{1/2}^{D}", "T_{1/2}^{G}"},  // sig==3
    };
    if (sig >= 1 && sig <= 3)
        for (int i = 0; i < nSig; i++)
            f->SetParName(i, kSigNames[sig-1][i]);

    if (bgType == 1) {
        f->SetParName(nSig, "BG");
    } else if (bgType == 2) {
        f->SetParName(nSig,   "BG_{flat}");
        f->SetParName(nSig+1, "A_{BG}");
        f->SetParName(nSig+2, "T_{BG}");
    } else {
        f->SetParName(nSig,   "A_{BG}");
        f->SetParName(nSig+1, "T_{BG}");
    }
    return f;
}


// Shared helper: build T1/2 row with saved label pointer
static void AddHalfLifeRow(TGCompositeFrame* par, const char* lbl,
                            TGLabel*& lblPtr, TGNumberEntry*& entry,
                            TGCheckButton*& fix, double defVal = 100.0)
{
    TGHorizontalFrame* row = new TGHorizontalFrame(par);
    par->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
    lblPtr = new TGLabel(row, lbl);
    lblPtr->SetWidth(95);
    row->AddFrame(lblPtr, new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
    entry = new TGNumberEntry(row, defVal, 10, -1,
                              TGNumberFormat::kNESRealFour,
                              TGNumberFormat::kNEAPositive);
    entry->SetWidth(90);
    row->AddFrame(entry, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
    fix = new TGCheckButton(row, "Fix");
    row->AddFrame(fix, new TGLayoutHints(kLHintsCenterY));
}

// Add a compact Lo/Hi bounds row indented under the matching T1/2 row
static void AddBoundsRow(TGCompositeFrame* par, TGNumberEntry*& lo, TGNumberEntry*& hi,
                         double defLo = 1e-6, double defHi = 1e12)
{
    TGHorizontalFrame* row = new TGHorizontalFrame(par);
    par->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 18, 2, 0, 2));
    row->AddFrame(new TGLabel(row, "Lo:"), new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
    lo = new TGNumberEntry(row, defLo, 8, -1,
                           TGNumberFormat::kNESRealFour, TGNumberFormat::kNEAAnyNumber);
    lo->SetWidth(72);
    row->AddFrame(lo, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));
    row->AddFrame(new TGLabel(row, "Hi:"), new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
    hi = new TGNumberEntry(row, defHi, 8, -1,
                           TGNumberFormat::kNESRealFour, TGNumberFormat::kNEAAnyNumber);
    hi->SetWidth(72);
    row->AddFrame(hi, new TGLayoutHints(kLHintsLeft));
}

void GammaFitGUI::BuildDecayTab(TGCompositeFrame* p)
{
    decaySubTabs_ = new TGTab(p, 305, 860);
    p->AddFrame(decaySubTabs_, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    // ═══════════════════════════════════════════════════════════════════════════
    // Sub-tab 1: Cuts  -  TH2 selection, asymmetric peak cut, peak list, rebin,
    //            peak counts vs time, decay model + T1/2 fitter, fit results
    // ═══════════════════════════════════════════════════════════════════════════
    {
        TGCompositeFrame* tf = decaySubTabs_->AddTab("Cuts");
        TGCanvas* sc = new TGCanvas(tf, 305, 820, kSunkenFrame);
        tf->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        TGCompositeFrame* cf = new TGCompositeFrame(sc->GetViewPort(), 293, 10, kVerticalFrame);
        sc->SetContainer(cf);
        TGCompositeFrame* p2 = cf;

        // ── TH2 selection ─────────────────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p2, "2D Histogram");
            p2->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

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
            axisRow->AddFrame(decayGammaAxisCombo_, new TGLayoutHints(kLHintsLeft));
            decayGammaAxisCombo_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                          "OnDecayTh2Changed(Int_t)");
        }

        // ── Fitted Peaks ──────────────────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p2, "Fitted Peaks");
            p2->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

            TGHorizontalFrame* cacheRow = new TGHorizontalFrame(grp);
            grp->AddFrame(cacheRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 2));
            cacheRow->AddFrame(new TGLabel(cacheRow, "Peaks cache:"),
                               new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            decayCacheCombo_ = new TGComboBox(cacheRow, 804);
            decayCacheCombo_->Resize(150, 22);
            cacheRow->AddFrame(decayCacheCombo_,
                               new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 0, 4, 0, 0));
            TGTextButton* scanBtn = new TGTextButton(cacheRow, "Scan");
            cacheRow->AddFrame(scanBtn, new TGLayoutHints(kLHintsCenterY));
            scanBtn->Connect("Clicked()", "GammaFitGUI", this, "OnDecayScanCaches()");
            scanBtn->SetToolTipText("Scan fit_caches/ and populate the dropdown.");

            // Asymmetric sigma window ─────────────────────────────────────────
            TGHorizontalFrame* sigRow = new TGHorizontalFrame(grp);
            grp->AddFrame(sigRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            sigRow->AddFrame(new TGLabel(sigRow, "Cut Lo sig:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
            decaySigLoEntry_ = new TGNumberEntry(sigRow, 1.0, 5, -1,
                                                 TGNumberFormat::kNESRealFour,
                                                 TGNumberFormat::kNEAPositive);
            decaySigLoEntry_->SetWidth(55);
            sigRow->AddFrame(decaySigLoEntry_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            decaySigLoEntry_->Connect("ValueSet(Long_t)",     "GammaFitGUI", this, "OnDecayCutChanged()");
            decaySigLoEntry_->Connect("ValueChanged(Long_t)", "GammaFitGUI", this, "OnDecayCutChanged()");
            sigRow->AddFrame(new TGLabel(sigRow, "Hi sig:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
            decaySigRangeEntry_ = new TGNumberEntry(sigRow, 1.0, 5, -1,
                                                    TGNumberFormat::kNESRealFour,
                                                    TGNumberFormat::kNEAPositive);
            decaySigRangeEntry_->SetWidth(55);
            sigRow->AddFrame(decaySigRangeEntry_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            decaySigRangeEntry_->Connect("ValueSet(Long_t)",     "GammaFitGUI", this, "OnDecayCutChanged()");
            decaySigRangeEntry_->Connect("ValueChanged(Long_t)", "GammaFitGUI", this, "OnDecayCutChanged()");
            TGTextButton* refreshBtn = new TGTextButton(sigRow, "Refresh");
            sigRow->AddFrame(refreshBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
            refreshBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRefreshDecayPeaks()");
            refreshBtn->SetToolTipText("Load peaks from the selected cache.");
            TGTextButton* previewBtn = new TGTextButton(sigRow, "Preview");
            sigRow->AddFrame(previewBtn, new TGLayoutHints(kLHintsLeft));
            previewBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPreviewDecay()");
            previewBtn->SetToolTipText("Show decay projection without fitting");

            // Live cut-edge readout: keV values + bin numbers
            decayBinInfoLabel_ = new TGLabel(grp, "(select a peak)");
            grp->AddFrame(decayBinInfoLabel_,
                          new TGLayoutHints(kLHintsExpandX, 4, 2, 0, 2));

            TGHorizontalFrame* rebinRow = new TGHorizontalFrame(grp);
            grp->AddFrame(rebinRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
            rebinRow->AddFrame(new TGLabel(rebinRow, "Rebin:"),
                               new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            decayRebinEntry_ = new TGNumberEntry(rebinRow, 1, 4, -1,
                                                 TGNumberFormat::kNESInteger,
                                                 TGNumberFormat::kNEAPositive,
                                                 TGNumberFormat::kNELLimitMinMax, 1, 1024);
            decayRebinEntry_->SetWidth(55);
            rebinRow->AddFrame(decayRebinEntry_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            TGTextButton* rebinApply = new TGTextButton(rebinRow, "Apply");
            rebinApply->Connect("Clicked()", "GammaFitGUI", this, "OnPreviewDecay()");
            rebinRow->AddFrame(rebinApply, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            TGTextButton* rebinReset = new TGTextButton(rebinRow, "Reset");
            rebinReset->Connect("Clicked()", "GammaFitGUI", this, "OnDecayRebinReset()");
            rebinRow->AddFrame(rebinReset, new TGLayoutHints(kLHintsLeft));

            decayPeakList_ = new TGListBox(grp, 802);
            decayPeakList_->Resize(280, 110);
            grp->AddFrame(decayPeakList_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            decayPeakList_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                    "OnDecayPeakSelected(Int_t)");

            TGTextButton* previewPeakBtn = new TGTextButton(grp, " Preview Gamma Peak ");
            grp->AddFrame(previewPeakBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
            previewPeakBtn->Connect("Clicked()", "GammaFitGUI", this, "OnDecayPreviewGammaPeak()");
            previewPeakBtn->SetToolTipText(
                "Show gamma spectrum zoomed to peak with asymmetric sigma-cut window marked");

            TGHorizontalFrame* lblRow = new TGHorizontalFrame(grp);
            grp->AddFrame(lblRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 0));
            lblRow->AddFrame(new TGLabel(lblRow, "Label:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            decayLabelEntry_ = new TGTextEntry(lblRow, "");
            decayLabelEntry_->SetWidth(120);
            lblRow->AddFrame(decayLabelEntry_, new TGLayoutHints(kLHintsLeft, 0, 8, 0, 0));

            TGTextButton* applyLblBtn = new TGTextButton(grp, "Apply Label to Peak");
            grp->AddFrame(applyLblBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            applyLblBtn->Connect("Clicked()", "GammaFitGUI", this, "OnDecayApplyLabel()");
            applyLblBtn->SetToolTipText("Save label to the cache entry for this peak");

            TGHorizontalFrame* cacheRow2 = new TGHorizontalFrame(grp);
            grp->AddFrame(cacheRow2, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));
            TGTextButton* loadBtn = new TGTextButton(cacheRow2, "Load Cache");
            cacheRow2->AddFrame(loadBtn, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 0, 2, 0, 0));
            loadBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLoadDecayCache()");
            TGTextButton* saveBtn = new TGTextButton(cacheRow2, "Save Cache");
            cacheRow2->AddFrame(saveBtn, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 2, 0, 0, 0));
            saveBtn->Connect("Clicked()", "GammaFitGUI", this, "OnSaveDecayCache()");
        }

        // ── Additional Gamma Cuts (summed into decay histogram) ───────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p2, "Additional Gamma Cuts");
            p2->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

            // Lo / Hi keV row
            TGHorizontalFrame* row1 = new TGHorizontalFrame(grp);
            grp->AddFrame(row1, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 2));
            row1->AddFrame(new TGLabel(row1, "Lo keV:"),
                           new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
            decayExtraCutLo_ = new TGNumberEntry(row1, 0.0, 8, -1,
                TGNumberFormat::kNESRealFour, TGNumberFormat::kNEAAnyNumber);
            decayExtraCutLo_->SetWidth(72);
            row1->AddFrame(decayExtraCutLo_, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));
            row1->AddFrame(new TGLabel(row1, "Hi keV:"),
                           new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
            decayExtraCutHi_ = new TGNumberEntry(row1, 0.0, 8, -1,
                TGNumberFormat::kNESRealFour, TGNumberFormat::kNEAAnyNumber);
            decayExtraCutHi_->SetWidth(72);
            row1->AddFrame(decayExtraCutHi_, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));
            TGTextButton* addCutBtn = new TGTextButton(row1, "Add");
            row1->AddFrame(addCutBtn, new TGLayoutHints(kLHintsCenterY));
            addCutBtn->Connect("Clicked()", "GammaFitGUI", this, "OnAddExtraDecayCut()");
            addCutBtn->SetToolTipText("Add this keV range to the decay projection sum");

            // List of added cuts
            decayExtraCutList_ = new TGListBox(grp, 870);
            decayExtraCutList_->Resize(280, 72);
            grp->AddFrame(decayExtraCutList_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

            TGHorizontalFrame* row2 = new TGHorizontalFrame(grp);
            grp->AddFrame(row2, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
            TGTextButton* remCutBtn = new TGTextButton(row2, " Remove ");
            row2->AddFrame(remCutBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            remCutBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRemoveExtraDecayCut()");
            remCutBtn->SetToolTipText("Remove selected extra cut");
            TGTextButton* clrCutBtn = new TGTextButton(row2, " Clear All ");
            row2->AddFrame(clrCutBtn, new TGLayoutHints(kLHintsLeft));
            clrCutBtn->Connect("Clicked()", "GammaFitGUI", this, "OnClearExtraDecayCuts()");
            clrCutBtn->SetToolTipText("Remove all extra cuts");
        }

        // ── Peak Counts vs Time ───────────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p2, "Peak Counts vs Time");
            p2->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));
            TGHorizontalFrame* row = new TGHorizontalFrame(grp);
            grp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            row->AddFrame(new TGLabel(row, "Slice width:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            decaySliceWidthEntry_ = new TGNumberEntry(row, 10.0, 8, -1,
                                                      TGNumberFormat::kNESRealFour,
                                                      TGNumberFormat::kNEAPositive);
            decaySliceWidthEntry_->SetWidth(80);
            row->AddFrame(decaySliceWidthEntry_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            row->AddFrame(new TGLabel(row, "ms"), new TGLayoutHints(kLHintsCenterY));
            TGTextButton* plotBtn = new TGTextButton(grp, "Plot Peak Counts vs Time");
            grp->AddFrame(plotBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            plotBtn->Connect("Clicked()", "GammaFitGUI", this, "OnMakePeakCountVsTime()");
            plotBtn->SetToolTipText(
                "For each time slice (x>0), project gamma axis and fit the selected peak; plot counts vs time");
        }

        // ── Decay Model (per-peak fitter, merged from old Fitter sub-tab) ─────
        {
            TGGroupFrame* grp = new TGGroupFrame(p2, "Decay Model");
            p2->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

            TGHorizontalFrame* modelRow = new TGHorizontalFrame(grp);
            grp->AddFrame(modelRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            modelRow->AddFrame(new TGLabel(modelRow, "Model:"),
                               new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
            decayModelCombo_ = new TGComboBox(modelRow, 803);
            decayModelCombo_->AddEntry("Parent + BG",              1);
            decayModelCombo_->AddEntry("Daughter (beta-) + BG",       2);
            decayModelCombo_->AddEntry("Granddaughter (beta-) + BG",  3);
            decayModelCombo_->AddEntry("Background only",          4);
            decayModelCombo_->AddEntry("Daughter (beta-n) + BG",      5);
            decayModelCombo_->AddEntry("Daughter (beta-2n) + BG",     6);
            decayModelCombo_->AddEntry("Granddaughter (beta-n)+BG",   7);
            decayModelCombo_->AddEntry("Granddaughter (beta-2n)+BG",  8);
            decayModelCombo_->Select(1, kFALSE);
            decayModelCombo_->Resize(200, 22);
            modelRow->AddFrame(decayModelCombo_, new TGLayoutHints(kLHintsLeft));
            decayModelCombo_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                      "OnDecayModelChanged(Int_t)");

            TGHorizontalFrame* bgRow = new TGHorizontalFrame(grp);
            grp->AddFrame(bgRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            bgRow->AddFrame(new TGLabel(bgRow, "BG type:"),
                            new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
            decayBGTypeCombo_ = new TGComboBox(bgRow, 805);
            decayBGTypeCombo_->AddEntry("Flat",         1);
            decayBGTypeCombo_->AddEntry("Flat + Exp",   2);
            decayBGTypeCombo_->AddEntry("Exp only",     3);
            decayBGTypeCombo_->Select(1, kFALSE);
            decayBGTypeCombo_->Resize(120, 22);
            bgRow->AddFrame(decayBGTypeCombo_, new TGLayoutHints(kLHintsLeft));
            decayBGTypeCombo_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                       "OnDecayBGTypeChanged(Int_t)");

            decayEquationLabel_ = new TGLabel(grp, "f(t) = A*exp(-ln2*t/T_P) + BG");
            decayEquationLabel_->SetTextJustify(kTextLeft);
            grp->AddFrame(decayEquationLabel_,
                          new TGLayoutHints(kLHintsExpandX, 6, 2, 2, 4));

            AddHalfLifeRow(grp, "T1/2 Parent:",    decayThalfPLbl_, decayThalfP_,    decayFixP_);
            AddBoundsRow(grp, decayThalfPLo_, decayThalfPHi_);
            AddHalfLifeRow(grp, "T1/2 Daughter:",  decayThalfDLbl_, decayThalfD_,    decayFixD_);
            AddBoundsRow(grp, decayThalfDLo_, decayThalfDHi_);
            AddHalfLifeRow(grp, "T1/2 GDaughter:", decayThalfGLbl_, decayThalfG_,    decayFixG_);
            AddBoundsRow(grp, decayThalfGLo_, decayThalfGHi_);
            AddHalfLifeRow(grp, "T1/2 Exp BG:",    decayThalfBGLbl_,decayThalfBGExp_,decayFixBGExp_, 1000.0);
            AddBoundsRow(grp, decayThalfBGExpLo_, decayThalfBGExpHi_, 1e-6, 1e15);
            decayThalfBGExp_->SetState(kFALSE);
            if (decayFixBGExp_) decayFixBGExp_->SetEnabled(kFALSE);

            TGHorizontalFrame* autoRow = new TGHorizontalFrame(grp);
            grp->AddFrame(autoRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 2));
            TGTextButton* autoBtn = new TGTextButton(autoRow, " Auto Model ");
            autoRow->AddFrame(autoBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            autoBtn->Connect("Clicked()", "GammaFitGUI", this, "OnDecayAutoModel()");
            autoBtn->SetToolTipText("Auto-select model from classification; seeds T1/2.");
            TGTextButton* seedBtn = new TGTextButton(autoRow, " Seed T1/2 ");
            autoRow->AddFrame(seedBtn, new TGLayoutHints(kLHintsLeft));
            seedBtn->Connect("Clicked()", "GammaFitGUI", this, "OnDecaySeedHalfLives()");
            seedBtn->SetToolTipText("Seed T1/2 from nuclear DB (traverses chain).");

            {
                TGHorizontalFrame* rr = new TGHorizontalFrame(grp);
                grp->AddFrame(rr, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 0));
                rr->AddFrame(new TGLabel(rr, "Fit from:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
                decayFitLo_ = new TGNumberEntry(rr, 0.0, 8, -1,
                                                TGNumberFormat::kNESRealFour,
                                                TGNumberFormat::kNEAAnyNumber);
                decayFitLo_->SetWidth(70);
                rr->AddFrame(decayFitLo_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
                rr->AddFrame(new TGLabel(rr, "to:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
                decayFitHi_ = new TGNumberEntry(rr, 0.0, 8, -1,
                                                TGNumberFormat::kNESRealFour,
                                                TGNumberFormat::kNEAAnyNumber);
                decayFitHi_->SetWidth(70);
                rr->AddFrame(decayFitHi_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
                decayFitFullRange_ = new TGCheckButton(rr, "Full");
                decayFitFullRange_->SetState(kButtonDown);
                rr->AddFrame(decayFitFullRange_, new TGLayoutHints(kLHintsCenterY));
            }

            {
                TGHorizontalFrame* rr = new TGHorizontalFrame(grp);
                grp->AddFrame(rr, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 2));
                decayFitMethod_ = new TGComboBox(rr, 840);
                decayFitMethod_->AddEntry("Chi2",             1);
                decayFitMethod_->AddEntry("Likelihood",       2);
                decayFitMethod_->AddEntry("Chi2+MINOS",       3);
                decayFitMethod_->AddEntry("Likelihood+MINOS", 4);
                decayFitMethod_->Select(2, kFALSE);
                decayFitMethod_->Resize(130, 22);
                rr->AddFrame(decayFitMethod_, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));
                decayShowResid_ = new TGCheckButton(rr, "Residuals");
                rr->AddFrame(decayShowResid_, new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
                decayErrBars_ = new TGCheckButton(rr, "Err bars");
                decayErrBars_->SetOn(kTRUE);
                rr->AddFrame(decayErrBars_, new TGLayoutHints(kLHintsCenterY));
                decayErrBars_->Connect("Toggled(Bool_t)", "GammaFitGUI", this, "OnDecayErrBarsToggled()");
            }
            {
                TGHorizontalFrame* lr = new TGHorizontalFrame(grp);
                grp->AddFrame(lr, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
                lr->AddFrame(new TGLabel(lr, "Leg X:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
                decayLegX_ = new TGNumberEntry(lr, 0.38, 5, -1,
                    TGNumberFormat::kNESRealFour, TGNumberFormat::kNEANonNegative,
                    TGNumberFormat::kNELLimitMinMax, 0.0, 0.95);
                decayLegX_->SetWidth(55);
                lr->AddFrame(decayLegX_, new TGLayoutHints(kLHintsLeft, 0, 8, 0, 0));
                lr->AddFrame(new TGLabel(lr, "Y:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
                decayLegY_ = new TGNumberEntry(lr, 0.72, 5, -1,
                    TGNumberFormat::kNESRealFour, TGNumberFormat::kNEANonNegative,
                    TGNumberFormat::kNELLimitMinMax, 0.0, 0.95);
                decayLegY_->SetWidth(55);
                lr->AddFrame(decayLegY_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
                lr->AddFrame(new TGLabel(lr, "(draggable in canvas)"),
                             new TGLayoutHints(kLHintsCenterY, 0, 0, 0, 0));
            }

            TGTextButton* fitBtn = new TGTextButton(grp, "  Fit Decay  ");
            grp->AddFrame(fitBtn, new TGLayoutHints(kLHintsCenterX, 0, 0, 4, 2));
            fitBtn->Connect("Clicked()", "GammaFitGUI", this, "OnFitDecay()");

            decayChiScanBtn_ = new TGTextButton(grp, " χ² vs Params ");
            decayChiScanBtn_->SetState(kButtonDisabled);
            grp->AddFrame(decayChiScanBtn_, new TGLayoutHints(kLHintsCenterX, 0, 0, 0, 4));
            decayChiScanBtn_->Connect("Clicked()", "GammaFitGUI", this, "OnChiScanDecay()");
        }

        // ── Fit Results ───────────────────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p2, "Fit Results");
            p2->AddFrame(grp, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 4, 4, 2, 4));
            decayResultView_ = new TGTextView(grp, 280, 180);
            grp->AddFrame(decayResultView_,
                          new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
        }

        // ── Cache Browser ─────────────────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p2, "Cache Browser");
            p2->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

            TGLabel* hint = new TGLabel(grp, "Click an entry to replay fit:");
            grp->AddFrame(hint, new TGLayoutHints(kLHintsLeft, 4, 2, 2, 0));

            decayCacheBrowserList_ = new TGListBox(grp, 871);
            decayCacheBrowserList_->Resize(280, 120);
            grp->AddFrame(decayCacheBrowserList_,
                          new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            decayCacheBrowserList_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                            "OnDecayCacheBrowserSelected(Int_t)");
        }
    } // end Cuts sub-tab

    // ═══════════════════════════════════════════════════════════════════════════
    // Sub-tab 2: Total Decay  -  sum multiple peaks, Bateman fit, population
    // ═══════════════════════════════════════════════════════════════════════════
    {
        TGCompositeFrame* tf = decaySubTabs_->AddTab("Total Decay");
        TGCanvas* sc = new TGCanvas(tf, 305, 820, kSunkenFrame);
        tf->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        TGCompositeFrame* cf = new TGCompositeFrame(sc->GetViewPort(), 293, 10, kVerticalFrame);
        sc->SetContainer(cf);
        TGCompositeFrame* p2 = cf;

        // ── Peak Selection ────────────────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p2, "Peak Selection (uses Cuts peak list)");
            p2->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));
            decayTdSumAll_ = new TGCheckButton(grp, "Sum all peaks in list");
            decayTdSumAll_->SetState(kButtonDown);
            grp->AddFrame(decayTdSumAll_, new TGLayoutHints(kLHintsLeft, 4, 2, 4, 2));
            decayTdSumAll_->SetToolTipText(
                "Checked: sum projections from every peak in the Cuts peak list.\n"
                "Unchecked: use only the selected peak in the Cuts list.");
        }

        // ── Decay Model ───────────────────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p2, "Decay Model");
            p2->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

            TGHorizontalFrame* modelRow = new TGHorizontalFrame(grp);
            grp->AddFrame(modelRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            modelRow->AddFrame(new TGLabel(modelRow, "Model:"),
                               new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
            decayTdModelCombo_ = new TGComboBox(modelRow, 810);
            decayTdModelCombo_->AddEntry("Parent + BG",              1);
            decayTdModelCombo_->AddEntry("Daughter (beta-) + BG",       2);
            decayTdModelCombo_->AddEntry("Granddaughter (beta-) + BG",  3);
            decayTdModelCombo_->AddEntry("Background only",          4);
            decayTdModelCombo_->AddEntry("Daughter (beta-n) + BG",      5);
            decayTdModelCombo_->AddEntry("Daughter (beta-2n) + BG",     6);
            decayTdModelCombo_->AddEntry("Granddaughter (beta-n)+BG",   7);
            decayTdModelCombo_->AddEntry("Granddaughter (beta-2n)+BG",  8);
            decayTdModelCombo_->Select(1, kFALSE);
            decayTdModelCombo_->Resize(200, 22);
            modelRow->AddFrame(decayTdModelCombo_, new TGLayoutHints(kLHintsLeft));
            decayTdModelCombo_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                        "OnDecayTdModelChanged(Int_t)");

            TGHorizontalFrame* bgRow = new TGHorizontalFrame(grp);
            grp->AddFrame(bgRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            bgRow->AddFrame(new TGLabel(bgRow, "BG type:"),
                            new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
            decayTdBGCombo_ = new TGComboBox(bgRow, 811);
            decayTdBGCombo_->AddEntry("Flat",         1);
            decayTdBGCombo_->AddEntry("Flat + Exp",   2);
            decayTdBGCombo_->AddEntry("Exp only",     3);
            decayTdBGCombo_->Select(1, kFALSE);
            decayTdBGCombo_->Resize(120, 22);
            bgRow->AddFrame(decayTdBGCombo_, new TGLayoutHints(kLHintsLeft));
            decayTdBGCombo_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                     "OnDecayTdBGTypeChanged(Int_t)");

            decayTdEquationLbl_ = new TGLabel(grp, "f(t) = A*exp(-ln2*t/T_P) + BG");
            decayTdEquationLbl_->SetTextJustify(kTextLeft);
            grp->AddFrame(decayTdEquationLbl_,
                          new TGLayoutHints(kLHintsExpandX, 6, 2, 2, 4));

            AddHalfLifeRow(grp, "T1/2 Parent:",    decayTdPLbl_, decayTdThalfP_,    decayTdFixP_);
            AddHalfLifeRow(grp, "T1/2 Daughter:",  decayTdDLbl_, decayTdThalfD_,    decayTdFixD_);
            AddHalfLifeRow(grp, "T1/2 GDaughter:", decayTdGLbl_, decayTdThalfG_,    decayTdFixG_);
            AddHalfLifeRow(grp, "T1/2 Exp BG:",    decayTdBGLbl_,decayTdThalfBGExp_,decayTdFixBGExp_, 1000.0);
            decayTdThalfBGExp_->SetState(kFALSE);
            if (decayTdFixBGExp_) decayTdFixBGExp_->SetEnabled(kFALSE);

            {
                TGHorizontalFrame* rr = new TGHorizontalFrame(grp);
                grp->AddFrame(rr, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 0));
                rr->AddFrame(new TGLabel(rr, "Fit from:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
                decayTdFitLo_ = new TGNumberEntry(rr, 0.0, 8, -1,
                                                  TGNumberFormat::kNESRealFour,
                                                  TGNumberFormat::kNEAAnyNumber);
                decayTdFitLo_->SetWidth(70);
                rr->AddFrame(decayTdFitLo_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
                rr->AddFrame(new TGLabel(rr, "to:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
                decayTdFitHi_ = new TGNumberEntry(rr, 0.0, 8, -1,
                                                  TGNumberFormat::kNESRealFour,
                                                  TGNumberFormat::kNEAAnyNumber);
                decayTdFitHi_->SetWidth(70);
                rr->AddFrame(decayTdFitHi_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
                decayTdFullRange_ = new TGCheckButton(rr, "Full");
                decayTdFullRange_->SetState(kButtonDown);
                rr->AddFrame(decayTdFullRange_, new TGLayoutHints(kLHintsCenterY));
            }

            // ── Rebin, method, residuals row ─────────────────────────────────
            {
                TGHorizontalFrame* rr = new TGHorizontalFrame(grp);
                grp->AddFrame(rr, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 2));

                rr->AddFrame(new TGLabel(rr, "Rebin:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
                decayTdRebinEntry_ = new TGNumberEntry(rr, 1, 4, -1,
                    TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive,
                    TGNumberFormat::kNELLimitMinMax, 1, 1024);
                decayTdRebinEntry_->SetWidth(50);
                rr->AddFrame(decayTdRebinEntry_, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));

                decayTdFitMethod_ = new TGComboBox(rr, 820);
                decayTdFitMethod_->AddEntry("Chi2",               1);
                decayTdFitMethod_->AddEntry("Likelihood",         2);
                decayTdFitMethod_->AddEntry("Chi2+MINOS",         3);
                decayTdFitMethod_->AddEntry("Likelihood+MINOS",   4);
                decayTdFitMethod_->Select(2, kFALSE);
                decayTdFitMethod_->Resize(130, 22);
                rr->AddFrame(decayTdFitMethod_, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));

                decayTdShowResid_ = new TGCheckButton(rr, "Residuals");
                rr->AddFrame(decayTdShowResid_, new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
                decayTdErrBars_ = new TGCheckButton(rr, "Err bars");
                decayTdErrBars_->SetOn(kTRUE);
                rr->AddFrame(decayTdErrBars_, new TGLayoutHints(kLHintsCenterY));
            }
            {
                TGHorizontalFrame* lr = new TGHorizontalFrame(grp);
                grp->AddFrame(lr, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
                lr->AddFrame(new TGLabel(lr, "Leg X:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
                decayTdLegX_ = new TGNumberEntry(lr, 0.38, 5, -1,
                    TGNumberFormat::kNESRealFour, TGNumberFormat::kNEANonNegative,
                    TGNumberFormat::kNELLimitMinMax, 0.0, 0.95);
                decayTdLegX_->SetWidth(55);
                lr->AddFrame(decayTdLegX_, new TGLayoutHints(kLHintsLeft, 0, 8, 0, 0));
                lr->AddFrame(new TGLabel(lr, "Y:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
                decayTdLegY_ = new TGNumberEntry(lr, 0.72, 5, -1,
                    TGNumberFormat::kNESRealFour, TGNumberFormat::kNEANonNegative,
                    TGNumberFormat::kNELLimitMinMax, 0.0, 0.95);
                decayTdLegY_->SetWidth(55);
                lr->AddFrame(decayTdLegY_, new TGLayoutHints(kLHintsLeft));
            }

            TGTextButton* fitTdBtn = new TGTextButton(grp, "  Fit Total Decay  ");
            grp->AddFrame(fitTdBtn, new TGLayoutHints(kLHintsCenterX, 0, 0, 4, 2));
            fitTdBtn->Connect("Clicked()", "GammaFitGUI", this, "OnFitTotalDecay()");

            decayTdChiScanBtn_ = new TGTextButton(grp, " χ² vs Params ");
            decayTdChiScanBtn_->SetState(kButtonDisabled);
            grp->AddFrame(decayTdChiScanBtn_, new TGLayoutHints(kLHintsCenterX, 0, 0, 0, 4));
            decayTdChiScanBtn_->Connect("Clicked()", "GammaFitGUI", this, "OnChiScanTotalDecay()");
        }

        // ── Fit Results ───────────────────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p2, "Fit Results");
            p2->AddFrame(grp, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 4, 4, 2, 4));
            decayTdResultView_ = new TGTextView(grp, 280, 200);
            grp->AddFrame(decayTdResultView_,
                          new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
        }
    } // end Total Decay sub-tab

    // ═══════════════════════════════════════════════════════════════════════════
    // Sub-tab 3: RK4 Full Chain  -  4-species ODE decay chain solved by RK4
    // ═══════════════════════════════════════════════════════════════════════════
    {
        TGCompositeFrame* tf = decaySubTabs_->AddTab("RK4 Chain");
        TGCanvas* sc = new TGCanvas(tf, 305, 820, kSunkenFrame);
        tf->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        TGCompositeFrame* cf = new TGCompositeFrame(sc->GetViewPort(), 293, 10, kVerticalFrame);
        sc->SetContainer(cf);
        TGCompositeFrame* p3 = cf;

        // ── Data Source ───────────────────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p3, "Data Source");
            p3->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

            // ── Mode A: TH2 peak projection sum ──────────────────────────────
            grp->AddFrame(new TGLabel(grp, "Mode A — sum TH2 projections (multi-select):"),
                          new TGLayoutHints(kLHintsLeft, 4, 2, 4, 2));

            rkPeakList_ = new TGListBox(grp, 841);
            rkPeakList_->SetMultipleSelections(kTRUE);
            rkPeakList_->Resize(280, 80);
            grp->AddFrame(rkPeakList_, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));

            TGHorizontalFrame* btnRow = new TGHorizontalFrame(grp);
            grp->AddFrame(btnRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
            TGTextButton* allBtn = new TGTextButton(btnRow, " All ");
            allBtn->SetToolTipText("Select all peaks");
            allBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRK4SelectAll()");
            btnRow->AddFrame(allBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            TGTextButton* noneBtn = new TGTextButton(btnRow, " None ");
            noneBtn->SetToolTipText("Deselect all peaks");
            noneBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRK4SelectNone()");
            btnRow->AddFrame(noneBtn, new TGLayoutHints(kLHintsLeft, 0, 8, 0, 0));
            TGTextButton* seedBtn = new TGTextButton(btnRow, " Seed T1/2 from Cuts ");
            seedBtn->SetToolTipText("Copy fitted T1/2 from the currently selected Cuts peak fit result");
            seedBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRK4SeedFromCuts()");
            btnRow->AddFrame(seedBtn, new TGLayoutHints(kLHintsLeft));

            // ── Mode B: direct 1D decay histogram ────────────────────────────
            TGHorizontalFrame* sepRow = new TGHorizontalFrame(grp);
            grp->AddFrame(sepRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 0));
            rkUse1DHist_ = new TGCheckButton(sepRow, "Mode B — use 1D histogram directly:");
            rkUse1DHist_->SetToolTipText(
                "When checked, ignore the peak list above and fit a pre-made\n"
                "decay histogram picked from the combo below.\n"
                "Scan finds TH1s in the open ROOT file and saved decay-slice files.");
            sepRow->AddFrame(rkUse1DHist_, new TGLayoutHints(kLHintsCenterY));

            TGHorizontalFrame* histRow = new TGHorizontalFrame(grp);
            grp->AddFrame(histRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));
            rkHistCombo_ = new TGComboBox(histRow, 843);
            rkHistCombo_->Resize(190, 22);
            histRow->AddFrame(rkHistCombo_, new TGLayoutHints(kLHintsLeft | kLHintsCenterY, 0, 4, 0, 0));
            TGTextButton* scanBtn = new TGTextButton(histRow, " Scan ");
            scanBtn->SetToolTipText("Scan open ROOT file and decay-slice caches for 1D histograms");
            scanBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRK4ScanHistograms()");
            histRow->AddFrame(scanBtn, new TGLayoutHints(kLHintsCenterY));
        }

        // ── Chain Model + equation display ────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p3, "Chain Model");
            p3->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

            TGHorizontalFrame* mrow = new TGHorizontalFrame(grp);
            grp->AddFrame(mrow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            mrow->AddFrame(new TGLabel(mrow, "Model:"),
                           new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            rkChainModel_ = new TGComboBox(mrow, 842);
            rkChainModel_->AddEntry("Full: P->D + P->Bn->Gd", 1);
            rkChainModel_->AddEntry("P->D only (no beta-n)",   2);
            rkChainModel_->AddEntry("Parent only",             3);
            rkChainModel_->Select(1, kFALSE);
            rkChainModel_->Resize(190, 22);
            mrow->AddFrame(rkChainModel_, new TGLayoutHints(kLHintsLeft));

            TGTextView* eqv = new TGTextView(grp, 280, 78);
            eqv->AddLine("dN_P/dt  = -lam_P*N_P");
            eqv->AddLine("dN_D/dt  = (1-Pn)*lam_P*N_P - lam_D*N_D");
            eqv->AddLine("dN_Bn/dt = Pn*lam_P*N_P    - lam_Bn*N_Bn");
            eqv->AddLine("dN_Gd/dt = lam_Bn*N_Bn     - lam_Gd*N_Gd");
            eqv->AddLine("f(t) = lam_P*N_P + lam_D*N_D + lam_Bn*N_Bn + lam_Gd*N_Gd + bg");
            grp->AddFrame(eqv, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        }

        // ── Half-lives ────────────────────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p3, "Half-lives (ms)");
            p3->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));
            AddHalfLifeRow(grp, "T1/2 Parent:",    rkPLbl_,  rkThalfP_,  rkFixP_,   100.0);
            AddHalfLifeRow(grp, "T1/2 Daughter:",  rkDLbl_,  rkThalfD_,  rkFixD_,   200.0);
            AddHalfLifeRow(grp, "T1/2 beta-n br:", rkBnLbl_, rkThalfBn_, rkFixBn_,  500.0);
            AddHalfLifeRow(grp, "T1/2 GDaughter:", rkGdLbl_, rkThalfGd_, rkFixGd_, 1000.0);
        }

        // ── Pn (beta-n probability) ───────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p3, "Beta-n Branching");
            p3->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));
            TGHorizontalFrame* row = new TGHorizontalFrame(grp);
            grp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 4));
            TGLabel* lb = new TGLabel(row, "Pn  [0,1]:");
            lb->SetWidth(78);
            row->AddFrame(lb, new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
            rkPn_ = new TGNumberEntry(row, 0.2, 8, -1,
                TGNumberFormat::kNESRealFour, TGNumberFormat::kNEANonNegative,
                TGNumberFormat::kNELLimitMinMax, 0.0, 1.0);
            rkPn_->SetWidth(80);
            row->AddFrame(rkPn_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            rkFixPn_ = new TGCheckButton(row, "Fix");
            row->AddFrame(rkFixPn_, new TGLayoutHints(kLHintsCenterY));
        }

        // ── Initial conditions ────────────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p3, "Initial Conditions");
            p3->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

            auto addIcRow = [&](const char* lbl, TGNumberEntry*& entry,
                                TGCheckButton*& fix, double defVal, bool posOnly) {
                TGHorizontalFrame* row = new TGHorizontalFrame(grp);
                grp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
                TGLabel* lb = new TGLabel(row, lbl);
                lb->SetWidth(78);
                row->AddFrame(lb, new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
                entry = new TGNumberEntry(row, defVal, 10, -1,
                    TGNumberFormat::kNESReal,
                    posOnly ? TGNumberFormat::kNEANonNegative : TGNumberFormat::kNEAAnyNumber);
                entry->SetWidth(90);
                row->AddFrame(entry, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
                fix = new TGCheckButton(row, "Fix");
                row->AddFrame(fix, new TGLayoutHints(kLHintsCenterY));
            };

            addIcRow("N0 (initial):", rkN0_, rkFixN0_, 220000.0, true);
            addIcRow("bg (flat):",    rkBg_, rkFixBg_,      0.0, false);
        }

        // ── Integration & rebin ───────────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p3, "Integration & Rebin");
            p3->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

            TGHorizontalFrame* row1 = new TGHorizontalFrame(grp);
            grp->AddFrame(row1, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            row1->AddFrame(new TGLabel(row1, "RK4 steps:"),
                           new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
            rkNsteps_ = new TGNumberEntry(row1, 200, 5, -1,
                TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive,
                TGNumberFormat::kNELLimitMinMax, 10, 10000);
            rkNsteps_->SetWidth(60);
            row1->AddFrame(rkNsteps_, new TGLayoutHints(kLHintsLeft, 0, 10, 0, 0));
            row1->AddFrame(new TGLabel(row1, "Rebin:"),
                           new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            rkRebinEntry_ = new TGNumberEntry(row1, 1, 4, -1,
                TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive,
                TGNumberFormat::kNELLimitMinMax, 1, 1024);
            rkRebinEntry_->SetWidth(50);
            row1->AddFrame(rkRebinEntry_, new TGLayoutHints(kLHintsLeft));
        }

        // ── Fit Range ─────────────────────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p3, "Fit Range (ms)");
            p3->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));
            TGHorizontalFrame* rr = new TGHorizontalFrame(grp);
            grp->AddFrame(rr, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            rr->AddFrame(new TGLabel(rr, "From:"), new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            rkFitLo_ = new TGNumberEntry(rr, 0.0, 8, -1,
                TGNumberFormat::kNESRealFour, TGNumberFormat::kNEAAnyNumber);
            rkFitLo_->SetWidth(65);
            rr->AddFrame(rkFitLo_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            rr->AddFrame(new TGLabel(rr, "to:"), new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            rkFitHi_ = new TGNumberEntry(rr, 0.0, 8, -1,
                TGNumberFormat::kNESRealFour, TGNumberFormat::kNEAAnyNumber);
            rkFitHi_->SetWidth(65);
            rr->AddFrame(rkFitHi_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            rkFitFull_ = new TGCheckButton(rr, "Full");
            rkFitFull_->SetState(kButtonDown);
            rr->AddFrame(rkFitFull_, new TGLayoutHints(kLHintsCenterY));
        }

        // ── Fit method, residuals, fit button ────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p3, "Fit Options");
            p3->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

            TGHorizontalFrame* optRow = new TGHorizontalFrame(grp);
            grp->AddFrame(optRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));

            optRow->AddFrame(new TGLabel(optRow, "Method:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            rkFitMethod_ = new TGComboBox(optRow, 830);
            rkFitMethod_->AddEntry("Chi2",               1);
            rkFitMethod_->AddEntry("Likelihood",         2);
            rkFitMethod_->AddEntry("Chi2+MINOS",         3);
            rkFitMethod_->AddEntry("Likelihood+MINOS",   4);
            rkFitMethod_->Select(2, kFALSE);
            rkFitMethod_->Resize(130, 22);
            optRow->AddFrame(rkFitMethod_, new TGLayoutHints(kLHintsLeft, 0, 8, 0, 0));

            rkShowResid_ = new TGCheckButton(optRow, "Residuals");
            optRow->AddFrame(rkShowResid_, new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
            rkErrBars_ = new TGCheckButton(optRow, "Err bars");
            rkErrBars_->SetOn(kTRUE);
            optRow->AddFrame(rkErrBars_, new TGLayoutHints(kLHintsCenterY));

            TGHorizontalFrame* lr = new TGHorizontalFrame(grp);
            grp->AddFrame(lr, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
            lr->AddFrame(new TGLabel(lr, "Leg X:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
            rkLegX_ = new TGNumberEntry(lr, 0.38, 5, -1,
                TGNumberFormat::kNESRealFour, TGNumberFormat::kNEANonNegative,
                TGNumberFormat::kNELLimitMinMax, 0.0, 0.95);
            rkLegX_->SetWidth(55);
            lr->AddFrame(rkLegX_, new TGLayoutHints(kLHintsLeft, 0, 8, 0, 0));
            lr->AddFrame(new TGLabel(lr, "Y:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
            rkLegY_ = new TGNumberEntry(lr, 0.66, 5, -1,
                TGNumberFormat::kNESRealFour, TGNumberFormat::kNEANonNegative,
                TGNumberFormat::kNELLimitMinMax, 0.0, 0.95);
            rkLegY_->SetWidth(55);
            lr->AddFrame(rkLegY_, new TGLayoutHints(kLHintsLeft));

            TGTextButton* btn = new TGTextButton(grp, "  Fit RK4 Decay Chain  ");
            grp->AddFrame(btn, new TGLayoutHints(kLHintsCenterX, 0, 0, 2, 2));
            btn->Connect("Clicked()", "GammaFitGUI", this, "OnFitRK4Decay()");

            rkChiScanBtn_ = new TGTextButton(grp, " χ² vs Params ");
            rkChiScanBtn_->SetState(kButtonDisabled);
            grp->AddFrame(rkChiScanBtn_, new TGLayoutHints(kLHintsCenterX, 0, 0, 0, 4));
            rkChiScanBtn_->Connect("Clicked()", "GammaFitGUI", this, "OnChiScanRK4()");
        }

        // ── Results ───────────────────────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p3, "Fit Results");
            p3->AddFrame(grp, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 4, 4, 2, 4));
            rkResultView_ = new TGTextView(grp, 280, 200);
            grp->AddFrame(rkResultView_,
                          new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
        }
    } // end RK4 Chain sub-tab
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
    decayPeakCacheNames_.clear();
    decayFitStore_.clear();
    decayGammaProjName_.clear();
    decayPeakList_->RemoveAll();
    decayPeakList_->MapSubwindows();
    decayPeakList_->Layout();
    if (decayLabelEntry_) decayLabelEntry_->SetText("");
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper: compose a human-readable equation for a given (modelId, bgType).
// The result fits in one TGLabel line.
// ─────────────────────────────────────────────────────────────────────────────
static std::string DecayEquationText(int modelId, int bgType)
{
    // Signal part
    static const char* kSig[] = {
        "",
        "A*exp(-ln2*t/T_P)",
        "A*lD/(lD-lP)*(e^-lP*t - e^-lD*t)",
        "Bateman(P->D->G) Ingrowth",
        "(no signal)",
        "Bateman-bn Daughter",
        "Bateman-b2n Daughter",
        "Bateman-bn GDaughter",
        "Bateman-b2n GDaughter",
    };
    std::string sig = (modelId >= 1 && modelId <= 8) ? kSig[modelId] : "";

    // BG part
    std::string bg;
    if (bgType == 1) bg = "bg";
    else if (bgType == 2) bg = "bg + A_bg*exp(-ln2*t/T_bg)";
    else                  bg = "A_bg*exp(-ln2*t/T_bg)";

    if (modelId == 4) return std::string("f(t) = ") + bg;
    return std::string("f(t) = ") + sig + " + " + bg;
}

void GammaFitGUI::UpdateCutsEquation()
{
    if (!decayEquationLabel_) return;
    int model  = decayModelCombo_   ? decayModelCombo_->GetSelected()   : 1;
    int bgType = decayBGTypeCombo_  ? decayBGTypeCombo_->GetSelected()  : 1;
    decayEquationLabel_->SetText(DecayEquationText(model, bgType).c_str());
    decayEquationLabel_->Layout();
}

void GammaFitGUI::UpdateTotalDecayEquation()
{
    if (!decayTdEquationLbl_) return;
    int model  = decayTdModelCombo_ ? decayTdModelCombo_->GetSelected() : 1;
    int bgType = decayTdBGCombo_    ? decayTdBGCombo_->GetSelected()    : 1;
    decayTdEquationLbl_->SetText(DecayEquationText(model, bgType).c_str());
    decayTdEquationLbl_->Layout();
}

void GammaFitGUI::OnDecayModelChanged(Int_t id)
{
    // Update T1/2 row labels to reflect the chain
    if (decayThalfPLbl_ && decayThalfDLbl_ && decayThalfGLbl_) {
        switch (id) {
            case 5:
                decayThalfPLbl_->SetText("T1/2 beta-n product:");
                decayThalfDLbl_->SetText("T1/2 Daughter:");
                break;
            case 6:
                decayThalfPLbl_->SetText("T1/2 beta-2n product:");
                decayThalfDLbl_->SetText("T1/2 Daughter:");
                break;
            case 7:
                decayThalfPLbl_->SetText("T1/2 Parent:");
                decayThalfDLbl_->SetText("T1/2 beta-n product:");
                decayThalfGLbl_->SetText("T1/2 GDaughter:");
                break;
            case 8:
                decayThalfPLbl_->SetText("T1/2 Parent:");
                decayThalfDLbl_->SetText("T1/2 beta-2n product:");
                decayThalfGLbl_->SetText("T1/2 GDaughter:");
                break;
            default:
                decayThalfPLbl_->SetText("T1/2 Parent:");
                decayThalfDLbl_->SetText("T1/2 Daughter:");
                decayThalfGLbl_->SetText("T1/2 GDaughter:");
                break;
        }
    }
    UpdateCutsEquation();
}

void GammaFitGUI::OnDecayBGTypeChanged(Int_t /*id*/)
{
    int bgType = decayBGTypeCombo_ ? decayBGTypeCombo_->GetSelected() : 1;
    bool hasExp = (bgType >= 2);
    if (decayThalfBGExp_) decayThalfBGExp_->SetState(hasExp ? kTRUE : kFALSE);
    if (decayFixBGExp_)   decayFixBGExp_->SetEnabled(hasExp ? kTRUE : kFALSE);
    UpdateCutsEquation();
}

void GammaFitGUI::OnDecayAutoModel()
{
    Int_t peakSel = decayPeakList_ ? decayPeakList_->GetSelected() : -1;
    if (peakSel < 1 || (size_t)peakSel > decayPeakEs_.size()) {
        AppendLog("[Decay] Select a peak first"); return;
    }
    // Get classification from gamma cache
    const std::string& cacheName = (peakSel <= (Int_t)decayPeakCacheNames_.size())
                                   ? decayPeakCacheNames_[peakSel-1] : decayGammaProjName_;
    FitDatabase fdb;
    fdb.Load(CacheFileFor(cacheName));
    std::string cls;
    auto it = fdb.GetEntries().find(decayPeakKeys_[peakSel-1]);
    if (it != fdb.GetEntries().end()) cls = it->second.classification;
    if (cls.empty() && decayLabelEntry_) {
        std::string lbl = decayLabelEntry_->GetText();
        if (labelClassMap_.count(lbl)) cls = labelClassMap_.at(lbl);
    }
    // Map classification string -> model ID
    static const struct { const char* key; int model; } kMap[] = {
        {"Granddaughter",       3},
        {"GDaughter",           3},
        {"Parent",              1},
        {"Daughter",            2},
        {"Background",          4},
        {"BG",                  4},
        {"bn Granddaughter",    7},
        {"b2n Granddaughter",   8},
        {"bn Daughter",         5},
        {"b2n Daughter",        6},
        {nullptr, 0}
    };
    int model = 0;
    for (int i = 0; kMap[i].key; ++i) {
        if (cls.find(kMap[i].key) != std::string::npos) { model = kMap[i].model; break; }
    }
    if (model == 0) {
        AppendLog(Form("[Decay] No model mapping for class '%s'", cls.c_str()));
        return;
    }
    if (decayModelCombo_) decayModelCombo_->Select(model, kTRUE);
    AppendLog(Form("[Decay] Auto model: '%s' -> model %d", cls.c_str(), model));
    OnDecaySeedHalfLives();
}

void GammaFitGUI::OnDecaySeedHalfLives()
{
    if (!decayLabelEntry_ || !decayThalfP_) return;
    std::string lbl = decayLabelEntry_->GetText();
    // Strip trailing spaces
    while (!lbl.empty() && lbl.back() == ' ') lbl.pop_back();
    if (lbl.empty()) { AppendLog("[Decay] Set a label first"); return; }

    int modelId  = decayModelCombo_ ? decayModelCombo_->GetSelected() : 1;
    int sigType  = DecaySignalType(modelId);

    auto findHl = [this](const std::string& id) -> double {
        // Try exact key
        auto it = nuclearDB_.find(id);
        if (it != nuclearDB_.end() && it->second.halflife_s > 0)
            return it->second.halflife_s * 1000.0;  // s -> ms
        // Try without spaces
        std::string clean = id;
        clean.erase(std::remove(clean.begin(), clean.end(), ' '), clean.end());
        it = nuclearDB_.find(clean);
        if (it != nuclearDB_.end() && it->second.halflife_s > 0)
            return it->second.halflife_s * 1000.0;
        return -1.0;
    };

    double T_this = findHl(lbl);
    if (T_this <= 0) {
        AppendLog(Form("[Decay] '%s' not found in nuclear DB  -  seed manually", lbl.c_str()));
        return;
    }

    switch (sigType) {
        case 1: // This IS the parent
            decayThalfP_->SetNumber(T_this);
            break;
        case 2: { // This IS the daughter
            if (decayThalfD_) decayThalfD_->SetNumber(T_this);
            auto pos = std::find(nucChainIsotopes_.begin(), nucChainIsotopes_.end(), lbl);
            if (pos != nucChainIsotopes_.end() && pos != nucChainIsotopes_.begin()) {
                double Tp = findHl(*std::prev(pos));
                if (Tp > 0) decayThalfP_->SetNumber(Tp);
            }
            break;
        }
        case 3: { // This IS the granddaughter
            if (decayThalfG_) decayThalfG_->SetNumber(T_this);
            auto pos = std::find(nucChainIsotopes_.begin(), nucChainIsotopes_.end(), lbl);
            if (pos != nucChainIsotopes_.end() && pos != nucChainIsotopes_.begin()) {
                auto ppos = std::prev(pos);
                double Td = findHl(*ppos);
                if (Td > 0 && decayThalfD_) decayThalfD_->SetNumber(Td);
                if (ppos != nucChainIsotopes_.begin()) {
                    double Tp = findHl(*std::prev(ppos));
                    if (Tp > 0) decayThalfP_->SetNumber(Tp);
                }
            }
            break;
        }
        default: break;
    }
    AppendLog(Form("[Decay] Seeded T1/2 for '%s' = %.4g ms", lbl.c_str(), T_this));
}

void GammaFitGUI::PopulateCacheBrowser()
{
    if (!decayCacheBrowserList_) return;
    decayCacheBrowserList_->RemoveAll();
    decayCacheBrowserEs_.clear();
    int id = 1;
    for (const auto& kv : decayFitStore_) {
        const DecayFitResult& r = kv.second;
        int sig = DecaySignalType(r.model);
        // Show: energy, T1/2(P), chi2/ndf, n extra cuts
        std::string entry = Form("%.2f keV", r.peakE);
        if (sig >= 1 && r.params.size() > 1)
            entry += Form("  T_{P}=%.4g ms", r.params[1]);
        if (r.chi2ndf >= 0)
            entry += Form("  χ²/ndf=%.3f", r.chi2ndf);
        if (!r.extraCuts.empty())
            entry += Form("  +%d cut(s)", (int)r.extraCuts.size());
        decayCacheBrowserList_->AddEntry(entry.c_str(), id);
        decayCacheBrowserEs_.push_back(r.peakE);
        ++id;
    }
    decayCacheBrowserList_->MapSubwindows();
    decayCacheBrowserList_->Layout();
}

void GammaFitGUI::OnDecayCacheBrowserSelected(Int_t id)
{
    if (id < 1 || (size_t)id > decayCacheBrowserEs_.size()) return;
    double E = decayCacheBrowserEs_[id - 1];
    auto it = decayFitStore_.find(E);
    if (it == decayFitStore_.end()) return;
    const DecayFitResult& r = it->second;

    // ── Show stats in result view ─────────────────────────────────────────────
    decayResultView_->Clear();
    int axisId = decayGammaAxisCombo_ ? decayGammaAxisCombo_->GetSelected() : 1;
    int bgType  = (r.bgType >= 1 && r.bgType <= 3) ? r.bgType : 1;
    int sig     = DecaySignalType(r.model);
    int nSig    = DecaySigNpar(sig);
    decayResultView_->AddLine(Form("Cache entry: E=%.2f keV", E));
    decayResultView_->AddLine(Form("  Window: [%.2f, %.2f keV]  Nsig: lo=%.2f hi=%.2f",
                                   r.eMin, r.eMax,
                                   r.NsigLo > 0 ? r.NsigLo : r.Nsig, r.Nsig));
    for (size_t ci = 0; ci < r.extraCuts.size(); ci++)
        decayResultView_->AddLine(Form("  ExtraCut[%zu]: [%.2f, %.2f keV]",
                                       ci+1, r.extraCuts[ci].first, r.extraCuts[ci].second));
    {
        const char* methodNames[] = {"","Chi2","Likelihood","Chi2+MINOS","Likelihood+MINOS"};
        int m = (r.fitMethod >= 1 && r.fitMethod <= 4) ? r.fitMethod : 2;
        decayResultView_->AddLine(Form("  Model: %d  bgType: %d  Method: %s  rebin: %d",
                                       r.model, bgType, methodNames[m], r.rebin));
        if (!r.fullRange)
            decayResultView_->AddLine(Form("  Fit range: [%.4g, %.4g] ms", r.fitLo, r.fitHi));
    }
    decayResultView_->AddLine(Form("  chi2/ndf=%.4g  status=%d", r.chi2ndf, r.status));
    // Parameter table
    TF1* fTmp = BuildDecayTF1("_cacheBrowse", r.model, bgType, r.fitLo > 0 ? r.fitLo : 0.0, r.fitHi > 0 ? r.fitHi : 1e6);
    if (fTmp) {
        int np = std::min((int)r.params.size(), fTmp->GetNpar());
        for (int i = 0; i < np; i++) {
            double err = (i < (int)r.errors.size()) ? r.errors[i] : 0.0;
            double plo = (i < (int)r.paramLo.size()) ? r.paramLo[i] : 0.0;
            double phi = (i < (int)r.paramHi.size()) ? r.paramHi[i] : 0.0;
            decayResultView_->AddLine(Form("  %-16s = %11.5g +/- %.4g  [%.4g, %.4g]",
                fTmp->GetParName(i), r.params[i], err, plo, phi));
        }
        delete fTmp;
    }
    if (!r.label.empty())
        decayResultView_->AddLine(Form("  Label: %s  (%s)", r.label.c_str(), r.classification.c_str()));
    decayResultView_->ShowBottom();

    // ── Replay on canvas ──────────────────────────────────────────────────────
    if (!inputFile_ || r.histName.empty()) return;
    TH2* h2 = (TH2*)inputFile_->Get(r.histName.c_str());
    if (!h2) return;

    TH1* hDecay = nullptr;
    if (axisId == 1) {
        int b1 = h2->GetXaxis()->FindBin(r.eMin), b2 = h2->GetXaxis()->FindBin(r.eMax);
        hDecay = h2->ProjectionY(Form("hCacheBrowse_%.1f", E), b1, b2);
    } else {
        int b1 = h2->GetYaxis()->FindBin(r.eMin), b2 = h2->GetYaxis()->FindBin(r.eMax);
        hDecay = h2->ProjectionX(Form("hCacheBrowse_%.1f", E), b1, b2);
    }
    if (!hDecay) return;
    hDecay->SetDirectory(nullptr);
    if (!r.extraCuts.empty()) ApplyExtraCuts(h2, axisId, hDecay, r.extraCuts);
    if (r.rebin > 1) hDecay->Rebin(r.rebin);

    double tLo = std::max(0.0, hDecay->GetXaxis()->GetXmin());
    double tHi = hDecay->GetXaxis()->GetXmax();
    double xlo = r.fullRange ? tLo : std::max(tLo, r.fitLo);
    double xhi = r.fullRange ? tHi : std::min(tHi, r.fitHi > 0 ? r.fitHi : tHi);
    hDecay->GetXaxis()->SetRangeUser(xlo, xhi);
    {
        std::string ttl = Form("E=%.2f keV  [%.2f, %.2f keV]", E, r.eMin, r.eMax);
        if (!r.extraCuts.empty())
            ttl += Form("  +%d extra cut(s)", (int)r.extraCuts.size());
        hDecay->SetTitle(ttl.c_str());
    }
    hDecay->GetXaxis()->SetTitle("Time (ms)");
    hDecay->GetYaxis()->SetTitle(r.rebin > 1
        ? Form("Counts / (%.4g ms)", hDecay->GetBinWidth(1)) : "Counts");
    hDecay->SetLineColor(kBlack); hDecay->SetMarkerSize(0);

    TF1* fReplay = BuildDecayTF1("fCacheBrowseReplay", r.model, bgType, xlo, xhi);
    if (!fReplay || (int)r.params.size() < fReplay->GetNpar()) { delete fReplay; return; }
    for (int i = 0; i < fReplay->GetNpar(); i++) fReplay->SetParameter(i, r.params[i]);
    fReplay->SetLineColor(kRed); fReplay->SetLineWidth(2);

    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();
    bool browseErrBars = !decayErrBars_ || decayErrBars_->IsOn();
    hDecay->Draw(browseErrBars ? "hist E" : "hist");
    fReplay->Draw("same");

    // BG component
    TF1* fBG = nullptr;
    if (sig != 4) {
        fBG = BuildDecayTF1("fCacheBrowseBG", 4, bgType, xlo, xhi);
        int nBGr = DecayBGNpar(bgType);
        for (int i = 0; i < nBGr; i++) fBG->SetParameter(i, fReplay->GetParameter(nSig + i));
        fBG->SetLineColor(kGreen+2); fBG->SetLineWidth(2); fBG->SetLineStyle(2);
        fBG->Draw("same");
    }

    // Legend
    TLegend* leg = new TLegend(0.38, 0.70, 0.62, 0.98);
    leg->SetBorderSize(1); leg->SetTextSize(0.016);
    leg->SetFillColor(kWhite); leg->SetFillStyle(1001);
    if (sig != 4 && !r.params.empty())
        leg->AddEntry((TObject*)nullptr,
            ("A = " + NNDCFormat(r.params[0], r.errors.empty() ? 0.0 : r.errors[0])).c_str(), "");
    if (sig >= 1 && r.params.size() > 1)
        leg->AddEntry((TObject*)nullptr,
            ("T_{1/2}(P) = " + NNDCFormat(r.params[1],
                r.errors.size() > 1 ? r.errors[1] : 0.0) + " ms").c_str(), "");
    if (sig >= 2 && r.params.size() > 2)
        leg->AddEntry((TObject*)nullptr,
            ("T_{1/2}(D) = " + NNDCFormat(r.params[2],
                r.errors.size() > 2 ? r.errors[2] : 0.0) + " ms").c_str(), "");
    if (sig >= 3 && r.params.size() > 3)
        leg->AddEntry((TObject*)nullptr,
            ("T_{1/2}(G) = " + NNDCFormat(r.params[3],
                r.errors.size() > 3 ? r.errors[3] : 0.0) + " ms").c_str(), "");
    {
        int bgP = nSig;
        if (bgType == 1 && (int)r.params.size() > bgP)
            leg->AddEntry((TObject*)nullptr,
                ("BG = " + NNDCFormat(r.params[bgP],
                    r.errors.size() > (size_t)bgP ? r.errors[bgP] : 0.0)).c_str(), "");
        else if (bgType == 2 && (int)r.params.size() > bgP+2)
            leg->AddEntry((TObject*)nullptr,
                ("BG_{flat} = " + NNDCFormat(r.params[bgP],
                    r.errors.size() > (size_t)bgP ? r.errors[bgP] : 0.0)).c_str(), "");
        else if (bgType == 3 && (int)r.params.size() > bgP+1)
            leg->AddEntry((TObject*)nullptr,
                ("T_{1/2}(BG) = " + NNDCFormat(r.params[bgP+1],
                    r.errors.size() > (size_t)bgP+1 ? r.errors[bgP+1] : 0.0) + " ms").c_str(), "");
    }
    if (r.chi2ndf >= 0)
        leg->AddEntry((TObject*)nullptr, Form("#chi^{2}/NDF=%.3f", r.chi2ndf), "");
    leg->AddEntry(hDecay, "Data", "l");
    leg->AddEntry(fReplay, "Fit", "l");
    if (fBG) leg->AddEntry(fBG, "BG", "l");
    leg->Draw();
    c->Modified(); c->Update();
    gSystem->ProcessEvents();
    SetStatus(Form("Cache browser: E=%.2f keV  chi2/ndf=%.3f", E, r.chi2ndf));
}

void GammaFitGUI::OnDecayErrBarsToggled()
{
    // Replay whichever peak is currently selected so the canvas updates immediately.
    if (decayPeakList_) {
        Int_t sel = decayPeakList_->GetSelected();
        if (sel >= 1) { OnDecayPeakSelected(sel); return; }
    }
    // No peak selected — just re-preview if possible.
    if (inputFile_ && !decayTh2Name_.empty()) OnPreviewDecay();
}

void GammaFitGUI::OnLoadDecayCache()
{
    if (decayTh2Name_.empty()) {
        AppendLog("[Decay] Run 'Refresh' first to set the active TH2");
        return;
    }
    LoadDecayFitCache();
    if (decayPeakList_) {
        Int_t sel = decayPeakList_->GetSelected();
        if (sel >= 1) OnDecayPeakSelected(sel);
    }
    SetStatus("Decay cache loaded: " + decayTh2Name_);
}

void GammaFitGUI::OnSaveDecayCache()
{
    if (decayTh2Name_.empty()) {
        AppendLog("[Decay] Run 'Refresh' first to set the active TH2");
        return;
    }
    if (decayFitStore_.empty()) {
        AppendLog("[Decay] No decay fits in memory to save");
        return;
    }
    SaveDecayFitCache();
    std::string path = DecayCacheFileFor();
    AppendLog(Form("[Decay] Saved %d fit(s) -> %s",
                   (int)decayFitStore_.size(), path.c_str()));
    SetStatus("Decay cache saved: " + decayTh2Name_);
}


void GammaFitGUI::OnDecayPeakSelected(Int_t id)
{
    if (id < 1 || (size_t)id > decayPeakKeys_.size()) return;
    const std::string& cacheName = (id <= (Int_t)decayPeakCacheNames_.size())
                                   ? decayPeakCacheNames_[id - 1]
                                   : decayGammaProjName_;
    if (cacheName.empty()) return;

    // Populate label from gamma projection cache
    FitDatabase fdb;
    fdb.Load(CacheFileFor(cacheName));
    const auto& entries = fdb.GetEntries();
    auto it = entries.find(decayPeakKeys_[id - 1]);
    if (it != entries.end()) {
        if (decayLabelEntry_)
            decayLabelEntry_->SetText(it->second.label.c_str());
    }

    double E = decayPeakEs_[id - 1];

    // Restore stored decay fit parameters
    auto fitIt = decayFitStore_.find(E);
    if (fitIt != decayFitStore_.end()) {
        const DecayFitResult& r = fitIt->second;
        if (decayModelCombo_) decayModelCombo_->Select(r.model, kFALSE);
        if (decayBGTypeCombo_) decayBGTypeCombo_->Select(r.bgType > 0 ? r.bgType : 1, kFALSE);
        OnDecayBGTypeChanged(0);
        int sig = DecaySignalType(r.model);
        if (sig >= 1 && r.params.size() >= 2 && decayThalfP_)
            decayThalfP_->SetNumber(r.params[1]);
        if (sig >= 2 && r.params.size() >= 3 && decayThalfD_)
            decayThalfD_->SetNumber(r.params[2]);
        if (sig >= 3 && r.params.size() >= 4 && decayThalfG_)
            decayThalfG_->SetNumber(r.params[3]);
        if (!r.label.empty() && decayLabelEntry_)
            decayLabelEntry_->SetText(r.label.c_str());
        if (r.Nsig > 0 && decaySigRangeEntry_)
            decaySigRangeEntry_->SetNumber(r.Nsig);
        if (decaySigLoEntry_)
            decaySigLoEntry_->SetNumber(r.NsigLo > 0 ? r.NsigLo : r.Nsig);
        if (r.rebin > 1 && decayRebinEntry_)
            decayRebinEntry_->SetNumber(r.rebin);
        if (decayFitMethod_) decayFitMethod_->Select(r.fitMethod, kFALSE);
        if (decayFitFullRange_) decayFitFullRange_->SetOn(r.fullRange, kFALSE);
        if (r.fitLo > 0 && decayFitLo_) decayFitLo_->SetNumber(r.fitLo);
        if (r.fitHi > 0 && decayFitHi_) decayFitHi_->SetNumber(r.fitHi);
        // Restore T1/2 bounds
        if (!r.paramLo.empty() && !r.paramHi.empty()) {
            int sig = DecaySignalType(r.model);
            if (sig >= 1) {
                if (decayThalfPLo_) decayThalfPLo_->SetNumber(r.paramLo.size() > 1 ? r.paramLo[1] : 1e-6);
                if (decayThalfPHi_) decayThalfPHi_->SetNumber(r.paramHi.size() > 1 ? r.paramHi[1] : 1e12);
            }
            if (sig >= 2) {
                if (decayThalfDLo_) decayThalfDLo_->SetNumber(r.paramLo.size() > 2 ? r.paramLo[2] : 1e-6);
                if (decayThalfDHi_) decayThalfDHi_->SetNumber(r.paramHi.size() > 2 ? r.paramHi[2] : 1e12);
            }
            if (sig >= 3) {
                if (decayThalfGLo_) decayThalfGLo_->SetNumber(r.paramLo.size() > 3 ? r.paramLo[3] : 1e-6);
                if (decayThalfGHi_) decayThalfGHi_->SetNumber(r.paramHi.size() > 3 ? r.paramHi[3] : 1e12);
            }
        }
        // Restore extra cuts
        extraDecayCuts_ = r.extraCuts;
        if (decayExtraCutList_) {
            decayExtraCutList_->RemoveAll();
            for (int i = 0; i < (int)extraDecayCuts_.size(); i++)
                decayExtraCutList_->AddEntry(
                    Form("%.2f - %.2f keV",
                         extraDecayCuts_[i].first, extraDecayCuts_[i].second), i + 1);
            decayExtraCutList_->MapSubwindows();
            decayExtraCutList_->Layout();
        }

        // ── Replay cached decay fit on canvas ─────────────────────────────────
        if (inputFile_ && !r.histName.empty()) {
            TH2* h2 = (TH2*)inputFile_->Get(r.histName.c_str());
            if (h2) {
                int axisId = decayGammaAxisCombo_ ? decayGammaAxisCombo_->GetSelected() : 1;
                TH1* hDecay = nullptr;
                if (axisId == 1) {
                    int b1 = h2->GetXaxis()->FindBin(r.eMin);
                    int b2 = h2->GetXaxis()->FindBin(r.eMax);
                    hDecay = h2->ProjectionY(Form("hDecayReplay_%.1f", E), b1, b2);
                } else {
                    int b1 = h2->GetYaxis()->FindBin(r.eMin);
                    int b2 = h2->GetYaxis()->FindBin(r.eMax);
                    hDecay = h2->ProjectionX(Form("hDecayReplay_%.1f", E), b1, b2);
                }

                if (hDecay && hDecay->GetEntries() > 0) {
                    hDecay->SetDirectory(nullptr);
                    // Replay extra cuts from cached result
                    if (!r.extraCuts.empty())
                        ApplyExtraCuts(h2, axisId, hDecay, r.extraCuts);
                    int rebin = decayRebinEntry_ ? (int)decayRebinEntry_->GetNumber() : 1;
                    if (rebin > 1) hDecay->Rebin(rebin);

                    double tLo = std::max(0.0, hDecay->GetXaxis()->GetXmin());
                    double tHi = hDecay->GetXaxis()->GetXmax();
                    double xlo = r.fullRange ? tLo : std::max(tLo, r.fitLo);
                    double xhi = r.fullRange ? tHi : (r.fitHi > 0 ? std::min(tHi, r.fitHi) : tHi);
                    hDecay->GetXaxis()->SetRangeUser(xlo, xhi);
                    {
                        std::string ttl = Form("E=%.2f keV  [%.2f, %.2f keV]",
                                               E, r.eMin, r.eMax);
                        if (!r.extraCuts.empty())
                            ttl += Form("  +%d extra cut(s)", (int)r.extraCuts.size());
                        hDecay->SetTitle(ttl.c_str());
                    }
                    hDecay->GetXaxis()->SetTitle("Time (ms)");
                    hDecay->GetYaxis()->SetTitle("Counts");
                    hDecay->SetLineColor(kBlack);
                    hDecay->SetMarkerSize(0);

                    int bgType  = (r.bgType >= 1 && r.bgType <= 3) ? r.bgType : 1;
                    int nSig    = DecaySigNpar(sig);
                    TF1* fStored = BuildDecayTF1("fDecayReplay", r.model, bgType, xlo, xhi);
                    if (fStored && (int)r.params.size() >= fStored->GetNpar()) {
                        for (int i = 0; i < fStored->GetNpar(); i++)
                            fStored->SetParameter(i, r.params[i]);
                        fStored->SetLineColor(kRed);
                        fStored->SetLineWidth(2);

                        TCanvas* c = canvas_->GetCanvas();
                        c->Clear(); c->cd();
                        bool replayErrBars = !decayErrBars_ || decayErrBars_->IsOn();
                        hDecay->Draw(replayErrBars ? "hist E" : "hist");
                        fStored->Draw("same");

                        // BG component from clean BGonly TF1
                        if (sig != 4) {
                            int replayBGType = (r.bgType >= 1 && r.bgType <= 3) ? r.bgType : 1;
                            TF1* fBG = BuildDecayTF1("fDecayReplayBG", 4, replayBGType, tLo, tHi);
                            int nBGr = DecayBGNpar(replayBGType);
                            for (int i = 0; i < nBGr; i++)
                                fBG->SetParameter(i, fStored->GetParameter(nSig + i));
                            fBG->SetLineColor(kGreen+2);
                            fBG->SetLineWidth(2);
                            fBG->SetLineStyle(2);
                            fBG->Draw("same");

                            TLegend* leg = new TLegend(0.38, 0.70, 0.62, 0.98);
                            leg->SetBorderSize(1);
                            leg->SetTextSize(0.016);
                            leg->SetFillColor(kWhite); leg->SetFillStyle(1001);
                            // Activity and half-lives from cached params
                            if (sig != 4 && r.params.size() > 0)
                                leg->AddEntry((TObject*)nullptr,
                                    ("A = " + NNDCFormat(r.params[0],
                                         r.errors.size() > 0 ? r.errors[0] : 0.0)).c_str(), "");
                            if (sig >= 1 && r.params.size() > 1)
                                leg->AddEntry((TObject*)nullptr,
                                    ("T_{1/2}(P) = " + NNDCFormat(r.params[1],
                                         r.errors.size() > 1 ? r.errors[1] : 0.0) + " ms").c_str(), "");
                            if (sig >= 2 && r.params.size() > 2)
                                leg->AddEntry((TObject*)nullptr,
                                    ("T_{1/2}(D) = " + NNDCFormat(r.params[2],
                                         r.errors.size() > 2 ? r.errors[2] : 0.0) + " ms").c_str(), "");
                            if (sig >= 3 && r.params.size() > 3)
                                leg->AddEntry((TObject*)nullptr,
                                    ("T_{1/2}(G) = " + NNDCFormat(r.params[3],
                                         r.errors.size() > 3 ? r.errors[3] : 0.0) + " ms").c_str(), "");
                            {   // BG from cached params
                                int bgP = nSig;
                                if (bgType == 1 && (int)r.params.size() > bgP)
                                    leg->AddEntry((TObject*)nullptr,
                                        ("BG = " + NNDCFormat(r.params[bgP],
                                             r.errors.size() > (size_t)bgP ? r.errors[bgP] : 0.0)).c_str(), "");
                                else if (bgType == 2 && (int)r.params.size() > bgP+2)
                                    leg->AddEntry((TObject*)nullptr,
                                        ("BG_{flat} = " + NNDCFormat(r.params[bgP],
                                             r.errors.size() > (size_t)bgP ? r.errors[bgP] : 0.0)).c_str(), "");
                                else if (bgType == 3 && (int)r.params.size() > bgP+1)
                                    leg->AddEntry((TObject*)nullptr,
                                        ("T_{1/2}(BG) = " + NNDCFormat(r.params[bgP+1],
                                             r.errors.size() > (size_t)bgP+1 ? r.errors[bgP+1] : 0.0) + " ms").c_str(), "");
                            }
                            if (r.chi2ndf >= 0)
                                leg->AddEntry((TObject*)nullptr,
                                    Form("#chi^{2}/NDF=%.3f", r.chi2ndf), "");
                            leg->AddEntry(hDecay,  "Data", "l");
                            leg->AddEntry(fStored, "Fit",  "l");
                            leg->AddEntry(fBG,     "BG",   "l");
                            leg->Draw();
                        }

                        c->Modified(); c->Update();
                        gSystem->ProcessEvents();
                    }

                    // Update result view
                    decayResultView_->Clear();
                    {
                        TAxis* gax2 = (axisId == 1) ? h2->GetXaxis() : h2->GetYaxis();
                        int cb1 = gax2->FindBin(r.eMin), cb2 = gax2->FindBin(r.eMax);
                        decayResultView_->AddLine(Form(
                            "Cached fit: E=%.2f keV  [%.2f, %.2f keV]  bins [%d, %d]",
                            E, r.eMin, r.eMax, cb1, cb2));
                        for (size_t ci = 0; ci < r.extraCuts.size(); ci++) {
                            int eb1 = gax2->FindBin(r.extraCuts[ci].first);
                            int eb2 = gax2->FindBin(r.extraCuts[ci].second);
                            decayResultView_->AddLine(Form(
                                "  ExtraCut[%zu]: [%.2f, %.2f keV]  bins [%d, %d]",
                                ci + 1, r.extraCuts[ci].first, r.extraCuts[ci].second, eb1, eb2));
                        }
                    }
                    decayResultView_->AddLine(
                        Form("Model: %d  bgType: %d  chi2/ndf=%.4g  status=%d",
                             r.model, bgType, r.chi2ndf, r.status));
                    if (fStored) {
                        for (int i = 0; i < fStored->GetNpar(); i++)
                            decayResultView_->AddLine(
                                Form("  %-14s = %11.5g +/- %.4g",
                                     fStored->GetParName(i), r.params[i],
                                     (i < (int)r.errors.size()) ? r.errors[i] : 0.0));
                    }
                    decayResultView_->ShowBottom();
                    SetStatus(Form("Cached decay fit: E=%.2f keV  chi2/ndf=%.3f", E, r.chi2ndf));
                } else {
                    delete hDecay;
                }
            }
        }
        OnDecayCutChanged();
        return;  // done
    }

    // ── No stored fit  -  default sigma window to where Gaussian meets background ──
    // Nsig = sqrt(2 * ln(A / bg_at_centroid)), clamped to [1.5, 5.0]
    double NsigDef = 2.5;
    if (it != entries.end()) {
        const auto& fe = it->second;
        FitLayout flay = DetectLayout((int)fe.params.size());
        if (flay.valid()) {
            // Find which peak in this (possibly multi-peak) entry is closest to E
            int bi = 0;
            double bestDiff = 1e9;
            for (int pi = 0; pi < flay.n; pi++) {
                double diff = std::abs(fe.params[3*pi+1] - E);
                if (diff < bestDiff) { bestDiff = diff; bi = pi; }
            }
            double A   = fe.params[3*bi];
            double sig = (3*bi+2 < (int)fe.params.size()) ? std::abs(fe.params[3*bi+2]) : 0.0;
            int bgB    = flay.bgBase();
            double bg  = fe.params[bgB] + fe.params[bgB+1]*E;

            // Compton step: erfc(0) = 1, so step contributes S_i to background at centroid
            if (flay.comptonStep) {
                int sIdx = StepParIdx(flay.n, flay.quadBG, bi);
                if (sIdx < (int)fe.params.size())
                    bg += fe.params[sIdx];
            }

            // Exponential tail: add tail value at centroid to the peak amplitude
            if (flay.expoTail && sig > 0.0) {
                int tIdx = TailAmplParIdx(flay.n, flay.quadBG, flay.comptonStep, bi);
                int bIdx = TailSlopeParIdx(flay.n, flay.quadBG, flay.comptonStep, bi);
                if (tIdx < (int)fe.params.size() && bIdx < (int)fe.params.size()) {
                    double T    = fe.params[tIdx];
                    double beta = fe.params[bIdx];
                    if (beta > 0.0)
                        A += T/2.0 * std::exp(sig*sig/(2.0*beta*beta))
                               * TMath::Erfc(sig / (beta * 1.41421356));
                }
            }

            if (A > 0.0 && bg > 0.0 && A > bg) {
                double nsig = std::sqrt(2.0 * std::log(A / bg));
                NsigDef = std::min(5.0, std::max(1.5, nsig));
            }
        }
    }
    if (decaySigRangeEntry_) decaySigRangeEntry_->SetNumber(NsigDef);
    if (decaySigLoEntry_)    decaySigLoEntry_->SetNumber(NsigDef);
    OnDecayCutChanged();
    OnDecayPreviewGammaPeak();
}

// ─────────────────────────────────────────────────────────────────────────────
// OnDecayCutChanged  -  update live cut-edge label (keV + bins)
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnDecayCutChanged()
{
    if (!decayBinInfoLabel_) return;

    Int_t peakSel = decayPeakList_ ? decayPeakList_->GetSelected() : -1;
    if (peakSel < 1 || (size_t)peakSel > decayPeakEs_.size()) {
        decayBinInfoLabel_->SetText("(select a peak)");
        return;
    }

    double E      = decayPeakEs_[peakSel - 1];
    double sig    = decayPeakSigs_[peakSel - 1];
    double NsigHi = decaySigRangeEntry_ ? decaySigRangeEntry_->GetNumber() : 1.0;
    double NsigLo = decaySigLoEntry_    ? decaySigLoEntry_->GetNumber()    : NsigHi;
    double eMin   = E - NsigLo * sig;
    double eMax   = E + NsigHi * sig;

    if (inputFile_ && !decayTh2Name_.empty()) {
        TH2* h2 = (TH2*)inputFile_->Get(decayTh2Name_.c_str());
        if (h2) {
            int    axisId = decayGammaAxisCombo_ ? decayGammaAxisCombo_->GetSelected() : 1;
            TAxis* ax     = (axisId == 1) ? h2->GetXaxis() : h2->GetYaxis();
            if (ax) {
                int b1 = ax->FindBin(eMin);
                int b2 = ax->FindBin(eMax);
                decayBinInfoLabel_->SetText(
                    Form("Lo: %.2f keV [bin %d]   Hi: %.2f keV [bin %d]", eMin, b1, eMax, b2));
                return;
            }
        }
    }

    // No TH2 available — show keV only
    decayBinInfoLabel_->SetText(
        Form("Lo: %.2f keV   Hi: %.2f keV", eMin, eMax));
}

// ─────────────────────────────────────────────────────────────────────────────
// OnDecayPreviewGammaPeak  -  show gamma spectrum zoomed to selected peak with cut
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnDecayPreviewGammaPeak()
{
    Int_t peakSel = decayPeakList_ ? decayPeakList_->GetSelected() : -1;
    if (peakSel < 1 || (size_t)peakSel > decayPeakEs_.size()) {
        AppendLog("[Decay] Select a peak first"); return;
    }
    if (!inputFile_ || decayTh2Name_.empty()) {
        AppendLog("[Decay] No TH2 loaded  -  run Refresh first"); return;
    }

    double E      = decayPeakEs_[peakSel - 1];
    double sig    = decayPeakSigs_[peakSel - 1];
    double NsigHi = decaySigRangeEntry_ ? decaySigRangeEntry_->GetNumber() : 1.0;
    double NsigLo = decaySigLoEntry_    ? decaySigLoEntry_->GetNumber()    : NsigHi;
    double eMin   = E - NsigLo * sig;
    double eMax   = E + NsigHi * sig;

    TH2* h2 = (TH2*)inputFile_->Get(decayTh2Name_.c_str());
    if (!h2) { AppendLog("[Decay] TH2 not found"); return; }

    int axisId = decayGammaAxisCombo_ ? decayGammaAxisCombo_->GetSelected() : 1;
    TH1* hGamma = (axisId == 1)
        ? h2->ProjectionX(Form("hGamPreview_%.1f", E))
        : h2->ProjectionY(Form("hGamPreview_%.1f", E));
    if (!hGamma || hGamma->GetEntries() == 0) { delete hGamma; return; }
    hGamma->SetDirectory(nullptr);

    hGamma->GetXaxis()->SetRangeUser(E - 6.0*sig, E + 6.0*sig);
    hGamma->GetXaxis()->SetTitle("Energy (keV)");
    hGamma->GetYaxis()->SetTitle("Counts");
    hGamma->SetLineColor(kBlack);
    hGamma->SetMarkerSize(0);

    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();
    hGamma->Draw("hist");
    c->Modified(); c->Update();

    double ylo = gPad->GetUymin();
    double yhi = gPad->GetUymax();

    TLine* lineLo = new TLine(eMin, ylo, eMin, yhi);
    lineLo->SetLineColor(kRed); lineLo->SetLineWidth(2); lineLo->SetLineStyle(2);
    lineLo->SetBit(kCanDelete); lineLo->Draw();

    TLine* lineHi = new TLine(eMax, ylo, eMax, yhi);
    lineHi->SetLineColor(kRed); lineHi->SetLineWidth(2); lineHi->SetLineStyle(2);
    lineHi->SetBit(kCanDelete); lineHi->Draw();

    // Draw extra gamma cuts in blue
    for (const auto& cut : extraDecayCuts_) {
        TLine* lLo = new TLine(cut.first,  ylo, cut.first,  yhi);
        lLo->SetLineColor(kBlue); lLo->SetLineWidth(2); lLo->SetLineStyle(2);
        lLo->SetBit(kCanDelete); lLo->Draw();
        TLine* lHi = new TLine(cut.second, ylo, cut.second, yhi);
        lHi->SetLineColor(kBlue); lHi->SetLineWidth(2); lHi->SetLineStyle(2);
        lHi->SetBit(kCanDelete); lHi->Draw();
    }

    c->Modified(); c->Update();
    gSystem->ProcessEvents();
    SetStatus(Form("Cut preview: %.4f keV  [%.4f, %.4f]  (lo=%.3f hi=%.3f sig)  +%d extra cut(s)",
                   E, eMin, eMax, NsigLo, NsigHi, (int)extraDecayCuts_.size()));
}

void GammaFitGUI::OnRefreshDecayPeaks()
{
    decayPeakEs_.clear();
    decayPeakSigs_.clear();
    decayPeakKeys_.clear();
    decayPeakCacheNames_.clear();
    decayPeakList_->RemoveAll();

    // Determine which cache to load peaks from.
    // Priority: explicitly selected cache combo entry > TH2-derived projection name.
    std::string peakCacheName;
    TGLBEntry* cacheSelE = decayCacheCombo_ ? decayCacheCombo_->GetSelectedEntry() : nullptr;
    if (cacheSelE) {
        peakCacheName = cacheSelE->GetTitle();
    } else {
        // Fall back to auto-derive from selected TH2
        TGLBEntry* th2SelE = decayTh2Combo_->GetSelectedEntry();
        if (!th2SelE) { AppendLog("[Decay] No TH2 or cache selected"); return; }
        decayTh2Name_ = th2SelE->GetTitle();
        int axisId = decayGammaAxisCombo_->GetSelected();
        peakCacheName = decayTh2Name_ + (axisId == 1 ? "_px" : "_py");
    }

    // Update active TH2 name if it's not already set (needed for projection later)
    if (decayTh2Name_.empty()) {
        TGLBEntry* th2SelE = decayTh2Combo_->GetSelectedEntry();
        if (th2SelE) decayTh2Name_ = th2SelE->GetTitle();
    }

    decayGammaProjName_ = peakCacheName;

    FitDatabase fdb;
    fdb.Load(CacheFileFor(peakCacheName));

    struct DPEntry { double E, sig; std::string key, cache, lbl; };
    std::vector<DPEntry> dpEntries;
    for (const auto& kv : fdb.GetEntries()) {
        const FitEntry& e = kv.second;
        if (e.params.size() < 5) continue;
        FitLayout lay = DetectLayout((int)e.params.size());
        if (!lay.valid()) continue;
        for (int i = 0; i < lay.n; i++) {
            double Ep  = e.params[3*i + 1];
            double sig = std::abs(e.params[3*i + 2]);
            if (sig <= 0 || sig > 100 || Ep <= 0) continue;
            std::string lbl = Form("%.2f keV  (sig=%.3f)", Ep, sig);
            if (!e.label.empty()) lbl = e.label + "  " + lbl;
            if (!e.classification.empty()) lbl += "  [" + e.classification + "]";
            if (e.needsRefit) lbl = "[R] " + lbl;
            dpEntries.push_back({Ep, sig, kv.first, peakCacheName, lbl});
        }
    }
    std::sort(dpEntries.begin(), dpEntries.end(),
              [](const DPEntry& a, const DPEntry& b){ return a.E < b.E; });

    int listIdx = 1;
    for (const auto& dp : dpEntries) {
        decayPeakEs_.push_back(dp.E);
        decayPeakSigs_.push_back(dp.sig);
        decayPeakKeys_.push_back(dp.key);
        decayPeakCacheNames_.push_back(dp.cache);
        decayPeakList_->AddEntry(dp.lbl.c_str(), listIdx++);
    }
    decayPeakList_->MapSubwindows();
    decayPeakList_->Layout();

    // Mirror to RK4 peak list
    if (rkPeakList_) {
        rkPeakList_->RemoveAll();
        int rkIdx = 1;
        for (const auto& dp : dpEntries)
            rkPeakList_->AddEntry(dp.lbl.c_str(), rkIdx++);
        rkPeakList_->MapSubwindows();
        rkPeakList_->Layout();
    }

    LoadDecayFitCache();
    AppendLog(Form("[Decay] %d peaks loaded from cache: %s",
                   (int)decayPeakEs_.size(), peakCacheName.c_str()));
}

// ─────────────────────────────────────────────────────────────────────────────
// Extra decay cut management
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnAddExtraDecayCut()
{
    if (!decayExtraCutLo_ || !decayExtraCutHi_ || !decayExtraCutList_) return;
    double lo = decayExtraCutLo_->GetNumber();
    double hi = decayExtraCutHi_->GetNumber();
    if (lo >= hi) { AppendLog("[Decay] Extra cut: Lo must be < Hi"); return; }
    extraDecayCuts_.push_back({lo, hi});
    int id = (int)extraDecayCuts_.size();
    decayExtraCutList_->AddEntry(Form("%.2f - %.2f keV", lo, hi), id);
    decayExtraCutList_->MapSubwindows();
    decayExtraCutList_->Layout();
}

void GammaFitGUI::OnRemoveExtraDecayCut()
{
    if (!decayExtraCutList_) return;
    Int_t sel = decayExtraCutList_->GetSelected();
    if (sel < 1 || (size_t)sel > extraDecayCuts_.size()) return;
    extraDecayCuts_.erase(extraDecayCuts_.begin() + sel - 1);
    // Rebuild list
    decayExtraCutList_->RemoveAll();
    for (int i = 0; i < (int)extraDecayCuts_.size(); i++)
        decayExtraCutList_->AddEntry(
            Form("%.2f - %.2f keV", extraDecayCuts_[i].first, extraDecayCuts_[i].second), i + 1);
    decayExtraCutList_->MapSubwindows();
    decayExtraCutList_->Layout();
}

void GammaFitGUI::OnClearExtraDecayCuts()
{
    extraDecayCuts_.clear();
    if (decayExtraCutList_) {
        decayExtraCutList_->RemoveAll();
        decayExtraCutList_->MapSubwindows();
        decayExtraCutList_->Layout();
    }
}

// Sum extra gamma cuts into hBase (in-place Add).
// axisId: 1 = gamma on X axis, 2 = gamma on Y axis.
static void ApplyExtraCuts(TH2* h2, int axisId, TH1* hBase,
                            const std::vector<std::pair<double,double>>& cuts)
{
    for (size_t i = 0; i < cuts.size(); i++) {
        double lo = cuts[i].first, hi = cuts[i].second;
        TH1* hExtra = nullptr;
        if (axisId == 1) {
            int b1 = h2->GetXaxis()->FindBin(lo);
            int b2 = h2->GetXaxis()->FindBin(hi);
            hExtra = h2->ProjectionY(Form("hExtraCut_%zu", i), b1, b2);
        } else {
            int b1 = h2->GetYaxis()->FindBin(lo);
            int b2 = h2->GetYaxis()->FindBin(hi);
            hExtra = h2->ProjectionX(Form("hExtraCut_%zu", i), b1, b2);
        }
        if (hExtra) {
            hBase->Add(hExtra);
            delete hExtra;
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnFitDecay()
{
    if (!inputFile_) { AppendLog("[Decay] No ROOT file loaded"); return; }
    if (decayTh2Name_.empty()) { AppendLog("[Decay] Run 'Refresh' first"); return; }

    Int_t peakSel = decayPeakList_->GetSelected();
    if (peakSel < 1 || (size_t)peakSel > decayPeakEs_.size()) {
        AppendLog("[Decay] Select a peak first"); return;
    }
    double E      = decayPeakEs_[peakSel - 1];
    double sig    = decayPeakSigs_[peakSel - 1];
    double NsigHi = decaySigRangeEntry_ ? decaySigRangeEntry_->GetNumber() : 1.0;
    double NsigLo = decaySigLoEntry_    ? decaySigLoEntry_->GetNumber()    : NsigHi;

    TH2* h2 = (TH2*)inputFile_->Get(decayTh2Name_.c_str());
    if (!h2) { AppendLog("[Decay] TH2 not found in file"); return; }

    int axisId = decayGammaAxisCombo_->GetSelected();
    TH1* hDecay = nullptr;
    double eMin = E - NsigLo * sig;
    double eMax = E + NsigHi * sig;

    if (axisId == 1) {
        // Gamma on X -> project onto Y axis
        int b1 = h2->GetXaxis()->FindBin(eMin);
        int b2 = h2->GetXaxis()->FindBin(eMax);
        hDecay = h2->ProjectionY(Form("hDecay_%.1f", E), b1, b2);
    } else {
        // Gamma on Y -> project onto X axis
        int b1 = h2->GetYaxis()->FindBin(eMin);
        int b2 = h2->GetYaxis()->FindBin(eMax);
        hDecay = h2->ProjectionX(Form("hDecay_%.1f", E), b1, b2);
    }

    if (!hDecay || hDecay->GetEntries() == 0) {
        AppendLog("[Decay] Projection is empty  -  check TH2 and peak selection");
        return;
    }
    hDecay->SetDirectory(nullptr);

    // Add any extra gamma cuts
    if (!extraDecayCuts_.empty())
        ApplyExtraCuts(h2, axisId, hDecay, extraDecayCuts_);

    // Rebin
    int rebin = decayRebinEntry_ ? (int)decayRebinEntry_->GetNumber() : 1;
    if (rebin > 1) hDecay->Rebin(rebin);

    // Restrict to t > 0
    double tAxisMin = std::max(0.0, hDecay->GetXaxis()->GetXmin());
    hDecay->GetXaxis()->SetRangeUser(tAxisMin, hDecay->GetXaxis()->GetXmax());

    TCanvas* c = canvas_->GetCanvas();
    c->Clear();
    c->cd();
    {
        std::string ttl = Form("E=%.2f keV  [%.2f, %.2f keV]  -%.2f/+%.2f#sigma",
                               E, eMin, eMax, NsigLo, NsigHi);
        if (!extraDecayCuts_.empty())
            ttl += Form("  +%d extra cut(s)", (int)extraDecayCuts_.size());
        hDecay->SetTitle(ttl.c_str());
    }
    hDecay->GetXaxis()->SetTitle("Time (ms)");
    if (rebin > 1) {
        double binW = hDecay->GetBinWidth(1);
        hDecay->GetYaxis()->SetTitle(Form("Counts / (%.4g ms)", binW));
    } else {
        hDecay->GetYaxis()->SetTitle("Counts");
    }
    hDecay->SetLineColor(kBlack);
    hDecay->SetMarkerSize(0);
    bool cutsErrBars = !decayErrBars_ || decayErrBars_->IsOn();
    Int_t savedStat = gStyle->GetOptStat();
    gStyle->SetOptStat(111111);  // entries, mean, RMS, underflow, overflow, integral
    hDecay->Draw(cutsErrBars ? "hist E" : "hist");
    c->Modified(); c->Update();
    gStyle->SetOptStat(savedStat);
    gSystem->ProcessEvents();

    // ── Build fit function ────────────────────────────────────────────────────
    int modelId = decayModelCombo_ ? decayModelCombo_->GetSelected() : 1;
    int bgType  = decayBGTypeCombo_ ? decayBGTypeCombo_->GetSelected() : 1;
    if (bgType < 1 || bgType > 3) bgType = 1;
    int sigType = DecaySignalType(modelId);
    int nSig    = DecaySigNpar(sigType);

    bool fullRange = !decayFitFullRange_ || decayFitFullRange_->IsOn();
    double xlo = fullRange ? tAxisMin : std::max(0.0, decayFitLo_->GetNumber());
    double xhi = fullRange ? hDecay->GetXaxis()->GetXmax() : decayFitHi_->GetNumber();
    if (xlo >= xhi) { xlo = tAxisMin; xhi = hDecay->GetXaxis()->GetXmax(); }

    // Restrict canvas view to the fit range so the fitted region fills the pad
    hDecay->GetXaxis()->SetRangeUser(xlo, xhi);

    TF1* fDecay = BuildDecayTF1("fDecay", modelId, bgType, xlo, xhi);

    // Seed signal parameters
    double initMax = std::max(1.0, hDecay->GetMaximum() - hDecay->GetMinimum());
    double initMin = std::max(0.0, hDecay->GetMinimum());
    double TP  = decayThalfP_  ? decayThalfP_->GetNumber()  : (xhi-xlo)/2.0;
    double TD  = decayThalfD_  ? decayThalfD_->GetNumber()  : (xhi-xlo)/4.0;
    double TG  = decayThalfG_  ? decayThalfG_->GetNumber()  : (xhi-xlo)/8.0;

    // Read user-specified bounds for each T1/2 parameter
    auto bndLo = [](TGNumberEntry* e, double def) { return e ? e->GetNumber() : def; };
    auto bndHi = [](TGNumberEntry* e, double def) { return e ? e->GetNumber() : def; };
    double pLoP = bndLo(decayThalfPLo_,      1e-6); double pHiP = bndHi(decayThalfPHi_,      1e12);
    double pLoD = bndLo(decayThalfDLo_,      1e-6); double pHiD = bndHi(decayThalfDHi_,      1e12);
    double pLoG = bndLo(decayThalfGLo_,      1e-6); double pHiG = bndHi(decayThalfGHi_,      1e12);
    double pLoBG= bndLo(decayThalfBGExpLo_,  1e-6); double pHiBG= bndHi(decayThalfBGExpHi_,  1e15);

    if (sigType == 1) {
        fDecay->SetParameter(0, initMax);
        fDecay->SetParameter(1, TP > 0 ? TP : (xhi-xlo)/2.0);
        fDecay->SetParLimits(1, pLoP, pHiP);
        if (decayFixP_ && decayFixP_->IsOn()) fDecay->FixParameter(1, TP);
    } else if (sigType == 2) {
        fDecay->SetParameter(0, initMax);
        fDecay->SetParameter(1, TP > 0 ? TP : (xhi-xlo)/4.0);
        fDecay->SetParameter(2, TD > 0 ? TD : (xhi-xlo)/2.0);
        fDecay->SetParLimits(1, pLoP, pHiP);
        fDecay->SetParLimits(2, pLoD, pHiD);
        if (decayFixP_ && decayFixP_->IsOn()) fDecay->FixParameter(1, TP);
        if (decayFixD_ && decayFixD_->IsOn()) fDecay->FixParameter(2, TD);
    } else if (sigType == 3) {
        fDecay->SetParameter(0, initMax);
        fDecay->SetParameter(1, TP > 0 ? TP : (xhi-xlo)/8.0);
        fDecay->SetParameter(2, TD > 0 ? TD : (xhi-xlo)/4.0);
        fDecay->SetParameter(3, TG > 0 ? TG : (xhi-xlo)/2.0);
        fDecay->SetParLimits(1, pLoP, pHiP);
        fDecay->SetParLimits(2, pLoD, pHiD);
        fDecay->SetParLimits(3, pLoG, pHiG);
        if (decayFixP_ && decayFixP_->IsOn()) fDecay->FixParameter(1, TP);
        if (decayFixD_ && decayFixD_->IsOn()) fDecay->FixParameter(2, TD);
        if (decayFixG_ && decayFixG_->IsOn()) fDecay->FixParameter(3, TG);
    }
    // Seed BG parameters
    double Tbg = decayThalfBGExp_ ? decayThalfBGExp_->GetNumber() : (xhi-xlo)*2.0;
    if (bgType == 1) {
        fDecay->SetParameter(nSig, initMin);
    } else if (bgType == 2) {
        fDecay->SetParameter(nSig,   initMin * 0.5);
        fDecay->SetParameter(nSig+1, initMin * 0.5);
        fDecay->SetParameter(nSig+2, Tbg > 0 ? Tbg : (xhi-xlo)*2.0);
        fDecay->SetParLimits(nSig+2, pLoBG, pHiBG);
        if (decayFixBGExp_ && decayFixBGExp_->IsOn()) fDecay->FixParameter(nSig+2, Tbg);
    } else {
        fDecay->SetParameter(nSig,   initMin);
        fDecay->SetParameter(nSig+1, Tbg > 0 ? Tbg : (xhi-xlo)*2.0);
        fDecay->SetParLimits(nSig+1, pLoBG, pHiBG);
        if (decayFixBGExp_ && decayFixBGExp_->IsOn()) fDecay->FixParameter(nSig+1, Tbg);
    }

    fDecay->SetLineColor(kRed);
    fDecay->SetLineWidth(2);

    int cutsMethod = decayFitMethod_ ? decayFitMethod_->GetSelected() : 2;
    std::string cutsOpts = "SR";
    if (cutsMethod == 2) cutsOpts = "SRL";
    if (cutsMethod == 3) cutsOpts = "SRE";
    if (cutsMethod == 4) cutsOpts = "SRLE";

    TFitResultPtr fitRes = hDecay->Fit(fDecay, cutsOpts.c_str(), "", xlo, xhi);

    double chi2val = fitRes.Get() ? fitRes->Chi2() : fDecay->GetChisquare();
    int    ndfVal  = fitRes.Get() ? fitRes->Ndf()  : fDecay->GetNDF();
    double chi2ndf = (ndfVal > 0) ? chi2val / ndfVal : -1.0;
    double pval    = (ndfVal > 0) ? TMath::Prob(chi2val, ndfVal) : -1.0;
    int    fitStatC = fitRes.Get() ? fitRes->Status() : -1;
    const char* statStrC = (fitStatC == 0) ? "OK" : Form("WARN(%d)", fitStatC);

    bool cutsShowRes = decayShowResid_ && decayShowResid_->IsOn();
    if (cutsShowRes) {
        DrawDecayResiduals(hDecay, fDecay, xlo, xhi);
    } else {
        c->cd();
        hDecay->Draw(cutsErrBars ? "hist E" : "hist");
        fDecay->Draw("same");

        // Build BG component from a clean BGonly TF1 (modelId=4) with copied BG params
        TF1* fBGComp = nullptr;
        if (sigType != 4) {
            fBGComp = BuildDecayTF1("fDecayBGComp", 4, bgType, xlo, xhi);
            int nBG = DecayBGNpar(bgType);
            for (int i = 0; i < nBG; i++)
                fBGComp->SetParameter(i, fDecay->GetParameter(nSig + i));
            fBGComp->SetLineColor(kGreen+2); fBGComp->SetLineWidth(2); fBGComp->SetLineStyle(2);
            fBGComp->Draw("same");
        }

        // ── Build legend ─────────────────────────────────────────────────────
        double lx = decayLegX_ ? decayLegX_->GetNumber() : 0.38;
        double ly = decayLegY_ ? decayLegY_->GetNumber() : 0.72;
        TLegend* leg = new TLegend(lx, ly, lx + 0.24, ly + 0.27);
        leg->SetBorderSize(1);
        leg->SetTextSize(0.016);
        leg->SetFillColor(kWhite); leg->SetFillStyle(1001);
        // Activity
        if (sigType != 4)
            leg->AddEntry((TObject*)nullptr,
                ("A = " + NNDCFormat(fDecay->GetParameter(0), fDecay->GetParError(0))).c_str(), "");
        // Half-lives
        if (sigType >= 1) leg->AddEntry((TObject*)nullptr,
            ("T_{1/2}(P) = " + NNDCFormat(fDecay->GetParameter(1), fDecay->GetParError(1)) + " ms").c_str(), "");
        if (sigType >= 2) leg->AddEntry((TObject*)nullptr,
            ("T_{1/2}(D) = " + NNDCFormat(fDecay->GetParameter(2), fDecay->GetParError(2)) + " ms").c_str(), "");
        if (sigType >= 3) leg->AddEntry((TObject*)nullptr,
            ("T_{1/2}(G) = " + NNDCFormat(fDecay->GetParameter(3), fDecay->GetParError(3)) + " ms").c_str(), "");
        // Background
        {
            int bgP = nSig;  // index of first BG parameter
            if (bgType == 1) {
                leg->AddEntry((TObject*)nullptr,
                    ("BG = " + NNDCFormat(fDecay->GetParameter(bgP), fDecay->GetParError(bgP))).c_str(), "");
            } else if (bgType == 2) {
                leg->AddEntry((TObject*)nullptr,
                    ("BG_{flat} = " + NNDCFormat(fDecay->GetParameter(bgP), fDecay->GetParError(bgP))).c_str(), "");
                leg->AddEntry((TObject*)nullptr,
                    ("T_{1/2}(BG) = " + NNDCFormat(fDecay->GetParameter(bgP+2), fDecay->GetParError(bgP+2)) + " ms").c_str(), "");
            } else if (bgType == 3) {
                leg->AddEntry((TObject*)nullptr,
                    ("T_{1/2}(BG) = " + NNDCFormat(fDecay->GetParameter(bgP+1), fDecay->GetParError(bgP+1)) + " ms").c_str(), "");
            }
        }
        leg->AddEntry((TObject*)nullptr,
            Form("#chi^{2}/NDF=%.3f  %s", chi2ndf, statStrC), "");
        leg->AddEntry(hDecay,  "Data", "l");
        leg->AddEntry(fDecay,  "Fit",  "l");
        if (fBGComp) leg->AddEntry(fBGComp, "BG", "l");
        leg->Draw();
        c->Modified(); c->Update();
    }

    // ── Report ────────────────────────────────────────────────────────────────
    decayResultView_->Clear();
    {
        // Compute bin numbers for main cut
        TAxis* gax = (axisId == 1) ? h2->GetXaxis() : h2->GetYaxis();
        int b1m = gax->FindBin(eMin), b2m = gax->FindBin(eMax);
        decayResultView_->AddLine(Form(
            "Peak: %.2f keV  [%.2f, %.2f keV]  bins [%d, %d]  cut: -%.2f/+%.2f sig  rebin:%d",
            E, eMin, eMax, b1m, b2m, NsigLo, NsigHi,
            decayRebinEntry_ ? (int)decayRebinEntry_->GetNumber() : 1));
        // Extra cuts
        for (size_t ci = 0; ci < extraDecayCuts_.size(); ci++) {
            int b1e = gax->FindBin(extraDecayCuts_[ci].first);
            int b2e = gax->FindBin(extraDecayCuts_[ci].second);
            decayResultView_->AddLine(Form(
                "  ExtraCut[%zu]: [%.2f, %.2f keV]  bins [%d, %d]",
                ci + 1, extraDecayCuts_[ci].first, extraDecayCuts_[ci].second, b1e, b2e));
        }
    }
    {
        const char* mn = cutsMethod==2?"Likelihood":cutsMethod==3?"Chi2+MINOS":
                          cutsMethod==4?"Likelihood+MINOS":"Chi2";
        decayResultView_->AddLine(Form("Method: %s   status: %d", mn,
                                       fitRes.Get() ? fitRes->Status() : -1));
    }
    decayResultView_->AddLine(Form("Chi2 = %.2f   NDF = %d   Chi2/NDF = %.4f", chi2val, ndfVal, chi2ndf));
    decayResultView_->AddLine(Form("p-value = %.6f%s", pval,
        pval >= 0 && pval < 0.01 ? "  *** poor fit ***" : pval < 0.05 ? "  * marginal *" : ""));
    for (int i = 0; i < fDecay->GetNpar(); i++) {
        decayResultView_->AddLine(Form("  %-16s = %11.5g +/- %.4g",
            fDecay->GetParName(i), fDecay->GetParameter(i), fDecay->GetParError(i)));
    }
    // Store fit result in memory + save to disk
    {
        DecayFitResult r;
        r.peakE          = E;
        r.model          = modelId;
        r.bgType         = bgType;
        r.chi2ndf        = chi2ndf;
        r.status         = fitRes.Get() ? fitRes->Status() : -1;
        r.histName       = decayTh2Name_;
        r.eMin           = eMin;
        r.eMax           = eMax;
        r.Nsig           = NsigHi;
        r.NsigLo         = NsigLo;
        r.rebin          = decayRebinEntry_ ? (int)decayRebinEntry_->GetNumber() : 1;
        r.fitMethod      = cutsMethod;
        r.fitLo          = xlo;
        r.fitHi          = xhi;
        r.fullRange      = fullRange;
        r.extraCuts      = extraDecayCuts_;
        // Preserve existing label/class if already set
        auto existing = decayFitStore_.find(E);
        if (existing != decayFitStore_.end()) {
            r.label          = existing->second.label;
            r.classification = existing->second.classification;
        }
        // Overwrite label from widget; class auto-populated from labelClassMap_
        if (decayLabelEntry_) {
            std::string lbl = decayLabelEntry_->GetText();
            if (!lbl.empty()) {
                r.label = lbl;
                if (labelClassMap_.count(lbl)) r.classification = labelClassMap_.at(lbl);
            }
        }
        int np = fDecay->GetNpar();
        r.params.resize(np); r.errors.resize(np);
        r.paramLo.resize(np); r.paramHi.resize(np);
        for (int i = 0; i < np; i++) {
            r.params[i] = fDecay->GetParameter(i);
            r.errors[i] = fDecay->GetParError(i);
            double plo = 0, phi = 0;
            fDecay->GetParLimits(i, plo, phi);
            r.paramLo[i] = plo; r.paramHi[i] = phi;
        }
        decayFitStore_[E] = r;
        SaveDecayFitCache();
    }

    // ── Peak count estimates ──────────────────────────────────────────────────
    if (sigType != 4) {  // signal-bearing models only
        int    npar = fDecay->GetNpar();
        double A    = fDecay->GetParameter(0);              // amplitude (counts/bin at t=0)
        double TP   = fDecay->GetParameter(1);              // parent T1/2 (ms)
        double BG   = fDecay->GetParameter(nSig);           // flat BG param
        double binW = hDecay->GetBinWidth(1);               // ms per bin
        const double ln2 = 0.6931471805599453;

        // ── (1) Counts in fit window ────────────────────────────────────────
        // ∫[xlo,xhi] (f(t) - BG) dt / binW
        // BG cancels analytically (dsignal/dBG = 0), so zero BG row/col in cov.
        double totalInteg   = fDecay->Integral(xlo, xhi);
        double windowCounts = (totalInteg - BG * (xhi - xlo)) / binW;

        double windowErr = -1.0;
        if (fitRes.Get() && fitRes->IsValid()) {
            const TMatrixDSym& cov = fitRes->GetCovarianceMatrix();
            if (cov.GetNrows() == npar) {
                std::vector<double> covMod(npar * npar, 0.0);
                for (int ii = 0; ii < npar - 1; ii++)
                    for (int jj = 0; jj < npar - 1; jj++)
                        covMod[ii * npar + jj] = cov(ii, jj);
                std::vector<double> pv(npar);
                for (int ii = 0; ii < npar; ii++) pv[ii] = fDecay->GetParameter(ii);
                double ie = fDecay->IntegralError(xlo, xhi,
                                                  pv.data(), covMod.data(), 1e-4);
                windowErr = ie / binW;
            }
        }

        // ── (2) Extrapolated total: N0 = A * tau_P / binW ────────────────────
        // For all models (Parent, Daughter, Granddaughter), integrating the
        // signal from 0 -> inf gives A * tau_P regardless of daughter half-lives.
        // This is the total number of peak counts the full decay would produce.
        double tau       = (TP > 0) ? TP / ln2 : 0.0;
        double totalN0   = (tau > 0 && binW > 0) ? A * tau / binW : -1.0;

        double totalErr = -1.0;
        if (totalN0 > 0 && fitRes.Get() && fitRes->IsValid()) {
            const TMatrixDSym& cov = fitRes->GetCovarianceMatrix();
            if (cov.GetNrows() == npar && A > 0 && TP > 0) {
                // sig^2(N0) = (tau/bw)^2*sigA^2 + (A/(ln2*bw))^2*sigTP^2
                //         + 2*(tau/bw)*(A/(ln2*bw))*cov(A,TP)
                double dA  = tau / binW;             // dN0/dA
                double dTP = A / (ln2 * binW);       // dN0/dTP
                double varN = dA*dA*cov(0,0)
                            + dTP*dTP*cov(1,1)
                            + 2.0*dA*dTP*cov(0,1);
                if (varN > 0) totalErr = std::sqrt(varN);
            }
        }

        decayResultView_->AddLine("─────────────────────────────────");
        if (windowErr > 0.0)
            decayResultView_->AddLine(
                Form("Counts in fit window:         %.0f +/- %.0f", windowCounts, windowErr));
        else
            decayResultView_->AddLine(
                Form("Counts in fit window:         %.0f", windowCounts));

        if (totalN0 > 0) {
            if (totalErr > 0.0)
                decayResultView_->AddLine(
                    Form("Total peak counts (0->inf):     %.0f +/- %.0f", totalN0, totalErr));
            else
                decayResultView_->AddLine(
                    Form("Total peak counts (0->inf):     %.0f", totalN0));
            decayResultView_->AddLine(
                Form("  [A=%.4g, T1/2=%.4g ms, tau=%.4g ms, bw=%.4g ms]",
                     A, TP, tau, binW));
        }
    }

    // Informational model note  -  does NOT write classification to gamma cache
    if (fitRes.Get() && fitRes->Status() == 0) {
        static const char* kModelCls[] = {
            "", "Parent", "Daughter(beta-)", "GDaughter(beta-)", "Background",
            "Daughter(beta-n)", "Daughter(beta-2n)", "GDaughter(beta-n)", "GDaughter(beta-2n)"
        };
        if (modelId >= 1 && modelId <= 8)
            decayResultView_->AddLine(
                Form("  -> Model: %s  bgType:%d  (Isotopes tab to set classification)",
                     kModelCls[modelId], bgType));
    }
    decayResultView_->ShowBottom();

    // Store snapshot for chi² scan popup
    delete lastDecayCutsHist_; lastDecayCutsHist_ = nullptr;
    delete lastDecayCutsTF1_;  lastDecayCutsTF1_  = nullptr;
    lastDecayCutsHist_ = (TH1*)hDecay->Clone("_lastDecayCutsH");
    lastDecayCutsHist_->SetDirectory(nullptr);
    lastDecayCutsTF1_  = (TF1*)fDecay->Clone("_lastDecayCutsF");
    lastDecayCutsXlo_  = xlo;
    lastDecayCutsXhi_  = xhi;
    if (decayChiScanBtn_) decayChiScanBtn_->SetState(kButtonUp);

    AppendLog(Form("[Decay] Fit done: peak=%.2f keV  chi2/ndf=%.3f", E, chi2ndf));
    SetStatus(Form("Decay fit: %.2f keV  chi2/ndf=%.3f", E, chi2ndf));
}

// ─────────────────────────────────────────────────────────────────────────────
// Decay fit cache helpers
// ─────────────────────────────────────────────────────────────────────────────

std::string GammaFitGUI::DecayCacheFileFor() const
{
    if (decayTh2Name_.empty()) return "";
    return CacheDirFor() + "/decay_fits_" + decayTh2Name_ + ".dat";
}

void GammaFitGUI::SaveDecayFitCache()
{
    std::string path = DecayCacheFileFor();
    if (path.empty()) return;
    EnsureCacheDir();
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
            << " " << r.bgType
            << " " << r.NsigLo
            << " " << r.rebin
            << " " << r.fitMethod
            << " " << r.fitLo
            << " " << r.fitHi
            << " " << (r.fullRange ? 1 : 0);
        int nlo = (int)r.paramLo.size();
        out << " " << nlo;
        for (int i = 0; i < nlo; i++) out << " " << r.paramLo[i] << " " << r.paramHi[i];
        int nec = (int)r.extraCuts.size();
        out << " " << nec;
        for (int i = 0; i < nec; i++) out << " " << r.extraCuts[i].first << " " << r.extraCuts[i].second;
        out << "\n";
    }
    PopulateCacheBrowser();
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
        int bgt = 1;
        if (ss >> bgt) r.bgType = bgt;   // optional; defaults 1 (Flat)
        double nsiglo = -1.0;
        if (ss >> nsiglo) r.NsigLo = nsiglo;
        int rebin = 1;
        if (ss >> rebin) r.rebin = rebin;
        int fitMethod = 2;
        if (ss >> fitMethod) r.fitMethod = fitMethod;
        double fitLo = 0.0, fitHi = 0.0;
        if (ss >> fitLo) r.fitLo = fitLo;
        if (ss >> fitHi) r.fitHi = fitHi;
        int fullRangeInt = 1;
        if (ss >> fullRangeInt) r.fullRange = (fullRangeInt != 0);
        int nBnds = 0;
        if (ss >> nBnds && nBnds > 0) {
            r.paramLo.resize(nBnds); r.paramHi.resize(nBnds);
            for (int i = 0; i < nBnds; i++) {
                double lo = 0, hi = 0;
                if (ss >> lo >> hi) { r.paramLo[i] = lo; r.paramHi[i] = hi; }
            }
        }
        int nec = 0;
        if (ss >> nec && nec > 0) {
            r.extraCuts.resize(nec);
            for (int i = 0; i < nec; i++) {
                double lo = 0, hi = 0;
                if (ss >> lo >> hi) r.extraCuts[i] = {lo, hi};
            }
        }
        decayFitStore_[r.peakE] = r;
        ++count;
    }
    AppendLog(Form("[Decay] Loaded %d cached decay fit(s) from %s", count, path.c_str()));
    PopulateCacheBrowser();
}

// ─────────────────────────────────────────────────────────────────────────────
// OnPreviewDecay  -  show decay projection without fitting
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnPreviewDecay()
{
    if (!inputFile_) { AppendLog("[Decay] No ROOT file loaded"); return; }
    if (decayTh2Name_.empty()) { AppendLog("[Decay] Run 'Refresh' first"); return; }

    Int_t peakSel = decayPeakList_->GetSelected();
    if (peakSel < 1 || (size_t)peakSel > decayPeakEs_.size()) {
        AppendLog("[Decay] Select a peak first"); return;
    }
    double E      = decayPeakEs_[peakSel - 1];
    double sig    = decayPeakSigs_[peakSel - 1];
    double NsigHi = decaySigRangeEntry_ ? decaySigRangeEntry_->GetNumber() : 1.0;
    double NsigLo = decaySigLoEntry_    ? decaySigLoEntry_->GetNumber()    : NsigHi;
    double eMin   = E - NsigLo * sig;
    double eMax   = E + NsigHi * sig;

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

    // Add any extra gamma cuts
    if (!extraDecayCuts_.empty())
        ApplyExtraCuts(h2, axisId, hDecay, extraDecayCuts_);

    // Rebin before display
    int rebin = decayRebinEntry_ ? (int)decayRebinEntry_->GetNumber() : 1;
    if (rebin > 1) hDecay->Rebin(rebin);

    // Restrict to t > 0
    double tLo = std::max(0.0, hDecay->GetXaxis()->GetXmin());
    hDecay->GetXaxis()->SetRangeUser(tLo, hDecay->GetXaxis()->GetXmax());

    {
        std::string ttl = Form("E=%.2f keV  [%.2f, %.2f keV]  -%.2f/+%.2f#sigma",
                               E, eMin, eMax, NsigLo, NsigHi);
        if (!extraDecayCuts_.empty())
            ttl += Form("  +%d extra cut(s)", (int)extraDecayCuts_.size());
        hDecay->SetTitle(ttl.c_str());
    }
    hDecay->GetXaxis()->SetTitle("Time (ms)");
    if (rebin > 1) {
        double binW = hDecay->GetBinWidth(1);
        hDecay->GetYaxis()->SetTitle(Form("Counts / (%.4g ms)", binW));
    } else {
        hDecay->GetYaxis()->SetTitle("Counts");
    }
    hDecay->SetLineColor(kBlack);
    hDecay->SetMarkerSize(0);

    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();
    bool previewErrBars = !decayErrBars_ || decayErrBars_->IsOn();
    hDecay->Draw(previewErrBars ? "hist E" : "hist");
    c->Modified(); c->Update();
    gSystem->ProcessEvents();

    SetStatus(Form("Preview: %.2f keV  [%.2f, %.2f]  cut: -%.2f/+%.2f sig", E, eMin, eMax, NsigLo, NsigHi));
}

// ─────────────────────────────────────────────────────────────────────────────
// OnDecayApplyLabel  -  save label/class from Decay tab fields to gamma cache
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnDecayApplyLabel()
{
    Int_t peakSel = decayPeakList_->GetSelected();
    if (peakSel < 1 || (size_t)peakSel > decayPeakKeys_.size()) {
        AppendLog("[Decay] Select a peak first"); return;
    }
    const std::string& cacheName = (peakSel <= (Int_t)decayPeakCacheNames_.size())
                                   ? decayPeakCacheNames_[peakSel - 1]
                                   : decayGammaProjName_;
    if (cacheName.empty()) { AppendLog("[Decay] No cache associated with this peak"); return; }

    FitDatabase fdb;
    fdb.Load(CacheFileFor(cacheName));
    const auto& entries = fdb.GetEntries();
    const std::string& key = decayPeakKeys_[peakSel - 1];
    auto it = entries.find(key);
    if (it == entries.end()) { AppendLog("[Decay] Cache entry not found"); return; }

    FitEntry e = it->second;
    if (decayLabelEntry_) {
        std::string lbl = decayLabelEntry_->GetText();
        if (!lbl.empty()) e.label = lbl;
    }
    // Auto-populate classification from labelClassMap_ (set via Isotopes tab)
    if (!e.label.empty() && labelClassMap_.count(e.label))
        e.classification = labelClassMap_.at(e.label);
    fdb.ForceStore(key, e);
    EnsureCacheDir();
    fdb.Save(CacheFileFor(cacheName));
    BackupCacheFile(CacheFileFor(cacheName));
    AppendLog("[Decay] Saved label=" + e.label + "  class=" + e.classification +
              "  -> " + key);
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
    double NsigHi = decaySigRangeEntry_ ? decaySigRangeEntry_->GetNumber() : 1.0;
    double NsigLo = decaySigLoEntry_    ? decaySigLoEntry_->GetNumber()    : NsigHi;
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

    // Gamma energy window (asymmetric cut)
    double eMin = E - NsigLo * sig;
    double eMax = E + NsigHi * sig;
    const double kSqrt2Pi = 2.5066282746310002;

    // Open output ROOT file for slice histograms (recreate each run)
    EnsureCacheDir();
    std::string sliceRootPath = CacheDirFor() + "/decay_slices_" +
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
        AppendLog("[Decay] No successful slice fits  -  try wider window or larger slice width");
        return;
    }

    // Plot
    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();

    int N = (int)tCenter.size();
    TGraphErrors* gr = new TGraphErrors(N);
    gr->SetName(Form("peakCounts_%.2fkeV", E));
    gr->SetTitle(Form("Peak counts vs time  (%.2f keV  -%.1f/+%.1f#sigma);Time (ms);Peak counts", E, NsigLo, NsigHi));
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

void GammaFitGUI::OnDecayRebinReset()
{
    if (decayRebinEntry_) decayRebinEntry_->SetNumber(1);
    OnPreviewDecay();
}

// ─────────────────────────────────────────────────────────────────────────────
// OnDecayScanCaches  -  populate cache picker with all non-decay .dat files
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnDecayScanCaches()
{
    if (!decayCacheCombo_) return;
    decayCacheCombo_->RemoveAll();

    // Collect histogram names from caches found in both the current per-file
    // subdirectory and the legacy flat directory.
    static const std::string kPrefix = "fit_cache_";
    static const std::string kSuffix = ".dat";

    std::vector<std::string> names;
    auto scanDir = [&](const std::string& dirPath) {
        void* d = gSystem->OpenDirectory(dirPath.c_str());
        if (!d) return;
        const char* ent;
        while ((ent = gSystem->GetDirEntry(d)) != nullptr) {
            std::string n = ent;
            if (n.size() <= kPrefix.size() + kSuffix.size()) continue;
            if (n.substr(0, kPrefix.size()) != kPrefix) continue;
            if (n.substr(n.size() - kSuffix.size()) != kSuffix) continue;
            // Strip prefix and suffix to get histogram name
            std::string hname = n.substr(kPrefix.size(),
                                         n.size() - kPrefix.size() - kSuffix.size());
            if (!hname.empty()) names.push_back(hname);
        }
        gSystem->FreeDirectory(d);
    };

    scanDir(CacheDirFor());  // current per-file cache dir

    // Also scan legacy flat dir if it differs from the current one
    std::string legacyDir = (launchDir_.empty() ? std::string(kCacheDir)
                                                 : launchDir_ + "/" + kCacheDir);
    if (legacyDir != CacheDirFor()) scanDir(legacyDir);

    // De-duplicate and sort
    std::sort(names.begin(), names.end());
    names.erase(std::unique(names.begin(), names.end()), names.end());

    int idx = 1;
    for (const auto& nm : names)
        decayCacheCombo_->AddEntry(nm.c_str(), idx++);
    decayCacheCombo_->MapSubwindows();
    decayCacheCombo_->Layout();

    if (names.empty())
        AppendLog("[Decay] No fit caches found  -  run AutoFit on a histogram first");
    else
        AppendLog(Form("[Decay] Found %d cache(s)", (int)names.size()));
}

// ─────────────────────────────────────────────────────────────────────────────
// Total Decay sub-tab slots
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnDecayTdModelChanged(Int_t id)
{
    if (decayTdPLbl_ && decayTdDLbl_ && decayTdGLbl_) {
        switch (id) {
            case 5:
                decayTdPLbl_->SetText("T1/2 beta-n product:");
                decayTdDLbl_->SetText("T1/2 Daughter:");
                break;
            case 6:
                decayTdPLbl_->SetText("T1/2 beta-2n product:");
                decayTdDLbl_->SetText("T1/2 Daughter:");
                break;
            case 7:
                decayTdPLbl_->SetText("T1/2 Parent:");
                decayTdDLbl_->SetText("T1/2 beta-n product:");
                decayTdGLbl_->SetText("T1/2 GDaughter:");
                break;
            case 8:
                decayTdPLbl_->SetText("T1/2 Parent:");
                decayTdDLbl_->SetText("T1/2 beta-2n product:");
                decayTdGLbl_->SetText("T1/2 GDaughter:");
                break;
            default:
                decayTdPLbl_->SetText("T1/2 Parent:");
                decayTdDLbl_->SetText("T1/2 Daughter:");
                decayTdGLbl_->SetText("T1/2 GDaughter:");
                break;
        }
    }
    UpdateTotalDecayEquation();
}

void GammaFitGUI::OnDecayTdBGTypeChanged(Int_t /*id*/)
{
    int bgType = decayTdBGCombo_ ? decayTdBGCombo_->GetSelected() : 1;
    bool hasExp = (bgType >= 2);
    if (decayTdThalfBGExp_) decayTdThalfBGExp_->SetState(hasExp ? kTRUE : kFALSE);
    if (decayTdFixBGExp_)   decayTdFixBGExp_->SetEnabled(hasExp ? kTRUE : kFALSE);
    UpdateTotalDecayEquation();
}

void GammaFitGUI::OnFitTotalDecay()
{
    if (!inputFile_) { AppendLog("[TotalDecay] No ROOT file loaded"); return; }
    if (decayTh2Name_.empty()) { AppendLog("[TotalDecay] Run 'Refresh' in Cuts tab first"); return; }
    if (decayPeakEs_.empty())  { AppendLog("[TotalDecay] No peaks loaded  -  refresh in Cuts tab"); return; }

    TH2* h2 = (TH2*)inputFile_->Get(decayTh2Name_.c_str());
    if (!h2) { AppendLog("[TotalDecay] TH2 not found"); return; }

    int axisId = decayGammaAxisCombo_ ? decayGammaAxisCombo_->GetSelected() : 1;
    double NsigHi = decaySigRangeEntry_ ? decaySigRangeEntry_->GetNumber() : 1.0;
    double NsigLo = decaySigLoEntry_    ? decaySigLoEntry_->GetNumber()    : NsigHi;

    // Determine which peaks to sum
    bool sumAll = !decayTdSumAll_ || decayTdSumAll_->IsOn();
    std::vector<int> peakIndices;
    if (sumAll) {
        for (int i = 0; i < (int)decayPeakEs_.size(); ++i) peakIndices.push_back(i);
    } else {
        Int_t sel = decayPeakList_ ? decayPeakList_->GetSelected() : -1;
        if (sel < 1 || (size_t)sel > decayPeakEs_.size()) {
            AppendLog("[TotalDecay] Select a peak in the Cuts list (or enable 'Sum all')");
            return;
        }
        peakIndices.push_back(sel - 1);
    }

    // Sum projections for each peak
    TH1* hTotal = nullptr;
    std::string summedPeaks;
    for (int idx : peakIndices) {
        double E      = decayPeakEs_[idx];
        double sig    = decayPeakSigs_[idx];
        // Use stored cut if available, otherwise current Lo/Hi entries
        double eMin, eMax;
        auto fitIt = decayFitStore_.find(E);
        if (fitIt != decayFitStore_.end()) {
            eMin = fitIt->second.eMin;
            eMax = fitIt->second.eMax;
        } else {
            eMin = E - NsigLo * sig;
            eMax = E + NsigHi * sig;
        }

        TH1* hPeak = nullptr;
        if (axisId == 1) {
            int b1 = h2->GetXaxis()->FindBin(eMin);
            int b2 = h2->GetXaxis()->FindBin(eMax);
            hPeak = h2->ProjectionY(Form("hTdPeak_%d", idx), b1, b2);
        } else {
            int b1 = h2->GetYaxis()->FindBin(eMin);
            int b2 = h2->GetYaxis()->FindBin(eMax);
            hPeak = h2->ProjectionX(Form("hTdPeak_%d", idx), b1, b2);
        }
        if (!hPeak || hPeak->GetEntries() == 0) { delete hPeak; continue; }
        hPeak->SetDirectory(nullptr);

        if (!hTotal) {
            hTotal = (TH1*)hPeak->Clone("hTotalDecay");
            hTotal->SetDirectory(nullptr);
            delete hPeak;
        } else {
            hTotal->Add(hPeak);
            delete hPeak;
        }
        summedPeaks += Form("%.2f ", E);
    }

    if (!hTotal || hTotal->GetEntries() == 0) {
        delete hTotal;
        AppendLog("[TotalDecay] No counts in summed projection");
        return;
    }

    // Rebin — use per-tab entry; fall back to Cuts tab entry
    int rebin = decayTdRebinEntry_ ? (int)decayTdRebinEntry_->GetNumber()
              : (decayRebinEntry_  ? (int)decayRebinEntry_->GetNumber() : 1);
    if (rebin > 1) hTotal->Rebin(rebin);

    double tAxisMin = std::max(0.0, hTotal->GetXaxis()->GetXmin());
    hTotal->GetXaxis()->SetRangeUser(tAxisMin, hTotal->GetXaxis()->GetXmax());
    hTotal->GetXaxis()->SetTitle("Time (ms)");
    if (rebin > 1) {
        double bw = hTotal->GetBinWidth(1);
        hTotal->GetYaxis()->SetTitle(Form("Counts / (%.4g ms)", bw));
    } else {
        hTotal->GetYaxis()->SetTitle("Counts");
    }
    hTotal->SetTitle(Form("Total decay  (%d peaks summed)  %s;Time (ms);Counts",
                         (int)peakIndices.size(), decayTh2Name_.c_str()));
    hTotal->SetLineColor(kBlack);
    hTotal->SetMarkerSize(0);

    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();
    hTotal->Draw(showErrorBars_ ? "hist E" : "hist");
    c->Modified(); c->Update();
    gSystem->ProcessEvents();

    // ── Build fit function ────────────────────────────────────────────────────
    int modelId = decayTdModelCombo_ ? decayTdModelCombo_->GetSelected() : 1;
    int bgType  = decayTdBGCombo_    ? decayTdBGCombo_->GetSelected()    : 1;
    if (bgType < 1 || bgType > 3) bgType = 1;
    int sigType = DecaySignalType(modelId);
    int nSig    = DecaySigNpar(sigType);

    bool fullRange = !decayTdFullRange_ || decayTdFullRange_->IsOn();
    double xlo = fullRange ? tAxisMin : std::max(0.0, decayTdFitLo_->GetNumber());
    double xhi = fullRange ? hTotal->GetXaxis()->GetXmax() : decayTdFitHi_->GetNumber();
    if (xlo >= xhi) { xlo = tAxisMin; xhi = hTotal->GetXaxis()->GetXmax(); }

    hTotal->GetXaxis()->SetRangeUser(xlo, xhi);

    TF1* fDecay = BuildDecayTF1("fTotalDecay", modelId, bgType, xlo, xhi);

    // Seed signal parameters
    double initMax = std::max(1.0, hTotal->GetMaximum() - hTotal->GetMinimum());
    double initMin = std::max(0.0, hTotal->GetMinimum());
    double TP  = decayTdThalfP_    ? decayTdThalfP_->GetNumber()    : (xhi-xlo)/2.0;
    double TD  = decayTdThalfD_    ? decayTdThalfD_->GetNumber()    : (xhi-xlo)/4.0;
    double TG  = decayTdThalfG_    ? decayTdThalfG_->GetNumber()    : (xhi-xlo)/8.0;
    double Tbg = decayTdThalfBGExp_? decayTdThalfBGExp_->GetNumber(): (xhi-xlo)*2.0;

    if (sigType == 1) {
        fDecay->SetParameter(0, initMax);
        fDecay->SetParameter(1, TP > 0 ? TP : (xhi-xlo)/2.0);
        fDecay->SetParLimits(1, 1e-6, 1e12);
        if (decayTdFixP_ && decayTdFixP_->IsOn()) fDecay->FixParameter(1, TP);
    } else if (sigType == 2) {
        fDecay->SetParameter(0, initMax);
        fDecay->SetParameter(1, TP > 0 ? TP : (xhi-xlo)/4.0);
        fDecay->SetParameter(2, TD > 0 ? TD : (xhi-xlo)/2.0);
        fDecay->SetParLimits(1, 1e-6, 1e12);
        fDecay->SetParLimits(2, 1e-6, 1e12);
        if (decayTdFixP_ && decayTdFixP_->IsOn()) fDecay->FixParameter(1, TP);
        if (decayTdFixD_ && decayTdFixD_->IsOn()) fDecay->FixParameter(2, TD);
    } else if (sigType == 3) {
        fDecay->SetParameter(0, initMax);
        fDecay->SetParameter(1, TP > 0 ? TP : (xhi-xlo)/8.0);
        fDecay->SetParameter(2, TD > 0 ? TD : (xhi-xlo)/4.0);
        fDecay->SetParameter(3, TG > 0 ? TG : (xhi-xlo)/2.0);
        fDecay->SetParLimits(1, 1e-6, 1e12);
        fDecay->SetParLimits(2, 1e-6, 1e12);
        fDecay->SetParLimits(3, 1e-6, 1e12);
        if (decayTdFixP_ && decayTdFixP_->IsOn()) fDecay->FixParameter(1, TP);
        if (decayTdFixD_ && decayTdFixD_->IsOn()) fDecay->FixParameter(2, TD);
        if (decayTdFixG_ && decayTdFixG_->IsOn()) fDecay->FixParameter(3, TG);
    }
    if (bgType == 1) {
        fDecay->SetParameter(nSig, initMin);
    } else if (bgType == 2) {
        fDecay->SetParameter(nSig,   initMin * 0.5);
        fDecay->SetParameter(nSig+1, initMin * 0.5);
        fDecay->SetParameter(nSig+2, Tbg > 0 ? Tbg : (xhi-xlo)*2.0);
        fDecay->SetParLimits(nSig+2, 1e-6, 1e12);
        if (decayTdFixBGExp_ && decayTdFixBGExp_->IsOn()) fDecay->FixParameter(nSig+2, Tbg);
    } else {
        fDecay->SetParameter(nSig,   initMin);
        fDecay->SetParameter(nSig+1, Tbg > 0 ? Tbg : (xhi-xlo)*2.0);
        fDecay->SetParLimits(nSig+1, 1e-6, 1e12);
        if (decayTdFixBGExp_ && decayTdFixBGExp_->IsOn()) fDecay->FixParameter(nSig+1, Tbg);
    }

    fDecay->SetLineColor(kRed);
    fDecay->SetLineWidth(2);

    // Build fit options from method combo
    int fitMethod = decayTdFitMethod_ ? decayTdFitMethod_->GetSelected() : 2;
    std::string fitOpts = "SR";
    if (fitMethod == 2) fitOpts = "SRL";       // Poisson log-likelihood
    if (fitMethod == 3) fitOpts = "SRE";       // Chi2 + MINOS errors
    if (fitMethod == 4) fitOpts = "SRLE";      // Likelihood + MINOS

    AppendLog(Form("[TotalDecay] Fitting with options \"%s\" ...", fitOpts.c_str()));
    Int_t savedStatTd = gStyle->GetOptStat();
    gStyle->SetOptStat(111111);
    TFitResultPtr fitRes = hTotal->Fit(fDecay, fitOpts.c_str());
    gStyle->SetOptStat(savedStatTd);

    double chi2    = fitRes.Get() ? fitRes->Chi2() : fDecay->GetChisquare();
    int    ndf     = fitRes.Get() ? fitRes->Ndf()  : fDecay->GetNDF();
    double chi2ndf = (ndf > 0) ? chi2 / ndf : -1.0;
    double pval    = (ndf > 0) ? TMath::Prob(chi2, ndf) : -1.0;
    int    fitStat = fitRes.Get() ? fitRes->Status() : -1;
    const char* statStr = (fitStat == 0) ? "OK" : Form("WARN(%d)", fitStat);

    bool tdErrBars = !decayTdErrBars_ || decayTdErrBars_->IsOn();
    bool showRes = decayTdShowResid_ && decayTdShowResid_->IsOn();
    if (showRes) {
        DrawDecayResiduals(hTotal, fDecay, xlo, xhi);
    } else {
        TCanvas* cv = canvas_->GetCanvas();
        cv->Clear(); cv->cd();
        hTotal->Draw(tdErrBars ? "hist E" : "hist");
        fDecay->Draw("same");

        int tdModelId = decayTdModelCombo_ ? decayTdModelCombo_->GetSelected() : 1;
        int tdSigType = DecaySignalType(tdModelId);
        TF1* fBGComp = nullptr;
        if (tdSigType != 4) {
            fBGComp = BuildDecayTF1("fTdBGComp", 4, bgType, xlo, xhi);
            int nBG = DecayBGNpar(bgType);
            for (int i = 0; i < nBG; i++)
                fBGComp->SetParameter(i, fDecay->GetParameter(nSig + i));
            fBGComp->SetLineColor(kGreen+2); fBGComp->SetLineWidth(2); fBGComp->SetLineStyle(2);
            fBGComp->Draw("same");
        }

        double lx = decayTdLegX_ ? decayTdLegX_->GetNumber() : 0.38;
        double ly = decayTdLegY_ ? decayTdLegY_->GetNumber() : 0.72;
        TLegend* leg = new TLegend(lx, ly, lx + 0.24, ly + 0.27);
        leg->SetBorderSize(1); leg->SetTextSize(0.016);
        leg->SetFillColor(kWhite); leg->SetFillStyle(1001);
        if (tdSigType != 4)
            leg->AddEntry((TObject*)nullptr,
                ("A = " + NNDCFormat(fDecay->GetParameter(0), fDecay->GetParError(0))).c_str(), "");
        if (tdSigType >= 1) leg->AddEntry((TObject*)nullptr,
            ("T_{1/2}(P) = " + NNDCFormat(fDecay->GetParameter(1), fDecay->GetParError(1)) + " ms").c_str(), "");
        if (tdSigType >= 2) leg->AddEntry((TObject*)nullptr,
            ("T_{1/2}(D) = " + NNDCFormat(fDecay->GetParameter(2), fDecay->GetParError(2)) + " ms").c_str(), "");
        if (tdSigType >= 3) leg->AddEntry((TObject*)nullptr,
            ("T_{1/2}(G) = " + NNDCFormat(fDecay->GetParameter(3), fDecay->GetParError(3)) + " ms").c_str(), "");
        {
            int bgP = nSig;
            if (bgType == 1) {
                leg->AddEntry((TObject*)nullptr,
                    ("BG = " + NNDCFormat(fDecay->GetParameter(bgP), fDecay->GetParError(bgP))).c_str(), "");
            } else if (bgType == 2) {
                leg->AddEntry((TObject*)nullptr,
                    ("BG_{flat} = " + NNDCFormat(fDecay->GetParameter(bgP), fDecay->GetParError(bgP))).c_str(), "");
                leg->AddEntry((TObject*)nullptr,
                    ("T_{1/2}(BG) = " + NNDCFormat(fDecay->GetParameter(bgP+2), fDecay->GetParError(bgP+2)) + " ms").c_str(), "");
            } else if (bgType == 3) {
                leg->AddEntry((TObject*)nullptr,
                    ("T_{1/2}(BG) = " + NNDCFormat(fDecay->GetParameter(bgP+1), fDecay->GetParError(bgP+1)) + " ms").c_str(), "");
            }
        }
        leg->AddEntry((TObject*)nullptr,
            Form("#chi^{2}/NDF=%.3f  %s", chi2ndf, statStr), "");
        leg->AddEntry(hTotal,  "Data", "l");
        leg->AddEntry(fDecay,  "Fit",  "l");
        if (fBGComp) leg->AddEntry(fBGComp, "BG", "l");
        leg->Draw();
        cv->Modified(); cv->Update();
    }

    // ── Report ────────────────────────────────────────────────────────────────
    if (decayTdResultView_) decayTdResultView_->Clear();
    if (decayTdResultView_) {
        decayTdResultView_->AddLine(Form("Total decay: %d peaks summed  TH2: %s",
                                        (int)peakIndices.size(), decayTh2Name_.c_str()));
        decayTdResultView_->AddLine(Form("Peaks: %s", summedPeaks.c_str()));
        const char* methodName = fitMethod==2?"Likelihood":fitMethod==3?"Chi2+MINOS":
                                  fitMethod==4?"Likelihood+MINOS":"Chi2";
        decayTdResultView_->AddLine(Form("Rebin: %d   method: %s   status: %d",
            rebin, methodName, fitRes.Get() ? fitRes->Status() : -1));
        decayTdResultView_->AddLine(Form("Chi2 = %.2f   NDF = %d   Chi2/NDF = %.4f",
                                        chi2, ndf, chi2ndf));
        decayTdResultView_->AddLine(Form("p-value = %.6f%s", pval,
            pval < 0.01 ? "  *** poor fit ***" : pval < 0.05 ? "  * marginal *" : ""));
        for (int i = 0; i < fDecay->GetNpar(); i++)
            decayTdResultView_->AddLine(Form("  %-16s = %11.5g +/- %.4g",
                fDecay->GetParName(i), fDecay->GetParameter(i), fDecay->GetParError(i)));

        // Population estimate
        const double ln2 = 0.6931471805599453;
        if (sigType != 4) {
            double A    = fDecay->GetParameter(0);
            double TP2  = fDecay->GetParameter(1);
            double binW = hTotal->GetBinWidth(1);
            double tau  = (TP2 > 0) ? TP2 / ln2 : 0.0;
            double Ntot = (tau > 0 && binW > 0) ? A * tau / binW : -1.0;
            if (Ntot > 0) {
                double totalErr = -1.0;
                if (fitRes.Get() && fitRes->IsValid()) {
                    const TMatrixDSym& cov = fitRes->GetCovarianceMatrix();
                    int np = fDecay->GetNpar();
                    if (cov.GetNrows() == np && A > 0 && TP2 > 0) {
                        double dA  = tau / binW;
                        double dTP = A / (ln2 * binW);
                        double varN = dA*dA*cov(0,0) + dTP*dTP*cov(1,1)
                                    + 2.0*dA*dTP*cov(0,1);
                        if (varN > 0) totalErr = std::sqrt(varN);
                    }
                }
                decayTdResultView_->AddLine("─────────────────────────────────");
                if (totalErr > 0)
                    decayTdResultView_->AddLine(
                        Form("Total peak counts (0->inf):  %.0f +/- %.0f", Ntot, totalErr));
                else
                    decayTdResultView_->AddLine(
                        Form("Total peak counts (0->inf):  %.0f", Ntot));
                decayTdResultView_->AddLine(
                    Form("  [A=%.4g  T1/2=%.4g ms  tau=%.4g ms  bw=%.4g ms]",
                         A, TP2, tau, binW));
                // Population = total parent decays N = A*tau/bw summed over peaks
                // (this already IS N for the total projection)
                decayTdResultView_->AddLine(
                    Form("Population (N_parent decays): %.4g", Ntot));
            }
        }
        decayTdResultView_->ShowBottom();
    }

    // ── Store and save ────────────────────────────────────────────────────────
    {
        decayTdFitResult_.model  = modelId;
        decayTdFitResult_.bgType = bgType;
        decayTdFitResult_.rebin  = rebin;
        decayTdFitResult_.chi2ndf= chi2ndf;
        decayTdFitResult_.status = fitRes.Get() ? fitRes->Status() : -1;
        decayTdFitResult_.histName = decayTh2Name_;
        decayTdFitResult_.eMin   = xlo;
        decayTdFitResult_.eMax   = xhi;
        decayTdFitResult_.Nsig   = (double)peakIndices.size();
        int np = fDecay->GetNpar();
        decayTdFitResult_.params.resize(np);
        decayTdFitResult_.errors.resize(np);
        for (int i = 0; i < np; i++) {
            decayTdFitResult_.params[i] = fDecay->GetParameter(i);
            decayTdFitResult_.errors[i] = fDecay->GetParError(i);
        }
        decayTdFitValid_ = true;
        SaveTotalDecayFitCache();
    }

    // Store snapshot for chi² scan popup
    delete lastDecayTdHist_; lastDecayTdHist_ = nullptr;
    delete lastDecayTdTF1_;  lastDecayTdTF1_  = nullptr;
    lastDecayTdHist_ = (TH1*)hTotal->Clone("_lastDecayTdH");
    lastDecayTdHist_->SetDirectory(nullptr);
    lastDecayTdTF1_  = (TF1*)fDecay->Clone("_lastDecayTdF");
    lastDecayTdXlo_  = xlo;
    lastDecayTdXhi_  = xhi;
    if (decayTdChiScanBtn_) decayTdChiScanBtn_->SetState(kButtonUp);

    AppendLog(Form("[TotalDecay] Fit done: %d peaks  chi2/ndf=%.3f", (int)peakIndices.size(), chi2ndf));
    SetStatus(Form("Total decay fit: %d peaks  chi2/ndf=%.3f", (int)peakIndices.size(), chi2ndf));
}

// ─────────────────────────────────────────────────────────────────────────────
// Total Decay cache helpers
// ─────────────────────────────────────────────────────────────────────────────

std::string GammaFitGUI::TotalDecayCacheFileFor() const
{
    if (decayTh2Name_.empty()) return "";
    return CacheDirFor() + "/decay_total_" + decayTh2Name_ + ".dat";
}

void GammaFitGUI::SaveTotalDecayFitCache()
{
    if (!decayTdFitValid_) return;
    std::string path = TotalDecayCacheFileFor();
    if (path.empty()) return;
    EnsureCacheDir();
    std::ofstream out(path);
    if (!out.is_open()) return;
    out << std::fixed << std::setprecision(8);
    const DecayFitResult& r = decayTdFitResult_;
    // Header: which peaks were summed (energies) and binning info
    out << "# Total decay fit for TH2: " << r.histName << "\n";
    out << "# nPeaksSummed " << (int)r.Nsig << "\n";
    out << "# rebin " << r.rebin << "\n";
    out << "# fitRange " << r.eMin << " " << r.eMax << "\n";
    // Per-peak energies from the summed projection
    if (!decayPeakEs_.empty()) {
        out << "# summedEnergies";
        for (double e : decayPeakEs_) out << " " << e;
        out << "\n";
    }
    // Fit result line: model bgType npar [pi ei]... chi2ndf status
    int np = (int)r.params.size();
    out << r.model << " " << r.bgType << " " << np;
    for (int i = 0; i < np; i++) out << " " << r.params[i] << " " << r.errors[i];
    out << " " << r.chi2ndf << " " << r.status << "\n";
    AppendLog("[TotalDecay] Cache saved -> " + path);
}

void GammaFitGUI::LoadTotalDecayFitCache()
{
    decayTdFitValid_ = false;
    std::string path = TotalDecayCacheFileFor();
    if (path.empty()) return;
    std::ifstream in(path);
    if (!in.is_open()) return;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        int model = 0, bgType = 1, np = 0;
        if (!(ss >> model >> bgType >> np)) continue;
        decayTdFitResult_.model  = model;
        decayTdFitResult_.bgType = bgType;
        decayTdFitResult_.params.resize(np);
        decayTdFitResult_.errors.resize(np);
        bool ok = true;
        for (int i = 0; i < np; i++)
            if (!(ss >> decayTdFitResult_.params[i] >> decayTdFitResult_.errors[i]))
                { ok = false; break; }
        if (!ok) continue;
        ss >> decayTdFitResult_.chi2ndf >> decayTdFitResult_.status;
        decayTdFitResult_.histName = decayTh2Name_;
        decayTdFitValid_ = true;
        // Seed T1/2 entries in the Total Decay tab
        int sig = DecaySignalType(model);
        if (sig >= 1 && np >= 2 && decayTdThalfP_) decayTdThalfP_->SetNumber(decayTdFitResult_.params[1]);
        if (sig >= 2 && np >= 3 && decayTdThalfD_) decayTdThalfD_->SetNumber(decayTdFitResult_.params[2]);
        if (sig >= 3 && np >= 4 && decayTdThalfG_) decayTdThalfG_->SetNumber(decayTdFitResult_.params[3]);
        if (decayTdModelCombo_) decayTdModelCombo_->Select(model, kTRUE);
        if (decayTdBGCombo_)    decayTdBGCombo_->Select(bgType, kTRUE);
        AppendLog("[TotalDecay] Loaded cached total decay fit from " + path);
        break;
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// DrawDecayResiduals — split-pad view: top=data+fit, bottom=pull (data-fit)/err
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::DrawDecayResiduals(TH1* h, TF1* fit, double xlo, double xhi)
{
    TCanvas* c = canvas_->GetCanvas();
    c->Clear();

    TPad* padTop = new TPad("dTop","",0.0,0.30,1.0,1.0);
    padTop->SetBottomMargin(0.12);
    padTop->SetLeftMargin(0.12);
    padTop->Draw();
    padTop->cd();

    if (h) {
        h->SetStats(0);
        h->Draw("hist E");
        if (fit) { fit->SetLineColor(kRed); fit->SetLineWidth(2); fit->Draw("same"); }
    }

    c->cd();
    TPad* padBot = new TPad("dBot","",0.0,0.0,1.0,0.30);
    padBot->SetTopMargin(0.02);
    padBot->SetBottomMargin(0.32);
    padBot->SetLeftMargin(0.12);
    padBot->Draw();
    padBot->cd();

    if (h && fit) {
        int b1 = h->FindBin(xlo), b2 = h->FindBin(xhi);
        int nb = b2 - b1 + 1;
        if (nb > 0) {
            TH1D* hPull = new TH1D("hDecayPull",";Time (ms);(data#minusfit)/#sigma",
                                    nb, h->GetBinLowEdge(b1), h->GetBinLowEdge(b2+1));
            hPull->SetBit(kCanDelete);
            for (int b = b1; b <= b2; b++) {
                double err = h->GetBinError(b);
                if (err > 0.0)
                    hPull->SetBinContent(b-b1+1,
                        (h->GetBinContent(b) - fit->Eval(h->GetBinCenter(b))) / err);
            }
            hPull->SetLineColor(kBlack);
            hPull->SetMarkerStyle(20);
            hPull->SetMarkerSize(0.6);
            hPull->GetXaxis()->SetLabelSize(0.10);
            hPull->GetXaxis()->SetTitleSize(0.10);
            hPull->GetXaxis()->SetTitleOffset(0.9);
            hPull->GetYaxis()->SetLabelSize(0.09);
            hPull->GetYaxis()->SetTitleSize(0.09);
            hPull->GetYaxis()->SetTitleOffset(0.5);
            hPull->GetYaxis()->SetNdivisions(505);
            hPull->Draw("P");

            TLine* zero = new TLine(xlo, 0, xhi, 0);
            zero->SetLineColor(kRed); zero->SetLineStyle(2);
            zero->Draw();
        }
    }
    c->Modified(); c->Update();
}

// ─────────────────────────────────────────────────────────────────────────────
// OnFitRK4Decay — 4-species decay chain via RK4 ODE solver
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnFitRK4Decay()
{
    if (!inputFile_) { AppendLog("[RK4] No ROOT file loaded"); return; }

    bool use1D = rkUse1DHist_ && rkUse1DHist_->IsOn();
    TH1* hTotal = nullptr;

    if (use1D) {
        // ── Mode B: direct 1D histogram ──────────────────────────────────────
        int selId = rkHistCombo_ ? rkHistCombo_->GetSelected() : -1;
        if (selId < 1 || (size_t)(selId - 1) >= rkHistSpecs_.size()) {
            AppendLog("[RK4] Select a histogram in Mode B (click Scan first)"); return;
        }
        std::string spec = rkHistSpecs_[selId - 1];
        auto sep = spec.find("::");
        std::string filePath = spec.substr(0, sep);
        std::string histName = spec.substr(sep + 2);

        TH1* hSrc = nullptr;
        if (filePath == "main") {
            TObject* obj = inputFile_->Get(histName.c_str());
            if (obj && obj->InheritsFrom(TH1::Class()))
                hSrc = (TH1*)((TH1*)obj)->Clone("hRK4Total_1D");
        } else {
            TFile* f = TFile::Open(filePath.c_str(), "READ");
            if (!f || f->IsZombie()) { AppendLog("[RK4] Cannot open: " + filePath); return; }
            TObject* obj = f->Get(histName.c_str());
            if (obj && obj->InheritsFrom(TH1::Class()))
                hSrc = (TH1*)((TH1*)obj)->Clone("hRK4Total_1D");
            f->Close(); delete f;
        }
        if (!hSrc) { AppendLog("[RK4] Histogram not found: " + histName); return; }
        hSrc->SetDirectory(nullptr);
        if (hSrc->GetEntries() == 0) { delete hSrc; AppendLog("[RK4] Histogram is empty"); return; }
        hTotal = hSrc;
        AppendLog(Form("[RK4] Mode B: '%s' (%.0f entries)", histName.c_str(), hTotal->GetEntries()));
    } else {
        // ── Mode A: sum TH2 projections ───────────────────────────────────────
        if (decayTh2Name_.empty()) { AppendLog("[RK4] Run 'Refresh' in Cuts tab first"); return; }
        if (decayPeakEs_.empty())  { AppendLog("[RK4] No peaks loaded — refresh in Cuts tab"); return; }

        TH2* h2 = (TH2*)inputFile_->Get(decayTh2Name_.c_str());
        if (!h2) { AppendLog("[RK4] TH2 not found: " + decayTh2Name_); return; }

        int axisId = decayGammaAxisCombo_ ? decayGammaAxisCombo_->GetSelected() : 1;
        double NsigHi = decaySigRangeEntry_ ? decaySigRangeEntry_->GetNumber() : 1.0;
        double NsigLo = decaySigLoEntry_    ? decaySigLoEntry_->GetNumber()    : NsigHi;

        std::vector<int> peakIdx;
        if (rkPeakList_) {
            TList selList;
            rkPeakList_->GetSelectedEntries(&selList);
            TIter next(&selList);
            while (TGLBEntry* entry = (TGLBEntry*)next()) {
                int id = entry->EntryId();
                if (id >= 1 && (size_t)id <= decayPeakEs_.size())
                    peakIdx.push_back(id - 1);
            }
        }
        if (peakIdx.empty()) {
            for (int i = 0; i < (int)decayPeakEs_.size(); ++i) peakIdx.push_back(i);
            AppendLog("[RK4] No peaks selected — using all peaks");
        }

        for (int idx : peakIdx) {
            double E   = decayPeakEs_[idx];
            double sig = decayPeakSigs_[idx];
            double eMin, eMax;
            auto fitIt = decayFitStore_.find(E);
            if (fitIt != decayFitStore_.end()) { eMin = fitIt->second.eMin; eMax = fitIt->second.eMax; }
            else { eMin = E - NsigLo * sig;   eMax = E + NsigHi * sig; }

            TH1* hp = nullptr;
            if (axisId == 1) {
                int b1 = h2->GetXaxis()->FindBin(eMin), b2 = h2->GetXaxis()->FindBin(eMax);
                hp = h2->ProjectionY(Form("hRK4_%d", idx), b1, b2);
            } else {
                int b1 = h2->GetYaxis()->FindBin(eMin), b2 = h2->GetYaxis()->FindBin(eMax);
                hp = h2->ProjectionX(Form("hRK4_%d", idx), b1, b2);
            }
            if (!hp || hp->GetEntries() == 0) { delete hp; continue; }
            hp->SetDirectory(nullptr);
            if (!hTotal) { hTotal = (TH1*)hp->Clone("hRK4Total"); hTotal->SetDirectory(nullptr); delete hp; }
            else         { hTotal->Add(hp); delete hp; }
        }
    }

    if (!hTotal || hTotal->GetEntries() == 0) {
        delete hTotal;
        AppendLog("[RK4] No counts in projection — check cuts / peak selection"); return;
    }

    int rebin = rkRebinEntry_ ? (int)rkRebinEntry_->GetNumber()
              : (decayRebinEntry_ ? (int)decayRebinEntry_->GetNumber() : 1);
    if (rebin > 1) hTotal->Rebin(rebin);

    double tAxisMin = std::max(0.0, hTotal->GetXaxis()->GetXmin());
    hTotal->GetXaxis()->SetRangeUser(tAxisMin, hTotal->GetXaxis()->GetXmax());
    hTotal->SetTitle(Form("RK4 Chain Fit  %s", decayTh2Name_.c_str()));
    hTotal->GetXaxis()->SetTitle("Time (ms)");
    double bw = hTotal->GetBinWidth(1);
    hTotal->GetYaxis()->SetTitle(Form("Counts / %.4g ms", bw));
    hTotal->SetLineColor(kBlack);

    double xlo = tAxisMin;
    double xhi = hTotal->GetXaxis()->GetXmax();
    if (rkFitFull_ && !rkFitFull_->IsOn()) {
        xlo = rkFitLo_ ? rkFitLo_->GetNumber() : xlo;
        xhi = rkFitHi_ ? rkFitHi_->GetNumber() : xhi;
    }
    hTotal->GetXaxis()->SetRangeUser(xlo, xhi);

    // ── Read parameters from UI ───────────────────────────────────────────────
    const double ln2 = 0.6931471805599453;
    auto getT = [&](TGNumberEntry* e, double def) { return e ? e->GetNumber() : def; };
    auto isFixed = [&](TGCheckButton* c) { return c && c->IsOn(); };

    double T_P  = getT(rkThalfP_,  100.0);
    double T_D  = getT(rkThalfD_,  200.0);
    double T_Bn = getT(rkThalfBn_, 500.0);
    double T_Gd = getT(rkThalfGd_,1000.0);
    double N0   = getT(rkN0_,  220000.0);
    double bg   = getT(rkBg_,       0.0);
    double Pn   = getT(rkPn_,       0.2);
    int    nSt  = rkNsteps_ ? (int)rkNsteps_->GetNumber() : 200;

    int chainMode = rkChainModel_ ? rkChainModel_->GetSelected() : 1;
    // chainMode 1 = full P+D+Bn+Gd, 2 = P+D only, 3 = parent only
    if (chainMode == 2) Pn = 0.0;   // no beta-n branch

    // ── Build RK4 TF1 ────────────────────────────────────────────────────────
    // par layout: [lam_P, lam_D, lam_Bn, lam_Gd, N0, bg, Pn]  (eff=1 for all)
    TF1* fRK4 = new TF1("fRK4Decay",
        [nSt, chainMode](double* x, double* p) -> double {
            double t = x[0];
            if (t <= 0.0) return p[5]; // bg only before t=0

            double lP  = p[0], lD  = p[1], lBn = p[2], lGd = p[3];
            double N0p = p[4];
            double Pn  = (chainMode >= 2) ? 0.0 : p[6];  // force Pn=0 for P+D / Parent only
            double dt  = t / nSt;

            double Np = N0p, Nd = 0.0, Nbn = 0.0, Ngd = 0.0;

            for (int i = 0; i < nSt; i++) {
                // RK4 slopes — k1
                double k1p  = -lP * Np;
                double k1d  =  (1.0 - Pn) * lP * Np - lD  * Nd;
                double k1bn =  Pn * lP * Np           - lBn * Nbn;
                double k1gd =  lBn * Nbn               - lGd * Ngd;
                // k2
                double Np2  = Np  + 0.5*dt*k1p;
                double Nd2  = Nd  + 0.5*dt*k1d;
                double Nbn2 = Nbn + 0.5*dt*k1bn;
                double Ngd2 = Ngd + 0.5*dt*k1gd;
                double k2p  = -lP * Np2;
                double k2d  =  (1.0-Pn)*lP*Np2 - lD*Nd2;
                double k2bn =  Pn*lP*Np2        - lBn*Nbn2;
                double k2gd =  lBn*Nbn2         - lGd*Ngd2;
                // k3
                double Np3  = Np  + 0.5*dt*k2p;
                double Nd3  = Nd  + 0.5*dt*k2d;
                double Nbn3 = Nbn + 0.5*dt*k2bn;
                double Ngd3 = Ngd + 0.5*dt*k2gd;
                double k3p  = -lP * Np3;
                double k3d  =  (1.0-Pn)*lP*Np3 - lD*Nd3;
                double k3bn =  Pn*lP*Np3        - lBn*Nbn3;
                double k3gd =  lBn*Nbn3         - lGd*Ngd3;
                // k4
                double Np4  = Np  + dt*k3p;
                double Nd4  = Nd  + dt*k3d;
                double Nbn4 = Nbn + dt*k3bn;
                double Ngd4 = Ngd + dt*k3gd;
                double k4p  = -lP * Np4;
                double k4d  =  (1.0-Pn)*lP*Np4 - lD*Nd4;
                double k4bn =  Pn*lP*Np4        - lBn*Nbn4;
                double k4gd =  lBn*Nbn4         - lGd*Ngd4;
                // update
                Np  += dt/6.0*(k1p  + 2*k2p  + 2*k3p  + k4p);
                Nd  += dt/6.0*(k1d  + 2*k2d  + 2*k3d  + k4d);
                Nbn += dt/6.0*(k1bn + 2*k2bn + 2*k3bn + k4bn);
                Ngd += dt/6.0*(k1gd + 2*k2gd + 2*k3gd + k4gd);
            }

            double counts = lP*Np;
            if (chainMode == 1) counts += lD*Nd + lBn*Nbn + lGd*Ngd;  // full chain
            if (chainMode == 2) counts += lD*Nd;                        // P+D only
            // chainMode==3: parent only — just lP*Np
            counts += p[5]; // bg
            return std::max(counts, 0.0);
        },
        xlo, xhi, 7);

    fRK4->SetParNames("lam_P","lam_D","lam_Bn","lam_Gd","N0","bg","Pn");

    // Convert T1/2 → lambda, guard against zero
    auto lam = [&](double T) { return (T > 0.0) ? ln2 / T : ln2 / 1.0; };

    fRK4->SetParameters(lam(T_P), lam(T_D), lam(T_Bn), lam(T_Gd), N0, bg, Pn);

    // Limits
    fRK4->SetParLimits(0, ln2/1e7, ln2/0.001); // lam_P
    fRK4->SetParLimits(1, ln2/1e7, ln2/0.001); // lam_D
    fRK4->SetParLimits(2, ln2/1e7, ln2/0.001); // lam_Bn
    fRK4->SetParLimits(3, ln2/1e7, ln2/0.001); // lam_Gd
    fRK4->SetParLimits(4, 0.0, 1e9);            // N0
    fRK4->SetParLimits(6, 0.0, 1.0);            // Pn

    // Fix parameters as requested
    auto fixPar = [&](int i, double v) { fRK4->FixParameter(i, v); };
    if (isFixed(rkFixP_))    fixPar(0, lam(T_P));
    if (isFixed(rkFixD_))    fixPar(1, lam(T_D));
    if (isFixed(rkFixBn_))   fixPar(2, lam(T_Bn));
    if (isFixed(rkFixGd_))   fixPar(3, lam(T_Gd));
    if (isFixed(rkFixN0_))   fixPar(4, N0);
    if (isFixed(rkFixBg_))   fixPar(5, bg);
    if (isFixed(rkFixPn_))   fixPar(6, Pn);
    // Chain model constraints
    if (chainMode == 2) { fixPar(6, 0.0); }                              // P+D: Pn=0
    if (chainMode == 3) { fixPar(1, lam(T_D)); fixPar(2, lam(T_Bn));    // Parent: fix D/Bn/Gd
                          fixPar(3, lam(T_Gd)); fixPar(6, 0.0); }

    fRK4->SetLineColor(kRed);
    fRK4->SetLineWidth(2);

    bool rkErrBars = !rkErrBars_ || rkErrBars_->IsOn();
    const char* rkDrawOpt = rkErrBars ? "hist E" : "hist";

    // ── Draw histogram first ──────────────────────────────────────────────────
    TCanvas* cv = canvas_->GetCanvas();
    cv->Clear(); cv->cd();
    Int_t savedStatRK = gStyle->GetOptStat();
    gStyle->SetOptStat(111111);
    hTotal->Draw(rkDrawOpt);
    cv->Modified(); cv->Update();
    gStyle->SetOptStat(savedStatRK);
    gSystem->ProcessEvents();

    // ── Build fit options ─────────────────────────────────────────────────────
    int rkMethod = rkFitMethod_ ? rkFitMethod_->GetSelected() : 2;
    std::string rkOpts = "RSN";
    if (rkMethod == 2) rkOpts = "RSNL";
    if (rkMethod == 3) rkOpts = "RSNE";
    if (rkMethod == 4) rkOpts = "RSNLE";

    AppendLog(Form("[RK4] Fitting range [%.1f,%.1f] ms, steps=%d, opts=%s ...",
                   xlo, xhi, nSt, rkOpts.c_str()));

    TFitResultPtr fitRes = hTotal->Fit(fRK4, rkOpts.c_str(), "", xlo, xhi);

    double chi2val = fRK4->GetChisquare();
    int    ndfVal  = fRK4->GetNDF();
    double chi2ndf = (ndfVal > 0) ? chi2val / ndfVal : -1.0;
    double pval    = (ndfVal > 0) ? TMath::Prob(chi2val, ndfVal) : -1.0;
    int    rkStat  = fitRes.Get() ? fitRes->Status() : -1;
    const char* rkStatStr = (rkStat == 0) ? "OK" : Form("WARN(%d)", rkStat);

    bool showRKRes = rkShowResid_ && rkShowResid_->IsOn();
    if (showRKRes) {
        DrawDecayResiduals(hTotal, fRK4, xlo, xhi);
    } else {
        cv->cd();
        hTotal->Draw(rkDrawOpt);
        fRK4->Draw("same");

        // ── Legend with equation, T1/2, chi2, status ────────────────────────
        double lx = rkLegX_ ? rkLegX_->GetNumber() : 0.38;
        double ly = rkLegY_ ? rkLegY_->GetNumber() : 0.66;
        TLegend* leg = new TLegend(lx, ly, lx + 0.24, ly + 0.32);
        leg->SetBorderSize(1); leg->SetTextSize(0.016);
        leg->SetFillColor(kWhite); leg->SetFillStyle(1001);
        auto tHalf = [&](int p) -> double {
            double lv = fRK4->GetParameter(p); return (lv > 0) ? ln2/lv : -1.0; };
        auto tErr  = [&](int p) -> double {
            double lv = fRK4->GetParameter(p), er = fRK4->GetParError(p);
            return (lv > 0 && er > 0) ? ln2*er/(lv*lv) : 0.0; };
        leg->AddEntry((TObject*)nullptr,
            ("T_{1/2}(P) = " + NNDCFormat(tHalf(0), tErr(0)) + " ms").c_str(), "");
        if (chainMode <= 2) leg->AddEntry((TObject*)nullptr,
            ("T_{1/2}(D) = " + NNDCFormat(tHalf(1), tErr(1)) + " ms").c_str(), "");
        if (chainMode == 1) {
            leg->AddEntry((TObject*)nullptr,
                ("T_{1/2}(Bn) = " + NNDCFormat(tHalf(2), tErr(2)) + " ms").c_str(), "");
            leg->AddEntry((TObject*)nullptr,
                ("T_{1/2}(Gd) = " + NNDCFormat(tHalf(3), tErr(3)) + " ms").c_str(), "");
            leg->AddEntry((TObject*)nullptr,
                ("Pn = " + NNDCFormat(fRK4->GetParameter(6), fRK4->GetParError(6))).c_str(), "");
        }
        leg->AddEntry((TObject*)nullptr,
            Form("#chi^{2}/NDF=%.3f  %s", chi2ndf, rkStatStr), "");
        leg->AddEntry(hTotal, "Data",    "l");
        leg->AddEntry(fRK4,   "RK4 fit", "l");
        leg->Draw();
        cv->Modified(); cv->Update();
    }

    // ── Report results ────────────────────────────────────────────────────────
    if (rkResultView_) {
        rkResultView_->Clear();
        const char* rkMethodName = rkMethod==2?"Likelihood":rkMethod==3?"Chi2+MINOS":
                                    rkMethod==4?"Likelihood+MINOS":"Chi2";
        rkResultView_->AddLine(Form("Method: %s   Rebin: %d", rkMethodName, rebin));
        rkResultView_->AddLine(Form("Fit status: %d  (%s)", rkStat, rkStatStr));
        rkResultView_->AddLine(Form("Chi2 = %.2f   NDF = %d   Chi2/NDF = %.4f",
                                    chi2val, ndfVal, chi2ndf));
        rkResultView_->AddLine(Form("p-value = %.6f%s", pval,
            pval < 0.01 ? "  *** poor fit ***" : pval < 0.05 ? "  * marginal *" : ""));
        rkResultView_->AddLine("---  Half-lives  ---");

        auto reportT = [&](const char* name, int parIdx) -> double {
            double lv  = fRK4->GetParameter(parIdx);
            double err = fRK4->GetParError(parIdx);
            double T   = (lv > 0) ? ln2 / lv : -1.0;
            double Terr= (lv > 0 && err > 0) ? ln2 * err / (lv * lv) : 0.0;
            rkResultView_->AddLine(Form("  %-12s  T1/2 = %.4f +/- %.4f ms",
                                        name, T, Terr));
            return T;
        };

        double tpFit  = reportT("Parent",    0);
        double tdFit  = reportT("Daughter",  1);
        double tbnFit = reportT("beta-n br.",2);
        double tgdFit = reportT("GDaughter", 3);

        rkResultView_->AddLine("---  Other  ---");
        rkResultView_->AddLine(Form("  N0   = %.4g +/- %.4g",
                                    fRK4->GetParameter(4), fRK4->GetParError(4)));
        rkResultView_->AddLine(Form("  bg   = %.4g +/- %.4g",
                                    fRK4->GetParameter(5), fRK4->GetParError(5)));
        rkResultView_->AddLine(Form("  Pn   = %.4f +/- %.4f",
                                    fRK4->GetParameter(6), fRK4->GetParError(6)));
        rkResultView_->Update();

        if (rkThalfP_  && tpFit  > 0) rkThalfP_ ->SetNumber(tpFit);
        if (rkThalfD_  && tdFit  > 0) rkThalfD_ ->SetNumber(tdFit);
        if (rkThalfBn_ && tbnFit > 0) rkThalfBn_->SetNumber(tbnFit);
        if (rkThalfGd_ && tgdFit > 0) rkThalfGd_->SetNumber(tgdFit);
    }

    // Store snapshot for chi² scan popup
    delete lastRK4Hist_; lastRK4Hist_ = nullptr;
    delete lastRK4TF1_;  lastRK4TF1_  = nullptr;
    lastRK4Hist_ = (TH1*)hTotal->Clone("_lastRK4H");
    lastRK4Hist_->SetDirectory(nullptr);
    lastRK4TF1_  = (TF1*)fRK4->Clone("_lastRK4F");
    lastRK4Xlo_  = xlo;
    lastRK4Xhi_  = xhi;
    if (rkChiScanBtn_) rkChiScanBtn_->SetState(kButtonUp);

    AppendLog(Form("[RK4] Done  chi2/NDF=%.4f  p=%.4f  T1/2_P=%.3f ms", chi2ndf, pval,
                   (fRK4->GetParameter(0) > 0) ? ln2 / fRK4->GetParameter(0) : -1.0));
    SetStatus(Form("RK4 chain fit: chi2/NDF=%.3f  p=%.4f", chi2ndf, pval));
    delete hTotal;
}

// ─────────────────────────────────────────────────────────────────────────────
// RK4 helper slots
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnRK4SelectAll()
{
    if (!rkPeakList_) return;
    for (int i = 1; i <= (int)decayPeakEs_.size(); ++i)
        rkPeakList_->Select(i, kTRUE);
    rkPeakList_->MapSubwindows();
    rkPeakList_->Layout();
}

void GammaFitGUI::OnRK4SelectNone()
{
    if (!rkPeakList_) return;
    for (int i = 1; i <= (int)decayPeakEs_.size(); ++i)
        rkPeakList_->Select(i, kFALSE);
    rkPeakList_->MapSubwindows();
    rkPeakList_->Layout();
}

void GammaFitGUI::OnRK4SeedFromCuts()
{
    // Copy the T1/2 values from the currently selected Cuts peak fit result
    Int_t sel = decayPeakList_ ? decayPeakList_->GetSelected() : -1;
    if (sel < 1 || (size_t)sel > decayPeakEs_.size()) {
        AppendLog("[RK4] Select a peak in the Cuts tab first to seed from");
        return;
    }
    double E = decayPeakEs_[sel - 1];
    auto it = decayFitStore_.find(E);
    if (it == decayFitStore_.end() || it->second.params.empty()) {
        AppendLog(Form("[RK4] No fit result stored for peak %.2f keV — fit it in Cuts first", E));
        return;
    }
    const DecayFitResult& r = it->second;
    int sigType = DecaySignalType(r.model);
    // params: [A, T_P, ...], T1/2 are in ms already
    if (r.params.size() >= 2 && rkThalfP_)
        rkThalfP_->SetNumber(r.params[1]);
    if (sigType >= 2 && r.params.size() >= 3 && rkThalfD_)
        rkThalfD_->SetNumber(r.params[2]);
    if (sigType >= 3 && r.params.size() >= 4 && rkThalfBn_)
        rkThalfBn_->SetNumber(r.params[3]);
    AppendLog(Form("[RK4] Seeded T1/2 from Cuts fit for peak %.2f keV (model %d)", E, r.model));
}

void GammaFitGUI::OnRK4ScanHistograms()
{
    if (!rkHistCombo_) return;
    rkHistCombo_->RemoveAll();
    rkHistSpecs_.clear();
    int comboIdx = 1;

    auto scanFile = [&](TFile* f, const std::string& filePath) {
        if (!f || f->IsZombie()) return;
        TList* keys = f->GetListOfKeys();
        if (!keys) return;
        TIter next(keys);
        while (TKey* k = (TKey*)next()) {
            TClass* cl = TClass::GetClass(k->GetClassName());
            if (!cl || !cl->InheritsFrom(TH1::Class()) || cl->InheritsFrom(TH2::Class())) continue;
            TH1* h = (TH1*)k->ReadObj();
            if (!h || h->GetNbinsX() < 10) { delete h; continue; }
            bool isFile = (filePath == "main");
            std::string label = Form("[%s] %s  (%d bins, %.1f-%.1f)",
                isFile ? "file" : "cache", k->GetName(), h->GetNbinsX(),
                h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
            rkHistCombo_->AddEntry(label.c_str(), comboIdx++);
            rkHistSpecs_.push_back(filePath + "::" + k->GetName());
            delete h;
        }
    };

    if (inputFile_) scanFile(inputFile_, "main");

    // Also scan decay-slice cache file if it exists
    std::string slicePath = CacheDirFor() + "/decay_slices_" + decayTh2Name_ + ".root";
    {
        TFile* sf = TFile::Open(slicePath.c_str(), "READ");
        if (sf && !sf->IsZombie()) { scanFile(sf, slicePath); sf->Close(); delete sf; }
    }

    rkHistCombo_->MapSubwindows();
    rkHistCombo_->Layout();
    int found = comboIdx - 1;
    if (found > 0) {
        rkHistCombo_->Select(1, kFALSE);
        AppendLog(Form("[RK4] Scan found %d 1D histograms", found));
    } else {
        AppendLog("[RK4] No 1D histograms found — check file and decay slice cache");
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Chi² scan popup  —  shared helper + three slots
// ─────────────────────────────────────────────────────────────────────────────

static void DoChiScanPopup(TH1* h, TF1* f,
                           double xlo, double xhi,
                           const std::string& title)
{
    if (!h || !f) return;

    int npar = f->GetNpar();
    // Collect free parameters: those whose stored error > 0
    std::vector<int>         freePars;
    std::vector<std::string> freeNames;
    for (int i = 0; i < npar; i++) {
        if (f->GetParError(i) > 0.0) {
            freePars.push_back(i);
            freeNames.push_back(f->GetParName(i));
        }
    }
    if (freePars.empty()) {
        printf("[ChiScan] No free parameters found (all errors = 0)\n");
        return;
    }

    int N     = (int)freePars.size();
    int ncols = (N <= 3) ? N : 3;
    int nrows = (N + ncols - 1) / ncols;
    int cw    = 310 * ncols;
    int ch    = 270 * nrows;

    static int scanIdx = 0;
    TCanvas* cScan = new TCanvas(Form("cChiScan_%d", scanIdx++),
                                  title.c_str(), cw, ch);
    cScan->Divide(ncols, nrows);

    // Clone TF1 for scanning so we never alter the drawn copy
    TF1* fs = (TF1*)f->Clone(Form("_fScanTmp_%p", (void*)f));

    const int nSteps = 50;

    for (int ip = 0; ip < N; ip++) {
        cScan->cd(ip + 1);
        gPad->SetLeftMargin(0.17);
        gPad->SetBottomMargin(0.17);

        int    pi  = freePars[ip];
        double val = f->GetParameter(pi);
        double err = f->GetParError(pi);
        double plo = val - 3.5 * err;
        double phi = val + 3.5 * err;

        // Compute chi² by holding param pi at each scan point; all others fixed
        // at their fitted values (simple surface slice — fast, no refit).
        std::vector<double> xv(nSteps), yv(nSteps);
        for (int j = 0; j < nSteps; j++) {
            double pv = plo + (phi - plo) * j / (nSteps - 1);
            // Reset all to fitted values then override scan parameter
            for (int k = 0; k < npar; k++) fs->SetParameter(k, f->GetParameter(k));
            fs->SetParameter(pi, pv);

            double chi2 = 0.0;
            for (int bin = 1; bin <= h->GetNbinsX(); bin++) {
                double x = h->GetBinCenter(bin);
                if (x < xlo || x > xhi) continue;
                double y = h->GetBinContent(bin);
                double e = h->GetBinError(bin);
                if (e <= 0.0) continue;
                double fy = fs->Eval(x);
                chi2 += (y - fy) * (y - fy) / (e * e);
            }
            xv[j] = pv;
            yv[j] = chi2;
        }

        // Find minimum chi² in scan
        double chi2min = *std::min_element(yv.begin(), yv.end());

        TGraph* gr = new TGraph(nSteps, xv.data(), yv.data());
        gr->SetName(Form("gScan_p%d", pi));
        gr->SetTitle(Form(";%s;#chi^{2}", freeNames[ip].c_str()));
        gr->SetLineColor(kBlue + 1);
        gr->SetLineWidth(2);
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(0.4);
        gr->Draw("ALP");
        gPad->Modified(); gPad->Update();

        double yaxis_lo = gr->GetHistogram()->GetMinimum();
        double yaxis_hi = gr->GetHistogram()->GetMaximum();

        // Vertical line at best-fit value
        TLine* vl = new TLine(val, yaxis_lo, val, yaxis_hi);
        vl->SetLineColor(kRed);
        vl->SetLineStyle(2);
        vl->SetLineWidth(2);
        vl->Draw();

        // 1σ and 2σ horizontal lines (Δχ² = 1 and 4)
        TLine* h1 = new TLine(plo, chi2min + 1.0, phi, chi2min + 1.0);
        h1->SetLineColor(kGreen + 2);
        h1->SetLineStyle(3);
        h1->SetLineWidth(2);
        h1->Draw();

        TLine* h2 = new TLine(plo, chi2min + 4.0, phi, chi2min + 4.0);
        h2->SetLineColor(kOrange + 1);
        h2->SetLineStyle(3);
        h2->SetLineWidth(2);
        h2->Draw();

        TLegend* leg = new TLegend(0.52, 0.74, 0.97, 0.97);
        leg->SetBorderSize(0);
        leg->SetTextSize(0.038);
        leg->SetFillStyle(0);
        leg->AddEntry(vl, Form("Best: %.4g", val), "l");
        leg->AddEntry(h1, "#Delta#chi^{2} = 1  (1#sigma)", "l");
        leg->AddEntry(h2, "#Delta#chi^{2} = 4  (2#sigma)", "l");
        leg->Draw();

        gPad->Modified(); gPad->Update();
    }

    delete fs;
    cScan->Update();
}

void GammaFitGUI::OnChiScanDecay()
{
    DoChiScanPopup(lastDecayCutsHist_, lastDecayCutsTF1_,
                   lastDecayCutsXlo_,  lastDecayCutsXhi_,
                   "chi2 vs Params — Cuts Decay");
}

void GammaFitGUI::OnChiScanTotalDecay()
{
    DoChiScanPopup(lastDecayTdHist_, lastDecayTdTF1_,
                   lastDecayTdXlo_,  lastDecayTdXhi_,
                   "chi2 vs Params — Total Decay");
}

void GammaFitGUI::OnChiScanRK4()
{
    DoChiScanPopup(lastRK4Hist_, lastRK4TF1_,
                   lastRK4Xlo_,  lastRK4Xhi_,
                   "chi2 vs Params — RK4 Chain");
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildIsotopesTab
// ─────────────────────────────────────────────────────────────────────────────

