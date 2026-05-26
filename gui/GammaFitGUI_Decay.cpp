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
            sigRow->AddFrame(new TGLabel(sigRow, "Hi sig:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
            decaySigRangeEntry_ = new TGNumberEntry(sigRow, 1.0, 5, -1,
                                                    TGNumberFormat::kNESRealFour,
                                                    TGNumberFormat::kNEAPositive);
            decaySigRangeEntry_->SetWidth(55);
            sigRow->AddFrame(decaySigRangeEntry_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            TGTextButton* refreshBtn = new TGTextButton(sigRow, "Refresh");
            sigRow->AddFrame(refreshBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
            refreshBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRefreshDecayPeaks()");
            refreshBtn->SetToolTipText("Load peaks from the selected cache.");
            TGTextButton* previewBtn = new TGTextButton(sigRow, "Preview");
            sigRow->AddFrame(previewBtn, new TGLayoutHints(kLHintsLeft));
            previewBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPreviewDecay()");
            previewBtn->SetToolTipText("Show decay projection without fitting");

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
            AddHalfLifeRow(grp, "T1/2 Daughter:",  decayThalfDLbl_, decayThalfD_,    decayFixD_);
            AddHalfLifeRow(grp, "T1/2 GDaughter:", decayThalfGLbl_, decayThalfG_,    decayFixG_);
            AddHalfLifeRow(grp, "T1/2 Exp BG:",    decayThalfBGLbl_,decayThalfBGExp_,decayFixBGExp_, 1000.0);
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

            TGTextButton* fitBtn = new TGTextButton(grp, "  Fit Decay  ");
            grp->AddFrame(fitBtn, new TGLayoutHints(kLHintsCenterX, 0, 0, 6, 4));
            fitBtn->Connect("Clicked()", "GammaFitGUI", this, "OnFitDecay()");
        }

        // ── Fit Results ───────────────────────────────────────────────────────
        {
            TGGroupFrame* grp = new TGGroupFrame(p2, "Fit Results");
            p2->AddFrame(grp, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 4, 4, 2, 4));
            decayResultView_ = new TGTextView(grp, 280, 180);
            grp->AddFrame(decayResultView_,
                          new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
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

            TGTextButton* fitTdBtn = new TGTextButton(grp, "  Fit Total Decay  ");
            grp->AddFrame(fitTdBtn, new TGLayoutHints(kLHintsCenterX, 0, 0, 6, 4));
            fitTdBtn->Connect("Clicked()", "GammaFitGUI", this, "OnFitTotalDecay()");
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

void GammaFitGUI::OnDecayModelChanged(Int_t id)
{
    if (!decayEquationLabel_) return;
    static const char* kEq[] = {
        "",
        "f(t) = A*exp(-ln2*t/T_P) + BG",                              // 1 Parent
        "f(t) = A*lD/(lD-lP)*(e^-lP*t - e^-lD*t) + BG",             // 2 Daughter beta-
        "f(t) = A*lP*lD*lG*Sum[e^-li*t/(prod)] + BG",                // 3 GDaughter beta-
        "f(t) = BG  (constant background)",                            // 4 BG only
        "f(t) = Daughter Bateman [beta-n chain] + BG",                   // 5 Daughter beta-n
        "f(t) = Daughter Bateman [beta-2n chain] + BG",                  // 6 Daughter beta-2n
        "f(t) = GDaughter Bateman [beta-n chain] + BG",                  // 7 GDaughter beta-n
        "f(t) = GDaughter Bateman [beta-2n chain] + BG",                 // 8 GDaughter beta-2n
    };
    if (id >= 1 && id <= 8) {
        decayEquationLabel_->SetText(kEq[id]);
        decayEquationLabel_->Layout();
    }

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
}

void GammaFitGUI::OnDecayBGTypeChanged(Int_t /*id*/)
{
    int bgType = decayBGTypeCombo_ ? decayBGTypeCombo_->GetSelected() : 1;
    bool hasExp = (bgType >= 2);
    if (decayThalfBGExp_) decayThalfBGExp_->SetState(hasExp ? kTRUE : kFALSE);
    if (decayFixBGExp_)   decayFixBGExp_->SetEnabled(hasExp ? kTRUE : kFALSE);
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
                    int rebin = decayRebinEntry_ ? (int)decayRebinEntry_->GetNumber() : 1;
                    if (rebin > 1) hDecay->Rebin(rebin);

                    double tLo = std::max(0.0, hDecay->GetXaxis()->GetXmin());
                    double tHi = hDecay->GetXaxis()->GetXmax();
                    hDecay->GetXaxis()->SetRangeUser(tLo, tHi);
                    hDecay->GetXaxis()->SetTitle("Time (ms)");
                    hDecay->GetYaxis()->SetTitle("Counts");
                    hDecay->SetLineColor(kBlack);
                    hDecay->SetMarkerSize(0);

                    int bgType  = (r.bgType >= 1 && r.bgType <= 3) ? r.bgType : 1;
                    int nSig    = DecaySigNpar(sig);
                    TF1* fStored = BuildDecayTF1("fDecayReplay", r.model, bgType, tLo, tHi);
                    if (fStored && (int)r.params.size() >= fStored->GetNpar()) {
                        for (int i = 0; i < fStored->GetNpar(); i++)
                            fStored->SetParameter(i, r.params[i]);
                        fStored->SetLineColor(kRed);
                        fStored->SetLineWidth(2);

                        TCanvas* c = canvas_->GetCanvas();
                        c->Clear(); c->cd();
                        hDecay->Draw(showErrorBars_ ? "hist E" : "hist");
                        fStored->Draw("same");

                        // Signal component (zero all BG params)
                        if (sig != 4) {
                            int npar = fStored->GetNpar();
                            TF1* fSig = (TF1*)fStored->Clone("fDecayReplaySig");
                            for (int i = nSig; i < npar; i++) fSig->FixParameter(i, 0.0);
                            fSig->SetLineColor(kBlue+1);
                            fSig->SetLineWidth(1);
                            fSig->SetLineStyle(2);
                            fSig->Draw("same");

                            // BG component
                            TF1* fBG = (TF1*)fStored->Clone("fDecayReplayBG");
                            if (sig >= 1) fBG->FixParameter(0, 0.0); // zero signal amplitude
                            fBG->SetLineColor(kGreen+2);
                            fBG->SetLineWidth(2);
                            fBG->SetLineStyle(2);
                            fBG->Draw("same");

                            TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.92);
                            leg->SetBorderSize(1);
                            leg->AddEntry(hDecay,  "Data",       "l");
                            leg->AddEntry(fStored, "Total fit",  "l");
                            leg->AddEntry(fSig,    "Signal",     "l");
                            leg->AddEntry(fBG,     "Background", "l");
                            leg->Draw();
                        }

                        c->Modified(); c->Update();
                        gSystem->ProcessEvents();
                    }

                    // Update result view
                    decayResultView_->Clear();
                    decayResultView_->AddLine(
                        Form("Cached fit: E=%.2f keV  window:[%.2f, %.2f]", E, r.eMin, r.eMax));
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
        return;  // done
    }

    // ── No stored fit  -  reset sigma windows to 1.0 and show gamma preview ──────
    if (decaySigRangeEntry_) decaySigRangeEntry_->SetNumber(1.0);
    if (decaySigLoEntry_)    decaySigLoEntry_->SetNumber(1.0);
    OnDecayPreviewGammaPeak();
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

    c->Modified(); c->Update();
    gSystem->ProcessEvents();
    SetStatus(Form("Cut preview: %.4f keV  [%.4f, %.4f]  (lo=%.3f hi=%.3f sig)",
                   E, eMin, eMax, NsigLo, NsigHi));
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
    LoadDecayFitCache();
    AppendLog(Form("[Decay] %d peaks loaded from cache: %s",
                   (int)decayPeakEs_.size(), peakCacheName.c_str()));
}

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
    if (rebin > 1) {
        double binW = hDecay->GetBinWidth(1);
        hDecay->GetYaxis()->SetTitle(Form("Counts / (%.4g ms)", binW));
    } else {
        hDecay->GetYaxis()->SetTitle("Counts");
    }
    hDecay->SetLineColor(kBlack);
    hDecay->SetMarkerSize(0);
    hDecay->Draw(showErrorBars_ ? "hist E" : "hist");
    c->Modified(); c->Update();
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

    TF1* fDecay = BuildDecayTF1("fDecay", modelId, bgType, xlo, xhi);

    // Seed signal parameters
    double initMax = std::max(1.0, hDecay->GetMaximum() - hDecay->GetMinimum());
    double initMin = std::max(0.0, hDecay->GetMinimum());
    double TP  = decayThalfP_  ? decayThalfP_->GetNumber()  : (xhi-xlo)/2.0;
    double TD  = decayThalfD_  ? decayThalfD_->GetNumber()  : (xhi-xlo)/4.0;
    double TG  = decayThalfG_  ? decayThalfG_->GetNumber()  : (xhi-xlo)/8.0;

    if (sigType == 1) {
        fDecay->SetParameter(0, initMax);
        fDecay->SetParameter(1, TP > 0 ? TP : (xhi-xlo)/2.0);
        fDecay->SetParLimits(1, 1e-6, 1e12);
        if (decayFixP_ && decayFixP_->IsOn()) fDecay->FixParameter(1, TP);
    } else if (sigType == 2) {
        fDecay->SetParameter(0, initMax);
        fDecay->SetParameter(1, TP > 0 ? TP : (xhi-xlo)/4.0);
        fDecay->SetParameter(2, TD > 0 ? TD : (xhi-xlo)/2.0);
        fDecay->SetParLimits(1, 1e-6, 1e12);
        fDecay->SetParLimits(2, 1e-6, 1e12);
        if (decayFixP_ && decayFixP_->IsOn()) fDecay->FixParameter(1, TP);
        if (decayFixD_ && decayFixD_->IsOn()) fDecay->FixParameter(2, TD);
    } else if (sigType == 3) {
        fDecay->SetParameter(0, initMax);
        fDecay->SetParameter(1, TP > 0 ? TP : (xhi-xlo)/8.0);
        fDecay->SetParameter(2, TD > 0 ? TD : (xhi-xlo)/4.0);
        fDecay->SetParameter(3, TG > 0 ? TG : (xhi-xlo)/2.0);
        fDecay->SetParLimits(1, 1e-6, 1e12);
        fDecay->SetParLimits(2, 1e-6, 1e12);
        fDecay->SetParLimits(3, 1e-6, 1e12);
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
        fDecay->SetParLimits(nSig+2, 1e-6, 1e12);
        if (decayFixBGExp_ && decayFixBGExp_->IsOn()) fDecay->FixParameter(nSig+2, Tbg);
    } else {
        fDecay->SetParameter(nSig,   initMin);
        fDecay->SetParameter(nSig+1, Tbg > 0 ? Tbg : (xhi-xlo)*2.0);
        fDecay->SetParLimits(nSig+1, 1e-6, 1e12);
        if (decayFixBGExp_ && decayFixBGExp_->IsOn()) fDecay->FixParameter(nSig+1, Tbg);
    }

    fDecay->SetLineColor(kRed);
    fDecay->SetLineWidth(2);

    TFitResultPtr fitRes = hDecay->Fit(fDecay, "SR");

    hDecay->Draw(showErrorBars_ ? "hist E" : "hist");
    fDecay->Draw("same");

    // ── Draw fit components ────────────────────────────────────────────────────
    if (sigType != 4) {
        int npar = fDecay->GetNpar();
        // Signal: zero all BG params
        TF1* fSig = (TF1*)fDecay->Clone("fDecaySigComp");
        for (int i = nSig; i < npar; i++) fSig->FixParameter(i, 0.0);
        fSig->SetLineColor(kBlue+1);
        fSig->SetLineWidth(1);
        fSig->SetLineStyle(2);
        fSig->Draw("same");

        // Background: zero signal amplitude (p[0])
        TF1* fBGComp = (TF1*)fDecay->Clone("fDecayBGComp");
        fBGComp->FixParameter(0, 0.0);
        fBGComp->SetLineColor(kGreen+2);
        fBGComp->SetLineWidth(2);
        fBGComp->SetLineStyle(2);
        fBGComp->Draw("same");

        TLegend* leg = new TLegend(0.60, 0.72, 0.92, 0.92);
        leg->SetBorderSize(1);
        leg->AddEntry(hDecay,   "Data",       "l");
        leg->AddEntry(fDecay,   "Total fit",  "l");
        leg->AddEntry(fSig,     "Signal",     "l");
        leg->AddEntry(fBGComp,  "Background", "l");
        leg->Draw();
    }

    c->Modified(); c->Update();

    // ── Report ────────────────────────────────────────────────────────────────
    decayResultView_->Clear();
    decayResultView_->AddLine(Form("Peak: %.2f keV  [%.2f, %.2f]  cut: -%.2f/+%.2f sig  rebin:%d",
                                   E, eMin, eMax, NsigLo, NsigHi,
                                   decayRebinEntry_ ? (int)decayRebinEntry_->GetNumber() : 1));
    double chi2ndf = (fitRes.Get() && fitRes->Ndf() > 0)
                     ? fitRes->Chi2() / fitRes->Ndf() : -1.0;
    decayResultView_->AddLine(Form("Chi2/NDF: %.4g  status: %d",
                                   chi2ndf, fitRes.Get() ? fitRes->Status() : -1));
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
        for (int i = 0; i < np; i++) {
            r.params[i] = fDecay->GetParameter(i);
            r.errors[i] = fDecay->GetParError(i);
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
        int bgt = 1;
        if (ss >> bgt) r.bgType = bgt;   // optional; defaults 1 (Flat)
        double nsiglo = -1.0;
        if (ss >> nsiglo) r.NsigLo = nsiglo;
        int rebin = 1;
        if (ss >> rebin) r.rebin = rebin;
        decayFitStore_[r.peakE] = r;
        ++count;
    }
    AppendLog(Form("[Decay] Loaded %d cached decay fit(s) from %s", count, path.c_str()));
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

    // Rebin before display
    int rebin = decayRebinEntry_ ? (int)decayRebinEntry_->GetNumber() : 1;
    if (rebin > 1) hDecay->Rebin(rebin);

    // Restrict to t > 0
    double tLo = std::max(0.0, hDecay->GetXaxis()->GetXmin());
    hDecay->GetXaxis()->SetRangeUser(tLo, hDecay->GetXaxis()->GetXmax());

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
    hDecay->Draw(showErrorBars_ ? "hist E" : "hist");
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
    if (!decayTdEquationLbl_) return;
    static const char* kEq[] = {
        "",
        "f(t) = A*exp(-ln2*t/T_P) + BG",
        "f(t) = A*lD/(lD-lP)*(e^-lP*t - e^-lD*t) + BG",
        "f(t) = A*lP*lD*lG*Sum[e^-li*t/(prod)] + BG",
        "f(t) = BG  (background only)",
        "f(t) = Daughter Bateman [beta-n chain] + BG",
        "f(t) = Daughter Bateman [beta-2n chain] + BG",
        "f(t) = GDaughter Bateman [beta-n chain] + BG",
        "f(t) = GDaughter Bateman [beta-2n chain] + BG",
    };
    if (id >= 1 && id <= 8) {
        decayTdEquationLbl_->SetText(kEq[id]);
        decayTdEquationLbl_->Layout();
    }
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
}

void GammaFitGUI::OnDecayTdBGTypeChanged(Int_t /*id*/)
{
    int bgType = decayTdBGCombo_ ? decayTdBGCombo_->GetSelected() : 1;
    bool hasExp = (bgType >= 2);
    if (decayTdThalfBGExp_) decayTdThalfBGExp_->SetState(hasExp ? kTRUE : kFALSE);
    if (decayTdFixBGExp_)   decayTdFixBGExp_->SetEnabled(hasExp ? kTRUE : kFALSE);
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

    // Rebin (use same rebin as Cuts tab)
    int rebin = decayRebinEntry_ ? (int)decayRebinEntry_->GetNumber() : 1;
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

    TFitResultPtr fitRes = hTotal->Fit(fDecay, "SR");

    hTotal->Draw(showErrorBars_ ? "hist E" : "hist");
    fDecay->Draw("same");

    // ── Draw components ───────────────────────────────────────────────────────
    TLegend* leg = new TLegend(0.58, 0.68, 0.92, 0.92);
    leg->SetBorderSize(1);
    leg->AddEntry(hTotal, "Data (sum)", "l");
    leg->AddEntry(fDecay, "Total fit",  "l");
    if (sigType != 4) {
        int npar = fDecay->GetNpar();
        TF1* fSig = (TF1*)fDecay->Clone("fTdSig");
        for (int i = nSig; i < npar; i++) fSig->FixParameter(i, 0.0);
        fSig->SetLineColor(kBlue+1); fSig->SetLineWidth(1); fSig->SetLineStyle(2);
        fSig->Draw("same");
        leg->AddEntry(fSig, "Signal", "l");

        TF1* fBGComp = (TF1*)fDecay->Clone("fTdBG");
        fBGComp->FixParameter(0, 0.0);
        fBGComp->SetLineColor(kGreen+2); fBGComp->SetLineWidth(2); fBGComp->SetLineStyle(2);
        fBGComp->Draw("same");
        leg->AddEntry(fBGComp, "Background", "l");
    }
    leg->Draw();
    c->Modified(); c->Update();

    // ── Report ────────────────────────────────────────────────────────────────
    double chi2ndf = (fitRes.Get() && fitRes->Ndf() > 0)
                     ? fitRes->Chi2() / fitRes->Ndf() : -1.0;
    if (decayTdResultView_) decayTdResultView_->Clear();
    if (decayTdResultView_) {
        decayTdResultView_->AddLine(Form("Total decay: %d peaks summed  TH2: %s",
                                        (int)peakIndices.size(), decayTh2Name_.c_str()));
        decayTdResultView_->AddLine(Form("Peaks: %s", summedPeaks.c_str()));
        decayTdResultView_->AddLine(Form("Rebin: %d   Chi2/NDF: %.4g   status: %d",
                                        rebin, chi2ndf, fitRes.Get() ? fitRes->Status() : -1));
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
// BuildIsotopesTab
// ─────────────────────────────────────────────────────────────────────────────

