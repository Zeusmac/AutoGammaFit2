#include "GammaFitGUI.h"
#include "GammaFitGUI_shared.h"

#include "TGTextEntry.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TGFileDialog.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLatex.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TMath.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <sys/stat.h>

// Ge Fano-limit constants: FWHM_stat = kFanoCoeff * sqrt(kFanoF * kFanoW * E)
static constexpr double kFanoF     = 0.12;       // Ge Fano factor
static constexpr double kFanoW     = 0.00296;    // keV  (= 2.96 eV per e-h pair)
static constexpr double kFanoCoeff = 2.344;
static constexpr double kFw        = kFanoF * kFanoW;  // 3.552e-4 keV

void GammaFitGUI::BuildFWHMTab(TGCompositeFrame* p)
{
    // ── Histogram selector ────────────────────────────────────────────────────
    TGGroupFrame* hg = new TGGroupFrame(p, "Histograms (combined)");
    p->AddFrame(hg, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

    {
        TGHorizontalFrame* comboRow = new TGHorizontalFrame(hg);
        hg->AddFrame(comboRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        fwhmCombo_ = new TGComboBox(comboRow, 800);
        fwhmCombo_->Resize(170, 22);
        comboRow->AddFrame(fwhmCombo_,
                           new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 0, 4, 0, 0));

        TGTextButton* addBtn = new TGTextButton(comboRow, " Add to Plot ");
        comboRow->AddFrame(addBtn, new TGLayoutHints(kLHintsCenterY, 0, 0, 0, 0));
        addBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLoadFWHM()");
        addBtn->SetToolTipText(
            "Append FWHM points from the selected histogram to the combined plot.\n"
            "You can add as many histograms as you like; each gets a distinct color.");
    }

    // Loaded-histogram list with Remove / Clear All
    fwhmHistList_ = new TGListBox(hg, -1);
    fwhmHistList_->Resize(290, 58);
    hg->AddFrame(fwhmHistList_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

    {
        TGHorizontalFrame* btnRow = new TGHorizontalFrame(hg);
        hg->AddFrame(btnRow, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 2));

        TGTextButton* remBtn = new TGTextButton(btnRow, " Remove ");
        btnRow->AddFrame(remBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        remBtn->Connect("Clicked()", "GammaFitGUI", this, "OnRemoveFWHMHist()");
        remBtn->SetToolTipText("Remove the selected histogram and its points from the combined plot.");

        TGTextButton* clrBtn = new TGTextButton(btnRow, " Clear All ");
        btnRow->AddFrame(clrBtn, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
        clrBtn->Connect("Clicked()", "GammaFitGUI", this, "OnClearFWHMHists()");
        clrBtn->SetToolTipText("Remove all histograms and data points from the combined plot.");
    }

    TGLabel* formulaLbl = new TGLabel(hg, "  FWHM^2 = a + b*E + c*E^2");
    hg->AddFrame(formulaLbl, new TGLayoutHints(kLHintsLeft, 4, 4, 2, 4));

    // ── Model parameters ──────────────────────────────────────────────────────
    TGGroupFrame* pg = new TGGroupFrame(p, "Resolution Model Parameters");
    p->AddFrame(pg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    // Helper: one row showing  "Label   [value]   lo:[lo]  to  [hi]"
    auto addParRow = [&](TGGroupFrame* grp, const char* lbl,
                         TGNumberEntry*& val, double vDef,
                         TGNumberEntry*& lo,  double loDef,
                         TGNumberEntry*& hi,  double hiDef)
    {
        TGHorizontalFrame* row = new TGHorizontalFrame(grp);
        grp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 1));
        TGLabel* l = new TGLabel(row, lbl);
        l->SetWidth(64);
        row->AddFrame(l, new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));

        val = new TGNumberEntry(row, vDef, 9, -1,
            TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELNoLimits);
        row->AddFrame(val, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));

        row->AddFrame(new TGLabel(row, "lo:"),
                      new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        lo = new TGNumberEntry(row, loDef, 7, -1,
            TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELNoLimits);
        row->AddFrame(lo, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));

        row->AddFrame(new TGLabel(row, "hi:"),
                      new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        hi = new TGNumberEntry(row, hiDef, 7, -1,
            TGNumberFormat::kNESReal, TGNumberFormat::kNEAAnyNumber,
            TGNumberFormat::kNELNoLimits);
        row->AddFrame(hi, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
    };

    addParRow(pg, "a (keV^2)", mFwhmA_, 0.5,   mFwhmAlo_, 1e-4,  mFwhmAhi_, 1000.0);
    addParRow(pg, "b (keV)",   mFwhmB_, 0.02,  mFwhmBlo_, 0.0,   mFwhmBhi_,   10.0);
    addParRow(pg, "c",         mFwhmC_, 0.0,   mFwhmClo_, 0.0,   mFwhmChi_,    0.1);

    TGTextButton* fitBtn = new TGTextButton(pg, "Fit Model to Data");
    pg->AddFrame(fitBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 6, 2));
    fitBtn->Connect("Clicked()", "GammaFitGUI", this, "OnFitFWHM()");
    fitBtn->SetToolTipText("Fit FWHM^2 = a + b*E + c*E^2 directly to the included (non-excluded) data points.\nFitting in FWHM^2 space is linear, giving more numerically stable convergence.");

    fwhmResultLbl_ = new TGLabel(pg, "No fit yet");
    pg->AddFrame(fwhmResultLbl_, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

    // ── Display options ───────────────────────────────────────────────────────
    TGGroupFrame* dg = new TGGroupFrame(p, "Display");
    p->AddFrame(dg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    fwhmStatLineChk_ = new TGCheckButton(dg,
        "Show Fano limit  [2.344*sqrt(F*w*E),  F=0.12,  w=2.96 eV (Ge)]");
    fwhmStatLineChk_->SetState(kButtonDown);
    dg->AddFrame(fwhmStatLineChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 1));
    fwhmStatLineChk_->Connect("Clicked()", "GammaFitGUI", this, "OnFWHMRedisplay()");
    fwhmStatLineChk_->SetToolTipText(
        "Draw the fundamental Fano statistical limit for a germanium detector:\n"
        "  FWHM_stat = 2.344 * sqrt(F * w * E)\n"
        "  F = 0.12  (Ge Fano factor)\n"
        "  w = 2.96 eV  (mean energy per electron-hole pair in Ge)\n"
        "This is the best resolution physically achievable — limited only by\n"
        "charge-carrier counting statistics.  Your actual FWHM exceeds this\n"
        "due to electronic noise (a term) and charge-trapping (c term).");

    fwhmShowResChk_ = new TGCheckButton(dg, "Show Resolution (%)  =  100 x FWHM / E");
    dg->AddFrame(fwhmShowResChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 1, 1));
    fwhmShowResChk_->Connect("Clicked()", "GammaFitGUI", this, "OnFWHMRedisplay()");
    fwhmShowResChk_->SetToolTipText(
        "Switch Y axis to Resolution (%) = 100 * FWHM / Energy.\n"
        "Statistical limit becomes 100*sqrt(b/E) = 2.355*sqrt(F*eps/E)*100");

    fwhmShowSigmaChk_ = new TGCheckButton(dg, "Show as sigma  (sqrt(FWHM^2) / 2.3548)");
    dg->AddFrame(fwhmShowSigmaChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 1, 4));
    fwhmShowSigmaChk_->Connect("Clicked()", "GammaFitGUI", this, "OnFWHMRedisplay()");
    fwhmShowSigmaChk_->SetToolTipText(
        "Switch Y axis from FWHM^2 (keV^2) to Gaussian sigma (keV).\n"
        "Has no effect when Resolution (%) mode is active.");

    // ── Point removal ─────────────────────────────────────────────────────────
    TGGroupFrame* rg = new TGGroupFrame(p, "Point Selection");
    p->AddFrame(rg, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    fwhmRemoveModeChk_ = new TGCheckButton(rg, "Click canvas to exclude/restore points");
    rg->AddFrame(fwhmRemoveModeChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));
    fwhmRemoveModeChk_->SetToolTipText(
        "When checked: click near a data point to toggle it excluded.\n"
        "Excluded points (gray circles) are ignored during fitting.");

    {
        TGHorizontalFrame* btnRow = new TGHorizontalFrame(rg);
        rg->AddFrame(btnRow, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 4));

        TGTextButton* restoreBtn = new TGTextButton(btnRow, " Restore All ");
        btnRow->AddFrame(restoreBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        restoreBtn->Connect("Clicked()", "GammaFitGUI", this, "OnFWHMRestoreAll()");
        restoreBtn->SetToolTipText("Re-include all excluded data points");
    }

    TGLabel* hint = new TGLabel(rg, "  Blue filled = included   Gray open = excluded");
    rg->AddFrame(hint, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 2));

    // ── Actions ───────────────────────────────────────────────────────────────
    TGGroupFrame* ag = new TGGroupFrame(p, "Actions");
    p->AddFrame(ag, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

    TGTextButton* accBtn = new TGTextButton(ag, "Accept & Save to Cache");
    ag->AddFrame(accBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    accBtn->Connect("Clicked()", "GammaFitGUI", this, "OnAcceptFWHM()");
    accBtn->SetToolTipText(
        "Save the fitted resolution model (and excluded point list) to cache.\n"
        "AutoFit will use the model for peak grouping and sigma seeding.");

    TGTextButton* saveCanvasBtn = new TGTextButton(ag, "Save Canvas to ROOT File");
    ag->AddFrame(saveCanvasBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
    saveCanvasBtn->Connect("Clicked()", "GammaFitGUI", this, "OnSaveFWHMCanvas()");
    saveCanvasBtn->SetToolTipText(
        "Save the current FWHM canvas to a ROOT file.\n"
        "Opens a file dialog; choose or type a .root filename.");
}

// ─────────────────────────────────────────────────────────────────────────────
// FWHM tab — load, fit, accept, display helpers
// ─────────────────────────────────────────────────────────────────────────────

// Color palette and marker styles for per-histogram coloring (up to 8 sources).
// Colors chosen for visibility on white background; marker styles vary so the
// plot is still readable in black-and-white print.
static const int kFWHMPalette[]    = { kBlue+1, kRed, kGreen+2, kMagenta+1,
                                        kOrange+7, kCyan+2, kViolet+2, kPink+1 };
static const int kFWHMMarkers[]    = { 20, 21, 22, 23, 33, 29, 34, 20 };
static const int kNFWHMPalette     = 8;

// Draw the FWHM plot onto any canvas.  All objects are heap-allocated and the
// canvas takes ownership via Draw() — do NOT delete them afterwards.
// fwhmTF1_ is NEVER drawn directly; fresh TF1 expressions are created here so
// that fwhmTF1_ stays fully under our control and safe to access later.
void GammaFitGUI::DrawFWHMToCanvas(TCanvas* c, bool showSigma, bool showStatLine,
                                    bool showResolution, bool showErrBars)
{
    if (fwhmAllX_.empty()) return;

    const double kSig = 2.3548200450309493;
    bool multiSource = (fwhmLoadedHists_.size() > 1);

    // Resolution (%) mode overrides sigma mode
    std::string yLabel, title;
    if (showResolution) {
        yLabel = "Resolution (%)";
        title  = "Resolution vs Energy  [" + fwhmHistName_ + "];Energy (keV);Resolution (%)";
    } else if (showSigma) {
        yLabel = "#sigma (keV)";
        title  = "#sigma vs Energy  [" + fwhmHistName_ + "];Energy (keV);#sigma (keV)";
    } else {
        yLabel = "FWHM^{2} (keV^{2})";
        title  = "FWHM^{2} vs Energy  [" + fwhmHistName_ + "];Energy (keV);FWHM^{2} (keV^{2})";
    }

    double xhi = *std::max_element(fwhmAllX_.begin(), fwhmAllX_.end()) * 1.15;

    auto toY = [&](double fwhm2, double E) -> double {
        if (showResolution) return 100.0 * std::sqrt(fwhm2) / E;
        if (showSigma)      return std::sqrt(fwhm2) / kSig;
        return fwhm2;
    };

    // ── Build per-source point sets (included) + global excluded/tied buckets ─
    // Excluded and tied points are always gray/orange regardless of source so
    // the included→excluded toggle is immediately obvious.
    struct SrcPoints { std::vector<double> x, y, ex, ey; };
    std::map<std::string, SrcPoints> srcIncl;
    std::vector<double> xEx, yEx;
    std::vector<double> xTied, yTied;

    for (size_t i = 0; i < fwhmAllX_.size(); i++) {
        double y = toY(fwhmAllY_[i], fwhmAllX_[i]);
        bool tied = (!fwhmTied_.empty() && fwhmTied_[i]);
        const std::string& src = fwhmHistSources_.empty() ? fwhmHistName_
                                                           : fwhmHistSources_[i];
        if (tied) {
            xTied.push_back(fwhmAllX_[i]); yTied.push_back(y);
        } else if (fwhmExcluded_[i]) {
            xEx.push_back(fwhmAllX_[i]); yEx.push_back(y);
        } else {
            auto& sp = srcIncl[src];
            sp.x.push_back(fwhmAllX_[i]); sp.y.push_back(y);
            sp.ex.push_back(0.0);          sp.ey.push_back(0.05 * y);
        }
    }

    // ── Compute y-axis range ──────────────────────────────────────────────────
    double yDataMax = 0.0;
    for (auto& [s, sp] : srcIncl) for (double v : sp.y) yDataMax = std::max(yDataMax, v);
    for (double v : yEx)    yDataMax = std::max(yDataMax, v);
    for (double v : yTied)  yDataMax = std::max(yDataMax, v);
    if (yDataMax <= 0) yDataMax = 1.0;
    double yAxisMax = yDataMax * 1.35;

    // ── Draw frame to establish axes ──────────────────────────────────────────
    {
        double xmin = *std::min_element(fwhmAllX_.begin(), fwhmAllX_.end());
        double frameXlo = showResolution ? 1.0 : std::max(0.0, xmin * 0.85);
        TH1F* frame = static_cast<TH1F*>(c->DrawFrame(frameXlo, 0.0, xhi, yAxisMax));
        frame->SetTitle(title.c_str());
        frame->GetXaxis()->SetTitle("Energy (keV)");
        frame->GetXaxis()->SetTitleSize(0.05);
        frame->GetXaxis()->SetLabelSize(0.04);
        frame->GetYaxis()->SetTitle(yLabel.c_str());
        frame->GetYaxis()->SetTitleSize(0.05);
        frame->GetYaxis()->SetLabelSize(0.04);
        frame->GetYaxis()->SetTitleOffset(1.2);
    }

    // ── Included points — one TGraphErrors per source histogram ──────────────
    // Order of drawing follows fwhmLoadedHists_ so colors are stable.
    const std::vector<std::string>& order =
        fwhmLoadedHists_.empty() ? std::vector<std::string>{fwhmHistName_}
                                 : fwhmLoadedHists_;
    for (size_t si = 0; si < order.size(); ++si) {
        const std::string& src = order[si];
        auto it = srcIncl.find(src);
        if (it == srcIncl.end() || it->second.x.empty()) continue;
        auto& sp = it->second;
        int col = kFWHMPalette[si % kNFWHMPalette];
        int mst = kFWHMMarkers[si % kNFWHMPalette];
        TGraphErrors* gr = new TGraphErrors((Int_t)sp.x.size(),
                                             sp.x.data(), sp.y.data(),
                                             sp.ex.data(), sp.ey.data());
        gr->SetMarkerStyle(mst); gr->SetMarkerSize(0.9);
        gr->SetMarkerColor(col); gr->SetLineColor(col);
        gr->Draw(showErrBars ? "EP same" : "P same");
    }

    // ── Excluded points (gray hollow circles, all sources) ───────────────────
    if (!xEx.empty()) {
        TGraph* grEx = new TGraph((Int_t)xEx.size(), xEx.data(), yEx.data());
        grEx->SetMarkerStyle(24); grEx->SetMarkerSize(0.9);
        grEx->SetMarkerColor(kGray + 1);
        grEx->Draw("P same");
    }

    // ── Tied-width points (orange hollow triangles, all sources) ─────────────
    if (!xTied.empty()) {
        TGraph* grTied = new TGraph((Int_t)xTied.size(), xTied.data(), yTied.data());
        grTied->SetMarkerStyle(26); grTied->SetMarkerSize(1.1);
        grTied->SetMarkerColor(kOrange + 7); grTied->SetLineColor(kOrange + 7);
        grTied->Draw("P same");
    }

    // ── Fano statistical limit ────────────────────────────────────────────────
    double xlo = showResolution ? 1.0 : 0.0;
    if (showStatLine) {
        TF1* statLine = nullptr;
        if (showResolution) {
            statLine = new TF1("fwhm_stat",
                Form("100.0*%.10g*sqrt(%.10g/x)", kFanoCoeff, kFw), xlo, xhi);
        } else if (showSigma) {
            statLine = new TF1("fwhm_stat",
                Form("%.10g*sqrt(%.10g*x)", kFanoCoeff / kSig, kFw), xlo, xhi);
        } else {
            statLine = new TF1("fwhm_stat",
                Form("%.10g*x", kFanoCoeff * kFanoCoeff * kFw), xlo, xhi);
        }
        statLine->SetLineColor(kGreen + 2); statLine->SetLineStyle(2);
        statLine->SetLineWidth(2); statLine->Draw("same");
    }

    // ── Model curve ───────────────────────────────────────────────────────────
    if (fwhmTF1_) {
        double a  = fwhmTF1_->GetParameter(0);
        double b  = fwhmTF1_->GetParameter(1);
        double cv = fwhmTF1_->GetParameter(2);

        TF1* modelDraw = nullptr;
        if (showResolution) {
            modelDraw = new TF1("fwhm_draw",
                Form("100.0*sqrt(%.10g+%.10g*x+%.10g*x*x)/x", a, b, cv), 1.0, xhi);
        } else if (showSigma) {
            modelDraw = new TF1("fwhm_draw",
                Form("sqrt((%.10g+%.10g*x+%.10g*x*x)/%.10g)",
                     a, b, cv, kSig * kSig), 0.0, xhi);
        } else {
            modelDraw = new TF1("fwhm_draw",
                Form("%.10g+%.10g*x+%.10g*x*x", a, b, cv), 0.0, xhi);
        }
        modelDraw->SetLineColor(kRed); modelDraw->SetLineWidth(2);
        modelDraw->Draw("same");

        // ── Fit info pave ─────────────────────────────────────────────────────
        TPaveText* pt = new TPaveText(0.58, 0.62, 0.93, 0.91, "NDC");
        pt->SetFillColor(0); pt->SetFillStyle(1001); pt->SetBorderSize(1);
        pt->SetTextSize(0.024); pt->SetTextAlign(12);
        if (showResolution)
            pt->AddText("R(%) = 100 #frac{#sqrt{a + bE + cE^{2}}}{E}");
        else if (showSigma)
            pt->AddText("#sigma = #frac{#sqrt{a + bE + cE^{2}}}{2.355}");
        else
            pt->AddText("FWHM^{2} = a + bE + cE^{2}");
        pt->AddText(Form("a = %.4g keV^{2}", a));
        pt->AddText(Form("b = %.4g keV",     b));
        pt->AddText(Form("c = %.4g",          cv));
        if (fwhmChi2Ndf_ >= 0)
            pt->AddText(Form("#chi^{2}/ndf = %.2f  (ndf = %d)", fwhmChi2Ndf_, fwhmNdf_));
        if (fwhmPValue_ >= 0)
            pt->AddText(Form("p-value = %.4f", fwhmPValue_));
        if (fwhmResidRMS_ >= 0)
            pt->AddText(Form("RMS_{res} = %.3f", fwhmResidRMS_));
        pt->Draw();
    }

    // ── Legend: always shown when multiple sources or tied points present ─────
    bool hasTied = !xTied.empty();
    bool hasExcl = !xEx.empty();
    if (multiSource || hasTied) {
        // Height: one row per loaded histogram + excluded row + tied row (if present)
        int nRows = (int)order.size() + (hasExcl ? 1 : 0) + (hasTied ? 1 : 0);
        double legH = 0.04 * nRows + 0.02;
        double legY2 = 0.91;
        TLegend* leg = new TLegend(0.11, legY2 - legH, 0.52, legY2);
        leg->SetBorderSize(1); leg->SetFillStyle(1001); leg->SetTextSize(0.026);
        for (size_t si = 0; si < order.size(); ++si) {
            int col = kFWHMPalette[si % kNFWHMPalette];
            int mst = kFWHMMarkers[si % kNFWHMPalette];
            TGraph* gL = new TGraph(1);
            gL->SetMarkerStyle(mst); gL->SetMarkerColor(col); gL->SetMarkerSize(0.9);
            leg->AddEntry(gL, order[si].c_str(), "p");
        }
        if (hasExcl) {
            TGraph* gEx = new TGraph(1);
            gEx->SetMarkerStyle(24); gEx->SetMarkerColor(kGray+1); gEx->SetMarkerSize(0.9);
            leg->AddEntry(gEx, "Excluded", "p");
        }
        if (hasTied) {
            TGraph* gTd = new TGraph(1);
            gTd->SetMarkerStyle(26); gTd->SetMarkerColor(kOrange+7); gTd->SetMarkerSize(1.1);
            leg->AddEntry(gTd, "Tied-width", "p");
        }
        leg->Draw();
    } else if (hasTied) {
        // Single-source with tied points — minimal legend
        TLegend* leg = new TLegend(0.12, 0.80, 0.45, 0.91);
        leg->SetBorderSize(1); leg->SetFillStyle(1001); leg->SetTextSize(0.030);
        TGraph* gL1 = new TGraph(1); gL1->SetMarkerStyle(20); gL1->SetMarkerColor(kBlue+1); gL1->SetMarkerSize(0.9);
        TGraph* gL2 = new TGraph(1); gL2->SetMarkerStyle(24); gL2->SetMarkerColor(kGray+1); gL2->SetMarkerSize(0.9);
        TGraph* gL3 = new TGraph(1); gL3->SetMarkerStyle(26); gL3->SetMarkerColor(kOrange+7); gL3->SetMarkerSize(1.1);
        leg->AddEntry(gL1, "Free-width",  "p");
        leg->AddEntry(gL2, "Excluded",    "p");
        leg->AddEntry(gL3, "Tied-width",  "p");
        leg->Draw();
    }
}

void GammaFitGUI::RedrawFWHM()
{
    if (fwhmAllX_.empty()) return;
    bool showSigma    = fwhmShowSigmaChk_ && fwhmShowSigmaChk_->IsOn();
    bool showStatLine = fwhmStatLineChk_  && fwhmStatLineChk_->IsOn();
    bool showRes      = fwhmShowResChk_   && fwhmShowResChk_->IsOn();
    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();
    DrawFWHMToCanvas(c, showSigma, showStatLine, showRes, showErrorBars_);
    c->Modified(); c->Update();
}

void GammaFitGUI::OnFWHMRedisplay()
{
    RedrawFWHM();
}

void GammaFitGUI::OnFWHMRestoreAll()
{
    if (fwhmExcluded_.empty()) { AppendLog("No FWHM data loaded."); return; }
    int n = 0;
    for (size_t i = 0; i < fwhmExcluded_.size(); i++) { if (fwhmExcluded_[i]) { fwhmExcluded_[i] = false; ++n; } }
    AppendLog("Restored " + std::to_string(n) + " excluded FWHM point(s).");
    RedrawFWHM();
}

void GammaFitGUI::OnLoadFWHM()
{
    Int_t id = fwhmCombo_->GetSelected();
    if (id < 1 || (size_t)id > histNames_.size()) {
        AppendLog("Select a histogram from the FWHM dropdown first.");
        return;
    }
    std::string hname = histNames_[id - 1];

    // Reject duplicates
    for (const auto& h : fwhmLoadedHists_) {
        if (h == hname) {
            AppendLog("FWHM: " + hname + " is already in the combined plot.");
            return;
        }
    }

    FitDatabase fitdb;
    if (!fitdb.Load(CacheFileFor(hname))) {
        AppendLog("No cache found for " + hname + " — run AutoFit first.");
        return;
    }

    // Remember where new points start so we can apply this histogram's exclusions
    size_t insertStart = fwhmAllX_.size();

    for (const auto& [key, entry] : fitdb.GetEntries()) {
        if (key.size() >= 2 && key[0] == '_' && key[1] == '_') continue;
        FitLayout lay = DetectLayout((int)entry.params.size());
        if (!lay.valid()) continue;
        for (int i = 0; i < lay.n; i++) {
            double E   = entry.params[3*i + 1];
            double sig = entry.params[3*i + 2];
            if (sig <= 0.0 || E <= 0.0) continue;
            double fwhm = 2.3548 * sig;
            fwhmAllX_.push_back(E);
            fwhmAllY_.push_back(fwhm * fwhm);
            fwhmTied_.push_back(entry.widthTied);
            fwhmExcluded_.push_back(entry.widthTied || entry.needsRefit);
            fwhmHistSources_.push_back(hname);
        }
    }

    if (fwhmAllX_.size() == insertStart) {
        AppendLog("No fitted peaks in cache for " + hname +
                  " — run AutoFit or accept manual fits first.");
        return;
    }

    // Restore saved per-histogram exclusions (search only the new block)
    auto exIt = fitdb.GetEntries().find(kExcludedFwhmKey);
    if (exIt != fitdb.GetEntries().end()) {
        for (const double exE : exIt->second.params) {
            for (size_t i = insertStart; i < fwhmAllX_.size(); i++) {
                if (std::abs(fwhmAllX_[i] - exE) < 0.5) { fwhmExcluded_[i] = true; break; }
            }
        }
    }
    auto inIt = fitdb.GetEntries().find(kFwhmInclKey);
    if (inIt != fitdb.GetEntries().end()) {
        for (const double inE : inIt->second.params) {
            for (size_t i = insertStart; i < fwhmAllX_.size(); i++) {
                if (std::abs(fwhmAllX_[i] - inE) < 0.5) { fwhmExcluded_[i] = false; break; }
            }
        }
    }

    // Seed resolution params from the first histogram's cache (only once)
    if (fwhmLoadedHists_.empty()) {
        delete fwhmTF1_; fwhmTF1_ = nullptr;
        fwhmChi2Ndf_ = -1.0; fwhmPValue_ = -1.0; fwhmResidRMS_ = -1.0; fwhmNdf_ = 0;
        auto rit = fitdb.GetEntries().find(kResolutionKey);
        if (rit != fitdb.GetEntries().end() && rit->second.params.size() == 3) {
            mFwhmA_->SetNumber(rit->second.params[0]);
            mFwhmB_->SetNumber(rit->second.params[1]);
            mFwhmC_->SetNumber(rit->second.params[2]);
            double xhi = *std::max_element(fwhmAllX_.begin(), fwhmAllX_.end()) * 1.15;
            fwhmTF1_ = new TF1("fwhm_model", "[0]+[1]*x+[2]*x*x", 0.0, xhi);
            fwhmTF1_->SetParameters(rit->second.params[0],
                                     rit->second.params[1],
                                     rit->second.params[2]);
            double cachedChi2 = rit->second.chi2ndf;
            if (cachedChi2 > 0 && cachedChi2 < 1e10) fwhmChi2Ndf_ = cachedChi2;
        }
    }

    fwhmLoadedHists_.push_back(hname);
    fwhmHistName_ = (fwhmLoadedHists_.size() == 1) ? hname : "Combined";

    // Update the listbox
    if (fwhmHistList_) {
        int nPts = 0;
        for (const auto& s : fwhmHistSources_) if (s == hname) ++nPts;
        std::string entry = hname + "  (" + std::to_string(nPts) + " pts)";
        fwhmHistList_->AddEntry(entry.c_str(), (Int_t)fwhmLoadedHists_.size());
        fwhmHistList_->MapSubwindows(); fwhmHistList_->Layout();
    }

    int nNew  = (int)(fwhmAllX_.size() - insertStart);
    int nExcl = 0;
    for (size_t i = insertStart; i < fwhmAllX_.size(); i++) if (fwhmExcluded_[i]) ++nExcl;
    AppendLog("FWHM: added " + std::to_string(nNew) + " pts from " + hname +
              "  (" + std::to_string(nExcl) + " excluded)  total=" +
              std::to_string(fwhmAllX_.size()));
    SetStatus("FWHM: " + fwhmHistName_);
    RedrawFWHM();
}

void GammaFitGUI::OnRemoveFWHMHist()
{
    if (!fwhmHistList_) return;
    Int_t sel = fwhmHistList_->GetSelected();
    if (sel < 1 || (size_t)sel > fwhmLoadedHists_.size()) {
        AppendLog("Select a histogram in the list first."); return;
    }
    std::string hname = fwhmLoadedHists_[sel - 1];

    // Erase all points belonging to this histogram
    int nRemoved = 0;
    for (int i = (int)fwhmAllX_.size() - 1; i >= 0; --i) {
        if (fwhmHistSources_[i] == hname) {
            fwhmAllX_      .erase(fwhmAllX_      .begin() + i);
            fwhmAllY_      .erase(fwhmAllY_      .begin() + i);
            fwhmExcluded_  .erase(fwhmExcluded_  .begin() + i);
            fwhmTied_      .erase(fwhmTied_      .begin() + i);
            fwhmHistSources_.erase(fwhmHistSources_.begin() + i);
            ++nRemoved;
        }
    }
    fwhmLoadedHists_.erase(fwhmLoadedHists_.begin() + (sel - 1));

    // Rebuild listbox
    fwhmHistList_->RemoveAll();
    for (size_t i = 0; i < fwhmLoadedHists_.size(); ++i) {
        const std::string& h = fwhmLoadedHists_[i];
        int nPts = 0;
        for (const auto& s : fwhmHistSources_) if (s == h) ++nPts;
        fwhmHistList_->AddEntry((h + "  (" + std::to_string(nPts) + " pts)").c_str(), (Int_t)(i + 1));
    }
    fwhmHistList_->MapSubwindows(); fwhmHistList_->Layout();

    if      (fwhmLoadedHists_.empty())          fwhmHistName_ = "";
    else if (fwhmLoadedHists_.size() == 1)      fwhmHistName_ = fwhmLoadedHists_[0];
    else                                         fwhmHistName_ = "Combined";

    AppendLog("FWHM: removed " + std::to_string(nRemoved) + " pts from " + hname +
              "  remaining=" + std::to_string(fwhmAllX_.size()));
    if (!fwhmAllX_.empty()) {
        RedrawFWHM();
    } else {
        TCanvas* c = canvas_->GetCanvas(); c->Clear(); c->Modified(); c->Update();
    }
}

void GammaFitGUI::OnClearFWHMHists()
{
    fwhmAllX_.clear(); fwhmAllY_.clear();
    fwhmExcluded_.clear(); fwhmTied_.clear();
    fwhmHistSources_.clear(); fwhmLoadedHists_.clear();
    fwhmHistName_ = "";
    delete fwhmTF1_; fwhmTF1_ = nullptr;
    fwhmChi2Ndf_ = -1.0; fwhmPValue_ = -1.0; fwhmResidRMS_ = -1.0; fwhmNdf_ = 0;
    if (fwhmHistList_) { fwhmHistList_->RemoveAll(); fwhmHistList_->MapSubwindows(); fwhmHistList_->Layout(); }
    if (fwhmResultLbl_) fwhmResultLbl_->SetText("No fit yet");
    TCanvas* c = canvas_->GetCanvas(); c->Clear(); c->Modified(); c->Update();
    AppendLog("FWHM: cleared all histograms.");
}

void GammaFitGUI::OnFitFWHM()
{
    if (fwhmAllX_.empty()) {
        AppendLog("Load FWHM data first.");
        return;
    }

    // Build a graph from non-excluded points only.
    // fwhmAllY_ holds FWHM² (keV²); error = 10% of FWHM² (propagated from 5% FWHM).
    std::vector<double> xIn, yIn, exIn, eyIn;
    for (size_t i = 0; i < fwhmAllX_.size(); i++) {
        if (fwhmExcluded_[i]) continue;
        double y = fwhmAllY_[i];  // FWHM²
        xIn.push_back(fwhmAllX_[i]); yIn.push_back(y);
        exIn.push_back(0.0); eyIn.push_back(0.10 * y);  // δ(FWHM²) = 2·FWHM·δFWHM ≈ 0.10·FWHM²
    }
    if ((int)xIn.size() < 3) {
        AppendLog("Need at least 3 included points to fit (currently " +
                  std::to_string(xIn.size()) + ").");
        return;
    }

    double xhi = *std::max_element(xIn.begin(), xIn.end()) * 1.15;

    // Auto-seed from data so the fit starts close to the correct answer.
    // The dominant term is b*E (charge-carrier statistics): b ≈ FWHM²/E.
    // Take the median across all included points so outliers don't skew it.
    // Noise floor: a ≈ FWHM²_lowest - b * E_lowest, clamped ≥ 0.
    // Charge-trapping term c starts at 0 (usually tiny for good detectors).
    double aSeed, bSeed, cSeed;
    {
        std::vector<std::pair<double,double>> pts;
        for (size_t i = 0; i < xIn.size(); i++)
            pts.push_back({xIn[i], yIn[i]});
        std::sort(pts.begin(), pts.end());

        // yIn values are FWHM² — b ≈ FWHM²/E at high energy (dominant term)
        std::vector<double> bEsts;
        for (const auto& [E, fwhm2] : pts)
            if (E > 10.0) bEsts.push_back(fwhm2 / E);
        std::sort(bEsts.begin(), bEsts.end());
        bSeed = bEsts.empty() ? mFwhmB_->GetNumber()
                              : bEsts[bEsts.size() / 2];

        double E0   = pts.front().first;
        double fw0  = pts.front().second;   // fw0 is FWHM²
        aSeed = std::max(1e-4, fw0 - bSeed * E0);
        cSeed = 0.0;
    }

    // Fit on a temporary graph owned here (not drawn to canvas)
    TGraphErrors grFit((Int_t)xIn.size(), xIn.data(), yIn.data(),
                        exIn.data(), eyIn.data());

    // Read user-editable bounds (fall back to defaults if widgets not ready)
    double aLo = mFwhmAlo_ ? mFwhmAlo_->GetNumber() : 1e-4;
    double aHi = mFwhmAhi_ ? mFwhmAhi_->GetNumber() : 1000.0;
    double bLo = mFwhmBlo_ ? mFwhmBlo_->GetNumber() : 0.0;
    double bHi = mFwhmBhi_ ? mFwhmBhi_->GetNumber() : 10.0;
    double cLo = mFwhmClo_ ? mFwhmClo_->GetNumber() : 0.0;
    double cHi = mFwhmChi_ ? mFwhmChi_->GetNumber() : 0.1;
    if (aLo >= aHi) { AppendLog("Invalid a bounds (lo >= hi)."); return; }
    if (bLo >= bHi) { AppendLog("Invalid b bounds (lo >= hi)."); return; }

    delete fwhmTF1_;
    fwhmTF1_ = new TF1("fwhm_model", "[0]+[1]*x+[2]*x*x", 0.0, xhi);
    fwhmTF1_->SetParameter(0, std::max(aSeed, aLo));
    fwhmTF1_->SetParameter(1, std::max(bSeed, bLo));
    fwhmTF1_->SetParameter(2, cSeed);
    fwhmTF1_->SetParLimits(0, aLo, aHi);
    fwhmTF1_->SetParLimits(1, bLo, bHi);
    if (cLo < cHi)
        fwhmTF1_->SetParLimits(2, cLo, cHi);
    else
        fwhmTF1_->FixParameter(2, cLo);

    grFit.Fit(fwhmTF1_, "R Q B");

    double a = fwhmTF1_->GetParameter(0);
    double b = fwhmTF1_->GetParameter(1);
    double c = fwhmTF1_->GetParameter(2);
    mFwhmA_->SetNumber(a); mFwhmB_->SetNumber(b); mFwhmC_->SetNumber(c);

    double chi2    = fwhmTF1_->GetChisquare();
    int    ndf     = fwhmTF1_->GetNDF();
    fwhmNdf_      = ndf;
    fwhmChi2Ndf_  = (ndf > 0) ? chi2 / ndf : -1.0;
    fwhmPValue_   = (ndf > 0) ? TMath::Prob(chi2, ndf) : -1.0;

    // Residual RMS: RMS of (data - model) / data  (both in FWHM² space)
    double rmsSum = 0.0; int rmsN = 0;
    for (size_t i = 0; i < xIn.size(); i++) {
        double pred = fwhmTF1_->Eval(xIn[i]);   // FWHM² predicted
        if (pred > 0 && yIn[i] > 0) {
            double pull = (yIn[i] - pred) / yIn[i];
            rmsSum += pull * pull;
            ++rmsN;
        }
    }
    fwhmResidRMS_ = (rmsN > 0) ? std::sqrt(rmsSum / rmsN) : -1.0;

    std::string result = Form("a=%.4g  b=%.4g  c=%.4g  chi2/ndf=%.2f  p=%.3f",
                               a, b, c, fwhmChi2Ndf_, fwhmPValue_);
    fwhmResultLbl_->SetText(result.c_str());
    AppendLog("FWHM fit (" + std::to_string(xIn.size()) + " pts): " + result);
    SetStatus("FWHM fit done  chi2/ndf=" + Fmt(fwhmChi2Ndf_, 2));

    RedrawFWHM();
}

void GammaFitGUI::OnAcceptFWHM()
{
    if (!fwhmTF1_) {
        AppendLog("Fit the resolution model first.");
        return;
    }
    if (fwhmLoadedHists_.empty() && fwhmHistName_.empty()) {
        AppendLog("No FWHM histogram loaded.");
        return;
    }

    double a = fwhmTF1_->GetParameter(0);
    double b = fwhmTF1_->GetParameter(1);
    double c = fwhmTF1_->GetParameter(2);

    double chi2ndf = std::numeric_limits<double>::max();
    int ndf = fwhmTF1_->GetNDF();
    if (ndf > 0) chi2ndf = fwhmTF1_->GetChisquare() / ndf;

    // Save the resolution model + per-histogram exclusions to every loaded cache.
    const std::vector<std::string>& targets =
        fwhmLoadedHists_.empty() ? std::vector<std::string>{fwhmHistName_}
                                 : fwhmLoadedHists_;

    for (const std::string& hname : targets) {
        std::string cacheFile = CacheFileFor(hname);
        FitDatabase fitdb;
        fitdb.Load(cacheFile);

        FitEntry eRes;
        eRes.key = kResolutionKey; eRes.params = { a, b, c };
        eRes.chi2ndf = chi2ndf; eRes.residualRMS = 0.0; eRes.maxPull = 0.0;
        fitdb.ForceStore(kResolutionKey, eRes);

        // Build exclusion/inclusion lists for this histogram's points only
        FitEntry eExcl, eIncl;
        eExcl.key = kExcludedFwhmKey; eExcl.chi2ndf = 0.0; eExcl.residualRMS = 0.0; eExcl.maxPull = 0.0;
        eIncl.key = kFwhmInclKey;     eIncl.chi2ndf = 0.0; eIncl.residualRMS = 0.0; eIncl.maxPull = 0.0;
        for (size_t i = 0; i < fwhmAllX_.size(); i++) {
            const std::string& src = fwhmHistSources_.empty() ? hname : fwhmHistSources_[i];
            if (src != hname) continue;
            if (fwhmExcluded_[i])
                eExcl.params.push_back(fwhmAllX_[i]);
            else if (!fwhmTied_.empty() && fwhmTied_[i])
                eIncl.params.push_back(fwhmAllX_[i]);
        }
        fitdb.ForceStore(kExcludedFwhmKey, eExcl);
        fitdb.ForceStore(kFwhmInclKey, eIncl);

        fitdb.rootFile = inputPath_;
        EnsureCacheDir();
        fitdb.Save(cacheFile);
        BackupCacheFile(cacheFile);

        // Register in Fit Results list
        std::string entryName = std::string(kFwhmPrefix) + hname;
        bool found = false;
        for (const auto& n : fittedHists_) if (n == entryName) { found = true; break; }
        if (!found) {
            fittedHists_.push_back(entryName);
            fitResultsList_->AddEntry(entryName.c_str(), (Int_t)fittedHists_.size());
            fitResultsList_->MapSubwindows(); fitResultsList_->Layout();
        }
    }

    // Update live resolution model
    res_.a = a; res_.b = b; res_.c = c;

    int nExcl = (int)std::count(fwhmExcluded_.begin(), fwhmExcluded_.end(), true);
    std::string saved = (targets.size() == 1) ? targets[0]
                      : std::to_string(targets.size()) + " histograms";
    AppendLog("Resolution model saved to " + saved + ":  a=" + Fmt(a, 5) +
              "  b=" + Fmt(b, 5) + "  c=" + Fmt(c, 8) +
              "  excluded=" + std::to_string(nExcl) + " pts");
    SetStatus("FWHM model saved: " + saved);
}

