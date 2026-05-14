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
    TGGroupFrame* hg = new TGGroupFrame(p, "Histogram");
    p->AddFrame(hg, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

    {
        TGHorizontalFrame* comboRow = new TGHorizontalFrame(hg);
        hg->AddFrame(comboRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        fwhmCombo_ = new TGComboBox(comboRow, 800);
        fwhmCombo_->Resize(215, 22);
        comboRow->AddFrame(fwhmCombo_,
                           new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 0, 4, 0, 0));

        TGTextButton* loadBtn = new TGTextButton(comboRow, " Load from Cache ");
        comboRow->AddFrame(loadBtn, new TGLayoutHints(kLHintsCenterY, 0, 0, 0, 0));
        loadBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLoadFWHM()");
        loadBtn->SetToolTipText(
            "Load FWHM data points, excluded-point list, and resolution model\n"
            "from the selected histogram's fit cache file.");
    }

    TGLabel* formulaLbl = new TGLabel(hg, "  FWHM = sqrt(a + b*E + c*E^2)");
    hg->AddFrame(formulaLbl, new TGLayoutHints(kLHintsLeft, 4, 4, 0, 4));

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
    fitBtn->SetToolTipText("Fit FWHM = sqrt(a + b*E + c*E^2) to the included (non-excluded) data points");

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

    fwhmShowSigmaChk_ = new TGCheckButton(dg, "Show as sigma (FWHM / 2.3548)");
    dg->AddFrame(fwhmShowSigmaChk_, new TGLayoutHints(kLHintsLeft, 2, 2, 1, 4));
    fwhmShowSigmaChk_->Connect("Clicked()", "GammaFitGUI", this, "OnFWHMRedisplay()");
    fwhmShowSigmaChk_->SetToolTipText(
        "Switch Y axis between FWHM and Gaussian sigma.\n"
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

// Draw the FWHM plot onto any canvas.  All objects are heap-allocated and the
// canvas takes ownership via Draw() — do NOT delete them afterwards.
// fwhmTF1_ is NEVER drawn directly; fresh TF1 expressions are created here so
// that fwhmTF1_ stays fully under our control and safe to access later.
void GammaFitGUI::DrawFWHMToCanvas(TCanvas* c, bool showSigma, bool showStatLine,
                                    bool showResolution)
{
    if (fwhmAllX_.empty()) return;

    const double kSig = 2.3548200450309493;

    // Resolution (%) mode overrides sigma mode
    std::string yLabel, title;
    if (showResolution) {
        yLabel = "Resolution (%)";
        title  = "Resolution vs Energy  [" + fwhmHistName_ +
                 "];Energy (keV);Resolution (%)";
    } else if (showSigma) {
        yLabel = "#sigma (keV)";
        title  = "FWHM vs Energy  [" + fwhmHistName_ +
                 "];Energy (keV);#sigma (keV)";
    } else {
        yLabel = "FWHM (keV)";
        title  = "FWHM vs Energy  [" + fwhmHistName_ +
                 "];Energy (keV);FWHM (keV)";
    }

    double xhi = *std::max_element(fwhmAllX_.begin(), fwhmAllX_.end()) * 1.15;

    // ── Split into included / excluded ───────────────────────────────────────
    auto toY = [&](double fwhm, double E) -> double {
        if (showResolution) return 100.0 * fwhm / E;
        if (showSigma)      return fwhm / kSig;
        return fwhm;
    };

    std::vector<double> xIn, yIn, exIn, eyIn;
    std::vector<double> xEx, yEx;
    for (size_t i = 0; i < fwhmAllX_.size(); i++) {
        double y = toY(fwhmAllY_[i], fwhmAllX_[i]);
        if (fwhmExcluded_[i]) {
            xEx.push_back(fwhmAllX_[i]);
            yEx.push_back(y);
        } else {
            xIn.push_back(fwhmAllX_[i]);
            yIn.push_back(y);
            exIn.push_back(0.0);
            eyIn.push_back(0.05 * y);
        }
    }

    // ── Compute data y-range for axis (clamp model curve to this scale) ──────
    double yDataMax = 0.0;
    for (double v : yIn) yDataMax = std::max(yDataMax, v);
    for (double v : yEx) yDataMax = std::max(yDataMax, v);
    if (yDataMax <= 0) yDataMax = 1.0;
    double yAxisMax = yDataMax * 1.35;

    // ── Included points (blue filled) — canvas takes ownership ───────────────
    TGraphErrors* grIn;
    if (!xIn.empty()) {
        grIn = new TGraphErrors((Int_t)xIn.size(), xIn.data(), yIn.data(),
                                 exIn.data(), eyIn.data());
    } else {
        std::vector<double> xA, yA, exA, eyA;
        for (size_t i = 0; i < fwhmAllX_.size(); i++) {
            xA.push_back(fwhmAllX_[i]);
            yA.push_back(toY(fwhmAllY_[i], fwhmAllX_[i]));
            exA.push_back(0.0); eyA.push_back(0.0);
        }
        grIn = new TGraphErrors((Int_t)xA.size(), xA.data(), yA.data(),
                                 exA.data(), eyA.data());
        grIn->SetMarkerSize(0);
    }
    grIn->SetTitle(title.c_str());
    grIn->SetMarkerStyle(20); grIn->SetMarkerSize(0.9);
    grIn->SetMarkerColor(kBlue + 1); grIn->SetLineColor(kBlue + 1);
    grIn->SetMinimum(0.0);
    grIn->SetMaximum(yAxisMax);
    grIn->Draw("AP");
    if (grIn->GetXaxis()) {
        grIn->GetXaxis()->SetTitle("Energy (keV)");
        grIn->GetXaxis()->SetTitleSize(0.05);
        grIn->GetXaxis()->SetLabelSize(0.04);
    }
    if (grIn->GetYaxis()) {
        grIn->GetYaxis()->SetTitle(yLabel.c_str());
        grIn->GetYaxis()->SetTitleSize(0.05);
        grIn->GetYaxis()->SetLabelSize(0.04);
        grIn->GetYaxis()->SetTitleOffset(1.2);
    }

    // ── Excluded points (gray hollow) ────────────────────────────────────────
    if (!xEx.empty()) {
        TGraph* grEx = new TGraph((Int_t)xEx.size(), xEx.data(), yEx.data());
        grEx->SetMarkerStyle(24); grEx->SetMarkerSize(0.9);
        grEx->SetMarkerColor(kGray + 1);
        grEx->Draw("P same");
    }

    // ── Fano statistical limit (independent of model fit) ────────────────────
    // FWHM_stat = 2.344 * sqrt(F*w*E), F=0.12, w=0.00296 keV
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
                Form("%.10g*sqrt(%.10g*x)", kFanoCoeff, kFw), xlo, xhi);
        }
        statLine->SetLineColor(kGreen + 2);
        statLine->SetLineStyle(2);
        statLine->SetLineWidth(2);
        statLine->Draw("same");
    }

    // ── Model curve ───────────────────────────────────────────────────────────
    if (fwhmTF1_) {
        double a  = fwhmTF1_->GetParameter(0);
        double b  = fwhmTF1_->GetParameter(1);
        double cv = fwhmTF1_->GetParameter(2);

        TF1* modelDraw = nullptr;
        if (showResolution) {
            modelDraw = new TF1("fwhm_draw",
                Form("100.0*sqrt(%.10g+%.10g*x+%.10g*x*x)/x", a, b, cv),
                1.0, xhi);
        } else if (showSigma) {
            modelDraw = new TF1("fwhm_draw",
                Form("sqrt((%.10g+%.10g*x+%.10g*x*x)/%.10g)",
                     a, b, cv, kSig * kSig), 0.0, xhi);
        } else {
            modelDraw = new TF1("fwhm_draw",
                Form("sqrt(%.10g+%.10g*x+%.10g*x*x)", a, b, cv), 0.0, xhi);
        }
        modelDraw->SetLineColor(kRed);
        modelDraw->SetLineWidth(2);
        modelDraw->Draw("same");

        // ── Fit info pave (upper right) ───────────────────────────────────────
        TPaveText* pt = new TPaveText(0.52, 0.60, 0.93, 0.93, "NDC");
        pt->SetFillColor(0);
        pt->SetFillStyle(1001);
        pt->SetBorderSize(1);
        pt->SetTextSize(0.030);
        pt->SetTextAlign(12);

        // Equation line — mode-dependent label
        if (showResolution) {
            pt->AddText("R(%) = 100 #cdot #sqrt{a + b#cdot E + c#cdot E^{2}} / E");
        } else if (showSigma) {
            pt->AddText("#sigma = #sqrt{a + b#cdot E + c#cdot E^{2}} / 2.355");
        } else {
            pt->AddText("FWHM = #sqrt{a + b#cdot E + c#cdot E^{2}}");
        }

        pt->AddText(Form("a = %.5g  keV^{2}", a));
        pt->AddText(Form("b = %.5g  keV",     b));
        pt->AddText(Form("c = %.5g",           cv));

        if (fwhmChi2Ndf_ >= 0) {
            pt->AddText(Form("#chi^{2}/ndf = %.2f  (ndf = %d)",
                             fwhmChi2Ndf_, fwhmNdf_));
        }
        if (fwhmPValue_ >= 0) {
            pt->AddText(Form("p-value = %.4f", fwhmPValue_));
        }
        if (fwhmResidRMS_ >= 0) {
            pt->AddText(Form("Residual RMS = %.3f", fwhmResidRMS_));
        }

        pt->Draw();
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
    DrawFWHMToCanvas(c, showSigma, showStatLine, showRes);
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
    fwhmHistName_ = histNames_[id - 1];

    FitDatabase fitdb;
    if (!fitdb.Load(CacheFileFor(fwhmHistName_))) {
        AppendLog("No cache found for " + fwhmHistName_ + " — run AutoFit first.");
        return;
    }

    // Collect all (energy, FWHM) pairs from NGauss entries
    fwhmAllX_.clear(); fwhmAllY_.clear(); fwhmExcluded_.clear();
    for (const auto& [key, entry] : fitdb.GetEntries()) {
        if (key.size() >= 2 && key[0] == '_' && key[1] == '_') continue;
        int npar = (int)entry.params.size();
        if (npar < 5 || (npar - 2) % 3 != 0) continue;
        int n = (npar - 2) / 3;
        for (int i = 0; i < n; i++) {
            double E   = entry.params[3*i + 1];
            double sig = entry.params[3*i + 2];
            if (sig <= 0.0 || E <= 0.0) continue;
            fwhmAllX_.push_back(E);
            fwhmAllY_.push_back(2.3548 * sig);
            fwhmExcluded_.push_back(false);
        }
    }

    if (fwhmAllX_.empty()) {
        AppendLog("No fitted peaks in cache for " + fwhmHistName_ +
                  " — run AutoFit or accept manual fits first.");
        return;
    }

    // Restore saved exclusions from cache (match by energy within 0.5 keV)
    auto exIt = fitdb.GetEntries().find(kExcludedFwhmKey);
    if (exIt != fitdb.GetEntries().end()) {
        for (const double exE : exIt->second.params) {
            for (size_t i = 0; i < fwhmAllX_.size(); i++) {
                if (std::abs(fwhmAllX_[i] - exE) < 0.5) {
                    fwhmExcluded_[i] = true; break;
                }
            }
        }
    }

    // Pre-populate parameter boxes from cached __RESOLUTION__
    delete fwhmTF1_; fwhmTF1_ = nullptr;
    fwhmChi2Ndf_ = -1.0; fwhmPValue_ = -1.0; fwhmResidRMS_ = -1.0; fwhmNdf_ = 0;
    auto rit = fitdb.GetEntries().find(kResolutionKey);
    if (rit != fitdb.GetEntries().end() && rit->second.params.size() == 3) {
        mFwhmA_->SetNumber(rit->second.params[0]);
        mFwhmB_->SetNumber(rit->second.params[1]);
        mFwhmC_->SetNumber(rit->second.params[2]);
        double xhi = *std::max_element(fwhmAllX_.begin(), fwhmAllX_.end()) * 1.15;
        fwhmTF1_ = new TF1("fwhm_model", "sqrt([0]+[1]*x+[2]*x*x)", 0.0, xhi);
        fwhmTF1_->SetParameters(rit->second.params[0],
                                 rit->second.params[1],
                                 rit->second.params[2]);
        // Restore fit quality from cache if available
        double cachedChi2 = rit->second.chi2ndf;
        if (cachedChi2 > 0 && cachedChi2 < 1e10) {
            fwhmChi2Ndf_ = cachedChi2;
            // p-value requires ndf; not stored separately, so skip pvalue here
        }
    }

    int nExcl = (int)std::count(fwhmExcluded_.begin(), fwhmExcluded_.end(), true);
    AppendLog("FWHM data loaded: " + std::to_string(fwhmAllX_.size()) +
              " points (" + std::to_string(nExcl) + " excluded) from " + fwhmHistName_);
    SetStatus("FWHM: " + fwhmHistName_);
    RedrawFWHM();
}

void GammaFitGUI::OnFitFWHM()
{
    if (fwhmAllX_.empty()) {
        AppendLog("Load FWHM data first.");
        return;
    }

    // Build a graph from non-excluded points only
    std::vector<double> xIn, yIn, exIn, eyIn;
    for (size_t i = 0; i < fwhmAllX_.size(); i++) {
        if (fwhmExcluded_[i]) continue;
        double y = fwhmAllY_[i];
        xIn.push_back(fwhmAllX_[i]); yIn.push_back(y);
        exIn.push_back(0.0); eyIn.push_back(0.05 * y);
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

        std::vector<double> bEsts;
        for (const auto& [E, fwhm] : pts)
            if (E > 10.0) bEsts.push_back(fwhm * fwhm / E);
        std::sort(bEsts.begin(), bEsts.end());
        bSeed = bEsts.empty() ? mFwhmB_->GetNumber()
                              : bEsts[bEsts.size() / 2];

        double E0   = pts.front().first;
        double fw0  = pts.front().second;
        aSeed = std::max(1e-4, fw0 * fw0 - bSeed * E0);
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
    fwhmTF1_ = new TF1("fwhm_model", "sqrt([0]+[1]*x+[2]*x*x)", 0.0, xhi);
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

    // Residual RMS: RMS of (data - model) / data
    double rmsSum = 0.0; int rmsN = 0;
    for (size_t i = 0; i < xIn.size(); i++) {
        double pred = fwhmTF1_->Eval(xIn[i]);
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
    if (fwhmHistName_.empty()) {
        AppendLog("No FWHM histogram loaded.");
        return;
    }

    double a = fwhmTF1_->GetParameter(0);
    double b = fwhmTF1_->GetParameter(1);
    double c = fwhmTF1_->GetParameter(2);

    std::string cacheFile = CacheFileFor(fwhmHistName_);
    FitDatabase fitdb;
    fitdb.Load(cacheFile);

    // Save resolution model
    FitEntry eRes;
    eRes.key    = kResolutionKey;
    eRes.params = { a, b, c };
    int ndf     = fwhmTF1_->GetNDF();
    eRes.chi2ndf = (ndf > 0) ? fwhmTF1_->GetChisquare() / ndf
                              : std::numeric_limits<double>::max();
    eRes.residualRMS = 0.0; eRes.maxPull = 0.0;
    fitdb.ForceStore(kResolutionKey, eRes);

    // Save excluded point list so it persists across sessions
    FitEntry eExcl;
    eExcl.key = kExcludedFwhmKey;
    eExcl.chi2ndf = 0.0; eExcl.residualRMS = 0.0; eExcl.maxPull = 0.0;
    for (size_t i = 0; i < fwhmAllX_.size(); i++)
        if (fwhmExcluded_[i]) eExcl.params.push_back(fwhmAllX_[i]);
    fitdb.ForceStore(kExcludedFwhmKey, eExcl);

    fitdb.rootFile = inputPath_;
    mkdir(kCacheDir, 0755);
    fitdb.Save(cacheFile);

    // Update live resolution model
    res_.a = a; res_.b = b; res_.c = c;

    // Register in Fit Results list
    std::string entryName = std::string(kFwhmPrefix) + fwhmHistName_;
    bool found = false;
    for (const auto& n : fittedHists_) if (n == entryName) { found = true; break; }
    if (!found) {
        fittedHists_.push_back(entryName);
        fitResultsList_->AddEntry(entryName.c_str(), (Int_t)fittedHists_.size());
        fitResultsList_->MapSubwindows(); fitResultsList_->Layout();
    }

    int nExcl = (int)std::count(fwhmExcluded_.begin(), fwhmExcluded_.end(), true);
    AppendLog("Resolution model saved: a=" + Fmt(a, 5) +
              "  b=" + Fmt(b, 5) + "  c=" + Fmt(c, 8) +
              "  excluded=" + std::to_string(nExcl) + " pts");
    SetStatus("FWHM model saved: " + fwhmHistName_);
}

