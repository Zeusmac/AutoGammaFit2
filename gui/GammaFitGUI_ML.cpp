#include "GammaFitGUI.h"
#include "PeakFitter.h"

#include "TGLabel.h"
#include "TGButton.h"
#include "TGListBox.h"
#include "TGComboBox.h"
#include "TGNumberEntry.h"
#include "TGTextEntry.h"
#include "TGTextView.h"
#include "TGLayout.h"
#include "TGFrame.h"
#include "TGTab.h"
#include "TGCanvas.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TH1F.h"
#include "TRandom3.h"
#include "TMath.h"

#include <algorithm>
#include <cstdio>
#include <ctime>
#include <fstream>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <vector>

// ════════════════════════════════════════════════════════════════════════════
// File-level helpers
// ════════════════════════════════════════════════════════════════════════════

static std::string MLDir(const std::string& launchDir)
{
    return launchDir + "/ml";
}
static std::string MLApprovedPath(const std::string& launchDir)
{
    return launchDir + "/ml/approved_spectra.txt";
}

static std::vector<std::string> ReadApprovedLines(const std::string& path)
{
    std::vector<std::string> out;
    std::ifstream f(path);
    if (!f) return out;
    std::string line;
    while (std::getline(f, line))
        if (!line.empty() && line[0] != '#') out.push_back(line);
    return out;
}

static bool ApprovedContains(const std::vector<std::string>& lines,
                             const std::string& histName)
{
    for (const auto& l : lines) {
        if (l.size() > histName.size() &&
            l.substr(0, histName.size()) == histName &&
            (l[histName.size()] == ' ' || l[histName.size()] == '\t'))
            return true;
    }
    return false;
}

// Returns "exists  size  date" for a file, or "not found".
static std::string FileStatus(const std::string& path)
{
    struct stat st{};
    if (::stat(path.c_str(), &st) != 0) return "not found";
    long long kb = st.st_size / 1024;
    std::time_t t = st.st_mtime;
    char buf[32];
    std::strftime(buf, sizeof(buf), "%Y-%m-%d %H:%M", std::localtime(&t));
    return std::string(buf) + "   " + std::to_string(kb) + " kB";
}

// Load a text file into a string (returns empty string if not found).
static std::string ReadTextFile(const std::string& path)
{
    std::ifstream f(path);
    if (!f) return "";
    return std::string(std::istreambuf_iterator<char>(f),
                       std::istreambuf_iterator<char>());
}

// ── convenience: add a section header label ──────────────────────────────────
static TGLabel* SectionLabel(TGCompositeFrame* parent, const char* text)
{
    auto* lbl = new TGLabel(parent, text);
    lbl->SetTextFont("helvetica-bold-r-*-*-24-*-*-*-*-*-*-*");
    parent->AddFrame(lbl, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, 6, 6, 14, 4));
    return lbl;
}

// ── convenience: add a mono-spaced line ──────────────────────────────────────
static void MonoLine(TGCompositeFrame* parent, const char* text, int padL = 16)
{
    auto* lbl = new TGLabel(parent, text);
    lbl->SetTextFont("courier-medium-r-*-*-20-*-*-*-*-*-*-*");
    parent->AddFrame(lbl, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, padL, 6, 3, 3));
}

// ── convenience: add a description label ─────────────────────────────────────
static void DescLabel(TGCompositeFrame* parent, const char* text, int padL = 8)
{
    auto* lbl = new TGLabel(parent, text);
    lbl->SetTextFont("helvetica-medium-r-*-*-20-*-*-*-*-*-*-*");
    parent->AddFrame(lbl, new TGLayoutHints(kLHintsLeft | kLHintsExpandX, padL, 6, 4, 4));
}

// ════════════════════════════════════════════════════════════════════════════
// BuildMLTab
// ════════════════════════════════════════════════════════════════════════════

void GammaFitGUI::BuildMLTab(TGCompositeFrame* tab)
{
    TGTab* mlTabs = new TGTab(tab, 308, 600);
    tab->AddFrame(mlTabs, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    // ── Sub-tab 1: Architecture ───────────────────────────────────────────────
    {
        TGCompositeFrame* p = mlTabs->AddTab("Architecture");

        TGCanvas* sc = new TGCanvas(p, 10, 600, kSunkenFrame);
        p->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        TGVerticalFrame* vf = new TGVerticalFrame(sc->GetViewPort(), 10, 10);
        sc->SetContainer(vf);

        // ── Background MLP ────────────────────────────────────────────────────
        TGGroupFrame* bgGrp = new TGGroupFrame(vf, "Background Estimation MLP");
        bgGrp->SetTitlePos(TGGroupFrame::kLeft);
        bgGrp->SetTextFont("helvetica-bold-r-*-*-22-*-*-*-*-*-*-*");
        vf->AddFrame(bgGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 8, 8));

        SectionLabel(bgGrp, "Architecture");
        MonoLine(bgGrp, "Input(101) -> FC(64) -> ReLU");
        MonoLine(bgGrp, "          -> FC(32) -> ReLU -> FC(1)");

        SectionLabel(bgGrp, "Sliding window");
        DescLabel(bgGrp, "101 bins centred on target bin b");
        DescLabel(bgGrp, "Zero-padded at spectrum edges");
        DescLabel(bgGrp, "Step: 1 bin  (one inference per bin)");

        SectionLabel(bgGrp, "Normalisation");
        MonoLine(bgGrp, "raw_max = max( counts[0..N] )");
        MonoLine(bgGrp, "win[k]  = counts[b-50+k] / raw_max");

        SectionLabel(bgGrp, "Output & scaling");
        MonoLine(bgGrp, "bg_norm = MLP( win )");
        MonoLine(bgGrp, "bkg[b]  = max( 0, bg_norm * raw_max )");
        DescLabel(bgGrp, "Subtracted from working histogram.");
        DescLabel(bgGrp, "Negative bins floored to 0.");

        SectionLabel(bgGrp, "Training loss");
        MonoLine(bgGrp, "MSELoss( bg_norm, true_bg_norm )");
        DescLabel(bgGrp, "Trained on 20 000 synthetic spectra.");

        // ── Peak Detector MLP ─────────────────────────────────────────────────
        TGGroupFrame* pkGrp = new TGGroupFrame(vf, "Peak Detector MLP");
        pkGrp->SetTitlePos(TGGroupFrame::kLeft);
        pkGrp->SetTextFont("helvetica-bold-r-*-*-22-*-*-*-*-*-*-*");
        vf->AddFrame(pkGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 8, 8));

        SectionLabel(pkGrp, "Architecture");
        MonoLine(pkGrp, "Input(51) -> FC(32) -> ReLU");
        MonoLine(pkGrp, "         -> FC(16) -> ReLU -> FC(1) -> Sigmoid");

        SectionLabel(pkGrp, "Sliding window");
        DescLabel(pkGrp, "51 bins centred on target bin b");
        DescLabel(pkGrp, "Applied to background-subtracted spectrum.");

        SectionLabel(pkGrp, "Output");
        MonoLine(pkGrp, "prob[b] = Peak_MLP( sub_norm_window )");
        DescLabel(pkGrp, "prob[b] in [0, 1]: probability that bin b");
        DescLabel(pkGrp, "contains a gamma peak.");

        SectionLabel(pkGrp, "Training loss (class imbalance fix)");
        MonoLine(pkGrp, "BCELoss  pos_weight = 20");
        DescLabel(pkGrp, "~1-5 % of bins are true peaks.");
        DescLabel(pkGrp, "pos_weight = 20 prevents the model");
        DescLabel(pkGrp, "from predicting 'no peak everywhere'.");

        // ── Non-Maximum Suppression ───────────────────────────────────────────
        TGGroupFrame* nmsGrp = new TGGroupFrame(vf, "Non-Maximum Suppression");
        nmsGrp->SetTitlePos(TGGroupFrame::kLeft);
        nmsGrp->SetTextFont("helvetica-bold-r-*-*-22-*-*-*-*-*-*-*");
        vf->AddFrame(nmsGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 8, 8));

        DescLabel(nmsGrp, "1.  Keep bins where prob[b] >= threshold");
        DescLabel(nmsGrp, "2.  In each region: keep the bin with");
        DescLabel(nmsGrp, "    the highest probability (local max");
        DescLabel(nmsGrp, "    within the ±25-bin window).");
        DescLabel(nmsGrp, "3.  Output bin centres -> AdaptiveFitter");
        DescLabel(nmsGrp, "    (same interface as TSpectrum output).");

        // ── Synthetic Training Data ───────────────────────────────────────────
        TGGroupFrame* synGrp = new TGGroupFrame(vf, "Synthetic Training Data");
        synGrp->SetTitlePos(TGGroupFrame::kLeft);
        synGrp->SetTextFont("helvetica-bold-r-*-*-22-*-*-*-*-*-*-*");
        vf->AddFrame(synGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 8, 8));

        SectionLabel(synGrp, "Energy axis");
        MonoLine(synGrp, "E in [100, 3000] keV  (4096 bins)");

        SectionLabel(synGrp, "Background shapes (random per spectrum)");
        MonoLine(synGrp, "Flat:        bg = c0");
        MonoLine(synGrp, "             c0 in [10, 2000]");
        MonoLine(synGrp, "Linear:      bg = c0 + c1*(E - Emin)");
        MonoLine(synGrp, "Exponential: bg = A * exp(-E / tau)");
        MonoLine(synGrp, "             tau in [300, 1500] keV");

        SectionLabel(synGrp, "Peaks");
        MonoLine(synGrp, "N ~ Poisson(5)");
        MonoLine(synGrp, "E_i ~ Uniform(150, 2800) keV");
        MonoLine(synGrp, "FWHM^2 = a + b*E_i  (±30 % per spectrum)");
        MonoLine(synGrp, "A_i ~ Uniform(2, 500) * bg(E_i)");

        SectionLabel(synGrp, "Noise");
        MonoLine(synGrp, "counts ~ Poisson( bg + sum_peaks )");

        SectionLabel(synGrp, "Background truth label");
        MonoLine(synGrp, "y_bg[b] = bg[b] / raw_max");

        SectionLabel(synGrp, "Peak truth label");
        MonoLine(synGrp, "y_peak[b] = 1 if bin b holds a peak");
        MonoLine(synGrp, "            centre, else 0");

        // ── Real Spectrum BG Truth ────────────────────────────────────────────
        TGGroupFrame* realGrp = new TGGroupFrame(vf, "Real-Spectrum BG Truth");
        realGrp->SetTitlePos(TGGroupFrame::kLeft);
        realGrp->SetTextFont("helvetica-bold-r-*-*-22-*-*-*-*-*-*-*");
        vf->AddFrame(realGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 8, 10));

        DescLabel(realGrp, "For each approved spectrum in");
        DescLabel(realGrp, "ml/approved_spectra.txt:");
        MonoLine(realGrp, "bg[b] = spectrum[b]");
        MonoLine(realGrp, "      - sum_i A_i*Gauss(E[b]; E_i, sig_i)");
        DescLabel(realGrp, "using fitted peaks from fit_caches/*.dat.");
        DescLabel(realGrp, "");
        DescLabel(realGrp, "Bins within 3*sigma of a known peak are");
        DescLabel(realGrp, "excluded from BG training windows");
        DescLabel(realGrp, "(BG truth unreliable there).");
        DescLabel(realGrp, "");
        DescLabel(realGrp, "Auto-skip rules:");
        MonoLine(realGrp, "> 20 % entries have needsRefit=1");
        MonoLine(realGrp, "any chi2/ndf > 5");

    } // end Architecture sub-tab

    // ── Sub-tab 2: Training Data ──────────────────────────────────────────────
    {
        TGCompositeFrame* p = mlTabs->AddTab("Training Data");

        TGCanvas* sc = new TGCanvas(p, 10, 600, kSunkenFrame);
        p->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        TGVerticalFrame* vf = new TGVerticalFrame(sc->GetViewPort(), 10, 10);
        sc->SetContainer(vf);

        // ── Approved Spectra ──────────────────────────────────────────────────
        TGGroupFrame* apGrp = new TGGroupFrame(vf, "Approved Spectra");
        vf->AddFrame(apGrp, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 2));

        DescLabel(apGrp, "Histograms whose fit cache is trusted");
        DescLabel(apGrp, "as ground truth for ML training.");

        mlApprovedList_ = new TGListBox(apGrp, -1, kSunkenFrame | kDoubleBorder);
        mlApprovedList_->SetMultipleSelections(kTRUE);
        mlApprovedList_->Resize(280, 130);
        apGrp->AddFrame(mlApprovedList_,
                        new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        {
            TGHorizontalFrame* r = new TGHorizontalFrame(apGrp);
            apGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            auto* addBtn  = new TGTextButton(r, " Add current ");
            auto* remBtn  = new TGTextButton(r, " Remove sel. ");
            auto* clrBtn  = new TGTextButton(r, " Clear all ");
            auto* refBtn  = new TGTextButton(r, " Refresh ");
            r->AddFrame(addBtn,  new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
            r->AddFrame(remBtn,  new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
            r->AddFrame(clrBtn,  new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
            r->AddFrame(refBtn,  new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
            addBtn->Connect("Clicked()", "GammaFitGUI", this, "OnMLTabAddCurrent()");
            remBtn->Connect("Clicked()", "GammaFitGUI", this, "OnMLTabRemoveApproved()");
            clrBtn->Connect("Clicked()", "GammaFitGUI", this, "OnMLTabClearApproved()");
            refBtn->Connect("Clicked()", "GammaFitGUI", this, "OnMLTabRefreshApproved()");
            addBtn->SetToolTipText("Append the currently selected histogram\n"
                                   "to ml/approved_spectra.txt");
            remBtn->SetToolTipText("Remove selected list entry from the file");
            clrBtn->SetToolTipText("Remove ALL entries from the file");
            refBtn->SetToolTipText("Re-read ml/approved_spectra.txt");
        }

        mlApprovedCntLbl_ = new TGLabel(apGrp, "0 spectra approved");
        apGrp->AddFrame(mlApprovedCntLbl_,
                        new TGLayoutHints(kLHintsLeft, 4, 4, 2, 4));

        // ── Dataset Files ─────────────────────────────────────────────────────
        TGGroupFrame* dsGrp = new TGGroupFrame(vf, "Dataset Files");
        vf->AddFrame(dsGrp, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));

        {
            auto* r = new TGHorizontalFrame(dsGrp);
            dsGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 1));
            r->AddFrame(new TGLabel(r, "synthetic.npz:"),
                        new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
            mlSynDataLbl_ = new TGLabel(r, "not found");
            r->AddFrame(mlSynDataLbl_, new TGLayoutHints(kLHintsCenterY));
        }
        {
            auto* r = new TGHorizontalFrame(dsGrp);
            dsGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 2));
            r->AddFrame(new TGLabel(r, "combined.npz:"),
                        new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
            mlCombDataLbl_ = new TGLabel(r, "not found");
            r->AddFrame(mlCombDataLbl_, new TGLayoutHints(kLHintsCenterY));
        }

        // ── Generation parameters ─────────────────────────────────────────
        {
            auto* r = new TGHorizontalFrame(dsGrp);
            dsGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 6, 1));
            r->AddFrame(new TGLabel(r, "Spectra:"),
                        new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            mlNSpectraEntry_ = new TGNumberEntry(r, 5000, 7, -1,
                TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive,
                TGNumberFormat::kNELLimitMinMax, 100, 200000);
            r->AddFrame(mlNSpectraEntry_,
                        new TGLayoutHints(kLHintsCenterY, 0, 16, 0, 0));
            r->AddFrame(new TGLabel(r, "Windows/spectrum:"),
                        new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            mlWinPerSpecEntry_ = new TGNumberEntry(r, 100, 5, -1,
                TGNumberFormat::kNESInteger, TGNumberFormat::kNEAPositive,
                TGNumberFormat::kNELLimitMinMax, 10, 2000);
            r->AddFrame(mlWinPerSpecEntry_,
                        new TGLayoutHints(kLHintsCenterY, 0, 0, 0, 0));
            mlNSpectraEntry_->Connect("ValueSet(Long_t)", "GammaFitGUI", this,
                                      "OnMLUpdateRAMEstimate()");
            mlNSpectraEntry_->GetNumberEntry()->Connect("ReturnPressed()", "GammaFitGUI", this,
                                                        "OnMLUpdateRAMEstimate()");
            mlWinPerSpecEntry_->Connect("ValueSet(Long_t)", "GammaFitGUI", this,
                                        "OnMLUpdateRAMEstimate()");
            mlWinPerSpecEntry_->GetNumberEntry()->Connect("ReturnPressed()", "GammaFitGUI", this,
                                                          "OnMLUpdateRAMEstimate()");
        }
        {
            auto* r = new TGHorizontalFrame(dsGrp);
            dsGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 4));
            mlRAMEstLbl_ = new TGLabel(r, "Est. RAM: —");
            r->AddFrame(mlRAMEstLbl_, new TGLayoutHints(kLHintsCenterY));
        }
        OnMLUpdateRAMEstimate();  // populate label with defaults on startup
        {
            auto* r = new TGHorizontalFrame(dsGrp);
            dsGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 1));
            auto* genBtn  = new TGTextButton(r, " Generate synthetic ");
            auto* loadBtn = new TGTextButton(r, " Load real spectra ");
            auto* refBtn  = new TGTextButton(r, " Refresh ");
            r->AddFrame(genBtn,  new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
            r->AddFrame(loadBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
            r->AddFrame(refBtn,  new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
            genBtn->Connect("Clicked()", "GammaFitGUI", this, "OnMLGenSynthetic()");
            loadBtn->Connect("Clicked()", "GammaFitGUI", this, "OnMLLoadReal()");
            refBtn->Connect("Clicked()", "GammaFitGUI", this, "OnMLRefreshDataStatus()");
            genBtn->SetToolTipText(
                "Launch:  python3 train/generate_spectra.py\n"
                "Uses the Spectra and Windows/spectrum values above.\n"
                "Output:  data/synthetic.npz\n"
                "Runs in background; stdout → ml/gen.log");
            loadBtn->SetToolTipText(
                "Launch:  python3 train/load_real_spectra.py\n"
                "Extracts windows from approved ROOT spectra.\n"
                "Output:  data/combined.npz\n"
                "Requires data/synthetic.npz to exist first.\n"
                "Runs in background; stdout → ml/load_real.log");
        }
        DescLabel(dsGrp,
            "Tip: only approve spectra where every visible\n"
            "feature is in fit_caches/ with chi2/ndf < 5\n"
            "and at most 1-2 needsRefit flags.");

    } // end Training Data sub-tab

    // ── Sub-tab 3: Train ──────────────────────────────────────────────────────
    {
        TGCompositeFrame* p = mlTabs->AddTab("Train");

        TGCanvas* sc = new TGCanvas(p, 10, 600, kSunkenFrame);
        p->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        TGVerticalFrame* vf = new TGVerticalFrame(sc->GetViewPort(), 10, 10);
        sc->SetContainer(vf);

        // ── Model Status ──────────────────────────────────────────────────────
        TGGroupFrame* modGrp = new TGGroupFrame(vf, "Model Files");
        vf->AddFrame(modGrp, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 2));

        {
            auto* r = new TGHorizontalFrame(modGrp);
            modGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 1));
            r->AddFrame(new TGLabel(r, "bg_model.onnx:"),
                        new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
            mlBgModelLbl_ = new TGLabel(r, "not found");
            r->AddFrame(mlBgModelLbl_, new TGLayoutHints(kLHintsCenterY));
        }
        {
            auto* r = new TGHorizontalFrame(modGrp);
            modGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 2));
            r->AddFrame(new TGLabel(r, "peak_model.onnx:"),
                        new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
            mlPeakModelLbl_ = new TGLabel(r, "not found");
            r->AddFrame(mlPeakModelLbl_, new TGLayoutHints(kLHintsCenterY));
        }
        {
            auto* r = new TGHorizontalFrame(modGrp);
            modGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 4));
            auto* refBtn = new TGTextButton(r, " Refresh status ");
            auto* relBtn = new TGTextButton(r, " Reload models ");
            r->AddFrame(refBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            r->AddFrame(relBtn, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
            refBtn->Connect("Clicked()", "GammaFitGUI", this, "OnMLRefreshModelStatus()");
            relBtn->Connect("Clicked()", "GammaFitGUI", this, "OnMLReloadModels()");
            refBtn->SetToolTipText("Refresh the model file status labels above");
            relBtn->SetToolTipText(
                "Clear cached model instances so the next AutoFit\n"
                "re-reads ml/*.onnx from disk.");
        }

        // ── Launch Training ───────────────────────────────────────────────────
        TGGroupFrame* trainGrp = new TGGroupFrame(vf, "Launch Training");
        vf->AddFrame(trainGrp, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        DescLabel(trainGrp, "Runs in background.  stdout -> ml/*.log");

        {
            auto* r = new TGHorizontalFrame(trainGrp);
            trainGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 2));
            auto* bgBtn    = new TGTextButton(r, " Train BG model ");
            auto* peakBtn  = new TGTextButton(r, " Train Peak model ");
            auto* sigmaBtn = new TGTextButton(r, " Train σ model ");
            r->AddFrame(bgBtn,    new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            r->AddFrame(peakBtn,  new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            r->AddFrame(sigmaBtn, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
            bgBtn->Connect("Clicked()",    "GammaFitGUI", this, "OnMLTrainBg()");
            peakBtn->Connect("Clicked()",  "GammaFitGUI", this, "OnMLTrainPeak()");
            sigmaBtn->Connect("Clicked()", "GammaFitGUI", this, "OnMLTrainSigma()");
            bgBtn->SetToolTipText(
                "Launch:  python3 train/train_bg.py\n"
                "Trains Linear(101→64→32→1) with MSELoss.\n"
                "Output:  ml/bg_model.onnx\n"
                "Log:     ml/train_bg.log");
            peakBtn->SetToolTipText(
                "Launch:  python3 train/train_peak.py\n"
                "Trains Linear(51→32→16→1→Sigmoid) with BCELoss.\n"
                "Output:  ml/peak_model.onnx\n"
                "Log:     ml/train_peak.log");
            sigmaBtn->SetToolTipText(
                "Launch:  python3 train/train_sigma.py\n"
                "Trains Linear(51→32→16→1→Softplus) with Huber loss.\n"
                "Output:  ml/sigma_model.onnx\n"
                "Log:     ml/train_sigma.log\n"
                "Constrains MIGRAD sigma to ±8% of prediction,\n"
                "reducing area uncertainty by breaking A-sigma correlation.");
        }

        // ── Training Log ──────────────────────────────────────────────────────
        TGGroupFrame* logGrp = new TGGroupFrame(vf, "Training Log");
        vf->AddFrame(logGrp, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));

        {
            auto* r = new TGHorizontalFrame(logGrp);
            logGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            r->AddFrame(new TGLabel(r, "Show log:"),
                        new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
            mlLogCombo_ = new TGComboBox(r, -1, kSunkenFrame | kDoubleBorder);
            mlLogCombo_->AddEntry("Background model  (train_bg.log)",  1);
            mlLogCombo_->AddEntry("Peak model  (train_peak.log)",       2);
            mlLogCombo_->AddEntry("Generate synthetic  (gen.log)",      3);
            mlLogCombo_->AddEntry("Load real spectra  (load_real.log)", 4);
            mlLogCombo_->AddEntry("Sigma model  (train_sigma.log)",     5);
            mlLogCombo_->Resize(180, 20);
            mlLogCombo_->Select(1, kFALSE);
            r->AddFrame(mlLogCombo_, new TGLayoutHints(kLHintsCenterY));
            mlLogCombo_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                  "OnMLLogComboChanged(Int_t)");

            auto* refBtn = new TGTextButton(r, " Refresh ");
            r->AddFrame(refBtn, new TGLayoutHints(kLHintsCenterY, 4, 0, 0, 0));
            refBtn->Connect("Clicked()", "GammaFitGUI", this, "OnMLRefreshLog()");
        }

        mlTrainLogView_ = new TGTextView(logGrp, 280, 260, kSunkenFrame | kDoubleBorder);
        logGrp->AddFrame(mlTrainLogView_,
                         new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));
        mlTrainLogView_->LoadBuffer("(no log yet — click Refresh after training starts)");

    } // end Train sub-tab

    // ── Sub-tab 4: Inspect ────────────────────────────────────────────────────
    {
        TGCompositeFrame* p = mlTabs->AddTab("Inspect");

        TGCanvas* sc = new TGCanvas(p, 10, 600, kSunkenFrame);
        p->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        TGVerticalFrame* vf = new TGVerticalFrame(sc->GetViewPort(), 10, 10);
        sc->SetContainer(vf);

        // ── Step-by-step inference ────────────────────────────────────────────
        TGGroupFrame* infGrp = new TGGroupFrame(vf, "Inference on Current Histogram");
        vf->AddFrame(infGrp, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 2));

        DescLabel(infGrp,
            "Run the ML pipeline step-by-step on the histogram\n"
            "currently loaded in the main canvas.\n"
            "Requires USE_ONNX=1 build + trained model files.");

        {
            TGHorizontalFrame* r = new TGHorizontalFrame(infGrp);
            infGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 6, 2));
            auto* b1 = new TGTextButton(r, " Raw + ML BG ");
            auto* b2 = new TGTextButton(r, " BG-subtracted ");
            r->AddFrame(b1, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            r->AddFrame(b2, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
            b1->Connect("Clicked()", "GammaFitGUI", this, "OnMLInspectBgOverlay()");
            b2->Connect("Clicked()", "GammaFitGUI", this, "OnMLInspectBgSub()");
            b1->SetToolTipText(
                "Draw the raw histogram on the main canvas with the ML\n"
                "background estimate overlaid in orange.");
            b2->SetToolTipText(
                "Show the ML background-subtracted histogram on the main canvas.");
        }
        {
            TGHorizontalFrame* r = new TGHorizontalFrame(infGrp);
            infGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            auto* b3 = new TGTextButton(r, " Peak probabilities ");
            auto* b4 = new TGTextButton(r, " All 4 steps ");
            r->AddFrame(b3, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            r->AddFrame(b4, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
            b3->Connect("Clicked()", "GammaFitGUI", this, "OnMLInspectPeakProb()");
            b4->Connect("Clicked()", "GammaFitGUI", this, "OnMLInspectAllSteps()");
            b3->SetToolTipText(
                "Open a new canvas showing the peak-detector MLP probability\n"
                "per bin on the BG-subtracted spectrum. The dashed line marks\n"
                "the current peak threshold.");
            b4->SetToolTipText(
                "Open a 4-pad canvas: raw+BG, BG-subtracted, peak\n"
                "probabilities, and BG-subtracted with NMS peak markers.");
        }
        {
            TGHorizontalFrame* r = new TGHorizontalFrame(infGrp);
            infGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 6));
            auto* b5 = new TGTextButton(r, " Compare TSpectrum vs ML BG ");
            r->AddFrame(b5, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
            b5->Connect("Clicked()", "GammaFitGUI", this, "OnMLInspectCompare()");
            b5->SetToolTipText(
                "Open a new canvas with both the TSpectrum background estimate\n"
                "and the ML BG estimate overlaid on the raw spectrum.");
        }

        // ── Training curves ───────────────────────────────────────────────────
        TGGroupFrame* curvGrp = new TGGroupFrame(vf, "Training Curves");
        vf->AddFrame(curvGrp, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        DescLabel(curvGrp,
            "Parse the training logs and plot loss / metrics vs epoch.\n"
            "Click a training button in the Train tab first, then refresh.");

        {
            TGHorizontalFrame* r = new TGHorizontalFrame(curvGrp);
            curvGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 6, 6));
            auto* b1 = new TGTextButton(r, " BG model curves ");
            auto* b2 = new TGTextButton(r, " Peak model curves ");
            r->AddFrame(b1, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            r->AddFrame(b2, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
            b1->Connect("Clicked()", "GammaFitGUI", this, "OnMLPlotBgCurves()");
            b2->Connect("Clicked()", "GammaFitGUI", this, "OnMLPlotPeakCurves()");
            b1->SetToolTipText(
                "Parse ml/train_bg.log and plot train MSE and val MSE\n"
                "vs epoch in a new ROOT canvas.");
            b2->SetToolTipText(
                "Parse ml/train_peak.log and plot train loss, val loss,\n"
                "precision, recall, and F1 vs epoch in a new ROOT canvas.");
        }

        // ── Synthetic preview ─────────────────────────────────────────────────
        TGGroupFrame* synGrp = new TGGroupFrame(vf, "Synthetic Spectrum Preview");
        vf->AddFrame(synGrp, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));

        DescLabel(synGrp,
            "Generate a random synthetic spectrum in C++ using the same\n"
            "algorithm as the Python training generator.\n"
            "Shows raw counts, the true background, and the true peak\n"
            "positions so you can see what the models are trained on.");

        {
            TGHorizontalFrame* r = new TGHorizontalFrame(synGrp);
            synGrp->AddFrame(r, new TGLayoutHints(kLHintsExpandX, 2, 2, 6, 6));
            auto* b1 = new TGTextButton(r, " New random spectrum ");
            r->AddFrame(b1, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
            b1->Connect("Clicked()", "GammaFitGUI", this, "OnMLPreviewSynSpectrum()");
            b1->SetToolTipText(
                "Generate a fresh random synthetic training spectrum and display\n"
                "it in a new ROOT canvas. The true background is shown in orange\n"
                "and each peak position is marked with a vertical dashed line.");
        }

    } // end Inspect sub-tab

    // Populate the approved list and status labels on startup
    OnMLTabRefreshApproved();
    OnMLRefreshDataStatus();
    OnMLRefreshModelStatus();
}

// ════════════════════════════════════════════════════════════════════════════
// Training Data sub-tab slots
// ════════════════════════════════════════════════════════════════════════════

void GammaFitGUI::OnMLTabRefreshApproved()
{
    if (!mlApprovedList_) return;
    mlApprovedList_->RemoveAll();

    const std::string path = MLApprovedPath(launchDir_);
    auto lines = ReadApprovedLines(path);

    int id = 1;
    for (const auto& l : lines) {
        // Show histogram name only; full line stored as tooltip on label
        std::string histName = l.substr(0, l.find_first_of(" \t"));
        mlApprovedList_->AddEntry(histName.c_str(), id++);
    }
    mlApprovedList_->MapSubwindows();
    mlApprovedList_->Layout();

    int n = (int)lines.size();
    if (mlApprovedCntLbl_) {
        std::string s = std::to_string(n) + " spectra approved";
        mlApprovedCntLbl_->SetText(s.c_str());
    }

    // Also sync the AutoFit tab checkbox/count if loaded
    RefreshMLApprovalState();
}

void GammaFitGUI::OnMLTabAddCurrent()
{
    if (currentHist_.empty()) {
        AppendLog("ML training: no histogram selected.");
        return;
    }
    const std::string filePath = inputPath_.empty() ? srcRootPath_ : inputPath_;
    if (filePath.empty()) {
        AppendLog("ML training: no ROOT file open.");
        return;
    }

    const std::string path = MLApprovedPath(launchDir_);
    auto lines = ReadApprovedLines(path);
    if (ApprovedContains(lines, currentHist_)) {
        AppendLog("ML training: " + currentHist_ + " already approved.");
        return;
    }

    gSystem->MakeDirectory(MLDir(launchDir_).c_str());
    std::ofstream f(path, std::ios::app);
    if (!f) { AppendLog("ML training: cannot write " + path); return; }
    f << currentHist_ << "\t" << filePath << "\n";
    AppendLog("ML training: approved  " + currentHist_);
    OnMLTabRefreshApproved();
}

void GammaFitGUI::OnMLTabRemoveApproved()
{
    if (!mlApprovedList_) return;

    // Collect selected IDs
    TList selected;
    mlApprovedList_->GetSelectedEntries(&selected);
    if (selected.IsEmpty()) {
        AppendLog("ML training: no entry selected.");
        return;
    }

    // Build set of histogram names to remove (1-based IDs into ReadApprovedLines vector)
    const std::string path = MLApprovedPath(launchDir_);
    auto lines = ReadApprovedLines(path);

    std::vector<bool> remove(lines.size(), false);
    TIter next(&selected);
    while (TGLBEntry* entry = (TGLBEntry*)next()) {
        int idx = entry->EntryId() - 1;
        if (idx >= 0 && idx < (int)lines.size())
            remove[idx] = true;
    }

    // Rewrite file preserving header comments
    std::ifstream fin(path);
    std::vector<std::string> comments;
    if (fin) {
        std::string line;
        while (std::getline(fin, line))
            if (!line.empty() && line[0] == '#') comments.push_back(line);
    }
    fin.close();

    std::ofstream fout(path, std::ios::trunc);
    for (const auto& c : comments) fout << c << "\n";
    for (int i = 0; i < (int)lines.size(); ++i)
        if (!remove[i]) fout << lines[i] << "\n";

    AppendLog("ML training: removed selected entries.");
    OnMLTabRefreshApproved();
}

void GammaFitGUI::OnMLTabClearApproved()
{
    const std::string path = MLApprovedPath(launchDir_);
    std::ifstream fin(path);
    std::vector<std::string> comments;
    if (fin) {
        std::string line;
        while (std::getline(fin, line))
            if (!line.empty() && line[0] == '#') comments.push_back(line);
    }
    fin.close();
    std::ofstream fout(path, std::ios::trunc);
    for (const auto& c : comments) fout << c << "\n";
    AppendLog("ML training: cleared all approved spectra.");
    OnMLTabRefreshApproved();
}

void GammaFitGUI::OnMLUpdateRAMEstimate()
{
    if (!mlNSpectraEntry_ || !mlWinPerSpecEntry_ || !mlRAMEstLbl_) return;

    long nSpec = (long)mlNSpectraEntry_->GetNumber();
    long winPS = (long)mlWinPerSpecEntry_->GetNumber();

    // bg_X:   nSpec × winPS × 101 floats × 4 bytes
    // peak_X: nSpec × winPS ×  51 floats × 4 bytes
    // + 20 % Python list / numpy conversion overhead
    double mb = (double)nSpec * winPS * (101 + 51) * 4 * 1.2 / (1024.0 * 1024.0);

    long availMB = -1;
    {
        std::ifstream f("/proc/meminfo");
        std::string key, unit;
        long val = 0;
        while (f >> key >> val >> unit)
            if (key == "MemAvailable:") { availMB = val / 1024; break; }
    }

    char buf[200];
    if (availMB >= 0) {
        const char* warn = (mb > availMB * 0.8) ? "  ** exceeds available RAM! **"
                         : (mb > availMB * 0.4) ? "  (large — consider fewer spectra)"
                         : "";
        snprintf(buf, sizeof(buf), "Est. RAM: ~%.0f MB  |  Available: %ld MB%s",
                 mb, availMB, warn);
    } else {
        snprintf(buf, sizeof(buf), "Est. RAM: ~%.0f MB", mb);
    }
    mlRAMEstLbl_->SetText(buf);
    mlRAMEstLbl_->Resize(mlRAMEstLbl_->GetDefaultSize());
}

void GammaFitGUI::OnMLRefreshDataStatus()
{
    const std::string base = launchDir_ + "/data";
    auto updateLbl = [&](TGLabel* lbl, const std::string& relPath) {
        if (!lbl) return;
        std::string s = FileStatus(launchDir_ + "/" + relPath);
        lbl->SetText(s.c_str());
        lbl->Resize(lbl->GetDefaultSize());
    };
    updateLbl(mlSynDataLbl_,  "data/synthetic.npz");
    updateLbl(mlCombDataLbl_, "data/combined.npz");
}

// ── script launchers (background via shell & operator) ────────────────────────
static void LaunchScript(const std::string& launchDir,
                         const std::string& script,
                         const std::string& logFile)
{
    // Ensure the ml/ directory exists for log output
    gSystem->MakeDirectory((launchDir + "/ml").c_str());
    std::string cmd =
        "bash -c 'cd " + launchDir +
        " && python3 " + script +
        " > " + launchDir + "/" + logFile + " 2>&1' &";
    gSystem->Exec(cmd.c_str());
}

void GammaFitGUI::OnMLGenSynthetic()
{
    int nSpec  = mlNSpectraEntry_  ? (int)mlNSpectraEntry_->GetNumber()  : 5000;
    int winPS  = mlWinPerSpecEntry_ ? (int)mlWinPerSpecEntry_->GetNumber() : 100;

    char scriptWithArgs[256];
    snprintf(scriptWithArgs, sizeof(scriptWithArgs),
             "train/generate_spectra.py --n-spectra %d --windows-per-spectrum %d",
             nSpec, winPS);

    AppendLog(std::string("Launching generate_spectra.py (") +
              std::to_string(nSpec) + " spectra, " +
              std::to_string(winPS) + " win/spec) → ml/gen.log");
    LaunchScript(launchDir_, scriptWithArgs, "ml/gen.log");
    SetStatus("generate_spectra.py launched — check log viewer for progress");
    if (mlLogCombo_) mlLogCombo_->Select(3);   // 3 = gen.log
    OnMLRefreshLog();
}

void GammaFitGUI::OnMLLoadReal()
{
    AppendLog("Launching load_real_spectra.py in background → ml/load_real.log");
    LaunchScript(launchDir_, "train/load_real_spectra.py", "ml/load_real.log");
    SetStatus("load_real_spectra.py launched (background)");
}

// ════════════════════════════════════════════════════════════════════════════
// Train sub-tab slots
// ════════════════════════════════════════════════════════════════════════════

void GammaFitGUI::OnMLTrainBg()
{
    AppendLog("Launching train_bg.py in background → ml/train_bg.log");
    LaunchScript(launchDir_, "train/train_bg.py", "ml/train_bg.log");
    SetStatus("train_bg.py launched (background)");
}

void GammaFitGUI::OnMLTrainPeak()
{
    AppendLog("Launching train_peak.py in background → ml/train_peak.log");
    LaunchScript(launchDir_, "train/train_peak.py", "ml/train_peak.log");
    SetStatus("train_peak.py launched (background)");
}

void GammaFitGUI::OnMLTrainSigma()
{
    AppendLog("Launching train_sigma.py in background → ml/train_sigma.log");
    LaunchScript(launchDir_, "train/train_sigma.py", "ml/train_sigma.log");
    SetStatus("train_sigma.py launched (background)");
    if (mlLogCombo_) mlLogCombo_->Select(5);   // 5 = train_sigma.log
    OnMLRefreshLog();
}

void GammaFitGUI::OnMLRefreshModelStatus()
{
    const std::string mldir = launchDir_ + "/ml";
    auto updateLbl = [&](TGLabel* lbl, const std::string& fname) {
        if (!lbl) return;
        std::string s = FileStatus(mldir + "/" + fname);
        lbl->SetText(s.c_str());
        lbl->Resize(lbl->GetDefaultSize());
    };
    updateLbl(mlBgModelLbl_,   "bg_model.onnx");
    updateLbl(mlPeakModelLbl_, "peak_model.onnx");
}

void GammaFitGUI::OnMLRefreshLog()
{
    if (!mlTrainLogView_ || !mlLogCombo_) return;
    int sel = mlLogCombo_->GetSelected();
    std::string logFile;
    switch (sel) {
        case 1: logFile = launchDir_ + "/ml/train_bg.log";    break;
        case 2: logFile = launchDir_ + "/ml/train_peak.log";  break;
        case 3: logFile = launchDir_ + "/ml/gen.log";         break;
        case 4: logFile = launchDir_ + "/ml/load_real.log";   break;
        case 5: logFile = launchDir_ + "/ml/train_sigma.log"; break;
        default: return;
    }
    std::string content = ReadTextFile(logFile);
    mlTrainLogView_->Clear();
    if (content.empty()) {
        mlTrainLogView_->LoadBuffer("(log file not found or empty)");
    } else {
        // Split into lines and load each
        std::istringstream ss(content);
        std::string line;
        std::vector<std::string> lines;
        while (std::getline(ss, line)) lines.push_back(line);
        // TGTextView::LoadBuffer for multi-line: use AddLine
        mlTrainLogView_->Clear();
        for (const auto& l : lines) mlTrainLogView_->AddLine(l.c_str());
        // Scroll to bottom to show latest output
        mlTrainLogView_->ShowBottom();
    }
}

void GammaFitGUI::OnMLLogComboChanged(Int_t /*id*/)
{
    OnMLRefreshLog();
}

// ════════════════════════════════════════════════════════════════════════════
// AutoFit tab approval slot (invoked from the checkbox in AutoFit tab)
// ════════════════════════════════════════════════════════════════════════════

void GammaFitGUI::RefreshMLApprovalState()
{
    if (!mlTrainingApprovedChk_ || !mlTrainingCountLbl_) return;

    auto lines = ReadApprovedLines(MLApprovedPath(launchDir_));

    std::string cnt = std::to_string((int)lines.size()) + " spectra approved";
    mlTrainingCountLbl_->SetText(cnt.c_str());
    mlTrainingCountLbl_->Resize(mlTrainingCountLbl_->GetDefaultSize());

    if (currentHist_.empty()) {
        mlTrainingApprovedChk_->SetState(kButtonUp);
        return;
    }
    bool approved = ApprovedContains(lines, currentHist_);
    mlTrainingApprovedChk_->SetState(approved ? kButtonDown : kButtonUp);
}

void GammaFitGUI::OnMLTrainingApprovalToggled()
{
    if (currentHist_.empty()) {
        AppendLog("ML training: no histogram selected.");
        if (mlTrainingApprovedChk_) mlTrainingApprovedChk_->SetState(kButtonUp);
        return;
    }
    const std::string filePath = inputPath_.empty() ? srcRootPath_ : inputPath_;
    if (filePath.empty()) {
        AppendLog("ML training: no ROOT file open.");
        if (mlTrainingApprovedChk_) mlTrainingApprovedChk_->SetState(kButtonUp);
        return;
    }

    const std::string path = MLApprovedPath(launchDir_);
    bool nowChecked = mlTrainingApprovedChk_ && mlTrainingApprovedChk_->IsOn();
    auto lines      = ReadApprovedLines(path);

    if (nowChecked) {
        if (!ApprovedContains(lines, currentHist_)) {
            gSystem->MakeDirectory(MLDir(launchDir_).c_str());
            std::ofstream f(path, std::ios::app);
            if (!f) { AppendLog("ML training: cannot write " + path); return; }
            f << currentHist_ << "\t" << filePath << "\n";
            AppendLog("ML training: approved  " + currentHist_);
        }
    } else {
        bool removed = false;
        std::vector<std::string> kept;
        for (const auto& l : lines) {
            bool match = (l.size() > currentHist_.size() &&
                          l.substr(0, currentHist_.size()) == currentHist_ &&
                          (l[currentHist_.size()] == ' ' ||
                           l[currentHist_.size()] == '\t'));
            if (match) { removed = true; continue; }
            kept.push_back(l);
        }
        if (removed) {
            std::ifstream fin(path);
            std::vector<std::string> comments;
            if (fin) {
                std::string line;
                while (std::getline(fin, line))
                    if (!line.empty() && line[0] == '#') comments.push_back(line);
            }
            fin.close();
            std::ofstream fout(path, std::ios::trunc);
            for (const auto& c : comments) fout << c << "\n";
            for (const auto& l : kept)     fout << l << "\n";
            AppendLog("ML training: removed   " + currentHist_);
        }
    }

    RefreshMLApprovalState();
    // Keep ML tab list in sync
    OnMLTabRefreshApproved();
}

// ════════════════════════════════════════════════════════════════════════════
// Reload models
// ════════════════════════════════════════════════════════════════════════════

void GammaFitGUI::OnMLReloadModels()
{
    PeakFitter::ReloadMLModels();
    SetStatus("ML models cleared — will reload from disk on next AutoFit.");
    AppendLog("ML models cleared (will reload from ml/*.onnx on next run).");
    OnMLRefreshModelStatus();
}

// ════════════════════════════════════════════════════════════════════════════
// Inspect sub-tab — inference visualisation
// ════════════════════════════════════════════════════════════════════════════

// Helper: return current ML model paths from the AutoFit tab widgets
static std::string BgModelPath()  { return "ml/bg_model.onnx"; }
static std::string PkModelPath()  { return "ml/peak_model.onnx"; }

// Draw raw histogram with ML background curve overlaid on the main canvas.
void GammaFitGUI::OnMLInspectBgOverlay()
{
    if (!rawHist_) { AppendLog("Inspect: no histogram loaded."); return; }
    TH1* bg = PeakFitter::GetMLBackground(rawHist_, BgModelPath());
    if (!bg) {
        AppendLog("Inspect: ML BG model not available "
                  "(need USE_ONNX=1 build + ml/bg_model.onnx).");
        return;
    }
    // Draw raw on main canvas, then overlay BG curve
    SafeDeleteHist(viewHist_);
    DrawOnCanvas(rawHist_);
    TCanvas* c = canvas_->GetCanvas();
    c->cd();
    bg->SetLineColor(kOrange + 1);
    bg->SetLineStyle(2);
    bg->SetLineWidth(2);
    bg->SetStats(0);
    bg->SetBit(kCanDelete);
    bg->Draw("hist same");
    TLegend* leg = new TLegend(0.55, 0.82, 0.98, 0.94, "", "NDC");
    leg->SetBorderSize(1);
    leg->SetFillColor(kWhite);
    leg->AddEntry(rawHist_, "Raw spectrum", "l");
    leg->AddEntry(bg, "ML background estimate", "l");
    leg->SetBit(kCanDelete);
    leg->Draw();
    c->Modified(); c->Update();
    SetStatus("Inspect: raw + ML BG overlay");
}

// Replace the main canvas view with the ML BG-subtracted histogram.
void GammaFitGUI::OnMLInspectBgSub()
{
    if (!rawHist_) { AppendLog("Inspect: no histogram loaded."); return; }
    TH1* h = MakeMLBgSubHist(rawHist_);
    if (!h) {
        AppendLog("Inspect: ML BG model not available "
                  "(need USE_ONNX=1 build + ml/bg_model.onnx).");
        return;
    }
    SafeDeleteHist(viewHist_);
    viewHist_ = h;
    RedrawView();
    SetStatus("Inspect: ML BG-subtracted view");
}

// Open a standalone TCanvas showing the peak-detector probability per bin.
void GammaFitGUI::OnMLInspectPeakProb()
{
    if (!rawHist_) { AppendLog("Inspect: no histogram loaded."); return; }
    TH1* bgSub = MakeMLBgSubHist(rawHist_);
    if (!bgSub) {
        AppendLog("Inspect: ML BG model not available.");
        return;
    }
    TH1* prob = PeakFitter::GetMLPeakProbabilities(bgSub, PkModelPath());
    delete bgSub;
    if (!prob) {
        AppendLog("Inspect: ML peak model not available "
                  "(need USE_ONNX=1 build + ml/peak_model.onnx).");
        return;
    }
    TCanvas* c = new TCanvas("ml_peak_prob",
                             Form("Peak probabilities — %s", rawHist_->GetName()),
                             1000, 400);
    c->SetLeftMargin(0.08); c->SetRightMargin(0.02);
    prob->SetTitle(Form("Peak-detector MLP probability — %s", rawHist_->GetName()));
    prob->GetXaxis()->SetTitle("Energy (keV)");
    prob->GetYaxis()->SetTitle("Peak probability");
    prob->SetLineColor(kBlue + 1);
    prob->SetFillColor(kBlue - 9);
    prob->SetFillStyle(1001);
    prob->SetStats(0);
    prob->SetMinimum(0.0);
    prob->SetMaximum(1.05);
    prob->SetBit(kCanDelete);
    prob->Draw("hist");
    // Draw threshold line
    float thr = mlThreshEntry_ ? (float)mlThreshEntry_->GetNumber() : 0.5f;
    TLine* tline = new TLine(prob->GetXaxis()->GetXmin(), thr,
                             prob->GetXaxis()->GetXmax(), thr);
    tline->SetLineColor(kRed);
    tline->SetLineStyle(2);
    tline->SetLineWidth(2);
    tline->SetBit(kCanDelete);
    tline->Draw();
    TLatex* lthr = new TLatex();
    lthr->SetNDC(kFALSE);
    lthr->SetTextSize(0.04);
    lthr->SetTextColor(kRed);
    lthr->SetBit(kCanDelete);
    lthr->DrawLatex(prob->GetXaxis()->GetXmin() + 30, thr + 0.02,
                    Form("threshold = %.2f", thr));
    c->Modified(); c->Update();
    SetStatus("Inspect: peak probabilities plotted");
}

// 4-pad TCanvas: raw+BG, BG-sub, peak-prob, BG-sub with peak markers.
void GammaFitGUI::OnMLInspectAllSteps()
{
    if (!rawHist_) { AppendLog("Inspect: no histogram loaded."); return; }

    TH1* bg = PeakFitter::GetMLBackground(rawHist_, BgModelPath());
    if (!bg) {
        AppendLog("Inspect: ML BG model not available.");
        return;
    }
    TH1* bgSub = MakeMLBgSubHist(rawHist_);
    TH1* prob  = bgSub ? PeakFitter::GetMLPeakProbabilities(bgSub, PkModelPath()) : nullptr;

    TCanvas* c = new TCanvas("ml_all_steps",
                             Form("ML inference steps — %s", rawHist_->GetName()),
                             1400, 900);
    c->Divide(2, 2, 0.005, 0.005);

    // ── Pad 1: raw + BG curve ─────────────────────────────────────────────────
    c->cd(1);
    gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.02);
    TH1* raw1 = (TH1*)rawHist_->Clone("_ml_raw1");
    raw1->SetTitle("Step 1 — Raw + ML background");
    raw1->SetStats(0);
    raw1->SetLineColor(kBlack);
    raw1->SetBit(kCanDelete);
    raw1->Draw("hist");
    bg->SetLineColor(kOrange + 1);
    bg->SetLineWidth(2);
    bg->SetLineStyle(2);
    bg->SetStats(0);
    bg->SetBit(kCanDelete);
    bg->Draw("hist same");
    TLegend* l1 = new TLegend(0.55, 0.80, 0.98, 0.94);
    l1->SetBorderSize(1); l1->SetFillColor(kWhite);
    l1->AddEntry(raw1, "Raw", "l");
    l1->AddEntry(bg,   "ML BG", "l");
    l1->SetBit(kCanDelete);
    l1->Draw();

    // ── Pad 2: BG-subtracted ──────────────────────────────────────────────────
    c->cd(2);
    gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.02);
    if (bgSub) {
        bgSub->SetTitle("Step 2 — BG-subtracted");
        bgSub->SetStats(0);
        bgSub->SetLineColor(kAzure + 1);
        bgSub->SetBit(kCanDelete);
        bgSub->Draw("hist");
    } else {
        TLatex* t = new TLatex(0.3, 0.5, "BG subtraction failed");
        t->SetNDC(); t->SetBit(kCanDelete); t->Draw();
    }

    // ── Pad 3: peak probability ───────────────────────────────────────────────
    c->cd(3);
    gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.02);
    float thr = mlThreshEntry_ ? (float)mlThreshEntry_->GetNumber() : 0.5f;
    if (prob) {
        prob->SetTitle("Step 3 — Peak probabilities");
        prob->GetYaxis()->SetTitle("Probability");
        prob->SetLineColor(kBlue + 1);
        prob->SetFillColor(kBlue - 9);
        prob->SetFillStyle(1001);
        prob->SetStats(0);
        prob->SetMinimum(0.0); prob->SetMaximum(1.05);
        prob->SetBit(kCanDelete);
        prob->Draw("hist");
        TLine* tl = new TLine(prob->GetXaxis()->GetXmin(), thr,
                              prob->GetXaxis()->GetXmax(), thr);
        tl->SetLineColor(kRed); tl->SetLineStyle(2); tl->SetLineWidth(2);
        tl->SetBit(kCanDelete); tl->Draw();
    } else {
        TLatex* t = new TLatex(0.3, 0.5, "Peak model not available");
        t->SetNDC(); t->SetBit(kCanDelete); t->Draw();
    }

    // ── Pad 4: BG-sub with NMS peak markers ───────────────────────────────────
    c->cd(4);
    gPad->SetLeftMargin(0.10); gPad->SetRightMargin(0.02);
    if (bgSub && prob) {
        TH1* bs4 = (TH1*)bgSub->Clone("_ml_bgsub4");
        bs4->SetTitle("Step 4 — Detected peaks (NMS)");
        bs4->SetStats(0);
        bs4->SetLineColor(kAzure + 1);
        bs4->SetBit(kCanDelete);
        bs4->Draw("hist");
        double yMax = bs4->GetMaximum();
        const int nBins = prob->GetNbinsX();
        // NMS: find local maxima above threshold
        for (int b = 1; b <= nBins; ++b) {
            float p = (float)prob->GetBinContent(b);
            if (p < thr) continue;
            bool isMax = true;
            for (int k = -25; k <= 25 && isMax; ++k) {
                int nb = b + k;
                if (nb < 1 || nb > nBins || nb == b) continue;
                if (prob->GetBinContent(nb) > p) isMax = false;
            }
            if (!isMax) continue;
            double E = prob->GetBinCenter(b);
            TLine* pl = new TLine(E, 0, E, yMax * 0.95);
            pl->SetLineColor(kRed); pl->SetLineStyle(2); pl->SetLineWidth(1);
            pl->SetBit(kCanDelete); pl->Draw();
        }
    } else {
        TLatex* t = new TLatex(0.3, 0.5, "Models not available");
        t->SetNDC(); t->SetBit(kCanDelete); t->Draw();
    }

    c->Modified(); c->Update();
    SetStatus("Inspect: all-steps canvas opened");
}

// Overlay TSpectrum BG and ML BG on the raw histogram in a new canvas.
void GammaFitGUI::OnMLInspectCompare()
{
    if (!rawHist_) { AppendLog("Inspect: no histogram loaded."); return; }

    TH1* mlBg   = PeakFitter::GetMLBackground(rawHist_, BgModelPath());
    int  iters  = bgIterEntry_ ? (int)bgIterEntry_->GetNumber() : 14;
    TH1* tsBg   = GetTSpectrumBg(rawHist_, iters);

    TCanvas* c = new TCanvas("ml_bg_compare",
                             Form("BG comparison — %s", rawHist_->GetName()),
                             1200, 500);
    c->SetLeftMargin(0.08); c->SetRightMargin(0.02);

    TH1* raw1 = (TH1*)rawHist_->Clone("_ml_cmp_raw");
    raw1->SetTitle(Form("Background comparison — %s", rawHist_->GetName()));
    raw1->SetStats(0);
    raw1->SetLineColor(kBlack);
    raw1->SetBit(kCanDelete);
    raw1->Draw("hist");

    TLegend* leg = new TLegend(0.55, 0.78, 0.98, 0.94);
    leg->SetBorderSize(1); leg->SetFillColor(kWhite);
    leg->AddEntry(raw1, "Raw spectrum", "l");

    if (tsBg) {
        tsBg->SetLineColor(kGreen + 2);
        tsBg->SetLineStyle(2);
        tsBg->SetLineWidth(2);
        tsBg->SetStats(0);
        tsBg->SetBit(kCanDelete);
        tsBg->Draw("hist same");
        leg->AddEntry(tsBg, Form("TSpectrum BG (%d iter)", iters), "l");
    } else {
        AppendLog("Inspect compare: TSpectrum BG failed.");
    }
    if (mlBg) {
        mlBg->SetLineColor(kOrange + 1);
        mlBg->SetLineStyle(2);
        mlBg->SetLineWidth(2);
        mlBg->SetStats(0);
        mlBg->SetBit(kCanDelete);
        mlBg->Draw("hist same");
        leg->AddEntry(mlBg, "ML background estimate", "l");
    } else {
        AppendLog("Inspect compare: ML BG model not available.");
    }
    leg->SetBit(kCanDelete);
    leg->Draw();
    c->Modified(); c->Update();
    SetStatus("Inspect: BG comparison canvas opened");
}

// ════════════════════════════════════════════════════════════════════════════
// Inspect sub-tab — training curve plots
// ════════════════════════════════════════════════════════════════════════════

// Parse "Epoch  X/Y  train MSE=A  val MSE=B" lines from train_bg.log
void GammaFitGUI::OnMLPlotBgCurves()
{
    const std::string logPath = launchDir_ + "/ml/train_bg.log";
    std::ifstream f(logPath);
    if (!f) {
        AppendLog("Plot BG curves: " + logPath + " not found. "
                  "Train the BG model first.");
        return;
    }

    std::vector<double> epochs, trainLoss, valLoss;
    std::string line;
    while (std::getline(f, line)) {
        // Format: "Epoch   X/Y  train MSE=A  val MSE=B"
        int ep, total;
        double tr, vl;
        if (std::sscanf(line.c_str(), "Epoch %d/%d  train MSE=%lf  val MSE=%lf",
                        &ep, &total, &tr, &vl) == 4) {
            epochs.push_back((double)ep);
            trainLoss.push_back(tr);
            valLoss.push_back(vl);
        }
    }
    if (epochs.empty()) {
        AppendLog("Plot BG curves: no epoch lines found in " + logPath +
                  "\nExpected format: 'Epoch X/Y  train MSE=N  val MSE=N'");
        return;
    }

    TCanvas* c = new TCanvas("ml_bg_curves", "BG Model Training Curves", 900, 500);
    c->SetLeftMargin(0.10); c->SetRightMargin(0.04);
    c->SetGridy();

    TGraph* gTr = new TGraph((int)epochs.size(), epochs.data(), trainLoss.data());
    TGraph* gVl = new TGraph((int)epochs.size(), epochs.data(), valLoss.data());

    TMultiGraph* mg = new TMultiGraph("mg_bg", "BG Model — MSE Loss vs Epoch");
    gTr->SetLineColor(kBlue + 1);  gTr->SetLineWidth(2);
    gTr->SetMarkerColor(kBlue + 1); gTr->SetMarkerStyle(20); gTr->SetMarkerSize(0.8);
    gVl->SetLineColor(kRed + 1);   gVl->SetLineWidth(2);
    gVl->SetMarkerColor(kRed + 1); gVl->SetMarkerStyle(21); gVl->SetMarkerSize(0.8);
    mg->Add(gTr, "lp"); mg->Add(gVl, "lp");
    mg->SetBit(kCanDelete);
    mg->Draw("a");
    mg->GetXaxis()->SetTitle("Epoch");
    mg->GetYaxis()->SetTitle("MSE Loss (normalised)");
    mg->GetYaxis()->SetTitleOffset(1.3);

    TLegend* leg = new TLegend(0.55, 0.75, 0.93, 0.90);
    leg->SetBorderSize(1); leg->SetFillColor(kWhite);
    leg->AddEntry(gTr, "Train MSE", "lp");
    leg->AddEntry(gVl, "Val MSE",   "lp");
    leg->SetBit(kCanDelete);
    leg->Draw();

    c->Modified(); c->Update();
    SetStatus(Form("BG training curves: %d epochs plotted", (int)epochs.size()));
}

// Parse "Epoch X/Y  train=A  val=B  P=C  R=D  F1=E" from train_peak.log
void GammaFitGUI::OnMLPlotPeakCurves()
{
    const std::string logPath = launchDir_ + "/ml/train_peak.log";
    std::ifstream f(logPath);
    if (!f) {
        AppendLog("Plot Peak curves: " + logPath + " not found. "
                  "Train the peak model first.");
        return;
    }

    std::vector<double> epochs, trainLoss, valLoss, prec, rec, f1;
    std::string line;
    while (std::getline(f, line)) {
        int ep, total;
        double tr, vl, p, r, fi;
        if (std::sscanf(line.c_str(),
                "Epoch %d/%d  train=%lf  val=%lf  P=%lf  R=%lf  F1=%lf",
                &ep, &total, &tr, &vl, &p, &r, &fi) == 7) {
            epochs.push_back((double)ep);
            trainLoss.push_back(tr);
            valLoss.push_back(vl);
            prec.push_back(p);
            rec.push_back(r);
            f1.push_back(fi);
        }
    }
    if (epochs.empty()) {
        AppendLog("Plot Peak curves: no epoch lines found in " + logPath);
        return;
    }

    TCanvas* c = new TCanvas("ml_peak_curves", "Peak Model Training Curves",
                             1200, 500);
    c->Divide(2, 1, 0.005, 0.005);

    // Left: loss curves
    c->cd(1);
    gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.03); gPad->SetGridy();
    TGraph* gTr = new TGraph((int)epochs.size(), epochs.data(), trainLoss.data());
    TGraph* gVl = new TGraph((int)epochs.size(), epochs.data(), valLoss.data());
    TMultiGraph* mgL = new TMultiGraph("mgL", "Loss vs Epoch");
    gTr->SetLineColor(kBlue + 1); gTr->SetLineWidth(2);
    gTr->SetMarkerColor(kBlue+1); gTr->SetMarkerStyle(20); gTr->SetMarkerSize(0.7);
    gVl->SetLineColor(kRed + 1);  gVl->SetLineWidth(2);
    gVl->SetMarkerColor(kRed+1);  gVl->SetMarkerStyle(21); gVl->SetMarkerSize(0.7);
    mgL->Add(gTr, "lp"); mgL->Add(gVl, "lp");
    mgL->SetBit(kCanDelete); mgL->Draw("a");
    mgL->GetXaxis()->SetTitle("Epoch");
    mgL->GetYaxis()->SetTitle("BCE Loss"); mgL->GetYaxis()->SetTitleOffset(1.5);
    TLegend* l1 = new TLegend(0.55, 0.75, 0.95, 0.90);
    l1->SetBorderSize(1); l1->SetFillColor(kWhite);
    l1->AddEntry(gTr, "Train loss", "lp"); l1->AddEntry(gVl, "Val loss", "lp");
    l1->SetBit(kCanDelete); l1->Draw();

    // Right: P/R/F1
    c->cd(2);
    gPad->SetLeftMargin(0.12); gPad->SetRightMargin(0.03); gPad->SetGridy();
    TGraph* gP  = new TGraph((int)epochs.size(), epochs.data(), prec.data());
    TGraph* gR  = new TGraph((int)epochs.size(), epochs.data(), rec.data());
    TGraph* gF1 = new TGraph((int)epochs.size(), epochs.data(), f1.data());
    TMultiGraph* mgR = new TMultiGraph("mgR", "Precision / Recall / F1 vs Epoch");
    gP->SetLineColor(kGreen + 2); gP->SetLineWidth(2);
    gP->SetMarkerColor(kGreen+2); gP->SetMarkerStyle(22); gP->SetMarkerSize(0.7);
    gR->SetLineColor(kOrange + 1); gR->SetLineWidth(2);
    gR->SetMarkerColor(kOrange+1); gR->SetMarkerStyle(23); gR->SetMarkerSize(0.7);
    gF1->SetLineColor(kMagenta + 1); gF1->SetLineWidth(3);
    gF1->SetMarkerColor(kMagenta+1); gF1->SetMarkerStyle(29); gF1->SetMarkerSize(0.9);
    mgR->Add(gP, "lp"); mgR->Add(gR, "lp"); mgR->Add(gF1, "lp");
    mgR->GetYaxis()->SetRangeUser(0.0, 1.05);
    mgR->SetBit(kCanDelete); mgR->Draw("a");
    mgR->GetXaxis()->SetTitle("Epoch");
    mgR->GetYaxis()->SetTitle("Score"); mgR->GetYaxis()->SetTitleOffset(1.5);
    TLegend* l2 = new TLegend(0.15, 0.15, 0.55, 0.35);
    l2->SetBorderSize(1); l2->SetFillColor(kWhite);
    l2->AddEntry(gP,  "Precision", "lp");
    l2->AddEntry(gR,  "Recall",    "lp");
    l2->AddEntry(gF1, "F1 score",  "lp");
    l2->SetBit(kCanDelete); l2->Draw();

    c->Modified(); c->Update();
    SetStatus(Form("Peak training curves: %d epochs plotted", (int)epochs.size()));
}

// ════════════════════════════════════════════════════════════════════════════
// Inspect sub-tab — synthetic spectrum preview (pure C++ generator)
// ════════════════════════════════════════════════════════════════════════════

void GammaFitGUI::OnMLPreviewSynSpectrum()
{
    // Use ROOT's TRandom3 (Mersenne Twister) for reproducible-looking random spectra.
    // Seed 0 = machine entropy → different every call.
    TRandom3 rng(0);

    const int    nBins = 1024;        // display at 1024 bins (full 4096 is slow)
    const double Emin  = 100.0;
    const double Emax  = 3000.0;
    const double dE    = (Emax - Emin) / nBins;

    // ── Background ────────────────────────────────────────────────────────────
    int bgType = (int)rng.Uniform(0, 3);  // 0=flat, 1=linear, 2=exp
    std::vector<double> bg(nBins);
    if (bgType == 0) {
        double c0 = rng.Uniform(10, 1000);
        for (int b = 0; b < nBins; b++) bg[b] = c0;
    } else if (bgType == 1) {
        double c0 = rng.Uniform(50, 1000);
        double c1 = rng.Uniform(-0.3, 0.3);
        for (int b = 0; b < nBins; b++) {
            double E = Emin + (b + 0.5) * dE;
            bg[b] = std::max(5.0, c0 + c1 * (E - Emin));
        }
    } else {
        double A   = rng.Uniform(200, 2000);
        double tau = rng.Uniform(300, 1500);
        for (int b = 0; b < nBins; b++) {
            double E = Emin + (b + 0.5) * dE;
            bg[b] = A * TMath::Exp(-E / tau);
        }
    }

    // ── Peaks ─────────────────────────────────────────────────────────────────
    int nPeaks = (int)rng.Poisson(5);
    nPeaks = std::max(1, std::min(nPeaks, 12));
    struct Peak { double E, sigma, amp; };
    std::vector<Peak> peaks(nPeaks);
    // FWHM^2 = a + b*E, representative HPGe values with ±30% scatter
    double a_fwhm = rng.Uniform(0.7, 1.3) * 1.0;      // keV^2 constant term
    double b_fwhm = rng.Uniform(0.7, 1.3) * 0.0025;   // keV slope
    for (int i = 0; i < nPeaks; i++) {
        peaks[i].E     = rng.Uniform(150, 2800);
        double fwhm2   = a_fwhm + b_fwhm * peaks[i].E;
        double sigma   = TMath::Sqrt(fwhm2) / 2.3548;
        peaks[i].sigma = sigma;
        // Find bin for this energy to get bg value there
        int bIdx = (int)((peaks[i].E - Emin) / dE);
        bIdx = std::max(0, std::min(bIdx, nBins - 1));
        peaks[i].amp = rng.Uniform(2, 500) * bg[bIdx];
    }

    // ── Combine and add Poisson noise ─────────────────────────────────────────
    static int ctr = 0; ++ctr;
    TH1F* hObs = new TH1F(Form("ml_syn_%d", ctr),
                           "Synthetic training spectrum preview",
                           nBins, Emin, Emax);
    TH1F* hBg  = new TH1F(Form("ml_syn_bg_%d", ctr), "", nBins, Emin, Emax);
    hObs->SetDirectory(nullptr);
    hBg->SetDirectory(nullptr);

    for (int b = 0; b < nBins; b++) {
        double E    = Emin + (b + 0.5) * dE;
        double sig  = bg[b];
        for (const auto& pk : peaks)
            sig += pk.amp * TMath::Gaus(E, pk.E, pk.sigma, kFALSE);
        double obs = rng.Poisson(std::max(0.01, sig));
        hObs->SetBinContent(b + 1, obs);
        hBg->SetBinContent(b + 1, bg[b]);
    }

    // ── Draw ──────────────────────────────────────────────────────────────────
    static const char* bgNames[] = {"Flat", "Linear", "Exponential"};
    TCanvas* c = new TCanvas(Form("ml_syn_canvas_%d", ctr),
                             Form("Synthetic spectrum — %s BG, %d peaks",
                                  bgNames[bgType], nPeaks),
                             1200, 500);
    c->SetLeftMargin(0.08); c->SetRightMargin(0.02);

    hObs->SetLineColor(kBlack);
    hObs->SetStats(0);
    hObs->GetXaxis()->SetTitle("Energy (keV)");
    hObs->GetYaxis()->SetTitle("Counts / bin");
    hObs->SetBit(kCanDelete);
    hObs->Draw("hist");

    hBg->SetLineColor(kOrange + 1);
    hBg->SetLineWidth(2);
    hBg->SetLineStyle(2);
    hBg->SetStats(0);
    hBg->SetBit(kCanDelete);
    hBg->Draw("hist same");

    double yMax = hObs->GetMaximum();
    for (const auto& pk : peaks) {
        TLine* l = new TLine(pk.E, 0, pk.E, yMax * 0.92);
        l->SetLineColor(kRed); l->SetLineStyle(3); l->SetLineWidth(1);
        l->SetBit(kCanDelete); l->Draw();
    }

    TLegend* leg = new TLegend(0.55, 0.80, 0.97, 0.95);
    leg->SetBorderSize(1); leg->SetFillColor(kWhite);
    leg->AddEntry(hObs, "Observed counts (Poisson noise)", "l");
    leg->AddEntry(hBg,  Form("True background (%s)", bgNames[bgType]), "l");
    TLine* dummy = new TLine(); dummy->SetLineColor(kRed); dummy->SetLineStyle(3);
    dummy->SetBit(kCanDelete);
    leg->AddEntry(dummy, "True peak positions", "l");
    leg->SetBit(kCanDelete);
    leg->Draw();

    TLatex* info = new TLatex();
    info->SetNDC(); info->SetTextSize(0.030);
    info->SetBit(kCanDelete);
    info->DrawLatex(0.10, 0.94,
        Form("%d peaks  |  FWHM(1 MeV) #approx %.1f keV  |  BG: %s",
             nPeaks,
             TMath::Sqrt(a_fwhm + b_fwhm * 1000.0) * 2.3548,
             bgNames[bgType]));

    c->Modified(); c->Update();
    SetStatus(Form("Synthetic preview: %s BG, %d peaks", bgNames[bgType], nPeaks));
}
