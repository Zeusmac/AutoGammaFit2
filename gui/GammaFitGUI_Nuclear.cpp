#include "GammaFitGUI.h"
#include "GammaFitGUI_shared.h"
#include "NNDCFetcher.h"

#include "TGTextEntry.h"
#include "TGFileDialog.h"
#include "TGMsgBox.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TLine.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TBox.h"
#include "TROOT.h"

#include <algorithm>
#include <deque>
#include <dirent.h>
#include <fstream>
#include <sstream>
#include <set>
#include <sys/stat.h>
#include <unistd.h>

// ─────────────────────────────────────────────────────────────────────────────
// NucCacheDirPath — returns the nuclear data cache directory path
// ─────────────────────────────────────────────────────────────────────────────
std::string GammaFitGUI::NucCacheDirPath() const
{
    return nucCacheDir_;
}

// ─────────────────────────────────────────────────────────────────────────────
// PopulateNucGammaRef — fill nucGammaRefView_ with gamma lines for chain
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::PopulateNucGammaRef()
{
    if (!nucGammaRefView_) return;
    nucGammaRefView_->Clear();

    for (const auto& isoID : nucChainIsotopes_) {
        auto it = nuclearDB_.find(isoID);
        if (it == nuclearDB_.end() || !it->second.valid()) continue;
        const NucIsotope& iso = it->second;

        // Header line
        std::string cls = labelClassMap_.count(isoID) ? labelClassMap_.at(isoID) : "?";
        std::string header = isoID + "  [" + cls + "]  T½=" + iso.hl_str;
        nucGammaRefView_->AddLine(header.c_str());

        // Gamma lines sorted by energy
        std::vector<NucGamma> sorted = iso.gammas;
        std::sort(sorted.begin(), sorted.end(),
            [](const NucGamma& a, const NucGamma& b){ return a.energy < b.energy; });

        for (const auto& gm : sorted) {
            if (gm.energy <= 0) continue;
            char buf[128];
            snprintf(buf, sizeof(buf), "  %.2f keV  I=%.1f%%", gm.energy, gm.intensity);
            nucGammaRefView_->AddLine(buf);
        }
        if (sorted.empty()) nucGammaRefView_->AddLine("  (no gamma data)");
        nucGammaRefView_->AddLine("");
    }
    nucGammaRefView_->Update();
}

// ─────────────────────────────────────────────────────────────────────────────
// RefreshIsoComboHelper — repopulate a TGComboBox with chain isotope names
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::RefreshIsoComboHelper(TGComboBox* combo,
                                         const std::vector<std::string>& chain,
                                         const std::map<std::string, NucIsotope>& /*db*/)
{
    if (!combo) return;
    // Remember current selection
    std::string curSel;
    if (TGLBEntry* e = combo->GetSelectedEntry()) curSel = e->GetTitle();

    combo->RemoveAll();
    combo->AddEntry("(none)", 1);
    int id = 2, restoreId = 1;
    for (const auto& isoID : chain) {
        if (isoID == curSel) restoreId = id;
        combo->AddEntry(isoID.c_str(), id++);
    }
    combo->Select(restoreId, kFALSE);
    combo->MapSubwindows();
    combo->Layout();
}

// ─────────────────────────────────────────────────────────────────────────────
// DrawLevelScheme — draw nuclear level scheme for a single isotope
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::DrawLevelScheme(const std::string& isoID)
{
    auto it = nuclearDB_.find(isoID);
    if (it == nuclearDB_.end() || !it->second.valid()) {
        AppendLog("Nuclear: no data for " + isoID); return;
    }
    const NucIsotope& iso = it->second;

    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();
    c->SetFillColor(kWhite);

    if (iso.levels.empty()) {
        AppendLog("Nuclear: no level data for " + isoID); return;
    }

    // Find energy range
    double eMax = 0;
    for (const auto& lv : iso.levels)
        if (lv.energy > eMax) eMax = lv.energy;
    if (eMax <= 0) eMax = 100.0;

    const double xL  = 0.15, xR = 0.75;
    const double yBot = 0.08, yTop = 0.92;

    TLatex tx; tx.SetNDC(kTRUE); tx.SetTextAlign(12); tx.SetTextSize(0.022);

    // Draw each level
    for (const auto& lv : iso.levels) {
        double y = yBot + (yTop - yBot) * (lv.energy / eMax);
        TLine* ln = new TLine(xL, y, xR, y);
        ln->SetLineColor(kBlue + 1);
        ln->SetLineWidth(1);
        ln->Draw();
        // label
        char buf[64];
        snprintf(buf, sizeof(buf), "%.1f  %s", lv.energy, lv.jpi.c_str());
        tx.DrawLatex(xR + 0.01, y, buf);
    }

    // Draw gamma transitions
    for (const auto& gm : iso.gammas) {
        if (gm.energy <= 0) continue;
        double yFrom = yBot + (yTop - yBot) * (gm.start_level / eMax);
        double eTo   = gm.start_level - gm.energy;
        double yTo   = yBot + (yTop - yBot) * (std::max(eTo, 0.0) / eMax);
        double xMid  = (xL + xR) / 2.0 + 0.05;
        TArrow* arr = new TArrow(xMid, yFrom, xMid, yTo, 0.010, "|>");
        arr->SetLineColor(kRed);
        arr->SetFillColor(kRed);
        arr->Draw();
    }

    // Title
    tx.SetTextSize(0.030); tx.SetTextAlign(22); tx.SetTextColor(kBlack);
    tx.DrawLatex(0.50, 0.96, (isoID + " Level Scheme").c_str());

    c->Modified(); c->Update();

    // If user level scheme is loaded for this isotope, delegate to enhanced drawing
    if (!lsData_.empty() && lsData_.isoID == isoID) {
        DrawEnhancedLevelScheme(isoID);
        return;
    }

    AppendLog("Level scheme drawn for " + isoID + " (NNDC data)");
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildNuclearTab — main Nuclear tab with 6 sub-tabs
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::BuildNuclearTab(TGCompositeFrame* p)
{
    TGTab* nucTab = new TGTab(p, 308, 860);
    p->AddFrame(nucTab, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));

    // ═══════════════════════════════════════════════════════════════════════
    // Sub-tab 1: Isotope Chain
    // ═══════════════════════════════════════════════════════════════════════
    {
        TGCompositeFrame* tab = nucTab->AddTab("Isotope Chain");
        TGCanvas* sc = new TGCanvas(tab, 306, 860, kSunkenFrame);
        tab->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        TGCompositeFrame* cf = new TGCompositeFrame(sc->GetViewPort(), 293, 10, kVerticalFrame);
        sc->SetContainer(cf);
        TGCompositeFrame* p1 = cf;

        // ── Group: Parent Nucleus ──────────────────────────────────────────
        TGGroupFrame* parentGrp = new TGGroupFrame(p1, "Parent Nucleus");
        p1->AddFrame(parentGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

        // Row 1: A, Symbol, Set Parent, Trace Chain
        {
            TGHorizontalFrame* row = new TGHorizontalFrame(parentGrp);
            parentGrp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 1));

            row->AddFrame(new TGLabel(row, "A:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
            nucAEntry_ = new TGNumberEntry(row, 44, 4, -1,
                TGNumberFormat::kNESInteger,
                TGNumberFormat::kNEAPositive,
                TGNumberFormat::kNELLimitMinMax, 1, 300);
            nucAEntry_->Resize(50, 22);
            row->AddFrame(nucAEntry_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));

            row->AddFrame(new TGLabel(row, "Sym:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
            nucSymbolEntry_ = new TGTextEntry(row, "S");
            nucSymbolEntry_->Resize(36, 22);
            row->AddFrame(nucSymbolEntry_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));

            TGTextButton* setParentBtn = new TGTextButton(row, "Set Parent");
            row->AddFrame(setParentBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
            setParentBtn->Connect("Clicked()", "GammaFitGUI", this,
                                  "OnNucSetParentFromChain()");

            nucTraceChainBtn_ = new TGTextButton(row, "Trace Chain");
            row->AddFrame(nucTraceChainBtn_, new TGLayoutHints(kLHintsLeft, 2, 0, 0, 0));
            nucTraceChainBtn_->Connect("Clicked()", "GammaFitGUI", this,
                                       "OnNucAutoTraceChain()");
            nucTraceChainBtn_->SetToolTipText(
                "Follow the beta-minus decay chain automatically from the parent nucleus.\n"
                "Fetches data from IAEA LiveChart (or cache) for each step.");
        }

        // Row 2: Z, N, Set Z/N, Name
        {
            TGHorizontalFrame* row = new TGHorizontalFrame(parentGrp);
            parentGrp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 2));

            row->AddFrame(new TGLabel(row, "Z:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
            isoParentZEntry_ = new TGNumberEntry(row, 0, 4, -1,
                TGNumberFormat::kNESInteger,
                TGNumberFormat::kNEANonNegative,
                TGNumberFormat::kNELLimitMinMax, 0, 120);
            isoParentZEntry_->Resize(45, 22);
            row->AddFrame(isoParentZEntry_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));

            row->AddFrame(new TGLabel(row, "N:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
            isoParentNEntry_ = new TGNumberEntry(row, 0, 4, -1,
                TGNumberFormat::kNESInteger,
                TGNumberFormat::kNEANonNegative,
                TGNumberFormat::kNELLimitMinMax, 0, 200);
            isoParentNEntry_->Resize(45, 22);
            row->AddFrame(isoParentNEntry_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));

            TGTextButton* setZNBtn = new TGTextButton(row, "Set Z/N");
            row->AddFrame(setZNBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            setZNBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoSetParent()");
        }

        // Row 3: Name entry (hidden / programmatic)
        {
            TGHorizontalFrame* row = new TGHorizontalFrame(parentGrp);
            parentGrp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 2));
            row->AddFrame(new TGLabel(row, "Name:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            isoParentNameEntry_ = new TGTextEntry(row, "");
            isoParentNameEntry_->Resize(100, 22);
            isoParentNameEntry_->SetToolTipText("Freeform isotope name, e.g. 44S");
            row->AddFrame(isoParentNameEntry_, new TGLayoutHints(kLHintsExpandX));
        }

        // ── Group: Decay Chain ─────────────────────────────────────────────
        TGGroupFrame* chainGrp = new TGGroupFrame(p1, "Decay Chain");
        p1->AddFrame(chainGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        nucChainList_ = new TGListBox(chainGrp, 1100);
        nucChainList_->Resize(285, 100);
        chainGrp->AddFrame(nucChainList_,
            new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        {
            TGHorizontalFrame* btnRow = new TGHorizontalFrame(chainGrp);
            chainGrp->AddFrame(btnRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
            TGTextButton* removeBtn = new TGTextButton(btnRow, "Remove");
            btnRow->AddFrame(removeBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            removeBtn->Connect("Clicked()", "GammaFitGUI", this, "OnNucRemoveFromChain()");
            TGTextButton* clearBtn = new TGTextButton(btnRow, "Clear All");
            btnRow->AddFrame(clearBtn, new TGLayoutHints(kLHintsLeft));
            clearBtn->Connect("Clicked()", "GammaFitGUI", this, "OnNucClearChain()");
        }

        chainGrp->AddFrame(new TGLabel(chainGrp, "Add manually:"),
                           new TGLayoutHints(kLHintsLeft, 2, 2, 2, 1));
        {
            TGHorizontalFrame* addRow = new TGHorizontalFrame(chainGrp);
            chainGrp->AddFrame(addRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));

            addRow->AddFrame(new TGLabel(addRow, "A:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
            nucAddAEntry_ = new TGNumberEntry(addRow, 44, 4, -1,
                TGNumberFormat::kNESInteger,
                TGNumberFormat::kNEAPositive,
                TGNumberFormat::kNELLimitMinMax, 1, 300);
            nucAddAEntry_->Resize(45, 22);
            addRow->AddFrame(nucAddAEntry_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));

            addRow->AddFrame(new TGLabel(addRow, "Sym:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
            nucAddSymEntry_ = new TGTextEntry(addRow, "");
            nucAddSymEntry_->Resize(36, 22);
            addRow->AddFrame(nucAddSymEntry_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));

            nucAddTypeCombo_ = new TGComboBox(addRow, 1110);
            nucAddTypeCombo_->AddEntry("(auto)", 1);
            nucAddTypeCombo_->AddEntry("Decay chain", 2);
            nucAddTypeCombo_->AddEntry("Background", 3);
            nucAddTypeCombo_->Select(1, kFALSE);
            nucAddTypeCombo_->Resize(90, 22);
            addRow->AddFrame(nucAddTypeCombo_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));

            TGTextButton* addBtn = new TGTextButton(addRow, "Add");
            addRow->AddFrame(addBtn, new TGLayoutHints(kLHintsLeft));
            addBtn->Connect("Clicked()", "GammaFitGUI", this, "OnNucAddToChain()");
        }

        // ── Group: Common Backgrounds ──────────────────────────────────────
        TGGroupFrame* bgGrp = new TGGroupFrame(p1, "Common Backgrounds");
        p1->AddFrame(bgGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        bgGrp->AddFrame(new TGLabel(bgGrp, "Quick add background isotopes:"),
                        new TGLayoutHints(kLHintsLeft, 2, 2, 2, 1));
        {
            TGHorizontalFrame* bgRow = new TGHorizontalFrame(bgGrp);
            bgGrp->AddFrame(bgRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 2));

            bgRow->AddFrame(new TGLabel(bgRow, "Isotope:"),
                            new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            nucBgCombo_ = new TGComboBox(bgRow, 1120);
            nucBgCombo_->AddEntry("40K",   1);
            nucBgCombo_->AddEntry("208Tl", 2);
            nucBgCombo_->AddEntry("214Bi", 3);
            nucBgCombo_->AddEntry("214Pb", 4);
            nucBgCombo_->AddEntry("137Cs", 5);
            nucBgCombo_->AddEntry("60Co",  6);
            nucBgCombo_->AddEntry("56Co",  7);
            nucBgCombo_->AddEntry("22Na",  8);
            nucBgCombo_->AddEntry("207Bi", 9);
            nucBgCombo_->AddEntry("241Am", 10);
            nucBgCombo_->AddEntry("152Eu", 11);
            nucBgCombo_->Select(1, kFALSE);
            nucBgCombo_->Resize(80, 22);
            bgRow->AddFrame(nucBgCombo_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));

            TGTextButton* addBgBtn = new TGTextButton(bgRow, "Add Background");
            bgRow->AddFrame(addBgBtn, new TGLayoutHints(kLHintsLeft));
            addBgBtn->Connect("Clicked()", "GammaFitGUI", this, "OnNucAddBackground()");
        }

        // ── Fetch / Load section ───────────────────────────────────────────
        TGGroupFrame* fetchGrp = new TGGroupFrame(p1, "Fetch / Load");
        p1->AddFrame(fetchGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        {
            TGHorizontalFrame* fetchRow = new TGHorizontalFrame(fetchGrp);
            fetchGrp->AddFrame(fetchRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 1));

            TGTextButton* fetchBtn = new TGTextButton(fetchRow, "Fetch NNDC");
            fetchRow->AddFrame(fetchBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
            fetchBtn->Connect("Clicked()", "GammaFitGUI", this, "OnNucFetchAll()");
            fetchBtn->SetToolTipText("Fetch data for all chain isotopes from IAEA LiveChart");

            TGTextButton* reloadBtn = new TGTextButton(fetchRow, "Reload Cache");
            fetchRow->AddFrame(reloadBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
            reloadBtn->Connect("Clicked()", "GammaFitGUI", this, "OnNucReloadCache()");

            TGTextButton* loadTxtBtn = new TGTextButton(fetchRow, "Load NNDC txt...");
            fetchRow->AddFrame(loadTxtBtn, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
            loadTxtBtn->Connect("Clicked()", "GammaFitGUI", this, "OnNucLoadNNDCTxt()");
            loadTxtBtn->SetToolTipText(
                "Load a locally-saved IAEA LiveChart CSV file for the isotope\n"
                "specified by the A / Sym entries above.");
        }
        nucStatusLbl_ = new TGLabel(fetchGrp, "");
        fetchGrp->AddFrame(nucStatusLbl_, new TGLayoutHints(kLHintsLeft, 2, 2, 1, 2));

        // ── Group: AME & NUBASE ────────────────────────────────────────────
        TGGroupFrame* ameGrp = new TGGroupFrame(p1, "AME & NUBASE");
        p1->AddFrame(ameGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        {
            TGHorizontalFrame* ameRow = new TGHorizontalFrame(ameGrp);
            ameGrp->AddFrame(ameRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 1));
            TGTextButton* ameBtn = new TGTextButton(ameRow, "Load AME Table");
            ameRow->AddFrame(ameBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            ameBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLoadAMETable()");
            ameLbl_ = new TGLabel(ameRow, "(not loaded)");
            ameRow->AddFrame(ameLbl_, new TGLayoutHints(kLHintsCenterY));
        }
        {
            TGHorizontalFrame* nbRow = new TGHorizontalFrame(ameGrp);
            ameGrp->AddFrame(nbRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 2));
            TGTextButton* nbBtn = new TGTextButton(nbRow, "Load NUBASE Table");
            nbRow->AddFrame(nbBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            nbBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLoadNubaseTable()");
            nubaseLbl_ = new TGLabel(nbRow, "(not loaded)");
            nbRow->AddFrame(nubaseLbl_, new TGLayoutHints(kLHintsCenterY));
        }

        // ── Group: Known Gamma Lines ───────────────────────────────────────
        TGGroupFrame* gammaRefGrp = new TGGroupFrame(p1, "Known Gamma Lines");
        p1->AddFrame(gammaRefGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

        nucGammaRefView_ = new TGTextView(gammaRefGrp, 285, 120);
        gammaRefGrp->AddFrame(nucGammaRefView_,
            new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        {
            TGTextButton* confirmBtn = new TGTextButton(gammaRefGrp, "Confirm Chain → Isotope DB");
            gammaRefGrp->AddFrame(confirmBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            confirmBtn->Connect("Clicked()", "GammaFitGUI", this, "OnNucConfirmChainToIsoDB()");
            confirmBtn->SetToolTipText(
                "Load all gamma lines for the current decay chain into the Isotope DB\n"
                "used by the Peak Matching tab.  Existing DB entries are preserved.");
        }
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Sub-tab 2: Peak Matching (previously the Isotopes tab)
    // ═══════════════════════════════════════════════════════════════════════
    {
        TGCompositeFrame* tab = nucTab->AddTab("Peak Matching");
        TGCanvas* sc = new TGCanvas(tab, 306, 860, kSunkenFrame);
        tab->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        TGCompositeFrame* cf = new TGCompositeFrame(sc->GetViewPort(), 293, 10, kVerticalFrame);
        sc->SetContainer(cf);
        TGCompositeFrame* p2 = cf;

        // ── Isotope DB load row ────────────────────────────────────────────
        {
            TGGroupFrame* isoDbLoadGrp = new TGGroupFrame(p2, "Isotope Database");
            p2->AddFrame(isoDbLoadGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

            TGHorizontalFrame* isoRow = new TGHorizontalFrame(isoDbLoadGrp);
            isoDbLoadGrp->AddFrame(isoRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 1));
            TGTextButton* isoBtn = new TGTextButton(isoRow, "Open Isotope DB...");
            isoRow->AddFrame(isoBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            isoBtn->Connect("Clicked()", "GammaFitGUI", this, "OnOpenIsotopeDB()");
            isoBtn->SetToolTipText("Browse for and load an isotope energy database text file");
            TGTextButton* reloadBtn = new TGTextButton(isoRow, "Reload");
            isoRow->AddFrame(reloadBtn, new TGLayoutHints(kLHintsLeft));
            reloadBtn->Connect("Clicked()", "GammaFitGUI", this, "OnReloadIsotopeDB()");
            reloadBtn->SetToolTipText("Reload the database from the current file");

            isotopeLbl_ = new TGLabel(isoDbLoadGrp, "(not loaded)");
            isoDbLoadGrp->AddFrame(isotopeLbl_, new TGLayoutHints(kLHintsLeft, 2, 2, 0, 2));
        }

        // ── Controls row ──────────────────────────────────────────────────
        {
            TGHorizontalFrame* topRow = new TGHorizontalFrame(p2);
            p2->AddFrame(topRow, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

            TGTextButton* refreshBtn = new TGTextButton(topRow, "Refresh");
            topRow->AddFrame(refreshBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            refreshBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoRefresh()");
            refreshBtn->SetToolTipText("Reload cache for the current histogram and repopulate the list");

            topRow->AddFrame(new TGLabel(topRow, "Filter:"),
                             new TGLayoutHints(kLHintsCenterY, 4, 4, 0, 0));
            isoFilterCombo_ = new TGComboBox(topRow, 950);
            isoFilterCombo_->AddEntry("All", 1);
            isoFilterCombo_->Select(1, kFALSE);
            isoFilterCombo_->Resize(110, 22);
            topRow->AddFrame(isoFilterCombo_, new TGLayoutHints(kLHintsExpandX));
            isoFilterCombo_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                     "OnIsoFilterChanged(Int_t)");
        }

        // ── Peak list ─────────────────────────────────────────────────────
        TGGroupFrame* listGrp = new TGGroupFrame(p2, "Peaks (sorted by energy)");
        p2->AddFrame(listGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        isoList_ = new TGListBox(listGrp, 960);
        isoList_->Resize(285, 200);
        listGrp->AddFrame(isoList_,
                          new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        isoList_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                          "OnIsoListSelected(Int_t)");

        {
            TGTextButton* previewBtn = new TGTextButton(listGrp, "Preview Selected Peak");
            listGrp->AddFrame(previewBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
            previewBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoPeakPreview()");
            previewBtn->SetToolTipText("Zoom canvas to the selected peak");
        }

        // ── Bulk action row ───────────────────────────────────────────────
        {
            TGHorizontalFrame* bulkRow = new TGHorizontalFrame(p2);
            p2->AddFrame(bulkRow, new TGLayoutHints(kLHintsExpandX, 4, 4, 0, 2));

            TGTextButton* rematchBtn = new TGTextButton(bulkRow, "Auto Re-match All");
            bulkRow->AddFrame(rematchBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            rematchBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoAutoMatchAll()");
            rematchBtn->SetToolTipText(
                "For every unlabeled peak in the cache, find the closest unused\n"
                "DB line (greedy unique match) and apply it as the label.");

            TGTextButton* clearAllBtn = new TGTextButton(bulkRow, "Clear All Labels");
            bulkRow->AddFrame(clearAllBtn, new TGLayoutHints(kLHintsLeft));
            clearAllBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoClearAll()");
            clearAllBtn->SetToolTipText("Remove label and classification from every peak in the cache");
        }

        // ── Match threshold ────────────────────────────────────────────────
        {
            TGHorizontalFrame* threshRow = new TGHorizontalFrame(p2);
            p2->AddFrame(threshRow, new TGLayoutHints(kLHintsExpandX, 4, 4, 0, 2));
            threshRow->AddFrame(new TGLabel(threshRow, "Match threshold (keV):"),
                                new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            isoMatchThreshEntry_ = new TGNumberEntry(threshRow, 3.0, 5, -1,
                TGNumberFormat::kNESRealOne,
                TGNumberFormat::kNEAPositive,
                TGNumberFormat::kNELLimitMinMax, 0.1, 50.0);
            threshRow->AddFrame(isoMatchThreshEntry_, new TGLayoutHints(kLHintsLeft));
        }

        {
            TGTextButton* applyMatchBtn = new TGTextButton(p2, "Load Matches to Histogram");
            p2->AddFrame(applyMatchBtn, new TGLayoutHints(kLHintsExpandX, 4, 4, 0, 2));
            applyMatchBtn->Connect("Clicked()", "GammaFitGUI", this,
                                   "OnIsoApplyAutoMatches()");
            applyMatchBtn->SetToolTipText("Apply all auto-matched labels to the cache");
        }

        // ── Edit group ────────────────────────────────────────────────────
        TGGroupFrame* editGrp = new TGGroupFrame(p2, "Edit Selected Peak");
        p2->AddFrame(editGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        {
            TGHorizontalFrame* lblRow = new TGHorizontalFrame(editGrp);
            editGrp->AddFrame(lblRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 0));
            lblRow->AddFrame(new TGLabel(lblRow, "Isotope:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            isoLabelCombo_ = new TGComboBox(lblRow, 961);
            isoLabelCombo_->AddEntry("(none)", 1);
            isoLabelCombo_->Select(1, kFALSE);
            isoLabelCombo_->Resize(200, 22);
            lblRow->AddFrame(isoLabelCombo_, new TGLayoutHints(kLHintsExpandX));
        }
        {
            TGHorizontalFrame* custLblRow = new TGHorizontalFrame(editGrp);
            editGrp->AddFrame(custLblRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
            custLblRow->AddFrame(new TGLabel(custLblRow, "Custom:"),
                                 new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            isoCustomLabelEntry_ = new TGTextEntry(custLblRow, "");
            isoCustomLabelEntry_->SetToolTipText(
                "Free-text label — overrides the combo above when non-empty");
            custLblRow->AddFrame(isoCustomLabelEntry_, new TGLayoutHints(kLHintsExpandX));
        }

        // Class row
        {
            TGHorizontalFrame* clsRow = new TGHorizontalFrame(editGrp);
            editGrp->AddFrame(clsRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 0));
            clsRow->AddFrame(new TGLabel(clsRow, "Class:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            isoClassCombo_ = new TGComboBox(clsRow, 970);
            isoClassCombo_->AddEntry("(none)",                         1);
            isoClassCombo_->AddEntry("Parent",                         2);
            isoClassCombo_->AddEntry("Daughter",                       3);
            isoClassCombo_->AddEntry("Granddaughter",                  4);
            isoClassCombo_->AddEntry("Great-Granddaughter",            5);
            isoClassCombo_->AddEntry("Great-Great-Granddaughter",      6);
            isoClassCombo_->AddEntry("Beta-n Daughter",                7);
            isoClassCombo_->AddEntry("Beta-2n Daughter",               8);
            isoClassCombo_->AddEntry("Beta-n Granddaughter",           9);
            isoClassCombo_->AddEntry("Beta-2n Granddaughter",          10);
            isoClassCombo_->AddEntry("Beta-n Great-Granddaughter",     11);
            isoClassCombo_->AddEntry("Beta-2n Great-Granddaughter",    12);
            isoClassCombo_->AddEntry("Beta-n Great-Great-Granddaughter",  13);
            isoClassCombo_->AddEntry("Beta-2n Great-Great-Granddaughter", 14);
            isoClassCombo_->AddEntry("Background",                     15);
            isoClassCombo_->AddEntry("Custom",                         16);
            isoClassCombo_->AddEntry("X-ray",                          17);
            isoClassCombo_->Select(1, kFALSE);
            isoClassCombo_->Resize(150, 22);
            clsRow->AddFrame(isoClassCombo_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            isoCustomEntry_ = new TGTextEntry(clsRow, "");
            isoCustomEntry_->SetToolTipText("Custom class name");
            clsRow->AddFrame(isoCustomEntry_, new TGLayoutHints(kLHintsExpandX));
        }

        TGTextButton* applyBtn = new TGTextButton(editGrp, "Apply Isotope to Selected Peak");
        editGrp->AddFrame(applyBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 1));
        applyBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoApply()");
        applyBtn->SetToolTipText("Save label + class to the selected peak's cache entry");

        TGTextButton* applyClsAllBtn = new TGTextButton(editGrp,
            "Apply Class to ALL peaks with this Label");
        editGrp->AddFrame(applyClsAllBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
        applyClsAllBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoApplyClassToAll()");
        applyClsAllBtn->SetToolTipText(
            "Apply the selected class to every peak in the cache that\n"
            "shares the currently selected label.");

        TGTextButton* applyAllBtn = new TGTextButton(editGrp, "Apply Label to All in Filter");
        editGrp->AddFrame(applyAllBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
        applyAllBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoApplyLabelAll()");
        applyAllBtn->SetToolTipText("Apply the current label to every peak currently shown in the list");

        TGTextButton* clearBtn = new TGTextButton(editGrp, "Clear Isotope");
        editGrp->AddFrame(clearBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 4));
        clearBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoClear()");
        clearBtn->SetToolTipText("Remove label and classification from the selected peak's cache entry");

        // ── Decay schematic ────────────────────────────────────────────────
        TGTextButton* schemBtn = new TGTextButton(p2, "Draw Decay Schematic on Canvas");
        p2->AddFrame(schemBtn, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));
        schemBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoDrawSchematic()");
        schemBtn->SetToolTipText(
            "Draw a decay chain diagram on the main canvas.\n"
            "Nodes are classes (Parent → Daughter → Granddaughter).\n"
            "Only classes with labeled peaks are shown.");

        // ── Isotope Database browser ───────────────────────────────────────
        TGGroupFrame* dbGrp = new TGGroupFrame(p2, "Isotope Database");
        p2->AddFrame(dbGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 4));

        {
            TGHorizontalFrame* searchRow = new TGHorizontalFrame(dbGrp);
            dbGrp->AddFrame(searchRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 2));
            searchRow->AddFrame(new TGLabel(searchRow, "Search:"),
                                new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            isoDbSearch_ = new TGTextEntry(searchRow, "");
            isoDbSearch_->SetToolTipText("Filter by isotope name or energy (partial match)");
            searchRow->AddFrame(isoDbSearch_, new TGLayoutHints(kLHintsExpandX, 0, 4, 0, 0));
            TGTextButton* srchBtn = new TGTextButton(searchRow, "Filter");
            searchRow->AddFrame(srchBtn, new TGLayoutHints(kLHintsLeft));
            srchBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoDbSearch()");
            TGTextButton* clearDbBtn = new TGTextButton(searchRow, "Clear");
            searchRow->AddFrame(clearDbBtn, new TGLayoutHints(kLHintsLeft, 4, 0, 0, 0));
            clearDbBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoDbClear()");
        }

        isoDbList_ = new TGListBox(dbGrp, 980);
        isoDbList_->Resize(285, 180);
        dbGrp->AddFrame(isoDbList_,
                        new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
        isoDbList_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                            "OnIsoDbLineSelected(Int_t)");

        {
            TGHorizontalFrame* clsRow = new TGHorizontalFrame(dbGrp);
            dbGrp->AddFrame(clsRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
            clsRow->AddFrame(new TGLabel(clsRow, "Class:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            isoDbClassCombo_ = new TGComboBox(clsRow, 985);
            isoDbClassCombo_->AddEntry("(none)",               1);
            isoDbClassCombo_->AddEntry("Parent",               2);
            isoDbClassCombo_->AddEntry("Daughter",             3);
            isoDbClassCombo_->AddEntry("Granddaughter",        4);
            isoDbClassCombo_->AddEntry("Beta-n Daughter",      5);
            isoDbClassCombo_->AddEntry("Beta-2n Daughter",     6);
            isoDbClassCombo_->AddEntry("Beta-n Granddaughter", 7);
            isoDbClassCombo_->AddEntry("Beta-2n Granddaughter",8);
            isoDbClassCombo_->AddEntry("Background",           9);
            isoDbClassCombo_->AddEntry("Custom",               10);
            isoDbClassCombo_->AddEntry("X-ray",                11);
            isoDbClassCombo_->Select(1, kFALSE);
            isoDbClassCombo_->Resize(150, 22);
            clsRow->AddFrame(isoDbClassCombo_, new TGLayoutHints(kLHintsLeft));
        }

        TGTextButton* dbApplyBtn = new TGTextButton(dbGrp,
            "Apply Selected DB Line to Fitted Peak");
        dbGrp->AddFrame(dbApplyBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));
        dbApplyBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoDbApply()");
        dbApplyBtn->SetToolTipText(
            "Apply the isotope name from the selected DB line as label,\n"
            "and the chosen class, to the peak selected in the Peaks list above.");
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Sub-tab 3: Level Scheme (editor + enhanced drawing)
    // ═══════════════════════════════════════════════════════════════════════
    {
        TGCompositeFrame* tab = nucTab->AddTab("Level Scheme");
        TGCanvas* sc = new TGCanvas(tab, 306, 860, kSunkenFrame);
        tab->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        TGCompositeFrame* cf = new TGCompositeFrame(sc->GetViewPort(), 293, 10, kVerticalFrame);
        sc->SetContainer(cf);
        BuildLSLevelSchemeTab(cf);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Sub-tab 4: Log ft & Branching Ratios
    // ═══════════════════════════════════════════════════════════════════════
    {
        TGCompositeFrame* tab = nucTab->AddTab("Log ft & BR");
        TGCanvas* sc = new TGCanvas(tab, 306, 860, kSunkenFrame);
        tab->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        TGCompositeFrame* cf = new TGCompositeFrame(sc->GetViewPort(), 293, 10, kVerticalFrame);
        sc->SetContainer(cf);
        BuildLSLogFtTab(cf);
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Sub-tab 5: gamma-gamma (TODO)
    // ═══════════════════════════════════════════════════════════════════════
    {
        TGCompositeFrame* tab = nucTab->AddTab("γ-γ (TODO)");
        tab->AddFrame(new TGLabel(tab, "Gamma-gamma coincidence analysis coming soon."),
                      new TGLayoutHints(kLHintsLeft, 4, 4, 8, 2));
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Sub-tab 6: Shell Model
    // ═══════════════════════════════════════════════════════════════════════
    {
        TGCompositeFrame* tab = nucTab->AddTab("Shell Model");
        tab->AddFrame(new TGLabel(tab, "Shell model comparisons coming soon."),
                      new TGLayoutHints(kLHintsLeft, 4, 4, 8, 2));
    }

    // ═══════════════════════════════════════════════════════════════════════
    // Sub-tab 7: References
    // ═══════════════════════════════════════════════════════════════════════
    {
        TGCompositeFrame* tab = nucTab->AddTab("References");
        TGCanvas* sc = new TGCanvas(tab, 306, 860, kSunkenFrame);
        tab->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
        TGCompositeFrame* cf = new TGCompositeFrame(sc->GetViewPort(), 293, 10, kVerticalFrame);
        sc->SetContainer(cf);

        // ── Filter row ────────────────────────────────────────────────────
        {
            TGGroupFrame* filterGrp = new TGGroupFrame(cf, "Filter");
            cf->AddFrame(filterGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));
            TGHorizontalFrame* row = new TGHorizontalFrame(filterGrp);
            filterGrp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

            row->AddFrame(new TGLabel(row, "Decay class:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            nucRefDecayFilter_ = new TGComboBox(row, 3700);
            nucRefDecayFilter_->AddEntry("All",                1);
            nucRefDecayFilter_->AddEntry("Parent",             2);
            nucRefDecayFilter_->AddEntry("Daughter",           3);
            nucRefDecayFilter_->AddEntry("Granddaughter",      4);
            nucRefDecayFilter_->AddEntry("Beta-n Daughter",    5);
            nucRefDecayFilter_->AddEntry("Beta-2n Daughter",   6);
            nucRefDecayFilter_->AddEntry("Other",              7);
            nucRefDecayFilter_->Resize(140, 22);
            nucRefDecayFilter_->Select(1, kFALSE);
            row->AddFrame(nucRefDecayFilter_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            nucRefDecayFilter_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                        "OnNucRefFilterChanged(Int_t)");
        }

        // ── Gamma reference list ──────────────────────────────────────────
        {
            TGGroupFrame* listGrp = new TGGroupFrame(cf, "Gamma Rays & References  (select → Open NNDC)");
            cf->AddFrame(listGrp, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 4, 4, 2, 2));

            nucRefList_ = new TGListBox(listGrp, 3710);
            nucRefList_->Resize(280, 420);
            nucRefList_->SetMultipleSelections(kFALSE);
            listGrp->AddFrame(nucRefList_,
                new TGLayoutHints(kLHintsExpandX | kLHintsExpandY, 2, 2, 2, 2));
            nucRefList_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                 "OnNucRefFilterChanged(Int_t)");

            // Detail label for selected row
            nucRefDetailLbl_ = new TGLabel(listGrp, "Select a gamma line above");
            nucRefDetailLbl_->SetTextJustify(kTextLeft);
            listGrp->AddFrame(nucRefDetailLbl_,
                new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

            // Button row
            TGHorizontalFrame* btnRow = new TGHorizontalFrame(listGrp);
            listGrp->AddFrame(btnRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));

            TGTextButton* openBtn = new TGTextButton(btnRow, "Open NNDC ↗");
            btnRow->AddFrame(openBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
            openBtn->Connect("Clicked()", "GammaFitGUI", this, "OnNucRefOpen()");
            openBtn->SetToolTipText("Open NuDat page for the selected isotope in your browser");

            TGTextButton* refreshBtn = new TGTextButton(btnRow, "Refresh");
            btnRow->AddFrame(refreshBtn, new TGLayoutHints(kLHintsLeft, 0, 0, 0, 0));
            refreshBtn->Connect("Clicked()", "GammaFitGUI", this, "PopulateNucRefTab()");
        }

        // ── Histogram → Parent nucleus ────────────────────────────────────
        {
            TGGroupFrame* histGrp = new TGGroupFrame(cf, "Assign Histogram to Parent Nucleus");
            cf->AddFrame(histGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

            TGHorizontalFrame* row = new TGHorizontalFrame(histGrp);
            histGrp->AddFrame(row, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

            row->AddFrame(new TGLabel(row, "Histogram:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            // Reuse nucCombo selector — we'll populate from the hist list at fill time
            // Use a plain text entry so any name can be typed
            TGTextEntry* histEntry = new TGTextEntry(row, "");
            histEntry->SetName("nucHistParentHistEntry");
            histEntry->Resize(110, 22);
            row->AddFrame(histEntry, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));

            row->AddFrame(new TGLabel(row, "Parent:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
            TGTextEntry* parentEntry = new TGTextEntry(row, "");
            parentEntry->SetName("nucHistParentEntry");
            parentEntry->Resize(60, 22);
            row->AddFrame(parentEntry, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));

            TGTextButton* setBtn = new TGTextButton(row, "Set");
            row->AddFrame(setBtn, new TGLayoutHints(kLHintsLeft));
            setBtn->Connect("Clicked()", "GammaFitGUI", this, "OnNucSetHistParent()");
            setBtn->SetToolTipText(
                "Associate the named histogram with a parent nucleus.\n"
                "Histograms sharing the same parent share the chain cache.");
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// OnNucSetParentFromChain
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnNucSetParentFromChain()
{
    if (!nucAEntry_ || !nucSymbolEntry_) return;
    int A = (int)nucAEntry_->GetNumber();
    std::string sym = nucSymbolEntry_->GetText();
    while (!sym.empty() && sym.front() == ' ') sym = sym.substr(1);
    while (!sym.empty() && sym.back()  == ' ') sym.pop_back();
    if (sym.empty() || A <= 0) { AppendLog("Nuclear: invalid parent entry"); return; }
    if (sym.size() >= 1) sym[0] = (char)toupper((unsigned char)sym[0]);
    if (sym.size() >= 2) sym[1] = (char)tolower((unsigned char)sym[1]);

    int Z = NucSymbolToZ(sym);
    if (Z <= 0) { AppendLog("Nuclear: unknown element symbol: " + sym); return; }
    int N = A - Z;

    isoParentZval_   = Z;
    isoParentNval_   = N;
    isoParentIsotope_ = std::to_string(A) + sym;
    if (isoParentNameEntry_) isoParentNameEntry_->SetText(isoParentIsotope_.c_str());
    if (isoParentZEntry_)    isoParentZEntry_->SetNumber(Z);
    if (isoParentNEntry_)    isoParentNEntry_->SetNumber(N);

    labelClassMap_[isoParentIsotope_] = "Parent";

    // Add parent to chain if not already there
    bool alreadyInChain = false;
    for (const auto& id : nucChainIsotopes_)
        if (id == isoParentIsotope_) { alreadyInChain = true; break; }
    if (!alreadyInChain) {
        nucChainIsotopes_.insert(nucChainIsotopes_.begin(), isoParentIsotope_);
        nucChainList_->RemoveAll();
        for (int i = 0; i < (int)nucChainIsotopes_.size(); i++) {
            const std::string& id = nucChainIsotopes_[i];
            auto dbIt = nuclearDB_.find(id);
            std::string entry = id;
            if (dbIt != nuclearDB_.end() && !dbIt->second.hl_str.empty())
                entry += "  T½=" + dbIt->second.hl_str;
            entry += "  [" + (labelClassMap_.count(id) ? labelClassMap_.at(id) : "?") + "]";
            nucChainList_->AddEntry(entry.c_str(), i + 1);
        }
        nucChainList_->Layout();
    }

    AppendLog("Parent set: " + isoParentIsotope_ +
              "  Z=" + std::to_string(Z) + "  N=" + std::to_string(N));
    if (nucStatusLbl_) nucStatusLbl_->SetText(("Parent: " + isoParentIsotope_).c_str());
    RefreshIsoComboHelper(nucIsoCombo_,      nucChainIsotopes_, nuclearDB_);
    RefreshIsoComboHelper(nucIsoLogftCombo_, nucChainIsotopes_, nuclearDB_);
    if (!isoHistName_.empty()) OnIsoSetParent();
}

// ─────────────────────────────────────────────────────────────────────────────
// OnNucAutoTraceChain
// BFS over all decay branches (β⁻, β⁻n, β⁻2n).  Each branch is traced until
// stable, long-lived, Q < 0, or not a β⁻ emitter.
// Generation naming: absStep encodes total steps from the Parent regardless of
// path.  prefix "" = main β⁻ chain; "Beta-n" / "Beta-2n" = neutron-emission
// sub-chains.
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnNucAutoTraceChain()
{
    if (isoParentZval_ <= 0 || isoParentNval_ <= 0) {
        OnNucSetParentFromChain();
        if (isoParentZval_ <= 0) { AppendLog("Nuclear: set parent first"); return; }
    }

    AppendLog("Nuclear: tracing decay chain from " + isoParentIsotope_ + " ...");
    if (nucStatusLbl_) nucStatusLbl_->SetText("Tracing chain...");

    const double kLongLivedThreshold = 3600.0;
    const int    kMaxAbsStep         = 12;
    const int    kMaxChainSize       = 50;
    const double kNeutronME          = 8071.3; // neutron mass excess (keV)

    // Generation name: 0=Parent, 1=Daughter, 2=Granddaughter, 3+=Great(s)-Granddaughter
    auto generationName = [](int step) -> std::string {
        if (step == 0) return "Parent";
        if (step == 1) return "Daughter";
        if (step == 2) return "Granddaughter";
        std::string s;
        for (int i = 3; i <= step; i++) s += "Great-";
        return s + "Granddaughter";
    };

    // Full class: prefix + " " + generation(absStep), or just generation for main chain.
    // Beta-n/2n products at absStep S encode "one step below" the emitter, matching
    // the same absolute distance from Parent as the main-chain nucleus at that step.
    auto makeClass = [&generationName](const std::string& prefix, int absStep) -> std::string {
        std::string gen = generationName(absStep);
        return prefix.empty() ? gen : (prefix + " " + gen);
    };

    // BFS queue entry
    struct QEntry { int Z, A, absStep; std::string prefix; };
    std::deque<QEntry>          queue;
    std::set<std::pair<int,int>> visited;

    queue.push_back({isoParentZval_, isoParentZval_ + isoParentNval_, 0, ""});

    while (!queue.empty() && (int)nucChainIsotopes_.size() < kMaxChainSize) {
        auto [Z, A, absStep, prefix] = queue.front();
        queue.pop_front();

        if (absStep > kMaxAbsStep) continue;
        if (visited.count(std::make_pair(Z, A))) continue;
        visited.insert(std::make_pair(Z, A));

        std::string sym = NucZToSymbol(Z);
        if (sym == "?" || sym.empty()) continue;
        std::string isoID = std::to_string(A) + sym;

        std::string cls = makeClass(prefix, absStep);
        labelClassMap_[isoID] = cls;

        bool found = false;
        for (const auto& id : nucChainIsotopes_) if (id == isoID) { found = true; break; }
        if (!found) nucChainIsotopes_.push_back(isoID);

        AppendLog("  Chain[" + std::to_string(absStep) +
                  (prefix.empty() ? "" : "/" + prefix) + "]: " + isoID + "  [" + cls + "]");

        // Fetch NNDC data if not cached
        auto dbIt = nuclearDB_.find(isoID);
        if (dbIt == nuclearDB_.end() || !dbIt->second.valid()) {
            NucIsotope iso;
            std::string cachePath = nucCacheDir_ + "/" + isoID + ".nucdat";
            if (!NNDCFetcher::LoadCache(cachePath, iso)) {
                if (nucStatusLbl_) nucStatusLbl_->SetText(("Fetching " + isoID + "...").c_str());
                gSystem->ProcessEvents();
                NNDCFetcher::Fetch(A, sym, iso, nucCacheDir_);
            }
            if (iso.valid()) nuclearDB_[isoID] = iso;
            dbIt = nuclearDB_.find(isoID);
        }

        bool hasData = (dbIt != nuclearDB_.end() && dbIt->second.valid());

        // Stop if stable or long-lived
        if (hasData) {
            const NucIsotope& iso = dbIt->second;
            if (iso.halflife_s == 0.0 &&
                (iso.hl_str == "stable" || iso.hl_str == "STABLE")) {
                AppendLog("  " + isoID + " is stable — branch ends.");
                continue;
            }
            if (absStep > 0 && iso.halflife_s > kLongLivedThreshold && iso.halflife_s > 0.0) {
                AppendLog("  " + isoID + "  T½=" + iso.hl_str + " — long-lived, stopping branch.");
                continue;
            }
        }

        // Determine if β⁻ emitter (also covers B-N/B-2N which are β⁻ + neutron emission)
        bool isBetaMinus = false;
        if (hasData) {
            const NucIsotope& iso = dbIt->second;
            auto isoBetaMode = [](const std::string& m) {
                return m == "B-" || m == "B-N" || m == "B-2N";
            };
            isBetaMinus = !iso.betaBranches.empty()
                          || isoBetaMode(iso.decayMode1)
                          || isoBetaMode(iso.decayMode2)
                          || isoBetaMode(iso.decayMode3);
        } else {
            isBetaMinus = ((A - Z) > Z); // neutron-rich heuristic
        }

        if (!isBetaMinus) {
            AppendLog("  " + isoID + " — not β⁻ emitter, stopping branch.");
            continue;
        }

        // ── β⁻ daughter: Q = ME(Z,A) - ME(Z+1,A) ────────────────────────────
        bool qBetaOk = true;
        if (ameLoaded_) {
            auto itP = ameTable_.find(std::make_pair(Z,   A));
            auto itD = ameTable_.find(std::make_pair(Z+1, A));
            if (itP != ameTable_.end() && itD != ameTable_.end()) {
                double Q = itP->second - itD->second;
                if (Q < 0.0) {
                    AppendLog("  " + isoID + "  Q_β⁻ = " + Fmt(Q, 1) + " keV < 0 — forbidden.");
                    qBetaOk = false;
                }
            }
        }
        if (qBetaOk && !visited.count(std::make_pair(Z+1, A)))
            queue.push_back({Z+1, A, absStep+1, prefix});

        // ── β⁻n / β⁻2n branching ratios ─────────────────────────────────────
        // Priority: NUBASE table (most accurate) → all three NNDC decay mode slots
        double brBetaN  = 0.0;
        double brBeta2N = 0.0;

        if (nubaseLoaded_) {
            auto nbit = nubaseTable_.find(std::make_pair(Z, A));
            if (nbit != nubaseTable_.end()) {
                if (nbit->second.brBetaN  > 0.0) brBetaN  = nbit->second.brBetaN;
                if (nbit->second.brBeta2N > 0.0) brBeta2N = nbit->second.brBeta2N;
            }
        }
        // NNDC fallback: scan all three decay mode slots
        if (hasData) {
            const NucIsotope& iso = dbIt->second;
            struct { const std::string& mode; double frac; } modes[] = {
                {iso.decayMode1, iso.decay1frac},
                {iso.decayMode2, iso.decay2frac},
                {iso.decayMode3, iso.decay3frac},
            };
            for (auto& m : modes) {
                if (brBetaN  < 1.0 && m.mode == "B-N"  && m.frac > 0.0)
                    brBetaN  = m.frac;
                if (brBeta2N < 1.0 && m.mode == "B-2N" && m.frac > 0.0)
                    brBeta2N = m.frac;
            }
            // If mode is present but fraction is zero/unknown, treat as significant
            if (brBetaN < 1.0) {
                for (auto& m : modes)
                    if (m.mode == "B-N")  { brBetaN  = 100.0; break; }
            }
            if (brBeta2N < 1.0) {
                for (auto& m : modes)
                    if (m.mode == "B-2N") { brBeta2N = 100.0; break; }
            }
        }

        // ── β⁻n daughter: (Z+1, A-1) ─────────────────────────────────────────
        // Q_βn = ME(Z,A) - ME(Z+1,A-1) - Δ_n
        if (brBetaN > 1.0) {
            bool qOk = true;
            if (ameLoaded_) {
                auto itP = ameTable_.find(std::make_pair(Z,   A  ));
                auto itD = ameTable_.find(std::make_pair(Z+1, A-1));
                if (itP != ameTable_.end() && itD != ameTable_.end()) {
                    double Q = itP->second - itD->second - kNeutronME;
                    if (Q < 0.0) {
                        AppendLog("  " + isoID + "  Q_β⁻n = " + Fmt(Q, 1) + " keV < 0 — forbidden.");
                        qOk = false;
                    }
                }
            }
            if (qOk && !visited.count(std::make_pair(Z+1, A-1)))
                queue.push_back({Z+1, A-1, absStep+1, "Beta-n"});
        }

        // ── β⁻2n daughter: (Z+1, A-2) ────────────────────────────────────────
        // Q_β2n = ME(Z,A) - ME(Z+1,A-2) - 2·Δ_n
        if (brBeta2N > 1.0) {
            bool qOk = true;
            if (ameLoaded_) {
                auto itP = ameTable_.find(std::make_pair(Z,   A  ));
                auto itD = ameTable_.find(std::make_pair(Z+1, A-2));
                if (itP != ameTable_.end() && itD != ameTable_.end()) {
                    double Q = itP->second - itD->second - 2.0 * kNeutronME;
                    if (Q < 0.0) {
                        AppendLog("  " + isoID + "  Q_β⁻2n = " + Fmt(Q, 1) + " keV < 0 — forbidden.");
                        qOk = false;
                    }
                }
            }
            if (qOk && !visited.count(std::make_pair(Z+1, A-2)))
                queue.push_back({Z+1, A-2, absStep+1, "Beta-2n"});
        }
    }

    // Rebuild chain list UI
    nucChainList_->RemoveAll();
    for (int i = 0; i < (int)nucChainIsotopes_.size(); i++) {
        const std::string& id = nucChainIsotopes_[i];
        auto dbIt2 = nuclearDB_.find(id);
        std::string entry = id;
        if (dbIt2 != nuclearDB_.end() && !dbIt2->second.hl_str.empty())
            entry += "  T½=" + dbIt2->second.hl_str;
        entry += "  [" + (labelClassMap_.count(id) ? labelClassMap_.at(id) : "?") + "]";
        nucChainList_->AddEntry(entry.c_str(), i + 1);
    }
    nucChainList_->Layout();

    PopulateNucGammaRef();
    RefreshIsoComboHelper(nucIsoCombo_,      nucChainIsotopes_, nuclearDB_);
    RefreshIsoComboHelper(nucIsoLogftCombo_, nucChainIsotopes_, nuclearDB_);

    std::string statusMsg = "Chain: " + std::to_string(nucChainIsotopes_.size()) + " isotopes";
    if (nucStatusLbl_) nucStatusLbl_->SetText(statusMsg.c_str());
    AppendLog("Chain trace complete: " + std::to_string(nucChainIsotopes_.size()) + " isotopes.");

    if (!isoHistName_.empty() && isoParentZval_ > 0) OnIsoSetParent();
}

// ─────────────────────────────────────────────────────────────────────────────
// OnNucFetchAll — fetch NNDC data for all chain isotopes
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnNucFetchAll()
{
    if (nucChainIsotopes_.empty()) {
        AppendLog("Nuclear: chain is empty — set parent and add isotopes first.");
        return;
    }
    NNDCFetcher::EnsureDir(nucCacheDir_);

    for (const auto& isoID : nucChainIsotopes_) {
        // Parse A and symbol from isoID
        size_t pos = 0;
        while (pos < isoID.size() && std::isdigit((unsigned char)isoID[pos])) pos++;
        if (pos == 0 || pos >= isoID.size()) continue;
        int A = std::stoi(isoID.substr(0, pos));
        std::string sym = isoID.substr(pos);

        if (nucStatusLbl_) nucStatusLbl_->SetText(("Fetching " + isoID + "...").c_str());
        gSystem->ProcessEvents();

        NucIsotope iso;
        bool ok = NNDCFetcher::Fetch(A, sym, iso, nucCacheDir_);
        if (ok) {
            nuclearDB_[isoID] = iso;
            AppendLog("Fetched: " + isoID + "  " + iso.hl_str +
                      "  (" + std::to_string(iso.gammas.size()) + " gammas)");
        } else {
            AppendLog("Fetch failed: " + isoID);
        }
    }

    // Rebuild chain list with half-life info
    nucChainList_->RemoveAll();
    for (int i = 0; i < (int)nucChainIsotopes_.size(); i++) {
        const std::string& id = nucChainIsotopes_[i];
        auto dbIt = nuclearDB_.find(id);
        std::string entry = id;
        if (dbIt != nuclearDB_.end() && !dbIt->second.hl_str.empty())
            entry += "  T½=" + dbIt->second.hl_str;
        entry += "  [" + (labelClassMap_.count(id) ? labelClassMap_.at(id) : "?") + "]";
        nucChainList_->AddEntry(entry.c_str(), i + 1);
    }
    nucChainList_->Layout();

    PopulateNucGammaRef();
    RefreshIsoComboHelper(nucIsoCombo_,      nucChainIsotopes_, nuclearDB_);
    RefreshIsoComboHelper(nucIsoLogftCombo_, nucChainIsotopes_, nuclearDB_);
    if (nucStatusLbl_) nucStatusLbl_->SetText("Fetch complete.");
}

// ─────────────────────────────────────────────────────────────────────────────
// OnNucReloadCache — reload nuclear data from local cache files
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnNucReloadCache()
{
    int n = 0;
    for (const auto& isoID : nucChainIsotopes_) {
        std::string cachePath = nucCacheDir_ + "/" + isoID + ".nucdat";
        NucIsotope iso;
        if (NNDCFetcher::LoadCache(cachePath, iso)) {
            nuclearDB_[isoID] = iso;
            ++n;
        }
    }
    // Rebuild chain list
    nucChainList_->RemoveAll();
    for (int i = 0; i < (int)nucChainIsotopes_.size(); i++) {
        const std::string& id = nucChainIsotopes_[i];
        auto dbIt = nuclearDB_.find(id);
        std::string entry = id;
        if (dbIt != nuclearDB_.end() && !dbIt->second.hl_str.empty())
            entry += "  T½=" + dbIt->second.hl_str;
        entry += "  [" + (labelClassMap_.count(id) ? labelClassMap_.at(id) : "?") + "]";
        nucChainList_->AddEntry(entry.c_str(), i + 1);
    }
    nucChainList_->Layout();
    PopulateNucGammaRef();
    RefreshIsoComboHelper(nucIsoCombo_,      nucChainIsotopes_, nuclearDB_);
    RefreshIsoComboHelper(nucIsoLogftCombo_, nucChainIsotopes_, nuclearDB_);
    AppendLog("Nuclear: reloaded " + std::to_string(n) + " isotopes from cache.");
    if (nucStatusLbl_) nucStatusLbl_->SetText(("Cache reloaded: " + std::to_string(n)).c_str());
}

// ─────────────────────────────────────────────────────────────────────────────
// OnNucAddToChain — manually add an isotope to the chain
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnNucAddToChain()
{
    if (!nucAddAEntry_ || !nucAddSymEntry_) return;
    int A = (int)nucAddAEntry_->GetNumber();
    std::string sym = nucAddSymEntry_->GetText();
    while (!sym.empty() && sym.front() == ' ') sym = sym.substr(1);
    while (!sym.empty() && sym.back()  == ' ') sym.pop_back();
    if (sym.empty() || A <= 0) { AppendLog("Nuclear: invalid add entry"); return; }
    if (sym.size() >= 1) sym[0] = (char)toupper((unsigned char)sym[0]);
    if (sym.size() >= 2) sym[1] = (char)tolower((unsigned char)sym[1]);

    std::string isoID = std::to_string(A) + sym;

    bool found = false;
    for (const auto& id : nucChainIsotopes_) if (id == isoID) { found = true; break; }
    if (found) { AppendLog("Nuclear: " + isoID + " already in chain"); return; }

    // Determine class from type combo
    std::string cls;
    if (nucAddTypeCombo_) {
        TGLBEntry* typeEntry = nucAddTypeCombo_->GetSelectedEntry();
        std::string typeStr = typeEntry ? typeEntry->GetTitle() : "(auto)";
        if (typeStr == "Background") cls = "Background";
        else if (typeStr == "Decay chain") {
            // auto-classify: if parent is set, make it Daughter
            cls = isoParentIsotope_.empty() ? "Decay chain" : "Daughter";
        }
    }
    if (!cls.empty()) labelClassMap_[isoID] = cls;

    nucChainIsotopes_.push_back(isoID);
    std::string entry = isoID;
    if (!cls.empty()) entry += "  [" + cls + "]";
    nucChainList_->AddEntry(entry.c_str(), (int)nucChainList_->GetNumberOfEntries() + 1);
    nucChainList_->Layout();
    AppendLog("Nuclear: added " + isoID + (cls.empty() ? "" : "  [" + cls + "]"));
    RefreshIsoComboHelper(nucIsoCombo_,      nucChainIsotopes_, nuclearDB_);
    RefreshIsoComboHelper(nucIsoLogftCombo_, nucChainIsotopes_, nuclearDB_);
}

// ─────────────────────────────────────────────────────────────────────────────
// OnNucRemoveFromChain
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnNucRemoveFromChain()
{
    Int_t sel = nucChainList_ ? nucChainList_->GetSelected() : -1;
    if (sel < 1 || (size_t)(sel - 1) >= nucChainIsotopes_.size()) {
        AppendLog("Nuclear: select an isotope from the chain list first."); return;
    }
    std::string removed = nucChainIsotopes_[sel - 1];
    nucChainIsotopes_.erase(nucChainIsotopes_.begin() + (sel - 1));
    nucChainList_->RemoveAll();
    for (int i = 0; i < (int)nucChainIsotopes_.size(); i++) {
        const std::string& id = nucChainIsotopes_[i];
        auto dbIt = nuclearDB_.find(id);
        std::string entry = id;
        if (dbIt != nuclearDB_.end() && !dbIt->second.hl_str.empty())
            entry += "  T½=" + dbIt->second.hl_str;
        entry += "  [" + (labelClassMap_.count(id) ? labelClassMap_.at(id) : "?") + "]";
        nucChainList_->AddEntry(entry.c_str(), i + 1);
    }
    nucChainList_->Layout();
    AppendLog("Nuclear: removed " + removed);
    RefreshIsoComboHelper(nucIsoCombo_,      nucChainIsotopes_, nuclearDB_);
    RefreshIsoComboHelper(nucIsoLogftCombo_, nucChainIsotopes_, nuclearDB_);
}

// ─────────────────────────────────────────────────────────────────────────────
// OnNucClearChain
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnNucClearChain()
{
    nucChainIsotopes_.clear();
    if (nucChainList_) { nucChainList_->RemoveAll(); nucChainList_->Layout(); }
    AppendLog("Nuclear: chain cleared.");
    RefreshIsoComboHelper(nucIsoCombo_,      nucChainIsotopes_, nuclearDB_);
    RefreshIsoComboHelper(nucIsoLogftCombo_, nucChainIsotopes_, nuclearDB_);
    if (nucStatusLbl_) nucStatusLbl_->SetText("Chain cleared.");
}

// ─────────────────────────────────────────────────────────────────────────────
// OnNucAddBackground — add a background isotope from the combo
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnNucAddBackground()
{
    if (!nucBgCombo_) return;
    TGLBEntry* sel = nucBgCombo_->GetSelectedEntry();
    if (!sel) return;
    std::string bgName = sel->GetTitle();

    // Parse A and symbol from bgName like "40K", "208Tl"
    size_t pos = 0;
    while (pos < bgName.size() && std::isdigit((unsigned char)bgName[pos])) pos++;
    if (pos == 0 || pos >= bgName.size()) return;
    int A = std::stoi(bgName.substr(0, pos));
    std::string sym = bgName.substr(pos);
    if (sym.empty() || A <= 0) return;

    std::string isoID = bgName;
    bool found = false;
    for (const auto& id : nucChainIsotopes_) if (id == isoID) { found = true; break; }
    if (found) { AppendLog("Nuclear: " + isoID + " already in chain"); return; }

    nucChainIsotopes_.push_back(isoID);
    labelClassMap_[isoID] = "Background";

    std::string entry = isoID + "  [Background]";
    nucChainList_->AddEntry(entry.c_str(), (int)nucChainList_->GetNumberOfEntries() + 1);
    nucChainList_->Layout();
    AppendLog("Nuclear: added background " + isoID);
    RefreshIsoComboHelper(nucIsoCombo_,      nucChainIsotopes_, nuclearDB_);
    RefreshIsoComboHelper(nucIsoLogftCombo_, nucChainIsotopes_, nuclearDB_);
}

// ─────────────────────────────────────────────────────────────────────────────
// OnNucLoadNNDCTxt — load a locally-saved IAEA CSV file
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnNucLoadNNDCTxt()
{
    TGFileInfo fi;
    static const char* kTypes[] = {
        "IAEA CSV files", "*.csv *.txt *.dat",
        "All files",      "*",
        nullptr,          nullptr
    };
    fi.fFileTypes = kTypes;
    fi.fIniDir    = StrDup(".");
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if (!fi.fFilename) return;

    std::string path = fi.fFilename;

    if (!nucAEntry_ || !nucSymbolEntry_) return;
    int A = (int)nucAEntry_->GetNumber();
    std::string sym = nucSymbolEntry_->GetText();
    while (!sym.empty() && sym.front() == ' ') sym = sym.substr(1);
    while (!sym.empty() && sym.back()  == ' ') sym.pop_back();
    if (sym.empty() || A <= 0) {
        AppendLog("Nuclear: set A and Sym first, then load txt file");
        return;
    }
    if (sym.size() >= 1) sym[0] = (char)toupper((unsigned char)sym[0]);
    if (sym.size() >= 2) sym[1] = (char)tolower((unsigned char)sym[1]);

    NucIsotope iso;
    bool ok = NNDCFetcher::LoadFromFile(path, A, sym, iso, nucCacheDir_);
    if (!ok) {
        AppendLog("Nuclear: failed to parse NNDC file: " + path);
        if (nucStatusLbl_) nucStatusLbl_->SetText("Load failed");
        return;
    }

    std::string isoID = std::to_string(A) + sym;
    auto it = nuclearDB_.find(isoID);
    if (it != nuclearDB_.end()) {
        NucIsotope& existing = it->second;
        if (!iso.gammas.empty())       existing.gammas       = iso.gammas;
        if (!iso.levels.empty())       existing.levels       = iso.levels;
        if (!iso.betaBranches.empty()) existing.betaBranches = iso.betaBranches;
        if (!iso.jpi.empty())          existing.jpi          = iso.jpi;
        if (iso.halflife_s > 0)        existing.halflife_s   = iso.halflife_s;
        if (!iso.hl_str.empty())       existing.hl_str       = iso.hl_str;
        existing._valid = true;
    } else {
        nuclearDB_[isoID] = iso;
        bool found = false;
        for (const auto& id : nucChainIsotopes_) if (id == isoID) { found = true; break; }
        if (!found) {
            nucChainIsotopes_.push_back(isoID);
            std::string entry = isoID + "  [loaded from file]";
            nucChainList_->AddEntry(entry.c_str(),
                (int)nucChainList_->GetNumberOfEntries() + 1);
            nucChainList_->Layout();
        }
    }

    PopulateNucGammaRef();
    RefreshIsoComboHelper(nucIsoCombo_,      nucChainIsotopes_, nuclearDB_);
    RefreshIsoComboHelper(nucIsoLogftCombo_, nucChainIsotopes_, nuclearDB_);
    AppendLog("Nuclear: loaded from file " + path + " → " + isoID +
              "  (" + std::to_string(iso.gammas.size()) + " gammas, " +
              std::to_string(iso.levels.size()) + " levels)");
    if (nucStatusLbl_) nucStatusLbl_->SetText(("Loaded: " + isoID).c_str());
}

// ─────────────────────────────────────────────────────────────────────────────
// OnNucDrawLevelScheme — draw on the embedded ROOT canvas
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnNucDrawLevelScheme()
{
    if (!nucIsoCombo_) return;
    TGLBEntry* e = nucIsoCombo_->GetSelectedEntry();
    if (!e) { AppendLog("[Nuclear] No isotope selected"); return; }
    std::string isoID = e->GetTitle();
    if (isoID == "(none)" || isoID.empty()) {
        AppendLog("[Nuclear] No isotope selected"); return;
    }
    if (nuclearDB_.find(isoID) == nuclearDB_.end()) {
        AppendLog("[Nuclear] No data for " + isoID + " — fetch first");
        return;
    }
    DrawLevelScheme(isoID);
}

// ─────────────────────────────────────────────────────────────────────────────
// OnNucOpenInteractive — generate Plotly HTML and open in browser
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnNucOpenInteractive()
{
    if (!nucIsoCombo_) return;
    TGLBEntry* e = nucIsoCombo_->GetSelectedEntry();
    if (!e) { AppendLog("[Nuclear] No isotope selected"); return; }
    std::string isoID = e->GetTitle();
    if (isoID == "(none)" || isoID.empty()) {
        AppendLog("[Nuclear] No isotope selected"); return;
    }
    if (nuclearDB_.find(isoID) == nuclearDB_.end()) {
        AppendLog("[Nuclear] No data for " + isoID + " — fetch first");
        return;
    }

    std::string html = GenerateLevelSchemePlotlyHTML(isoID);
    if (html.empty()) {
        AppendLog("[Nuclear] Could not generate level scheme for " + isoID);
        return;
    }

    std::string outPath = "/tmp/agf_level_scheme_" + isoID + ".html";
    std::ofstream f(outPath);
    if (!f.is_open()) {
        AppendLog("[Nuclear] Cannot write to " + outPath);
        return;
    }
    f << html;
    f.close();

    std::string cmd = "xdg-open \"" + outPath + "\" &";
    int ret = system(cmd.c_str());
    if (ret != 0)
        AppendLog("[Nuclear] Could not open browser — file saved to " + outPath);
    else
        AppendLog("[Nuclear] Opened interactive level scheme: " + outPath);
}

// ─────────────────────────────────────────────────────────────────────────────
// GenerateLevelSchemePlotlyHTML — build a self-contained Plotly HTML string
// ─────────────────────────────────────────────────────────────────────────────
std::string GammaFitGUI::GenerateLevelSchemePlotlyHTML(const std::string& isoID) const
{
    auto it = nuclearDB_.find(isoID);
    if (it == nuclearDB_.end()) return "";
    const NucIsotope& iso = it->second;
    if (iso.levels.empty() && iso.gammas.empty()) return "";

    // Sort levels by energy
    std::vector<NucLevel> levels = iso.levels;
    std::sort(levels.begin(), levels.end(),
        [](const NucLevel& a, const NucLevel& b){ return a.energy < b.energy; });

    // Ensure ground state is present
    if (levels.empty() || levels[0].energy > 1.0) {
        NucLevel gs; gs.energy = 0.0; gs.jpi = iso.jpi;
        gs.hl_str = iso.hl_str; gs.halflife_s = iso.halflife_s;
        levels.insert(levels.begin(), gs);
    }

    // Cap at 80 levels for readability
    if (levels.size() > 80) levels.resize(80);

    // Max energy for axis
    double maxE = 200.0;
    for (auto& lv : levels) maxE = std::max(maxE, lv.energy);
    maxE *= 1.12;

    // Helper: find closest level energy within 20 keV tolerance
    auto closestLevelE = [&](double E) -> double {
        double best = E; double bestD = 20.0;
        for (auto& lv : levels) {
            double d = std::abs(lv.energy - E);
            if (d < bestD) { bestD = d; best = lv.energy; }
        }
        return best;
    };

    // Escape string for JSON
    auto jsStr = [](const std::string& s) -> std::string {
        std::string out;
        for (char c : s) {
            if      (c == '"')  out += "\\\"";
            else if (c == '\\') out += "\\\\";
            else if (c == '\n') out += "<br>";
            else                out += c;
        }
        return out;
    };

    std::ostringstream shapes, annots, hx, hy, htxt;
    shapes << "[";
    annots << "[";
    bool firstS = true, firstA = true, firstH = true;

    // ── Level lines ──────────────────────────────────────────────────────────
    for (auto& lv : levels) {
        double y = lv.energy;

        // Line shape
        if (!firstS) shapes << ",";
        shapes << "{\"type\":\"line\","
               << "\"x0\":0.04,\"x1\":0.72,"
               << "\"y0\":" << y << ",\"y1\":" << y << ","
               << "\"line\":{\"color\":\"#7788cc\",\"width\":2}}";
        firstS = false;

        // Energy label (left)
        std::string eStr = (y == 0.0) ? "0 g.s." : Form("%.1f", y);
        if (!firstA) annots << ",";
        annots << "{\"x\":0.03,\"y\":" << y << ","
               << "\"text\":\"" << jsStr(eStr) << "\","
               << "\"showarrow\":false,\"xanchor\":\"right\","
               << "\"font\":{\"color\":\"#99aade\",\"size\":9}}";
        firstA = false;

        // Jpi label (right)
        std::string jpi = lv.jpi.empty() ? "?" : lv.jpi;
        std::string hl  = lv.hl_str.empty() ? (y == 0.0 ? iso.hl_str : "?") : lv.hl_str;
        annots << ",{\"x\":0.74,\"y\":" << y << ","
               << "\"text\":\"" << jsStr(jpi) << "\","
               << "\"showarrow\":false,\"xanchor\":\"left\","
               << "\"font\":{\"color\":\"#ccd4ff\",\"size\":9}}";

        // Hover scatter point
        if (!firstH) { hx << ","; hy << ","; htxt << ","; }
        hx   << "0.38";
        hy   << y;
        htxt << "\"E=" << Form("%.3f", y) << " keV"
             << "<br>Jπ=" << jsStr(jpi)
             << "<br>T½=" << jsStr(hl) << "\"";
        firstH = false;
    }

    // ── Gamma arrows ─────────────────────────────────────────────────────────
    // Group by start level so we can spread x positions
    std::map<double, std::vector<size_t>> byStart;
    for (size_t i = 0; i < iso.gammas.size(); i++) {
        double sl = closestLevelE(iso.gammas[i].start_level);
        byStart[sl].push_back(i);
    }

    // Limit total arrows to 120 (keep highest intensity)
    std::vector<const NucGamma*> sorted;
    for (auto& gm : iso.gammas) sorted.push_back(&gm);
    std::sort(sorted.begin(), sorted.end(),
        [](const NucGamma* a, const NucGamma* b){ return a->intensity > b->intensity; });
    std::set<size_t> topIdx;
    for (size_t i = 0; i < std::min((size_t)120, sorted.size()); i++) {
        topIdx.insert(sorted[i] - iso.gammas.data());
    }

    for (auto& kv : byStart) {
        const auto& indices = kv.second;
        int n = (int)indices.size();
        for (int j = 0; j < n; j++) {
            size_t gi = indices[j];
            if (topIdx.find(gi) == topIdx.end()) continue;
            const NucGamma& gm = iso.gammas[gi];

            double startE = closestLevelE(gm.start_level);
            double endE   = closestLevelE(startE - gm.energy);
            if (startE <= endE + 0.5) continue;

            double xPos  = 0.10 + (n > 1 ? (double)j / (n - 1) * 0.56 : 0.28);
            double width = std::max(0.8, std::min(5.0, gm.intensity / 20.0));
            double alpha = std::max(0.25, std::min(1.0, 0.25 + gm.intensity / 70.0));

            std::string color = Form("rgba(255,210,30,%.2f)", alpha);
            std::string label;
            if (gm.intensity > 2.0)
                label = Form("%.1f keV", gm.energy);

            annots << ",{\"x\":" << xPos << ",\"y\":" << endE << ","
                   << "\"ax\":" << xPos << ",\"ay\":" << startE << ","
                   << "\"axref\":\"x\",\"ayref\":\"y\","
                   << "\"arrowhead\":2,\"arrowsize\":1,"
                   << "\"arrowwidth\":" << Form("%.1f", width) << ","
                   << "\"arrowcolor\":\"" << color << "\","
                   << "\"text\":\"" << jsStr(label) << "\","
                   << "\"font\":{\"color\":\"#ffe060\",\"size\":8},"
                   << "\"showarrow\":true,"
                   << "\"hovertext\":\"Eγ=" << Form("%.3f", gm.energy)
                   << " keV<br>I=" << Form("%.1f", gm.intensity) << "%"
                   << "<br>from " << Form("%.1f", startE) << " keV\","
                   << "\"bgcolor\":\"rgba(20,20,40,0.85)\","
                   << "\"borderpad\":2}";
        }
    }
    shapes << "]";
    annots << "]";

    std::string title = Form("<sup>%d</sup>%s Level Scheme", iso.A, iso.symbol.c_str());
    std::string sub   = "Jπ=" + iso.jpi + "  T½=" + iso.hl_str
                      + Form("  (%d levels, %d γs shown)",
                             (int)levels.size(),
                             (int)std::min((size_t)120, iso.gammas.size()));

    std::ostringstream html;
    html << R"(<!DOCTYPE html>
<html><head>
<meta charset="utf-8">
<title>)" << isoID << R"( Level Scheme</title>
<script src="https://cdn.plot.ly/plotly-2.27.0.min.js" charset="utf-8"></script>
<style>body{margin:0;background:#13131f;font-family:sans-serif;}
#info{position:fixed;bottom:8px;left:8px;color:#556;font-size:11px;}
</style>
</head>
<body>
<div id="plot" style="width:100%;height:100vh;"></div>
<div id="info">AutoGammaFit — hover over levels and arrows for details</div>
<script>
var trace = {
  x:[)" << hx.str() << R"(],
  y:[)" << hy.str() << R"(],
  text:[)" << htxt.str() << R"(],
  mode:'markers', type:'scatter',
  marker:{size:18,opacity:0},
  hoverinfo:'text',
  hoverlabel:{bgcolor:'#1e1e3e',bordercolor:'#445',font:{color:'#ddeeff',size:12}},
  showlegend:false
};
var layout = {
  paper_bgcolor:'#13131f', plot_bgcolor:'#13131f',
  title:{text:')" << title << "<br><sup>" << jsStr(sub) << R"(</sup>',
         font:{color:'#c8d8ff',size:15},x:0.5},
  xaxis:{visible:false,range:[0,1],fixedrange:true},
  yaxis:{
    title:{text:'Energy (keV)',font:{color:'#8899bb',size:12},standoff:8},
    tickfont:{color:'#8899bb',size:10},
    gridcolor:'#222233', gridwidth:1,
    zeroline:true, zerolinecolor:'#3344aa', zerolinewidth:1.5,
    range:[-)" << maxE * 0.05 << "," << maxE << R"(]
  },
  shapes:)" << shapes.str() << R"(,
  annotations:)" << annots.str() << R"(,
  margin:{l:65,r:100,t:75,b:40},
  hovermode:'closest'
};
Plotly.newPlot('plot',[trace],layout,
  {responsive:true,displayModeBar:true,
   modeBarButtonsToRemove:['lasso2d','select2d']});
</script>
</body></html>)";

    return html.str();
}

// ─────────────────────────────────────────────────────────────────────────────
// OnNucConfirmChainToIsoDB
// Load all gamma lines for the current decay chain into the Isotope DB (db_).
// Existing entries are preserved; duplicates (same isotope + energy within
// 0.1 keV) are skipped so re-clicking is safe.
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnNucConfirmChainToIsoDB()
{
    if (nucChainIsotopes_.empty()) {
        AppendLog("Nuclear: chain is empty — trace a chain first.");
        return;
    }

    int nAdded = 0;
    for (const auto& isoID : nucChainIsotopes_) {
        auto it = nuclearDB_.find(isoID);
        if (it == nuclearDB_.end() || !it->second.valid()) continue;
        const NucIsotope& iso = it->second;
        if (iso.gammas.empty()) continue;

        for (const auto& gm : iso.gammas) {
            if (gm.energy <= 0.0) continue;
            // Skip if already in db_ for this isotope within 0.1 keV
            bool dup = false;
            for (const auto& line : db_.db) {
                if (line.isotope == isoID &&
                    std::abs(line.energy - gm.energy) < 0.1) { dup = true; break; }
            }
            if (dup) continue;
            GammaLine gl;
            gl.isotope      = isoID;
            gl.energy       = gm.energy;
            gl.intensity    = gm.intensity;
            gl.hasIntensity = (gm.intensity > 0.0);
            db_.db.push_back(gl);
            ++nAdded;
        }
    }

    dbLoaded_ = !db_.db.empty();
    if (isotopeLbl_)
        isotopeLbl_->SetText(
            (std::string("Chain DB  (") + std::to_string(db_.db.size()) + " lines)").c_str());
    PopulateIsoDbList();
    AppendLog("Confirm Chain: added " + std::to_string(nAdded) +
              " lines from " + std::to_string(nucChainIsotopes_.size()) +
              " isotopes to Isotope DB.  Total: " + std::to_string(db_.db.size()) + " lines.");
    SetStatus("Isotope DB updated from chain.");
}

// ─────────────────────────────────────────────────────────────────────────────
// PopulateNucRefTab — fill nucRefList_ with gammas classified by decay type,
// each row carrying ENSDF author/cutoff reference info.
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::PopulateNucRefTab()
{
    if (!nucRefList_) return;
    nucRefList_->RemoveAll();

    // Which decay class filter is active?
    std::string filterClass;
    if (nucRefDecayFilter_) {
        TGLBEntry* fe = nucRefDecayFilter_->GetSelectedEntry();
        if (fe) {
            filterClass = fe->GetTitle();
            if (filterClass == "All") filterClass.clear();
        }
    }

    int rowId = 1;
    for (const auto& isoID : nucChainIsotopes_) {
        auto it = nuclearDB_.find(isoID);
        if (it == nuclearDB_.end() || !it->second.valid()) continue;
        const NucIsotope& iso = it->second;

        std::string cls = labelClassMap_.count(isoID) ? labelClassMap_.at(isoID) : "Other";
        if (!filterClass.empty() && cls != filterClass) continue;

        std::string decayMode = iso.decayMode1.empty() ? "?" : iso.decayMode1;

        // Section header row (not selectable — we mark it with negative id)
        std::string hdr = "── " + isoID + "  [" + cls + "]  decay=" + decayMode
                        + "  T½=" + iso.hl_str + " ──";
        nucRefList_->AddEntry(hdr.c_str(), -rowId);
        ++rowId;

        // Gamma rows sorted by energy
        std::vector<NucGamma> sorted = iso.gammas;
        std::sort(sorted.begin(), sorted.end(),
            [](const NucGamma& a, const NucGamma& b){ return a.energy < b.energy; });

        for (const auto& gm : sorted) {
            if (gm.energy <= 0) continue;
            char buf[320];
            std::string ref = gm.ensdf_authors.empty() ? "NNDC"
                            : gm.ensdf_authors + "  (" + gm.ensdf_cutoff + ")";
            snprintf(buf, sizeof(buf),
                "  [%s] %s  %.2f keV  I=%.1f%%  %s  | %s",
                cls.c_str(), isoID.c_str(),
                gm.energy, gm.intensity,
                gm.multipolarity.empty() ? "--" : gm.multipolarity.c_str(),
                ref.c_str());
            nucRefList_->AddEntry(buf, rowId);
            ++rowId;
        }
        if (sorted.empty()) {
            nucRefList_->AddEntry("  (no gamma data)", rowId++);
        }
        // Blank separator
        nucRefList_->AddEntry("", rowId++);
    }
    nucRefList_->Layout();
    if (nucRefDetailLbl_)
        nucRefDetailLbl_->SetText(
            nucChainIsotopes_.empty() ? "Build a decay chain first (Isotope Chain tab)"
                                      : "Select a gamma line to see details");
}

// ─────────────────────────────────────────────────────────────────────────────
// OnNucRefFilterChanged — re-filter the reference list
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnNucRefFilterChanged(Int_t /*id*/)
{
    PopulateNucRefTab();

    // If a row was actually selected, update the detail label
    if (!nucRefList_) return;
    TGLBEntry* sel = nucRefList_->GetSelectedEntry();
    if (!sel || !nucRefDetailLbl_) return;
    std::string txt = sel->GetTitle();
    // Skip header/blank rows
    if (txt.empty() || txt.substr(0,2) == "──") return;
    nucRefDetailLbl_->SetText(txt.c_str());
}

// ─────────────────────────────────────────────────────────────────────────────
// OnNucRefOpen — open the NuDat page for the selected gamma's isotope
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnNucRefOpen()
{
    if (!nucRefList_) return;
    TGLBEntry* sel = nucRefList_->GetSelectedEntry();
    if (!sel) { AppendLog("References: select a gamma line first."); return; }

    std::string txt = sel->GetTitle();
    if (txt.empty() || txt.substr(0,2) == "──") {
        AppendLog("References: select a gamma line row (not a header)."); return;
    }

    // Extract isotope ID from the row — format: "  [Class] ISOID  E keV ..."
    // Find the isotope ID: first non-space token after ']'
    std::string isoID;
    auto bracket = txt.find(']');
    if (bracket != std::string::npos) {
        std::istringstream ss(txt.substr(bracket + 1));
        ss >> isoID;
    }
    if (isoID.empty()) { AppendLog("References: could not parse isotope from row."); return; }

    // Build NuDat URL: https://www.nndc.bnl.gov/nudat3/getdatasetClassic.jsp?nucleus=44S
    std::string url = "https://www.nndc.bnl.gov/nudat3/getdatasetClassic.jsp?nucleus=" + isoID;
    std::string cmd = "xdg-open '" + url + "' &";
    gSystem->Exec(cmd.c_str());
    AppendLog("References: opening " + url);
    if (nucRefDetailLbl_) nucRefDetailLbl_->SetText(("NNDC: " + url).c_str());
}

// ─────────────────────────────────────────────────────────────────────────────
// OnNucSetHistParent — assign the typed histogram to the typed parent nucleus
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnNucSetHistParent()
{
    // Find the two text entries by name in the References sub-tab
    TGTextEntry* histEntry   = nullptr;
    TGTextEntry* parentEntry = nullptr;
    // Walk the widget tree looking for our named entries
    auto findEntry = [&](TGCompositeFrame* frame, const char* name) -> TGTextEntry* {
        if (!frame) return nullptr;
        TList* list = frame->GetList();
        if (!list) return nullptr;
        TIter next(list);
        while (TObject* obj = next()) {
            if (TGTextEntry* te = dynamic_cast<TGTextEntry*>(obj)) {
                if (strcmp(te->GetName(), name) == 0) return te;
            }
        }
        return nullptr;
    };
    // Search in the main frame
    histEntry   = findEntry(this, "nucHistParentHistEntry");
    parentEntry = findEntry(this, "nucHistParentEntry");

    std::string histName, parentID;
    if (histEntry)   histName = histEntry->GetText();
    if (parentEntry) parentID = parentEntry->GetText();

    // Trim
    auto trim = [](std::string s) {
        while (!s.empty() && s.front() == ' ') s = s.substr(1);
        while (!s.empty() && s.back()  == ' ') s.pop_back();
        return s;
    };
    histName = trim(histName);
    parentID = trim(parentID);

    if (histName.empty() || parentID.empty()) {
        AppendLog("References: enter both histogram name and parent nucleus ID.");
        return;
    }

    histParent_[histName] = parentID;
    SaveChainCache();
    AppendLog("References: '" + histName + "' → parent '" + parentID + "'  (saved to chain cache)");
    SetStatus("Histogram parent set: " + histName + " → " + parentID);
}

// ─────────────────────────────────────────────────────────────────────────────
// SaveChainCache — persist decay chain + histParent_ to a .chaindat file
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::SaveChainCache()
{
    if (nucChainIsotopes_.empty()) return;
    std::string parentID = nucChainIsotopes_.empty() ? "unknown" : nucChainIsotopes_.front();
    std::string dir = nucCacheDir_.empty()
                    ? (launchDir_ + "/nuclear_cache")
                    : nucCacheDir_;
    // Ensure directory exists
    struct stat st{};
    if (stat(dir.c_str(), &st) != 0) mkdir(dir.c_str(), 0755);

    std::string path = dir + "/" + parentID + ".chaindat";
    std::ofstream f(path);
    if (!f.is_open()) { AppendLog("SaveChainCache: cannot write " + path); return; }

    f << "# AutoGammaFit chain cache — do not edit manually\n";
    f << "PARENT " << parentID << "\n";
    f << "CHAIN_SIZE " << nucChainIsotopes_.size() << "\n";
    for (size_t i = 0; i < nucChainIsotopes_.size(); i++)
        f << "CHAIN_" << i << " " << nucChainIsotopes_[i] << "\n";

    // Histogram → parent mappings
    for (const auto& kv : histParent_)
        f << "HIST_PARENT " << kv.first << "|" << kv.second << "\n";

    // Per-isotope gamma lines with references
    for (const auto& isoID : nucChainIsotopes_) {
        auto it = nuclearDB_.find(isoID);
        if (it == nuclearDB_.end() || !it->second.valid()) continue;
        f << "ISO_BEGIN " << isoID << "\n";
        for (const auto& gm : it->second.gammas) {
            if (gm.energy <= 0) continue;
            f << "REF_GAMMA " << gm.energy << "|" << gm.intensity << "|"
              << gm.start_level << "|" << gm.end_level << "|"
              << gm.multipolarity << "|" << gm.ensdf_authors << "|"
              << gm.ensdf_cutoff << "\n";
        }
        f << "ISO_END\n";
    }
    AppendLog("Chain cache saved: " + path);
}

// ─────────────────────────────────────────────────────────────────────────────
// LoadChainCache — restore chain + histParent_ from .chaindat file
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::LoadChainCache()
{
    std::string dir = nucCacheDir_.empty()
                    ? (launchDir_ + "/nuclear_cache")
                    : nucCacheDir_;

    // Find any .chaindat files and load the most recently modified one
    std::string bestFile;
    double bestTime = -1;
    DIR* d = opendir(dir.c_str());
    if (!d) return;
    struct dirent* ent;
    while ((ent = readdir(d)) != nullptr) {
        std::string name = ent->d_name;
        if (name.size() > 9 && name.substr(name.size()-9) == ".chaindat") {
            std::string fp = dir + "/" + name;
            struct stat st{};
            if (stat(fp.c_str(), &st) == 0 && (double)st.st_mtime > bestTime) {
                bestTime = (double)st.st_mtime;
                bestFile = fp;
            }
        }
    }
    closedir(d);
    if (bestFile.empty()) return;

    std::ifstream f(bestFile);
    if (!f.is_open()) return;

    nucChainIsotopes_.clear();
    histParent_.clear();

    std::string line, curISO;
    size_t chainSize = 0;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string tag; ss >> tag;

        if (tag == "CHAIN_SIZE") { ss >> chainSize; nucChainIsotopes_.resize(chainSize); }
        else if (tag.substr(0, 6) == "CHAIN_") {
            size_t idx = 0;
            try { idx = std::stoul(tag.substr(6)); } catch (...) {}
            std::string iso; ss >> iso;
            if (idx < nucChainIsotopes_.size()) nucChainIsotopes_[idx] = iso;
        }
        else if (tag == "HIST_PARENT") {
            std::string rec; std::getline(ss, rec);
            auto pipe = rec.find('|');
            if (pipe != std::string::npos) {
                std::string hist   = rec.substr(0, pipe);
                std::string parent = rec.substr(pipe + 1);
                auto trim = [](std::string s) {
                    while (!s.empty() && s.front()==' ') s=s.substr(1);
                    while (!s.empty() && s.back()==' ') s.pop_back();
                    return s; };
                hist = trim(hist); parent = trim(parent);
                if (!hist.empty() && !parent.empty())
                    histParent_[hist] = parent;
            }
        }
        else if (tag == "ISO_BEGIN") { ss >> curISO; }
        else if (tag == "REF_GAMMA" && !curISO.empty()) {
            std::string rec; std::getline(ss, rec);
            std::istringstream ps(rec);
            std::string field;
            std::vector<std::string> fields;
            while (std::getline(ps, field, '|')) fields.push_back(field);
            if (fields.size() >= 2) {
                NucGamma gm;
                try { gm.energy      = std::stod(fields[0]); } catch (...) {}
                try { gm.intensity   = std::stod(fields[1]); } catch (...) {}
                if (fields.size() > 2) try { gm.start_level = std::stod(fields[2]); } catch (...) {}
                if (fields.size() > 3) try { gm.end_level   = std::stod(fields[3]); } catch (...) {}
                if (fields.size() > 4) gm.multipolarity  = fields[4];
                if (fields.size() > 5) gm.ensdf_authors  = fields[5];
                if (fields.size() > 6) gm.ensdf_cutoff   = fields[6];
                if (gm.energy > 0) {
                    auto& iso = nuclearDB_[curISO];
                    iso._valid = true;
                    iso.symbol = curISO;
                    iso.gammas.push_back(gm);
                }
            }
        }
        else if (tag == "ISO_END") { curISO.clear(); }
    }

    if (!nucChainIsotopes_.empty())
        AppendLog("Chain cache loaded: " + bestFile + " ("
                  + std::to_string(nucChainIsotopes_.size()) + " isotopes)");
}
