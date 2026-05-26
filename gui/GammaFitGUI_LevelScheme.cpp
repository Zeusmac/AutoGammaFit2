#include "GammaFitGUI.h"
#include "GammaFitGUI_shared.h"
#include "LevelScheme.h"
#include "NNDCFetcher.h"

#include "TGTextEntry.h"
#include "TGFileDialog.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TLine.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TBox.h"
#include "TLegend.h"

#include <algorithm>
#include <cstdio>
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>

// ─────────────────────────────────────────────────────────────────────────────
// Helpers
// ─────────────────────────────────────────────────────────────────────────────

static std::string LSCachePath(const std::string& launchDir, const std::string& isoID)
{
    return launchDir + "/level_schemes/" + isoID + ".lsdat";
}

static void EnsureLSDir(const std::string& launchDir)
{
    std::string d = launchDir + "/level_schemes";
    ::mkdir(d.c_str(), 0755);
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildLSLevelSchemeTab  -  sub-tab 3 content
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::BuildLSLevelSchemeTab(TGCompositeFrame* p)
{
    // ── Top action row ────────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Isotope");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

        // Isotope combo (shared with OnNucDrawLevelScheme)
        TGHorizontalFrame* row1 = new TGHorizontalFrame(grp);
        grp->AddFrame(row1, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        row1->AddFrame(new TGLabel(row1, "Isotope:"),
                       new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        nucIsoCombo_ = new TGComboBox(row1, 1200);
        nucIsoCombo_->AddEntry("(none)", 1);
        nucIsoCombo_->Select(1, kFALSE);
        nucIsoCombo_->Resize(160, 22);
        row1->AddFrame(nucIsoCombo_, new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 0, 4, 0, 0));

        TGTextButton* seedBtn = new TGTextButton(row1, "Seed NNDC");
        row1->AddFrame(seedBtn, new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        seedBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLSSeedFromNNDC()");
        seedBtn->SetToolTipText("Populate the level scheme from NNDC nuclear data for the selected isotope");

        TGHorizontalFrame* row2 = new TGHorizontalFrame(grp);
        grp->AddFrame(row2, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));

        TGTextButton* loadBtn = new TGTextButton(row2, "Load .lsdat");
        row2->AddFrame(loadBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        loadBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLSLoad()");

        TGTextButton* saveBtn = new TGTextButton(row2, "Save .lsdat");
        row2->AddFrame(saveBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        saveBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLSSave()");

        TGTextButton* clearBtn = new TGTextButton(row2, "Clear");
        row2->AddFrame(clearBtn, new TGLayoutHints(kLHintsLeft));
        clearBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLSDrawEnhanced()");
        clearBtn->SetToolTipText("Draw the enhanced level scheme on the main canvas");
        // Re-use slot; "Clear" button is actually "Draw Enhanced" here
        // (label set by separate button below)
    }

    // ── Level editor ─────────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Nuclear Levels");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        lsLevelList_ = new TGListBox(grp, 3100);
        lsLevelList_->Resize(280, 120);
        grp->AddFrame(lsLevelList_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        lsLevelList_->AddEntry("  (load or seed from NNDC)", 1);
        lsLevelList_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                              "OnLSLevelSelected(Int_t)");

        // Add-level row
        TGHorizontalFrame* addRow = new TGHorizontalFrame(grp);
        grp->AddFrame(addRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        addRow->AddFrame(new TGLabel(addRow, "E(keV):"),
                         new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        lsLevelEnergyEntry_ = new TGNumberEntry(addRow, 0.0, 7, -1,
                                                TGNumberFormat::kNESRealFour,
                                                TGNumberFormat::kNEANonNegative);
        lsLevelEnergyEntry_->SetWidth(65);
        addRow->AddFrame(lsLevelEnergyEntry_, new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));

        addRow->AddFrame(new TGLabel(addRow, "Jpi:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        lsLevelJpiEntry_ = new TGTextEntry(addRow, "0+");
        lsLevelJpiEntry_->Resize(40, 22);
        addRow->AddFrame(lsLevelJpiEntry_, new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));

        addRow->AddFrame(new TGLabel(addRow, "beta%:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        lsBetaFeedEntry_ = new TGNumberEntry(addRow, 0.0, 5, -1,
                                             TGNumberFormat::kNESRealFour,
                                             TGNumberFormat::kNEANonNegative);
        lsBetaFeedEntry_->SetWidth(50);
        addRow->AddFrame(lsBetaFeedEntry_, new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));

        // Beta type entry
        addRow->AddFrame(new TGLabel(addRow, "type:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        lsBetaTypeEntry_ = new TGTextEntry(addRow, "GT");
        lsBetaTypeEntry_->Resize(38, 22);
        lsBetaTypeEntry_->SetToolTipText("GT, Fermi, Mixed, 1F, Forbidden, Unknown");
        addRow->AddFrame(lsBetaTypeEntry_, new TGLayoutHints(kLHintsCenterY));

        TGHorizontalFrame* btnRow = new TGHorizontalFrame(grp);
        grp->AddFrame(btnRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));

        TGTextButton* addLvlBtn = new TGTextButton(btnRow, " Add Level ");
        btnRow->AddFrame(addLvlBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        addLvlBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLSAddLevel()");

        TGTextButton* remLvlBtn = new TGTextButton(btnRow, " Remove ");
        btnRow->AddFrame(remLvlBtn, new TGLayoutHints(kLHintsLeft));
        remLvlBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLSRemoveLevel()");
    }

    // ── Transition editor ─────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Gamma Transitions");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        lsTransList_ = new TGListBox(grp, 3200);
        lsTransList_->Resize(280, 100);
        grp->AddFrame(lsTransList_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        lsTransList_->AddEntry("  (no transitions)", 1);
        lsTransList_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                              "OnLSTransSelected(Int_t)");

        // From / To combos
        TGHorizontalFrame* ftRow = new TGHorizontalFrame(grp);
        grp->AddFrame(ftRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        ftRow->AddFrame(new TGLabel(ftRow, "From:"),
                        new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        lsTransFromCombo_ = new TGComboBox(ftRow, 3210);
        lsTransFromCombo_->AddEntry("(none)", 1);
        lsTransFromCombo_->Select(1, kFALSE);
        lsTransFromCombo_->Resize(95, 22);
        ftRow->AddFrame(lsTransFromCombo_,
                        new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));

        ftRow->AddFrame(new TGLabel(ftRow, "To:"),
                        new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        lsTransToCombo_ = new TGComboBox(ftRow, 3211);
        lsTransToCombo_->AddEntry("(none)", 1);
        lsTransToCombo_->Select(1, kFALSE);
        lsTransToCombo_->Resize(95, 22);
        ftRow->AddFrame(lsTransToCombo_, new TGLayoutHints(kLHintsCenterY));

        // Energy / multipolarity / intensity row
        TGHorizontalFrame* emiRow = new TGHorizontalFrame(grp);
        grp->AddFrame(emiRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));

        emiRow->AddFrame(new TGLabel(emiRow, "E(keV):"),
                         new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        lsTransEnergyEntry_ = new TGNumberEntry(emiRow, 0.0, 7, -1,
                                                TGNumberFormat::kNESRealFour,
                                                TGNumberFormat::kNEANonNegative);
        lsTransEnergyEntry_->SetWidth(65);
        emiRow->AddFrame(lsTransEnergyEntry_, new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));

        emiRow->AddFrame(new TGLabel(emiRow, "Mult:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        lsTransMultEntry_ = new TGTextEntry(emiRow, "E2");
        lsTransMultEntry_->Resize(42, 22);
        lsTransMultEntry_->SetToolTipText("E1, E2, M1, M1+E2, etc.");
        emiRow->AddFrame(lsTransMultEntry_, new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));

        emiRow->AddFrame(new TGLabel(emiRow, "I:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        lsTransIntEntry_ = new TGNumberEntry(emiRow, 0.0, 6, -1,
                                             TGNumberFormat::kNESRealFour,
                                             TGNumberFormat::kNEANonNegative);
        lsTransIntEntry_->SetWidth(55);
        emiRow->AddFrame(lsTransIntEntry_, new TGLayoutHints(kLHintsCenterY));

        // Buttons
        TGHorizontalFrame* tBtnRow = new TGHorizontalFrame(grp);
        grp->AddFrame(tBtnRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));

        TGTextButton* addTrBtn = new TGTextButton(tBtnRow, " Add Transition ");
        tBtnRow->AddFrame(addTrBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        addTrBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLSAddTransition()");

        TGTextButton* remTrBtn = new TGTextButton(tBtnRow, " Remove ");
        tBtnRow->AddFrame(remTrBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        remTrBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLSRemoveTransition()");

        TGTextButton* linkBtn = new TGTextButton(tBtnRow, " Link Peak ");
        tBtnRow->AddFrame(linkBtn, new TGLayoutHints(kLHintsLeft));
        linkBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLSLinkPeak()");
        linkBtn->SetToolTipText(
            "Link the currently selected Peak Table row to the selected transition.\n"
            "Sets the intensity and cache link automatically.");
    }

    // ── Drawing ───────────────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Draw");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

        TGHorizontalFrame* dRow = new TGHorizontalFrame(grp);
        grp->AddFrame(dRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        TGTextButton* drawNNDCBtn = new TGTextButton(dRow, "NNDC Only");
        dRow->AddFrame(drawNNDCBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        drawNNDCBtn->Connect("Clicked()", "GammaFitGUI", this, "OnNucDrawLevelScheme()");
        drawNNDCBtn->SetToolTipText("Draw basic level scheme from NNDC data only");

        TGTextButton* drawEnhBtn = new TGTextButton(dRow, "Enhanced (user data)");
        dRow->AddFrame(drawEnhBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        drawEnhBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLSDrawEnhanced()");
        drawEnhBtn->SetToolTipText(
            "Draw the level scheme using user-edited level/transition data.\n"
            "Arrow thickness scales with intensity; beta-feeding shown on the left.");

        TGTextButton* interBtn = new TGTextButton(dRow, "Browser");
        dRow->AddFrame(interBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        interBtn->Connect("Clicked()", "GammaFitGUI", this, "OnNucOpenInteractive()");
        interBtn->SetToolTipText("Open interactive Plotly level scheme in browser");

        TGTextButton* schemBtn = new TGTextButton(dRow, "Decay Schematic");
        dRow->AddFrame(schemBtn, new TGLayoutHints(kLHintsLeft));
        schemBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoDrawSchematic()");
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildLSLogFtTab  -  sub-tab 4 content
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::BuildLSLogFtTab(TGCompositeFrame* p)
{
    // ── Parent parameters (can override NNDC values) ──────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Parent Parameters");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

        // Isotope selector (shared nucIsoLogftCombo_)
        TGHorizontalFrame* row0 = new TGHorizontalFrame(grp);
        grp->AddFrame(row0, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        row0->AddFrame(new TGLabel(row0, "Daughter:"),
                       new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        nucIsoLogftCombo_ = new TGComboBox(row0, 1210);
        nucIsoLogftCombo_->AddEntry("(none)", 1);
        nucIsoLogftCombo_->Select(1, kFALSE);
        nucIsoLogftCombo_->Resize(160, 22);
        row0->AddFrame(nucIsoLogftCombo_,
                       new TGLayoutHints(kLHintsExpandX | kLHintsCenterY));

        TGHorizontalFrame* row1 = new TGHorizontalFrame(grp);
        grp->AddFrame(row1, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));

        row1->AddFrame(new TGLabel(row1, "T1/2(s):"),
                       new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        lsParentHlEntry_ = new TGNumberEntry(row1, -1.0, 8, -1,
                                             TGNumberFormat::kNESRealFour,
                                             TGNumberFormat::kNEAAnyNumber);
        lsParentHlEntry_->SetWidth(70);
        row1->AddFrame(lsParentHlEntry_, new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));

        row1->AddFrame(new TGLabel(row1, "Q_beta(keV):"),
                       new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        lsParentQbEntry_ = new TGNumberEntry(row1, -1.0, 8, -1,
                                             TGNumberFormat::kNESRealFour,
                                             TGNumberFormat::kNEAAnyNumber);
        lsParentQbEntry_->SetWidth(70);
        row1->AddFrame(lsParentQbEntry_, new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));

        row1->AddFrame(new TGLabel(row1, "Z:"),
                       new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        lsParentZEntry_ = new TGNumberEntry(row1, -1.0, 4, -1,
                                            TGNumberFormat::kNESInteger,
                                            TGNumberFormat::kNEAAnyNumber);
        lsParentZEntry_->SetWidth(40);
        row1->AddFrame(lsParentZEntry_, new TGLayoutHints(kLHintsCenterY));
    }

    // ── Branching Ratios ──────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Branching Ratios");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        TGHorizontalFrame* brRow = new TGHorizontalFrame(grp);
        grp->AddFrame(brRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        brRow->AddFrame(new TGLabel(brRow, "From level:"),
                        new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        lsBrFromCombo_ = new TGComboBox(brRow, 3300);
        lsBrFromCombo_->AddEntry("(none)", 1);
        lsBrFromCombo_->Select(1, kFALSE);
        lsBrFromCombo_->Resize(120, 22);
        brRow->AddFrame(lsBrFromCombo_,
                        new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 0, 4, 0, 0));
        lsBrFromCombo_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                                "OnLSBrLevelSelected(Int_t)");

        TGTextButton* calcBRBtn = new TGTextButton(brRow, " Calculate ");
        brRow->AddFrame(calcBRBtn, new TGLayoutHints(kLHintsCenterY));
        calcBRBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLSCalcBR()");

        lsBrList_ = new TGListBox(grp, 3301);
        lsBrList_->Resize(280, 90);
        grp->AddFrame(lsBrList_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        lsBrList_->AddEntry("  (select a level and click Calculate)", 1);
    }

    // ── Log ft & Transition Strength ─────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Log ft  /  B(GT)  /  B(F)");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        TGTextButton* calcLogFtBtn = new TGTextButton(grp, " Calculate Log ft & Strengths ");
        grp->AddFrame(calcLogFtBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        calcLogFtBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLSCalcLogFt()");
        calcLogFtBtn->SetToolTipText(
            "Compute log ft (Feenberg-Trigg approx, +/-0.3) and B(GT)/B(F)\n"
            "for every level with non-zero beta-feeding.\n"
            "Requires parent T1/2, Q_beta, and Z to be set.");

        lsLogFtView_ = new TGTextView(grp, 280, 140);
        grp->AddFrame(lsLogFtView_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        lsLogFtView_->AddLine("  (click Calculate Log ft)");
    }

    // ── Cascade Balance ───────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Cascade Intensity Balance");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

        TGTextButton* calcBalBtn = new TGTextButton(grp, " Calculate Cascade Balance ");
        grp->AddFrame(calcBalBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        calcBalBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLSCalcBalance()");
        calcBalBtn->SetToolTipText(
            "For each level: check that beta-feed + cascade-in ~= gamma-out.\n"
            "Large residuals indicate missing feedings or mis-linked peaks.\n"
            "(Infrastructure for the future cascade-balance solver.)");

        lsBalanceList_ = new TGListBox(grp, 3400);
        lsBalanceList_->Resize(280, 100);
        grp->AddFrame(lsBalanceList_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        lsBalanceList_->AddEntry("  (click Calculate Cascade Balance)", 1);
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper: refresh the level list box from lsData_
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::RefreshLSLevelList()
{
    if (!lsLevelList_) return;
    lsLevelList_->RemoveAll();

    if (lsData_.empty()) {
        lsLevelList_->AddEntry("  (no levels)", 1);
        lsLevelList_->MapSubwindows(); lsLevelList_->Layout();
        return;
    }

    char buf[256];
    for (const auto& lv : lsData_.levels) {
        std::string tag = (lv.id == 0) ? " GS" : "";
        std::string typeStr = (lv.betaFeedPct > 0.0)
            ? Form("  beta=%.1f%% %s", lv.betaFeedPct, LSBetaTypeStr(lv.betaType).c_str())
            : "";
        std::string logftStr = (lv.logFt > 0.0)
            ? Form("  logft=%.2f  B(GT)=%.3g", lv.logFt, lv.B_GT) : "";

        std::snprintf(buf, sizeof(buf), "[%d] E=%7.1f  Jpi=%s%s%s%s",
                      lv.id, lv.energy,
                      lv.jpi.empty() ? "?" : lv.jpi.c_str(),
                      tag.c_str(), typeStr.c_str(), logftStr.c_str());
        lsLevelList_->AddEntry(buf, lv.id + 1);  // TGListBox id offset by 1
    }
    lsLevelList_->MapSubwindows();
    lsLevelList_->Layout();
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper: refresh the transition list box from lsData_
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::RefreshLSTransList()
{
    if (!lsTransList_) return;
    lsTransList_->RemoveAll();

    if (lsData_.transitions.empty()) {
        lsTransList_->AddEntry("  (no transitions)", 1);
        lsTransList_->MapSubwindows(); lsTransList_->Layout();
        return;
    }

    char buf[256];
    for (const auto& tr : lsData_.transitions) {
        const LSLevel* from = lsData_.FindLevel(tr.from_id);
        const LSLevel* to   = lsData_.FindLevel(tr.to_id);
        std::string fromStr = from ? Form("%.1f", from->energy) : "?";
        std::string toStr   = to   ? Form("%.1f", to->energy)   : "?";
        std::string linkStr = tr.linkedCache.empty() ? "" : " [linked]";
        std::string multStr = tr.multipolStr.empty() ? "?" : tr.multipolStr;

        std::snprintf(buf, sizeof(buf),
                      "[%d] %s->%s keV  E=%.2f  %s  I=%.1f%s",
                      tr.id, fromStr.c_str(), toStr.c_str(),
                      tr.energy, multStr.c_str(), tr.intensity, linkStr.c_str());
        lsTransList_->AddEntry(buf, tr.id + 1);
    }
    lsTransList_->MapSubwindows();
    lsTransList_->Layout();
}

// ─────────────────────────────────────────────────────────────────────────────
// Helper: rebuild From/To/BR combos from lsData_.levels
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::RefreshLSLevelCombos()
{
    auto rebuild = [&](TGComboBox* cb) {
        if (!cb) return;
        int curSel = cb->GetSelected();
        cb->RemoveAll();
        cb->AddEntry("(none)", 1);
        for (const auto& lv : lsData_.levels) {
            char buf[64];
            std::snprintf(buf, sizeof(buf), "[%d] %.1f %s",
                          lv.id, lv.energy,
                          lv.jpi.empty() ? "?" : lv.jpi.c_str());
            cb->AddEntry(buf, lv.id + 2);  // offset: 1=none, 2+=level id+2
        }
        cb->Select(curSel > 1 ? curSel : 1, kFALSE);
        cb->MapSubwindows(); cb->Layout();
    };
    rebuild(lsTransFromCombo_);
    rebuild(lsTransToCombo_);
    rebuild(lsBrFromCombo_);
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnLSSeedFromNNDC
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnLSSeedFromNNDC()
{
    if (!nucIsoCombo_) return;
    TGLBEntry* e = nucIsoCombo_->GetSelectedEntry();
    if (!e || std::string(e->GetTitle()) == "(none)") {
        AppendLog("[LevelScheme] Select an isotope first (Isotope Chain tab -> Fetch All).");
        return;
    }
    std::string isoID = e->GetTitle();

    auto it = nuclearDB_.find(isoID);
    if (it == nuclearDB_.end() || !it->second.valid()) {
        AppendLog("[LevelScheme] No NNDC data for " + isoID + ". Use 'Fetch All' in Isotope Chain.");
        return;
    }

    // Try to find parent in chain (the isotope listed just before this one)
    const NucIsotope* parent = nullptr;
    int parentZ = 0;
    for (int i = 1; i < (int)nucChainIsotopes_.size(); i++) {
        if (nucChainIsotopes_[i] == isoID) {
            auto pit = nuclearDB_.find(nucChainIsotopes_[i-1]);
            if (pit != nuclearDB_.end() && pit->second.valid()) {
                parent  = &pit->second;
                parentZ = parent->Z;
            }
            break;
        }
    }

    lsData_.SeedFromNNDC(it->second, parent, parentZ);

    // Update parent parameter entries
    if (lsParentHlEntry_)  lsParentHlEntry_->SetNumber(lsData_.parentHL_s);
    if (lsParentQbEntry_)  lsParentQbEntry_->SetNumber(lsData_.parentQbeta);
    if (lsParentZEntry_)   lsParentZEntry_->SetNumber(lsData_.parentZ);

    RefreshLSLevelList();
    RefreshLSTransList();
    RefreshLSLevelCombos();

    AppendLog(Form("[LevelScheme] Seeded %s: %d levels, %d transitions",
                   isoID.c_str(), (int)lsData_.levels.size(),
                   (int)lsData_.transitions.size()));
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnLSSave
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnLSSave()
{
    if (lsData_.empty()) { AppendLog("[LevelScheme] Nothing to save."); return; }
    EnsureLSDir(launchDir_);
    std::string path = LSCachePath(launchDir_, lsData_.isoID);
    if (lsData_.Save(path))
        AppendLog("[LevelScheme] Saved -> " + path);
    else
        AppendLog("[LevelScheme] ERROR: could not write " + path);
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnLSLoad
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnLSLoad()
{
    static const char* types[] = {
        "Level Scheme files", "*.lsdat",
        "All files",          "*",
        nullptr, nullptr
    };
    TGFileInfo fi;
    fi.fFileTypes = types;
    EnsureLSDir(launchDir_);
    std::string initDir = launchDir_ + "/level_schemes";
    fi.fIniDir = StrDup(initDir.c_str());
    OpenFileDialog(this, kFDOpen, &fi);
    if (!fi.fFilename) return;

    if (!lsData_.Load(fi.fFilename)) {
        AppendLog(std::string("[LevelScheme] Could not load ") + fi.fFilename);
        return;
    }

    if (lsParentHlEntry_)  lsParentHlEntry_->SetNumber(lsData_.parentHL_s);
    if (lsParentQbEntry_)  lsParentQbEntry_->SetNumber(lsData_.parentQbeta);
    if (lsParentZEntry_)   lsParentZEntry_->SetNumber(lsData_.parentZ);

    RefreshLSLevelList();
    RefreshLSTransList();
    RefreshLSLevelCombos();
    AppendLog(Form("[LevelScheme] Loaded %s: %d levels, %d transitions",
                   lsData_.isoID.c_str(), (int)lsData_.levels.size(),
                   (int)lsData_.transitions.size()));
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnLSAddLevel
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnLSAddLevel()
{
    double E   = lsLevelEnergyEntry_ ? lsLevelEnergyEntry_->GetNumber() : 0.0;
    std::string jpi  = lsLevelJpiEntry_   ? lsLevelJpiEntry_->GetText()   : "?";
    double feed = lsBetaFeedEntry_  ? lsBetaFeedEntry_->GetNumber()  : 0.0;
    std::string btype = lsBetaTypeEntry_  ? lsBetaTypeEntry_->GetText()   : "Unknown";

    int id = lsData_.AddLevel(E, jpi);
    LSLevel* lv = lsData_.FindLevel(id);
    if (lv) {
        lv->betaFeedPct = feed;
        lv->betaType    = LSBetaTypeFromStr(btype);
    }
    RefreshLSLevelList();
    RefreshLSLevelCombos();
    AppendLog(Form("[LevelScheme] Added level [%d] E=%.1f keV  Jpi=%s  beta=%.1f%% %s",
                   id, E, jpi.c_str(), feed, btype.c_str()));
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnLSRemoveLevel
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnLSRemoveLevel()
{
    if (!lsLevelList_) return;
    TGLBEntry* sel = lsLevelList_->GetSelectedEntry();
    if (!sel) { AppendLog("[LevelScheme] Select a level first."); return; }
    int lvId = sel->EntryId() - 1;  // undo offset
    lsData_.RemoveLevel(lvId);
    RefreshLSLevelList();
    RefreshLSTransList();
    RefreshLSLevelCombos();
    AppendLog(Form("[LevelScheme] Removed level [%d]", lvId));
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnLSLevelSelected  -  fill entry fields for editing
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnLSLevelSelected(Int_t id)
{
    int lvId = id - 1;
    const LSLevel* lv = lsData_.FindLevel(lvId);
    if (!lv) return;
    if (lsLevelEnergyEntry_) lsLevelEnergyEntry_->SetNumber(lv->energy);
    if (lsLevelJpiEntry_)    lsLevelJpiEntry_->SetText(lv->jpi.c_str());
    if (lsBetaFeedEntry_)    lsBetaFeedEntry_->SetNumber(lv->betaFeedPct);
    if (lsBetaTypeEntry_)    lsBetaTypeEntry_->SetText(LSBetaTypeStr(lv->betaType).c_str());
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnLSAddTransition
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnLSAddTransition()
{
    if (!lsTransFromCombo_ || !lsTransToCombo_) return;

    int fromSel = lsTransFromCombo_->GetSelected();
    int toSel   = lsTransToCombo_->GetSelected();
    if (fromSel <= 1 || toSel <= 1) {
        AppendLog("[LevelScheme] Select From and To levels first."); return;
    }
    int fromId = fromSel - 2;
    int toId   = toSel   - 2;

    if (fromId == toId) { AppendLog("[LevelScheme] From and To levels must differ."); return; }

    // Ensure from > to (higher energy first)
    const LSLevel* fromLv = lsData_.FindLevel(fromId);
    const LSLevel* toLv   = lsData_.FindLevel(toId);
    if (!fromLv || !toLv) return;
    if (fromLv->energy < toLv->energy) std::swap(fromId, toId);

    double E     = lsTransEnergyEntry_ ? lsTransEnergyEntry_->GetNumber() : 0.0;
    std::string mult  = lsTransMultEntry_   ? lsTransMultEntry_->GetText()     : "";
    double I     = lsTransIntEntry_    ? lsTransIntEntry_->GetNumber()    : 0.0;

    int tid = lsData_.AddTransition(fromId, toId, E);
    LSTransition* tr = lsData_.FindTransition(tid);
    if (tr) {
        tr->multipolStr = mult;
        tr->multipol    = LSMultipolFromStr(mult);
        tr->intensity   = I;
    }
    RefreshLSTransList();
    AppendLog(Form("[LevelScheme] Added transition [%d] %d->%d  E=%.2f  %s  I=%.1f",
                   tid, fromId, toId, E, mult.c_str(), I));
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnLSRemoveTransition
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnLSRemoveTransition()
{
    if (!lsTransList_) return;
    TGLBEntry* sel = lsTransList_->GetSelectedEntry();
    if (!sel) { AppendLog("[LevelScheme] Select a transition first."); return; }
    int tid = sel->EntryId() - 1;
    lsData_.RemoveTransition(tid);
    RefreshLSTransList();
    AppendLog(Form("[LevelScheme] Removed transition [%d]", tid));
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnLSTransSelected  -  fill entry fields for editing
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnLSTransSelected(Int_t id)
{
    int tid = id - 1;
    const LSTransition* tr = lsData_.FindTransition(tid);
    if (!tr) return;
    if (lsTransEnergyEntry_) lsTransEnergyEntry_->SetNumber(tr->energy);
    if (lsTransMultEntry_)   lsTransMultEntry_->SetText(tr->multipolStr.c_str());
    if (lsTransIntEntry_)    lsTransIntEntry_->SetNumber(tr->intensity);
    // Set combos to correct levels
    if (lsTransFromCombo_) lsTransFromCombo_->Select(tr->from_id + 2, kFALSE);
    if (lsTransToCombo_)   lsTransToCombo_->Select(tr->to_id   + 2, kFALSE);
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnLSLinkPeak  -  link selected Peak Table row to selected transition
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnLSLinkPeak()
{
    if (ptSelectedRow_ < 0 || ptSelectedRow_ >= (int)ptRows_.size()) {
        AppendLog("[LevelScheme] Select a peak in the Peak Table first."); return;
    }
    if (!lsTransList_) return;
    TGLBEntry* sel = lsTransList_->GetSelectedEntry();
    if (!sel) { AppendLog("[LevelScheme] Select a transition first."); return; }
    int tid = sel->EntryId() - 1;
    LSTransition* tr = lsData_.FindTransition(tid);
    if (!tr) return;

    const PeakTableRow& r = ptRows_[ptSelectedRow_];
    tr->linkedCache  = r.cacheFile;
    tr->linkedEnergy = r.energy;
    tr->intensity    = r.area;
    tr->intensityErr = r.areaErr;
    // Auto-set energy if not set
    if (tr->energy <= 0.0) tr->energy = r.energy;

    RefreshLSTransList();
    AppendLog(Form("[LevelScheme] Linked transition [%d] -> E=%.3f keV  I=%.1f+/-%.1f  [%s]",
                   tid, r.energy, r.area, r.areaErr, r.cacheFile.c_str()));
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnLSCalcBR  -  compute branching ratios for the selected level
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnLSCalcBR()
{
    if (!lsBrFromCombo_ || !lsBrList_) return;
    int sel = lsBrFromCombo_->GetSelected();
    if (sel <= 1) { AppendLog("[LevelScheme] Select a level."); return; }
    int lvId = sel - 2;

    auto brs = lsData_.BranchingRatios(lvId);
    lsBrList_->RemoveAll();

    if (brs.empty()) {
        lsBrList_->AddEntry("  (no transitions from this level with intensities)", 1);
        lsBrList_->MapSubwindows(); lsBrList_->Layout();
        return;
    }

    const LSLevel* fromLv = lsData_.FindLevel(lvId);
    std::string hdr = Form("  From [%d] E=%.1f %s:",
                           lvId,
                           fromLv ? fromLv->energy : 0.0,
                           fromLv && !fromLv->jpi.empty() ? fromLv->jpi.c_str() : "?");
    lsBrList_->AddEntry(hdr.c_str(), 1);

    char buf[256];
    int row = 2;
    for (const auto& br : brs) {
        const LSLevel* toLv = lsData_.FindLevel(br.to_id);
        std::string toStr = toLv
            ? Form("E=%.1f %s", toLv->energy, toLv->jpi.c_str())
            : Form("[%d]", br.to_id);
        std::snprintf(buf, sizeof(buf),
                      "  -> %s  I=%.2g  BR=%.3f (%.1f%%)  RelBR=%.4f",
                      toStr.c_str(), br.absIntensity,
                      br.branchingRatio, br.branchingRatio * 100.0,
                      br.relativeBR);
        lsBrList_->AddEntry(buf, row++);
    }
    lsBrList_->MapSubwindows(); lsBrList_->Layout();
    AppendLog(Form("[LevelScheme] Branching ratios from [%d]: %d branches", lvId, (int)brs.size()));
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnLSBrLevelSelected  -  recalculate BRs when combo changes
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnLSBrLevelSelected(Int_t /*id*/)
{
    OnLSCalcBR();
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnLSCalcLogFt
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnLSCalcLogFt()
{
    if (!lsLogFtView_) return;
    lsLogFtView_->Clear();

    // Apply any user overrides from the entry fields
    if (lsParentHlEntry_ && lsParentHlEntry_->GetNumber() > 0.0)
        lsData_.parentHL_s  = lsParentHlEntry_->GetNumber();
    if (lsParentQbEntry_ && lsParentQbEntry_->GetNumber() > 0.0)
        lsData_.parentQbeta = lsParentQbEntry_->GetNumber();
    if (lsParentZEntry_  && lsParentZEntry_->GetNumber()  >= 0.0)
        lsData_.parentZ     = (int)lsParentZEntry_->GetNumber();

    if (lsData_.parentHL_s <= 0.0 || lsData_.parentQbeta <= 0.0) {
        lsLogFtView_->AddLine("  ERROR: parent T1/2 and Q_beta must be > 0.");
        lsLogFtView_->Update(); return;
    }

    lsData_.ComputeAllLogFt();
    RefreshLSLevelList();  // update list with new logFt values

    lsLogFtView_->AddLine(Form("Parent: %s  T1/2=%.4g s  Q_beta=%.1f keV  Z=%d",
                                lsData_.parentID.c_str(), lsData_.parentHL_s,
                                lsData_.parentQbeta, lsData_.parentZ));
    lsLogFtView_->AddLine("──────────────────────────────────────────────────────");
    lsLogFtView_->AddLine("  Level       Jpi      beta-feed%  Type     log ft   B(GT)    B(F)");

    bool anyFed = false;
    for (const auto& lv : lsData_.levels) {
        if (lv.betaFeedPct <= 0.0) continue;
        anyFed = true;
        char buf[200];
        if (lv.logFt > 0.0) {
            std::snprintf(buf, sizeof(buf),
                "  E=%7.1f  %-6s  %6.2f%%   %-8s  %5.2f    %.3g    %.3g",
                lv.energy,
                lv.jpi.empty() ? "?" : lv.jpi.c_str(),
                lv.betaFeedPct,
                LSBetaTypeStr(lv.betaType).c_str(),
                lv.logFt, lv.B_GT, lv.B_Fermi);
        } else {
            std::snprintf(buf, sizeof(buf),
                "  E=%7.1f  %-6s  %6.2f%%   %-8s  (Q_beta too low or endpoint <= 0)",
                lv.energy,
                lv.jpi.empty() ? "?" : lv.jpi.c_str(),
                lv.betaFeedPct,
                LSBetaTypeStr(lv.betaType).c_str());
        }
        lsLogFtView_->AddLine(buf);
    }
    if (!anyFed)
        lsLogFtView_->AddLine("  (no levels with beta-feeding set)");

    lsLogFtView_->AddLine("");
    lsLogFtView_->AddLine("  Note: log ft uses Feenberg-Trigg approximation (+/-0.3).");
    lsLogFtView_->AddLine("  B(GT) = 6146.5 / (ft * g_A^2), g_A = 1.2754");
    lsLogFtView_->AddLine("  B(F)  = 6146.5 / ft (pure Fermi limit)");
    lsLogFtView_->Update();

    AppendLog("[LevelScheme] Log ft calculated for " + lsData_.isoID);
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnLSCalcBalance
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnLSCalcBalance()
{
    if (!lsBalanceList_) return;
    lsBalanceList_->RemoveAll();

    if (lsData_.empty()) {
        lsBalanceList_->AddEntry("  (no level scheme loaded)", 1);
        lsBalanceList_->MapSubwindows(); lsBalanceList_->Layout();
        return;
    }

    auto rows = lsData_.CascadeBalance();
    lsBalanceList_->AddEntry(
        "  Level(keV)    beta-in    gamma-in    gamma-out   residual%", 1);
    int id = 2;
    bool ok = true;
    for (const auto& row : rows) {
        const LSLevel* lv = lsData_.FindLevel(row.level_id);
        char buf[200];
        double resPct = row.residual * 100.0;
        std::snprintf(buf, sizeof(buf),
                      "  [%d] E=%7.1f  %6.1f  %6.1f  %6.1f   %.1f%%%s",
                      row.level_id,
                      lv ? lv->energy : 0.0,
                      row.betaIn, row.gammaIn, row.gammaOut,
                      resPct,
                      resPct > 20.0 ? "  <- imbalanced" : "");
        lsBalanceList_->AddEntry(buf, id++);
        if (resPct > 20.0) ok = false;
    }
    lsBalanceList_->MapSubwindows(); lsBalanceList_->Layout();
    AppendLog(Form("[LevelScheme] Cascade balance: %d levels%s",
                   (int)rows.size(), ok ? "  -  all balanced" : "  -  check flagged levels"));
}

// ─────────────────────────────────────────────────────────────────────────────
// DrawEnhancedLevelScheme  -  level scheme with intensity-scaled arrows,
// multipolarity labels, beta-feeding arrows, and Jpi annotations.
// Falls back to NNDC-only draw if lsData_ is empty.
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::DrawEnhancedLevelScheme(const std::string& isoID)
{
    // If no user data, fall through to NNDC draw (already called from DrawLevelScheme)
    if (lsData_.empty() || (!isoID.empty() && lsData_.isoID != isoID)) {
        DrawLevelScheme(isoID.empty() ? lsData_.isoID : isoID);
        return;
    }

    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();
    c->SetFillColor(kWhite);

    if (lsData_.levels.empty()) {
        AppendLog("[LevelScheme] No levels to draw."); return;
    }

    // Energy range
    double eMax = 0.0;
    for (const auto& lv : lsData_.levels) if (lv.energy > eMax) eMax = lv.energy;
    if (eMax <= 0.0) eMax = 100.0;
    eMax *= 1.08;  // headroom

    // Layout constants
    const double xBeta = 0.04;   // beta-feeding column left edge
    const double xL    = 0.22;   // level lines left edge
    const double xR    = 0.70;   // level lines right edge
    const double xLbl  = 0.72;   // label column
    const double yBot  = 0.08;
    const double yTop  = 0.92;

    auto eToY = [&](double e) {
        return yBot + (yTop - yBot) * (e / eMax);
    };

    TLatex tx;
    tx.SetNDC(kTRUE);
    tx.SetTextAlign(12);
    tx.SetTextSize(0.018);

    // Maximum intensity for arrow-width scaling
    double maxI = 0.0;
    for (const auto& tr : lsData_.transitions)
        if (tr.intensity > maxI) maxI = tr.intensity;

    // ── Draw level lines ──────────────────────────────────────────────────────
    for (const auto& lv : lsData_.levels) {
        double y = eToY(lv.energy);
        int col = (lv.id == 0) ? kBlack : kBlue + 1;
        TLine* ln = new TLine(xL, y, xR, y);
        ln->SetLineColor(col);
        ln->SetLineWidth(lv.id == 0 ? 2 : 1);
        ln->Draw();

        // Right label: energy + Jpi
        char buf[80];
        std::snprintf(buf, sizeof(buf), "%.1f  %s",
                      lv.energy, lv.jpi.empty() ? "?" : lv.jpi.c_str());
        tx.DrawLatex(xLbl, y, buf);

        // beta-feeding: horizontal bar on the left proportional to betaFeedPct
        if (lv.betaFeedPct > 0.0) {
            double barW = 0.12 * lv.betaFeedPct / 100.0;
            TBox* box = new TBox(xBeta, y - 0.005, xBeta + barW, y + 0.005);
            int bcol = (lv.betaType == LSBetaType::GamowTeller)  ? kMagenta + 1 :
                       (lv.betaType == LSBetaType::Fermi)         ? kGreen   + 2 :
                       (lv.betaType == LSBetaType::FirstForbidden)? kOrange  + 1 :
                                                                     kGray    + 1;
            box->SetFillColor(bcol); box->SetLineColor(bcol);
            box->Draw("f");

            char bfBuf[60];
            std::snprintf(bfBuf, sizeof(bfBuf), "%.1f%% %s",
                          lv.betaFeedPct, LSBetaTypeStr(lv.betaType).c_str());
            tx.SetTextAlign(12); tx.SetTextSize(0.015);
            tx.DrawLatex(xBeta + 0.01, y + 0.008, bfBuf);
            tx.SetTextSize(0.018);
        }
    }

    // ── Draw gamma transition arrows ──────────────────────────────────────────
    // Stagger x-positions to reduce overlap for close energies
    int nTrans = (int)lsData_.transitions.size();
    for (int ti = 0; ti < nTrans; ti++) {
        const auto& tr = lsData_.transitions[ti];
        const LSLevel* fromLv = lsData_.FindLevel(tr.from_id);
        const LSLevel* toLv   = lsData_.FindLevel(tr.to_id);
        if (!fromLv || !toLv) continue;

        double yFrom = eToY(fromLv->energy);
        double yTo   = eToY(toLv->energy);

        // Stagger: spread x across the level-line region
        double xFrac = (nTrans > 1) ? (double)ti / (double)(nTrans - 1) : 0.5;
        double xArr  = xL + (xR - xL) * (0.15 + 0.70 * xFrac);

        // Line width scaled by intensity (1..4)
        int lw = 1;
        if (maxI > 0.0 && tr.intensity > 0.0)
            lw = 1 + (int)(3.0 * tr.intensity / maxI + 0.5);

        // Color by multipolarity
        int col = kRed;
        if (tr.multipol == LSMultipol::E1 || tr.multipol == LSMultipol::E1M2)
            col = kRed;
        else if (tr.multipol == LSMultipol::E2 || tr.multipol == LSMultipol::E3)
            col = kBlue + 1;
        else if (tr.multipol == LSMultipol::M1 || tr.multipol == LSMultipol::M1E2)
            col = kGreen + 2;
        else if (tr.multipol == LSMultipol::M2 || tr.multipol == LSMultipol::M3)
            col = kOrange + 1;

        TArrow* arr = new TArrow(xArr, yFrom, xArr, yTo, 0.007, "|>");
        arr->SetLineColor(col); arr->SetFillColor(col); arr->SetLineWidth(lw);
        arr->Draw();

        // Label: energy + multipolarity
        if (tr.energy > 0.0 || !tr.multipolStr.empty()) {
            char lbuf[60];
            std::snprintf(lbuf, sizeof(lbuf), "%.1f %s",
                          tr.energy, tr.multipolStr.c_str());
            double yMid = 0.5 * (yFrom + yTo);
            tx.SetTextAlign(12); tx.SetTextSize(0.014); tx.SetTextColor(col);
            tx.DrawLatex(xArr + 0.005, yMid, lbuf);
            tx.SetTextColor(kBlack); tx.SetTextSize(0.018);
        }
    }

    // ── Title and legend ──────────────────────────────────────────────────────
    tx.SetTextSize(0.026); tx.SetTextAlign(22); tx.SetTextColor(kBlack);
    std::string title = lsData_.isoID + " Level Scheme";
    if (!lsData_.parentID.empty()) title += "  (beta- from " + lsData_.parentID + ")";
    tx.DrawLatex(0.48, 0.97, title.c_str());

    // Small legend for beta-feeding colors
    tx.SetTextSize(0.016); tx.SetTextAlign(12);
    double lx = 0.04, ly = 0.30;
    TBox* bGT = new TBox(lx, ly, lx+0.015, ly+0.010);
    bGT->SetFillColor(kMagenta+1); bGT->Draw("f");
    tx.DrawLatex(lx + 0.018, ly + 0.005, "GT");
    ly -= 0.030;
    TBox* bF = new TBox(lx, ly, lx+0.015, ly+0.010);
    bF->SetFillColor(kGreen+2); bF->Draw("f");
    tx.DrawLatex(lx + 0.018, ly + 0.005, "Fermi");
    ly -= 0.030;
    TBox* b1F = new TBox(lx, ly, lx+0.015, ly+0.010);
    b1F->SetFillColor(kOrange+1); b1F->Draw("f");
    tx.DrawLatex(lx + 0.018, ly + 0.005, "1st Forb.");

    c->Modified(); c->Update();
    AppendLog("[LevelScheme] Enhanced drawing: " + lsData_.isoID +
              Form(" (%d levels, %d transitions)", (int)lsData_.levels.size(),
                   (int)lsData_.transitions.size()));
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnLSDrawEnhanced
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnLSDrawEnhanced()
{
    if (lsData_.empty()) {
        AppendLog("[LevelScheme] No level scheme loaded. Use Seed NNDC or Load .lsdat.");
        return;
    }
    DrawEnhancedLevelScheme(lsData_.isoID);
}
