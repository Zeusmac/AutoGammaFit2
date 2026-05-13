#include "GammaFitGUI.h"
#include "GammaFitGUI_shared.h"

#include "TGTextEntry.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TSystem.h"
#include "TBox.h"
#include "TLine.h"
#include "TArrow.h"
#include "TLatex.h"
#include "TPave.h"

#include <set>
#include <map>
#include <algorithm>
#include <sstream>
#include <limits>
#include <sys/stat.h>

// Key for persisting the label→class map inside the cache file
static constexpr const char* kLabelClassesKey = "__LABEL_CLASSES__";

// Serialize / deserialize labelClassMap_ as "iso:cls;iso:cls;..."
static std::string SerializeLCMap(const std::map<std::string,std::string>& m)
{
    std::string s;
    for (auto& [k,v] : m) {
        if (k.empty() || v.empty()) continue;
        s += k + ':' + v + ';';
    }
    return s;
}
static std::map<std::string,std::string> DeserializeLCMap(const std::string& s)
{
    std::map<std::string,std::string> m;
    std::istringstream ss(s);
    std::string tok;
    while (std::getline(ss, tok, ';')) {
        auto pos = tok.find(':');
        if (pos == std::string::npos || pos == 0) continue;
        m[tok.substr(0,pos)] = tok.substr(pos+1);
    }
    return m;
}

void GammaFitGUI::LoadLabelClassMap(FitDatabase& fitdb)
{
    labelClassMap_.clear();
    auto it = fitdb.GetEntries().find(kLabelClassesKey);
    if (it != fitdb.GetEntries().end())
        labelClassMap_ = DeserializeLCMap(it->second.label);
    // Also build from per-peak classifications for any labels not in the map
    for (const auto& kv : fitdb.GetEntries()) {
        if (kv.first.empty() || kv.first[0] == '_') continue;
        const auto& e = kv.second;
        if (!e.label.empty() && !e.classification.empty()
                && !labelClassMap_.count(e.label))
            labelClassMap_[e.label] = e.classification;
    }
}

void GammaFitGUI::SaveLabelClassMap(FitDatabase& fitdb)
{
    FitEntry e;
    e.key         = kLabelClassesKey;
    e.label       = SerializeLCMap(labelClassMap_);
    e.chi2ndf     = 0;
    e.residualRMS = 0;
    e.maxPull     = 0;
    fitdb.ForceStore(kLabelClassesKey, e);
}

// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::BuildIsotopesTab(TGCompositeFrame* p)
{
    TGCanvas* sc = new TGCanvas(p, 308, 860, kSunkenFrame);
    p->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    TGCompositeFrame* cf = new TGCompositeFrame(sc->GetViewPort(), 295, 10, kVerticalFrame);
    sc->SetContainer(cf);
    p = cf;

    // ── Controls row ──────────────────────────────────────────────────────────
    {
        TGHorizontalFrame* topRow = new TGHorizontalFrame(p);
        p->AddFrame(topRow, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

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

    // ── Peak list ─────────────────────────────────────────────────────────────
    TGGroupFrame* listGrp = new TGGroupFrame(p, "Peaks (sorted by label)");
    p->AddFrame(listGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    isoList_ = new TGListBox(listGrp, 960);
    isoList_->Resize(285, 200);
    listGrp->AddFrame(isoList_,
                      new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    isoList_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                      "OnIsoListSelected(Int_t)");

    // ── Bulk action row ───────────────────────────────────────────────────────
    {
        TGHorizontalFrame* bulkRow = new TGHorizontalFrame(p);
        p->AddFrame(bulkRow, new TGLayoutHints(kLHintsExpandX, 4, 4, 0, 2));

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

    // ── Edit group ────────────────────────────────────────────────────────────
    TGGroupFrame* editGrp = new TGGroupFrame(p, "Edit Selected Peak");
    p->AddFrame(editGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    {
        TGHorizontalFrame* lblRow = new TGHorizontalFrame(editGrp);
        editGrp->AddFrame(lblRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 0));
        lblRow->AddFrame(new TGLabel(lblRow, "Label:"),
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

    // Class row — class is now the LABEL's class (applies to all peaks with same label)
    {
        TGHorizontalFrame* clsRow = new TGHorizontalFrame(editGrp);
        editGrp->AddFrame(clsRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 0));
        clsRow->AddFrame(new TGLabel(clsRow, "Class:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        isoClassCombo_ = new TGComboBox(clsRow, 970);
        isoClassCombo_->AddEntry("(none)",               1);
        isoClassCombo_->AddEntry("Parent",               2);
        isoClassCombo_->AddEntry("Daughter",             3);
        isoClassCombo_->AddEntry("Granddaughter",        4);
        isoClassCombo_->AddEntry("Beta-n Daughter",      5);
        isoClassCombo_->AddEntry("Beta-2n Daughter",     6);
        isoClassCombo_->AddEntry("Beta-n Granddaughter", 7);
        isoClassCombo_->AddEntry("Beta-2n Granddaughter",8);
        isoClassCombo_->AddEntry("Background",           9);
        isoClassCombo_->AddEntry("Custom",               10);
        isoClassCombo_->AddEntry("X-ray",                11);
        isoClassCombo_->Select(1, kFALSE);
        isoClassCombo_->Resize(150, 22);
        clsRow->AddFrame(isoClassCombo_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        isoCustomEntry_ = new TGTextEntry(clsRow, "");
        isoCustomEntry_->SetToolTipText("Custom class name");
        clsRow->AddFrame(isoCustomEntry_, new TGLayoutHints(kLHintsExpandX));
    }

    TGTextButton* applyBtn = new TGTextButton(editGrp, "Apply to Selected Peak");
    editGrp->AddFrame(applyBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 1));
    applyBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoApply()");
    applyBtn->SetToolTipText("Save label + class to the selected peak's cache entry");

    TGTextButton* applyClsAllBtn = new TGTextButton(editGrp,
        "Set Class for ALL peaks with this Label");
    editGrp->AddFrame(applyClsAllBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
    applyClsAllBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoApplyClassToAll()");
    applyClsAllBtn->SetToolTipText(
        "Apply the selected class to every peak in the cache that\n"
        "shares the currently selected label.\n"
        "This defines the class for the LABEL, not just this peak.");

    TGTextButton* applyAllBtn = new TGTextButton(editGrp, "Apply Label to All in Filter");
    editGrp->AddFrame(applyAllBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
    applyAllBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoApplyLabelAll()");
    applyAllBtn->SetToolTipText("Apply the current label to every peak currently shown in the list");

    TGTextButton* clearBtn = new TGTextButton(editGrp, "Clear Label & Class");
    editGrp->AddFrame(clearBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 4));
    clearBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoClear()");
    clearBtn->SetToolTipText("Remove label and classification from the selected peak's cache entry");

    // ── Decay schematic ───────────────────────────────────────────────────────
    TGTextButton* schemBtn = new TGTextButton(p, "Draw Decay Schematic on Canvas");
    p->AddFrame(schemBtn, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));
    schemBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoDrawSchematic()");
    schemBtn->SetToolTipText(
        "Draw a decay chain diagram on the main canvas.\n"
        "Nodes are classes (Parent → Daughter → Granddaughter).\n"
        "Only classes with labeled peaks are shown.\n"
        "Assign classes via 'Set Class for ALL peaks with this Label'.");

    // ── Isotope Database browser ──────────────────────────────────────────────
    TGGroupFrame* dbGrp = new TGGroupFrame(p, "Isotope Database");
    p->AddFrame(dbGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 4));

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
        clearDbBtn->SetToolTipText("Show all database lines");
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

// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::PopulateIsoList(const std::string& filterLabel)
{
    isoList_->RemoveAll();
    isoListKeys_.clear();
    isoListAutoMatch_.clear();
    isoListDbEnergy_.clear();

    if (isoHistName_.empty()) return;

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));

    struct Row {
        std::string key;
        double energy;
        std::string label;
        std::string classification;
    };
    std::vector<Row> rows;
    for (const auto& kv : fitdb.GetEntries()) {
        if (kv.first.empty() || kv.first[0] == '_') continue;
        const FitEntry& fe = kv.second;
        if (!filterLabel.empty() && filterLabel != "All" && fe.label != filterLabel) continue;

        int npar = (int)fe.params.size();
        if (npar >= 5 && (npar - 2) % 3 == 0) {
            int n = (npar - 2) / 3;
            for (int i = 0; i < n; i++) {
                double E = fe.params[3*i + 1];
                rows.push_back({kv.first, E, fe.label, fe.classification});
            }
        } else {
            double E = 0.0;
            try { E = std::stod(kv.first.substr(0, kv.first.find('_'))); } catch (...) {}
            rows.push_back({kv.first, E, fe.label, fe.classification});
        }
    }

    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b){
        if (a.label != b.label) return a.label < b.label;
        return a.energy < b.energy;
    });

    // ── Unique auto-matching for unlabeled peaks ──────────────────────────────
    // Greedy: sort unlabeled peaks by energy, match each to the best unclaimed DB line.
    std::set<std::string> claimedDbLines;  // "isotope_energy" keys already used

    // First claim DB lines for labeled peaks (don't let unlabeled steal them)
    if (dbLoaded_) {
        for (const auto& r : rows) {
            if (r.label.empty()) continue;
            double best = std::numeric_limits<double>::max();
            std::string bestKey;
            for (const auto& gl : db_.db) {
                if (gl.isotope != r.label) continue;
                double d = std::abs(gl.energy - r.energy);
                if (d < best) {
                    best = d;
                    bestKey = gl.isotope + Form("_%.1f", gl.energy);
                }
            }
            if (!bestKey.empty()) claimedDbLines.insert(bestKey);
        }
    }

    // Pre-compute auto-match for each unlabeled peak (unique)
    std::vector<std::string> autoMatchVec(rows.size());
    std::vector<double>      dbEnergyVec(rows.size(), 0.0);
    if (dbLoaded_) {
        for (size_t i = 0; i < rows.size(); i++) {
            if (!rows[i].label.empty() || rows[i].energy <= 0) continue;
            double fwhm = res_.FWHM(rows[i].energy);
            auto dbM = db_.Match(rows[i].energy, fwhm);
            for (const auto& m : dbM) {
                std::string key = m.isotope + Form("_%.1f", m.energy);
                if (claimedDbLines.count(key)) continue;
                autoMatchVec[i] = m.isotope;
                dbEnergyVec[i]  = m.energy;
                claimedDbLines.insert(key);
                break;
            }
        }
    }

    int id = 1;
    for (size_t idx = 0; idx < rows.size(); idx++) {
        const auto& r = rows[idx];

        // Effective label for display: confirmed > (unlabeled)
        // No "?" suffix — auto-match only shown in DB column
        std::string displayLabel;
        double      dbEnergy = 0.0;
        double      deltaE   = 0.0;

        if (!r.label.empty()) {
            displayLabel = r.label;
            // Find closest DB line for this confirmed label
            if (dbLoaded_) {
                double best = std::numeric_limits<double>::max();
                for (const auto& gl : db_.db) {
                    if (gl.isotope != r.label) continue;
                    double d = std::abs(gl.energy - r.energy);
                    if (d < best) { best = d; dbEnergy = gl.energy; }
                }
                deltaE = (dbEnergy > 0) ? r.energy - dbEnergy : 0.0;
            }
        } else {
            displayLabel = "(unlabeled)";
            // Use the unique auto-match for the DB column
            if (dbEnergyVec[idx] > 0) {
                dbEnergy = dbEnergyVec[idx];
                deltaE   = r.energy - dbEnergy;
            }
        }

        // DB column: "ISO 1173.24(+0.03)" or "  ---"
        std::string dbCol;
        if (dbEnergy > 0) {
            std::string isoName = !r.label.empty() ? r.label : autoMatchVec[idx];
            dbCol = Form("%-10s %8.2f(%+.2f)", isoName.c_str(), dbEnergy, deltaE);
        } else {
            dbCol = "---";
        }

        // Wider spacing: label(16) | energy(10.2f) | db column
        std::string entry = Form("%-16s  %8.2f    %s",
                                 displayLabel.c_str(), r.energy, dbCol.c_str());
        isoList_->AddEntry(entry.c_str(), id++);
        isoListKeys_.push_back(r.key);
        isoListAutoMatch_.push_back(autoMatchVec[idx]);
        isoListDbEnergy_.push_back(dbEnergy);
    }
    isoList_->MapSubwindows(); isoList_->Layout();
}

void GammaFitGUI::RefreshIsoDisplay()
{
    if (isoHistName_.empty()) return;

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));
    LoadLabelClassMap(fitdb);

    std::set<std::string> labels;
    for (const auto& kv : fitdb.GetEntries()) {
        if (!kv.first.empty() && kv.first[0] != '_' && !kv.second.label.empty())
            labels.insert(kv.second.label);
    }

    if (isoFilterCombo_) {
        std::string curFilter;
        if (TGLBEntry* ent = isoFilterCombo_->GetSelectedEntry())
            curFilter = ent->GetTitle();

        isoFilterCombo_->RemoveAll();
        isoFilterCombo_->AddEntry("All", 1);
        int fid = 2, restoreId = 1;
        for (const auto& l : labels) {
            if (l == curFilter) restoreId = fid;
            isoFilterCombo_->AddEntry(l.c_str(), fid++);
        }
        isoFilterCombo_->Select(restoreId, kFALSE);
    }

    PopulateIsoList();
}

void GammaFitGUI::OnIsoRefresh()
{
    isoHistName_ = currentHist_;
    if (isoHistName_.empty()) { AppendLog("Load a histogram first."); return; }
    RefreshIsoDisplay();
    AppendLog("Isotopes: loaded cache for " + isoHistName_ +
              "  (" + std::to_string(isoListKeys_.size()) + " peaks)");
}

void GammaFitGUI::OnIsoFilterChanged(Int_t /*id*/)
{
    if (!isoFilterCombo_) return;
    TGLBEntry* ent = isoFilterCombo_->GetSelectedEntry();
    std::string filter = ent ? ent->GetTitle() : "All";
    PopulateIsoList(filter == "All" ? "" : filter);
}

void GammaFitGUI::OnIsoListSelected(Int_t id)
{
    if (id < 1 || (size_t)id > isoListKeys_.size()) return;
    const std::string& key = isoListKeys_[id - 1];
    const std::string autoMatch = ((size_t)(id-1) < isoListAutoMatch_.size())
                                  ? isoListAutoMatch_[id - 1] : "";

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));
    const auto& entries = fitdb.GetEntries();
    auto it = entries.find(key);
    if (it == entries.end()) return;

    if (isoLabelCombo_) {
        isoLabelCombo_->RemoveAll();
        isoLabelCombo_->AddEntry("(none)", 1);
        int cid = 2;
        const std::string& curLabel = it->second.label;

        if (!curLabel.empty())
            isoLabelCombo_->AddEntry(curLabel.c_str(), cid++);

        int autoMatchIdx = -1;
        if (!autoMatch.empty() && autoMatch != curLabel) {
            autoMatchIdx = cid;
            isoLabelCombo_->AddEntry(autoMatch.c_str(), cid++);
        }

        if (dbLoaded_) {
            std::set<std::string> seen;
            if (!curLabel.empty())  seen.insert(curLabel);
            if (!autoMatch.empty()) seen.insert(autoMatch);
            std::vector<std::string> names;
            for (const auto& gl : db_.db)
                if (seen.insert(gl.isotope).second)
                    names.push_back(gl.isotope);
            std::sort(names.begin(), names.end());
            for (const auto& nm : names)
                isoLabelCombo_->AddEntry(nm.c_str(), cid++);
        }

        int selIdx = 1;
        if (!curLabel.empty())     selIdx = 2;
        else if (autoMatchIdx > 0) selIdx = autoMatchIdx;
        isoLabelCombo_->Select(selIdx, kFALSE);
        isoLabelCombo_->MapSubwindows(); isoLabelCombo_->Layout();
    }

    if (isoCustomLabelEntry_) {
        const std::string& curLabel = it->second.label;
        bool inDB = false;
        if (dbLoaded_ && !curLabel.empty())
            for (const auto& gl : db_.db)
                if (gl.isotope == curLabel) { inDB = true; break; }
        isoCustomLabelEntry_->SetText((!curLabel.empty() && !inDB) ? curLabel.c_str() : "");
    }

    if (isoClassCombo_) {
        // Class comes from labelClassMap_ for the peak's label, else stored classification
        const std::string& lbl = it->second.label;
        std::string cls;
        if (!lbl.empty() && labelClassMap_.count(lbl))
            cls = labelClassMap_.at(lbl);
        else
            cls = it->second.classification;
        // Auto-suggest X-ray for peaks below 100 keV
        double peakE = (it->second.params.size() >= 2) ? it->second.params[1] : 0.0;
        if (cls.empty() && peakE > 0 && peakE < 100.0) cls = "X-ray";
        isoClassCombo_->Select(ClassToComboIndex(cls), kFALSE);
    }
    if (isoCustomEntry_) {
        const std::string& cls = it->second.classification;
        if (cls.size() > 7 && cls.substr(0,7) == "Custom:")
            isoCustomEntry_->SetText(cls.substr(7).c_str());
        else
            isoCustomEntry_->SetText("");
    }
}

void GammaFitGUI::OnIsoApply()
{
    if (isoHistName_.empty()) { AppendLog("Refresh cache first."); return; }
    Int_t sel = isoList_->GetSelected();
    if (sel < 1 || (size_t)sel > isoListKeys_.size()) {
        AppendLog("Select a peak from the list first."); return;
    }
    const std::string& key = isoListKeys_[sel - 1];

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));
    const auto& entries = fitdb.GetEntries();
    auto it = entries.find(key);
    if (it == entries.end()) { AppendLog("Cache entry not found."); return; }

    FitEntry e = it->second;
    {
        std::string customLbl = isoCustomLabelEntry_ ? isoCustomLabelEntry_->GetText() : "";
        if (!customLbl.empty()) {
            e.label = customLbl;
        } else if (isoLabelCombo_) {
            TGLBEntry* le = isoLabelCombo_->GetSelectedEntry();
            std::string lbl = le ? le->GetTitle() : "";
            e.label = (lbl != "(none)") ? lbl : "";
        }
    }
    if (isoClassCombo_) {
        std::string cust = isoCustomEntry_ ? isoCustomEntry_->GetText() : "";
        e.classification = ClassToString(isoClassCombo_->GetSelected(), cust);
    }

    fitdb.ForceStore(key, e);
    mkdir(kCacheDir, 0755);
    fitdb.Save(CacheFileFor(isoHistName_));
    AppendLog("Updated: " + key + "  label=" + e.label + "  class=" + e.classification);
    RefreshIsoDisplay();
}

void GammaFitGUI::OnIsoApplyClassToAll()
{
    if (isoHistName_.empty()) { AppendLog("Refresh cache first."); return; }

    // Determine the label to apply the class to
    std::string targetLabel;
    {
        std::string customLbl = isoCustomLabelEntry_ ? isoCustomLabelEntry_->GetText() : "";
        if (!customLbl.empty()) {
            targetLabel = customLbl;
        } else if (isoLabelCombo_) {
            TGLBEntry* le = isoLabelCombo_->GetSelectedEntry();
            std::string lbl = le ? le->GetTitle() : "";
            if (lbl != "(none)") targetLabel = lbl;
        }
    }
    if (targetLabel.empty()) {
        AppendLog("Select a label first."); return;
    }

    std::string newCls;
    if (isoClassCombo_) {
        std::string cust = isoCustomEntry_ ? isoCustomEntry_->GetText() : "";
        newCls = ClassToString(isoClassCombo_->GetSelected(), cust);
    }

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));

    // Update labelClassMap_ for this label
    if (!newCls.empty()) labelClassMap_[targetLabel] = newCls;

    // Apply class to ALL peaks with this label
    int nUpdated = 0;
    for (const auto& kv : fitdb.GetEntries()) {
        if (kv.first.empty() || kv.first[0] == '_') continue;
        if (kv.second.label != targetLabel) continue;
        FitEntry e = kv.second;
        e.classification = newCls;
        fitdb.ForceStore(kv.first, e);
        ++nUpdated;
    }

    SaveLabelClassMap(fitdb);
    mkdir(kCacheDir, 0755);
    fitdb.Save(CacheFileFor(isoHistName_));
    AppendLog("Set class '" + newCls + "' for all " + std::to_string(nUpdated) +
              " peaks labeled '" + targetLabel + "'");
    RefreshIsoDisplay();
}

void GammaFitGUI::OnIsoApplyLabelAll()
{
    if (isoHistName_.empty() || isoListKeys_.empty()) {
        AppendLog("Refresh cache first."); return;
    }
    std::string newLabel;
    if (isoLabelCombo_) {
        TGLBEntry* le = isoLabelCombo_->GetSelectedEntry();
        std::string t = le ? le->GetTitle() : "";
        if (t != "(none)") newLabel = t;
    }
    std::string cust   = isoCustomEntry_ ? std::string(isoCustomEntry_->GetText()) : "";
    std::string newCls = isoClassCombo_ ? ClassToString(isoClassCombo_->GetSelected(), cust) : "";

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));
    for (const auto& key : isoListKeys_) {
        auto it = fitdb.GetEntries().find(key);
        if (it == fitdb.GetEntries().end()) continue;
        FitEntry e = it->second;
        if (!newLabel.empty()) e.label = newLabel;
        if (!newCls.empty())   e.classification = newCls;
        fitdb.ForceStore(key, e);
    }
    if (!newLabel.empty() && !newCls.empty())
        labelClassMap_[newLabel] = newCls;
    SaveLabelClassMap(fitdb);
    mkdir(kCacheDir, 0755);
    fitdb.Save(CacheFileFor(isoHistName_));
    AppendLog("Applied label/class to " + std::to_string(isoListKeys_.size()) + " peaks.");
    RefreshIsoDisplay();
}

void GammaFitGUI::OnIsoClear()
{
    if (isoHistName_.empty()) { AppendLog("Refresh cache first."); return; }
    Int_t sel = isoList_->GetSelected();
    if (sel < 1 || (size_t)sel > isoListKeys_.size()) {
        AppendLog("Select a peak from the list first."); return;
    }
    const std::string& key = isoListKeys_[sel - 1];
    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));
    auto it = fitdb.GetEntries().find(key);
    if (it == fitdb.GetEntries().end()) { AppendLog("Cache entry not found."); return; }
    FitEntry e = it->second;
    e.label          = "";
    e.classification = "";
    fitdb.ForceStore(key, e);
    mkdir(kCacheDir, 0755);
    fitdb.Save(CacheFileFor(isoHistName_));
    if (isoLabelCombo_)       isoLabelCombo_->Select(1, kFALSE);
    if (isoCustomLabelEntry_) isoCustomLabelEntry_->SetText("");
    if (isoClassCombo_)       isoClassCombo_->Select(1, kFALSE);
    if (isoCustomEntry_)      isoCustomEntry_->SetText("");
    AppendLog("Cleared label/class for: " + key);
    RefreshIsoDisplay();
}

void GammaFitGUI::OnIsoClearAll()
{
    if (isoHistName_.empty()) { AppendLog("Refresh cache first."); return; }
    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));
    int n = 0;
    for (const auto& kv : fitdb.GetEntries()) {
        if (kv.first.empty() || kv.first[0] == '_') continue;
        if (kv.second.label.empty() && kv.second.classification.empty()) continue;
        FitEntry e = kv.second;
        e.label = ""; e.classification = "";
        fitdb.ForceStore(kv.first, e);
        ++n;
    }
    labelClassMap_.clear();
    SaveLabelClassMap(fitdb);
    mkdir(kCacheDir, 0755);
    fitdb.Save(CacheFileFor(isoHistName_));
    AppendLog("Cleared labels/classes for " + std::to_string(n) + " peaks.");
    RefreshIsoDisplay();
}

void GammaFitGUI::OnIsoAutoMatchAll()
{
    if (isoHistName_.empty()) { AppendLog("Refresh cache first."); return; }
    if (!dbLoaded_) { AppendLog("Load isotope DB first."); return; }

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));

    // Collect all peaks and their energies
    struct PeakInfo { std::string key; double energy; bool hasLabel; };
    std::vector<PeakInfo> allPeaks;
    for (const auto& kv : fitdb.GetEntries()) {
        if (kv.first.empty() || kv.first[0] == '_') continue;
        const FitEntry& fe = kv.second;
        int npar = (int)fe.params.size();
        if (npar >= 5 && (npar - 2) % 3 == 0) {
            int n = (npar - 2) / 3;
            for (int i = 0; i < n; i++)
                allPeaks.push_back({kv.first, fe.params[3*i+1], !fe.label.empty()});
        } else {
            double E = 0.0;
            try { E = std::stod(kv.first.substr(0, kv.first.find('_'))); } catch (...) {}
            allPeaks.push_back({kv.first, E, !fe.label.empty()});
        }
    }

    // Claim DB lines for already-labeled peaks first
    std::set<std::string> claimedLines;
    for (const auto& pk : allPeaks) {
        if (!pk.hasLabel) continue;
        const auto& e = fitdb.GetEntries().at(pk.key);
        if (e.label.empty()) continue;
        double best = std::numeric_limits<double>::max();
        std::string bestKey;
        for (const auto& gl : db_.db) {
            if (gl.isotope != e.label) continue;
            double d = std::abs(gl.energy - pk.energy);
            if (d < best) { best = d; bestKey = gl.isotope + Form("_%.1f", gl.energy); }
        }
        if (!bestKey.empty()) claimedLines.insert(bestKey);
    }

    // Sort unlabeled peaks by energy for deterministic matching
    std::vector<PeakInfo*> unlabeled;
    for (auto& pk : allPeaks)
        if (!pk.hasLabel && pk.energy > 0) unlabeled.push_back(&pk);
    std::sort(unlabeled.begin(), unlabeled.end(),
              [](const PeakInfo* a, const PeakInfo* b){ return a->energy < b->energy; });

    int nMatched = 0;
    for (auto* pk : unlabeled) {
        double fwhm = res_.FWHM(pk->energy);
        auto dbM = db_.Match(pk->energy, fwhm);
        for (const auto& m : dbM) {
            std::string dbKey = m.isotope + Form("_%.1f", m.energy);
            if (claimedLines.count(dbKey)) continue;
            // Apply this match
            auto it = fitdb.GetEntries().find(pk->key);
            if (it == fitdb.GetEntries().end()) break;
            FitEntry e = it->second;
            e.label = m.isotope;
            fitdb.ForceStore(pk->key, e);
            claimedLines.insert(dbKey);
            ++nMatched;
            break;
        }
    }

    mkdir(kCacheDir, 0755);
    fitdb.Save(CacheFileFor(isoHistName_));
    AppendLog("Auto-matched " + std::to_string(nMatched) + " unlabeled peaks (unique).");
    RefreshIsoDisplay();
}

// ─────────────────────────────────────────────────────────────────────────────
// Decay chain schematic
// ─────────────────────────────────────────────────────────────────────────────

// Class hierarchy columns:
//   col 0 (left):   Parent
//   col 1 (middle): Daughter, Beta-n Daughter, Beta-2n Daughter
//   col 2 (right):  Granddaughter, Beta-n Granddaughter, Beta-2n Granddaughter
// Arrow connections (col0→col1, col1→col2 per chain):
//   Parent → Daughter → Granddaughter
//   Parent → Beta-n Daughter → Beta-n Granddaughter
//   Parent → Beta-2n Daughter → Beta-2n Granddaughter

static const std::vector<std::tuple<std::string,int>> kClassLayout = {
    {"Parent",               0},
    {"Daughter",             1},
    {"Beta-n Daughter",      1},
    {"Beta-2n Daughter",     1},
    {"Granddaughter",        2},
    {"Beta-n Granddaughter", 2},
    {"Beta-2n Granddaughter",2},
};
// Chain linkage: {parent_class, child_class}
static const std::vector<std::pair<std::string,std::string>> kChainLinks = {
    {"Parent",          "Daughter"},
    {"Parent",          "Beta-n Daughter"},
    {"Parent",          "Beta-2n Daughter"},
    {"Daughter",        "Granddaughter"},
    {"Beta-n Daughter", "Beta-n Granddaughter"},
    {"Beta-2n Daughter","Beta-2n Granddaughter"},
};

void GammaFitGUI::DrawDecaySchematic(TCanvas* c)
{
    if (isoHistName_.empty()) return;

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));
    LoadLabelClassMap(fitdb);

    // Build: class → list of {isotope, peak_count}
    std::map<std::string, std::map<std::string,int>> classIso;  // class→{iso→count}
    for (const auto& kv : fitdb.GetEntries()) {
        if (kv.first.empty() || kv.first[0] == '_') continue;
        const FitEntry& fe = kv.second;
        if (fe.label.empty()) continue;
        std::string cls;
        if (labelClassMap_.count(fe.label))     cls = labelClassMap_.at(fe.label);
        else if (!fe.classification.empty())    cls = fe.classification;
        if (cls.empty()) continue;

        int nPeaks = 1;
        int npar = (int)fe.params.size();
        if (npar >= 5 && (npar - 2) % 3 == 0) nPeaks = (npar - 2) / 3;
        classIso[cls][fe.label] += nPeaks;
    }

    // Collect nodes that actually have data
    struct Node {
        std::string cls;
        int col;
        std::vector<std::pair<std::string,int>> isos;  // {name, peak count}
        double xc = 0, yc = 0, h = 0;  // center x/y and height (NDC)
    };
    std::vector<Node> nodes;
    for (auto& [cls, col] : kClassLayout) {
        auto it = classIso.find(cls);
        if (it == classIso.end()) continue;
        Node nd;
        nd.cls = cls; nd.col = col;
        for (auto& [iso, cnt] : it->second)
            nd.isos.push_back({iso, cnt});
        std::sort(nd.isos.begin(), nd.isos.end());
        nodes.push_back(nd);
    }

    if (nodes.empty()) {
        AppendLog("No labeled peaks with class assignments. "
                  "Use 'Set Class for ALL peaks with this Label' to assign classes.");
        return;
    }

    c->Clear(); c->cd();
    c->SetFillColor(kWhite);
    gPad->SetFillColor(kWhite);

    // Layout: x positions per column
    const double xCol[3] = {0.12, 0.50, 0.88};
    const double boxW    = 0.22;   // half-width
    const double lineH   = 0.040; // per text line
    const double padV    = 0.020; // vertical padding inside box

    // Assign y centers per column
    std::map<int,std::vector<Node*>> byCols;
    for (auto& nd : nodes) byCols[nd.col].push_back(&nd);

    for (auto& [col, ndlist] : byCols) {
        // Height of each node depends on number of isotope lines
        for (auto* nd : ndlist)
            nd->h = (nd->isos.size() + 1) * lineH + 2 * padV;
        // Total height needed
        double totalH = 0;
        for (auto* nd : ndlist) totalH += nd->h;
        double spacing = (ndlist.size() > 1)
            ? (0.85 - 0.10 - totalH) / (ndlist.size() - 1) : 0.0;
        double y = 0.90;
        for (auto* nd : ndlist) {
            nd->xc = xCol[col];
            nd->yc = y - nd->h / 2;
            y -= nd->h + spacing;
        }
    }

    // Draw boxes and text
    TLatex tx; tx.SetNDC(kTRUE); tx.SetTextAlign(22);

    // Box colors by column
    const int kBoxFill[3] = {kOrange-9, kAzure-9, kGreen-9};

    for (const auto& nd : nodes) {
        double x1 = nd.xc - boxW / 2, x2 = nd.xc + boxW / 2;
        double y1 = nd.yc - nd.h / 2, y2 = nd.yc + nd.h / 2;

        TBox* box = new TBox(x1, y1, x2, y2);
        box->SetFillColor(kBoxFill[nd.col]);
        box->SetLineColor(kBlack);
        box->SetLineWidth(1);
        box->Draw();

        // Class name (bold header)
        tx.SetTextSize(0.030); tx.SetTextColor(kBlack);
        double ty = y2 - padV - lineH * 0.5;
        tx.DrawLatex(nd.xc, ty, nd.cls.c_str());

        // Isotope lines
        tx.SetTextSize(0.025); tx.SetTextColor(kBlue + 1);
        for (const auto& [iso, cnt] : nd.isos) {
            ty -= lineH;
            tx.DrawLatex(nd.xc, ty,
                Form("%s (%d peak%s)", iso.c_str(), cnt, cnt == 1 ? "" : "s"));
        }
    }

    // Draw chain arrows
    // Build a map: class name → node pointer for quick lookup
    std::map<std::string, const Node*> nodeMap;
    for (const auto& nd : nodes) nodeMap[nd.cls] = &nd;

    for (const auto& [fromCls, toCls] : kChainLinks) {
        auto fi = nodeMap.find(fromCls);
        auto ti = nodeMap.find(toCls);
        if (fi == nodeMap.end() || ti == nodeMap.end()) continue;
        const Node& from = *fi->second;
        const Node& to   = *ti->second;

        // Arrow from right edge of 'from' to left edge of 'to'
        double ax1 = from.xc + boxW / 2;
        double ax2 = to.xc   - boxW / 2;
        double ay  = (from.yc + to.yc) / 2;
        // Clamp y to both boxes
        ay = std::max(ay, std::max(from.yc - from.h/2 + padV, to.yc - to.h/2 + padV));
        ay = std::min(ay, std::min(from.yc + from.h/2 - padV, to.yc + to.h/2 - padV));

        TArrow* arr = new TArrow(ax1, ay, ax2, ay, 0.012, "|>");
        arr->SetLineColor(kBlack);
        arr->SetFillColor(kBlack);
        arr->SetLineWidth(2);
        arr->Draw();
    }

    // Title
    tx.SetTextSize(0.035); tx.SetTextColor(kBlack); tx.SetTextAlign(22);
    tx.DrawLatex(0.50, 0.97,
        Form("Decay Schematic: %s", isoHistName_.c_str()));

    c->Modified(); c->Update();
}

void GammaFitGUI::OnIsoDrawSchematic()
{
    if (isoHistName_.empty()) {
        AppendLog("Click 'Refresh' to load a histogram's cache first."); return;
    }
    DrawDecaySchematic(canvas_->GetCanvas());
    AppendLog("Decay schematic drawn for " + isoHistName_);
}

// ─────────────────────────────────────────────────────────────────────────────
// Isotope DB browser
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::PopulateIsoDbList(const std::string& filter)
{
    if (!isoDbList_) return;
    isoDbList_->RemoveAll();
    isoDbIndices_.clear();

    std::string lf = filter;
    std::transform(lf.begin(), lf.end(), lf.begin(), ::tolower);

    int id = 1;
    for (int i = 0; i < (int)db_.db.size(); ++i) {
        const GammaLine& gl = db_.db[i];
        std::string iso_lc = gl.isotope;
        std::transform(iso_lc.begin(), iso_lc.end(), iso_lc.begin(), ::tolower);
        std::string estr = Form("%.1f", gl.energy);
        if (!lf.empty() && iso_lc.find(lf) == std::string::npos
                        && estr.find(filter) == std::string::npos)
            continue;
        std::string intStr = gl.hasIntensity ? Form("%.1f%%", gl.intensity) : "N/A";
        std::string entry  = Form("%-14s  %8.2f keV  I=%-7s",
                                  gl.isotope.c_str(), gl.energy, intStr.c_str());
        isoDbList_->AddEntry(entry.c_str(), id++);
        isoDbIndices_.push_back(i);
    }
    isoDbList_->MapSubwindows();
    isoDbList_->Layout();
}

void GammaFitGUI::OnIsoDbSearch()
{
    std::string filter = isoDbSearch_ ? std::string(isoDbSearch_->GetText()) : "";
    PopulateIsoDbList(filter);
}

void GammaFitGUI::OnIsoDbClear()
{
    if (isoDbSearch_) isoDbSearch_->SetText("");
    PopulateIsoDbList();
}

void GammaFitGUI::OnIsoDbApply()
{
    Int_t dbId = isoDbList_ ? isoDbList_->GetSelected() : -1;
    if (dbId < 1 || (size_t)dbId > isoDbIndices_.size()) {
        AppendLog("Select a line from the Isotope Database list first."); return;
    }
    Int_t peakId = isoList_ ? isoList_->GetSelected() : -1;
    if (peakId < 1 || (size_t)peakId > isoListKeys_.size()) {
        AppendLog("Select a fitted peak from the Peaks list first."); return;
    }
    if (isoHistName_.empty()) {
        AppendLog("Refresh the Isotopes tab for a histogram first."); return;
    }

    const GammaLine& gl = db_.db[isoDbIndices_[dbId - 1]];
    const std::string& key = isoListKeys_[peakId - 1];

    std::string newCls;
    if (isoDbClassCombo_) newCls = ClassToString(isoDbClassCombo_->GetSelected(), "");

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));
    auto it = fitdb.GetEntries().find(key);
    if (it == fitdb.GetEntries().end()) {
        AppendLog("Cache entry not found for key: " + key); return;
    }

    FitEntry e = it->second;
    e.label = gl.isotope;
    if (!newCls.empty()) {
        e.classification = newCls;
        labelClassMap_[gl.isotope] = newCls;
    }
    fitdb.ForceStore(key, e);
    SaveLabelClassMap(fitdb);
    mkdir(kCacheDir, 0755);
    fitdb.Save(CacheFileFor(isoHistName_));

    double fittedE = (e.params.size() >= 2) ? e.params[1] : 0.0;
    AppendLog("Matched: " + key + "  →  " + gl.isotope +
              Form("  DB: %.2f keV  dE=%+.2f keV", gl.energy, fittedE - gl.energy) +
              (newCls.empty() ? "" : "  [" + newCls + "]"));
    RefreshIsoDisplay();
}

void GammaFitGUI::OnIsoDbLineSelected(Int_t id)
{
    if (!isoDbList_ || id < 1 || (size_t)id > isoDbIndices_.size()) return;
    const GammaLine& gl = db_.db[isoDbIndices_[id - 1]];

    if (isoLabelCombo_) {
        int found = -1;
        TGListBox* lb = isoLabelCombo_->GetListBox();
        for (int i = 1; lb && i <= lb->GetNumberOfEntries(); ++i) {
            TGLBEntry* e = lb->GetEntry(i);
            if (!e) continue;
            if (std::string(e->GetTitle()) == gl.isotope) { found = i; break; }
        }
        if (found < 0) {
            int nxt = isoLabelCombo_->GetNumberOfEntries() + 1;
            isoLabelCombo_->AddEntry(gl.isotope.c_str(), nxt);
            found = nxt;
        }
        isoLabelCombo_->Select(found, kFALSE);
        isoLabelCombo_->MapSubwindows(); isoLabelCombo_->Layout();
    }
    AppendLog(Form("DB line selected: %s  %.2f keV  I=%.1f%%",
                   gl.isotope.c_str(), gl.energy, gl.intensity));
}
