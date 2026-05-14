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

// Keys for persisting data inside the cache file
static constexpr const char* kLabelClassesKey = "__LABEL_CLASSES__";
static constexpr const char* kParentInfoKey   = "__PARENT_INFO__";

// Element symbols Z=1..118
static const char* const kElementSymbols[] = {
    "H",  "He", "Li", "Be", "B",  "C",  "N",  "O",  "F",  "Ne",
    "Na", "Mg", "Al", "Si", "P",  "S",  "Cl", "Ar", "K",  "Ca",
    "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Ga", "Ge", "As", "Se", "Br", "Kr", "Rb", "Sr", "Y",  "Zr",
    "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In", "Sn",
    "Sb", "Te", "I",  "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd",
    "Pm", "Sm", "Eu", "Gd", "Tb", "Dy", "Ho", "Er", "Tm", "Yb",
    "Lu", "Hf", "Ta", "W",  "Re", "Os", "Ir", "Pt", "Au", "Hg",
    "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th",
    "Pa", "U",  "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es", "Fm",
    "Md", "No", "Lr", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds",
    "Rg", "Cn", "Nh", "Fl", "Mc", "Lv", "Ts", "Og"
};
static constexpr int kNElements = 118;

static const char* ElementSymbol(int Z) {
    if (Z < 1 || Z > kNElements) return "?";
    return kElementSymbols[Z - 1];
}

static int SymbolToZ(const std::string& sym) {
    for (int z = 1; z <= kNElements; z++)
        if (sym == kElementSymbols[z - 1]) return z;
    return 0;
}

// Beta-minus decay chain: given parent (Z0,N0), compute (Z,N) for a class.
// Each step: Z+1, N-1 (beta-minus emission).
// Returns false if parent not set or class is not in the beta-minus chain.
static bool ClassZN(const std::string& cls, int Z0, int N0, int& Z, int& N)
{
    if (Z0 <= 0 || N0 <= 0) return false;
    if      (cls == "Parent")                { Z = Z0;   N = N0;   }
    else if (cls == "Daughter")              { Z = Z0+1; N = N0-1; }
    else if (cls == "Granddaughter"         ||
             cls == "Beta-n Daughter")       { Z = Z0+2; N = N0-2; }
    else if (cls == "Beta-n Granddaughter"  ||
             cls == "Beta-2n Daughter")      { Z = Z0+3; N = N0-3; }
    else if (cls == "Beta-2n Granddaughter") { Z = Z0+4; N = N0-4; }
    else return false;
    return (Z > 0 && N > 0);
}

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

    // Load parent nucleus info from __PARENT_INFO__ entry
    isoParentIsotope_.clear(); isoParentZval_ = 0; isoParentNval_ = 0;
    auto pit = fitdb.GetEntries().find(kParentInfoKey);
    if (pit != fitdb.GetEntries().end()) {
        // Format: "name=Kr-90;Z=36;N=54"
        std::istringstream ss(pit->second.label);
        std::string tok;
        while (std::getline(ss, tok, ';')) {
            auto eq = tok.find('=');
            if (eq == std::string::npos) continue;
            std::string k = tok.substr(0, eq), v = tok.substr(eq + 1);
            if      (k == "name") isoParentIsotope_ = v;
            else if (k == "Z")    try { isoParentZval_ = std::stoi(v); } catch (...) {}
            else if (k == "N")    try { isoParentNval_ = std::stoi(v); } catch (...) {}
        }
        if (isoParentNameEntry_) isoParentNameEntry_->SetText(isoParentIsotope_.c_str());
        if (isoParentZEntry_)    isoParentZEntry_->SetNumber(isoParentZval_);
        if (isoParentNEntry_)    isoParentNEntry_->SetNumber(isoParentNval_);
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
    TGGroupFrame* listGrp = new TGGroupFrame(p, "Peaks (sorted by energy)");
    p->AddFrame(listGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    isoList_ = new TGListBox(listGrp, 960);
    isoList_->Resize(285, 200);
    listGrp->AddFrame(isoList_,
                      new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    isoList_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                      "OnIsoListSelected(Int_t)");

    {
        TGTextButton* previewBtn = new TGTextButton(listGrp, "Preview Selected Peak");
        listGrp->AddFrame(previewBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        previewBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoPeakPreview()");
        previewBtn->SetToolTipText(
            "Zoom the main canvas to the selected peak and overlay its cached fit.");
    }

    // ── Bulk action row ───────────────────────────────────────────────────────
    {
        TGHorizontalFrame* bulkRow = new TGHorizontalFrame(p);
        p->AddFrame(bulkRow, new TGLayoutHints(kLHintsExpandX, 4, 4, 0, 2));

        TGTextButton* rematchBtn = new TGTextButton(bulkRow, "Auto Re-match All");
        bulkRow->AddFrame(rematchBtn, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        rematchBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoAutoMatchAll()");
        rematchBtn->SetToolTipText(
            "For every unlabeled peak, find the closest unused DB line within\n"
            "the match distance threshold and apply it as the label.\n"
            "Peaks > 4 keV apart in a multi-Gaussian fit are matched individually.");

        TGTextButton* clearAllBtn = new TGTextButton(bulkRow, "Clear All Isotopes");
        bulkRow->AddFrame(clearAllBtn, new TGLayoutHints(kLHintsLeft));
        clearAllBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoClearAll()");
        clearAllBtn->SetToolTipText("Remove isotope label and classification from every peak in the cache");
    }
    // ── Match threshold ───────────────────────────────────────────────────────
    {
        TGHorizontalFrame* thrRow = new TGHorizontalFrame(p);
        p->AddFrame(thrRow, new TGLayoutHints(kLHintsExpandX, 4, 4, 0, 4));
        thrRow->AddFrame(new TGLabel(thrRow, "Match threshold (keV):"),
                         new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        isoMatchThreshEntry_ = new TGNumberEntry(thrRow, 10.0, 6, -1,
            TGNumberFormat::kNESRealOne,
            TGNumberFormat::kNEAPositive,
            TGNumberFormat::kNELLimitMin, 0.1);
        thrRow->AddFrame(isoMatchThreshEntry_, new TGLayoutHints(kLHintsLeft));
    }
    {
        TGTextButton* loadMatchBtn = new TGTextButton(p, "Load Matches to Histogram");
        p->AddFrame(loadMatchBtn, new TGLayoutHints(kLHintsExpandX, 4, 4, 0, 4));
        loadMatchBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoApplyAutoMatches()");
        loadMatchBtn->SetToolTipText(
            "Auto-match unlabeled peaks within threshold, save labels to cache,\n"
            "then redraw the histogram with isotope names shown above each peak energy.");
    }

    // ── Edit group ────────────────────────────────────────────────────────────
    TGGroupFrame* editGrp = new TGGroupFrame(p, "Edit Selected Peak");
    p->AddFrame(editGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

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

    TGTextButton* applyBtn = new TGTextButton(editGrp, "Apply Isotope to Selected Peak");
    editGrp->AddFrame(applyBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 1));
    applyBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoApply()");
    applyBtn->SetToolTipText("Save the chosen isotope label to the selected peak's cache entry");

    TGTextButton* setLblDecayBtn = new TGTextButton(editGrp, "Set Isotope & Decay Type...");
    editGrp->AddFrame(setLblDecayBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 1));
    setLblDecayBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoSetLabelDecay()");
    setLblDecayBtn->SetToolTipText(
        "Open a dialog to assign a decay type (class) to any DB isotope.\n"
        "No peak selection needed — assigns the class globally for that isotope\n"
        "and propagates it to ALL peaks in the cache labeled with it.");

    TGTextButton* clearBtn = new TGTextButton(editGrp, "Clear Isotope");
    editGrp->AddFrame(clearBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 1, 4));
    clearBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoClear()");
    clearBtn->SetToolTipText("Remove isotope label and classification from the selected peak's cache entry");

    // ── Decay schematic ───────────────────────────────────────────────────────
    TGTextButton* schemBtn = new TGTextButton(p, "Draw Decay Schematic on Canvas");
    p->AddFrame(schemBtn, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));
    schemBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoDrawSchematic()");
    schemBtn->SetToolTipText(
        "Draw a decay chain diagram on the main canvas.\n"
        "Nodes are classes (Parent → Daughter → Granddaughter).\n"
        "Only classes with labeled peaks are shown.\n"
        "Assign classes via 'Set Class for ALL peaks with this Label'.");

    // ── Parent Nucleus ────────────────────────────────────────────────────────
    TGGroupFrame* parentGrp = new TGGroupFrame(p, "Parent Nucleus (#beta^{-} chain)");
    p->AddFrame(parentGrp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

    {
        TGHorizontalFrame* nameRow = new TGHorizontalFrame(parentGrp);
        parentGrp->AddFrame(nameRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 4, 2));
        nameRow->AddFrame(new TGLabel(nameRow, "Isotope:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        isoParentNameEntry_ = new TGTextEntry(nameRow, "");
        isoParentNameEntry_->SetToolTipText("Parent isotope name, e.g. Kr-90");
        nameRow->AddFrame(isoParentNameEntry_, new TGLayoutHints(kLHintsExpandX));
    }
    {
        TGHorizontalFrame* znRow = new TGHorizontalFrame(parentGrp);
        parentGrp->AddFrame(znRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 2));
        znRow->AddFrame(new TGLabel(znRow, "Z:"),
                        new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        isoParentZEntry_ = new TGNumberEntry(znRow, 0, 4, -1,
            TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
            TGNumberFormat::kNELLimitMin, 0.0, 118.0);
        znRow->AddFrame(isoParentZEntry_, new TGLayoutHints(kLHintsLeft, 0, 10, 0, 0));
        znRow->AddFrame(new TGLabel(znRow, "N:"),
                        new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        isoParentNEntry_ = new TGNumberEntry(znRow, 0, 4, -1,
            TGNumberFormat::kNESInteger, TGNumberFormat::kNEANonNegative,
            TGNumberFormat::kNELLimitMin, 0.0, 300.0);
        znRow->AddFrame(isoParentNEntry_, new TGLayoutHints(kLHintsLeft));
    }
    {
        TGTextButton* setParentBtn = new TGTextButton(parentGrp, "Set Parent Info");
        parentGrp->AddFrame(setParentBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));
        setParentBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoSetParent()");
        setParentBtn->SetToolTipText(
            "Save parent isotope name, Z, and N to the cache.\n"
            "Daughter/granddaughter Z and N are auto-calculated\n"
            "via beta-minus decay steps when drawing the schematic.");
    }

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
    Int_t savedScroll = 0;
    if (TGScrollBar* vsb = isoList_->GetVScrollbar())
        savedScroll = vsb->GetPosition();

    isoList_->RemoveAll();
    isoListKeys_.clear();
    isoListGaussIdx_.clear();
    isoListAutoMatch_.clear();
    isoListDbEnergy_.clear();

    if (isoHistName_.empty()) return;

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));

    struct Row {
        std::string key;
        double      energy;
        std::string label;         // effective label (per-Gaussian if set, else entry)
        std::string classification;
        int         gaussIdx;      // -1 = whole entry; 0..n-1 = specific Gaussian
    };
    std::vector<Row> rows;
    for (const auto& kv : fitdb.GetEntries()) {
        if (kv.first.empty() || kv.first[0] == '_') continue;
        const FitEntry& fe = kv.second;

        int npar = (int)fe.params.size();
        if (npar >= 5 && (npar - 2) % 3 == 0) {
            int n = (npar - 2) / 3;
            for (int i = 0; i < n; i++) {
                double E        = fe.params[3*i + 1];
                std::string lbl = fe.PeakLabel(i);
                std::string cls = fe.PeakClass(i);
                if (!filterLabel.empty() && filterLabel != "All" && lbl != filterLabel) continue;
                rows.push_back({kv.first, E, lbl, cls, i});
            }
        } else {
            double E = 0.0;
            try { E = std::stod(kv.first.substr(0, kv.first.find('_'))); } catch (...) {}
            if (!filterLabel.empty() && filterLabel != "All" && fe.label != filterLabel) continue;
            rows.push_back({kv.first, E, fe.label, fe.classification, -1});
        }
    }

    std::sort(rows.begin(), rows.end(), [](const Row& a, const Row& b){
        return a.energy < b.energy;
    });

    // ── Unique auto-matching for unlabeled rows ───────────────────────────────
    // 1. Claim DB lines for already-labeled rows (prevents unlabeled stealing them).
    // 2. For remaining unlabeled rows, assign the closest unclaimed DB line;
    //    process in distance-first order so exact matches are claimed before
    //    nearby peaks can steal them.
    std::set<std::string> claimedDbLines;

    if (dbLoaded_) {
        for (const auto& r : rows) {
            if (r.label.empty()) continue;
            double best = std::numeric_limits<double>::max();
            std::string bestKey;
            for (const auto& gl : db_.db) {
                if (gl.isotope != r.label) continue;
                double d = std::abs(gl.energy - r.energy);
                if (d < best) { best = d; bestKey = gl.isotope + Form("_%.1f", gl.energy); }
            }
            if (!bestKey.empty()) claimedDbLines.insert(bestKey);
        }
    }

    // Collect unlabeled rows ordered by distance to nearest DB line.
    // No window cutoff here — the display column always shows the nearest line
    // regardless of the match threshold.  Distance-first ordering ensures that
    // exact-match rows claim their DB lines before nearby rows steal them.
    std::vector<std::pair<double,size_t>> unlabeledOrder; // {minDist, rowIndex}
    if (dbLoaded_) {
        for (size_t i = 0; i < rows.size(); i++) {
            if (!rows[i].label.empty() || rows[i].energy <= 0) continue;
            double minD = std::numeric_limits<double>::max();
            for (const auto& gl : db_.db)
                minD = std::min(minD, std::abs(gl.energy - rows[i].energy));
            unlabeledOrder.push_back({minD, i});
        }
        std::sort(unlabeledOrder.begin(), unlabeledOrder.end(),
                  [](const auto& a, const auto& b){ return a.first < b.first; });
    }

    // For each unlabeled row, assign the closest unclaimed DB line — no distance
    // cutoff so the display always shows a suggestion even for distant matches.
    std::vector<std::string> autoMatchVec(rows.size());
    std::vector<double>      dbEnergyVec(rows.size(), 0.0);
    for (const auto& [dist, i] : unlabeledOrder) {
        double bestD  = std::numeric_limits<double>::max();
        std::string bestIso; double bestE = 0.0;
        for (const auto& gl : db_.db) {
            double d = std::abs(gl.energy - rows[i].energy);
            std::string k = gl.isotope + Form("_%.1f", gl.energy);
            if (claimedDbLines.count(k)) continue;
            if (d < bestD) { bestD = d; bestIso = gl.isotope; bestE = gl.energy; }
        }
        if (!bestIso.empty()) {
            autoMatchVec[i] = bestIso;
            dbEnergyVec[i]  = bestE;
            claimedDbLines.insert(bestIso + Form("_%.1f", bestE));
        }
    }

    int id = 1;
    for (size_t idx = 0; idx < rows.size(); idx++) {
        const auto& r = rows[idx];

        std::string displayLabel;
        double      dbEnergy = 0.0;
        double      deltaE   = 0.0;

        if (!r.label.empty()) {
            displayLabel = r.label;
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
            if (dbEnergyVec[idx] > 0) {
                dbEnergy = dbEnergyVec[idx];
                deltaE   = r.energy - dbEnergy;
            }
        }

        std::string dbCol;
        if (dbEnergy > 0) {
            std::string isoName = !r.label.empty() ? r.label : autoMatchVec[idx];
            dbCol = Form("%-10s %8.2f(%+.2f)", isoName.c_str(), dbEnergy, deltaE);
        } else {
            dbCol = "---";
        }

        std::string entry = Form("%-16s  %8.2f    %s",
                                 displayLabel.c_str(), r.energy, dbCol.c_str());
        isoList_->AddEntry(entry.c_str(), id++);
        isoListKeys_.push_back(r.key);
        isoListGaussIdx_.push_back(r.gaussIdx);
        isoListAutoMatch_.push_back(autoMatchVec[idx]);
        isoListDbEnergy_.push_back(dbEnergy);
    }
    isoList_->MapSubwindows(); isoList_->Layout();
    if (savedScroll > 0) {
        if (TGScrollBar* vsb = isoList_->GetVScrollbar())
            vsb->SetPosition(savedScroll);
    }
}

void GammaFitGUI::RefreshIsoDisplay()
{
    if (isoHistName_.empty()) return;

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));
    LoadLabelClassMap(fitdb);

    std::set<std::string> labels;
    for (const auto& kv : fitdb.GetEntries()) {
        if (kv.first.empty() || kv.first[0] == '_') continue;
        if (!kv.second.label.empty()) labels.insert(kv.second.label);
        for (const auto& pl : kv.second.peakLabels)
            if (!pl.empty()) labels.insert(pl);
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
    isoLabelDecayPeakSel_ = id;  // keep saved selection in sync
    const std::string& key = isoListKeys_[id - 1];
    int gaussIdx = ((size_t)(id-1) < isoListGaussIdx_.size()) ? isoListGaussIdx_[id-1] : -1;
    const std::string autoMatch = ((size_t)(id-1) < isoListAutoMatch_.size())
                                  ? isoListAutoMatch_[id - 1] : "";

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));
    const auto& entries = fitdb.GetEntries();
    auto it = entries.find(key);
    if (it == entries.end()) return;

    // Effective label for this row: per-Gaussian if available, else entry label
    const std::string curLabel = it->second.PeakLabel(gaussIdx);

    if (isoLabelCombo_) {
        isoLabelCombo_->RemoveAll();
        isoLabelCombo_->AddEntry("(none)", 1);
        int cid = 2;
        // curLabel is already per-Gaussian or entry-level

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
        bool inDB = false;
        if (dbLoaded_ && !curLabel.empty())
            for (const auto& gl : db_.db)
                if (gl.isotope == curLabel) { inDB = true; break; }
        isoCustomLabelEntry_->SetText((!curLabel.empty() && !inDB) ? curLabel.c_str() : "");
    }

    if (isoClassCombo_) {
        // Class: per-Gaussian if available, else from labelClassMap_ / entry classification
        std::string cls = it->second.PeakClass(gaussIdx);
        if (cls.empty() && !curLabel.empty() && labelClassMap_.count(curLabel))
            cls = labelClassMap_.at(curLabel);
        // Auto-suggest X-ray for peaks below 100 keV
        int npar = (int)it->second.params.size();
        double peakE = 0.0;
        if (gaussIdx >= 0 && npar >= 5 && (npar-2)%3==0 && gaussIdx < (npar-2)/3)
            peakE = it->second.params[3*gaussIdx + 1];
        else if (npar >= 2)
            peakE = it->second.params[1];
        if (cls.empty() && peakE > 0 && peakE < 100.0) cls = "X-ray";
        isoClassCombo_->Select(ClassToComboIndex(cls), kFALSE);
    }
    if (isoCustomEntry_) {
        std::string cls = it->second.PeakClass(gaussIdx);
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
    const std::string& key    = isoListKeys_[sel - 1];
    int                gaussIdx = ((size_t)(sel-1) < isoListGaussIdx_.size())
                                  ? isoListGaussIdx_[sel-1] : -1;

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));
    const auto& entries = fitdb.GetEntries();
    auto it = entries.find(key);
    if (it == entries.end()) { AppendLog("Cache entry not found."); return; }

    FitEntry e = it->second;

    // Determine new label
    std::string newLabel;
    {
        std::string customLbl = isoCustomLabelEntry_ ? isoCustomLabelEntry_->GetText() : "";
        if (!customLbl.empty()) {
            newLabel = customLbl;
        } else if (isoLabelCombo_) {
            TGLBEntry* le = isoLabelCombo_->GetSelectedEntry();
            std::string lbl = le ? le->GetTitle() : "";
            newLabel = (lbl != "(none)") ? lbl : "";
        }
    }

    // Apply label: per-Gaussian for wide entries, entry-level for tight/single
    if (gaussIdx >= 0) {
        int npar = (int)e.params.size();
        int n    = (npar >= 5 && (npar-2)%3==0) ? (npar-2)/3 : 1;
        if ((int)e.peakLabels.size() < n)          e.peakLabels.resize(n);
        if ((int)e.peakClassifications.size() < n) e.peakClassifications.resize(n);
        e.peakLabels[gaussIdx] = newLabel;
    } else {
        e.label = newLabel;
        // Whole-entry assignment clears per-Gaussian labels (they'd be stale)
        e.peakLabels.clear();
        e.peakClassifications.clear();
    }

    // Auto-classify
    const std::string& appliedLabel = (gaussIdx >= 0) ? e.peakLabels[gaussIdx] : e.label;
    if (!appliedLabel.empty()) {
        std::string cls = AutoClassFromParent(appliedLabel);
        if (cls.empty()) {
            auto mapIt = labelClassMap_.find(appliedLabel);
            if (mapIt != labelClassMap_.end()) cls = mapIt->second;
        }
        if (!cls.empty()) {
            if (gaussIdx >= 0)
                e.peakClassifications[gaussIdx] = cls;
            else
                e.classification = cls;
            labelClassMap_[appliedLabel] = cls;
        }
    }

    fitdb.ForceStore(key, e);

    // Propagate entry-level classification to peers with same entry label
    if (gaussIdx < 0 && !e.classification.empty() && !e.label.empty()) {
        for (const auto& kv : fitdb.GetEntries()) {
            if (kv.first == key || kv.first.empty() || kv.first[0] == '_') continue;
            if (kv.second.label != e.label) continue;
            FitEntry pe = kv.second;
            pe.classification = e.classification;
            fitdb.ForceStore(kv.first, pe);
        }
        SaveLabelClassMap(fitdb);
    } else if (gaussIdx >= 0 && !appliedLabel.empty()) {
        SaveLabelClassMap(fitdb);
    }

    mkdir(kCacheDir, 0755);
    fitdb.Save(CacheFileFor(isoHistName_));
    std::string lbl = (gaussIdx >= 0) ? e.peakLabels[gaussIdx] : e.label;
    std::string cls = (gaussIdx >= 0) ? e.peakClassifications[gaussIdx] : e.classification;
    AppendLog("Updated: " + key + (gaussIdx >= 0 ? Form("[G%d]", gaussIdx) : "") +
              "  label=" + lbl + "  class=" + cls);
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
    const std::string& key    = isoListKeys_[sel - 1];
    int                gaussIdx = ((size_t)(sel-1) < isoListGaussIdx_.size())
                                  ? isoListGaussIdx_[sel-1] : -1;
    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));
    auto it = fitdb.GetEntries().find(key);
    if (it == fitdb.GetEntries().end()) { AppendLog("Cache entry not found."); return; }
    FitEntry e = it->second;
    if (gaussIdx >= 0 && gaussIdx < (int)e.peakLabels.size()) {
        e.peakLabels[gaussIdx]          = "";
        if (gaussIdx < (int)e.peakClassifications.size())
            e.peakClassifications[gaussIdx] = "";
    } else {
        e.label          = "";
        e.classification = "";
        e.peakLabels.clear();
        e.peakClassifications.clear();
    }
    fitdb.ForceStore(key, e);
    mkdir(kCacheDir, 0755);
    fitdb.Save(CacheFileFor(isoHistName_));
    if (isoLabelCombo_)       isoLabelCombo_->Select(1, kFALSE);
    if (isoCustomLabelEntry_) isoCustomLabelEntry_->SetText("");
    if (isoClassCombo_)       isoClassCombo_->Select(1, kFALSE);
    if (isoCustomEntry_)      isoCustomEntry_->SetText("");
    AppendLog("Cleared label/class for: " + key +
              (gaussIdx >= 0 ? Form("[G%d]", gaussIdx) : ""));
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
        bool hasAny = !kv.second.label.empty() || !kv.second.classification.empty()
                      || !kv.second.peakLabels.empty() || !kv.second.peakClassifications.empty();
        if (!hasAny) continue;
        FitEntry e = kv.second;
        e.label = ""; e.classification = "";
        e.peakLabels.clear(); e.peakClassifications.clear();
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

    // Each "match unit" is one claimable peak.  For multi-Gaussian entries:
    //   - if all Gaussians span ≤ 4 keV → one unit (highest-amplitude Gaussian)
    //   - if span > 4 keV               → one unit per Gaussian (each gets its own label)
    // gaussIdx == -1 means the label goes to the entry .label field.
    // gaussIdx >= 0 means the label goes to entry .peakLabels[gaussIdx].
    struct MatchUnit {
        std::string entryKey;
        int         gaussIdx;  // -1 = whole entry, 0..n-1 = specific Gaussian
        double      energy;
        bool        hasLabel;
    };

    std::vector<MatchUnit> units;
    for (const auto& kv : fitdb.GetEntries()) {
        if (kv.first.empty() || kv.first[0] == '_') continue;
        const FitEntry& fe = kv.second;
        int npar = (int)fe.params.size();

        // Extract per-Gaussian energies and amplitudes
        std::vector<std::pair<double,double>> gaussians;  // {energy, amplitude}
        if (npar >= 5 && (npar - 2) % 3 == 0) {
            int n = (npar - 2) / 3;
            for (int i = 0; i < n; i++)
                gaussians.push_back({fe.params[3*i + 1], fe.params[3*i]});
        } else if (npar >= 2) {
            double repE = fe.params[1];
            if (repE <= 0)
                try { repE = std::stod(kv.first.substr(0, kv.first.find('_'))); } catch (...) {}
            if (repE > 0) gaussians.push_back({repE, npar >= 1 ? fe.params[0] : 1.0});
        }
        if (gaussians.empty()) continue;

        // Determine spread
        double minE = gaussians[0].first, maxE = gaussians[0].first;
        for (const auto& [e, a] : gaussians) { minE = std::min(minE, e); maxE = std::max(maxE, e); }
        bool tightCluster = (maxE - minE) <= 4.0;

        if (tightCluster || gaussians.size() == 1) {
            // One unit for the whole entry — use highest-amplitude Gaussian as representative
            double repE = gaussians[0].first, bestAmp = gaussians[0].second;
            for (const auto& [e, a] : gaussians) if (a > bestAmp) { bestAmp = a; repE = e; }
            bool hasLbl = !fe.label.empty();
            units.push_back({kv.first, -1, repE, hasLbl});
        } else {
            // One unit per Gaussian — each gets matched and labeled independently.
            // Only per-Gaussian labels count here.  A legacy entry-level label (set
            // by old single-unit matching) is intentionally ignored so that all
            // Gaussians in the wide entry are re-matched individually.
            int n = (int)gaussians.size();
            for (int i = 0; i < n; i++) {
                bool hasLbl = (i < (int)fe.peakLabels.size() && !fe.peakLabels[i].empty());
                units.push_back({kv.first, i, gaussians[i].first, hasLbl});
            }
        }
    }

    // Threshold from UI (default 10 keV)
    const double matchThresh = (isoMatchThreshEntry_
                                ? isoMatchThreshEntry_->GetNumber() : 10.0);

    // Distance-only unclaimed line lookup — rejects lines beyond matchThresh
    auto findClosestUnclaimed = [&](double energy,
                                    const std::set<std::string>& claimed)
        -> std::pair<std::string, double>
    {
        double bestDist = std::numeric_limits<double>::max();
        std::string bestIso; double bestE = 0.0;
        for (const auto& gl : db_.db) {
            double dE = std::fabs(gl.energy - energy);
            if (dE > matchThresh) continue;
            std::string k = gl.isotope + Form("_%.1f", gl.energy);
            if (claimed.count(k)) continue;
            if (dE < bestDist) { bestDist = dE; bestIso = gl.isotope; bestE = gl.energy; }
        }
        return {bestIso, bestE};
    };

    // Pass 1: claim DB lines for already-labeled units (prevents them being stolen)
    std::set<std::string> claimedLines;
    for (const auto& u : units) {
        if (!u.hasLabel) continue;
        const FitEntry& fe = fitdb.GetEntries().at(u.entryKey);
        std::string lbl = fe.PeakLabel(u.gaussIdx);
        if (lbl.empty()) continue;
        double best = std::numeric_limits<double>::max();
        std::string bestKey;
        for (const auto& gl : db_.db) {
            if (gl.isotope != lbl) continue;
            double d = std::fabs(gl.energy - u.energy);
            if (d < best) { best = d; bestKey = gl.isotope + Form("_%.1f", gl.energy); }
        }
        if (!bestKey.empty()) claimedLines.insert(bestKey);
    }

    // Pass 2: sort unlabeled units within matchThresh by distance (distance-first
    // ordering ensures exact-match units claim their lines before nearby units steal them)
    std::vector<std::pair<double, size_t>> unlabeled;  // {minDist, unitIndex}
    for (size_t i = 0; i < units.size(); i++) {
        if (units[i].hasLabel) continue;
        double minD = std::numeric_limits<double>::max();
        for (const auto& gl : db_.db) {
            double d = std::fabs(gl.energy - units[i].energy);
            if (d <= matchThresh) minD = std::min(minD, d);
        }
        if (minD <= matchThresh)
            unlabeled.push_back({minD, i});
    }
    std::sort(unlabeled.begin(), unlabeled.end(),
              [](const auto& a, const auto& b){ return a.first < b.first; });

    int nMatched = 0;
    for (const auto& [dist, ui] : unlabeled) {
        const MatchUnit& u = units[ui];
        auto [iso, dbE] = findClosestUnclaimed(u.energy, claimedLines);
        if (iso.empty()) continue;

        auto it = fitdb.GetEntries().find(u.entryKey);
        if (it == fitdb.GetEntries().end()) continue;
        FitEntry e = it->second;

        // Auto-classify
        std::string cls = AutoClassFromParent(iso);
        if (cls.empty()) {
            auto mapIt = labelClassMap_.find(iso);
            if (mapIt != labelClassMap_.end()) cls = mapIt->second;
        }

        if (u.gaussIdx < 0) {
            // Tight cluster or single peak — label the whole entry
            e.label = iso;
            if (!cls.empty()) { e.classification = cls; }
        } else {
            // Wide multi-Gaussian — label only this Gaussian
            int npar = (int)e.params.size();
            int n    = (npar >= 5 && (npar-2)%3==0) ? (npar-2)/3 : 1;
            if ((int)e.peakLabels.size() < n)          e.peakLabels.resize(n);
            if ((int)e.peakClassifications.size() < n) e.peakClassifications.resize(n);
            e.peakLabels[u.gaussIdx] = iso;
            if (!cls.empty()) e.peakClassifications[u.gaussIdx] = cls;
        }
        if (!cls.empty()) labelClassMap_[iso] = cls;

        fitdb.ForceStore(u.entryKey, e);
        claimedLines.insert(iso + Form("_%.1f", dbE));
        ++nMatched;
    }

    SaveLabelClassMap(fitdb);
    mkdir(kCacheDir, 0755);
    fitdb.Save(CacheFileFor(isoHistName_));
    AppendLog("Auto-matched " + std::to_string(nMatched) + " peaks (≤4 keV clusters = one unit).");
    RefreshIsoDisplay();
}

void GammaFitGUI::OnIsoApplyAutoMatches()
{
    // Match labels and save to cache
    OnIsoAutoMatchAll();

    if (isoHistName_.empty() || !canvas_) return;

    // Switch to the matched histogram if not already showing it
    if (currentHist_ != isoHistName_) {
        currentHist_ = isoHistName_;
        bool owned = false;
        if (rawHistOwned_ && rawHist_) { delete rawHist_; rawHist_ = nullptr; }
        rawHist_      = LoadHistFromFile(isoHistName_, owned);
        rawHistOwned_ = owned;
    }
    if (!rawHist_) return;

    TCanvas* c = canvas_->GetCanvas();
    if (histViewCombo_) OnHistViewChanged(histViewCombo_->GetSelected());
    else                RedrawCurrent();
    OverlayFitPeaks(isoHistName_, c);
    AppendLog("Histogram redrawn with isotope labels for " + isoHistName_);
}

void GammaFitGUI::OnIsoSetParent()
{
    if (isoHistName_.empty()) { AppendLog("Refresh cache first."); return; }

    isoParentIsotope_ = isoParentNameEntry_ ? std::string(isoParentNameEntry_->GetText()) : "";
    isoParentZval_    = isoParentZEntry_    ? (int)isoParentZEntry_->GetNumber()  : 0;
    isoParentNval_    = isoParentNEntry_    ? (int)isoParentNEntry_->GetNumber()  : 0;

    // Parent is always "Parent" class
    if (!isoParentIsotope_.empty()) labelClassMap_[isoParentIsotope_] = "Parent";

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));

    // Persist parent info
    FitEntry pe;
    pe.key     = kParentInfoKey;
    pe.label   = Form("name=%s;Z=%d;N=%d",
                      isoParentIsotope_.c_str(), isoParentZval_, isoParentNval_);
    pe.chi2ndf = 0; pe.residualRMS = 0; pe.maxPull = 0;
    fitdb.ForceStore(kParentInfoKey, pe);

    // Auto-classify every labeled peak using the new parent chain
    int nClassified = 0;
    for (const auto& kv : fitdb.GetEntries()) {
        if (kv.first.empty() || kv.first[0] == '_') continue;
        const FitEntry& fe = kv.second;
        bool changed = false;
        FitEntry e = fe;

        // Entry-level label
        if (!e.label.empty()) {
            std::string cls = AutoClassFromParent(e.label);
            if (!cls.empty() && cls != e.classification) {
                e.classification = cls;
                labelClassMap_[e.label] = cls;
                changed = true;
            }
        }
        // Per-Gaussian labels (wide multi-peak entries)
        for (int i = 0; i < (int)e.peakLabels.size(); i++) {
            if (e.peakLabels[i].empty()) continue;
            std::string cls = AutoClassFromParent(e.peakLabels[i]);
            if (cls.empty()) continue;
            if ((int)e.peakClassifications.size() <= i)
                e.peakClassifications.resize(i + 1);
            if (cls != e.peakClassifications[i]) {
                e.peakClassifications[i] = cls;
                labelClassMap_[e.peakLabels[i]] = cls;
                changed = true;
            }
        }
        if (changed) { fitdb.ForceStore(kv.first, e); ++nClassified; }
    }

    SaveLabelClassMap(fitdb);
    mkdir(kCacheDir, 0755);
    fitdb.Save(CacheFileFor(isoHistName_));

    const char* sym = ElementSymbol(isoParentZval_);
    int A = isoParentZval_ + isoParentNval_;
    AppendLog(Form("Parent set: %s  Z=%d  N=%d  A=%d  (%s)  — auto-classified %d peaks",
                   isoParentIsotope_.c_str(), isoParentZval_, isoParentNval_, A, sym,
                   nClassified));
    RefreshIsoDisplay();
}

// ─────────────────────────────────────────────────────────────────────────────
// Set Label & Decay dialog
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnIsoSetLabelDecay()
{
    if (!dbLoaded_) { AppendLog("Load isotope DB first."); return; }
    if (isoHistName_.empty()) { AppendLog("Refresh cache first."); return; }

    // If already open, raise to front
    if (isoLabelDecayDlg_) {
        isoLabelDecayDlg_->RaiseWindow();
        return;
    }

    // Save the currently selected peak — the modeless dialog steals focus and
    // causes isoList_->GetSelected() to return -1 when Apply is clicked.
    isoLabelDecayPeakSel_ = isoList_ ? isoList_->GetSelected() : -1;

    // Create modeless popup — no WaitFor, stays open across multiple Apply clicks
    isoLabelDecayDlg_ = new TGTransientFrame(gClient->GetRoot(), this, 390, 530);
    isoLabelDecayDlg_->SetWindowName("Set Isotope & Decay Type");
    isoLabelDecayDlg_->SetCleanup(kDeepCleanup);
    // Null member pointers when closed via X or our Close button
    isoLabelDecayDlg_->Connect("CloseWindow()", "GammaFitGUI", this,
                               "OnIsoLabelDecayDlgClosed()");

    // Search row
    {
        TGHorizontalFrame* sRow = new TGHorizontalFrame(isoLabelDecayDlg_);
        isoLabelDecayDlg_->AddFrame(sRow,
            new TGLayoutHints(kLHintsExpandX, 8, 8, 8, 2));
        sRow->AddFrame(new TGLabel(sRow, "Search:"),
            new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        isoLabelDecaySearch_ = new TGTextEntry(sRow, "");
        sRow->AddFrame(isoLabelDecaySearch_, new TGLayoutHints(kLHintsExpandX, 0, 4, 0, 0));
        TGTextButton* fBtn = new TGTextButton(sRow, "Filter");
        sRow->AddFrame(fBtn, new TGLayoutHints(kLHintsLeft));
        fBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoLabelDecaySearch()");
    }

    // Isotope list
    isoLabelDecayList_ = new TGListBox(isoLabelDecayDlg_, 1100);
    isoLabelDecayList_->Resize(370, 310);
    isoLabelDecayDlg_->AddFrame(isoLabelDecayList_,
        new TGLayoutHints(kLHintsExpandX, 8, 8, 2, 4));

    // Populate with unique sorted DB names (append [class] when known)
    {
        std::set<std::string> seen;
        std::vector<std::string> names;
        for (const auto& gl : db_.db)
            if (seen.insert(gl.isotope).second)
                names.push_back(gl.isotope);
        std::sort(names.begin(), names.end());
        int idx = 1;
        for (const auto& nm : names) {
            std::string display = nm;
            auto mapIt = labelClassMap_.find(nm);
            if (mapIt != labelClassMap_.end() && !mapIt->second.empty())
                display += "  [" + mapIt->second + "]";
            isoLabelDecayList_->AddEntry(display.c_str(), idx++);
        }
        isoLabelDecayList_->MapSubwindows(); isoLabelDecayList_->Layout();
    }

    // Decay type row
    {
        TGHorizontalFrame* dtRow = new TGHorizontalFrame(isoLabelDecayDlg_);
        isoLabelDecayDlg_->AddFrame(dtRow,
            new TGLayoutHints(kLHintsExpandX, 8, 8, 4, 4));
        dtRow->AddFrame(new TGLabel(dtRow, "Decay type:"),
            new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));
        isoDecayTypeCombo_ = new TGComboBox(dtRow, 1101);
        isoDecayTypeCombo_->AddEntry("(none)",                1);
        isoDecayTypeCombo_->AddEntry("Daughter",              2);
        isoDecayTypeCombo_->AddEntry("Granddaughter",         3);
        isoDecayTypeCombo_->AddEntry("Beta-n Daughter",       4);
        isoDecayTypeCombo_->AddEntry("Beta-2n Daughter",      5);
        isoDecayTypeCombo_->AddEntry("Beta-n Granddaughter",  6);
        isoDecayTypeCombo_->AddEntry("Beta-2n Granddaughter", 7);
        isoDecayTypeCombo_->AddEntry("Background",            8);
        isoDecayTypeCombo_->AddEntry("X-ray",                 9);
        isoDecayTypeCombo_->Select(1, kFALSE);
        isoDecayTypeCombo_->Resize(180, 22);
        dtRow->AddFrame(isoDecayTypeCombo_, new TGLayoutHints(kLHintsLeft));
    }

    // Buttons row — Apply keeps dialog open; Close dismisses it
    {
        TGHorizontalFrame* bRow = new TGHorizontalFrame(isoLabelDecayDlg_);
        isoLabelDecayDlg_->AddFrame(bRow,
            new TGLayoutHints(kLHintsCenterX, 8, 8, 6, 8));
        TGTextButton* applyBtn = new TGTextButton(bRow, "  Apply  ");
        bRow->AddFrame(applyBtn, new TGLayoutHints(kLHintsLeft, 0, 12, 0, 0));
        applyBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoLabelDecayApply()");
        TGTextButton* closeBtn = new TGTextButton(bRow, "  Close  ");
        bRow->AddFrame(closeBtn, new TGLayoutHints(kLHintsLeft));
        closeBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoLabelDecayClose()");
    }

    isoLabelDecayDlg_->MapSubwindows();
    isoLabelDecayDlg_->Layout();
    isoLabelDecayDlg_->MapWindow();
    isoLabelDecayDlg_->CenterOnParent();
    // Dialog is now open and modeless — control returns immediately
}

void GammaFitGUI::OnIsoLabelDecayApply()
{
    if (!isoLabelDecayList_ || !isoDecayTypeCombo_) return;
    if (isoHistName_.empty()) return;

    TGLBEntry* le = isoLabelDecayList_->GetSelectedEntry();
    std::string isotope = le ? le->GetTitle() : "";
    // Strip " [class]" annotation added for display
    auto bracketPos = isotope.find("  [");
    if (bracketPos != std::string::npos) isotope = isotope.substr(0, bracketPos);
    if (isotope.empty()) { AppendLog("Select an isotope from the list first."); return; }

    TGLBEntry* de = isoDecayTypeCombo_->GetSelectedEntry();
    std::string dt = de ? de->GetTitle() : "";
    std::string decayType = (dt == "(none)") ? "" : dt;

    if (decayType.empty()) {
        AppendLog("Select a decay type first."); return;
    }

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));

    // Assign class to this isotope in the global map
    labelClassMap_[isotope] = decayType;

    // Propagate to ALL peaks (entry-level and per-Gaussian) labeled with this isotope
    int nUpdated = 0;
    for (const auto& kv : fitdb.GetEntries()) {
        if (kv.first.empty() || kv.first[0] == '_') continue;
        FitEntry e = kv.second;
        bool changed = false;

        if (e.label == isotope) {
            e.classification = decayType;
            changed = true;
        }
        for (int i = 0; i < (int)e.peakLabels.size(); i++) {
            if (e.peakLabels[i] == isotope) {
                if ((int)e.peakClassifications.size() <= i)
                    e.peakClassifications.resize(i + 1);
                e.peakClassifications[i] = decayType;
                changed = true;
            }
        }
        if (changed) {
            fitdb.ForceStore(kv.first, e);
            ++nUpdated;
        }
    }

    SaveLabelClassMap(fitdb);
    mkdir(kCacheDir, 0755);
    fitdb.Save(CacheFileFor(isoHistName_));
    AppendLog("Assigned decay type '" + decayType + "' to isotope '" + isotope +
              "' — updated " + std::to_string(nUpdated) + " peaks.");
    RefreshIsoDisplay();
    // Dialog stays open for the next assignment
}

void GammaFitGUI::OnIsoLabelDecayClose()
{
    // Triggered by the Close button — let CloseWindow handle deletion
    if (isoLabelDecayDlg_) isoLabelDecayDlg_->CloseWindow();
}

void GammaFitGUI::OnIsoLabelDecayDlgClosed()
{
    // Triggered by the CloseWindow signal (X button or our Close button) —
    // null member pointers before the widget tree is destroyed
    isoLabelDecayDlg_    = nullptr;
    isoLabelDecayList_   = nullptr;
    isoDecayTypeCombo_   = nullptr;
    isoLabelDecaySearch_ = nullptr;
}

void GammaFitGUI::OnIsoLabelDecaySearch()
{
    if (!isoLabelDecayList_ || !isoLabelDecaySearch_) return;
    std::string filter = isoLabelDecaySearch_->GetText();
    std::string lf = filter;
    std::transform(lf.begin(), lf.end(), lf.begin(), ::tolower);

    isoLabelDecayList_->RemoveAll();
    std::set<std::string> seen;
    std::vector<std::string> names;
    for (const auto& gl : db_.db)
        if (seen.insert(gl.isotope).second) {
            std::string nl = gl.isotope;
            std::transform(nl.begin(), nl.end(), nl.begin(), ::tolower);
            if (lf.empty() || nl.find(lf) != std::string::npos)
                names.push_back(gl.isotope);
        }
    std::sort(names.begin(), names.end());
    int idx = 1;
    for (const auto& nm : names) {
        std::string display = nm;
        auto mapIt = labelClassMap_.find(nm);
        if (mapIt != labelClassMap_.end() && !mapIt->second.empty())
            display += "  [" + mapIt->second + "]";
        isoLabelDecayList_->AddEntry(display.c_str(), idx++);
    }
    isoLabelDecayList_->MapSubwindows();
    isoLabelDecayList_->Layout();
}

// ─────────────────────────────────────────────────────────────────────────────
// Auto-classify from parent decay chain
// ─────────────────────────────────────────────────────────────────────────────

// Parse "44Cl", "Bi214", "K40" etc. into (A, Z, N). Returns false if unrecognised.
static bool ParseIsotopeName(const std::string& name, int& A, int& Z, int& N)
{
    if (name.empty()) return false;

    // Leading digits: "44Cl", "43Ar", "44S"
    {
        size_t i = 0;
        while (i < name.size() && std::isdigit(name[i])) ++i;
        if (i > 0 && i < name.size()) {
            try {
                int a = std::stoi(name.substr(0, i));
                int z = SymbolToZ(name.substr(i));
                if (z > 0 && a > z) { A=a; Z=z; N=a-z; return true; }
            } catch (...) {}
        }
    }

    // Trailing digits: "Bi214", "K40", "Ti208"
    {
        size_t j = name.size();
        while (j > 0 && std::isdigit(name[j-1])) --j;
        if (j < name.size() && j > 0) {
            try {
                int a = std::stoi(name.substr(j));
                int z = SymbolToZ(name.substr(0, j));
                if (z > 0 && a > z) { A=a; Z=z; N=a-z; return true; }
            } catch (...) {}
        }
    }
    return false;
}

std::string GammaFitGUI::AutoClassFromParent(const std::string& label)
{
    if (isoParentZval_ <= 0 || isoParentNval_ <= 0) return "";
    int A, Z, N;
    if (!ParseIsotopeName(label, A, Z, N)) return "";

    int dZ = Z - isoParentZval_;
    int dN = N - isoParentNval_;

    // Pure beta-minus chain: each step dZ=+1, dN=-1
    if (dZ >= 0 && dZ <= 4 && dZ == -dN) {
        switch (dZ) {
        case 0: return "Parent";
        case 1: return "Daughter";
        case 2: return "Granddaughter";
        case 3: return "Beta-n Granddaughter";
        case 4: return "Beta-2n Granddaughter";
        }
    }
    return "";
}

// ─────────────────────────────────────────────────────────────────────────────
// Preview selected peak on canvas
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnIsoPeakPreview()
{
    if (!isoList_) return;
    Int_t sel = isoList_->GetSelected();
    if (sel < 1 || (size_t)sel > isoListKeys_.size()) {
        AppendLog("Select a peak first."); return;
    }
    if (isoHistName_.empty()) return;

    const std::string& key = isoListKeys_[sel - 1];
    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));
    auto it = fitdb.GetEntries().find(key);
    if (it == fitdb.GetEntries().end()) return;

    // Find the highest-amplitude peak energy in this entry
    const FitEntry& fe = it->second;
    double peakE = 0.0;
    int npar = (int)fe.params.size();
    if (npar >= 5 && (npar - 2) % 3 == 0) {
        int n = (npar - 2) / 3;
        double bestAmp = -1.0;
        for (int i = 0; i < n; i++) {
            if (fe.params[3*i] > bestAmp) {
                bestAmp = fe.params[3*i];
                peakE   = fe.params[3*i + 1];
            }
        }
    }
    if (peakE <= 0.0 && npar >= 2) peakE = fe.params[1];
    if (peakE <= 0.0) {
        try { peakE = std::stod(key.substr(0, key.find('_'))); } catch (...) {}
    }
    if (peakE <= 0.0) { AppendLog("Cannot determine peak energy."); return; }

    // Switch histogram if needed
    if (currentHist_ != isoHistName_) {
        currentHist_ = isoHistName_;
        bool owned = false;
        if (rawHistOwned_ && rawHist_) { delete rawHist_; rawHist_ = nullptr; }
        rawHist_      = LoadHistFromFile(isoHistName_, owned);
        rawHistOwned_ = owned;
    }
    if (!rawHist_) { AppendLog("Could not load histogram: " + isoHistName_); return; }

    // Zoom to ±6σ (at least 20 keV) around the peak
    double sigma = res_.Sigma(peakE);
    double halfW = std::max(6.0 * sigma, 20.0);
    viewXmin_ = std::max(peakE - halfW, rawHist_->GetXaxis()->GetXmin());
    viewXmax_ = std::min(peakE + halfW, rawHist_->GetXaxis()->GetXmax());

    if (histViewCombo_) OnHistViewChanged(histViewCombo_->GetSelected());
    else                RedrawCurrent();

    OverlayFitPeaks(currentHist_, canvas_->GetCanvas());
    canvas_->GetCanvas()->Modified();
    canvas_->GetCanvas()->Update();
    AppendLog(Form("Preview: %.2f keV  [%s]", peakE, key.c_str()));
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
        std::string nucInfo;  // e.g. "^{90}Kr  (Z=36, N=54)"
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
        // Compute nuclear identity from beta-minus chain
        int ZZ = 0, NN = 0;
        if (ClassZN(cls, isoParentZval_, isoParentNval_, ZZ, NN) && ZZ > 0 && NN > 0) {
            int A = ZZ + NN;
            nd.nucInfo = Form("^{%d}%s  (Z=%d, N=%d)", A, ElementSymbol(ZZ), ZZ, NN);
        }
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
        // Height of each node: class header + optional nucInfo line + isotope lines
        for (auto* nd : ndlist)
            nd->h = ((int)nd->isos.size() + 1 + (nd->nucInfo.empty() ? 0 : 1)) * lineH
                    + 2 * padV;
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

        // Nuclear identity line (^{A}El, Z, N) if parent is set
        if (!nd.nucInfo.empty()) {
            ty -= lineH;
            tx.SetTextSize(0.026); tx.SetTextColor(kRed + 1);
            tx.DrawLatex(nd.xc, ty, nd.nucInfo.c_str());
        }

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

        // β⁻ label above the arrow midpoint
        TLatex* bmlbl = new TLatex((ax1 + ax2) / 2.0, ay + 0.018, "#beta^{-}");
        bmlbl->SetNDC(kTRUE);
        bmlbl->SetTextSize(0.022);
        bmlbl->SetTextAlign(22);
        bmlbl->SetTextColor(kMagenta + 1);
        bmlbl->Draw();
    }

    // Title
    tx.SetTextSize(0.032); tx.SetTextColor(kBlack); tx.SetTextAlign(22);
    if (!isoParentIsotope_.empty() && isoParentZval_ > 0) {
        int A = isoParentZval_ + isoParentNval_;
        tx.DrawLatex(0.50, 0.97,
            Form("Decay Schematic: %s   |   Parent: %s (^{%d}%s, Z=%d, N=%d)",
                 isoHistName_.c_str(),
                 isoParentIsotope_.c_str(), A,
                 ElementSymbol(isoParentZval_),
                 isoParentZval_, isoParentNval_));
    } else {
        tx.DrawLatex(0.50, 0.97,
            Form("Decay Schematic: %s   |   Parent not set",
                 isoHistName_.c_str()));
    }

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
    AppendLog("Matched: " + key + "  ->  " + gl.isotope +
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
