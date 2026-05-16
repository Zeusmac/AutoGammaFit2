#include "GammaFitGUI.h"
#include "GammaFitGUI_shared.h"

#include "TGTextEntry.h"
#include "TGFileDialog.h"
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
#include <fstream>
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

    // ── Active histogram indicator ────────────────────────────────────────────
    {
        TGHorizontalFrame* histRow = new TGHorizontalFrame(p);
        p->AddFrame(histRow, new TGLayoutHints(kLHintsExpandX, 4, 4, 0, 2));
        histRow->AddFrame(new TGLabel(histRow, "Histogram:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        isoHistLabel_ = new TGLabel(histRow, "(none — click Refresh)");
        histRow->AddFrame(isoHistLabel_, new TGLayoutHints(kLHintsLeft));
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
        parentGrp->AddFrame(setParentBtn, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        setParentBtn->Connect("Clicked()", "GammaFitGUI", this, "OnIsoSetParent()");
        setParentBtn->SetToolTipText(
            "Save parent isotope name, Z, and N to the cache.\n"
            "Daughter/granddaughter Z and N are auto-calculated\n"
            "via beta-minus decay steps when drawing the schematic.");
    }
    {
        // AME mass table for Q-value computation
        TGHorizontalFrame* ameRow = new TGHorizontalFrame(parentGrp);
        parentGrp->AddFrame(ameRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        TGTextButton* ameBtn = new TGTextButton(ameRow, "Load AME Table");
        ameRow->AddFrame(ameBtn, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));
        ameBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLoadAMETable()");
        ameBtn->SetToolTipText(
            "Load the NNDC/IAEA AME2020 mass table (mass_1.mas) to compute\n"
            "Q values for each decay step shown in the decay schematic.\n"
            "Download from: https://www-nds.iaea.org/amdc/ame2020/mass_1.mas");
        ameLbl_ = new TGLabel(ameRow, "(not loaded)");
        ameRow->AddFrame(ameLbl_, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 0, 0, 0, 0));
    }
    {
        // NUBASE2020 for half-lives, branching ratios
        TGHorizontalFrame* nbRow = new TGHorizontalFrame(parentGrp);
        parentGrp->AddFrame(nbRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 0, 4));
        TGTextButton* nbBtn = new TGTextButton(nbRow, "Load NUBASE Table");
        nbRow->AddFrame(nbBtn, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));
        nbBtn->Connect("Clicked()", "GammaFitGUI", this, "OnLoadNubaseTable()");
        nbBtn->SetToolTipText(
            "Load the IAEA NUBASE2020 table (nubase_2020.txt) for half-lives\n"
            "and beta-decay branching ratios shown in the decay schematic.\n"
            "Download from: https://www-nds.iaea.org/amdc/ame2020/nubase_2020.txt");
        nubaseLbl_ = new TGLabel(nbRow, "(not loaded)");
        nbRow->AddFrame(nubaseLbl_, new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 0, 0, 0, 0));
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
    if (isoHistLabel_) isoHistLabel_->SetText(isoHistName_.c_str());
    RefreshIsoDisplay();
    AppendLog("Isotopes: loaded cache for " + isoHistName_ +
              "  (" + std::to_string(isoListKeys_.size()) + " peaks)");
    if (schematicDrawn_)
        DrawDecaySchematic(canvas_->GetCanvas());
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

// ─────────────────────────────────────────────────────────────────────────────
// AME mass-table loader  (NNDC/IAEA AME2020 mass_1.mas format)
// Q_β- = Δ(parent) - Δ(daughter) - n_emitted × Δ(neutron)
// ─────────────────────────────────────────────────────────────────────────────

bool GammaFitGUI::LoadAMETable(const std::string& path)
{
    std::ifstream f(path);
    if (!f.is_open()) return false;

    ameTable_.clear();
    std::string line;
    int nLoaded = 0;

    while (std::getline(f, line)) {
        if (line.empty()) continue;
        // Skip page-header lines that begin with "1"
        if (!line.empty() && line[0] == '1') continue;

        // Tokenise the line
        std::istringstream ss(line);
        std::vector<std::string> tok;
        std::string t;
        while (ss >> t) tok.push_back(t);

        if (tok.size() < 7) continue;

        // Find the element-symbol token: shortest alphabetic string (1-3 chars)
        // that sits between purely-numeric tokens.  It separates (NZ,N,Z,A) from
        // the mass-excess column.
        int symIdx = -1;
        for (int i = 1; i < (int)tok.size() - 1; i++) {
            if (tok[i].size() >= 1 && tok[i].size() <= 3 &&
                std::isalpha((unsigned char)tok[i][0])) {
                symIdx = i;
                break;
            }
        }
        // Need at least three numeric tokens before the symbol (N, Z, A)
        if (symIdx < 3) continue;

        try {
            int N = std::stoi(tok[symIdx - 3]);
            int Z = std::stoi(tok[symIdx - 2]);
            int A = std::stoi(tok[symIdx - 1]);
            if (Z < 0 || N < 0 || A != Z + N) continue;

            // Mass excess is the first token after the symbol.
            // Strip '#' (extrapolated) and '*' (estimated) markers.
            std::string mexStr = tok[symIdx + 1];
            mexStr.erase(std::remove_if(mexStr.begin(), mexStr.end(),
                         [](char c){ return c == '#' || c == '*'; }), mexStr.end());
            if (mexStr.empty()) continue;

            ameTable_[{Z, A}] = std::stod(mexStr);
            ++nLoaded;
        } catch (...) {
            continue;
        }
    }

    ameLoaded_ = (nLoaded > 100);
    return ameLoaded_;
}

void GammaFitGUI::OnLoadAMETable()
{
    TGFileInfo fi;
    static const char* kFileTypes[] = {
        "AME mass files", "*.mas *.txt *.dat",
        "All files",      "*",
        nullptr,          nullptr
    };
    fi.fFileTypes = kFileTypes;
    fi.fIniDir    = StrDup(".");
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if (!fi.fFilename) return;

    std::string path = fi.fFilename;
    bool ok = LoadAMETable(path);
    if (ok) {
        std::string disp = path;
        // Show only filename (strip path)
        auto slash = disp.rfind('/');
        if (slash != std::string::npos) disp = disp.substr(slash + 1);
        if (ameLbl_) ameLbl_->SetText(disp.c_str());
        AppendLog("AME table loaded: " + path + "  (" +
                  std::to_string(ameTable_.size()) + " nuclides)");
    } else {
        if (ameLbl_) ameLbl_->SetText("(load failed)");
        AppendLog("ERROR: Could not parse AME table from " + path);
    }
}

// ── NUBASE2020 parser ─────────────────────────────────────────────────────────
bool GammaFitGUI::LoadNubaseTable(const std::string& path)
{
    std::ifstream in(path);
    if (!in.is_open()) return false;
    nubaseTable_.clear();

    // Unit → seconds conversion
    auto unitToSec = [](const std::string& u) -> double {
        if (u == "ps") return 1e-12;
        if (u == "ns") return 1e-9;
        if (u == "us") return 1e-6;
        if (u == "ms") return 1e-3;
        if (u == "s")  return 1.0;
        if (u == "m")  return 60.0;
        if (u == "h")  return 3600.0;
        if (u == "d")  return 86400.0;
        if (u == "y")  return 3.156e7;
        if (u == "ky") return 3.156e10;
        if (u == "My") return 3.156e13;
        if (u == "Gy") return 3.156e16;
        if (u == "Ty") return 3.156e19;
        if (u == "stbl") return 0.0;    // stable
        return -1.0;                     // unknown
    };

    // Extract branching ratio value after "KEY=" in a decay mode field
    auto parseBR = [](const std::string& field, const std::string& key) -> double {
        auto pos = field.find(key + "=");
        if (pos == std::string::npos) return -1.0;
        pos += key.size() + 1;
        // skip leading spaces
        while (pos < field.size() && field[pos] == ' ') ++pos;
        std::string tok;
        while (pos < field.size() && (std::isdigit(field[pos]) || field[pos] == '.' || field[pos] == '-'))
            tok += field[pos++];
        if (tok.empty()) return -1.0;
        try { return std::stod(tok); } catch (...) { return -1.0; }
    };

    int nLoaded = 0;
    std::string line;
    while (std::getline(in, line)) {
        if (line.size() < 10) continue;
        if (line[0] == '#') continue;

        // Cols 0-2: A, cols 4-6: Z (1-indexed in spec → 0-indexed here)
        std::string aStr = line.substr(0, 3);
        std::string zStr = (line.size() >= 7) ? line.substr(4, 3) : "";
        int A = 0, Z = 0;
        try { A = std::stoi(aStr); Z = std::stoi(zStr); } catch (...) { continue; }
        if (A <= 0 || Z < 0 || Z > 118 || Z > A) continue;

        // Isomer indicator at cols 7-10 — skip isomers (keep only ground state)
        if (line.size() >= 8) {
            char iso = line[7];
            if (iso == 'm' || iso == 'n') continue;  // isomeric state
        }

        NubaseEntry nb;

        // Half-life: value at cols 57-68, unit at cols 69-77 (approximate)
        // Use token search: find a time unit token
        if (line.size() > 69) {
            // The half-life field spans ~cols 57-77; extract it
            std::string hlField = line.substr(57, std::min((int)line.size() - 57, 22));
            std::istringstream hls(hlField);
            std::string tok1, tok2;
            if (hls >> tok1 >> tok2) {
                // tok1 might be the value, tok2 the unit
                double factor = unitToSec(tok2);
                if (factor >= 0.0) {
                    // strip # and < > from value
                    for (char c : {'#','<','>',' '}) tok1.erase(std::remove(tok1.begin(), tok1.end(), c), tok1.end());
                    try {
                        nb.halflifeSec  = std::stod(tok1) * factor;
                        // Build human-readable string
                        nb.halflifeStr  = tok1 + " " + tok2;
                    } catch (...) {}
                } else if (tok1 == "stbl" || tok2 == "stbl") {
                    nb.halflifeSec = 0.0;
                    nb.halflifeStr = "stable";
                }
            } else if (tok1 == "stbl") {
                nb.halflifeSec = 0.0;
                nb.halflifeStr = "stable";
            }
        }

        // Decay modes: from col 119 onward (if line is long enough)
        if (line.size() > 119) {
            std::string decayField = line.substr(119);
            // Replace semicolons with spaces for tokenizing, but keep = intact
            // e.g. "B-=85.7 4;B-N=14.0 4;B-2N=0.30 10"
            nb.brBetaMinus = parseBR(decayField, "B-");
            // Make sure "B-N" match doesn't steal from "B-"
            // parseBR for "B-" finds "B-=" but not "B-N=" — correct since B- comes before =
            nb.brBetaN     = parseBR(decayField, "B-N");
            nb.brBeta2N    = parseBR(decayField, "B-2N");
        }

        nubaseTable_[{Z, A}] = nb;
        ++nLoaded;
    }

    nubaseLoaded_ = (nLoaded > 100);
    return nubaseLoaded_;
}

void GammaFitGUI::OnLoadNubaseTable()
{
    TGFileInfo fi;
    static const char* kFileTypes[] = {
        "NUBASE files", "*.txt *.dat",
        "All files",    "*",
        nullptr,        nullptr
    };
    fi.fFileTypes = kFileTypes;
    fi.fIniDir    = StrDup(".");
    new TGFileDialog(gClient->GetRoot(), this, kFDOpen, &fi);
    if (!fi.fFilename) return;

    std::string path = fi.fFilename;
    bool ok = LoadNubaseTable(path);
    if (ok) {
        std::string disp = path;
        auto slash = disp.rfind('/');
        if (slash != std::string::npos) disp = disp.substr(slash + 1);
        if (nubaseLbl_) nubaseLbl_->SetText(disp.c_str());
        AppendLog("NUBASE table loaded: " + path + "  (" +
                  std::to_string(nubaseTable_.size()) + " nuclides)");
    } else {
        if (nubaseLbl_) nubaseLbl_->SetText("(load failed)");
        AppendLog("ERROR: Could not parse NUBASE table from " + path);
    }
}

// Fixed grid layout: col 0 = Parent, col 1 = Daughters, col 2 = Granddaughters.
// All 7 nodes are always drawn, even if no peaks are matched to them.

struct SchNodeDef { const char* cls; int col; int row; };
static const SchNodeDef kSchNodes[] = {
    {"Parent",               0, 1},
    {"Daughter",             1, 0},
    {"Beta-n Daughter",      1, 1},
    {"Beta-2n Daughter",     1, 2},
    {"Granddaughter",        2, 0},
    {"Beta-n Granddaughter", 2, 1},
    {"Beta-2n Granddaughter",2, 2},
};
static const int kNSchNodes = 7;

struct SchArrowDef { const char* from; const char* to; int neutrons; };
static const SchArrowDef kSchArrows[] = {
    {"Parent",          "Daughter",             0},
    {"Parent",          "Beta-n Daughter",      1},
    {"Parent",          "Beta-2n Daughter",     2},
    {"Daughter",        "Granddaughter",        0},
    {"Beta-n Daughter", "Beta-n Granddaughter", 0},
    {"Beta-2n Daughter","Beta-2n Granddaughter",0},
};
static const int kNSchArrows = 6;

void GammaFitGUI::DrawDecaySchematic(TCanvas* c)
{
    if (isoHistName_.empty()) return;

    FitDatabase fitdb;
    fitdb.Load(CacheFileFor(isoHistName_));
    LoadLabelClassMap(fitdb);

    // Collect peak counts and primary isotope name per class
    std::map<std::string, int>         classPeaks;
    std::map<std::string, std::string> classIsoName;

    for (const auto& kv : fitdb.GetEntries()) {
        if (kv.first.empty() || kv.first[0] == '_') continue;
        const FitEntry& fe = kv.second;

        auto addPeak = [&](const std::string& lbl, const std::string& cls) {
            std::string resolvedCls = cls;
            if (resolvedCls.empty() && !lbl.empty() && labelClassMap_.count(lbl))
                resolvedCls = labelClassMap_.at(lbl);
            if (resolvedCls.empty() && !lbl.empty())
                resolvedCls = AutoClassFromParent(lbl);
            if (resolvedCls.empty()) return;
            classPeaks[resolvedCls]++;
            if (!classIsoName.count(resolvedCls) && !lbl.empty())
                classIsoName[resolvedCls] = lbl;
        };

        int npar = (int)fe.params.size();
        if (npar >= 5 && (npar-2)%3==0) {
            int n = (npar-2)/3;
            for (int i = 0; i < n; i++)
                addPeak(fe.PeakLabel(i), fe.PeakClass(i));
        } else {
            addPeak(fe.label, fe.classification);
        }
    }

    c->Clear(); c->cd();
    c->SetFillColor(kWhite);
    gPad->SetFillColor(kWhite);
    gPad->Range(0, 0, 1, 1);   // TBox/TArrow use user coords; reset to NDC-space

    // ── NDC grid ──────────────────────────────────────────────────────────────
    const double xCol[3]  = {0.13, 0.47, 0.83};
    const double yRows[3] = {0.80, 0.50, 0.20};
    const double bW       = 0.090;  // box half-width
    const double bH       = 0.095;  // box half-height

    const int kFill[3] = {kOrange-9, kAzure-9, kGreen-9};

    // NDC center map — all 7 nodes always present
    std::map<std::string, std::pair<double,double>> pos;
    for (int i = 0; i < kNSchNodes; i++)
        pos[kSchNodes[i].cls] = {xCol[kSchNodes[i].col], yRows[kSchNodes[i].row]};

    TLatex ltx; ltx.SetNDC(kTRUE);

    // Q-value helper: Q_β- = Δ(from) - Δ(to) - neutrons × Δ(n) [keV]
    // Returns -9999.0 when mass data are unavailable (AME not loaded / parent not set).
    const double kNeutronMex = 8071.3171;  // keV
    auto computeQ = [&](const char* fromCls, const char* toCls, int neutrons) -> double {
        if (!ameLoaded_ || isoParentZval_ <= 0 || isoParentNval_ <= 0) return -9999.0;
        int Zf = 0, Nf = 0, Zt = 0, Nt = 0;
        if (!ClassZN(fromCls, isoParentZval_, isoParentNval_, Zf, Nf)) return -9999.0;
        if (!ClassZN(toCls,   isoParentZval_, isoParentNval_, Zt, Nt)) return -9999.0;
        auto ipf = ameTable_.find({Zf, Zf+Nf});
        auto ipt = ameTable_.find({Zt, Zt+Nt});
        if (ipf == ameTable_.end() || ipt == ameTable_.end()) return -9999.0;
        return ipf->second - ipt->second - neutrons * kNeutronMex;
    };

    // An arrow is drawn only when Q ≥ 0 (or when the AME table isn't loaded so
    // we have no way to judge).
    auto arrowEnabled = [&](const SchArrowDef& ad) -> bool {
        double Q = computeQ(ad.from, ad.to, ad.neutrons);
        return Q >= 0.0 || Q < -9000.0;  // enabled when Q>0 OR Q unknown
    };

    // ── Pass 1: colored box fills for ALL 7 nodes ─────────────────────────────
    for (int i = 0; i < kNSchNodes; i++) {
        const SchNodeDef& nd = kSchNodes[i];
        double xc = pos[nd.cls].first, yc = pos[nd.cls].second;
        TBox* fill = new TBox(xc-bW, yc-bH, xc+bW, yc+bH);
        fill->SetFillColor(kFill[nd.col]);
        fill->SetLineColor(kBlack);
        fill->SetLineWidth(2);
        fill->Draw();
    }

    // ── Pass 2: arrows (drawn on top of fills so arrowheads are always visible)
    for (int ai = 0; ai < kNSchArrows; ai++) {
        const SchArrowDef& ad = kSchArrows[ai];
        if (!arrowEnabled(ad)) continue;   // skip impossible decays (Q < 0)

        double fx = pos[ad.from].first,  fy = pos[ad.from].second;
        double tx = pos[ad.to].first,    ty = pos[ad.to].second;

        double ax1 = fx + bW, ay1 = fy;
        double ax2 = tx - bW, ay2 = ty;

        TArrow* arr = new TArrow(ax1, ay1, ax2, ay2, 0.016, "|>");
        arr->SetLineColor(kBlack); arr->SetFillColor(kBlack);
        arr->SetLineWidth(3); arr->Draw();

        // Perpendicular unit vector (always pointing toward top of canvas)
        double dx = ax2-ax1, dy = ay2-ay1;
        double len = std::sqrt(dx*dx + dy*dy);
        double px = (len > 0) ? -dy/len : 0.0;
        double py = (len > 0) ?  dx/len : 1.0;
        if (py < 0) { px = -px; py = -py; }

        double midX = (ax1+ax2)/2.0, midY = (ay1+ay2)/2.0;

        std::string decayLbl = "#beta^{-}";
        if (ad.neutrons == 1) decayLbl += "+n";
        if (ad.neutrons == 2) decayLbl += "+2n";

        ltx.SetTextSize(0.021); ltx.SetTextColor(kMagenta+1); ltx.SetTextAlign(22);
        ltx.DrawLatex(midX + px*0.038, midY + py*0.038, decayLbl.c_str());

        // Q value (green, slightly further offset)
        double Q = computeQ(ad.from, ad.to, ad.neutrons);
        if (Q > -9000.0) {
            ltx.SetTextSize(0.018); ltx.SetTextColor(kGreen+2); ltx.SetTextAlign(22);
            ltx.DrawLatex(midX + px*0.060, midY + py*0.060,
                          Form("Q=%.2f MeV", Q / 1000.0));
        }

        // Neutron symbol box for beta-n / beta-2n arrows
        if (ad.neutrons > 0) {
            double nox = midX + px*0.098, noy = midY + py*0.098;
            const double nW = 0.030, nH = 0.026;

            TArrow* narr = new TArrow(midX + px*0.066, midY + py*0.066,
                                      nox, noy - nH, 0.008, "|>");
            narr->SetLineColor(kRed+1); narr->SetFillColor(kRed+1);
            narr->SetLineStyle(2); narr->SetLineWidth(1); narr->Draw();

            TBox* nbox = new TBox(nox-nW, noy-nH, nox+nW, noy+nH);
            nbox->SetFillColor(kRed-9); nbox->SetLineColor(kBlack);
            nbox->SetLineWidth(1); nbox->Draw();

            std::string nlbl = (ad.neutrons == 1) ? "^{1}_{0}n" : "2 ^{1}_{0}n";
            ltx.SetTextSize(0.018); ltx.SetTextColor(kRed+2); ltx.SetTextAlign(22);
            ltx.DrawLatex(nox, noy - 0.010, nlbl.c_str());
        }
    }

    // ── Pass 3: box text (element symbol, Z, N, A, isotope, T½, BRs, Sn) ───────
    // Drawn last so text is never obscured by arrows.
    for (int i = 0; i < kNSchNodes; i++) {
        const SchNodeDef& nd = kSchNodes[i];
        double xc = pos[nd.cls].first, yc = pos[nd.cls].second;
        double x1 = xc-bW, x2 = xc+bW, y1 = yc-bH, y2 = yc+bH;

        // Nuclear identity from beta-minus chain
        int Z = 0, N = 0;
        bool hasZ = ClassZN(nd.cls, isoParentZval_, isoParentNval_, Z, N) && Z > 0 && N > 0;
        int A = Z + N;

        if (hasZ) {
            // Mass number — top left
            ltx.SetTextSize(0.018); ltx.SetTextColor(kBlack); ltx.SetTextAlign(11);
            ltx.DrawLatex(x1+0.005, y2-0.026, Form("%d", A));
            // Z — bottom left
            ltx.SetTextAlign(11);
            ltx.DrawLatex(x1+0.005, y1+0.006, Form("Z=%d", Z));
            // N — bottom right
            ltx.SetTextAlign(31);
            ltx.DrawLatex(x2-0.005, y1+0.006, Form("N=%d", N));
            // Element symbol — large, upper center
            ltx.SetTextSize(0.035); ltx.SetTextColor(kBlack); ltx.SetTextAlign(22);
            ltx.DrawLatex(xc, yc+0.018, ElementSymbol(Z));
        }

        // Class label (abbreviated)
        std::string scls = nd.cls;
        if      (scls == "Beta-n Daughter")        scls = "#beta^{-}n Dau.";
        else if (scls == "Beta-2n Daughter")       scls = "#beta^{-}2n Dau.";
        else if (scls == "Beta-n Granddaughter")   scls = "#beta^{-}n GD";
        else if (scls == "Beta-2n Granddaughter")  scls = "#beta^{-}2n GD";
        ltx.SetTextSize(0.013); ltx.SetTextColor(kGray+2); ltx.SetTextAlign(22);
        ltx.DrawLatex(xc, yc-(hasZ ? 0.003 : 0.008), scls.c_str());

        // Matched isotope name (blue)
        std::string isoName = classIsoName.count(nd.cls) ? classIsoName.at(nd.cls) : "";
        if (!isoName.empty()) {
            ltx.SetTextSize(0.014); ltx.SetTextColor(kBlue+1); ltx.SetTextAlign(22);
            ltx.DrawLatex(xc, yc-0.020, isoName.c_str());
        }

        // NUBASE: half-life and branching ratios
        double textY = yc - 0.034;
        if (nubaseLoaded_ && hasZ) {
            auto nbit = nubaseTable_.find({Z, A});
            if (nbit != nubaseTable_.end()) {
                const NubaseEntry& nb = nbit->second;
                if (!nb.halflifeStr.empty()) {
                    ltx.SetTextSize(0.012); ltx.SetTextColor(kGreen+2); ltx.SetTextAlign(22);
                    ltx.DrawLatex(xc, textY, Form("T_{1/2}=%s", nb.halflifeStr.c_str()));
                    textY -= 0.014;
                }
                // Beta-minus branching ratios
                std::string brStr;
                if (nb.brBetaMinus >= 0) brStr += Form("#beta^{-}=%.1f%%", nb.brBetaMinus);
                if (nb.brBetaN     >= 0) {
                    if (!brStr.empty()) brStr += " ";
                    brStr += Form("#beta^{-}n=%.1f%%", nb.brBetaN);
                }
                if (!brStr.empty()) {
                    ltx.SetTextSize(0.011); ltx.SetTextColor(kMagenta+1); ltx.SetTextAlign(22);
                    ltx.DrawLatex(xc, textY, brStr.c_str());
                    textY -= 0.013;
                }
                if (nb.brBeta2N >= 0) {
                    ltx.SetTextSize(0.011); ltx.SetTextColor(kMagenta+1); ltx.SetTextAlign(22);
                    ltx.DrawLatex(xc, textY, Form("#beta^{-}2n=%.2f%%", nb.brBeta2N));
                    textY -= 0.013;
                }
            }
        }

        // AME: neutron separation energy Sn = Δ(Z,A-1) + Δ_n - Δ(Z,A)
        if (ameLoaded_ && hasZ && A > 1) {
            auto it0 = ameTable_.find({Z, A});
            auto it1 = ameTable_.find({Z, A-1});
            if (it0 != ameTable_.end() && it1 != ameTable_.end()) {
                double Sn = it1->second + kNeutronMex - it0->second;  // keV
                ltx.SetTextSize(0.012); ltx.SetTextColor(kOrange+2); ltx.SetTextAlign(22);
                ltx.DrawLatex(xc, textY, Form("S_{n}=%.2f MeV", Sn * 1e-3));
                textY -= 0.013;
            }
        }

        // Peak count — only shown when > 0
        int nPks = classPeaks.count(nd.cls) ? classPeaks.at(nd.cls) : 0;
        if (nPks > 0) {
            ltx.SetTextSize(0.013); ltx.SetTextColor(kBlack); ltx.SetTextAlign(22);
            ltx.DrawLatex(xc, y1+0.012, Form("%d peak%s", nPks, nPks==1?"":"s"));
        }
    }

    // ── Title ─────────────────────────────────────────────────────────────────
    ltx.SetTextSize(0.026); ltx.SetTextColor(kBlack); ltx.SetTextAlign(22);
    if (!isoParentIsotope_.empty() && isoParentZval_ > 0) {
        int A = isoParentZval_ + isoParentNval_;
        ltx.DrawLatex(0.50, 0.96,
            Form("Decay Chain: %s   |   Parent: %s (^{%d}%s, Z=%d, N=%d)",
                 isoHistName_.c_str(), isoParentIsotope_.c_str(), A,
                 ElementSymbol(isoParentZval_), isoParentZval_, isoParentNval_));
    } else {
        ltx.DrawLatex(0.50, 0.96,
            Form("Decay Chain: %s   |   Set parent nucleus (Z, N) to show element symbols",
                 isoHistName_.c_str()));
    }

    c->Modified(); c->Update();
}

void GammaFitGUI::OnIsoDrawSchematic()
{
    if (isoHistName_.empty()) {
        AppendLog("Click 'Refresh' to load a histogram's cache first."); return;
    }
    schematicDrawn_ = true;
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
