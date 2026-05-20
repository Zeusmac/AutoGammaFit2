#include "GammaFitGUI.h"
#include "GammaFitGUI_shared.h"
#include "FitDatabase.h"

#include "TGTextEntry.h"
#include "TGFileDialog.h"
#include "TH1.h"
#include "TFile.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TLine.h"

#include <algorithm>
#include <cstring>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <set>
#include <dirent.h>

// ln(ε) = a - b·ln(E) + c·ln(E)² - d/E²
static double EvalEfficiency(double E_keV, double a, double b, double c, double d) {
    if (E_keV <= 0.0) return 0.0;
    double lnE = std::log(E_keV);
    return std::exp(a - b*lnE + c*lnE*lnE - d/(E_keV*E_keV));
}

// ─────────────────────────────────────────────────────────────────────────────
// LoadPeakCacheIntoTable — load all peaks from a cache file
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::LoadPeakCacheIntoTable(const std::string& cachePath,
                                          const std::string& histName)
{
    FitDatabase fitdb;
    fitdb.cacheOnly = true;
    if (!fitdb.Load(cachePath)) return;

    bool alreadyLoaded = false;
    for (const auto& p : ptLoadedCaches_)
        if (p == cachePath) { alreadyLoaded = true; break; }
    if (!alreadyLoaded) {
        ptLoadedCaches_.push_back(cachePath);
        if (ptCacheList_) {
            int id = (int)ptCacheList_->GetNumberOfEntries() + 1;
            std::string fname = cachePath;
            size_t sl = fname.rfind('/');
            if (sl != std::string::npos) fname = fname.substr(sl+1);
            ptCacheList_->AddEntry(fname.c_str(), id);
            ptCacheList_->MapSubwindows();
            ptCacheList_->Layout();
        }
    }

    int histIdx = (int)ptLoadedCaches_.size() - 1;

    for (const auto& kv : fitdb.GetEntries()) {
        if (kv.first.empty() || kv.first[0] == '_') continue;
        const FitEntry& fe = kv.second;

        FitLayout lay = DetectLayout((int)fe.params.size());
        if (!lay.valid()) continue;

        double bw   = (fe.binWidth > 0.0) ? fe.binWidth : 1.0;
        double chi2 = fe.chi2ndf;

        for (int gi = 0; gi < lay.n; gi++) {
            double A   = fe.params[3*gi];
            double E   = fe.params[3*gi + 1];
            double sig = fe.params[3*gi + 2];

            if (sig <= 0.0 || A <= 0.0) continue;

            double area    = A * sig * std::sqrt(2.0 * M_PI) / bw;
            double areaErr = 0.0;
            if (!fe.paramErrors.empty() && (int)fe.paramErrors.size() > 3*gi) {
                double dA   = (3*gi   < (int)fe.paramErrors.size()) ? fe.paramErrors[3*gi]   : 0.0;
                double dSig = (3*gi+2 < (int)fe.paramErrors.size()) ? fe.paramErrors[3*gi+2] : 0.0;
                double relA   = (A   > 0.0) ? dA/A     : 0.0;
                double relSig = (sig > 0.0) ? dSig/sig : 0.0;
                areaErr = area * std::sqrt(relA*relA + relSig*relSig);
            }

            double fwhm      = 2.3548200 * sig;
            double energyErr = 0.0;
            if ((int)fe.paramErrors.size() > 3*gi+1)
                energyErr = fe.paramErrors[3*gi+1];

            PeakTableRow row;
            row.histName       = histName;
            row.histIdx        = histIdx;
            row.energy         = E;
            row.energyErr      = energyErr;
            row.sigma          = sig;
            row.fwhm           = fwhm;
            row.area           = area;
            row.areaErr        = areaErr;
            row.chi2ndf        = chi2;
            row.label          = fe.PeakLabel(gi);
            row.classification = fe.PeakClass(gi);
            row.cacheFile      = cachePath;

            ptRows_.push_back(row);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildPeakTableTab
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::BuildPeakTableTab(TGCompositeFrame* p) {
    TGCanvas* sc = new TGCanvas(p, 308, 860, kSunkenFrame);
    p->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    TGCompositeFrame* cf = new TGCompositeFrame(sc->GetViewPort(), 295, 10, kVerticalFrame);
    sc->SetContainer(cf);
    p = cf;

    // ── Cache Sources ─────────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Cache Sources");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

        TGHorizontalFrame* btnRow = new TGHorizontalFrame(grp);
        grp->AddFrame(btnRow, new TGLayoutHints(kLHintsLeft, 2, 2, 2, 2));

        TGTextButton* scanBtn = new TGTextButton(btnRow, " Scan & Load All ");
        btnRow->AddFrame(scanBtn, new TGLayoutHints(kLHintsLeft, 0, 3, 0, 0));
        scanBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPTScanAll()");
        scanBtn->SetToolTipText("Scan the fit_caches directory and load all .dat files");

        TGTextButton* addBtn = new TGTextButton(btnRow, " Add Cache ");
        btnRow->AddFrame(addBtn, new TGLayoutHints(kLHintsLeft, 0, 3, 0, 0));
        addBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPTAddCache()");
        addBtn->SetToolTipText("Browse for a specific cache file to add");

        TGTextButton* removeBtn = new TGTextButton(btnRow, " Remove ");
        btnRow->AddFrame(removeBtn, new TGLayoutHints(kLHintsLeft, 0, 3, 0, 0));
        removeBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPTRemoveCache()");

        TGTextButton* clearBtn = new TGTextButton(btnRow, " Clear ");
        btnRow->AddFrame(clearBtn, new TGLayoutHints(kLHintsLeft));
        clearBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPTClearCaches()");

        ptCacheList_ = new TGListBox(grp, -1);
        ptCacheList_->Resize(280, 60);
        grp->AddFrame(ptCacheList_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    }

    // ── Table ─────────────────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Peak Table");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        // Filter row
        TGHorizontalFrame* filterRow = new TGHorizontalFrame(grp);
        grp->AddFrame(filterRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        filterRow->AddFrame(new TGLabel(filterRow, "Filter:"),
                            new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptFilterLabel_ = new TGTextEntry(filterRow, "");
        ptFilterLabel_->Resize(70, 22);
        filterRow->AddFrame(ptFilterLabel_, new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));

        filterRow->AddFrame(new TGLabel(filterRow, "Class:"),
                            new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptClassFilter_ = new TGComboBox(filterRow, 5200);
        ptClassFilter_->AddEntry("All",          1);
        ptClassFilter_->AddEntry("Parent",       2);
        ptClassFilter_->AddEntry("Daughter",     3);
        ptClassFilter_->AddEntry("Granddaughter",4);
        ptClassFilter_->AddEntry("Background",   5);
        ptClassFilter_->AddEntry("Custom",       6);
        ptClassFilter_->AddEntry("(unlabeled)",  7);
        ptClassFilter_->Select(1, kFALSE);
        ptClassFilter_->Resize(90, 22);
        filterRow->AddFrame(ptClassFilter_, new TGLayoutHints(kLHintsCenterY));

        // Sort + rebuild row
        TGHorizontalFrame* sortRow = new TGHorizontalFrame(grp);
        grp->AddFrame(sortRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        sortRow->AddFrame(new TGLabel(sortRow, "Sort:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptSortCombo_ = new TGComboBox(sortRow, 5201);
        ptSortCombo_->AddEntry("Energy",    1);
        ptSortCombo_->AddEntry("Area",      2);
        ptSortCombo_->AddEntry("FWHM",      3);
        ptSortCombo_->AddEntry("Chi2/NDF",  4);
        ptSortCombo_->AddEntry("Label",     5);
        ptSortCombo_->AddEntry("Histogram", 6);
        ptSortCombo_->Select(1, kFALSE);
        ptSortCombo_->Resize(90, 22);
        sortRow->AddFrame(ptSortCombo_, new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));

        TGTextButton* rebuildBtn = new TGTextButton(sortRow, " Rebuild ");
        sortRow->AddFrame(rebuildBtn, new TGLayoutHints(kLHintsCenterY, 0, 3, 0, 0));
        rebuildBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPTRebuildTable()");

        // Selectable list box — click a row to preview the peak
        ptTableList_ = new TGListBox(grp, 5210);
        ptTableList_->Resize(280, 220);
        grp->AddFrame(ptTableList_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        ptTableList_->AddEntry("  (load caches and click Rebuild)", 1);
        ptTableList_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                              "OnPTRowSelected(Int_t)");
    }

    // ── Graph ─────────────────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Graph");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

        TGHorizontalFrame* axisRow = new TGHorizontalFrame(grp);
        grp->AddFrame(axisRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        axisRow->AddFrame(new TGLabel(axisRow, "X:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptXAxisCombo_ = new TGComboBox(axisRow, 5202);
        ptXAxisCombo_->AddEntry("Histogram Index", 1);
        ptXAxisCombo_->AddEntry("Energy",          2);
        ptXAxisCombo_->AddEntry("Area",            3);
        ptXAxisCombo_->AddEntry("FWHM",            4);
        ptXAxisCombo_->AddEntry("Chi2/NDF",        5);
        ptXAxisCombo_->Select(2, kFALSE);
        ptXAxisCombo_->Resize(110, 22);
        axisRow->AddFrame(ptXAxisCombo_, new TGLayoutHints(kLHintsCenterY, 0, 8, 0, 0));

        axisRow->AddFrame(new TGLabel(axisRow, "Y:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptYAxisCombo_ = new TGComboBox(axisRow, 5203);
        ptYAxisCombo_->AddEntry("Area",     1);
        ptYAxisCombo_->AddEntry("FWHM",     2);
        ptYAxisCombo_->AddEntry("Energy",   3);
        ptYAxisCombo_->AddEntry("Chi2/NDF", 4);
        ptYAxisCombo_->AddEntry("Sigma",    5);
        ptYAxisCombo_->Select(1, kFALSE);
        ptYAxisCombo_->Resize(110, 22);
        axisRow->AddFrame(ptYAxisCombo_, new TGLayoutHints(kLHintsCenterY));

        TGHorizontalFrame* plotRow = new TGHorizontalFrame(grp);
        grp->AddFrame(plotRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 4));

        plotRow->AddFrame(new TGLabel(plotRow, "Filter label:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptGraphLabel_ = new TGTextEntry(plotRow, "");
        ptGraphLabel_->Resize(70, 22);
        plotRow->AddFrame(ptGraphLabel_, new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));

        TGTextButton* plotBtn = new TGTextButton(plotRow, " Plot ");
        plotRow->AddFrame(plotBtn, new TGLayoutHints(kLHintsCenterY, 0, 3, 0, 0));
        plotBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPTPlot()");
        plotBtn->SetToolTipText("Draw selected columns as a TGraph on the main canvas");

        TGTextButton* csvBtn = new TGTextButton(plotRow, " Export CSV ");
        plotRow->AddFrame(csvBtn, new TGLayoutHints(kLHintsCenterY));
        csvBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPTExportCSV()");
        csvBtn->SetToolTipText("Export the peak table to a CSV file");
    }

    // ── Intensity Calculation ─────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Intensity Calculation");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

        // Efficiency model row
        TGHorizontalFrame* effRow = new TGHorizontalFrame(grp);
        grp->AddFrame(effRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        effRow->AddFrame(new TGLabel(effRow, "Eff model:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptEffCombo_ = new TGComboBox(effRow, 5220);
        ptEffCombo_->AddEntry("Manual", 1);
        ptEffCombo_->Select(1, kFALSE);
        ptEffCombo_->Resize(100, 22);
        effRow->AddFrame(ptEffCombo_,
                         new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 0, 3, 0, 0));
        ptEffCombo_->Connect("Selected(Int_t)", "GammaFitGUI", this,
                             "OnPTEffSelected(Int_t)");

        TGTextButton* scanEffBtn = new TGTextButton(effRow, "Scan");
        effRow->AddFrame(scanEffBtn, new TGLayoutHints(kLHintsCenterY, 0, 2, 0, 0));
        scanEffBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPTScanEffCaches()");
        scanEffBtn->SetToolTipText("Load saved efficiency fits from eff_caches.dat");

        // Save-current-eff row
        TGHorizontalFrame* saveEffRow = new TGHorizontalFrame(grp);
        grp->AddFrame(saveEffRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        saveEffRow->AddFrame(new TGLabel(saveEffRow, "Save as:"),
                             new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptEffNameEntry_ = new TGTextEntry(saveEffRow, "");
        ptEffNameEntry_->Resize(110, 22);
        ptEffNameEntry_->SetToolTipText("Name for this efficiency fit (used when saving)");
        saveEffRow->AddFrame(ptEffNameEntry_,
                             new TGLayoutHints(kLHintsExpandX | kLHintsCenterY, 0, 3, 0, 0));
        TGTextButton* saveEffBtn = new TGTextButton(saveEffRow, "Save eff.");
        saveEffRow->AddFrame(saveEffBtn, new TGLayoutHints(kLHintsCenterY));
        saveEffBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPTSaveEffCache()");
        saveEffBtn->SetToolTipText(
            "Save the current AutoFit efficiency parameters (a,b,c,d) as a named cache");

        // Activity + time row
        TGHorizontalFrame* actRow = new TGHorizontalFrame(grp);
        grp->AddFrame(actRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        actRow->AddFrame(new TGLabel(actRow, "Activity:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptActivityEntry_ = new TGNumberEntry(actRow, 1.0, 10, -1,
                                             TGNumberFormat::kNESRealFour,
                                             TGNumberFormat::kNEAPositive);
        ptActivityEntry_->SetWidth(85);
        actRow->AddFrame(ptActivityEntry_, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
        actRow->AddFrame(new TGLabel(actRow, "Bq"),
                         new TGLayoutHints(kLHintsCenterY, 0, 8, 0, 0));
        actRow->AddFrame(new TGLabel(actRow, "Time:"),
                         new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptTimeEntry_ = new TGNumberEntry(actRow, 1.0, 10, -1,
                                         TGNumberFormat::kNESRealFour,
                                         TGNumberFormat::kNEAPositive);
        ptTimeEntry_->SetWidth(75);
        actRow->AddFrame(ptTimeEntry_, new TGLayoutHints(kLHintsLeft, 0, 2, 0, 0));
        actRow->AddFrame(new TGLabel(actRow, "s"),
                         new TGLayoutHints(kLHintsCenterY));

        // Energy + efficiency row (auto-filled on peak select)
        TGHorizontalFrame* eRow = new TGHorizontalFrame(grp);
        grp->AddFrame(eRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        eRow->AddFrame(new TGLabel(eRow, "E (keV):"),
                       new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptEnergyEntry_ = new TGNumberEntry(eRow, 0.0, 8, -1,
                                           TGNumberFormat::kNESRealFour,
                                           TGNumberFormat::kNEAPositive);
        ptEnergyEntry_->SetWidth(75);
        eRow->AddFrame(ptEnergyEntry_, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));
        eRow->AddFrame(new TGLabel(eRow, "Eff:"),
                       new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptEffValEntry_ = new TGNumberEntry(eRow, 1.0, 8, -1,
                                           TGNumberFormat::kNESReal,
                                           TGNumberFormat::kNEAPositive);
        ptEffValEntry_->SetWidth(80);
        eRow->AddFrame(ptEffValEntry_, new TGLayoutHints(kLHintsLeft));

        // Area row (auto-filled on peak select)
        TGHorizontalFrame* areaRow = new TGHorizontalFrame(grp);
        grp->AddFrame(areaRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        areaRow->AddFrame(new TGLabel(areaRow, "Area:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptAreaEntry_ = new TGNumberEntry(areaRow, 0.0, 10, -1,
                                         TGNumberFormat::kNESRealFour,
                                         TGNumberFormat::kNEAAnyNumber);
        ptAreaEntry_->SetWidth(90);
        areaRow->AddFrame(ptAreaEntry_, new TGLayoutHints(kLHintsLeft, 0, 4, 0, 0));
        areaRow->AddFrame(new TGLabel(areaRow, "+/-"),
                          new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptAreaErrEntry_ = new TGNumberEntry(areaRow, 0.0, 8, -1,
                                            TGNumberFormat::kNESRealFour,
                                            TGNumberFormat::kNEAAnyNumber);
        ptAreaErrEntry_->SetWidth(75);
        areaRow->AddFrame(ptAreaErrEntry_, new TGLayoutHints(kLHintsLeft, 0, 6, 0, 0));
        TGTextButton* calcBtn = new TGTextButton(areaRow, "Calculate");
        areaRow->AddFrame(calcBtn, new TGLayoutHints(kLHintsCenterY));
        calcBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPTCalculateIntensity()");

        // Result
        ptIntensityLbl_ = new TGLabel(grp, "  I = ---");
        ptIntensityLbl_->SetTextJustify(kTextLeft);
        grp->AddFrame(ptIntensityLbl_,
                      new TGLayoutHints(kLHintsExpandX, 4, 2, 2, 4));
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnPTScanAll
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnPTScanAll() {
    std::string cacheDir = CacheDirFor();
    DIR* dir = opendir(cacheDir.c_str());
    if (!dir) {
        AppendLog("PeakTable: cannot open cache dir: " + cacheDir);
        return;
    }

    std::set<std::string> alreadySet(ptLoadedCaches_.begin(), ptLoadedCaches_.end());
    int nLoaded = 0;

    struct dirent* entry;
    while ((entry = readdir(dir)) != nullptr) {
        std::string fname = entry->d_name;
        if (fname.size() < 5) continue;
        if (fname.substr(fname.size()-4) != ".dat") continue;
        if (fname[0] == '_') continue;

        std::string fullPath = cacheDir + "/" + fname;
        if (alreadySet.count(fullPath)) continue;

        std::string histName = fname;
        const std::string prefix = "fit_cache_";
        if (histName.size() > prefix.size() &&
            histName.substr(0, prefix.size()) == prefix)
            histName = histName.substr(prefix.size(),
                                       histName.size() - prefix.size() - 4);
        else
            histName = fname;

        LoadPeakCacheIntoTable(fullPath, histName);
        alreadySet.insert(fullPath);
        nLoaded++;
    }
    closedir(dir);

    AppendLog("PeakTable: loaded " + std::to_string(nLoaded) +
              " cache files (" + std::to_string(ptRows_.size()) + " peaks total)");
    OnPTRebuildTable();
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnPTAddCache
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnPTAddCache() {
    static const char* types[] = {
        "Fit cache files", "*.dat",
        "All files",       "*",
        nullptr, nullptr
    };
    TGFileInfo fi;
    fi.fFileTypes = types;
    std::string initDir = CacheDirFor();
    fi.fIniDir = const_cast<char*>(initDir.c_str());
    OpenFileDialog(this, kFDOpen, &fi);
    if (!fi.fFilename) return;

    std::string path = fi.fFilename;
    std::string fname = path;
    size_t sl = fname.rfind('/');
    if (sl != std::string::npos) fname = fname.substr(sl+1);

    std::string histName = fname;
    const std::string prefix = "fit_cache_";
    if (histName.size() > prefix.size() &&
        histName.substr(0, prefix.size()) == prefix)
        histName = histName.substr(prefix.size(),
                                   histName.size() - prefix.size() - 4);

    LoadPeakCacheIntoTable(path, histName);
    AppendLog("PeakTable: added " + fname);
    OnPTRebuildTable();
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnPTRemoveCache
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnPTRemoveCache() {
    if (!ptCacheList_) return;
    TGLBEntry* sel = ptCacheList_->GetSelectedEntry();
    if (!sel) { AppendLog("PeakTable: no cache selected"); return; }
    int selID = sel->EntryId() - 1;
    if (selID < 0 || selID >= (int)ptLoadedCaches_.size()) return;

    std::string removedPath = ptLoadedCaches_[selID];
    ptLoadedCaches_.erase(ptLoadedCaches_.begin() + selID);

    ptRows_.erase(
        std::remove_if(ptRows_.begin(), ptRows_.end(),
            [&](const PeakTableRow& r){ return r.cacheFile == removedPath; }),
        ptRows_.end());

    ptCacheList_->RemoveAll();
    for (int i = 0; i < (int)ptLoadedCaches_.size(); i++) {
        std::string fn = ptLoadedCaches_[i];
        size_t sl = fn.rfind('/'); if (sl!=std::string::npos) fn=fn.substr(sl+1);
        ptCacheList_->AddEntry(fn.c_str(), i+1);
    }
    ptCacheList_->MapSubwindows();
    ptCacheList_->Layout();

    for (auto& row : ptRows_) {
        for (int i = 0; i < (int)ptLoadedCaches_.size(); i++) {
            if (row.cacheFile == ptLoadedCaches_[i]) { row.histIdx = i; break; }
        }
    }

    AppendLog("PeakTable: removed cache and its peaks");
    OnPTRebuildTable();
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnPTClearCaches
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnPTClearCaches() {
    ptLoadedCaches_.clear();
    ptRows_.clear();
    ptFilteredRows_.clear();
    if (ptCacheList_) { ptCacheList_->RemoveAll(); ptCacheList_->MapSubwindows(); ptCacheList_->Layout(); }
    if (ptTableList_) {
        ptTableList_->RemoveAll();
        ptTableList_->AddEntry("  (cleared)", 1);
        ptTableList_->MapSubwindows();
        ptTableList_->Layout();
    }
    AppendLog("PeakTable: cleared all caches");
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnPTRebuildTable
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnPTRebuildTable() {
    if (!ptTableList_) return;

    std::string filterLabel;
    if (ptFilterLabel_) filterLabel = ptFilterLabel_->GetText();
    int classFilter = ptClassFilter_ ? ptClassFilter_->GetSelected() : 1;
    int sortSel     = ptSortCombo_   ? ptSortCombo_->GetSelected()   : 1;

    // Build filtered list
    std::vector<size_t> indices;
    for (size_t i = 0; i < ptRows_.size(); i++) {
        const PeakTableRow& r = ptRows_[i];
        if (!filterLabel.empty() && r.label.find(filterLabel) == std::string::npos)
            continue;
        if (classFilter == 2 && r.classification != "Parent")        continue;
        if (classFilter == 3 && r.classification != "Daughter")      continue;
        if (classFilter == 4 && r.classification != "Granddaughter") continue;
        if (classFilter == 5 && r.classification != "Background")    continue;
        if (classFilter == 6 && r.classification.substr(0,6) != "Custom") continue;
        if (classFilter == 7 && !r.label.empty())                    continue;
        indices.push_back(i);
    }

    // Sort
    std::sort(indices.begin(), indices.end(), [&](size_t a, size_t b) {
        const PeakTableRow& ra = ptRows_[a]; const PeakTableRow& rb = ptRows_[b];
        switch (sortSel) {
            case 2: return ra.area    < rb.area;
            case 3: return ra.fwhm    < rb.fwhm;
            case 4: return ra.chi2ndf < rb.chi2ndf;
            case 5: return ra.label   < rb.label;
            case 6: return ra.histIdx < rb.histIdx;
            default: return ra.energy < rb.energy;
        }
    });

    ptFilteredRows_ = indices;

    ptTableList_->RemoveAll();

    if (indices.empty()) {
        ptTableList_->AddEntry("  (no rows match filter)", 1);
        ptTableList_->MapSubwindows();
        ptTableList_->Layout();
        return;
    }

    // Each row: E=xxxx.xxx +/-x.xxx  FWHM=x.xxx  Area=xxxxxxx +/-xxxxx  chi2=x.xx  label [class]  {hist}
    char buf[512];
    for (int i = 0; i < (int)indices.size(); i++) {
        const PeakTableRow& r = ptRows_[indices[i]];

        std::string lbl = r.label.empty() ? "(unlabeled)" : r.label;
        std::string cls = r.classification.empty() ? "" : "[" + r.classification + "]";

        // Truncate histogram name to keep line readable
        std::string hname = r.histName;
        if (hname.size() > 20) hname = hname.substr(0, 17) + "...";

        std::snprintf(buf, sizeof(buf),
            "  E=%9.3f +/-%6.3f  FWHM=%6.3f  Area=%8.1f +/-%7.1f  chi2=%5.2f  %-18s %-16s {%s}",
            r.energy, r.energyErr, r.fwhm,
            r.area, r.areaErr, r.chi2ndf,
            lbl.c_str(), cls.c_str(), hname.c_str());

        ptTableList_->AddEntry(buf, i + 1);
    }

    ptTableList_->MapSubwindows();
    ptTableList_->Layout();
    AppendLog(Form("PeakTable: %d peaks shown", (int)indices.size()));
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnPTRowSelected — preview the selected peak on the canvas
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnPTRowSelected(Int_t id) {
    if (id < 1 || (size_t)id > ptFilteredRows_.size()) return;
    const PeakTableRow& r = ptRows_[ptFilteredRows_[id - 1]];

    // Try to find the histogram in the open file
    TH1* h = nullptr;
    if (inputFile_) {
        h = dynamic_cast<TH1*>(inputFile_->Get(r.histName.c_str()));
    }
    if (!h) {
        // Histogram not in open file — just log the info
        AppendLog(Form("PeakTable: E=%.3f keV  FWHM=%.3f  Area=%.1f +/-%.1f  %s  {%s}",
                       r.energy, r.fwhm, r.area, r.areaErr,
                       r.label.c_str(), r.histName.c_str()));
        return;
    }

    double sig = r.sigma;
    if (sig <= 0) sig = r.fwhm / 2.3548;
    if (sig <= 0) sig = 1.0;

    // Zoom to +/-8 sigma around the peak
    double lo = r.energy - 8.0 * sig;
    double hi = r.energy + 8.0 * sig;
    h->GetXaxis()->SetRangeUser(lo, hi);
    h->GetXaxis()->SetTitle("Energy (keV)");
    h->GetYaxis()->SetTitle("Counts");
    h->SetLineColor(kBlack);
    h->SetMarkerSize(0);

    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();
    h->Draw(showErrorBars_ ? "hist E" : "hist");
    c->Modified(); c->Update();

    // Draw vertical line at peak centre
    double ylo = gPad->GetUymin();
    double yhi = gPad->GetUymax();
    TLine* lc = new TLine(r.energy, ylo, r.energy, yhi);
    lc->SetLineColor(kRed); lc->SetLineWidth(2); lc->Draw();

    // Draw +/-sigma boundary lines
    TLine* ll = new TLine(r.energy - sig, ylo, r.energy - sig, yhi);
    ll->SetLineColor(kBlue+1); ll->SetLineWidth(1); ll->SetLineStyle(2); ll->Draw();
    TLine* lr = new TLine(r.energy + sig, ylo, r.energy + sig, yhi);
    lr->SetLineColor(kBlue+1); lr->SetLineWidth(1); lr->SetLineStyle(2); lr->Draw();

    c->Modified(); c->Update();
    gSystem->ProcessEvents();

    SetStatus(Form("PeakTable preview: E=%.3f keV  FWHM=%.3f  Area=%.1f+/-%.1f  %s",
                   r.energy, r.fwhm, r.area, r.areaErr, r.label.c_str()));

    // ── Auto-fill intensity calculation fields ────────────────────────────────
    if (ptEnergyEntry_)  ptEnergyEntry_->SetNumber(r.energy);
    if (ptAreaEntry_)    ptAreaEntry_->SetNumber(r.area);
    if (ptAreaErrEntry_) ptAreaErrEntry_->SetNumber(r.areaErr);

    // Compute efficiency from selected model
    if (ptEffCombo_ && ptEffValEntry_) {
        Int_t effSel = ptEffCombo_->GetSelected();
        // effSel == 1 → Manual; effSel >= 2 → ptEffCaches_[effSel-2]
        if (effSel >= 2) {
            int ci = effSel - 2;
            if (ci < (int)ptEffCaches_.size()) {
                const EfficiencyCache& ec = ptEffCaches_[ci];
                double eff = EvalEfficiency(r.energy, ec.a, ec.b, ec.c, ec.d);
                ptEffValEntry_->SetNumber(eff);
            }
        }
    }
    // Reset intensity display so stale result isn't shown for new peak
    if (ptIntensityLbl_) ptIntensityLbl_->SetText("  I = --- (click Calculate)");
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnPTPlot
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnPTPlot() {
    if (!ptXAxisCombo_ || !ptYAxisCombo_) return;

    int xSel = ptXAxisCombo_->GetSelected();
    int ySel = ptYAxisCombo_->GetSelected();

    std::string labelFilter;
    if (ptGraphLabel_) labelFilter = ptGraphLabel_->GetText();

    static const char* xNames[] = {"", "Histogram Index", "Energy (keV)", "Area (counts)", "FWHM (keV)", "Chi2/NDF"};
    static const char* yNames[] = {"", "Area (counts)", "FWHM (keV)", "Energy (keV)", "Chi2/NDF", "Sigma (keV)"};

    auto getX = [&](const PeakTableRow& r) -> double {
        switch (xSel) {
            case 1: return (double)r.histIdx;
            case 2: return r.energy;
            case 3: return r.area;
            case 4: return r.fwhm;
            case 5: return r.chi2ndf;
            default: return r.energy;
        }
    };
    auto getY = [&](const PeakTableRow& r) -> double {
        switch (ySel) {
            case 1: return r.area;
            case 2: return r.fwhm;
            case 3: return r.energy;
            case 4: return r.chi2ndf;
            case 5: return r.sigma;
            default: return r.area;
        }
    };
    auto getYErr = [&](const PeakTableRow& r) -> double {
        if (ySel == 1) return r.areaErr;
        if (ySel == 3) return r.energyErr;
        return 0.0;
    };
    auto getXErr = [&](const PeakTableRow& r) -> double {
        if (xSel == 2) return r.energyErr;
        return 0.0;
    };

    std::vector<const PeakTableRow*> rows;
    for (const auto& r : ptRows_) {
        if (!labelFilter.empty() && r.label.find(labelFilter) == std::string::npos)
            continue;
        rows.push_back(&r);
    }

    if (rows.empty()) {
        AppendLog("PeakTable: no data to plot (check filter)");
        return;
    }

    int n = (int)rows.size();
    std::vector<double> vx(n), vy(n), vxe(n), vye(n);
    for (int i = 0; i < n; i++) {
        vx[i]  = getX(*rows[i]);
        vy[i]  = getY(*rows[i]);
        vxe[i] = getXErr(*rows[i]);
        vye[i] = getYErr(*rows[i]);
    }

    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();
    c->SetLeftMargin(0.14); c->SetBottomMargin(0.12);

    TGraphErrors* g = new TGraphErrors(n, vx.data(), vy.data(), vxe.data(), vye.data());
    g->SetBit(kCanDelete);
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.8);
    g->SetMarkerColor(kBlue+1);
    g->SetLineColor(kBlue+1);

    std::string title = std::string(xSel <= 5 ? xNames[xSel] : "?") +
                        " vs " +
                        std::string(ySel <= 5 ? yNames[ySel] : "?");
    if (!labelFilter.empty()) title += "  (" + labelFilter + ")";
    g->SetTitle(title.c_str());
    g->GetXaxis()->SetTitle(xSel <= 5 ? xNames[xSel] : "");
    g->GetYaxis()->SetTitle(ySel <= 5 ? yNames[ySel] : "");

    g->Draw("AP");
    c->Update();
    AppendLog("PeakTable: plotted " + std::to_string(n) + " points");
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnPTExportCSV
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnPTExportCSV() {
    static const char* types[] = {
        "CSV files", "*.csv",
        "All files", "*",
        nullptr, nullptr
    };
    TGFileInfo fi;
    fi.fFileTypes = types;
    OpenFileDialog(this, kFDSave, &fi);
    if (!fi.fFilename) return;

    std::string path = fi.fFilename;
    if (path.size() < 4 || path.substr(path.size()-4) != ".csv")
        path += ".csv";

    std::ofstream out(path);
    if (!out.is_open()) {
        AppendLog("PeakTable: cannot open for writing: " + path);
        return;
    }

    out << "Histogram,Energy_keV,EnergyErr_keV,Sigma_keV,FWHM_keV,Area,AreaErr,Chi2NDF,Label,Classification,CacheFile\n";

    for (const auto& r : ptRows_) {
        out << r.histName << ","
            << std::fixed << std::setprecision(4)
            << r.energy   << "," << r.energyErr << ","
            << r.sigma    << "," << r.fwhm      << ","
            << std::setprecision(2)
            << r.area     << "," << r.areaErr   << ","
            << r.chi2ndf  << ","
            << "\"" << r.label << "\","
            << "\"" << r.classification << "\","
            << "\"" << r.cacheFile << "\"\n";
    }
    out.close();
    AppendLog("PeakTable: exported " + std::to_string(ptRows_.size()) + " peaks to " + path);
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnPTScanEffCaches — load eff_caches.dat from the launch directory
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnPTScanEffCaches() {
    ptEffCaches_.clear();
    if (ptEffCombo_) {
        ptEffCombo_->RemoveAll();
        ptEffCombo_->AddEntry("Manual", 1);
        ptEffCombo_->Select(1, kFALSE);
    }

    std::string path = launchDir_ + "/eff_caches.dat";
    std::ifstream in(path);
    if (!in.is_open()) {
        AppendLog("PeakTable: no eff_caches.dat found at " + path);
        if (ptEffCombo_) ptEffCombo_->Layout();
        return;
    }

    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        EfficiencyCache ec;
        if (!(ss >> ec.name >> ec.a >> ec.b >> ec.c >> ec.d)) continue;
        ptEffCaches_.push_back(ec);
        if (ptEffCombo_)
            ptEffCombo_->AddEntry(ec.name.c_str(), (int)ptEffCaches_.size() + 1);
    }

    if (ptEffCombo_) { ptEffCombo_->MapSubwindows(); ptEffCombo_->Layout(); }
    AppendLog(Form("PeakTable: loaded %d efficiency fit(s)", (int)ptEffCaches_.size()));
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnPTSaveEffCache — save current effA..D params under a user-chosen name
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnPTSaveEffCache() {
    std::string name = ptEffNameEntry_ ? ptEffNameEntry_->GetText() : "";
    if (name.empty()) { AppendLog("PeakTable: enter a name before saving"); return; }

    double a = effA_ ? effA_->GetNumber() : 0.0;
    double b = effB_ ? effB_->GetNumber() : 0.0;
    double c = effC_ ? effC_->GetNumber() : 0.0;
    double d = effD_ ? effD_->GetNumber() : 0.0;

    // Append to eff_caches.dat
    std::string path = launchDir_ + "/eff_caches.dat";
    std::ofstream out(path, std::ios::app);
    if (!out.is_open()) {
        AppendLog("PeakTable: cannot write eff_caches.dat at " + path);
        return;
    }
    out << std::setprecision(10) << name << " " << a << " " << b
        << " " << c << " " << d << "\n";
    out.close();

    // Add to in-memory list and combo
    EfficiencyCache ec{ name, a, b, c, d };
    ptEffCaches_.push_back(ec);
    if (ptEffCombo_) {
        ptEffCombo_->AddEntry(name.c_str(), (int)ptEffCaches_.size() + 1);
        ptEffCombo_->MapSubwindows(); ptEffCombo_->Layout();
        ptEffCombo_->Select((int)ptEffCaches_.size() + 1, kFALSE);
    }
    AppendLog(Form("PeakTable: saved eff cache '%s'  a=%.4g b=%.4g c=%.4g d=%.4g",
                   name.c_str(), a, b, c, d));
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnPTEffSelected — recompute efficiency entry when model changes
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnPTEffSelected(Int_t id) {
    if (id <= 1 || !ptEffValEntry_ || !ptEnergyEntry_) return;  // "Manual" → user edits
    int ci = id - 2;
    if (ci < 0 || ci >= (int)ptEffCaches_.size()) return;

    double E = ptEnergyEntry_->GetNumber();
    if (E <= 0) return;

    double eff = EvalEfficiency(E, ptEffCaches_[ci].a, ptEffCaches_[ci].b,
                                   ptEffCaches_[ci].c, ptEffCaches_[ci].d);
    ptEffValEntry_->SetNumber(eff);
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnPTCalculateIntensity — I = Area / (ε × A₀ × t)
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnPTCalculateIntensity() {
    if (!ptIntensityLbl_) return;

    double area    = ptAreaEntry_     ? ptAreaEntry_->GetNumber()     : 0.0;
    double areaErr = ptAreaErrEntry_  ? ptAreaErrEntry_->GetNumber()  : 0.0;
    double eff     = ptEffValEntry_   ? ptEffValEntry_->GetNumber()   : 0.0;
    double actBq   = ptActivityEntry_ ? ptActivityEntry_->GetNumber() : 0.0;
    double time_s  = ptTimeEntry_     ? ptTimeEntry_->GetNumber()     : 0.0;

    if (eff <= 0 || actBq <= 0 || time_s <= 0) {
        ptIntensityLbl_->SetText("  I = --- (need eff, activity, and time > 0)");
        ptIntensityLbl_->Layout();
        return;
    }

    // Absolute intensity (gamma emission probability per decay)
    double denom  = eff * actBq * time_s;
    double I      = area    / denom;
    double I_err  = areaErr / denom;

    // Also compute as percentage
    char buf[256];
    if (I_err > 0.0)
        std::snprintf(buf, sizeof(buf),
            "  I = %.4g +/- %.4g  (%.3f +/- %.3f %%)",
            I, I_err, I * 100.0, I_err * 100.0);
    else
        std::snprintf(buf, sizeof(buf),
            "  I = %.4g  (%.3f %%)", I, I * 100.0);

    ptIntensityLbl_->SetText(buf);
    ptIntensityLbl_->Layout();

    AppendLog(Form("PeakTable: intensity E=%.3f keV  I=%.4g  eff=%.4g  A=%.4g Bq  t=%.4g s",
                   ptEnergyEntry_ ? ptEnergyEntry_->GetNumber() : 0.0,
                   I, eff, actBq, time_s));
}
