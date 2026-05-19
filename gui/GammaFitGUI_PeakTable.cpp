#include "GammaFitGUI.h"
#include "GammaFitGUI_shared.h"
#include "FitDatabase.h"

#include "TGTextEntry.h"
#include "TGFileDialog.h"
#include "TH1.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TAxis.h"
#include "TMath.h"
#include "TROOT.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <set>
#include <dirent.h>

// ─────────────────────────────────────────────────────────────────────────────
// LoadPeakCacheIntoTable — load all peaks from a cache file
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::LoadPeakCacheIntoTable(const std::string& cachePath,
                                          const std::string& histName)
{
    FitDatabase fitdb;
    fitdb.cacheOnly = true;
    if (!fitdb.Load(cachePath)) return;

    // Check if histName already in ptCacheList_
    bool alreadyLoaded = false;
    for (const auto& p : ptLoadedCaches_)
        if (p == cachePath) { alreadyLoaded = true; break; }
    if (!alreadyLoaded) {
        ptLoadedCaches_.push_back(cachePath);
        if (ptCacheList_) {
            int id = (int)ptCacheList_->GetNumberOfEntries() + 1;
            // Show just the file name
            std::string fname = cachePath;
            size_t sl = fname.rfind('/');
            if (sl != std::string::npos) fname = fname.substr(sl+1);
            ptCacheList_->AddEntry(fname.c_str(), id);
            ptCacheList_->Layout();
        }
    }

    int histIdx = (int)ptLoadedCaches_.size() - 1;

    for (const auto& kv : fitdb.GetEntries()) {
        if (kv.first.empty() || kv.first[0] == '_') continue;
        const FitEntry& fe = kv.second;

        FitLayout lay = DetectLayout((int)fe.params.size());
        if (!lay.valid()) continue;

        double bw = (fe.binWidth > 0.0) ? fe.binWidth : 1.0;
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
                // Error propagation: δ(area) = area * sqrt((δA/A)^2 + (δσ/σ)^2)
                double relA   = (A   > 0.0) ? dA/A     : 0.0;
                double relSig = (sig > 0.0) ? dSig/sig : 0.0;
                areaErr = area * std::sqrt(relA*relA + relSig*relSig);
            }

            double fwhm   = 2.3548200 * sig;
            double energyErr = 0.0;
            if ((int)fe.paramErrors.size() > 3*gi+1)
                energyErr = fe.paramErrors[3*gi+1];

            PeakTableRow row;
            row.histName   = histName;
            row.histIdx    = histIdx;
            row.energy     = E;
            row.energyErr  = energyErr;
            row.sigma      = sig;
            row.fwhm       = fwhm;
            row.area       = area;
            row.areaErr    = areaErr;
            row.chi2ndf    = chi2;
            row.label      = fe.PeakLabel(gi);
            row.classification = fe.PeakClass(gi);
            row.cacheFile  = cachePath;

            ptRows_.push_back(row);
        }
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// BuildPeakTableTab
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::BuildPeakTableTab(TGCompositeFrame* p) {
    // Scrollable wrapper
    TGCanvas* sc = new TGCanvas(p, 308, 860, kSunkenFrame);
    p->AddFrame(sc, new TGLayoutHints(kLHintsExpandX | kLHintsExpandY));
    TGCompositeFrame* cf = new TGCompositeFrame(sc->GetViewPort(), 295, 10, kVerticalFrame);
    sc->SetContainer(cf);
    p = cf;

    // ── Cache Sources ─────────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Cache Sources");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 4, 2));

        // Button row
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

        // Cache list box
        ptCacheList_ = new TGListBox(grp, -1);
        ptCacheList_->Resize(280, 70);
        grp->AddFrame(ptCacheList_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
    }

    // ── Table ─────────────────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Table");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 2));

        // Filter row
        TGHorizontalFrame* filterRow = new TGHorizontalFrame(grp);
        grp->AddFrame(filterRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        filterRow->AddFrame(new TGLabel(filterRow, "Filter label:"),
                            new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptFilterLabel_ = new TGTextEntry(filterRow, "");
        ptFilterLabel_->Resize(80, 22);
        filterRow->AddFrame(ptFilterLabel_, new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));

        filterRow->AddFrame(new TGLabel(filterRow, "Class:"),
                            new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptClassFilter_ = new TGComboBox(filterRow, 5200);
        ptClassFilter_->AddEntry("All",                  1);
        ptClassFilter_->AddEntry("Parent",               2);
        ptClassFilter_->AddEntry("Daughter",             3);
        ptClassFilter_->AddEntry("Granddaughter",        4);
        ptClassFilter_->AddEntry("Background",           5);
        ptClassFilter_->AddEntry("Custom",               6);
        ptClassFilter_->AddEntry("(unlabeled)",          7);
        ptClassFilter_->Select(1, kFALSE);
        ptClassFilter_->Resize(100, 22);
        filterRow->AddFrame(ptClassFilter_, new TGLayoutHints(kLHintsCenterY));

        // Sort row
        TGHorizontalFrame* sortRow = new TGHorizontalFrame(grp);
        grp->AddFrame(sortRow, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));

        sortRow->AddFrame(new TGLabel(sortRow, "Sort by:"),
                          new TGLayoutHints(kLHintsCenterY, 0, 4, 0, 0));
        ptSortCombo_ = new TGComboBox(sortRow, 5201);
        ptSortCombo_->AddEntry("Energy",      1);
        ptSortCombo_->AddEntry("Area",        2);
        ptSortCombo_->AddEntry("FWHM",        3);
        ptSortCombo_->AddEntry("Chi2/NDF",    4);
        ptSortCombo_->AddEntry("Label",       5);
        ptSortCombo_->AddEntry("Histogram",   6);
        ptSortCombo_->Select(1, kFALSE);
        ptSortCombo_->Resize(100, 22);
        sortRow->AddFrame(ptSortCombo_, new TGLayoutHints(kLHintsCenterY, 0, 6, 0, 0));

        TGTextButton* rebuildBtn = new TGTextButton(sortRow, " Rebuild Table ");
        sortRow->AddFrame(rebuildBtn, new TGLayoutHints(kLHintsCenterY));
        rebuildBtn->Connect("Clicked()", "GammaFitGUI", this, "OnPTRebuildTable()");

        // Table text view
        ptTableView_ = new TGTextView(grp, 280, 250);
        grp->AddFrame(ptTableView_, new TGLayoutHints(kLHintsExpandX, 2, 2, 2, 2));
        ptTableView_->AddLine("  (load caches and click Rebuild Table)");
    }

    // ── Graph ─────────────────────────────────────────────────────────────────
    {
        TGGroupFrame* grp = new TGGroupFrame(p, "Graph");
        p->AddFrame(grp, new TGLayoutHints(kLHintsExpandX, 4, 4, 2, 4));

        // X axis
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
        ptYAxisCombo_->AddEntry("Area",        1);
        ptYAxisCombo_->AddEntry("FWHM",        2);
        ptYAxisCombo_->AddEntry("Energy",      3);
        ptYAxisCombo_->AddEntry("Chi2/NDF",    4);
        ptYAxisCombo_->AddEntry("Sigma",       5);
        ptYAxisCombo_->Select(1, kFALSE);
        ptYAxisCombo_->Resize(110, 22);
        axisRow->AddFrame(ptYAxisCombo_, new TGLayoutHints(kLHintsCenterY));

        // Label filter + buttons
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
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnPTScanAll — scan CacheDirFor() for .dat files and load all
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
        // Accept .dat files that look like fit cache files
        if (fname.size() < 5) continue;
        if (fname.substr(fname.size()-4) != ".dat") continue;
        if (fname[0] == '_') continue;

        std::string fullPath = cacheDir + "/" + fname;
        if (alreadySet.count(fullPath)) continue;

        // Derive histogram name from filename: fit_cache_HISTNAME.dat
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
// Slot: OnPTAddCache — browse for a cache file
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
    int selID = sel->EntryId() - 1;  // 0-based
    if (selID < 0 || selID >= (int)ptLoadedCaches_.size()) return;

    std::string removedPath = ptLoadedCaches_[selID];
    ptLoadedCaches_.erase(ptLoadedCaches_.begin() + selID);

    // Remove all rows from that cache
    ptRows_.erase(
        std::remove_if(ptRows_.begin(), ptRows_.end(),
            [&](const PeakTableRow& r){ return r.cacheFile == removedPath; }),
        ptRows_.end());

    // Rebuild list box
    ptCacheList_->RemoveAll();
    for (int i = 0; i < (int)ptLoadedCaches_.size(); i++) {
        std::string fname = ptLoadedCaches_[i];
        size_t sl = fname.rfind('/'); if (sl!=std::string::npos) fname=fname.substr(sl+1);
        ptCacheList_->AddEntry(fname.c_str(), i+1);
    }
    ptCacheList_->Layout();

    // Update histIdx for remaining rows
    for (auto& row : ptRows_) {
        // Recompute histIdx based on current ptLoadedCaches_
        for (int i = 0; i < (int)ptLoadedCaches_.size(); i++) {
            if (row.cacheFile == ptLoadedCaches_[i]) {
                row.histIdx = i; break;
            }
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
    if (ptCacheList_) { ptCacheList_->RemoveAll(); ptCacheList_->Layout(); }
    if (ptTableView_) {
        ptTableView_->Clear();
        ptTableView_->AddLine("  (cleared)");
        ptTableView_->Update();
    }
    AppendLog("PeakTable: cleared all caches");
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnPTRebuildTable — rebuild the table view with filter/sort
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnPTRebuildTable() {
    if (!ptTableView_) return;

    // Get filter settings
    std::string filterLabel;
    if (ptFilterLabel_) filterLabel = ptFilterLabel_->GetText();

    int classFilter = ptClassFilter_ ? ptClassFilter_->GetSelected() : 1;
    int sortSel     = ptSortCombo_   ? ptSortCombo_->GetSelected()   : 1;

    // Build working copy of rows
    std::vector<const PeakTableRow*> rows;
    for (const auto& r : ptRows_) {
        // Label filter
        if (!filterLabel.empty() && r.label.find(filterLabel) == std::string::npos)
            continue;
        // Class filter
        if (classFilter == 2 && r.classification != "Parent") continue;
        if (classFilter == 3 && r.classification != "Daughter") continue;
        if (classFilter == 4 && r.classification != "Granddaughter") continue;
        if (classFilter == 5 && r.classification != "Background") continue;
        if (classFilter == 6 && r.classification.substr(0,6) != "Custom") continue;
        if (classFilter == 7 && !r.label.empty()) continue;

        rows.push_back(&r);
    }

    // Sort
    auto sortFn = [&](const PeakTableRow* a, const PeakTableRow* b) -> bool {
        switch (sortSel) {
            case 2: return a->area    < b->area;
            case 3: return a->fwhm    < b->fwhm;
            case 4: return a->chi2ndf < b->chi2ndf;
            case 5: return a->label   < b->label;
            case 6: return a->histIdx < b->histIdx;
            default: return a->energy < b->energy;
        }
    };
    std::sort(rows.begin(), rows.end(), sortFn);

    ptTableView_->Clear();
    // Header
    ptTableView_->AddLine("  Histogram              |  Energy  | ±σ_E |  FWHM  |   Area    | ±σ_A  | Chi2  | Label               | Class");
    ptTableView_->AddLine("  -----------------------|----------|------|--------|-----------|-------|-------|---------------------|----------");

    char buf[320];
    for (const PeakTableRow* r : rows) {
        // Truncate histName to 22 chars
        std::string hname = r->histName;
        if (hname.size() > 22) hname = hname.substr(0, 19) + "...";

        std::string lbl = r->label;
        if (lbl.size() > 20) lbl = lbl.substr(0, 17) + "...";
        std::string cls = r->classification;
        if (cls.size() > 10) cls = cls.substr(0, 10);

        std::snprintf(buf, sizeof(buf),
            "  %-22s | %8.3f | %4.3f | %6.4f | %9.1f | %5.1f | %5.2f | %-20s| %s",
            hname.c_str(),
            r->energy, r->energyErr, r->fwhm,
            r->area, r->areaErr, r->chi2ndf,
            lbl.c_str(), cls.c_str());
        ptTableView_->AddLine(buf);
    }

    if (rows.empty())
        ptTableView_->AddLine("  (no rows match filter)");
    else {
        std::snprintf(buf, sizeof(buf), "  --- %d peaks shown ---", (int)rows.size());
        ptTableView_->AddLine(buf);
    }

    ptTableView_->Update();
}

// ─────────────────────────────────────────────────────────────────────────────
// Slot: OnPTPlot — draw TGraph on main canvas
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnPTPlot() {
    if (!ptXAxisCombo_ || !ptYAxisCombo_) return;

    int xSel = ptXAxisCombo_->GetSelected();
    int ySel = ptYAxisCombo_->GetSelected();

    std::string labelFilter;
    if (ptGraphLabel_) labelFilter = ptGraphLabel_->GetText();

    // Axis labels
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
        switch (ySel) {
            case 1: return r.areaErr;
            case 3: return r.energyErr;
            default: return 0.0;
        }
    };
    auto getXErr = [&](const PeakTableRow& r) -> double {
        if (xSel == 2) return r.energyErr;
        return 0.0;
    };

    // Filter rows
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
// Slot: OnPTExportCSV — export table to CSV
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

    // Header
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
