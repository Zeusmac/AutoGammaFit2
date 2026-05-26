#include "GammaFitGUI.h"
#include "GammaFitGUI_shared.h"
#include "database/FitDatabase.h"

#include "TGTextEntry.h"
#include "TGFileDialog.h"
#include "TGMsgBox.h"
#include "TH1.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TLatex.h"
#include "TAxis.h"
#include "TROOT.h"
#include "TSystem.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <dirent.h>

#include "TFile.h"
#include "TKey.h"
#include "TClass.h"
#include "NNDCFetcher.h"

// ─────────────────────────────────────────────────────────────────────────────
// Path helpers
// ─────────────────────────────────────────────────────────────────────────────

static std::string ECalibPath(const std::string& cacheDir, const std::string& name)
{
    std::string safe = name;
    for (char& c : safe)
        if (!std::isalnum((unsigned char)c) && c != '_' && c != '-') c = '_';
    return cacheDir + "/energy_cal_" + safe + ".txt";
}

// ─────────────────────────────────────────────────────────────────────────────
// RefreshECalCalibs — populate both apply combos with saved calibrations
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnRefreshECalCalibs()
{
    std::string cdir = CacheDirFor();
    std::vector<std::pair<std::string,std::string>> found; // (name, path)

    DIR* d = opendir(cdir.c_str());
    if (d) {
        struct dirent* de;
        while ((de = readdir(d))) {
            std::string fn = de->d_name;
            static const std::string pfx = "energy_cal_", sfx = ".txt";
            if (fn.size() > pfx.size() + sfx.size()
                    && fn.substr(0, pfx.size()) == pfx
                    && fn.substr(fn.size() - sfx.size()) == sfx) {
                std::string nm = fn.substr(pfx.size(), fn.size() - pfx.size() - sfx.size());
                found.push_back({nm, cdir + "/" + fn});
            }
        }
        closedir(d);
    }
    std::sort(found.begin(), found.end());

    auto fillCombo = [&](TGComboBox* cb) {
        if (!cb) return;
        cb->RemoveAll();
        if (found.empty()) {
            cb->AddEntry("(none saved)", 1);
            cb->Select(1, kFALSE);
        } else {
            int id = 1;
            for (const auto& [nm, path] : found)
                cb->AddEntry(nm.c_str(), id++);
            cb->Select(1, kFALSE);
        }
        cb->MapSubwindows(); cb->Layout();
    };
    fillCombo(ecalApplyCombo_);
    fillCombo(autoFitEcalCombo_);
    fillCombo(effEcalCombo_);
    AppendLog(Form("Energy Cal: found %d saved calibration(s).", (int)found.size()));
}

// ─────────────────────────────────────────────────────────────────────────────
// Load calibration parameters from a named file
// ─────────────────────────────────────────────────────────────────────────────

static bool LoadECalib(const std::string& path, double& a, double& b, double& c)
{
    std::ifstream f(path);
    if (!f.is_open()) return false;
    a = 0.0; b = 1.0; c = 0.0;
    std::string line;
    while (std::getline(f, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string key; double val;
        if (!(ss >> key >> val)) continue;
        if      (key == "a") a = val;
        else if (key == "b") b = val;
        else if (key == "c") c = val;
    }
    return true;
}

// ─────────────────────────────────────────────────────────────────────────────
// GetSelectedECalib — read the selected combo entry and load calibration
// ─────────────────────────────────────────────────────────────────────────────

static bool GetSelectedECalib(TGComboBox* combo, const std::string& cacheDir,
                               double& a, double& b, double& c, std::string& name)
{
    if (!combo) return false;
    TGLBEntry* le = combo->GetSelectedEntry();
    if (!le) return false;
    name = le->GetTitle();
    return LoadECalib(ECalibPath(cacheDir, name), a, b, c);
}

// ─────────────────────────────────────────────────────────────────────────────
// Apply calibration E_new(ch) = a + b*ch + c*ch^2 to a TH1
// Clones the histogram with a rescaled x-axis and injects it into histNames_
// ─────────────────────────────────────────────────────────────────────────────

static TH1* ApplyECalibToHist(TH1* h, double a, double b, double c,
                               const std::string& suffix)
{
    if (!h) return nullptr;
    int nbins = h->GetNbinsX();
    // Build calibrated bin edges
    std::vector<double> edges(nbins + 1);
    for (int i = 0; i <= nbins; i++) {
        double ch = h->GetBinLowEdge(i + 1);
        edges[i] = a + b * ch + c * ch * ch;
    }
    std::string newName  = std::string(h->GetName()) + suffix;
    std::string newTitle = std::string(h->GetTitle()) + "  [cal]";
    TH1F* hcal = new TH1F(newName.c_str(), newTitle.c_str(), nbins, edges.data());
    hcal->SetDirectory(nullptr);
    for (int i = 1; i <= nbins; i++) {
        hcal->SetBinContent(i, h->GetBinContent(i));
        hcal->SetBinError(i, h->GetBinError(i));
    }
    hcal->GetXaxis()->SetTitle("Energy (keV)");
    return hcal;
}

// ─────────────────────────────────────────────────────────────────────────────
// OnAddECalHist — take current srcLines_ assignments and push into ecal data
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnAddECalHist()
{
    if (srcLines_.empty() || srcPeakEs_.empty()) {
        AppendLog("ECal: load source lines and run AutoFit / Load Cache first.");
        return;
    }

    std::vector<std::pair<double,double>> pts; // (fittedE, refE)
    for (const auto& sl : srcLines_) {
        if (sl.assigned < 0 || sl.assigned >= (int)srcPeakEs_.size()) continue;
        pts.push_back({srcPeakEs_[sl.assigned], sl.energy});
    }
    if (pts.empty()) {
        AppendLog("ECal: no assigned peaks — use Auto Identify or Manual Assign first.");
        return;
    }

    const std::string& hname = srcHist_.empty() ? "?" : srcHist_;

    // Reject duplicate source
    for (const auto& s : ecalLoadedHists_) {
        if (s == hname) {
            AppendLog("ECal: " + hname + " is already in the calibration set.");
            return;
        }
    }

    for (const auto& [fE, rE] : pts) {
        ecalAllX_.push_back(fE);
        ecalAllY_.push_back(rE);
        ecalHistSources_.push_back(hname);
    }
    ecalLoadedHists_.push_back(hname);

    // Update point list
    if (ecalPtList_) {
        ecalPtList_->RemoveAll();
        for (size_t i = 0; i < ecalAllX_.size(); i++) {
            std::string e = Form("%.3f -> %.3f keV  [%s]",
                                 ecalAllX_[i], ecalAllY_[i], ecalHistSources_[i].c_str());
            ecalPtList_->AddEntry(e.c_str(), (Int_t)i + 1);
        }
        ecalPtList_->MapSubwindows(); ecalPtList_->Layout();
    }
    if (ecalHistList_) {
        ecalHistList_->RemoveAll();
        int idx = 1;
        for (const auto& h : ecalLoadedHists_)
            ecalHistList_->AddEntry(h.c_str(), idx++);
        ecalHistList_->MapSubwindows(); ecalHistList_->Layout();
    }
    AppendLog(Form("ECal: added %d point(s) from %s  (total %d)",
                   (int)pts.size(), hname.c_str(), (int)ecalAllX_.size()));
}

// ─────────────────────────────────────────────────────────────────────────────
// OnRemoveECalHist — remove selected point from list
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnRemoveECalHist()
{
    if (!ecalPtList_) return;
    Int_t sel = ecalPtList_->GetSelected();
    if (sel < 1 || (size_t)sel > ecalAllX_.size()) return;
    size_t idx = (size_t)(sel - 1);
    std::string src = ecalHistSources_[idx];
    ecalAllX_.erase(ecalAllX_.begin() + idx);
    ecalAllY_.erase(ecalAllY_.begin() + idx);
    ecalHistSources_.erase(ecalHistSources_.begin() + idx);
    // Remove from loaded list if no more points from that source
    if (std::find(ecalHistSources_.begin(), ecalHistSources_.end(), src) == ecalHistSources_.end())
        ecalLoadedHists_.erase(std::find(ecalLoadedHists_.begin(), ecalLoadedHists_.end(), src));

    ecalPtList_->RemoveAll();
    for (size_t i = 0; i < ecalAllX_.size(); i++) {
        std::string e = Form("%.3f -> %.3f keV  [%s]",
                             ecalAllX_[i], ecalAllY_[i], ecalHistSources_[i].c_str());
        ecalPtList_->AddEntry(e.c_str(), (Int_t)i + 1);
    }
    ecalPtList_->MapSubwindows(); ecalPtList_->Layout();
    if (ecalHistList_) {
        ecalHistList_->RemoveAll();
        int id = 1;
        for (const auto& h : ecalLoadedHists_) ecalHistList_->AddEntry(h.c_str(), id++);
        ecalHistList_->MapSubwindows(); ecalHistList_->Layout();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// OnClearECalHists — clear all calibration data
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnClearECalHists()
{
    ecalAllX_.clear(); ecalAllY_.clear();
    ecalHistSources_.clear(); ecalLoadedHists_.clear();
    if (ecalPtList_) { ecalPtList_->RemoveAll(); ecalPtList_->MapSubwindows(); ecalPtList_->Layout(); }
    if (ecalHistList_) { ecalHistList_->RemoveAll(); ecalHistList_->MapSubwindows(); ecalHistList_->Layout(); }
    if (ecalResultLbl_) ecalResultLbl_->SetText("  No fit yet");
    AppendLog("ECal: cleared all calibration points.");
}

// ─────────────────────────────────────────────────────────────────────────────
// OnFitECal — linear fit E_ref = a + b * E_fit
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnFitECal()
{
    int n = (int)ecalAllX_.size();
    if (n < 2) {
        AppendLog("ECal: need at least 2 calibration points to fit.");
        return;
    }

    // Weighted least squares: E_ref = a + b * E_fit
    double sx = 0, sy = 0, sxx = 0, sxy = 0;
    for (int i = 0; i < n; i++) {
        sx  += ecalAllX_[i];
        sy  += ecalAllY_[i];
        sxx += ecalAllX_[i] * ecalAllX_[i];
        sxy += ecalAllX_[i] * ecalAllY_[i];
    }
    double det = n * sxx - sx * sx;
    if (std::abs(det) < 1e-30) {
        AppendLog("ECal: singular matrix — all fitted E values identical?"); return;
    }
    double a = (sy * sxx - sx * sxy) / det;
    double b = (n * sxy - sx * sy) / det;

    // Residuals
    double rms = 0.0;
    for (int i = 0; i < n; i++) {
        double r = ecalAllY_[i] - (a + b * ecalAllX_[i]);
        rms += r * r;
    }
    rms = std::sqrt(rms / (n - 2));

    ecalFitA_ = a; ecalFitB_ = b; ecalFitC_ = 0.0;
    if (ecalA_) ecalA_->SetNumber(a);
    if (ecalB_) ecalB_->SetNumber(b);
    if (ecalC_) ecalC_->SetNumber(0.0);

    std::string res = Form("  a = %.4f   b = %.6f   residual RMS = %.4f keV  (%d pts)",
                           a, b, rms, n);
    if (ecalResultLbl_) ecalResultLbl_->SetText(res.c_str());
    AppendLog("ECal fit: " + res);

    // Draw E_ref vs E_fit with fit line on canvas
    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();
    TGraph* gr = new TGraph(n, ecalAllX_.data(), ecalAllY_.data());
    gr->SetTitle("Energy Calibration;Fitted position;Reference energy (keV)");
    gr->SetMarkerStyle(20); gr->SetMarkerSize(1.0);
    gr->Draw("AP");
    double xlo = *std::min_element(ecalAllX_.begin(), ecalAllX_.end());
    double xhi = *std::max_element(ecalAllX_.begin(), ecalAllX_.end());
    TF1* fl = new TF1("ecal_fit", "[0]+[1]*x", xlo * 0.9, xhi * 1.1);
    fl->SetParameters(a, b);
    fl->SetLineColor(kRed); fl->SetLineWidth(2);
    fl->Draw("same");
    c->Modified(); c->Update();
}

// ─────────────────────────────────────────────────────────────────────────────
// OnAcceptECal — save calibration to file
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnAcceptECal()
{
    if (!ecalNameEntry_) return;
    std::string name = ecalNameEntry_->GetText();
    while (!name.empty() && name.front() == ' ') name = name.substr(1);
    while (!name.empty() && name.back()  == ' ') name.pop_back();
    if (name.empty()) name = "default";

    double a = ecalA_ ? ecalA_->GetNumber() : ecalFitA_;
    double b = ecalB_ ? ecalB_->GetNumber() : ecalFitB_;
    double c = ecalC_ ? ecalC_->GetNumber() : ecalFitC_;

    EnsureCacheDir();
    std::string path = ECalibPath(CacheDirFor(), name);
    std::ofstream f(path);
    if (!f.is_open()) {
        AppendLog("ECal: cannot write " + path); return;
    }
    f << "# Energy calibration: E_ref = a + b*E_fit + c*E_fit^2\n";
    f << "a " << a << "\n";
    f << "b " << b << "\n";
    f << "c " << c << "\n";
    f << "name " << name << "\n";
    f << "points " << ecalAllX_.size() << "\n";
    f.close();

    AppendLog(Form("ECal: saved '%s'  a=%.4f  b=%.6f  c=%.6f  -> %s",
                   name.c_str(), a, b, c, path.c_str()));
    OnRefreshECalCalibs();
}

// ─────────────────────────────────────────────────────────────────────────────
// Apply helpers shared between AutoFit and Source tabs
// ─────────────────────────────────────────────────────────────────────────────

static void InjectCalibratedHist(TH1* hcal, std::vector<std::string>& histNames,
                                  TGComboBox* manualCombo)
{
    if (!hcal) return;
    std::string newName = hcal->GetName();
    if (std::find(histNames.begin(), histNames.end(), newName) == histNames.end()) {
        histNames.push_back(newName);
        if (manualCombo)
            manualCombo->AddEntry(newName.c_str(), (Int_t)histNames.size());
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// OnApplyECalFromAutoFit — apply to selected histogram in AutoFit histList_
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnApplyECalFromAutoFit()
{
    double a, b, c; std::string calName;
    if (!GetSelectedECalib(autoFitEcalCombo_, CacheDirFor(), a, b, c, calName)) {
        AppendLog("ECal: select a saved calibration and click Refresh first."); return;
    }
    Int_t id = histList_ ? histList_->GetSelected() : 0;
    if (id < 1 || (size_t)id > histNames_.size()) {
        AppendLog("ECal: select a histogram in the AutoFit list first."); return;
    }
    const std::string& hname = histNames_[id - 1];
    bool owned = false;
    TH1* h = LoadHistFromFile(hname, owned);
    if (!h) { AppendLog("ECal: cannot load histogram " + hname); return; }

    TH1* hcal = ApplyECalibToHist(h, a, b, c, "_ecal_" + calName);
    if (owned) delete h;
    if (!hcal) return;

    // Store in the open ROOT file so it persists
    if (inputFile_ && !inputFile_->IsZombie()) {
        inputFile_->cd();
        hcal->SetDirectory(static_cast<TDirectory*>(inputFile_));
        hcal->Write("", TObject::kOverwrite);
    } else {
        hcal->SetDirectory(nullptr);
    }
    std::string newName = hcal->GetName();
    InjectCalibratedHist(hcal, histNames_, manualCombo_);
    PopulateHistWidgets();
    AppendLog(Form("ECal: applied '%s' to %s -> %s  (a=%.4f, b=%.6f)",
                   calName.c_str(), hname.c_str(), newName.c_str(), a, b));
    SetStatus("Calibrated: " + newName);
}

// ─────────────────────────────────────────────────────────────────────────────
// OnApplyECalFromSource — apply to currently selected source histogram
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnApplyECalFromSource()
{
    double a, b, c; std::string calName;
    // Try ecalApplyCombo_ first (source tab), fall back to autoFitEcalCombo_
    TGComboBox* cb = ecalApplyCombo_ ? ecalApplyCombo_ : autoFitEcalCombo_;
    if (!GetSelectedECalib(cb, CacheDirFor(), a, b, c, calName)) {
        AppendLog("ECal: select a saved calibration and click Refresh first."); return;
    }
    if (srcHist_.empty()) {
        AppendLog("ECal: select a source histogram first."); return;
    }
    bool owned = false;
    TH1* h = GetSrcHistogram(srcHist_, owned);
    if (!h) { AppendLog("ECal: cannot load source histogram " + srcHist_); return; }

    TH1* hcal = ApplyECalibToHist(h, a, b, c, "_ecal_" + calName);
    if (owned) delete h;
    if (!hcal) return;

    // Save to a new ROOT file alongside the source file
    std::string outPath = CacheDirFor() + "/ecal_" + srcHist_ + "_" + calName + ".root";
    // sanitise
    for (char& ch : outPath)
        if (ch == '/' && &ch != outPath.c_str() + outPath.rfind('/')) {}  // keep last /
    TFile* fout = TFile::Open(outPath.c_str(), "RECREATE");
    if (fout && !fout->IsZombie()) {
        hcal->SetDirectory(static_cast<TDirectory*>(fout));
        hcal->Write("", TObject::kOverwrite);
        fout->Close(); delete fout;
    } else {
        hcal->SetDirectory(nullptr);
    }

    // Inject into source histogram list
    std::string newName = hcal->GetName();
    if (std::find(srcHistNames_.begin(), srcHistNames_.end(), newName) == srcHistNames_.end()) {
        srcHistNames_.push_back(newName);
        SourceHistMeta& m = srcHistMeta_[newName];
        m.isotope   = srcHistMeta_.count(srcHist_) ? srcHistMeta_[srcHist_].isotope : "";
        m.externalFile = outPath;
        m.pathInFile   = newName;
        m.sourceRootFile = srcRootPath_;
        SaveSrcHistMeta();
    }
    PopulateSrcHistCombo();
    PopulateHistWidgets();
    AppendLog(Form("ECal: applied '%s' to source %s -> %s",
                   calName.c_str(), srcHist_.c_str(), newName.c_str()));
    SetStatus("Source calibrated: " + newName);
}

// ─────────────────────────────────────────────────────────────────────────────
// OnApplyECalFromEfficiency — apply cached calibration to source histogram,
// called from the Efficiency sub-tab's "Apply Cached Calibration" button.
// Identical logic to OnApplyECalFromSource but reads from effEcalCombo_.
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnApplyECalFromEfficiency()
{
    double a, b, c; std::string calName;
    TGComboBox* cb = effEcalCombo_ ? effEcalCombo_
                   : (ecalApplyCombo_ ? ecalApplyCombo_ : autoFitEcalCombo_);
    if (!GetSelectedECalib(cb, CacheDirFor(), a, b, c, calName)) {
        AppendLog("ECal: select a saved calibration and click Refresh first."); return;
    }
    if (srcHist_.empty()) {
        AppendLog("ECal: select a source histogram first (Source tab > Spectrum)."); return;
    }
    bool owned = false;
    TH1* h = GetSrcHistogram(srcHist_, owned);
    if (!h) { AppendLog("ECal: cannot load source histogram " + srcHist_); return; }

    TH1* hcal = ApplyECalibToHist(h, a, b, c, "_ecal_" + calName);
    if (owned) delete h;
    if (!hcal) return;

    std::string outPath = CacheDirFor() + "/ecal_" + srcHist_ + "_" + calName + ".root";
    TFile* fout = TFile::Open(outPath.c_str(), "RECREATE");
    if (fout && !fout->IsZombie()) {
        hcal->SetDirectory(static_cast<TDirectory*>(fout));
        hcal->Write("", TObject::kOverwrite);
        fout->Close(); delete fout;
    } else {
        hcal->SetDirectory(nullptr);
    }

    std::string newName = hcal->GetName();
    if (std::find(srcHistNames_.begin(), srcHistNames_.end(), newName) == srcHistNames_.end()) {
        srcHistNames_.push_back(newName);
        SourceHistMeta& m = srcHistMeta_[newName];
        m.isotope        = srcHistMeta_.count(srcHist_) ? srcHistMeta_[srcHist_].isotope : "";
        m.externalFile   = outPath;
        m.pathInFile     = newName;
        m.sourceRootFile = srcRootPath_;
        SaveSrcHistMeta();
    }
    PopulateSrcHistCombo();
    PopulateHistWidgets();
    AppendLog(Form("ECal: applied '%s' to source %s -> %s  (a=%.4f, b=%.6f)",
                   calName.c_str(), srcHist_.c_str(), newName.c_str(), a, b));
    SetStatus("Source calibrated: " + newName);
}

// ─────────────────────────────────────────────────────────────────────────────
// OnLoadFWHMFromSource — add source histogram cache to FWHM combined plot
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnLoadFWHMFromSource()
{
    if (!fwhmSrcCombo_ || srcHistNames_.empty()) {
        AppendLog("FWHM: open a source ROOT file first."); return;
    }
    Int_t id = fwhmSrcCombo_->GetSelected();
    if (id < 1) { AppendLog("FWHM: select a source histogram."); return; }

    // Map combo id to srcHistNames_ (non-TH2 histograms in order)
    std::vector<std::string> nonTh2;
    for (const auto& n : srcHistNames_)
        if (!srcTh2Names_.count(n)) nonTh2.push_back(n);
    if ((size_t)id > nonTh2.size()) { AppendLog("FWHM: invalid selection."); return; }
    const std::string& hname = nonTh2[id - 1];

    // Reject duplicates
    for (const auto& h : fwhmLoadedHists_)
        if (h == hname) { AppendLog("FWHM: " + hname + " already in plot."); return; }

    FitDatabase fitdb;
    if (!fitdb.Load(CacheFileFor(hname))) {
        AppendLog("FWHM: no cache for " + hname + " — run AutoFit or load cache first.");
        return;
    }

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
        AppendLog("FWHM: no fitted peaks in cache for " + hname); return;
    }

    // Restore exclusions saved in this cache
    auto exIt = fitdb.GetEntries().find(kExcludedFwhmKey);
    if (exIt != fitdb.GetEntries().end())
        for (const double exE : exIt->second.params)
            for (size_t i = insertStart; i < fwhmAllX_.size(); i++)
                if (std::abs(fwhmAllX_[i] - exE) < 0.5) { fwhmExcluded_[i] = true; break; }

    fwhmLoadedHists_.push_back(hname);
    fwhmHistName_ = (fwhmLoadedHists_.size() == 1) ? hname : "Combined";

    if (fwhmHistList_) {
        int nPts = 0;
        for (const auto& s : fwhmHistSources_) if (s == hname) ++nPts;
        std::string entry = "[SRC] " + hname + "  (" + std::to_string(nPts) + " pts)";
        fwhmHistList_->AddEntry(entry.c_str(), (Int_t)fwhmLoadedHists_.size());
        fwhmHistList_->MapSubwindows(); fwhmHistList_->Layout();
    }

    int nNew = (int)(fwhmAllX_.size() - insertStart);
    AppendLog(Form("FWHM: added %d points from source %s", nNew, hname.c_str()));
    OnFWHMRedisplay();
}

// ─────────────────────────────────────────────────────────────────────────────
// ECalMatchAndAdd — shared helper: match fitted peaks from a cache to NNDC
// lines for a given isotope, and push matched pairs into the ECal data.
// Returns number of pairs added.
// ─────────────────────────────────────────────────────────────────────────────

static int ECalMatchAndAdd(const std::string& hname,
                            const std::string& isoLabel,
                            const std::string& cacheFile,
                            const std::string& nucCacheDir,
                            std::vector<double>& allX,
                            std::vector<double>& allY,
                            std::vector<std::string>& histSources,
                            double matchTol = 20.0)
{
    if (isoLabel.empty()) return 0;

    FitDatabase fitdb;
    if (!fitdb.Load(cacheFile)) return 0;

    // Collect fitted peak energies from all cache entries
    std::vector<double> fittedEs;
    for (const auto& [key, entry] : fitdb.GetEntries()) {
        if (key.size() >= 2 && key[0] == '_' && key[1] == '_') continue;
        FitLayout lay = DetectLayout((int)entry.params.size());
        if (!lay.valid()) continue;
        for (int i = 0; i < lay.n; i++) {
            double E = entry.params[3*i + 1];
            if (E > 0.0) fittedEs.push_back(E);
        }
    }
    if (fittedEs.empty()) return 0;

    // Load NNDC data for this isotope
    size_t pos = 0;
    while (pos < isoLabel.size() && std::isdigit((unsigned char)isoLabel[pos])) pos++;
    if (pos == 0 || pos >= isoLabel.size()) return 0;
    int A = std::stoi(isoLabel.substr(0, pos));
    std::string sym = isoLabel.substr(pos);

    NucIsotope iso;
    std::string cachePath = nucCacheDir + "/" + std::to_string(A) + sym + ".nucdat";
    if (!NNDCFetcher::LoadCache(cachePath, iso) || iso.gammas.empty()) return 0;

    // Match each NNDC gamma to the nearest fitted peak within tolerance
    int nAdded = 0;
    std::vector<bool> used(fittedEs.size(), false);
    // Sort gammas by intensity descending so bright lines get first pick
    std::vector<const NucGamma*> sorted;
    for (const auto& g : iso.gammas) sorted.push_back(&g);
    std::sort(sorted.begin(), sorted.end(),
              [](const NucGamma* a, const NucGamma* b){ return a->intensity > b->intensity; });

    for (const NucGamma* gp : sorted) {
        double refE = gp->energy;
        if (refE <= 0.0) continue;
        int bestIdx = -1;
        double bestDist = matchTol;
        for (size_t i = 0; i < fittedEs.size(); i++) {
            if (used[i]) continue;
            double d = std::abs(fittedEs[i] - refE);
            if (d < bestDist) { bestDist = d; bestIdx = (int)i; }
        }
        if (bestIdx >= 0) {
            used[bestIdx] = true;
            allX.push_back(fittedEs[bestIdx]);
            allY.push_back(refE);
            histSources.push_back(hname);
            ++nAdded;
        }
    }
    return nAdded;
}

// ─────────────────────────────────────────────────────────────────────────────
// OnAddECalFromSelected — add the source histogram selected in ecalSrcCombo_
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnAddECalFromSelected()
{
    if (!ecalSrcCombo_) return;

    // Map combo id to non-TH2 source histograms (same as PopulateSrcHistCombo)
    Int_t id = ecalSrcCombo_->GetSelected();
    if (id < 1) { AppendLog("ECal: select a source histogram first."); return; }

    std::vector<std::string> nonTh2;
    for (const auto& n : srcHistNames_)
        if (!srcTh2Names_.count(n)) nonTh2.push_back(n);
    if ((size_t)id > nonTh2.size()) { AppendLog("ECal: invalid selection."); return; }
    const std::string& hname = nonTh2[id - 1];

    // Reject duplicates
    for (const auto& h : ecalLoadedHists_)
        if (h == hname) { AppendLog("ECal: " + hname + " already in set."); return; }

    auto mit = srcHistMeta_.find(hname);
    std::string iso = (mit != srcHistMeta_.end()) ? mit->second.isotope : "";

    std::string nucDir = launchDir_ + "/nuclear_cache";
    size_t before = ecalAllX_.size();
    int n = ECalMatchAndAdd(hname, iso, CacheFileFor(hname), nucDir,
                            ecalAllX_, ecalAllY_, ecalHistSources_);

    if (n == 0) {
        AppendLog("ECal: no matched peaks for " + hname +
                  (iso.empty() ? " (no isotope label in metadata)" : " [" + iso + "]"));
        // Roll back any partial additions
        ecalAllX_.resize(before); ecalAllY_.resize(before);
        ecalHistSources_.resize(before);
        return;
    }

    ecalLoadedHists_.push_back(hname);

    // Rebuild point list
    if (ecalPtList_) {
        ecalPtList_->RemoveAll();
        for (size_t i = 0; i < ecalAllX_.size(); i++) {
            std::string e = Form("%.3f -> %.3f keV  [%s]",
                                 ecalAllX_[i], ecalAllY_[i], ecalHistSources_[i].c_str());
            ecalPtList_->AddEntry(e.c_str(), (Int_t)i + 1);
        }
        ecalPtList_->MapSubwindows(); ecalPtList_->Layout();
    }
    if (ecalHistList_) {
        ecalHistList_->RemoveAll();
        int idx = 1;
        for (const auto& h : ecalLoadedHists_)
            ecalHistList_->AddEntry(h.c_str(), idx++);
        ecalHistList_->MapSubwindows(); ecalHistList_->Layout();
    }
    AppendLog(Form("ECal: added %d point(s) from %s [%s]  (total %d)",
                   n, hname.c_str(), iso.c_str(), (int)ecalAllX_.size()));
}

// ─────────────────────────────────────────────────────────────────────────────
// OnAddECalFromRootFile — add all projections from a user-selected .root file
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnAddECalFromRootFile()
{
    static const char* kTypes[] = {
        "ROOT files", "*.root", "All files", "*", nullptr, nullptr
    };
    TGFileInfo fi;
    fi.fFileTypes = kTypes;
    fi.fIniDir    = StrDup(launchDir_.c_str());
    OpenFileDialog(this, kFDOpen, &fi);
    if (!fi.fFilename) return;
    std::string rootPath = fi.fFilename;

    TFile* f = TFile::Open(rootPath.c_str(), "READ");
    if (!f || f->IsZombie()) {
        AppendLog("ECal: cannot open " + rootPath); delete f; return;
    }

    std::string nucDir = launchDir_ + "/nuclear_cache";

    // Collect all TH1 objects recursively (walk TDirectories)
    std::function<void(TDirectory*)> walk = [&](TDirectory* dir) {
        TIter next(dir->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)next())) {
            TClass* cl = gROOT->GetClass(key->GetClassName());
            if (!cl) continue;
            if (cl->InheritsFrom("TDirectory")) {
                TDirectory* sub = (TDirectory*)key->ReadObj();
                if (sub) walk(sub);
            } else if (cl->InheritsFrom("TH1") && !cl->InheritsFrom("TH2")) {
                std::string hname = key->GetName();

                // Look up isotope from metadata (keyed by bare name, ignoring dir prefix)
                std::string iso;
                auto mit = srcHistMeta_.find(hname);
                if (mit != srcHistMeta_.end()) iso = mit->second.isotope;

                if (iso.empty()) {
                    // Try to infer from the directory name (label directory)
                    std::string dirName = dir->GetName();
                    // srcHistMeta_ might key this as "dirName/hname"
                    std::string fullKey = std::string(dir->GetName()) + "/" + hname;
                    auto mit2 = srcHistMeta_.find(fullKey);
                    if (mit2 != srcHistMeta_.end()) iso = mit2->second.isotope;
                }

                if (iso.empty()) {
                    AppendLog("ECal: skipping " + hname + " — no isotope label");
                    continue;
                }

                // Already in set?
                if (std::find(ecalLoadedHists_.begin(), ecalLoadedHists_.end(), hname)
                        != ecalLoadedHists_.end()) continue;

                size_t before = ecalAllX_.size();
                int n = ECalMatchAndAdd(hname, iso, CacheFileFor(hname), nucDir,
                                        ecalAllX_, ecalAllY_, ecalHistSources_);
                if (n > 0) {
                    ecalLoadedHists_.push_back(hname);
                    AppendLog(Form("ECal: +%d pts from %s [%s]", n, hname.c_str(), iso.c_str()));
                } else {
                    ecalAllX_.resize(before); ecalAllY_.resize(before);
                    ecalHistSources_.resize(before);
                }
            }
        }
    };
    walk(f);
    f->Close(); delete f;

    // Rebuild UI lists
    if (ecalPtList_) {
        ecalPtList_->RemoveAll();
        for (size_t i = 0; i < ecalAllX_.size(); i++) {
            std::string e = Form("%.3f -> %.3f keV  [%s]",
                                 ecalAllX_[i], ecalAllY_[i], ecalHistSources_[i].c_str());
            ecalPtList_->AddEntry(e.c_str(), (Int_t)i + 1);
        }
        ecalPtList_->MapSubwindows(); ecalPtList_->Layout();
    }
    if (ecalHistList_) {
        ecalHistList_->RemoveAll();
        int idx = 1;
        for (const auto& h : ecalLoadedHists_)
            ecalHistList_->AddEntry(h.c_str(), idx++);
        ecalHistList_->MapSubwindows(); ecalHistList_->Layout();
    }
    AppendLog(Form("ECal: ROOT file scan done — %d total calibration points from %d histograms.",
                   (int)ecalAllX_.size(), (int)ecalLoadedHists_.size()));
}

// ─────────────────────────────────────────────────────────────────────────────
// OnPlotCalibrationPoints — draw energy vs channel for all enabled calib points,
// skipping any whose source histogram fit entry has needsRefit=true.
// ─────────────────────────────────────────────────────────────────────────────

void GammaFitGUI::OnPlotCalibrationPoints()
{
    const auto& pts = calibBuilder_.GetPoints();
    if (pts.empty()) {
        AppendLog("Calib Plot: no calibration points — add some with 'Add Point' first.");
        return;
    }

    // Cache FitDatabases per source histogram so we only load each file once
    std::map<std::string, FitDatabase> dbCache;
    auto getDB = [&](const std::string& hname) -> FitDatabase* {
        auto it = dbCache.find(hname);
        if (it != dbCache.end()) return &it->second;
        FitDatabase db;
        db.Load(CacheFileFor(hname));
        dbCache[hname] = std::move(db);
        return &dbCache[hname];
    };

    // Collect points to plot
    std::vector<double> xs, ys, exs, eys;
    std::vector<std::string> labels;
    int nSkippedRefit = 0;

    for (const auto& pt : pts) {
        if (!pt.enabled) continue;

        // Check if the source peak is marked for refit
        bool markedRefit = false;
        if (!pt.sourceHistogram.empty()) {
            FitDatabase* db = getDB(pt.sourceHistogram);
            // Scan all entries for the nearest peak to pt.channel
            double bestDist = 5.0; // within 5 keV
            for (const auto& [key, entry] : db->GetEntries()) {
                if (key.size() >= 2 && key[0] == '_' && key[1] == '_') continue;
                if (!entry.needsRefit) continue;
                FitLayout lay = DetectLayout((int)entry.params.size());
                if (!lay.valid()) continue;
                for (int i = 0; i < lay.n; i++) {
                    double peakE = entry.params[3 * i + 1];
                    double d = std::abs(peakE - pt.channel);
                    if (d < bestDist) { bestDist = d; markedRefit = true; }
                }
            }
        }

        if (markedRefit) { ++nSkippedRefit; continue; }

        xs.push_back(pt.channel);
        ys.push_back(pt.energy);
        exs.push_back(pt.channelErr);
        eys.push_back(pt.energyErr);
        labels.push_back(pt.isotope.empty() ? pt.sourceHistogram : pt.isotope);
    }

    if (xs.empty()) {
        AppendLog(Form("Calib Plot: all %d point(s) were excluded (refit or disabled).", (int)pts.size()));
        return;
    }

    // Build TGraphErrors
    TGraphErrors* gr = new TGraphErrors((Int_t)xs.size(),
                                        xs.data(), ys.data(),
                                        exs.data(), eys.data());
    gr->SetName("calibPtGraph");
    gr->SetTitle("Energy Calibration Points;Fitted Position;Reference Energy (keV)");
    gr->SetMarkerStyle(20);
    gr->SetMarkerSize(1.1);
    gr->SetMarkerColor(kBlue + 1);
    gr->SetLineColor(kBlue + 1);

    canvas_->GetCanvas()->cd();
    canvas_->GetCanvas()->Clear();
    gr->Draw("AP");

    // Overlay fitted line if available
    if (calibrationFitted_) {
        TF1* fitLine = currentCalibration_.GetFitFunction("calib_overlay");
        if (fitLine) {
            fitLine->SetLineColor(kRed);
            fitLine->SetLineWidth(2);
            fitLine->Draw("same");
        }
    }

    // Annotate each point with its isotope label
    for (size_t i = 0; i < xs.size(); i++) {
        TLatex* lbl = new TLatex(xs[i], ys[i], labels[i].c_str());
        lbl->SetTextSize(0.022);
        lbl->SetTextColor(kBlue + 1);
        lbl->Draw("same");
    }

    canvas_->GetCanvas()->Update();

    std::string msg = Form("Calib Plot: drew %d point(s)", (int)xs.size());
    if (nSkippedRefit) msg += Form("  (%d skipped — marked for Refit)", nSkippedRefit);
    AppendLog(msg);
}
