#include "GammaFitGUI.h"
#include "GammaFitGUI_shared.h"

#include "TCanvas.h"
#include "TF1.h"
#include "TFitResult.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TAxis.h"
#include "TMath.h"
#include "TSystem.h"
#include "TROOT.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <dirent.h>

// ─────────────────────────────────────────────────────────────────────────────
// EfficiencyCache::eval
// ─────────────────────────────────────────────────────────────────────────────
double EfficiencyCache::eval(double E_keV) const
{
    if (E_keV <= 0.0) return 0.0;
    switch (type) {
    case kLog4: {
        double lnE = std::log(E_keV);
        return std::exp(a - b*lnE + c*lnE*lnE - d/(E_keV*E_keV));
    }
    case kG3LogPoly: {
        if (g3params.size() < 9) return 0.0;
        double lnE = std::log(E_keV);
        double s = 0.0;
        for (int i = 0; i < 9; i++) s += g3params[i] * std::pow(lnE, i);
        return std::exp(s);
    }
    case kStep: {
        if (points.empty()) return 0.0;
        if (E_keV <= points.front().first) return points.front().second;
        if (E_keV >= points.back().first)  return points.back().second;
        auto it = std::lower_bound(points.begin(), points.end(),
                                   std::make_pair(E_keV, 0.0));
        if (it == points.end())    return points.back().second;
        if (it == points.begin())  return points.front().second;
        auto prev = std::prev(it);
        double x0 = prev->first, y0 = prev->second;
        double x1 = it->first,   y1 = it->second;
        double t  = (E_keV - x0) / (x1 - x0);
        return y0 + t*(y1-y0);
    }
    }
    return 0.0;
}

// ─────────────────────────────────────────────────────────────────────────────
// File I/O helpers (file-local)
// ─────────────────────────────────────────────────────────────────────────────

static std::string EffCurvesDir(const std::string& launchDir)
{ return launchDir + "/eff_curves"; }

static std::string EffCurveFilePath(const std::string& launchDir,
                                    const std::string& name)
{ return EffCurvesDir(launchDir) + "/" + name + ".eff"; }

static bool ParseEffFile(const std::string& path, EfficiencyCache& ec)
{
    std::ifstream in(path);
    if (!in.is_open()) return false;
    std::string line;
    bool gotType = false, gotParams = false;
    ec.points.clear(); ec.g3params.clear();
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string key; ss >> key;
        if (key == "type") {
            std::string val; ss >> val;
            if      (val == "log4")      ec.type = EfficiencyCache::kLog4;
            else if (val == "g3logpoly") ec.type = EfficiencyCache::kG3LogPoly;
            else if (val == "step")      ec.type = EfficiencyCache::kStep;
            else return false;
            gotType = true;
        } else if (key == "params") {
            if (ec.type == EfficiencyCache::kLog4) {
                ss >> ec.a >> ec.b >> ec.c >> ec.d; gotParams = true;
            } else if (ec.type == EfficiencyCache::kG3LogPoly) {
                ec.g3params.resize(9, 0.0);
                for (int i = 0; i < 9; i++) ss >> ec.g3params[i];
                gotParams = true;
            }
        } else if (key == "point") {
            double E = 0, eff = 0; ss >> E >> eff;
            ec.points.push_back({E, eff});
        }
    }
    bool ok = gotType && (gotParams || ec.type == EfficiencyCache::kStep);
    if (ok && ec.type == EfficiencyCache::kStep && !ec.points.empty())
        std::sort(ec.points.begin(), ec.points.end());
    return ok;
}

static bool WriteEffFile(const std::string& path, const EfficiencyCache& ec)
{
    std::ofstream out(path);
    if (!out.is_open()) return false;
    out << "# AutoGammaFit2 efficiency curve\n# name " << ec.name << "\n";
    if (ec.type == EfficiencyCache::kLog4) {
        out << "type log4\nparams " << std::setprecision(12)
            << ec.a << " " << ec.b << " " << ec.c << " " << ec.d << "\n";
    } else if (ec.type == EfficiencyCache::kG3LogPoly) {
        out << "type g3logpoly\nparams";
        for (double v : ec.g3params) out << " " << std::setprecision(12) << v;
        out << "\n";
    } else {
        out << "type step\n";
        for (const auto& pt : ec.points)
            out << "point " << std::setprecision(6) << pt.first
                << " " << std::setprecision(10) << pt.second << "\n";
    }
    return true;
}

// Collect source efficiency points (same logic as OnShowEfficiency).
// Returns false if no points found.
static bool CollectEffPoints(
    double A_meas, double liveTime,
    const std::vector<SourceLine>& srcLines,
    const std::vector<double>& srcPeakEs,
    const std::vector<double>& srcPeakCounts,
    const std::vector<double>& srcPeakCountsErr,
    std::vector<double>& outEs,
    std::vector<double>& outEffs,
    std::vector<double>& outErrs)
{
    outEs.clear(); outEffs.clear(); outErrs.clear();
    bool absolute = (A_meas > 0.0 && liveTime > 0.0);
    for (const auto& sl : srcLines) {
        if (sl.assigned < 0 || sl.assigned >= (int)srcPeakEs.size()) continue;
        if (sl.intensity <= 0.0) continue;
        double counts    = srcPeakCounts[sl.assigned];
        double countsErr = srcPeakCountsErr[sl.assigned];
        double denom     = absolute ? (A_meas * sl.intensity * liveTime) : sl.intensity;
        if (denom <= 0.0) continue;
        outEs.push_back(srcPeakEs[sl.assigned]);
        outEffs.push_back(counts / denom);
        outErrs.push_back(countsErr / denom);
    }
    return !outEs.empty();
}

// ─────────────────────────────────────────────────────────────────────────────
// OnEffFitCurve
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnEffFitCurve()
{
    double A_meas   = ComputeDecayedActivity();
    double liveTime = srcLiveTime_ ? srcLiveTime_->GetNumber() : 1.0;

    std::vector<double> Es, effs, errs;
    if (!CollectEffPoints(A_meas, liveTime,
                          srcLines_, srcPeakEs_, srcPeakCounts_, srcPeakCountsErr_,
                          Es, effs, errs)) {
        AppendLog("EffCurve: no assigned source peaks — run AutoIdentify or ManualAssign first.");
        return;
    }
    int n = (int)Es.size();

    if (!effCurveType_) return;
    int typeId = effCurveType_->GetSelected();   // 1=Log4, 2=G3, 3=Step

    effPendingG3_.clear();
    effPendingPoints_.clear();

    if (typeId == 3) {
        for (int i = 0; i < n; i++)
            effPendingPoints_.push_back({Es[i], effs[i]});
        std::sort(effPendingPoints_.begin(), effPendingPoints_.end());
        if (effCurveFitLbl_) {
            effCurveFitLbl_->SetText(Form("  Step: %d point(s) stored", n));
            effCurveFitLbl_->Layout();
        }
        AppendLog(Form("EffCurve: stored %d source points for Step model.", n));
        return;
    }

    std::vector<double> exErr(n, 0.0);
    TGraphErrors* gr = new TGraphErrors(n, Es.data(), effs.data(),
                                         exErr.data(), errs.data());
    gr->SetBit(kCanDelete);

    TF1* tf = nullptr;
    double Elo = Es.front()*0.5, Ehi = Es.back()*2.0;

    if (typeId == 1) {
        // Log4: ln(eff) = a - b*lnE + c*lnE^2 - d/E^2
        tf = new TF1("effLog4",
            "exp([0] - [1]*log(x) + [2]*pow(log(x),2) - [3]/(x*x))",
            Elo, Ehi);
        tf->SetParameters(-3.0, 0.5, 0.05, 100.0);
    } else {
        // G3LogPoly: eff = exp([0]+[1]*lnE+...+[8]*lnE^8)
        tf = new TF1("effG3",
            "exp([0]+[1]*log(x)+[2]*pow(log(x),2)+[3]*pow(log(x),3)"
            "+[4]*pow(log(x),4)+[5]*pow(log(x),5)+[6]*pow(log(x),6)"
            "+[7]*pow(log(x),7)+[8]*pow(log(x),8))",
            Elo, Ehi);
        // Linear seed in log-log space
        double lnEmean = 0, lnEffMean = 0;
        for (int i = 0; i < n; i++) {
            lnEmean   += std::log(Es[i]);
            lnEffMean += std::log(std::max(effs[i], 1e-30));
        }
        lnEmean /= n; lnEffMean /= n;
        double sxy = 0, sxx = 0;
        for (int i = 0; i < n; i++) {
            double dx = std::log(Es[i]) - lnEmean;
            double dy = std::log(std::max(effs[i], 1e-30)) - lnEffMean;
            sxy += dx*dy; sxx += dx*dx;
        }
        double slope = (sxx > 0) ? sxy/sxx : -1.0;
        double intercept = lnEffMean - slope*lnEmean;
        tf->SetParameter(0, intercept);
        tf->SetParameter(1, slope);
        for (int i = 2; i <= 8; i++) tf->SetParameter(i, 0.0);
    }
    tf->SetBit(kCanDelete);

    TFitResultPtr res = gr->Fit(tf, "RQNS", "", Es.front()*0.9, Es.back()*1.1);

    double chi2ndf = -1.0;
    if (res.Get()) {
        double ndf = (double)res->Ndf();
        if (ndf > 0) chi2ndf = res->Chi2() / ndf;
    }

    if (typeId == 1) {
        effPendingG3_ = { tf->GetParameter(0), tf->GetParameter(1),
                          tf->GetParameter(2), tf->GetParameter(3) };
    } else {
        effPendingG3_.resize(9);
        for (int i = 0; i < 9; i++) effPendingG3_[i] = tf->GetParameter(i);
    }

    std::string chi2str = (chi2ndf >= 0) ? Form("  chi2/ndf=%.3f", chi2ndf) : "";
    std::string msg = Form("EffCurve: fit done (%d pts)%s", n, chi2str.c_str());
    if (effCurveFitLbl_) {
        effCurveFitLbl_->SetText(msg.c_str());
        effCurveFitLbl_->Layout();
    }
    AppendLog(msg);

    delete gr;
    delete tf;
}

// ─────────────────────────────────────────────────────────────────────────────
// OnEffSaveCurve
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnEffSaveCurve()
{
    if (!effCurveName_) return;
    std::string name = effCurveName_->GetText();
    if (name.empty()) { AppendLog("EffCurve: enter a name before saving."); return; }

    if (!effCurveType_) return;
    int typeId = effCurveType_->GetSelected();

    if (effPendingG3_.empty() && effPendingPoints_.empty()) {
        AppendLog("EffCurve: run 'Fit to Source Data' first.");
        return;
    }

    std::string dir  = EffCurvesDir(launchDir_);
    gSystem->mkdir(dir.c_str(), true);
    std::string path = EffCurveFilePath(launchDir_, name);

    EfficiencyCache ec;
    ec.name = name;
    if (typeId == 3) {
        ec.type   = EfficiencyCache::kStep;
        ec.points = effPendingPoints_;
    } else if (typeId == 1) {
        ec.type = EfficiencyCache::kLog4;
        if (effPendingG3_.size() >= 4) {
            ec.a = effPendingG3_[0]; ec.b = effPendingG3_[1];
            ec.c = effPendingG3_[2]; ec.d = effPendingG3_[3];
        }
    } else {
        ec.type     = EfficiencyCache::kG3LogPoly;
        ec.g3params = effPendingG3_;
        if ((int)ec.g3params.size() < 9) ec.g3params.resize(9, 0.0);
    }

    if (!WriteEffFile(path, ec)) {
        AppendLog("EffCurve: failed to write " + path); return;
    }
    AppendLog("EffCurve: saved to " + path);
    OnEffScanCurves();
}

// ─────────────────────────────────────────────────────────────────────────────
// OnEffScanCurves
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnEffScanCurves()
{
    effCurves_.clear();
    effCurveSelected_ = -1;
    if (effCurveList_) { effCurveList_->RemoveAll(); effCurveList_->Layout(); }

    std::string dir = EffCurvesDir(launchDir_);
    DIR* dp = opendir(dir.c_str());
    if (!dp) {
        AppendLog("EffCurve: no eff_curves/ directory (will be created on first save).");
        return;
    }

    std::vector<std::string> files;
    struct dirent* ep;
    while ((ep = readdir(dp))) {
        std::string fn = ep->d_name;
        if (fn.size() > 4 && fn.substr(fn.size()-4) == ".eff")
            files.push_back(fn);
    }
    closedir(dp);
    std::sort(files.begin(), files.end());

    int idx = 1;
    for (const auto& fn : files) {
        EfficiencyCache ec;
        ec.name = fn.substr(0, fn.size()-4);
        if (!ParseEffFile(dir + "/" + fn, ec)) {
            AppendLog("EffCurve: skipping malformed " + fn); continue;
        }
        ec.name = fn.substr(0, fn.size()-4);
        effCurves_.push_back(ec);
        static const char* ts[] = { "log4", "g3logpoly", "step" };
        std::string entry = ec.name + "  [" + ts[(int)ec.type] + "]";
        if (ec.type == EfficiencyCache::kStep)
            entry += Form("  %d pts", (int)ec.points.size());
        if (effCurveList_) effCurveList_->AddEntry(entry.c_str(), idx++);
    }
    if (effCurveList_) { effCurveList_->MapSubwindows(); effCurveList_->Layout(); }
    AppendLog(Form("EffCurve: %d curve(s) in %s", (int)effCurves_.size(), dir.c_str()));
}

// ─────────────────────────────────────────────────────────────────────────────
// OnEffCurveSelected
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnEffCurveSelected(Int_t id)
{
    int ci = id - 1;
    if (ci < 0 || ci >= (int)effCurves_.size()) return;
    effCurveSelected_ = ci;
    const EfficiencyCache& ec = effCurves_[ci];

    static const char* ts[] = { "Log4", "G3LogPoly", "Step" };
    std::string info = ec.name + "  [" + ts[(int)ec.type] + "]";
    if (ec.type == EfficiencyCache::kStep) {
        info += Form("  %d pts", (int)ec.points.size());
        if (!ec.points.empty())
            info += Form("  E=[%.0f,%.0f] keV",
                         ec.points.front().first, ec.points.back().first);
    } else if (ec.type == EfficiencyCache::kLog4) {
        info += Form("  a=%.3g b=%.3g c=%.3g d=%.3g", ec.a, ec.b, ec.c, ec.d);
        if (effA_) effA_->SetNumber(ec.a);
        if (effB_) effB_->SetNumber(ec.b);
        if (effC_) effC_->SetNumber(ec.c);
        if (effD_) effD_->SetNumber(ec.d);
    }
    if (effCurveFitLbl_) {
        effCurveFitLbl_->SetText(info.c_str());
        effCurveFitLbl_->Layout();
    }
    AppendLog("EffCurve: selected '" + ec.name + "'");
}

// ─────────────────────────────────────────────────────────────────────────────
// OnEffPlotCurve
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnEffPlotCurve()
{
    // Determine which curve to use
    EfficiencyCache tmpCurve;
    const EfficiencyCache* curvePtr = nullptr;

    if (effCurveSelected_ >= 0 && effCurveSelected_ < (int)effCurves_.size()) {
        curvePtr = &effCurves_[effCurveSelected_];
    } else if (!effPendingG3_.empty() || !effPendingPoints_.empty()) {
        tmpCurve = EfficiencyCache{};
        tmpCurve.name = "(pending)";
        int typeId = effCurveType_ ? effCurveType_->GetSelected() : 2;
        if (!effPendingPoints_.empty()) {
            tmpCurve.type   = EfficiencyCache::kStep;
            tmpCurve.points = effPendingPoints_;
        } else if (typeId == 1) {
            tmpCurve.type = EfficiencyCache::kLog4;
            if ((int)effPendingG3_.size() >= 4) {
                tmpCurve.a = effPendingG3_[0]; tmpCurve.b = effPendingG3_[1];
                tmpCurve.c = effPendingG3_[2]; tmpCurve.d = effPendingG3_[3];
            }
        } else {
            tmpCurve.type     = EfficiencyCache::kG3LogPoly;
            tmpCurve.g3params = effPendingG3_;
            if ((int)tmpCurve.g3params.size() < 9) tmpCurve.g3params.resize(9, 0.0);
        }
        curvePtr = &tmpCurve;
    } else {
        AppendLog("EffCurve: nothing to plot — select a curve or run Fit first.");
        return;
    }

    // Collect source data for overlay
    double A_meas   = ComputeDecayedActivity();
    double liveTime = srcLiveTime_ ? srcLiveTime_->GetNumber() : 1.0;
    std::vector<double> Es, effs, errs, exErr;
    CollectEffPoints(A_meas, liveTime,
                     srcLines_, srcPeakEs_, srcPeakCounts_, srcPeakCountsErr_,
                     Es, effs, errs);
    int npts = (int)Es.size();

    // Plot range
    double E_lo = 50.0, E_hi = 3000.0;
    if (curvePtr->type == EfficiencyCache::kStep && !curvePtr->points.empty()) {
        E_lo = curvePtr->points.front().first * 0.8;
        E_hi = curvePtr->points.back().first  * 1.2;
    }
    if (npts > 0) {
        E_lo = std::min(E_lo, *std::min_element(Es.begin(), Es.end()) * 0.7);
        E_hi = std::max(E_hi, *std::max_element(Es.begin(), Es.end()) * 1.3);
    }

    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();
    c->SetLogx(1); c->SetLogy(1);

    bool drawnFirst = false;
    if (npts > 0) {
        exErr.assign(npts, 0.0);
        TGraphErrors* gr = new TGraphErrors(npts, Es.data(), effs.data(),
                                             exErr.data(), errs.data());
        gr->SetBit(kCanDelete);
        gr->SetTitle(Form("Efficiency  -  %s;Energy (keV);Efficiency",
                          curvePtr->name.c_str()));
        gr->SetMarkerStyle(20);
        gr->SetMarkerSize(1.1);
        gr->SetMarkerColor(kBlue+1);
        gr->SetLineColor(kBlue+1);
        gr->Draw("AP");
        drawnFirst = true;
    }

    // Sample curve
    int nCurve = 400;
    std::vector<double> cE(nCurve), cEff(nCurve);
    double logLo = std::log10(std::max(E_lo, 1.0));
    double logHi = std::log10(std::max(E_hi, 10.0));
    bool allPos = true;
    for (int i = 0; i < nCurve; i++) {
        cE[i]   = std::pow(10.0, logLo + i*(logHi-logLo)/(nCurve-1));
        cEff[i] = curvePtr->eval(cE[i]);
        if (cEff[i] <= 0) allPos = false;
    }
    TGraph* gc = new TGraph(nCurve, cE.data(), cEff.data());
    gc->SetBit(kCanDelete);
    gc->SetLineColor(kRed+1);
    gc->SetLineWidth(2);
    if (!drawnFirst) {
        gc->SetTitle(Form("Efficiency  -  %s;Energy (keV);Efficiency",
                          curvePtr->name.c_str()));
        gc->Draw("AL");
    } else {
        gc->Draw("L same");
    }
    c->Modified(); c->Update();

    AppendLog("EffCurve: plotted '" + curvePtr->name + "'");
    if (!allPos)
        AppendLog("  (warning: curve has zero/negative values — check parameters)");
}

// ─────────────────────────────────────────────────────────────────────────────
// OnEffDeleteCurve
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnEffDeleteCurve()
{
    if (effCurveSelected_ < 0 || effCurveSelected_ >= (int)effCurves_.size()) {
        AppendLog("EffCurve: select a curve first."); return;
    }
    std::string name = effCurves_[effCurveSelected_].name;
    std::string path = EffCurveFilePath(launchDir_, name);
    if (gSystem->Unlink(path.c_str()) != 0) {
        AppendLog("EffCurve: could not delete " + path); return;
    }
    AppendLog("EffCurve: deleted " + path);
    OnEffScanCurves();
}

// ─────────────────────────────────────────────────────────────────────────────
// OnEffApplyToHist
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnEffApplyToHist()
{
    if (effCurveSelected_ < 0 || effCurveSelected_ >= (int)effCurves_.size()) {
        AppendLog("EffCurve: select a curve in the list first."); return;
    }
    const EfficiencyCache& ec = effCurves_[effCurveSelected_];

    TH1* src = rawHist_;
    if (!src) {
        AppendLog("EffCurve: no histogram loaded — open a file and select a histogram.");
        return;
    }

    std::string newName = std::string(src->GetName()) + "_effcorr";
    TObject* old = gROOT->FindObject(newName.c_str());
    if (old) delete old;

    TH1* h = (TH1*)src->Clone(newName.c_str());
    h->SetTitle(Form("%s  (eff-corrected: %s)", src->GetTitle(), ec.name.c_str()));
    h->SetDirectory(nullptr);

    int zeroCount = 0;
    for (int b = 1; b <= h->GetNbinsX(); b++) {
        double E   = h->GetXaxis()->GetBinCenter(b);
        double eff = ec.eval(E);
        if (eff <= 0.0) { ++zeroCount; continue; }
        h->SetBinContent(b, h->GetBinContent(b) / eff);
        h->SetBinError  (b, h->GetBinError(b)   / eff);
    }

    TCanvas* c = canvas_->GetCanvas();
    c->Clear(); c->cd();
    h->SetBit(kCanDelete);
    h->Draw("hist E");
    c->Modified(); c->Update();

    std::string msg = "EffCurve: applied '" + ec.name + "' -> '" + newName + "'";
    if (zeroCount > 0)
        msg += Form("  (%d bins skipped, eff=0)", zeroCount);
    AppendLog(msg);
}

// ─────────────────────────────────────────────────────────────────────────────
// OnEffAddStepPoint
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnEffAddStepPoint()
{
    if (!effManualStepE_ || !effManualStepV_) return;
    double E   = effManualStepE_->GetNumber();
    double eff = effManualStepV_->GetNumber();
    if (E <= 0) { AppendLog("EffCurve: E must be > 0 keV."); return; }
    if (eff < 0){ AppendLog("EffCurve: efficiency must be >= 0."); return; }

    // Remove any existing entry at the same E (±0.001 keV)
    effManualPoints_.erase(
        std::remove_if(effManualPoints_.begin(), effManualPoints_.end(),
                       [E](const std::pair<double,double>& p){
                           return std::abs(p.first - E) < 0.001; }),
        effManualPoints_.end());
    effManualPoints_.push_back({E, eff});
    std::sort(effManualPoints_.begin(), effManualPoints_.end());

    // Rebuild listbox
    if (effManualStepList_) {
        effManualStepList_->RemoveAll();
        for (int i = 0; i < (int)effManualPoints_.size(); i++) {
            std::string entry = Form("%.4g keV  eff=%.6g",
                                     effManualPoints_[i].first,
                                     effManualPoints_[i].second);
            effManualStepList_->AddEntry(entry.c_str(), i+1);
        }
        effManualStepList_->MapSubwindows();
        effManualStepList_->Layout();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// OnEffRemoveStepPoint
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnEffRemoveStepPoint()
{
    if (!effManualStepList_) return;
    Int_t sel = effManualStepList_->GetSelected();
    if (sel < 1 || sel > (int)effManualPoints_.size()) return;
    effManualPoints_.erase(effManualPoints_.begin() + (sel - 1));

    effManualStepList_->RemoveAll();
    for (int i = 0; i < (int)effManualPoints_.size(); i++) {
        std::string entry = Form("%.4g keV  eff=%.6g",
                                  effManualPoints_[i].first,
                                  effManualPoints_[i].second);
        effManualStepList_->AddEntry(entry.c_str(), i+1);
    }
    effManualStepList_->MapSubwindows();
    effManualStepList_->Layout();
}

// ─────────────────────────────────────────────────────────────────────────────
// OnEffClearStepPoints
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnEffClearStepPoints()
{
    effManualPoints_.clear();
    if (effManualStepList_) {
        effManualStepList_->RemoveAll();
        effManualStepList_->Layout();
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// OnEffSaveManual
// ─────────────────────────────────────────────────────────────────────────────
void GammaFitGUI::OnEffSaveManual()
{
    if (!effCurveName_ || !effCurveType_) return;
    std::string name = effCurveName_->GetText();
    if (name.empty()) { AppendLog("EffCurve: enter a name in the Name field first."); return; }

    int typeId = effCurveType_->GetSelected();  // 1=Log4, 2=G3, 3=Step

    EfficiencyCache ec;
    ec.name = name;

    if (typeId == 1) {
        // Log4: read a, b, c, d from manual entries
        ec.type = EfficiencyCache::kLog4;
        ec.a = effManualA_ ? effManualA_->GetNumber() : 0.0;
        ec.b = effManualB_ ? effManualB_->GetNumber() : 0.0;
        ec.c = effManualC_ ? effManualC_->GetNumber() : 0.0;
        ec.d = effManualD_ ? effManualD_->GetNumber() : 0.0;
    } else if (typeId == 2) {
        // G3: parse the text field
        ec.type = EfficiencyCache::kG3LogPoly;
        ec.g3params.assign(9, 0.0);
        if (effManualG3Entry_) {
            std::istringstream ss(effManualG3Entry_->GetText());
            for (int i = 0; i < 9; i++) {
                double v = 0.0;
                if (!(ss >> v)) break;
                ec.g3params[i] = v;
            }
        }
        // Quick sanity check: at least one non-zero param
        bool anyNonZero = false;
        for (double v : ec.g3params) if (v != 0.0) { anyNonZero = true; break; }
        if (!anyNonZero) {
            AppendLog("EffCurve: all G3 params are zero — enter values in the g0..g8 field.");
            return;
        }
    } else {
        // Step
        ec.type = EfficiencyCache::kStep;
        if (effManualPoints_.empty()) {
            AppendLog("EffCurve: no step points entered — use Add to build the list.");
            return;
        }
        ec.points = effManualPoints_;
    }

    std::string dir  = EffCurvesDir(launchDir_);
    gSystem->mkdir(dir.c_str(), true);
    std::string path = EffCurveFilePath(launchDir_, name);

    if (!WriteEffFile(path, ec)) {
        AppendLog("EffCurve: failed to write " + path); return;
    }
    AppendLog("EffCurve: saved manual entry to " + path);
    OnEffScanCurves();
}
