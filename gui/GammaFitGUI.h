#ifndef GAMMAFITGUI_H
#define GAMMAFITGUI_H

// ROOT GUI
#include "TGFrame.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGNumberEntry.h"
#include "TGListBox.h"
#include "TGComboBox.h"
#include "TGTab.h"
#include "TGTextView.h"
#include "TGStatusBar.h"
#include "TRootEmbeddedCanvas.h"

// Standard
#include <map>
#include <set>
#include <string>
#include <vector>

// Project
#include "GammaDB.h"
#include "ResolutionModel.h"
#include "FitStorage.h"
#include "PeakTracker.h"
#include "FitDatabase.h"
#include "Debug.h"

class TH1;
class TH2;
class TFile;
class TCanvas;
class TF1;
class TGraph;
class TGraphErrors;
class TGTextEntry;

struct CustomProjDef {
    std::string th2Name;
    bool        projX;   // true = ProjectionX (cut on Y); false = ProjectionY (cut on X)
    double      lo, hi;  // cut range on the non-projected axis
};

struct SourceLine {
    double energy;    // keV (reference value)
    double intensity; // branching ratio [0,1]
    int    assigned;  // index into srcPeakEs_; -1 = unassigned
};

// ─────────────────────────────────────────────────────────────────────────────
// GammaFitGUI — Main interactive analysis window.
//
// Left panel (310 px): TGTab with "AutoFit", "Manual Fit", "Fit Results" tabs
// Right panel         : shared spectrum canvas + log strip
// ─────────────────────────────────────────────────────────────────────────────
class GammaFitGUI : public TGMainFrame {
public:
    explicit GammaFitGUI(const TGWindow* p,
                         UInt_t w = 1400, UInt_t h = 920);
    ~GammaFitGUI() override;

    // ── File & histogram ──────────────────────────────────────────────────────
    void OnOpenFile();
    void OnOpenIsotopeDB();
    void OnHistogramSelected(Int_t id);

    // ── AutoFit tab ───────────────────────────────────────────────────────────
    void OnRunSelected();
    void OnRunAll();
    void OnLoadCacheSelected();
    void OnLoadCacheAll();
    void OnApplyBgSelected();
    void OnApplyBgAll();
    void OnResetBgSub();
    void OnDebugAllOn();
    void OnDebugAllOff();

    // ── Source tab ────────────────────────────────────────────────────────────
    void OnOpenSourceRootFile();
    void OnLoadSourceFile();
    void OnRunSourceAutoFit();
    void OnLoadSourceCache();
    void OnAutoIdentify();
    void OnManualAssign();
    void OnShowEnergyCalib();
    void OnShowEfficiency();
    void OnShowSourceFWHM();
    void OnActivityUnitChanged(Int_t id);

    // ── Fit Results tab ───────────────────────────────────────────────────────
    void OnFitResultSelected(Int_t id);
    void OnSaveSelected();
    void OnSaveAll();
    void OnExportCacheCSV();

    // ── Decay tab ─────────────────────────────────────────────────────────────
    void OnDecayTh2Changed(Int_t id);
    void OnDecayModelChanged(Int_t id);
    void OnDecayPeakSelected(Int_t id);
    void OnFitDecay();
    void OnRefreshDecayPeaks();
    void OnPreviewDecay();
    void OnDecayApplyLabel();
    void OnMakePeakCountVsTime();
    void OnLoadDecayCache();

    // ── Histogram classification ──────────────────────────────────────────────
    void OnHistClassSet();

    // ── Custom projections ────────────────────────────────────────────────────
    void OnAddCustomProjection();

    // ── Isotopes tab ──────────────────────────────────────────────────────────
    void OnIsoRefresh();
    void OnIsoFilterChanged(Int_t id);
    void OnIsoListSelected(Int_t id);
    void OnIsoApply();
    void OnIsoApplyLabelAll();
    void OnIsoClear();
    void OnIsoClearAll();
    void OnIsoAutoMatchAll();
    void OnIsoDrawSchematic();
    void OnIsoApplyClassToAll();
    void OnIsoDbSearch();
    void OnIsoDbClear();
    void OnIsoDbLineSelected(Int_t id);
    void OnIsoDbApply();

    // ── FWHM tab ──────────────────────────────────────────────────────────────
    void OnLoadFWHM();
    void OnFitFWHM();
    void OnAcceptFWHM();
    void OnFWHMRestoreAll();
    void OnFWHMRedisplay();
    void OnSaveFWHMCanvas();

    // ── Manual Fit tab ────────────────────────────────────────────────────────
    void OnHistViewChanged(Int_t id);
    void OnLoadManual();
    void OnLoadCache();
    void OnCanvasEvent(Int_t event, Int_t px, Int_t py, TObject* obj);
    void OnSeedParams();
    void OnPreview();
    void OnManualFit();
    void OnAcceptFit();
    void OnRejectFit();
    void OnClearPeaks();
    void OnRemovePeak();
    void OnSetBgFromCanvas();
    void OnFitBackground();
    void OnClearBackground();
    void OnSetRangeFromCanvas();
    void OnClearFitRange();
    void OnParameterScan();

    // ── Peak navigation ───────────────────────────────────────────────────────
    void OnPrevPeak();
    void OnNextPeak();
    void OnZoomIn();
    void OnZoomOut();
    void OnDeleteCacheEntry();

    // ── Residuals ─────────────────────────────────────────────────────────────
    void OnToggleResiduals();
    void OnSelectResidualFit(Int_t id);

    // ── 2D histogram ─────────────────────────────────────────────────────────
    void OnApplyTh2Labels();
    void OnApplySrcTh2Labels();

    // ── Utilities ─────────────────────────────────────────────────────────────
    void AppendLog(const std::string& msg);
    void SetStatus(const std::string& msg);
    void OnSavePlot();
    void OnClearHistCache();
    void OnApplyCanvasAnnotations();

    ClassDefOverride(GammaFitGUI, 0)

private:
    // ── Analysis objects ──────────────────────────────────────────────────────
    GammaDB       db_;
    ResolutionModel res_;
    FitStorage    storage_;
    PeakTracker   tracker_;
    bool          dbLoaded_ = false;

    // ── Run-time state ────────────────────────────────────────────────────────
    TFile*       inputFile_   = nullptr;
    std::string  inputPath_;
    std::vector<std::string> histNames_;
    std::string  currentHist_;
    TH1*         rawHist_    = nullptr;
    bool         rawHistOwned_ = false;   // true when rawHist_ is a projection we created
    TH1*         viewHist_   = nullptr;   // owned bg-sub copy for view modes 2/3
    TF1*         manualTF1_       = nullptr;
    std::vector<TF1*> fitComponents_;  // bg + individual Gaussian components, owned here
    double       clickEnergy_     = 500.0;
    double       lastManualChi2ndf_ = -1.0;  // chi2/ndf from last OnManualFit run

    // Histogram classification: "Gamma" | "Decay" | "2D"
    std::map<std::string, std::string> histClass_;

    // Custom projections — created by the user with explicit cut ranges
    std::map<std::string, CustomProjDef> customProjDefs_;

    // TH2 support — main file
    std::set<std::string>              th2Names_;    // TH2 names in inputFile_
    std::map<std::string, std::string> projParent_;  // "name_px"/"name_py" -> parent TH2 name

    // TH2 support — source file
    std::set<std::string>              srcTh2Names_;
    std::map<std::string, std::string> srcProjParent_;

    // ── Shared widgets ────────────────────────────────────────────────────────
    TGStatusBar*         statusBar_ = nullptr;
    TGTextView*          logView_   = nullptr;
    TRootEmbeddedCanvas* canvas_    = nullptr;

    // ── AutoFit tab widgets ───────────────────────────────────────────────────
    std::string    isotopePath_;
    TGLabel*       fileLbl_          = nullptr;
    TGLabel*       isotopeLbl_       = nullptr;
    TGCheckButton* bgSubtractChk_    = nullptr;
    TGNumberEntry* bgIterEntry_      = nullptr;
    TGNumberEntry* tspecSigmaEntry_  = nullptr;
    TGNumberEntry* tspecThreshEntry_ = nullptr;
    TGListBox*     histList_         = nullptr;
    TGCheckButton* useSeedsChk_      = nullptr;
    TGCheckButton* autoLogLikChk_    = nullptr;
    TGCheckButton* autoImprovChk_    = nullptr;
    TGCheckButton* debugChk_[8]      = {};
    TGComboBox*    histClassCombo_   = nullptr;  // Gamma / Decay / 2D

    // Custom projection widgets
    TGComboBox*    custProjTh2Combo_  = nullptr;
    TGComboBox*    custProjAxisCombo_ = nullptr;  // 1=ProjX(cut Y), 2=ProjY(cut X)
    TGNumberEntry* custProjLo_        = nullptr;
    TGNumberEntry* custProjHi_        = nullptr;
    TGTextEntry*   custProjName_      = nullptr;
    TGTextEntry*   th2XLabelEntry_ = nullptr;  // axis label for TH2 X axis (main file)
    TGTextEntry*   th2YLabelEntry_ = nullptr;  // axis label for TH2 Y axis (main file)

    // ── Isotopes tab widgets ──────────────────────────────────────────────────
    TGListBox*     isoList_        = nullptr;
    TGComboBox*    isoFilterCombo_ = nullptr;
    TGComboBox*    isoLabelCombo_       = nullptr;  // DB-match combo for label
    TGTextEntry*   isoCustomLabelEntry_ = nullptr;  // free-text custom label (overrides combo)
    TGComboBox*    isoClassCombo_       = nullptr;
    TGTextEntry*   isoCustomEntry_      = nullptr;
    // parallel arrays: one entry per isoList_ item
    std::vector<std::string> isoListKeys_;       // FitEntry keys matching each row
    std::vector<std::string> isoListAutoMatch_;  // best DB match isotope (may be empty)
    std::vector<double>      isoListDbEnergy_;   // matched DB line energy (0 if none)
    std::string              isoHistName_;  // histogram whose cache is shown
    // Label → class mapping (label name → class string, e.g. "Co-60" → "Parent")
    std::map<std::string, std::string> labelClassMap_;
    // Isotope DB browser
    TGListBox*     isoDbList_      = nullptr;
    TGTextEntry*   isoDbSearch_    = nullptr;
    TGComboBox*    isoDbClassCombo_ = nullptr;
    std::vector<int> isoDbIndices_;  // index into db_.db for each visible row

    // ── Fit Results tab widgets ───────────────────────────────────────────────
    TGListBox*           fitResultsList_ = nullptr;
    std::vector<std::string> fittedHists_;   // parallel to fitResultsList_ entries

    // Canvas annotation widgets
    TGTextEntry*   frTitleEntry_  = nullptr;
    TGTextEntry*   frXLabelEntry_ = nullptr;
    TGTextEntry*   frYLabelEntry_ = nullptr;
    TGTextEntry*   frAnnotText_   = nullptr;
    TGNumberEntry* frAnnotX_      = nullptr;
    TGNumberEntry* frAnnotY_      = nullptr;
    TGNumberEntry* frAnnotSize_   = nullptr;

    // ── Manual Fit tab widgets ────────────────────────────────────────────────
    TGComboBox*    manualCombo_  = nullptr;
    TGComboBox*    histViewCombo_ = nullptr;  // Raw / BG Sub / BG Sub + AutoFit
    TGNumberEntry* mEnergy_      = nullptr;
    TGNumberEntry* mSigma_       = nullptr;
    TGNumberEntry* mAmp_         = nullptr;
    TGNumberEntry* mBg0_         = nullptr;
    TGNumberEntry* mBg1_         = nullptr;
    TGNumberEntry* mRange_       = nullptr;
    TGLabel*       mResultLbl_   = nullptr;
    TGTextView*    mFitParamsView_ = nullptr;  // per-parameter values + uncertainties
    TGTextView*    peakStatsView_  = nullptr;  // derived statistics per fitted peak

    // Multi-peak
    std::vector<double>  manualPeaks_;
    TGListBox*           peakListBox_  = nullptr;
    TGCheckButton*       addPeakChk_   = nullptr;

    // Peak label & classification (saved with cache entry in OnAcceptFit)
    TGComboBox*          mPeakLabelCombo_ = nullptr;  // populated with DB matches after fit
    TGComboBox*          mPeakClass_      = nullptr;
    TGTextEntry*         mPeakCustom_     = nullptr;

    // Fit options
    TGCheckButton*       mFitLogLikChk_ = nullptr;
    TGCheckButton*       mFitImprovChk_ = nullptr;
    TGCheckButton*       mFitMinosChk_  = nullptr;

    // Background region
    TGNumberEntry*       mBgLo_        = nullptr;
    TGNumberEntry*       mBgHi_        = nullptr;
    bool                 bgClickMode_  = false;
    int                  bgClickCount_ = 0;
    TF1*                 bgTF1_        = nullptr;

    // Fit range (2-click)
    TGNumberEntry*       mFitLo_          = nullptr;
    TGNumberEntry*       mFitHi_          = nullptr;
    bool                 rangeClickMode_  = false;
    int                  rangeClickCount_ = 0;

    // ── Source tab widgets ────────────────────────────────────────────────────
    TGLabel*       srcRootFileLbl_  = nullptr;
    TGComboBox*    srcHistCombo_    = nullptr;
    TGCheckButton* srcBgSubChk_     = nullptr;
    TGTextEntry*   srcTh2XLabelEntry_ = nullptr;  // axis label for TH2 X axis (source)
    TGTextEntry*   srcTh2YLabelEntry_ = nullptr;  // axis label for TH2 Y axis (source)
    TGNumberEntry* srcBgIterEntry_  = nullptr;
    TGLabel*       srcFileLbl_      = nullptr;
    TGLabel*       srcInfoLbl_      = nullptr;
    TGComboBox*    srcActivityUnit_ = nullptr;  // 1=Bq, 2=µCi
    TGNumberEntry* srcLiveTime_     = nullptr;
    TGListBox*     srcLineList_     = nullptr;
    TGNumberEntry* srcManualE_      = nullptr;

    // Source ROOT file (separate from the main inputFile_)
    TFile*                   srcRootFile_  = nullptr;
    std::string              srcRootPath_;
    std::vector<std::string> srcHistNames_;

    // Source description data
    std::string    srcIsotope_;
    double         srcActivityRaw_  = 0.0;
    double         srcActivity_     = 0.0;   // always Bq
    double         srcHalflifeDays_ = 0.0;
    std::string    srcCalDate_;
    std::string    srcMeasDate_;
    std::vector<SourceLine> srcLines_;

    // Peaks extracted from cache for source histogram
    std::string         srcHist_;
    std::vector<double> srcPeakEs_;
    std::vector<double> srcPeakCounts_;
    std::vector<double> srcPeakCountsErr_;

    // ── FWHM tab widgets ─────────────────────────────────────────────────────
    TGComboBox*    fwhmCombo_         = nullptr;
    TGNumberEntry* mFwhmA_            = nullptr;
    TGNumberEntry* mFwhmB_            = nullptr;
    TGNumberEntry* mFwhmC_            = nullptr;
    TGLabel*       fwhmResultLbl_     = nullptr;
    TGCheckButton* fwhmShowSigmaChk_  = nullptr;
    TGCheckButton* fwhmStatLineChk_   = nullptr;
    TGCheckButton* fwhmShowResChk_    = nullptr;
    TGCheckButton* fwhmRemoveModeChk_ = nullptr;
    TF1*           fwhmTF1_           = nullptr;   // owned; clones drawn to canvas
    std::string    fwhmHistName_;
    std::vector<double> fwhmAllX_;       // all FWHM data points (energy)
    std::vector<double> fwhmAllY_;       // all FWHM data points (FWHM value)
    std::vector<bool>   fwhmExcluded_;   // parallel exclusion flags
    double         fwhmChi2Ndf_   = -1.0;
    double         fwhmPValue_    = -1.0;
    double         fwhmResidRMS_  = -1.0;
    int            fwhmNdf_       = 0;

    // ── Decay tab widgets ─────────────────────────────────────────────────────
    TGComboBox*    decayTh2Combo_       = nullptr;  // TH2 selection
    TGComboBox*    decayGammaAxisCombo_ = nullptr;  // 1=X is gamma, 2=Y is gamma
    TGNumberEntry* decaySigRangeEntry_  = nullptr;  // ±N*σ window
    TGListBox*     decayPeakList_       = nullptr;  // fitted peaks from gamma projection
    TGComboBox*    decayModelCombo_     = nullptr;  // 1=Parent, 2=Daughter, 3=GDaughter, 4=BG
    TGNumberEntry* decayThalfP_         = nullptr;
    TGCheckButton* decayFixP_           = nullptr;
    TGNumberEntry* decayThalfD_         = nullptr;
    TGCheckButton* decayFixD_           = nullptr;
    TGNumberEntry* decayThalfG_         = nullptr;
    TGCheckButton* decayFixG_           = nullptr;
    TGTextView*    decayResultView_     = nullptr;
    TGTextEntry*   decayLabelEntry_     = nullptr;  // label for selected peak
    TGComboBox*    decayLabelClassCombo_= nullptr;  // class for selected peak
    TGNumberEntry* decaySliceWidthEntry_= nullptr;  // time-slice width for count vs time plot
    TGNumberEntry* decayFitLo_          = nullptr;  // decay fit range low
    TGNumberEntry* decayFitHi_          = nullptr;  // decay fit range high
    TGCheckButton* decayFitFullRange_   = nullptr;  // use full axis range when checked
    TGNumberEntry* decayRebinEntry_     = nullptr;  // rebin factor (1 = no rebin)
    TGLabel*       decayEquationLabel_  = nullptr;  // shows formula for selected model

    // Decay fit result — stored per peak energy for the active TH2
    struct DecayFitResult {
        double peakE  = 0.0;
        int    model  = 0;
        std::vector<double> params;
        std::vector<double> errors;
        double chi2ndf = -1.0;
        int    status  = -1;
        // provenance / annotation
        std::string histName;        // which TH2 this came from
        double      eMin   = 0.0;   // exact gamma energy cut low edge
        double      eMax   = 0.0;   // exact gamma energy cut high edge
        double      Nsig   = 1.0;   // sigma window used
        std::string label;
        std::string classification;
    };

    // Decay state
    std::string                        decayTh2Name_;
    std::string                        decayGammaProjName_;
    std::vector<double>                decayPeakEs_;
    std::vector<double>                decayPeakSigs_;
    std::vector<std::string>           decayPeakKeys_;
    std::map<double, DecayFitResult>   decayFitStore_;  // keyed by peak energy

    // ── Peak navigation ───────────────────────────────────────────────────────
    std::vector<std::string> peakNavKeys_;
    int              peakNavIdx_  = 0;
    TGLabel*         peakNavLbl_  = nullptr;

    // ── Residuals ─────────────────────────────────────────────────────────────
    TGCheckButton*   residualChk_   = nullptr;
    TGComboBox*      residualCombo_ = nullptr;
    bool             residualsOn_   = false;

    // Zoom window (keV); 0,0 = full range
    double viewXmin_ = 0.0;
    double viewXmax_ = 0.0;

    // ── Build helpers ─────────────────────────────────────────────────────────
    void BuildAutoFitTab    (TGCompositeFrame* parent);
    void BuildSourceTab     (TGCompositeFrame* parent);
    void BuildManualFitTab  (TGCompositeFrame* parent);
    void BuildFWHMTab       (TGCompositeFrame* parent);
    void BuildDecayTab      (TGCompositeFrame* parent);
    void BuildFitResultsTab (TGCompositeFrame* parent);
    void BuildIsotopesTab   (TGCompositeFrame* parent);
    void PopulateHistWidgets();
    void SyncDebugToggles();
    void DrawOnCanvas(TH1* h, TF1* fit = nullptr);
    void DrawWithResiduals(TH1* h, TF1* fit, double xlo, double xhi);
    void DrawFitComponents(TCanvas* c, TF1* f);
    void RedrawCurrent();
    void RunFitOnHistogram(const std::string& hname,
                           TFile* overrideFile = nullptr,
                           const std::vector<double>& forcedSeeds = {});
    std::string CacheFileFor(const std::string& hname) const;
    std::string DecayCacheFileFor() const;
    void        SaveDecayFitCache();
    void        LoadDecayFitCache();

    // Fit Results helpers
    void ShowFitResult(const std::string& hname);
    TH1* MakeBgSubHist(TH1* raw, bool doSubtract = false, int iterations = 14);
    void OverlayFitPeaks(const std::string& hname, TCanvas* c);

    // FWHM helpers
    void RedrawFWHM();
    void DrawFWHMToCanvas(TCanvas* c, bool showSigma, bool showStatLine,
                          bool showResolution = false);

    // Source tab helpers
    double      ComputeDecayedActivity() const;
    void        PopulateSourceList();
    void        ExtractPeaksFromCache(const std::string& hname);
    void        UpdateSourceInfoLabel();
    void        PopulateSrcHistCombo();
    std::string SourceAnalysisFileFor(const std::string& hname) const;
    void        SaveSourceAnalysis();

    // Peak-nav helpers
    void NavigateToPeak(int idx);
    void PopulateNavAndResidual();
    TF1* BuildFromCacheKey(const std::string& key);
    TF1* RebuildFromEntry(const FitEntry& entry, double xlo, double xhi);
    void DrawPeakLabels(TF1* f);
    void SaveFitResultToFile(const std::string& hname, TFile* fout);
    void UpdatePeakStats(TF1* f, TH1* h, double xlo, double xhi);
    void PopulateIsoList(const std::string& filterLabel = "");
    void PopulateIsoDbList(const std::string& filter = "");
    void RefreshIsoDisplay();
    void DrawDecaySchematic(TCanvas* c);
    void LoadLabelClassMap(FitDatabase& fitdb);
    void SaveLabelClassMap(FitDatabase& fitdb);
    void PopulateDecayTh2Combo();
    void PopulateCustProjTh2Combo();

    // Histogram classification helpers
    std::string ClassOf(const std::string& name) const;
    std::string MetadataFileFor() const;
    void        SaveMetadata() const;
    void        LoadMetadata();

    // Returns histogram for hname from inputFile_; creates projection if needed.
    // Sets owned=true when caller must delete the returned pointer.
    TH1* LoadHistFromFile(const std::string& hname, bool& owned) const;
    TH1* LoadProjection(TFile* f, const std::string& projName,
                        const std::map<std::string, std::string>& projParent,
                        const TGTextEntry* xLblEntry,
                        const TGTextEntry* yLblEntry) const;
};

#endif
