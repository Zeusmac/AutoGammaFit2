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
#include "SourceFitDatabase.h"
#include "CalibrationModel.h"
#include "EnergyCalibration.h"
#include "Debug.h"
#include "NuclearData.h"
#include "OptionalDeps.h"
#include "LevelScheme.h"

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

// Keyboard shortcut descriptor  -  one per configurable action.
struct KbShortcut {
    std::string   id;           // unique action ID, e.g. "choose_peaks"
    std::string   label;        // human-readable label shown in dialog
    UInt_t        keysym     = 0;  // current binding (0 = unbound)
    UInt_t        modifier   = 0;  // current modifier mask (kKeyControlMask etc.)
    UInt_t        defKeysym  = 0;  // factory default
    UInt_t        defModifier = 0;
    TGTextButton* dlgBtn   = nullptr;  // binding button in the shortcuts dialog
    TGTextButton* clearBtn = nullptr;  // clear button in the shortcuts dialog
};

struct BgSubtractDef {
    std::string srcName;  // histogram being subtracted from
    std::string bgName;   // background histogram to subtract
    double      scale;    // scale factor applied to the background (default 1.0)
};

struct SourceLine {
    double energy;    // keV (reference value)
    double intensity; // branching ratio [0,1]
    int    assigned;  // index into srcPeakEs_; -1 = unassigned
};

struct SrcMatchedLine {
    double refE;       // NNDC reference energy (keV)
    double intensity;  // relative intensity [0,1]
    double fittedE;    // fitted peak energy (keV), or -1 if unmatched
    double counts;
    double countsErr;
};

struct PeakTableRow {
    std::string histName;
    int         histIdx      = 0;
    double      energy       = 0.0;
    double      energyErr    = 0.0;
    double      sigma        = 0.0;
    double      fwhm         = 0.0;
    double      area         = 0.0;
    double      areaErr      = 0.0;
    double      chi2ndf      = -1.0;
    double      intensity    = 0.0;  // absolute emission probability (0 = not computed)
    double      intensityErr = 0.0;
    std::string label;
    std::string classification;
    std::string cacheFile;
    bool        needsRefit     = false;
};

// Named efficiency fit / curve.
// type kLog4:      ln(eff) = a - b*ln(E) + c*ln(E)^2 - d/E^2  (legacy 4-param model)
// type kG3LogPoly: eff = exp(sum(g3[i]*lnE^i, i=0..8))  (log-log polynomial, 9 params)
// type kStep:      piecewise-linear interpolation on (E, eff) pairs
struct EfficiencyCache {
    enum Type { kLog4 = 0, kG3LogPoly = 1, kStep = 2 };
    std::string name;
    Type type = kLog4;
    // kLog4 fields (kept for backward compat with PeakTable)
    double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
    // kG3LogPoly: 9 coefficients, g3[0]+g3[1]*lnE+...+g3[8]*lnE^8 in exponent
    std::vector<double> g3params;
    // kStep: sorted (E_keV, efficiency) pairs
    std::vector<std::pair<double,double>> points;

    double eval(double E_keV) const;
};

// ─────────────────────────────────────────────────────────────────────────────
// GammaFitGUI  -  Main interactive analysis window.
//
// Left panel (310 px): TGTab with "AutoFit", "Manual Fit", "Fit Results" tabs
// Right panel         : shared spectrum canvas + log strip
// ─────────────────────────────────────────────────────────────────────────────
class GammaFitGUI : public TGMainFrame {
public:
    explicit GammaFitGUI(const TGWindow* p,
                         UInt_t w = 1400, UInt_t h = 920);
    ~GammaFitGUI() override;

    Bool_t HandleKey(Event_t* event) override;

    // ── File & histogram ──────────────────────────────────────────────────────
    void OnOpenFile();
    void OnOpenIsotopeDB();
    void OnReloadIsotopeDB();
    void OnHistogramSelected(Int_t id);

    // ── AutoFit tab ───────────────────────────────────────────────────────────
    void OnRunSelected();
    void OnRunAll();
    void OnLoadCacheSelected();
    void OnLoadCacheAll();
    void OnApplyBgSelected();
    void OnApplyBgAll();
    void OnResetBgSub();
    void OnApplyBgToCurrent();
    void OnToggleShowBgLine();
    void OnToggleNegBinsRed();
    void OnTogglePeakClickZoom();
    void OnDeleteHistogram();
    void OnDebugAllOn();
    void OnDebugAllOff();

    // ── Source tab ────────────────────────────────────────────────────────────
    void OnOpenSourceRootFile();
    void OnLoadSourceFile();
    void OnRunSourceAutoFit();
    void OnLoadSourceCache();
    void OnAutoIdentify();
    void OnManualAssign();
    void RefreshSrcPeakCombo();
    void OnShowEnergyCalib();
    // ── Energy calibration tab (Source) ───────────────────────────────────────
    void OnAddECalHist();           // legacy helper — still used by OnAddECalFromSelected
    void OnAddECalFromSelected();   // add single projection via ecalSrcCombo_
    void OnAddECalFromRootFile();   // add all projections from a .root file
    void OnRemoveECalHist();
    void OnClearECalHists();
    void OnFitECal();
    void OnAcceptECal();
    void OnRefreshECalCalibs();
    void OnApplyECalFromAutoFit();
    void OnApplyECalFromSource();
    void OnApplyECalFromEfficiency();   // apply cache cal to source hist, from Efficiency sub-tab
    // ── FWHM source cache ──────────────────────────────────────────────────────
    void OnLoadFWHMFromSource();
    void OnShowEfficiency();
    void OnShowSourceFWHM();
    void OnActivityUnitChanged(Int_t id);
    void OnSrcDetTh2Changed(Int_t id);
    void OnSrcExtractDetector();
    void OnSrcTh2ListSelected(Int_t id);
    void OnSrcSetTh2Label();
    void OnSrcAutoDetectLabels();
    void OnSrcFetchNNDC();
    void OnSrcLoadLinesForTh2();
    void OnSrcIsoComboChanged(Int_t id);
    void OnSrcLineSearch();
    void PopulateSourceIsoCombo();
    void PopulateSrcIsoLines(const std::string& isoID, const std::string& filterKeV = "");
    void OnSrcAutoMatchAll();
    void ApplyAutoMatchToCurrentHist();
    void OnSrcClearHistCache();
    void OnSrcClearFileCache();
    void OnSrcAutoFitAllProjections();
    void OnSrcAutoProjectAll();
    void OnSrcAutoProjectInPlace();   // project all TH2 detectors into TDirs in the source file itself
    void OnSrcSendToFitResults();      // add current source histogram to the Fit Results list
    void OnSrcOpenCacheProjections();  // open auto-generated cache projections file
    void OnOpenRecentSrc();                 // open a source file from the recent-files combo
    void OnSrcHistSearch();                 // filter srcHistCombo_ by search text
    void OnSrcHistComboChanged(Int_t id);   // histogram selected in srcHistCombo_
    void OnSrcHistMetaSave();               // persist metadata for current histogram
    void OnSrcPreviewTh2FromList();         // preview TH2 selected in srcTh2List_
    void OnSrcCalcExpected();               // compute & display expected counts with equations

    // ── Source 2D histogram preview and calibration ────────────────────────────
    void OnSrcPreviewTh2();          // display/preview selected TH2
    void OnSrcRefreshTh2Preview();   // refresh preview after selection change
    void OnSrcToggleTh2Projection(); // toggle between 2D and 1D projection view
    void OnSrcBuildCalibrationPlot();     // build energy vs channel from selected sources
    void OnSrcFitCalibration();           // fit linear calibration to selected points
    void OnSrcApplyCalibrationToHist();   // apply fitted calibration to any histogram
    void OnSrcClearCalibration();         // clear current calibration data
    void OnSrcSaveCalibration();          // save calibration to file
    void OnSrcLoadCalibration();          // load calibration from file
    void OnSrcCalibrateHist();            // apply calibration interactively to selected hist
    void OnSrcRemoveCalibPoint(Int_t id); // remove point from calibration
    void OnSrcToggleCalibPoint(Int_t id); // enable/disable calibration point
    void OnPlotCalibrationPoints();        // plot energy vs channel graph

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
    void OnDecayPreviewGammaPeak();
    void OnPreviewDecay();
    void OnDecayAutoModel();
    void OnDecaySeedHalfLives();
    void OnDecayBGTypeChanged(Int_t id);
    void OnDecayRebinReset();
    void OnDecayApplyLabel();
    void OnMakePeakCountVsTime();
    void OnLoadDecayCache();
    void OnSaveDecayCache();
    void OnDecayScanCaches();
    void OnDecayCutChanged();
    void OnDecayErrBarsToggled();
    void OnDecayCacheBrowserSelected(Int_t id);
    void PopulateCacheBrowser();
    void OnAddExtraDecayCut();
    void OnRemoveExtraDecayCut();
    void OnClearExtraDecayCuts();
    // ── Total Decay sub-tab ───────────────────────────────────────────────────
    void OnDecayTdModelChanged(Int_t id);
    void OnDecayTdBGTypeChanged(Int_t id);
    void OnFitTotalDecay();
    void OnFitRK4Decay();
    void OnRK4SelectAll();
    void OnRK4SelectNone();
    void OnRK4SeedFromCuts();
    void OnRK4ScanHistograms();
    void OnChiScanDecay();
    void OnChiScanTotalDecay();
    void OnChiScanRK4();
    void UpdateCutsEquation();
    void UpdateTotalDecayEquation();

    // ── Histogram classification ──────────────────────────────────────────────
    void OnHistClassSet();

    // ── Custom projections ────────────────────────────────────────────────────
    void OnAddCustomProjection();
    void OnCustProjRangeChanged();

    // ── Background histogram subtraction ─────────────────────────────────────
    void OnSubtractHistogram();

    // ── Rebin histogram ───────────────────────────────────────────────────────
    void OnRebinPreview();
    void OnRebinApply();
    void OnRebinReset();

    // ── Peak label display ────────────────────────────────────────────────────
    void OnToggleIsoLabels();
    void OnTogglePeakFWHMArea();

    // ── gnuScope export ───────────────────────────────────────────────────────
    void OnExportGnuScope();
    void OnExportGnuScopeAll();

    // ── Extended fit options (Manual Fit) ─────────────────────────────────────
    void OnApplyBgAnchors();
    void OnClearBgAnchors();
    void OnBgIterChanged();
    void OnImportAllCachesFromFile();
    void OnToggleErrorBars();
    void OnResModelComboChanged();
    void OnLoadResFromHist();

    // ── Channel->keV calibration ───────────────────────────────────────────────
    void OnApplyCalibration();

    // ── Efficiency correction ─────────────────────────────────────────────────
    void OnApplyEfficiency();

    // ── Efficiency curve cache (Efficiency sub-tab) ───────────────────────────
    void OnEffFitCurve();            // fit G3 log-poly to source efficiency points
    void OnEffSaveCurve();           // save current fit/points to eff_curves/<name>.eff
    void OnEffScanCurves();          // scan eff_curves/ and populate listbox
    void OnEffCurveSelected(Int_t id);  // load selected curve into display
    void OnEffPlotCurve();           // draw selected curve on canvas
    void OnEffDeleteCurve();         // delete selected curve file
    void OnEffApplyToHist();         // divide selected histogram by efficiency
    // manual entry
    void OnEffAddStepPoint();        // add (E, eff) row to the manual step list
    void OnEffRemoveStepPoint();     // remove selected row from the manual step list
    void OnEffClearStepPoints();     // clear all manual step points
    void OnEffSaveManual();          // save manually-entered params to eff_curves/

    // ── Recent files ─────────────────────────────────────────────────────────
    void OnOpenRecent();

    // ── Isotopes tab ──────────────────────────────────────────────────────────
    void OnIsoRefresh();
    void OnIsoFilterChanged(Int_t id);
    void OnIsoListSelected(Int_t id);
    void OnIsoApply();
    void OnIsoApplyLabelAll();
    void OnIsoClear();
    void OnIsoClearAll();
    void OnIsoAutoMatchAll();
    void OnIsoApplyAutoMatches();
    void OnIsoDrawSchematic();
    void OnIsoApplyClassToAll();
    void OnIsoSetParent();
    void OnLoadAMETable();
    void OnLoadNubaseTable();
    void OnIsoSetLabelDecay();
    void OnIsoLabelDecayApply();
    void OnIsoLabelDecaySearch();
    void OnIsoLabelDecayClose();
    void OnIsoLabelDecayDlgClosed();
    std::string AutoClassFromParent(const std::string& label);  // beta-minus chain classification
    void OnIsoPeakPreview();
    void OnIsoPopulateIntensity();
    void OnIsoDbSearch();
    void OnIsoDbClear();
    void OnIsoDbLineSelected(Int_t id);
    void OnIsoDbApply();

    // ── Nuclear tab ───────────────────────────────────────────────────────────
    void OnNucDrawLevelScheme();
    void OnNucOpenInteractive();
    void OnNucSetParentFromChain();
    void OnNucAutoTraceChain();
    void OnNucFetchAll();
    void OnNucReloadCache();
    void OnNucAddToChain();
    void OnNucRemoveFromChain();
    void OnNucClearChain();
    void OnNucAddBackground();
    void OnNucLoadNNDCTxt();
    void OnNucConfirmChainToIsoDB();

    // ── Level Scheme tab (Nuclear sub-tab 3 & 4) ──────────────────────────────
    void OnLSSeedFromNNDC();      // populate lsData_ from nuclearDB_ for selected isotope
    void OnLSSave();              // save lsData_ to level_schemes/<isoID>.lsdat
    void OnLSLoad();              // load .lsdat file via file browser
    void OnLSAddLevel();          // add level from entry fields
    void OnLSRemoveLevel();       // remove selected level
    void OnLSLevelSelected(Int_t id);  // fill detail fields when level clicked
    void OnLSAddTransition();     // add transition from combo+entry fields
    void OnLSRemoveTransition();  // remove selected transition
    void OnLSTransSelected(Int_t id);  // fill detail fields when transition clicked
    void OnLSLinkPeak();          // link selected ptRows_ peak to selected transition
    void OnLSDrawEnhanced();      // draw enhanced level scheme using lsData_ + NNDC
    void OnLSCalcBR();            // compute branching ratios for selected level
    void OnLSCalcLogFt();         // compute log ft + B(GT)/B(F) for all fed levels
    void OnLSCalcBalance();       // compute cascade balance for all levels
    void OnLSBrLevelSelected(Int_t id); // update BR display when level combo changes
    // Internal helpers (not slots)
    void BuildLSLevelSchemeTab(TGCompositeFrame* p);
    void BuildLSLogFtTab(TGCompositeFrame* p);
    void RefreshLSLevelList();
    void RefreshLSTransList();
    void RefreshLSLevelCombos();
    void DrawEnhancedLevelScheme(const std::string& isoID);

    // ── Peak Table tab ────────────────────────────────────────────────────────
    void OnPTScanAll();
    void OnPTRestoreMissing();
    void OnPTAddCache();
    void OnPTRemoveCache();
    void OnPTClearCaches();
    void OnPTRebuildTable();
    void OnPTRowSelected(Int_t id);
    void OnPTPlot();
    void OnPTExportCSV();
    void LoadPeakCacheIntoTable(const std::string& cachePath,
                                const std::string& histName);
    // Intensity / efficiency
    void OnPTScanEffCaches();
    void OnPTSaveEffCache();
    void OnPTEffSelected(Int_t id);
    void OnPTCalculateIntensity();
    void OnPTSaveIntensity();
    void OnPTPopulate();
    void OnPTPreviewPeak();
    // DB comparison
    void OnPTDbCompare();

    // ── FWHM tab ──────────────────────────────────────────────────────────────
    TGComboBox*    fwhmSrcCombo_   = nullptr;   // source histogram selector for FWHM
    void OnLoadFWHM();
    void OnRemoveFWHMHist();
    void OnClearFWHMHists();
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
    void UpdateFitEquation();
    void OnShowCompToggled();
    void OnAcceptFit();
    void OnRejectFit();
    void OnClearPeaks();
    void OnRemovePeak();
    void OnSetBgFromCanvas();
    void OnFitBackground();
    void OnClearBackground();
    void OnSetAnch1FromCanvas();
    void OnSetAnch2FromCanvas();
    void OnSetRangeFromCanvas();
    void OnClearFitRange();
    void OnParameterScan();
    void OnShowFitParams();
    void OnFitParamClose();
    void OnSeedBoundsFromModel();
    void OnTransferCache();
    void OnTogglePeakPlaceMode();

    // ── Peak navigation ───────────────────────────────────────────────────────
    void OnPrevPeak();
    void OnNextPeak();
    void OnZoomIn();
    void OnZoomOut();
    void OnNavXRangeGo();
    void OnDeleteCacheEntry();
    void OnToggleMarkRefit();
    void OnAddPeakNoFit();

    void OnApplyPickedLabel();

    // ── Residuals ─────────────────────────────────────────────────────────────
    void OnToggleResiduals();
    void OnSelectResidualFit(Int_t id);

    // ── 2D histogram ─────────────────────────────────────────────────────────
    void OnApplyTh2Labels();
    void OnApplySrcTh2Labels();

    // ── Keyboard shortcuts configuration ─────────────────────────────────────
    void OnConfigureShortcuts();
    void OnShortcutCaptureNext();
    void OnShortcutClear();
    void OnShortcutResetDefaults();
    void OnShortcutSave();
    void OnShortcutDialogClose();

    // ── References sub-tab ───────────────────────────────────────────────────
    void PopulateNucRefTab();
    void OnNucRefOpen();
    void OnNucRefFilterChanged(Int_t id);
    void OnNucSetHistParent();
    void SaveChainCache();
    void LoadChainCache();

    // ── AutoFit tab: histogram search bar ────────────────────────────────────
    void OnHistSearch();

    // ── AutoFit tab: parent nucleus quick-assign ──────────────────────────────
    void OnAutoFitSetHistParent();

    // ── Cache backup ──────────────────────────────────────────────────────────
    void BackupCacheFile(const std::string& srcPath);  // copy to <cacheDir>_backup/

    // ── Utilities ─────────────────────────────────────────────────────────────
    void AppendLog(const std::string& msg);
    void SetStatus(const std::string& msg);
    void OnSavePlot();
    void OnClearHistCache();
    void OnSaveHistCache();
    void OnArchiveHistCache();
    void OnRestoreArchivedCache();
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
    std::string  launchDir_;    // CWD at construction  -  all cache paths anchored here
    OptionalDeps optDeps_;      // probed once at construction
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
    double       lastManualEdm_     = -1.0;  // EDM from last OnManualFit run

    // Histogram classification: "Gamma" | "Decay" | "2D"
    std::map<std::string, std::string> histClass_;

    // Custom projections  -  created by the user with explicit cut ranges
    std::map<std::string, CustomProjDef> customProjDefs_;

    // TH2 support  -  main file
    std::set<std::string>              th2Names_;    // TH2 names in inputFile_
    std::map<std::string, std::string> projParent_;  // "name_px"/"name_py" -> parent TH2 name

    // TH2 support  -  source file
    std::set<std::string>              srcTh2Names_;
    std::map<std::string, std::string> srcProjParent_;
    std::vector<std::string>           srcDirNames_;   // directories found in source file scan

    // ── Shared widgets ────────────────────────────────────────────────────────
    TGStatusBar*         statusBar_ = nullptr;
    TGTextView*          logView_   = nullptr;
    TRootEmbeddedCanvas* canvas_    = nullptr;

    // ── AutoFit tab widgets ───────────────────────────────────────────────────
    std::string    isotopePath_;
    TGLabel*       fileLbl_          = nullptr;
    TGLabel*       isotopeLbl_       = nullptr;
    TGCheckButton* bgSubtractChk_    = nullptr;
    TGCheckButton* showBgLineChk_    = nullptr;
    TGCheckButton* showNegBinsChk_   = nullptr;
    bool           showBgLine_       = false;
    bool           showNegBinsRed_   = true;
    bool           showErrorBars_    = true;
    TGTextButton*  errorBarsBtn_     = nullptr;
    TGNumberEntry* bgIterEntry_      = nullptr;
    TGNumberEntry* tspecSigmaEntry_  = nullptr;
    TGNumberEntry* tspecThreshEntry_ = nullptr;
    TGTextEntry*   histSearchEntry_  = nullptr;  // filter for histList_
    TGListBox*     histList_         = nullptr;
    TGCheckButton* useSeedsChk_      = nullptr;
    TGCheckButton* autoLogLikChk_    = nullptr;
    TGCheckButton* autoImprovChk_    = nullptr;
    TGCheckButton* debugChk_[8]      = {};
    TGComboBox*    histClassCombo_   = nullptr;  // Gamma / Decay / 2D / Background
    TGComboBox*    recentCombo_      = nullptr;  // recent ROOT files
    std::vector<std::string> recentFiles_;       // most-recent-first
    TGTextButton*  peakZoomBtn_       = nullptr;
    bool           peakClickZoom_     = false;

    // Parent nucleus quick-assign (AutoFit tab)
    TGTextEntry*   autoParentEntry_   = nullptr;

    // Custom projection widgets
    TGComboBox*    custProjTh2Combo_  = nullptr;
    TGComboBox*    custProjAxisCombo_ = nullptr;  // 1=ProjX(cut Y), 2=ProjY(cut X)
    TGNumberEntry* custProjLo_        = nullptr;
    TGNumberEntry* custProjHi_        = nullptr;
    TGTextEntry*   custProjName_      = nullptr;
    TGLabel*       custProjBinInfoLabel_ = nullptr;  // live cut-bin / proj-range readout

    // Background subtraction widgets
    TGComboBox*    bgSubSrcCombo_     = nullptr;
    TGComboBox*    bgSubBgCombo_      = nullptr;
    TGNumberEntry* bgSubScaleEntry_   = nullptr;
    TGTextEntry*   bgSubNameEntry_    = nullptr;

    // Rebin widget
    TGNumberEntry* rebinEntry_        = nullptr;   // rebin factor (1 = no rebin)

    // S/N pre-screen threshold for AutoFit (0 = disabled)
    TGNumberEntry* snThreshEntry_     = nullptr;

    // Efficiency correction parameters (ln(eff) = a - b*ln(E) + c*(ln(E))^2 - d/E^2)
    TGNumberEntry* effA_              = nullptr;
    TGNumberEntry* effB_              = nullptr;
    TGNumberEntry* effC_              = nullptr;
    TGNumberEntry* effD_              = nullptr;

    // Efficiency curve cache (Source → Efficiency sub-tab)
    TGListBox*     effCurveList_      = nullptr;  // saved curves
    TGTextEntry*   effCurveName_      = nullptr;  // name for new curve
    TGComboBox*    effCurveType_      = nullptr;  // kLog4/kG3LogPoly/kStep
    TGLabel*       effCurveFitLbl_    = nullptr;  // chi2/ndf of last fit
    std::vector<EfficiencyCache> effCurves_;      // in-memory loaded curves
    int            effCurveSelected_  = -1;       // index into effCurves_
    // last-fitted G3 params (9) and step points for pending save
    std::vector<double>                  effPendingG3_;
    std::vector<std::pair<double,double>> effPendingPoints_;
    // Manual entry widgets
    TGNumberEntry* effManualA_        = nullptr;  // Log4 a
    TGNumberEntry* effManualB_        = nullptr;  // Log4 b
    TGNumberEntry* effManualC_        = nullptr;  // Log4 c
    TGNumberEntry* effManualD_        = nullptr;  // Log4 d
    TGTextEntry*   effManualG3Entry_  = nullptr;  // G3: 9 space-separated values
    TGNumberEntry* effManualStepE_    = nullptr;  // Step: E (keV) to add
    TGNumberEntry* effManualStepV_    = nullptr;  // Step: eff value to add
    TGListBox*     effManualStepList_ = nullptr;  // display of manual step points
    std::vector<std::pair<double,double>> effManualPoints_;  // pending step data

    // Channel->keV calibration parameters (E = a + b*ch + c*ch^2)
    TGNumberEntry* calibA_            = nullptr;
    TGNumberEntry* calibB_            = nullptr;
    TGNumberEntry* calibC_            = nullptr;

    // Background subtraction definitions (persistent via metadata)
    std::map<std::string, BgSubtractDef> bgSubtractDefs_;

    // Stored rebin factors per histogram name (persistent via metadata)
    std::map<std::string, int>           rebinFactors_;
    TGTextEntry*   th2XLabelEntry_ = nullptr;  // axis label for TH2 X axis (main file)
    TGTextEntry*   th2YLabelEntry_ = nullptr;  // axis label for TH2 Y axis (main file)

    // ── AME mass table (NNDC/IAEA AME2020 mass_1.mas) ────────────────────────
    std::map<std::pair<int,int>, double> ameTable_;  // (Z,A) -> mass excess (keV)
    bool           ameLoaded_  = false;
    TGLabel*       ameLbl_     = nullptr;
    bool LoadAMETable(const std::string& path);

    struct NubaseEntry {
        std::string halflifeStr;       // human-readable, e.g. "125 ms"
        double halflifeSec  = -1.0;    // seconds; 0=stable, -1=unknown
        double brBetaMinus  = -1.0;    // % beta-   (-1=unknown)
        double brBetaN      = -1.0;    // % beta-n
        double brBeta2N     = -1.0;    // % beta-2n
    };
    std::map<std::pair<int,int>, NubaseEntry> nubaseTable_;  // (Z,A) -> entry
    bool           nubaseLoaded_ = false;
    TGLabel*       nubaseLbl_    = nullptr;
    bool LoadNubaseTable(const std::string& path);

    // ── Isotopes tab widgets ──────────────────────────────────────────────────
    TGLabel*       isoHistLabel_   = nullptr;  // shows isoHistName_ in the tab
    TGListBox*     isoList_        = nullptr;
    TGComboBox*    isoFilterCombo_ = nullptr;
    TGComboBox*    isoLabelCombo_       = nullptr;  // DB-match combo for label
    TGTextEntry*   isoCustomLabelEntry_ = nullptr;  // free-text custom label (overrides combo)
    TGComboBox*    isoClassCombo_       = nullptr;
    TGTextEntry*   isoCustomEntry_      = nullptr;
    TGNumberEntry* isoMatchThreshEntry_ = nullptr;  // auto-match max distance (keV)
    // parallel arrays: one entry per isoList_ item
    std::vector<std::string> isoListKeys_;       // FitEntry keys matching each row
    std::vector<int>         isoListGaussIdx_;   // Gaussian index in entry (-1 = whole entry)
    std::vector<std::string> isoListAutoMatch_;  // best DB match isotope (may be empty)
    std::vector<double>      isoListDbEnergy_;   // matched DB line energy (0 if none)
    std::string              isoHistName_;  // histogram whose cache is shown
    bool                     schematicDrawn_ = false;  // true after first DrawDecaySchematic call
    // Label -> class mapping (label name -> class string, e.g. "Co-60" -> "Parent")
    std::map<std::string, std::string> labelClassMap_;
    // Parent nucleus info (set by user in Isotopes tab)
    std::string isoParentIsotope_;
    int         isoParentZval_ = 0;
    int         isoParentNval_ = 0;
    // Parent entry widgets
    TGTextEntry*   isoParentNameEntry_ = nullptr;
    TGNumberEntry* isoParentZEntry_    = nullptr;
    TGNumberEntry* isoParentNEntry_    = nullptr;
    // Set Label & Decay dialog (created lazily, nullptr when closed)
    TGTransientFrame* isoLabelDecayDlg_    = nullptr;
    TGListBox*        isoLabelDecayList_   = nullptr;
    TGComboBox*       isoDecayTypeCombo_   = nullptr;
    TGTextEntry*      isoLabelDecaySearch_ = nullptr;
    std::string       isoLabelDecaySelected_;  // isotope name chosen in dialog
    std::string       isoDecayTypeSelected_;   // decay type chosen in dialog
    Int_t             isoLabelDecayPeakSel_ = -1;  // isoList_ selection saved at dialog open
    // Isotope DB browser
    TGListBox*     isoDbList_      = nullptr;
    TGTextEntry*   isoDbSearch_    = nullptr;
    TGComboBox*    isoDbClassCombo_ = nullptr;
    std::vector<int> isoDbIndices_;  // index into db_.db for each visible row
    TGLabel*       isoIntensityInfoLbl_ = nullptr;  // shows matched decay fit A0/T1/2

    // ── Nuclear tab widgets ───────────────────────────────────────────────────
    TGListBox*    nucChainList_     = nullptr;  // isotopes in the decay chain
    TGComboBox*   nucBgCombo_       = nullptr;  // common background selector
    TGNumberEntry* nucAddAEntry_    = nullptr;  // A for manual add
    TGTextEntry*  nucAddSymEntry_   = nullptr;  // symbol for manual add
    TGComboBox*   nucAddTypeCombo_  = nullptr;  // type: auto/decay/background
    TGTextView*   nucGammaRefView_  = nullptr;  // known gammas display
    TGLabel*      nucStatusLbl_     = nullptr;  // status line
    TGTextButton* nucTraceChainBtn_ = nullptr;
    TGNumberEntry* nucAEntry_       = nullptr;  // parent A entry
    TGTextEntry*  nucSymbolEntry_   = nullptr;  // parent symbol entry
    TGComboBox*   nucIsoCombo_      = nullptr;  // isotope combo for Level Scheme
    TGComboBox*   nucIsoLogftCombo_ = nullptr;  // isotope combo for Log ft tab
    std::string   nucCacheDir_;                 // path to nuclear data cache
    std::vector<std::string> nucChainIsotopes_; // isotope IDs in chain (e.g. "44S")
    std::map<std::string, NucIsotope> nuclearDB_; // fetched nuclear data by ID

    // ── References sub-tab widgets ────────────────────────────────────────────
    TGListBox*    nucRefList_        = nullptr;  // gamma reference list
    TGLabel*      nucRefDetailLbl_   = nullptr;  // selected entry detail
    TGComboBox*   nucRefDecayFilter_ = nullptr;  // filter by decay class
    // histogram -> parent nucleus mapping (persisted in cache)
    std::map<std::string, std::string> histParent_; // histName -> parentNucleusID (e.g. "44S")

    // ── Level Scheme data + widgets ───────────────────────────────────────────
    LevelSchemeData lsData_;              // current user-edited level scheme
    // Sub-tab 3: Level Scheme editor
    TGListBox*     lsLevelList_       = nullptr;  // shows all levels
    TGListBox*     lsTransList_       = nullptr;  // shows all transitions
    TGNumberEntry* lsLevelEnergyEntry_= nullptr;  // E for new level (keV)
    TGTextEntry*   lsLevelJpiEntry_   = nullptr;  // Jpi for new level
    TGNumberEntry* lsBetaFeedEntry_   = nullptr;  // beta-feeding % for selected level
    TGTextEntry*   lsBetaTypeEntry_   = nullptr;  // "GT", "Fermi", "Mixed", "1F"
    TGComboBox*    lsTransFromCombo_  = nullptr;  // from-level combo for new transition
    TGComboBox*    lsTransToCombo_    = nullptr;  // to-level combo for new transition
    TGNumberEntry* lsTransEnergyEntry_= nullptr;  // E for new transition (keV)
    TGTextEntry*   lsTransMultEntry_  = nullptr;  // multipolarity ("E2", "M1+E2")
    TGNumberEntry* lsTransIntEntry_   = nullptr;  // intensity for new transition
    // Sub-tab 4: Log ft & BR
    TGComboBox*    lsBrFromCombo_     = nullptr;  // level whose BRs to display
    TGListBox*     lsBrList_          = nullptr;  // branching ratio results
    TGTextView*    lsLogFtView_       = nullptr;  // log ft / B(GT) table
    TGListBox*     lsBalanceList_     = nullptr;  // cascade balance results
    TGNumberEntry* lsParentHlEntry_   = nullptr;  // parent T1/2 (s) override
    TGNumberEntry* lsParentQbEntry_   = nullptr;  // parent Q_beta (keV) override
    TGNumberEntry* lsParentZEntry_    = nullptr;  // parent Z override

    // ── Peak Table tab widgets ────────────────────────────────────────────────
    TGListBox*   ptCacheList_   = nullptr;
    TGTextEntry* ptFilterLabel_ = nullptr;
    TGComboBox*  ptClassFilter_ = nullptr;
    TGComboBox*  ptSortCombo_   = nullptr;
    TGListBox*   ptTableList_   = nullptr;
    TGComboBox*  ptXAxisCombo_  = nullptr;
    TGComboBox*  ptYAxisCombo_  = nullptr;
    TGTextEntry* ptGraphLabel_  = nullptr;
    std::vector<PeakTableRow>   ptRows_;
    std::vector<std::string>    ptLoadedCaches_;
    std::vector<size_t>         ptFilteredRows_;

    // Intensity calculation widgets
    TGComboBox*    ptEffCombo_      = nullptr;  // select efficiency cache
    TGTextEntry*   ptEffNameEntry_  = nullptr;  // name when saving a new cache
    TGNumberEntry* ptPopulationEntry_ = nullptr; // total population (decays)
    TGNumberEntry* ptEnergyEntry_   = nullptr;  // peak energy (auto-filled)
    TGNumberEntry* ptEffValEntry_   = nullptr;  // efficiency at E (auto-filled or manual)
    TGNumberEntry* ptAreaEntry_     = nullptr;  // peak area (auto-filled)
    TGNumberEntry* ptAreaErrEntry_  = nullptr;  // peak area error (auto-filled)
    TGLabel*       ptIntensityLbl_  = nullptr;  // computed intensity result
    TGLabel*       ptDecayInfoLbl_  = nullptr;  // shows A0/T1/2 after Populate
    std::vector<EfficiencyCache> ptEffCaches_;  // loaded efficiency fits
    int            ptSelectedRow_   = -1;       // index into ptRows_ of current selection
    TGCheckButton* ptShowRefitOnly_ = nullptr;  // when on: hide refit-marked peaks (show fitted only)

    // DB comparison widgets
    TGTextEntry*   ptDbIsoEntry_    = nullptr;
    TGNumberEntry* ptDbTolEntry_    = nullptr;
    TGListBox*     ptDbResultList_  = nullptr;

    // ── Fit Results tab widgets ───────────────────────────────────────────────
    TGListBox*           fitResultsList_ = nullptr;
    std::vector<std::string> fittedHists_;     // parallel to fitResultsList_ entries
    std::set<std::string>    fittedSrcHists_;  // subset of fittedHists_ that came from the source tab

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
    // Fit bounds (editable via Fit Parameters popup)
    TGNumberEntry* mEnergyWin_   = nullptr;  // +/- keV window for E
    TGNumberEntry* mAmpLoFrac_   = nullptr;  // A lower = val * frac
    TGNumberEntry* mAmpHiFrac_   = nullptr;  // A upper = val * frac
    TGNumberEntry* mSigLoFrac_   = nullptr;  // sigma lower = model * frac
    TGNumberEntry* mSigHiFrac_   = nullptr;  // sigma upper = model * frac
    TGMainFrame*   fitParamDlg_  = nullptr;
    TGLabel*       mResultLbl_   = nullptr;
    TGTextView*    mFitParamsView_ = nullptr;  // per-parameter values + uncertainties
    TGTextView*    peakStatsView_  = nullptr;  // derived statistics per fitted peak
    TGTextView*    fitEqView_      = nullptr;  // current fit equation display

    // Multi-peak
    std::vector<double>  manualPeaks_;
    TGListBox*           peakListBox_  = nullptr;
    TGCheckButton*       addPeakChk_   = nullptr;

    // Peak label & classification (saved with cache entry in OnAcceptFit)
    TGComboBox*          mPeakLabelCombo_ = nullptr;  // populated with DB matches after fit
    TGComboBox*          mPeakClass_      = nullptr;
    TGTextEntry*         mPeakCustom_     = nullptr;

    // Choose Peaks mode  -  when ON, canvas clicks place peak seeds; default is label-pick
    bool                 peakPlaceMode_   = false;
    TGTextButton*        choosePeakBtn_   = nullptr;

    // Label pick (always-on default)  -  click on canvas to select a cached peak and relabel it
    std::string          labelPickKey_;          // cache key of the selected entry
    int                  labelPickGaussIdx_ = -1; // which Gaussian in that entry (-1 = whole entry)
    TGLabel*             labelPickInfo_   = nullptr;
    TGTextButton*        applyLabelBtn_   = nullptr;

    // Fit options
    TGCheckButton*       mFitLogLikChk_ = nullptr;
    TGCheckButton*       mFitImprovChk_ = nullptr;
    TGCheckButton*       mFitMinosChk_  = nullptr;
    TGCheckButton*       mBgFlatChk_    = nullptr;  // flat (constant) background
    TGCheckButton*       mShowManualFitChk_ = nullptr;  // toggle manual fit overlay
    TGCheckButton*       mShowCacheFitsChk_ = nullptr;  // toggle cached fit overlay
    TGTextButton*        markRefitBtn_       = nullptr;  // mark/unmark current peak for refit
    TGCheckButton*       mShowCompChk_     = nullptr;  // show BG + individual Gaussian components
    TGCheckButton*       showIsoLabelsChk_     = nullptr;
    bool                 showIsoLabels_        = false;
    TGCheckButton*       showPeakFWHMAreaChk_  = nullptr;  // toggle FWHM+area on peak labels
    bool                 showPeakFWHMArea_     = false;
    TGCheckButton*       mBgQuadChk_      = nullptr;  // quadratic background term
    TGCheckButton*       mComptonStepChk_ = nullptr;  // Compton step (Erfc term per peak)
    TGCheckButton*       mExpoTailChk_    = nullptr;  // low-energy exponential tail (HYPERMET)
    TGCheckButton*       mTieWidthsChk_    = nullptr;  // tie sigma to resolution model
    TGCheckButton*       mTieSameSigmaChk_ = nullptr;  // force same sigma for all Gaussians
    TGComboBox*          mResModelCombo_  = nullptr;  // "Auto" | "Custom"
    TGNumberEntry*       mResParA_        = nullptr;  // FWHM^2=a+b*E+c*E^2  -  a
    TGNumberEntry*       mResParB_        = nullptr;  // b
    TGNumberEntry*       mResParC_        = nullptr;  // c

    // BG anchor markers  -  two regions used to seed the linear BG
    TGNumberEntry*       mBgAnch1Lo_      = nullptr;
    TGNumberEntry*       mBgAnch1Hi_      = nullptr;
    TGNumberEntry*       mBgAnch2Lo_      = nullptr;
    TGNumberEntry*       mBgAnch2Hi_      = nullptr;

    // Background region
    TGNumberEntry*       mBgLo_        = nullptr;
    TGNumberEntry*       mBgHi_        = nullptr;
    bool                 bgClickMode_  = false;
    int                  bgClickCount_ = 0;
    bool                 anchClickMode_  = false;
    int                  anchClickCount_ = 0;
    int                  anchTarget_     = 0;  // 0=anchor1, 1=anchor2
    TF1*                 bgTF1_        = nullptr;

    // Fit range (2-click)
    TGNumberEntry*       mFitLo_          = nullptr;
    TGNumberEntry*       mFitHi_          = nullptr;
    bool                 rangeClickMode_  = false;
    int                  rangeClickCount_ = 0;

    // ── Source tab widgets ────────────────────────────────────────────────────
    TGLabel*       srcRootFileLbl_      = nullptr;
    TGComboBox*    recentSrcCombo_      = nullptr;  // recent source ROOT files
    TGTextEntry*   srcHistSearchEntry_  = nullptr;  // filter for srcHistCombo_
    TGComboBox*    srcHistCombo_    = nullptr;
    // Detector-array (TH2) sub-group
    TGComboBox*    srcDetTh2Combo_  = nullptr;  // source TH2 selector
    TGComboBox*    srcDetAxisCombo_ = nullptr;  // which axis holds detector number
    TGNumberEntry* srcDetLoEntry_   = nullptr;  // first detector bin (0 = all)
    TGNumberEntry* srcDetHiEntry_   = nullptr;  // last detector bin (>= lo)
    TGLabel*       srcDetInfoLbl_   = nullptr;
    TH1*           srcDetHist_      = nullptr;  // owned extracted 1D slice
    // Multi-source manager
    TGListBox*     srcTh2List_       = nullptr;  // TH2 names + labels
    TGTextEntry*   srcTh2LabelEntry_ = nullptr;  // label entry for selected TH2
    TGLabel*       srcMultiSrcLbl_   = nullptr;  // status label
    std::map<std::string, std::string> srcTh2Labels_; // th2name -> isotope ID
    std::map<std::string, NucIsotope>  srcNucDB_;     // isoID -> nuclear data
    TGCheckButton* srcBgSubChk_          = nullptr;
    TGTextEntry*   srcTh2XLabelEntry_    = nullptr;
    TGTextEntry*   srcTh2YLabelEntry_    = nullptr;
    TGNumberEntry* srcBgIterEntry_       = nullptr;
    // Extended source-tab fit options (mirror the AutoFit tab options)
    TGNumberEntry* srcTspecSigmaEntry_   = nullptr;  // TSpectrum sigma (bins)
    TGNumberEntry* srcTspecThreshEntry_  = nullptr;  // TSpectrum threshold
    TGCheckButton* srcLogLikChk_         = nullptr;  // log-likelihood fit
    TGCheckButton* srcImprovChk_         = nullptr;  // IMPROVE after MIGRAD
    TGNumberEntry* srcSnThreshEntry_     = nullptr;  // min S/N ratio pre-screen
    TGLabel*       srcFileLbl_        = nullptr;
    TGLabel*       srcInfoLbl_        = nullptr;
    TGComboBox*    srcActivityUnit_   = nullptr;  // 1=Bq, 2=uCi
    TGNumberEntry* srcLiveTime_       = nullptr;
    // Editable source description fields (per-histogram cache)
    TGTextEntry*   srcIsotopeEntry_   = nullptr;  // isotope / source label
    TGTextEntry*   srcCalDateEntry_   = nullptr;  // calibration date yyyy-mm-dd
    TGTextEntry*   srcMeasDateEntry_  = nullptr;  // measurement date yyyy-mm-dd
    // Equations / expected-counts display
    TGTextView*    srcPhysicsView_    = nullptr;
    TGComboBox*    srcIsoCombo_          = nullptr;  // isotope selector for lines browser (ID 957)
    TGTextEntry*   srcLineSearchEntry_   = nullptr;  // energy filter for the lines list
    TGListBox*     srcLineList_          = nullptr;
    TGNumberEntry* srcManualE_           = nullptr;  // kept for legacy; hidden when srcPeakCombo_ present
    TGComboBox*    srcPeakCombo_         = nullptr;  // dropdown of fitted peaks for manual assign

    // ── Energy Calibration tab widgets ────────────────────────────────────────
    TGComboBox*    ecalSrcCombo_      = nullptr;  // source histogram selector
    TGListBox*     ecalHistList_      = nullptr;  // list of added histograms
    TGListBox*     ecalPtList_        = nullptr;  // calibration point list
    TGLabel*       ecalResultLbl_     = nullptr;
    TGNumberEntry* ecalA_             = nullptr;
    TGNumberEntry* ecalB_             = nullptr;
    TGNumberEntry* ecalC_             = nullptr;
    TGTextEntry*   ecalNameEntry_     = nullptr;
    TGComboBox*    ecalApplyCombo_    = nullptr;  // saved calibration selector (src tab)
    TGComboBox*    autoFitEcalCombo_  = nullptr;  // saved calibration selector (autofit tab)
    TGComboBox*    effEcalCombo_      = nullptr;  // saved calibration selector (efficiency tab)

    // Energy calibration data
    std::vector<double>      ecalAllX_;           // fitted peak positions
    std::vector<double>      ecalAllY_;           // reference energies from NNDC
    std::vector<std::string> ecalHistSources_;    // per-point source histogram name
    std::vector<std::string> ecalLoadedHists_;    // histograms contributing
    double ecalFitA_ = 0.0, ecalFitB_ = 1.0, ecalFitC_ = 0.0;

    // Source ROOT file (separate from the main inputFile_)
    TFile*                   srcRootFile_  = nullptr;
    std::string              srcRootPath_;
    std::vector<std::string> recentSrcFiles_;  // most-recent-first
    std::vector<std::string> srcHistNames_;

    // Per-histogram source description metadata
    struct SourceHistMeta {
        std::string isotope;
        std::string calDate;
        std::string measDate;
        double      livetime        = 1.0;
        double      activity        = 0.0;   // Bq at calDate
        std::string th2Parent;               // TH2 this 1D was projected from
        std::string th2Label;                // source/isotope label of parent TH2
        std::string sourceRootFile;          // source ROOT file the projection was derived from
        std::string externalFile;            // if set, load histogram from this file
        std::string pathInFile;              // path within externalFile (e.g. "Co60/Co60_det01")
    };
    std::map<std::string, SourceHistMeta> srcHistMeta_;  // histname -> metadata

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

    // Per-histogram auto-match results (populated by OnSrcAutoMatchAll)
    std::map<std::string, std::vector<SrcMatchedLine>> srcAutoMatches_;
    
    // ── Source Fit Database (separate cache for source spectrum fits) ─────────
    SourceFitDatabase   srcFitDB_;
    std::string         srcFitCachePath_ = "fit_caches/source_fits.txt";
    
    // ── 2D Histogram Preview ─────────────────────────────────────────────────
    TH2*                srcTh2Preview_  = nullptr;  // current TH2 being previewed
    bool                srcShowingTh2_  = false;    // whether 2D view is active
    TRootEmbeddedCanvas* srcPreviewCanvas_ = nullptr;
    
    // ── Calibration Infrastructure ──────────────────────────────────────────
    CalibrationBuilder  calibBuilder_;          // aggregates calibration points
    LinearCalibrationModel currentCalibration_;  // fitted linear model
    EnergyCalibration   energyCalibManager_;    // manages multiple calibrations
    std::string         activeCalibName_;       // name of active calibration
    bool                calibrationFitted_ = false;
    
    // ── Calibration UI Widgets ─────────────────────────────────────────────
    TGListBox*     calibPointList_     = nullptr;  // display calibration points
    TGNumberEntry* calibPointEnergyEntry_ = nullptr;
    TGNumberEntry* calibPointChannelEntry_ = nullptr;
    TGTextEntry*   calibPointSourceEntry_ = nullptr;
    TGLabel*       calibResultsLbl_    = nullptr;  // displays fitted A, B, chi^2
    TGCheckButton* calibUsePointChk_   = nullptr;  // enable/disable point
    TGComboBox*    calibNameCombo_     = nullptr;  // select which calibration to use

    // ── FWHM tab widgets ─────────────────────────────────────────────────────
    TGComboBox*    fwhmCombo_         = nullptr;
    TGListBox*     fwhmHistList_      = nullptr;   // loaded histograms for combined plot
    TGNumberEntry* mFwhmA_            = nullptr;
    TGNumberEntry* mFwhmB_            = nullptr;
    TGNumberEntry* mFwhmC_            = nullptr;
    TGNumberEntry* mFwhmAlo_          = nullptr;
    TGNumberEntry* mFwhmAhi_          = nullptr;
    TGNumberEntry* mFwhmBlo_          = nullptr;
    TGNumberEntry* mFwhmBhi_          = nullptr;
    TGNumberEntry* mFwhmClo_          = nullptr;
    TGNumberEntry* mFwhmChi_          = nullptr;
    TGLabel*       fwhmResultLbl_     = nullptr;
    TGCheckButton* fwhmShowSigmaChk_  = nullptr;
    TGCheckButton* fwhmStatLineChk_   = nullptr;
    TGCheckButton* fwhmShowResChk_    = nullptr;
    TGCheckButton* fwhmRemoveModeChk_ = nullptr;
    TF1*           fwhmTF1_           = nullptr;   // owned; clones drawn to canvas
    std::string    fwhmHistName_;                  // last loaded / "Combined"
    std::vector<std::string> fwhmLoadedHists_;    // names of loaded histograms in order
    std::vector<std::string> fwhmHistSources_;    // per-point: which histogram it came from
    std::vector<double> fwhmAllX_;       // all FWHM data points (energy)
    std::vector<double> fwhmAllY_;       // all FWHM data points (FWHM^2 value, keV^2)
    std::vector<bool>   fwhmExcluded_;   // parallel exclusion flags
    std::vector<bool>   fwhmTied_;       // true = sigma was fixed to res model
    double         fwhmChi2Ndf_   = -1.0;
    double         fwhmPValue_    = -1.0;
    double         fwhmResidRMS_  = -1.0;
    int            fwhmNdf_       = 0;

    // ── Decay tab widgets ─────────────────────────────────────────────────────
    TGTab*         decaySubTabs_        = nullptr;  // "Cuts" / "Fitter" sub-tabs
    TGComboBox*    decayTh2Combo_       = nullptr;  // TH2 selection
    TGComboBox*    decayGammaAxisCombo_ = nullptr;  // 1=X is gamma, 2=Y is gamma
    TGNumberEntry* decaySigLoEntry_     = nullptr;  // Lo-side sigma window (below centroid)
    TGNumberEntry* decaySigRangeEntry_  = nullptr;  // Hi-side sigma window (above centroid)
    TGLabel*       decayBinInfoLabel_   = nullptr;  // live readout of cut edges in keV + bins
    TGListBox*     decayPeakList_       = nullptr;  // fitted peaks from gamma projection
    TGComboBox*    decayModelCombo_     = nullptr;
    TGComboBox*    decayBGTypeCombo_    = nullptr;  // 1=Flat, 2=Flat+Exp, 3=Exp only
    // Half-life rows  -  labels are updated dynamically per model
    TGLabel*       decayThalfPLbl_      = nullptr;
    TGNumberEntry* decayThalfP_         = nullptr;
    TGCheckButton* decayFixP_           = nullptr;
    TGLabel*       decayThalfDLbl_      = nullptr;
    TGNumberEntry* decayThalfD_         = nullptr;
    TGCheckButton* decayFixD_           = nullptr;
    TGLabel*       decayThalfGLbl_      = nullptr;
    TGNumberEntry* decayThalfG_         = nullptr;
    TGCheckButton* decayFixG_           = nullptr;
    TGLabel*       decayThalfBGLbl_     = nullptr;
    TGNumberEntry* decayThalfBGExp_     = nullptr;  // T1/2 of exponential BG component (ms)
    TGCheckButton* decayFixBGExp_       = nullptr;
    // Per-parameter bounds (lo/hi) for each T1/2 — used to call SetParLimits
    TGNumberEntry* decayThalfPLo_      = nullptr;
    TGNumberEntry* decayThalfPHi_      = nullptr;
    TGNumberEntry* decayThalfDLo_      = nullptr;
    TGNumberEntry* decayThalfDHi_      = nullptr;
    TGNumberEntry* decayThalfGLo_      = nullptr;
    TGNumberEntry* decayThalfGHi_      = nullptr;
    TGNumberEntry* decayThalfBGExpLo_  = nullptr;
    TGNumberEntry* decayThalfBGExpHi_  = nullptr;
    TGTextView*    decayResultView_     = nullptr;
    TGTextEntry*   decayLabelEntry_     = nullptr;
    TGComboBox*    decayLabelClassCombo_= nullptr;
    TGNumberEntry* decaySliceWidthEntry_= nullptr;
    TGNumberEntry* decayFitLo_          = nullptr;
    TGNumberEntry* decayFitHi_          = nullptr;
    TGCheckButton* decayFitFullRange_   = nullptr;
    TGNumberEntry* decayRebinEntry_     = nullptr;
    TGLabel*       decayEquationLabel_  = nullptr;
    TGComboBox*    decayFitMethod_      = nullptr;
    TGCheckButton* decayShowResid_      = nullptr;
    TGCheckButton* decayErrBars_        = nullptr;
    TGNumberEntry* decayLegX_           = nullptr;
    TGNumberEntry* decayLegY_           = nullptr;
    TGListBox*     decayCacheBrowserList_ = nullptr;
    std::vector<double> decayCacheBrowserEs_; // parallel to list entries
    // Extra gamma-cut additions (summed into decay histogram)
    TGNumberEntry* decayExtraCutLo_     = nullptr;
    TGNumberEntry* decayExtraCutHi_     = nullptr;
    TGListBox*     decayExtraCutList_   = nullptr;
    std::vector<std::pair<double,double>> extraDecayCuts_;

    // Decay fit result  -  stored per peak energy for the active TH2
    struct DecayFitResult {
        double peakE  = 0.0;
        int    model  = 0;
        int    bgType = 1;  // 1=Flat, 2=Flat+Exp, 3=Exp only
        int    rebin  = 1;  // time-axis rebin factor used for the fit
        std::vector<double> params;
        std::vector<double> errors;
        std::vector<double> paramLo;   // lower bounds passed to SetParLimits
        std::vector<double> paramHi;   // upper bounds passed to SetParLimits
        double chi2ndf = -1.0;
        int    status  = -1;
        int    fitMethod = 2;     // 1=Chi2 2=Likelihood 3=Chi2+MINOS 4=Likelihood+MINOS
        double fitLo   = 0.0;    // fit range low edge (ms)
        double fitHi   = 0.0;    // fit range high edge (ms)
        bool   fullRange = true; // true = fit full histogram range
        // provenance / annotation
        std::string histName;        // which TH2 this came from
        double      eMin   = 0.0;   // exact gamma energy cut low edge
        double      eMax   = 0.0;   // exact gamma energy cut high edge
        double      Nsig   = 1.0;   // Hi-side sigma window
        double      NsigLo = -1.0;  // Lo-side sigma window (-1 = same as Nsig)
        std::string label;
        std::string classification;
        std::vector<std::pair<double,double>> extraCuts; // additional gamma-energy cuts summed in
    };

    // Decay state
    std::string                        decayTh2Name_;
    std::string                        decayGammaProjName_;
    std::vector<double>                decayPeakEs_;
    std::vector<double>                decayPeakSigs_;
    std::vector<std::string>           decayPeakKeys_;
    std::map<double, DecayFitResult>   decayFitStore_;  // keyed by peak energy

    // Decay cache picker (which gamma histogram's cache to use for tie widths)
    TGComboBox*              decayCacheCombo_      = nullptr;
    std::vector<std::string> decayPeakCacheNames_; // cache name for each peak in decayPeakList_

    // ── Total Decay sub-tab widgets ───────────────────────────────────────────
    TGComboBox*    decayTdModelCombo_   = nullptr;
    TGComboBox*    decayTdBGCombo_      = nullptr;
    TGLabel*       decayTdEquationLbl_  = nullptr;
    TGLabel*       decayTdPLbl_         = nullptr;
    TGNumberEntry* decayTdThalfP_       = nullptr;
    TGCheckButton* decayTdFixP_         = nullptr;
    TGLabel*       decayTdDLbl_         = nullptr;
    TGNumberEntry* decayTdThalfD_       = nullptr;
    TGCheckButton* decayTdFixD_         = nullptr;
    TGLabel*       decayTdGLbl_         = nullptr;
    TGNumberEntry* decayTdThalfG_       = nullptr;
    TGCheckButton* decayTdFixG_         = nullptr;
    TGLabel*       decayTdBGLbl_        = nullptr;
    TGNumberEntry* decayTdThalfBGExp_   = nullptr;
    TGCheckButton* decayTdFixBGExp_     = nullptr;
    TGNumberEntry* decayTdFitLo_        = nullptr;
    TGNumberEntry* decayTdFitHi_        = nullptr;
    TGCheckButton* decayTdFullRange_    = nullptr;
    TGCheckButton* decayTdSumAll_       = nullptr;  // sum all peaks (vs selected only)
    TGTextView*    decayTdResultView_   = nullptr;
    DecayFitResult decayTdFitResult_;
    bool           decayTdFitValid_     = false;
    TGNumberEntry* decayTdRebinEntry_   = nullptr;
    TGCheckButton* decayTdShowResid_    = nullptr;
    TGCheckButton* decayTdErrBars_      = nullptr;
    TGComboBox*    decayTdFitMethod_    = nullptr;
    TGNumberEntry* decayTdLegX_         = nullptr;
    TGNumberEntry* decayTdLegY_         = nullptr;

    // ── RK4 Full Chain sub-tab ────────────────────────────────────────────────
    TGCheckButton* rkSumAll_      = nullptr;
    TGLabel*       rkPLbl_        = nullptr;
    TGNumberEntry* rkThalfP_      = nullptr;
    TGCheckButton* rkFixP_        = nullptr;
    TGLabel*       rkDLbl_        = nullptr;
    TGNumberEntry* rkThalfD_      = nullptr;
    TGCheckButton* rkFixD_        = nullptr;
    TGLabel*       rkBnLbl_       = nullptr;
    TGNumberEntry* rkThalfBn_     = nullptr;
    TGCheckButton* rkFixBn_       = nullptr;
    TGLabel*       rkGdLbl_       = nullptr;
    TGNumberEntry* rkThalfGd_     = nullptr;
    TGCheckButton* rkFixGd_       = nullptr;
    TGNumberEntry* rkEffP_        = nullptr;
    TGCheckButton* rkFixEffP_     = nullptr;
    TGNumberEntry* rkEffD_        = nullptr;
    TGCheckButton* rkFixEffD_     = nullptr;
    TGNumberEntry* rkEffBn_       = nullptr;
    TGCheckButton* rkFixEffBn_    = nullptr;
    TGNumberEntry* rkEffGd_       = nullptr;
    TGCheckButton* rkFixEffGd_    = nullptr;
    TGNumberEntry* rkPn_          = nullptr;
    TGCheckButton* rkFixPn_       = nullptr;
    TGNumberEntry* rkN0_          = nullptr;
    TGCheckButton* rkFixN0_       = nullptr;
    TGNumberEntry* rkBg_          = nullptr;
    TGCheckButton* rkFixBg_       = nullptr;
    TGNumberEntry* rkNsteps_      = nullptr;
    TGListBox*     rkPeakList_    = nullptr;
    TGComboBox*    rkChainModel_  = nullptr;
    TGNumberEntry* rkLegX_        = nullptr;
    TGNumberEntry* rkLegY_        = nullptr;
    TGCheckButton* rkUse1DHist_   = nullptr;  // use direct 1D histogram instead of TH2 projection
    TGComboBox*    rkHistCombo_   = nullptr;  // available 1D histograms
    std::vector<std::string> rkHistSpecs_;    // "filePath::histName" parallel to rkHistCombo_ entries
    TGNumberEntry* rkRebinEntry_  = nullptr;
    TGCheckButton* rkShowResid_   = nullptr;
    TGCheckButton* rkErrBars_     = nullptr;
    TGComboBox*    rkFitMethod_   = nullptr;
    TGNumberEntry* rkFitLo_       = nullptr;
    TGNumberEntry* rkFitHi_       = nullptr;
    TGCheckButton* rkFitFull_     = nullptr;
    TGTextView*    rkResultView_  = nullptr;

    // ── Chi² scan popup buttons ───────────────────────────────────────────────
    TGTextButton*  decayChiScanBtn_   = nullptr;
    TGTextButton*  decayTdChiScanBtn_ = nullptr;
    TGTextButton*  rkChiScanBtn_      = nullptr;

    // ── Stored last-fit state for chi² scan ───────────────────────────────────
    TH1*   lastDecayCutsHist_  = nullptr;
    TF1*   lastDecayCutsTF1_   = nullptr;
    double lastDecayCutsXlo_   = 0.0;
    double lastDecayCutsXhi_   = 0.0;
    TH1*   lastDecayTdHist_    = nullptr;
    TF1*   lastDecayTdTF1_     = nullptr;
    double lastDecayTdXlo_     = 0.0;
    double lastDecayTdXhi_     = 0.0;
    TH1*   lastRK4Hist_        = nullptr;
    TF1*   lastRK4TF1_         = nullptr;
    double lastRK4Xlo_         = 0.0;
    double lastRK4Xhi_         = 0.0;

    // ── Peak navigation ───────────────────────────────────────────────────────
    std::vector<std::string> peakNavKeys_;
    int              peakNavIdx_  = 0;
    TGLabel*         peakNavLbl_  = nullptr;
    TGNumberEntry*   navXMinEntry_ = nullptr;
    TGNumberEntry*   navXMaxEntry_ = nullptr;

    // ── Residuals ─────────────────────────────────────────────────────────────
    TGCheckButton*   residualChk_   = nullptr;
    TGComboBox*      residualCombo_ = nullptr;
    bool             residualsOn_   = false;

    // Zoom window (keV); 0,0 = full range
    double viewXmin_ = 0.0;
    double viewXmax_ = 0.0;

    // ── Keyboard shortcuts state ──────────────────────────────────────────────
    std::vector<KbShortcut>              shortcuts_;
    int                                  capturingShortcutIdx_ = -1;
    TGTransientFrame*                    shortcutDlg_          = nullptr;
    std::vector<std::pair<Int_t,Int_t>>  boundKeys_;   // (keycode, modifier) registered

    // ── Build helpers ─────────────────────────────────────────────────────────
    void BuildAutoFitTab    (TGCompositeFrame* parent);
    void BuildSourceTab     (TGCompositeFrame* parent);
    void BuildManualFitTab  (TGCompositeFrame* parent);
    void BuildFWHMTab       (TGCompositeFrame* parent);
    void BuildDecayTab      (TGCompositeFrame* parent);
    void BuildFitResultsTab (TGCompositeFrame* parent);
    void BuildPeakTableTab  (TGCompositeFrame* parent);
    void BuildNuclearTab    (TGCompositeFrame* parent);
    void PopulateHistWidgets();
    std::string NucCacheDirPath() const;
    std::string GenerateLevelSchemePlotlyHTML(const std::string& isoID) const;
    void PopulateNucGammaRef();
    void RefreshIsoComboHelper(TGComboBox* combo,
                               const std::vector<std::string>& chain,
                               const std::map<std::string, NucIsotope>& db);
    void DrawLevelScheme(const std::string& isoID);
    void SyncDebugToggles();
    void DrawOnCanvas(TH1* h, TF1* fit = nullptr);
    void DrawWithResiduals(TH1* h, TF1* fit, double xlo, double xhi);
    void DrawDecayResiduals(TH1* h, TF1* fit, double xlo, double xhi);
    void DrawFitComponents(TCanvas* c, TF1* f);
    void RedrawCurrent();
    void RunFitOnHistogram(const std::string& hname,
                           TFile* overrideFile = nullptr,
                           const std::vector<double>& forcedSeeds = {});
    std::string CacheDirFor() const;
    std::string ArchiveDirFor() const;
    void        EnsureCacheDir() const;
    std::string CacheFileFor(const std::string& hname) const;
    std::string DecayCacheFileFor() const;
    void        SaveDecayFitCache();
    void        LoadDecayFitCache();
    std::string TotalDecayCacheFileFor() const;
    void        SaveTotalDecayFitCache();
    void        LoadTotalDecayFitCache();

    // Fit Results helpers
    void ShowFitResult(const std::string& hname);
    TH1* MakeBgSubHist(TH1* raw, bool doSubtract = false, int iterations = 14);
    TH1* GetTSpectrumBg(TH1* raw, int iterations = 14);
    void SyncResParFields();   // push res_.a/b/c to display fields (Auto mode only)
    ResolutionModel GetTieResModel() const;  // effective model for Tie widths
    void OverlayFitPeaks(const std::string& hname, TCanvas* c);

    // Peak stats snapshot for before/after Run Fit comparison
    std::vector<std::string> peakStatsCurrent_;
    std::vector<std::string> peakStatsOld_;

    // FWHM helpers
    void RedrawFWHM();
    void DrawFWHMToCanvas(TCanvas* c, bool showSigma, bool showStatLine,
                          bool showResolution = false, bool showErrBars = false);

    // Source tab helpers
    double      ComputeDecayedActivity() const;
    void        PopulateSourceList();
    void        ExtractPeaksFromCache(const std::string& hname);
    void        UpdateSourceInfoLabel();
    void        PopulateSrcHistCombo();
    void        PopulateSrcTh2List();
    std::string SourceAnalysisFileFor(const std::string& hname) const;
    void        SaveSourceAnalysis();
    // Source histogram metadata helpers
    std::string SrcHistMetaPath() const;
    void        SaveSrcHistMeta();
    void        LoadSrcHistMeta();
    void        ApplySrcMetaToUI(const std::string& histName);
    void        CollectSrcMetaFromUI(const std::string& histName);
    // TH2 label persistence (Multi-Source Manager)
    std::string SrcTh2LabelsPath() const;
    void        SaveSrcTh2Labels();
    void        LoadSrcTh2Labels();
    // Source recent files
    void        LoadRecentSrcFiles();
    void        SaveRecentSrcFiles() const;
    void        AddToRecentSrcFiles(const std::string& path);
    void        RefreshRecentSrcCombo();
    // Load a source histogram  -  checks externalFile in metadata first
    TH1*        GetSrcHistogram(const std::string& histName, bool& owned) const;

    // Peak-nav helpers
    void NavigateToPeak(int idx);
    void PopulateNavAndResidual();
    TF1* BuildFromCacheKey(const std::string& key);
    TF1* RebuildFromEntry(const FitEntry& entry, double xlo, double xhi);
    void DrawPeakLabels(TF1* f);
    void RedrawView();
    void SaveFitResultToFile(const std::string& hname, TFile* fout);
    void UpdatePeakStats(TF1* f, TH1* h, double xlo, double xhi,
                         TF1* cachedF = nullptr, TH1* cachedH = nullptr,
                         double cxlo = 0.0, double cxhi = 0.0);
    void PopulateIsoList(const std::string& filterLabel = "");
    void PopulateIsoDbList(const std::string& filter = "");
    void RefreshIsoDisplay();
    void DrawDecaySchematic(TCanvas* c);
    void LoadLabelClassMap(FitDatabase& fitdb);
    void SaveLabelClassMap(FitDatabase& fitdb);
    void PopulateDecayTh2Combo();
    void PopulateCustProjTh2Combo();
    void PopulateBgSubCombos();
    void AddToRecentFiles(const std::string& path);
    void LoadRecentFiles();
    void SaveRecentFiles() const;
    void RefreshRecentCombo();

    // gnuScope export helper  -  returns false on I/O error
    bool ExportGnuScopeFile(TH1* h, const std::string& outPath) const;
    bool ExportGnuScopeFit(TH1* h, const std::string& outPath,
                           const std::string& histName) const;

    // Histogram classification helpers
    std::string ClassOf(const std::string& name) const;
    std::string MetadataFileFor() const;
    void        SaveMetadata() const;
    void        LoadMetadata();

    // Keyboard shortcut helpers
    void        InitShortcuts();
    void        ApplyShortcuts();
    void        LoadShortcuts();
    void        SaveShortcuts() const;
    void        DispatchShortcut(const std::string& id);
    std::string ShortcutsFilePath() const;
    static std::string KeybindingText(UInt_t keysym, UInt_t modifier);

    void SafeDeleteHist(TH1*& h);  // remove from canvas primitives then delete

    // Returns histogram for hname from inputFile_; creates projection if needed.
    // Sets owned=true when caller must delete the returned pointer.
    TH1* LoadHistFromFile(const std::string& hname, bool& owned) const;
    TH1* LoadProjection(TFile* f, const std::string& projName,
                        const std::map<std::string, std::string>& projParent,
                        const TGTextEntry* xLblEntry,
                        const TGTextEntry* yLblEntry) const;
};

#endif
