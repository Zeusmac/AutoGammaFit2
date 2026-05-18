#ifndef ANALYSISTAB_H
#define ANALYSISTAB_H

#include <string>
#include "UISettings.h"

// ─────────────────────────────────────────────────────────────────────────────
// AnalysisTab — lightweight interface that all tab implementations satisfy.
//
// Each tab owns a ROOT widget tree built in BuildUI(). The GammaFitGUI
// main class calls OnHistogramChanged() when the active histogram switches,
// SaveState()/RestoreState() on shutdown/startup, and Reset() when the
// user clears the session.
//
// Not a TObject — tabs communicate back to the host via direct slot calls
// on GammaFitGUI (passed as a pointer) rather than through ROOT signals, to
// keep the interface thin and avoid rootcling overhead for the base class.
// ─────────────────────────────────────────────────────────────────────────────
class AnalysisTab {
public:
    virtual ~AnalysisTab() = default;

    // Display name shown in the tab strip.
    virtual std::string Name() const = 0;

    // Called after the user selects a different histogram.
    virtual void OnHistogramChanged(const std::string& histName) {}

    // Persist tab-specific UI state into settings (called on shutdown).
    virtual void SaveState(UISettings& /*s*/) {}

    // Restore tab-specific UI state from settings (called after BuildUI).
    virtual void RestoreState(const UISettings& /*s*/) {}

    // Clear any cached results / displayed data (user-triggered reset).
    virtual void Reset() {}
};

#endif // ANALYSISTAB_H
