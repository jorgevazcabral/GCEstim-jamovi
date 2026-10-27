# GCEstimjmv v1.1.0

This release introduces a standalone Ridge Trace analysis and graphical improvements across GCEstimjmv for jamovi.

## New: Ridge Trace

A dedicated analysis based on `GCEstim::ridgetrace()`, with:

- Automatic geometric and custom ridge penalty grids.
- Cross-validation with configurable folds and random seed.
- Minimum-error and module-derived 1-SE lambda selection.
- Manual selection from the evaluated lambda grid.
- Coefficient estimates, coefficient ranges and model summaries.
- Seven plots for coefficient trajectories, regularization diagnostics and cross-validation stability.
- Cached numerical-result reuse when changing graphical options or lambda selection.

## Improvements to existing analyses

WGCE, TARWGCE and TASWGCE retain their interfaces and gain consistent coefficient colours, improved coefficient/factor labels, corrected selected-point identities, improved legends, corrected normalized entropy trajectories and saved-state compatibility improvements where applicable.

## Bug fixes

- Fixed incorrect coefficient-colour associations in selected support points.
- Fixed omission of valid entropy observations following failed support configurations.
- Fixed graphical theme overrides affecting coefficient colours.
- Corrected standalone Ridge Trace selected-lambda reference-line presentation and preserved coefficient-specific triangular extrema.
- Fixed coefficient-path and ridge legend clipping using adaptive wrapping.

## Compatibility and known limitations

Existing GCEstim estimation procedures are preserved. The underlying GCEstim R package and pinned dependency are unchanged.

The deferred upstream GCEstim support-selection 1-SE concern remains unresolved; returned TARWGCE/TASWGCE selections are preserved. Standardized ridge cross-validation currently uses full-data preprocessing before folds. Ridge coefficient ranges and stability displays are descriptive, not confidence intervals. Unsupported engine configurations produce informative errors; MAPE and MASE are not exposed.

Older embedded images may require toggling OLS comparison on and off, then saving. Legacy entropy states without sufficient verified recovery data display an explicit notice.

## Installation

Download the attached `GCEstimjmv_1.1.0.jmo` and install it through **Modules → Side-load** in jamovi. This is a Windows x64 build validated with jamovi 2.7.41; it is not a macOS/Linux build. Compatible prebuilt JMO users do not need a separate R installation.

This GitHub release does not imply that v1.1.0 is available in the jamovi Library. No Library submission is part of this release.

## Validation

The complete automated suite passed **1,999 assertions**, with zero failures, test warnings or skips. The Windows build uses jamovi 2.7.41.

Actual jamovi checks confirmed module version 1.1.0, all four analysis entries, complete wrapped coefficient-path legends, entropy plots, all seven Ridge Trace plots, PNG export and saved-file reopening. Fitted caches remained identical after graphical refreshes; TARWGCE and TASWGCE retained all 20 and 25 valid entropy trajectory observations, respectively.
