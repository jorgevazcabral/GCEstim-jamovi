# Standalone Ridge Trace

Current release: **1.1.0**. The Windows JMO and source are distributed through the GitHub v1.1.0 release. Historical local-validation records below describe earlier 1.0.5 development builds. No jamovi Library submission is part of this release.

The registered analyses are **WGCE**, **TARWGCE**, **TASWGCE**, and **RidgeTrace** (menu title **Ridge Trace**). LMGCE is not a registered analysis and is not restored by this version. References to GCEstim::lmgce() or lmgce objects describe the entropy analyses’ engine and object classes, not an LMGCE module analysis.

## Estimation and controls

The analysis appears under GCEstim / Regression / Ridge Trace. Estimation is exclusively GCEstim::ridgetrace(), using the existing pinned GCEstim dependency. GCEstim itself and the existing analyses are unchanged.

Select one numeric dependent variable, numeric covariates and/or factors. Complete cases are used, unused factor levels are dropped, and coefficient names follow the engine's model matrix. An intercept and standardization are enabled initially. Intercept penalization applies only without standardization and when an intercept is present.

The default grid contains 100 geometric penalties from 0.001 to 1000. A custom comma-separated vector must contain at least two distinct finite positive values. Custom values are sorted ascending. Matrix columns are aligned by position with the numeric lambda vector; rounded column names are discarded.

CV defaults to five folds and seed 230676. Available error measures are RMSE, MSE, MAE and sMAPE. Returned error values are never replaced by a separate error calculation.

## Lambda selection

For K folds, the CV mean at lambda l is mean(E[,l]); its descriptive standard error is sd(E[,l])/sqrt(K). All folds must be finite for that lambda's CV mean and SE to be available. Automatic selection excludes unavailable errors; it does not discard individual invalid folds and average the remainder.

Minimum error selects the smallest evaluated lambda attaining the minimum mean CV error (ascending grid, first tie). Without CV, it selects minimum returned in-sample error and is explicitly labelled descriptive, without claims of optimal predictive performance.

The 1-SE rule selects the largest evaluated lambda with mean CV error no greater than the minimum mean CV error plus the SE at that minimum. The 1-SE rule is calculated by GCEstimjmv using the fold-specific CV errors returned by GCEstim::ridgetrace(). It requires CV. Fold SE bars are not confidence intervals for generalization error.

Manual selection uses an integer grid position (1 = smallest lambda). The optional Lambda Results table shows these positions. This avoids ambiguity from rounded displayed values and never interpolates coefficients or refits a model. The summary and coefficient table title identify the selected numeric lambda.

If no finite automatic selection exists, paths and descriptive grid ranges remain available. Selected estimates are unavailable until a manual position is chosen.

## Tables

- Model Summary: formula, complete cases, number of predictor columns, intercept and standardization details, grid, CV settings, error measure, selection and GCEstim version.
- Selected Coefficients: original-unit coefficients from the selected engine path column. No SEs, confidence intervals or p-values are invented.
- Coefficient Ranges: returned minimum, maximum and maximum absolute coefficient plus selected estimate. Ranges describe the evaluated grid, not confidence intervals.
- Lambda Results (optional): grid position, lambda, returned in-sample error, CV mean/SE and minimum, 1-SE and selected indicators. Ascending-order rows are limited by a display control; the limit never changes estimation or selection.

## Plots

Default: coefficient paths and CV mean error. The in-sample plot defaults to automatic display when CV is disabled; its Show/Hide setting can override automatic visibility.

There are seven plots: Coefficient Paths, Coefficient Ranges, Local Coefficient Sensitivity, Cross-Validation Error, In-Sample Error, Fold Error Paths, and Fold-Training Coefficient Stability. Ranges, sensitivity, fold errors and stability are optional. The in-sample plot can also be shown with CV enabled.

Coefficient paths support log10 or linear lambda axes, a zero reference, intercept display (enabled by default), extrema markers and facets. Displayed coefficient names can be restricted using a comma-separated list; quote names containing commas with double quotes. Names must match the Coefficients table exactly. This filter affects coefficient plots only.

Coefficient extrema follow TARWGCE: maxima use upward triangles (shape 24), minima downward triangles (shape 25), size 1 and coefficient colours. CV references use black dotted lines for Minimum error and black dashed lines for 1-SE rule. Coefficient Ranges marks the selected estimates as points and labels the selected lambda; the ranges are not confidence intervals.

The heatmap and standalone selected-lambda coefficient plot were removed from local v1.0.5 after interface review. Their controls, result objects and renderers are removed; the Selected Coefficients table remains available.

Local sensitivity is (b[j,l+1]-b[j,l])/(log10(lambda[l+1])-log10(lambda[l])), plotted at log-interval midpoints. It is descriptive and grid-dependent, not significance or an automatic optimum criterion.

Fold error paths optionally overlay their mean. Stability facets show fold-training paths and the full-data coefficient path. They are not bootstrap distributions or coefficient confidence intervals. Fold-based coefficient SD was optional and is not added; the returned trajectories are shown directly.

Plots use ggplot2 render callbacks with jamovi themes and colours. Plot data come from numerical engine fields; plot.ridgetrace() and its internal ggplot layers are not used. Tables and coefficient paths survive failures of optional plots.

## Preserved GCEstim limitations

When standardization and CV are enabled, a results note explains that the current engine preprocesses the complete dataset before making folds. No fold-specific preprocessing is substituted.

The backend rejects standardized models without an intercept and standardized CV models with an intercept and one expanded predictor column. Messages explicitly identify an implementation limitation rather than a statistical restriction. The second guard also covers a binary factor producing the same matrix shape.

The existing CV accuracy-call argument reversal affects MAPE and MASE, which are not exposed. The four exposed measures are symmetric under argument reversal. Undefined/nonfinite returned errors generate warnings and unavailable entries rather than replacement statistics.

Unexpected engine errors and incompatible returned dimensions produce coded, informative analysis errors. These issues remain future GCEstim work; its source is not patched.

## State and updates

The implementation reuses the existing hidden-table state cache pattern. A versioned plain list stores numerical matrices/vectors once, engine version, fit warnings and elapsed time. Image states contain selected indices and presentation settings, not another coefficient path or a ggplot object. CV summaries are derived from the canonical fold-error matrix.

Data/model/grid/CV/error/seed options invalidate the fit cache. Plot settings, table row limits and selection settings do not. Compatible saved state reconstructs tables and plots after reopening without estimation. An incompatible cache version refits once. The current refinement retains cache version 1. Obsolete plotHeatmap, plotSelected and normalize options from older saved analyses are retired using jamovi NONE option-deletion markers (the host retains their names as inactive tombstones, with no active values); obsolete image settings are removed even when a retained plot is hidden. No estimation is needed for this migration. The seed field uses the supported largest UI width; its default, validation and cache dependencies are unchanged. Optional-plot input is validated only when its plot is requested.

## Validation

Run tools/validate-ridge.R using an R development library with the package dependencies. The script supports prepare, core, ridge, test, plots and build arguments. On Windows it discovers the current user R library and an installed jamovi directory. GCESTIMJMV_R_LIBRARY and JAMOVI_HOME can override those paths. Build uses the installed jamovi compiler in build-only mode with dependencies/remotes unchanged, writes a separately named development JMO, and restores pre-existing compiler scratch files with hash checks.

Automated tests cover direct numerical agreement, matrix orientation, tiny custom penalties with rounded-name collisions, numeric/factor models, supported intercept/standardization combinations, CV on/off, CV means/SEs, selection rules and manual lookup, unsupported configurations, nonfinite errors, all plot renderers, optional-plot failures, state serialization/replay, precise invalidation and existing analyses.

## Final manual jamovi checklist

1. Use GCEstimjmv_1.1.0.jmo in a test jamovi session. Confirm WGCE, TARWGCE, TASWGCE and Ridge Trace under Regression. Smoke-test all three existing entropy analyses.
2. Fit y on X001 and X002 with the bundled dataset, then add X005. Check coefficient names, summary, all tables and plot exports.
3. Switch minimum/1-SE/manual grid position. Check the selected coefficients and markers change while the result says Previous fit reused.
4. Toggle every plot, axes, intercept, facets, extrema and display filters. Check legends, captions and larger-model sizing. Invalid display names should skip coefficient plots while preserving tables and performance plots.
5. Disable CV: check the automatic in-sample plot, its Show/Hide override and descriptive selection label. Confirm 1-SE is unavailable and an active stale 1-SE selection gives an informative error.
6. Try custom unsorted, duplicated, zero, negative, nonfinite and malformed grids. Check grid positions and original lambda values; changing grid/model settings must refit.
7. Try both guarded engine configurations. Check that settings stay unchanged and messages identify implementation limitations.
8. Save an .omv with all plots, close/reopen, toggle settings and export plots. Confirm cached reuse, unchanged coefficients and no duplicated table rows.
9. Change a model variable, factor, intercept setting, standardization, CV folds, seed or error measure. Check obsolete results clear and computation reruns. Confirm the full-data preprocessing note appears only with standardization and CV enabled.

## Prior development validation record (2026-10-27)

- Core milestone: 115 assertions passed before extended plot integration.
- Final complete suite: 1,056 assertions passed (267 Ridge Trace plus 789 existing); zero failures, test warnings or skips.
- All nine plot types rendered through jamovi-compatible callbacks, including both axis scales, original/normalized heatmaps, factors and intercept display. PNG exports were visually inspected.
- jamovi 2.7.41 on Windows: definitions/UI compiled and the package installed into the compiler's build library and passed its loading checks. A separately named development JMO was built without installing into the user's jamovi profile or publishing.
- Development artifact: build/ridge-validation/GCEstimjmv-ridge-development.jmo. That prior development artifact has version metadata 1.0.4; it is not the v1.0.5 local test build and the published root artifacts are unchanged.
- Full-suite transcript and per-test CSV are under build/ridge-validation/. The R environment emits pre-existing locale startup warnings and a testthat R-version warning outside the test results.
- GUI interactions, actual .omv save/reopen and user-triggered plot export remain manual checks. State-codec roundtrip and no-refit replay are automated.
- Compiler-only whitespace changes to existing analysis headers were restored to preserve their original source. No existing analysis backend or test was changed.

Run the compiled-package smoke check as a separate fresh process: Rscript --vanilla tools/validate-ridge.R bundled.

## Changed-file inventory

Created: R/ridgetrace.b.R, R/ridgetrace.h.R (generated), R/ridgetrace-state.R; jamovi/ridgetrace.a.yaml, jamovi/ridgetrace.r.yaml, jamovi/ridgetrace.u.yaml; tests/testthat/test-ridgetrace.R; tools/validate-ridge.R; RIDGE-TRACE.md.

Modified: jamovi/0000.yaml (menu integration), DESCRIPTION (analysis description and standard utils namespace declaration), README.md and NEWS.md (development documentation).

## Historical local v1.0.5 preparation

At that development stage, module and R package metadata were 1.0.5. Analysis-interface versions in the individual .a.yaml files are independent schema versions and remain unchanged. The local testing package is GCEstimjmv_1.0.5.jmo in the repository root. Version-specific test, archive and installation evidence is recorded in build/local-1.0.5/. Historical public changelog entries and the earlier development validation record are retained.

### Initial local v1.0.5 validation and installation result (before refinement)

The complete suite passed 1,056 assertions (267 Ridge Trace and 789 existing) without failures, test warnings or skips. The Windows JMO compiled with jamovi 2.7.41, and its embedded jamovi.yaml, jamovi-full.yaml and R DESCRIPTION each report 1.0.5. All four compiled analysis interfaces are included. A fresh-process compiled-package check loaded GCEstim 1.1.0 and matched Ridge Trace selected coefficients to the direct engine.

Local installation completed through the verified jamovi-compiler installer, which invokes jamovi.exe --install. The jamovi GUI subsequently displayed GCEstimjmv 1.0.5 and the four registered analyses. Ridge Trace output for the 75-row dataThesis example was observed with coefficient plots, the full-data preprocessing note and Previous fit reused.

All existing analysis backends, interfaces, tests and the 18 bundled GCEstim files retain their pre-preparation hashes. Detailed evidence is in build/local-1.0.5/POST-INSTALLATION.md and manifest.json. Remaining manual checks include full GUI runs of the three existing entropy analyses, .omv save/reopen, plot export and exhaustive plot/option combinations. No public release, commit or push was performed.

## Local v1.0.5 interface refinement

The refined analysis contains seven plots, with no heatmap or standalone selected-coefficient image result. The Selected Coefficients table is preserved. New analyses display the intercept in coefficient plots by default; saved display preferences remain respected. The random seed input uses jamovi’s supported largest width, with default 230676 and unchanged validation/RNG behaviour. CV references use Minimum error (dotted) and 1-SE rule (dashed); extrema use TARWGCE’s upward/downward triangles.

The complete automated suite passed 1,076 assertions: 287 Ridge Trace and all 789 existing assertions, with zero failures, test warnings or skips. Removal of plot-specific checks reduced the adapted original Ridge tests from 267 to 232; 55 new regression assertions bring Ridge Trace to 287. New coverage verifies reference styles, extrema positions/colours, seven registered plots, intercept/seed defaults, UI width and retirement of obsolete saved options and hidden image settings without refitting. All seven callbacks render on both axis scales with factors. The numerical estimator and selection functions are unchanged.

Validation evidence for this refinement is under build/ridge-refinement-1.0.5/; complete-suite logs and CSV are also under build/local-1.0.5/. Pre-existing R locale startup warnings and the testthat build-version warning occur outside the test results.

The refined JMO compiled with jamovi 2.7.41 and passed embedded-version and fresh-process numerical/loading checks. Local reinstallation completed through jamovi-compiler’s supported installer (jamovi.exe --install). In the actual jamovi GUI, the installed module reports 1.0.5, the four-analysis menu remains available, both removed plot controls are absent, and all seven remaining plots rendered. CV reference styles, extrema, selected estimates in ranges, concise 1-SE labels, the checked intercept-display default and an unclipped ten-digit seed were observed. The default seed was restored after the width check.

A pre-refinement .omv was reopened, refreshed through a display change without fitting, and saved as a separate migrated copy. The structured saved-state audit confirms seven completed PNG results, no obsolete image results, three inactive NONE tombstones instead of old option values, and a byte-for-byte unchanged fit-cache state. Run Rscript --vanilla tools/refine-ridge.R saved-refine to repeat the audit against the two local evidence files. Remaining manual checks are user-triggered plot export, other display scaling, and exhaustive datasets/option combinations.

Refinement file inventory and automated/GUI evidence: build/ridge-refinement-1.0.5/REPORT.md. The JMO SHA-256 is 3f49f263992b9f5c4e7b444ead507a48e1b5010e91eea767ee0b0072917774ce.

## Compatibility and methodological limitations

Coefficient colours use internal identifiers; display labels remove formula backticks and disambiguate collisions. Coefficient-path and entropy legends wrap at existing dimensions. Compatible saved analyses reuse numerical results; old embedded images can be refreshed by toggling **OLS comparison** on and off, then saving. An explicit notice appears if legacy entropy data are insufficient for complete recovery.

Ridge coefficient ranges, sensitivity and fold-training stability are descriptive, not confidence intervals. Standardized ridge CV currently uses full-data preprocessing before folds; unsupported engine configurations produce informative errors. The deferred upstream GCEstim support-selection 1-SE concern remains unresolved: TARWGCE/TASWGCE preserve the engine’s returned selections. The standalone Ridge Trace 1-SE lambda rule is derived separately from returned fold errors.

This GitHub release is available by sideloading; it does not imply approval or availability of v1.1.0 in the jamovi Library.
