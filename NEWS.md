# GCEstimjmv 1.1.2

- Removed redundant explanatory notes from the Ridge Trace Selected Lambda and Selected Coefficients tables.
- No changes to statistical estimation or analysis functionality.

# GCEstimjmv 1.1.1

## Audit-driven bug fixes and usability improvements

- Give every standalone Ridge Trace coefficient its own jamovi theme colour, including models with seven or more coefficients. Keep colours consistent across plots, filtering and intercept visibility.
- Accept displayed coefficient labels and complete internal identifiers in the coefficient filter, including quoted names containing commas. Reject unknown or ambiguous names once, listing valid displayed names before plot construction.
- Keep requested plots visible after build failures with an error using the plot title. Report unexpected rendering failures instead of silently returning a blank image.
- Reuse shared coefficient label and palette helpers. Preserve statistical estimation, lambda selection, cross-validation, numerical cache format and all other analyses.
- Restrict each Ridge Trace plot's invalidation and saved presentation state to the options it reads; retain unaffected images without rebuilding them.
- Disable standardization without an intercept and pass the effective value FALSE while retaining the stored preference. With CV disabled, retain a stored 1-SE preference but use descriptive minimum in-sample error with an INFO notice.
- Show minimum, maximum and selected lambda in numeric result cells, with a concise coefficient-table title and full internal precision.
- Keep statistical explanations in table footnotes and the full-data CV preprocessing notice. Remove general captions from exported Ridge Trace images; move genuine legacy entropy-recovery failures to WARNING notices.
- Remove development-only Ridge Trace migration helpers and their internal jmvcore API dependency; preserve released numerical caches and saved image compatibility.
- Use title-case Ridge Trace plot titles, “CV mean ± 1 SE” and “Manual (grid position)” labels consistently.

# GCEstimjmv 1.1.0

## New features

- Add standalone Ridge Trace using `GCEstim::ridgetrace()`, automatic geometric or custom penalty grids, configurable cross-validation, minimum-error and module-derived 1-SE lambda selection, and manual selection by grid position.
- Add model summaries, selected coefficients, coefficient ranges and optional lambda results, plus seven coefficient, prediction-error and stability plots.
- Reuse evaluated ridge results when changing lambda selection or presentation.

## Improvements

- Use deterministic named coefficient colours from jamovi’s categorical palette across each fitted model’s plots.
- Format coefficient and factor labels for presentation only, safely disambiguating collisions while retaining complete internal identifiers for joins and filtering.
- Wrap coefficient-path, TARWGCE ridge and entropy legends adaptively at existing dimensions; retain bottom legends without unnecessary titles.
- Improve current and legacy saved-state reconstruction and graphical refreshes without unnecessary estimation when compatible cached data are available.

## Bug fixes

- Preserve coefficient identities on selected support-estimate and entropy points in TARWGCE/TASWGCE, including reordered predictors and factors.
- Restore valid coefficient entropy observations omitted after failed support configurations using verified same-stage cached data; show an explicit notice when legacy recovery data are insufficient.
- Apply named coefficient scales after jamovi themes to prevent palette overrides.
- Match standalone Ridge Trace selected-lambda line types to the active rule: minimum error dotted, 1-SE dashed and manual solid. Preserve simultaneous CV references.
- Preserve coefficient-specific triangular extrema and correct selected estimates.
- Fix coefficient-path legend clipping for multiple or long coefficient names.

## Compatibility

- Preserve WGCE, TARWGCE and TASWGCE interfaces, numerical results, estimation, cross-validation, optimization and existing caching architecture.
- Keep the GCEstim dependency pinned to v1.1.0; no GCEstim package changes.
- Retain informative Ridge Trace errors for unsupported engine configurations, document full-data preprocessing during standardized CV, and exclude MAPE/MASE.
- The deferred upstream GCEstim support-selection 1-SE concern remains unresolved. Returned entropy-analysis support selections are preserved; this release does not claim to correct that upstream matter.

Windows x64 release asset: `GCEstimjmv_1.1.0.jmo` (jamovi 2.7.41).

# GCEstimjmv 1.0.4

- Cache W-GCE's compact results and confidence intervals, including bootstrap
  intervals and OLS comparison data. Plot and OLS toggles reuse the fit, and
  saved columns and reopened analyses can be reconstructed from plain state.
- Standardize Model Summary CI types as Asymptotic normal, signal probability
  headings as p (point 1), p (point 2), etc., and adaptive support-point hints.
- Validate true coefficients only for selected, available precision-error plots.
  Invalid optional input produces a warning and skips those plots while retaining
  tables and other plots, including after factor or intercept changes.
- Store adaptive fit caches as inspectable lists (format v4), without nested
  serialization, gzip packing, or duplicate copies of ordinary image data.
  Older cache formats refit once; current caches replay saved outputs without
  refitting. Precision plots use compact coefficient data so changing true
  coefficients or plot options never requires cross-validation.
- Invalidate W-GCE model results when signal or noise prior options change.

Validation: 789 automated assertions passed, including compact cache serialization,
replay, bootstrap reuse, saved-output alignment, optional precision input, and
result dependencies. Windows x64 build for jamovi 2.7.33.0.

Release asset: `GCEstimjmv_1.0.4.jmo`.

# GCEstimjmv 1.0.3

- TARW-GCE and TASW-GCE now cache compact results and plot data so fitted
  results can be reused after saving and reopening an `.omv` file. Serialized
  fit caches contain no `lmgce` or `ggplot` objects, functions, or environments.
- Plot visibility, OLS comparisons, and normalized-entropy options can be
  updated using cached results without repeating cross-validation.
- Saved fitted values, residuals, and noise probabilities (`w`) remain
  available from the cache, with correct alignment for incomplete cases.
- Restore OLS reference lines and normalized-entropy confidence limits in
  diagnostic plots, and retain plot warnings when cached results are reused.
- Fix TARW-GCE precision plots when a candidate support fit is unsuccessful.
- Keep TARW-GCE coefficient and signal-probability tables when presentation
  options change; avoid clearing CV results and ridge traces for unrelated
  options.
- Prevent duplicate TARW-GCE table rows and probability columns when retained
  results are reused.
- Use consistent SE labels, result headings, and optimization-method labels.
  Remove duplicate in-plot titles, describe CI methods in words on the estimate
  axis, clarify the true-coefficient label, and group cross-validation controls
  with the model settings.
- Cache format v3 replaces earlier caches; an older saved analysis may need
  one refit before its new cache can be reused.

Validation: `devtools::test()` passed 225 assertions with no failures, test
warnings, or skips. The Windows x64 module built and installed successfully
with jamovi 2.7.33.0. Manual jamovi GUI checks passed for `.omv` saving and
reopening, cache reuse, plot and OLS toggles, normalized entropy, and saved
columns. The rebuilt module also passed manual GUI checks after the label,
panel-order, result-invalidation, and retained-table fixes.

Release asset: `GCEstimjmv_1.0.3.jmo` (Windows x64, jamovi 2.7.33.0).
SHA-256: `f0115479ed04f6b4cb61214c40272e7ecf6bcb7798f966df8bacaa8d741267d4`.
