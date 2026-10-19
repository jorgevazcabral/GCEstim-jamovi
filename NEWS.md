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
