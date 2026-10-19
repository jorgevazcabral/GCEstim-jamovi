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
- Cache format v3 replaces earlier caches; an older saved analysis may need
  one refit before its new cache can be reused.

Validation: `devtools::test()` passed 208 assertions with no failures, test
warnings, or skips. The Windows x64 module built and installed successfully
with jamovi 2.7.33.0. Manual jamovi GUI checks passed for `.omv` saving and
reopening, cache reuse, plot and OLS toggles, normalized entropy, and saved
columns.

Release asset: `GCEstimjmv_1.0.3.jmo` (Windows x64, jamovi 2.7.33.0).
SHA-256: `c13adfdaf80c7ddad9bf9ec5d679597be0b79ca85e9f04d72ba1ca4695233aca`.
