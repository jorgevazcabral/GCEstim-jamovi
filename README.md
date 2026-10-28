# GCEstim for jamovi

**GCEstimjmv** brings entropy-based linear regression from the [GCEstim R package](https://CRAN.R-project.org/package=GCEstim) to jamovi. It provides four analyses under **Analyses → GCEstim → Regression**, including a standalone Ridge Trace analysis for exploring ridge penalties and coefficient paths.

Current published version: **1.1.2**. See the [GitHub release](https://github.com/jorgevazcabral/GCEstim-jamovi/releases/tag/v1.1.2) and [release notes](NEWS.md).

## Main features

- Entropy-based linear models with user-specified, ridge-informed or standardized-coefficient-based signal supports.
- Numeric covariates and factors, with complete-case analysis for the selected variables.
- Cross-validation and diagnostic plots for the adaptive entropy analyses and standalone Ridge Trace.
- Coefficient estimates, signal support probabilities (*p*), confidence intervals and optional saved fitted values, residuals and noise support probabilities (*w*) for the entropy analyses.
- Reuse of compatible numerical results when changing presentation or reopening saved analyses.

## Available analyses

| Analysis | jamovi menu label | Purpose |
| --- | --- | --- |
| **WGCE — Weighted Generalized Cross Entropy** | **W-GCE** | Fit a weighted GCE linear model with user-specified signal support limits. |
| **TARWGCE — Two-Step Adaptive Ridge-Weighted GCE** | **TARW-GCE** | Use ridge-informed symmetric or asymmetric supports and compare candidate support sizes and weights by cross-validation. |
| **TASWGCE — Two-Step Adaptive Standardized-Weighted GCE** | **TASW-GCE** | Construct supports from standardized coefficients and compare candidate support sizes and weights by cross-validation. |
| **Ridge Trace** | **Ridge Trace** | Explore ridge coefficient paths, evaluated penalties, prediction errors and fold-training coefficient stability. |

The registered analysis identifiers are `WGCE`, `TARWGCE`, `TASWGCE` and `RidgeTrace`. W-GCE uses `GCEstim::lmgce()`; the adaptive entropy analyses use `GCEstim::cv.lmgce()`. These R function names do not identify additional module analyses.

For W-GCE, **Lower limits** and **Upper limits** accept one number for all coefficients or comma-separated values for individual coefficients. Include the intercept, if selected, and one value for each non-reference factor level when supplying individual values. TARW-GCE and TASW-GCE accept comma-separated candidate signal points, noise points and weights, such as `3, 5, 7, 9` and `0.1, 0.3, 0.5, 0.7, 0.9`. Support point counts must be odd integers of at least 3; a list may contain a single candidate. Larger candidate grids and bootstrap settings can increase computation time substantially.

Saving noise support probabilities (*w*) as data columns is available for the data-based GLM case; moment-based cases do not provide one probability row per observation.

## Installation

### Preferred: jamovi Library, when available

In jamovi, open **Modules (+) → jamovi library**, find **GCEstimjmv** and select **Install** when it is listed for your platform and jamovi version. Availability of **v1.1.2** in the Library is not confirmed here; the GitHub release does not imply Library approval or availability. If the module or required version is unavailable, use the compatible GitHub build below.

### Alternative: GitHub release

1. Open the official [v1.1.2 release](https://github.com/jorgevazcabral/GCEstim-jamovi/releases/tag/v1.1.2) and download [`GCEstimjmv_1.1.2.jmo`](https://github.com/jorgevazcabral/GCEstim-jamovi/releases/download/v1.1.2/GCEstimjmv_1.1.2.jmo).
2. In jamovi, select **Modules (+) → Side-load**.
3. **Click the upload icon** and choose the `.jmo` file. Dragging it into the main window may attempt to open it as a data file instead.
4. Find **GCEstim** in the Analyses ribbon.

The published JMO is a **Windows x64** build for the tested jamovi environment described below. A compatible prebuilt module includes its required R dependencies through jamovi's module distribution mechanism: ordinary jamovi users do **not** need to install R or GCEstim separately. Builds are specific to an operating system, architecture and jamovi series; this Windows build cannot be used on macOS or Linux.

## Standalone Ridge Trace

Ridge Trace uses **`GCEstim::ridgetrace()`** to explore how coefficients and errors change over an evaluated lambda grid. It complements the entropy-based estimators by providing ridge diagnostics and penalty selection; it does not replace WGCE, TARWGCE or TASWGCE.

Its controls and results include:

- Geometric lambda grids or custom comma-separated positive lambda values.
- Cross-validation with configurable folds, random seed and RMSE, MSE, MAE or sMAPE error measures.
- Minimum-error and 1-SE selection rules, plus **Manual (grid position)** selection from the evaluated grid. The module derives the 1-SE choice from returned fold errors.
- Model Summary and numeric Selected Lambda summaries, selected coefficient estimates, coefficient ranges and an optional Lambda Results table.
- Coefficient filtering using displayed labels or internal identifiers, including quoted names containing commas; coefficient colours remain consistent across plots and filtering.

Seven graphical outputs are available:

1. Ridge Coefficient Paths.
2. Coefficient Ranges.
3. Local Coefficient Sensitivity.
4. Cross-Validation Error.
5. In-Sample Error.
6. Fold-Specific CV Errors.
7. Coefficient Stability Across CV Folds.

Coefficient paths and CV error are shown by default when CV is enabled; other plots can be enabled in **Plots**. The in-sample plot appears automatically when CV is disabled unless its display setting is overridden. Changing the selected evaluated lambda or presentation reuses the existing numerical path without refitting.

See [Ridge Trace controls, methodology and validation history](RIDGE-TRACE.md) for details. Historical version and development-validation sections in that document describe earlier builds; the current release is v1.1.2.

## Quick start with dataThesis

The bundled [`dataThesis.csv`](data/dataThesis.csv) contains 75 observations, numeric predictors `X001`–`X004`, a three-level factor `X005` (A, B, C) and a continuous outcome `y`.

1. After installing the module, open jamovi's **main menu → Open → Data Library** and select **dataThesis**.
2. Open **Analyses → GCEstim → Regression → Ridge Trace**.
3. Set **Dependent Variable** to `y`, add `X001`–`X004` as **Covariates**, and optionally add `X005` under **Factors**.
4. Under **Cross-Validation and Selection**, keep **Cross-validation** enabled. The defaults are five folds, seed `230676` and minimum-error selection. Review the **Ridge Penalty** grid, then select **Fit and update model**.
5. Inspect **Selected Lambda**, the coefficient estimates and **Ridge Coefficient Paths**. Try the **1-SE rule** or **Manual (grid position)** to compare evaluated penalties; enable additional outputs under **Plots**.

The same dataset can be used with any entropy analysis. For a quicker first TARW-GCE or TASW-GCE run, use `3, 5` signal support points, `3` noise support points, `0.5` noise weight and `2` signal support spaces.

## Documentation

- [Release notes](NEWS.md).
- [Ridge Trace controls, methodology and validation history](RIDGE-TRACE.md).
- [Package metadata and dependencies](DESCRIPTION).
- [v1.1.1 audit corrections and validation report](AUDIT-LOW-v1.1.1.txt); v1.1.2 subsequently removed two redundant table notes without changing estimation.

### Build from source (developers)

Building requires a compatible installation of jamovi, R and [`jmvtools`](https://github.com/jamovi/jmvtools). With this repository open as an RStudio project, run in the R console, adjusting the jamovi path:

```r
options(jamovi_home = "C:/Program Files/jamovi 2.7.41.0")

jmvtools::prepare()
jmvtools::install()
```

`DESCRIPTION` declares `GCEstim (>= 1.1.0)` and pins the remote dependency to `GCEstim@v1.1.0`. `jmvtools::install()` downloads and bundles dependencies into a platform-specific `.jmo`; a fresh build may take time. The source does not rely on a personal R library path. These development requirements do not apply to users installing a compatible prebuilt JMO.

## Compatibility and limitations

- **Module:** v1.1.2. **Tested environment:** jamovi 2.7.41 on Windows. The published asset is Windows x64; testing on macOS or Linux is not documented.
- **Dependency:** GCEstim v1.1.0 is used by the distributed build. The module's maintenance releases do not modify the underlying GCEstim R package.
- **Standardized Ridge Trace CV:** GCEstim preprocesses the complete dataset before constructing CV folds. Results retain that behaviour rather than independent fold-specific preprocessing; an INFO notice explains the limitation.
- **Unsupported configuration:** standardized CV with an intercept and only one expanded predictor column is rejected, including a binary factor yielding that shape. Without an intercept, standardization is disabled and the effective value is FALSE; its stored preference is retained.
- **CV disabled:** automatic selection uses minimum returned in-sample error and is labelled descriptive. A stored 1-SE preference falls back to that rule with an INFO notice and is restored when CV is re-enabled.
- **Interpretation:** coefficient ranges, local sensitivity and fold-training stability are descriptive diagnostics, not coefficient confidence intervals. CV fold SE bars are not confidence intervals for generalization error. Unavailable or nonfinite numerical results generate informative warnings.
- **Saved analyses:** compatible numerical results are reused. Legacy entropy images can be refreshed by toggling **OLS comparison** and saving; genuine recovery failures produce a warning while available graphical results are retained.

## Citation, references and support

- Cabral, J. (2026). *GCEstim: Regression Coefficients Estimation Using the Generalized Cross Entropy*. [R package](https://CRAN.R-project.org/package=GCEstim).
- [GCEstim source](https://github.com/jorgevazcabral/GCEstim).
- [Report a module issue](https://github.com/jorgevazcabral/GCEstim-jamovi/issues).

Author and maintainer: Jorge Cabral, Department of Mathematics, University of Aveiro. License: **GPL-3**. See [`DESCRIPTION`](DESCRIPTION) for package metadata.
