# GCEstim for jamovi

**GCEstimjmv** brings entropy-based linear regression from the [GCEstim R package](https://CRAN.R-project.org/package=GCEstim) to jamovi. It provides four analyses in **Analyses → GCEstim → Regression**:

| Analysis | Signal support | Main use |
| --- | --- | --- |
| **W-GCE** | Limits entered by the user | Fit a weighted GCE linear model with specified signal supports. |
| **TARW-GCE** | Ridge-informed, with symmetric or asymmetric supports | Compare candidate support sizes and weights by cross-validation. |
| **TASW-GCE** | Based on standardized coefficients | Compare candidate support sizes and weights by cross-validation. |
| **Ridge Trace** | Standalone ridge coefficient paths | Inspect penalties, prediction errors and coefficient stability using GCEstim::ridgetrace(). |

The adaptive analyses use `GCEstim::cv.lmgce()`; W-GCE uses `GCEstim::lmgce()`. The three entropy analyses accept numeric covariates and factors, report coefficient estimates and estimated signal support probabilities (*p*), and offer confidence intervals and optional saved columns for fitted values, residuals, and estimated noise support probabilities (*w*). TARW-GCE and TASW-GCE additionally provide a cross-validation results table and diagnostic plots.

## Version 1.1.0

The registered analyses are **WGCE**, **TARWGCE**, **TASWGCE**, and **RidgeTrace** (displayed as **Ridge Trace**). This release includes standalone Ridge Trace and graphical/state compatibility corrections. See [release notes](NEWS.md).

## Install on Windows

1. Download [`GCEstimjmv_1.1.0.jmo`](https://github.com/jorgevazcabral/GCEstim-jamovi/releases/download/v1.1.0/GCEstimjmv_1.1.0.jmo) from the [v1.1.0 release](https://github.com/jorgevazcabral/GCEstim-jamovi/releases/tag/v1.1.0). This Windows x64 build is validated with jamovi 2.7.41.
2. In jamovi, select **Modules (+) → Side-load**.
3. **Click the upload icon**, then choose the `.jmo` file in the file picker. Dragging the file into the jamovi window may attempt to open it as a data file instead.
4. Find **GCEstim** in the Analyses ribbon.

Users of a compatible prebuilt `.jmo` do **not** need to install R or the GCEstim R package separately. The `.jmo` bundles the R packages used by the module. Builds are specific to an operating system, architecture, and jamovi series; this Windows build cannot be used on macOS or Linux.

## Run an analysis

1. Open a data set in jamovi, select **GCEstim** in the Analyses ribbon, and choose **W-GCE**, **TARW-GCE**, **TASW-GCE**, or **Ridge Trace**.
2. Set the dependent variable and add covariates or factors. The model uses complete cases for the selected variables.
3. For an entropy analysis, review signal and noise supports. For Ridge Trace, review the penalty grid and CV options. Select **Fit and update model**.
4. Inspect the model summary and coefficient estimates. The entropy analyses also report signal probabilities (*p*) and offer **Save** options. Ridge Trace reports coefficient ranges and optional lambda results, with seven graphical outputs.

For W-GCE, **Lower limits** and **Upper limits** accept one number for all coefficients or comma-separated values for individual coefficients. When supplying individual values, include the intercept if selected and one value for each non-reference factor level. For TARW-GCE and TASW-GCE, enter comma-separated candidate signal points, noise points, and weights, for example `3, 5, 7, 9` and `0.1, 0.3, 0.5, 0.7, 0.9`. Each signal or noise support point count must be an odd integer of at least 3; a list can contain just one candidate value. The adaptive analyses offer minimum-error and one-standard-error rules for selecting the cross-validation result. Larger candidate grids and bootstrap settings can increase computation time substantially.

Saving noise support probabilities (*w*) as data columns is available for the data-based GLM case; the moment-based cases do not provide one probability row per observation.

## Example data

The module includes [`dataThesis.csv`](data/dataThesis.csv). It contains 75 observations, numeric predictors `X001`–`X004`, a three-level factor `X005` (A, B, C), and a continuous outcome `y`.

After installing the module, open jamovi's **main menu → Open → Data Library** and select **dataThesis**. To try any of the four analyses, set **Dependent Variable** to `y`, add `X001`–`X004` as **Covariates**, and optionally add `X005` as a **Factor**. Select **Fit and update model**. For a quicker first run of TARW-GCE or TASW-GCE, use `3, 5` for signal support points, `3` for noise support points, `0.5` for the noise weight, and `2` signal support spaces.

## Build from source

Building requires a compatible installation of jamovi, R, and [`jmvtools`](https://github.com/jamovi/jmvtools). With this repository open as an RStudio project, run in the R console (adjust the jamovi path to your installation):

```r
options(jamovi_home = "C:/Program Files/jamovi 2.7.41.0")  # adjust if needed

jmvtools::prepare()
jmvtools::install()
```

The package dependency is declared in `DESCRIPTION` and pinned to `GCEstim` v1.1.0. `jmvtools::install()` downloads and bundles dependencies into a platform-specific `.jmo` in the project directory. It may take time on a fresh build. The source does not rely on a personal R library path.

## References and support

- Cabral, J. (2026). *GCEstim: Regression Coefficients Estimation Using the Generalized Cross Entropy*. [R package](https://CRAN.R-project.org/package=GCEstim).
- [GCEstim source](https://github.com/jorgevazcabral/GCEstim)
- [Report a module issue](https://github.com/jorgevazcabral/GCEstim-jamovi/issues)

License: GPL-3. See [`DESCRIPTION`](DESCRIPTION) for package metadata.

## Standalone Ridge Trace

The source now includes **GCEstim → Regression → Ridge Trace**, using only `GCEstim::ridgetrace()` for estimation. It provides selected coefficients, coefficient ranges, optional lambda results and seven plots. Minimum CV error, a GCEstimjmv-derived 1-SE rule, and manual selection by ascending grid position reuse the evaluated coefficient path without refitting. Only RMSE, MSE, MAE and sMAPE are exposed.

See [Ridge Trace implementation and validation](RIDGE-TRACE.md) for controls, statistical limitations, cache behaviour and manual checks. Ridge Trace is included in v1.1.0.

## Compatibility and methodological limitations

Coefficient colours use internal identifiers; display labels remove formula backticks and disambiguate collisions. Coefficient-path and entropy legends wrap at existing dimensions. Compatible saved analyses reuse numerical results; old embedded images can be refreshed by toggling **OLS comparison** on and off, then saving. An explicit notice appears if legacy entropy data are insufficient for complete recovery.

Ridge coefficient ranges, sensitivity and fold-training stability are descriptive, not confidence intervals. Standardized ridge CV currently uses full-data preprocessing before folds; unsupported engine configurations produce informative errors. The deferred upstream GCEstim support-selection 1-SE concern remains unresolved: TARWGCE/TASWGCE preserve the engine’s returned selections. The standalone Ridge Trace 1-SE lambda rule is derived separately from returned fold errors.

This GitHub release is available by sideloading; it does not imply approval or availability of v1.1.0 in the jamovi Library.
