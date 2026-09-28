# GCEstim for jamovi

**GCEstimjmv** brings entropy-based linear regression from the [GCEstim R package](https://CRAN.R-project.org/package=GCEstim) to jamovi. It provides three analyses in **Analyses → GCEstim → Regression**:

| Analysis | Signal support | Main use |
| --- | --- | --- |
| **W-GCE** | Limits entered by the user | Fit a weighted GCE linear model with specified signal supports. |
| **TARW-GCE** | Ridge-informed, with symmetric or asymmetric supports | Compare candidate support sizes and weights by cross-validation. |
| **TASW-GCE** | Based on standardized coefficients | Compare candidate support sizes and weights by cross-validation. |

The adaptive analyses use `GCEstim::cv.lmgce()`; W-GCE uses `GCEstim::lmgce()`. All three accept numeric covariates and factors, report coefficient estimates and estimated signal support probabilities (*p*), and offer confidence intervals and optional saved columns for fitted values, residuals, and estimated noise support probabilities (*w*). TARW-GCE and TASW-GCE additionally provide a cross-validation results table and diagnostic plots.

## Install on Windows

1. Download the Windows `.jmo` file from [Releases](https://github.com/jorgevazcabral/GCEstim-jamovi/releases), if a build for your jamovi version is available. The initial tested build is named `GCEstimjmv_0.1.0.jmo` and was built with jamovi 2.7.33.0 for Windows x64.
2. In jamovi, select **Modules (+) → Side-load**.
3. **Click the upload icon**, then choose the `.jmo` file in the file picker. Dragging the file into the jamovi window may attempt to open it as a data file instead.
4. Find **GCEstim** in the Analyses ribbon.

Users of a compatible prebuilt `.jmo` do **not** need to install R or the GCEstim R package separately. The `.jmo` bundles the R packages used by the module. Builds are specific to an operating system, architecture, and jamovi series; a Windows build cannot be used on macOS or Linux.

## Run an analysis

1. Open a data set in jamovi, select **GCEstim** in the Analyses ribbon, and choose **W-GCE**, **TARW-GCE**, or **TASW-GCE**.
2. Set the dependent variable and add covariates or factors. The model uses complete cases for the selected variables.
3. Review the signal and noise support options, then click **Run model**.
4. Inspect the model summary, coefficient estimates, and estimated signal probabilities (*p*). Select the **Save** options to add fitted values, residuals, or noise probabilities (*w*) to the data set.

For W-GCE, **Lower limits** and **Upper limits** accept one number for all coefficients or comma-separated values for individual coefficients. For TARW-GCE and TASW-GCE, enter comma-separated candidate signal points, noise points, and weights, for example `3, 5, 7, 9` and `0.1, 0.3, 0.5, 0.7, 0.9`. Signal and noise point counts must be odd integers of at least 3. The adaptive analyses offer `min` and `1se` rules for selecting the cross-validation result. Larger candidate grids and bootstrap settings can increase computation time substantially.

## Build from source

Building requires a compatible installation of jamovi, R, and [`jmvtools`](https://github.com/jamovi/jmvtools). With this repository open as an RStudio project, run in the R console (adjust the jamovi path to your installation):

```r
options(jamovi_home = "C:/Program Files/jamovi 2.7.33.0")

jmvtools::prepare()
jmvtools::install()
```

The package dependency is declared in `DESCRIPTION` and pinned to `GCEstim` v1.1.0. `jmvtools::install()` downloads and bundles dependencies into a platform-specific `.jmo` in the project directory. It may take time on a fresh build. The source does not rely on a personal R library path.

## References and support

- Cabral, J. (2026). *GCEstim: Regression Coefficients Estimation Using the Generalized Cross Entropy*. [R package](https://CRAN.R-project.org/package=GCEstim).
- [GCEstim source](https://github.com/jorgevazcabral/GCEstim)
- [Report a module issue](https://github.com/jorgevazcabral/GCEstim-jamovi/issues)

License: GPL-3. See [`DESCRIPTION`](DESCRIPTION) for package metadata.
