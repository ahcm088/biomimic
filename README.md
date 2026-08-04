# biomimic

<!-- badges: start -->
[![R-CMD-check](https://github.com/ahcm088/biomimic/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ahcm088/biomimic/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R >= 4.1](https://img.shields.io/badge/R-%3E%3D%204.1-blue.svg)](https://cran.r-project.org/)
<!-- badges: end -->

Latent differential analysis for high-dimensional biological data: screen
thousands of candidate variables, then fit a structural equation model to the
ones that matter — without needing to write SEM syntax.

The usual approach to antibodyome and proteomics data borrows from
transcriptomics: test each protein independently, then correct for multiple
testing. That works for large marginal effects but discards the correlation
structure, which is often where the biology is. A MIMIC model asks a different
question — *does the latent immunological state differ between groups?* —
and answers it with a single structural coefficient, so no multiple-testing
correction is needed.

The obstacle is practical: general-purpose SEM software takes the indicator set
as given, and choosing 15 indicators out of 8,000 is not something it does.
`biomimic` closes that gap with a three-stage pipeline.

## Installation

`biomimic` is not on CRAN. Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("ahcm088/biomimic")
```

Or with `pak`:

```r
# install.packages("pak")
pak::pak("ahcm088/biomimic")
```

To build the vignette locally, add `build_vignettes = TRUE`:

```r
remotes::install_github("ahcm088/biomimic", build_vignettes = TRUE)
```

Requires R >= 4.1.0 and [`lavaan`](https://lavaan.ugent.be) (>= 0.6-17), which
is installed automatically.

## Quick start

The package ships a subset of the neurodegenerative antibodyome data used in
the paper (Alzheimer's disease vs. control), already model-ready.

```r
library(biomimic)
data(neuro_antibodyome)

fit <- biomimic(neuro_antibodyome$proteins,  # 137 samples x 400 features
                neuro_antibodyome$X,         # intercept + group + age + sex
                K = 1,                       # latent factors
                k = 15)                      # indicators to retain
summary(fit)
```

```
=== biomimic SEM Results ===
Model: MIMIC | Estimator: ML | n = 137

--- Model Fit ---
  Chi-sq = 380.7 (df = 132, p = < 0.001)
  RMSEA  = 0.117 [0.104, 0.131]  [Poor fit]
  CFI    = 0.938  [Acceptable]
  SRMR   = 0.029  [Excellent fit]

--- Group Effects on Latent Pathways ---
  group -> Latent pathway: est = -0.849, SE = 0.196, z = -4.33, p = < 0.001 ***  [-1.233, -0.465]
  age -> Latent pathway: est = 0.362, SE = 0.088, z = 4.11, p = < 0.001 ***  [0.190, 0.534]
  sex -> Latent pathway: est = -0.489, SE = 0.149, z = -3.29, p = 0.001 **  [-0.781, -0.197]

--- Protein-Pathway Associations (Factor Loadings) ---
  SSR4 <- Latent pathway: 1.000 (std = 0.960)
  ZDHHC7 <- Latent pathway: 1.038 (std = 0.939) ***
  RAMP1 <- Latent pathway: 0.923 (std = 0.952) ***
  ...
```

One number carries the group comparison: the latent immunological state is
lower in AD than in controls (est = -0.849, p < 0.001), adjusted for age and
sex. The whole run takes under two seconds. The full `lavaan` object remains
available via `fit$lavaan_fit` for anything the wrapper does not cover.

## What the pipeline does

| Stage | Function | Purpose |
|---|---|---|
| 1. Screening | `vb_ard()` | Ranks all *p* variables by posterior signal-to-noise ratio using a variational Bayes factor model with ARD priors. The Gamma extension simultaneously flags variables whose covariate association bypasses the latent factor. |
| 2. Selection | `top_k()`, `rank_variables()`, `discover_factors()` | Retains the top *k* variables and, for multi-factor models, discovers block structure via SVD with varimax rotation. |
| 3. Estimation | `fit_sem()`, `build_lavaan_model()` | Generates and fits the `lavaan` model: MIMIC, multi-group CFA, MIMIC with direct effects, or multi-factor. |

`biomimic()` chains all three stages in one call. The result is a
`biomimic_fit` object with `print()`, `summary()`, `coef()`, `confint()` and
`plot()` methods that report in biological language rather than SEM jargon.

Post-estimation tools: `goodness_of_fit()`, `test_hypotheses()`, `lr_test()`,
`test_invariance()` (measurement invariance across groups), `diagnose()`,
`residual_analysis()`, `mahalanobis_distance()`, `predict_scores()`,
`predict_observed()` and `network_analysis()`.

Plotting helpers include `plot_screening()`, `plot_loadings()`, `plot_path()`,
`plot_scores()`, `plot_structural_effects()`, `plot_group_effects()`,
`plot_mediation()`, `plot_residuals()`, `plot_mahalanobis()` and
`plot_selection_recurrence()`, themed by `theme_biomimic()` and
`palette_biomimic()`.

Simulation utilities: `simulate_mimic()`, `run_simulation()`, `sim_metrics()`
and `sim_screening_metrics()`.

## Documentation

```r
?biomimic                  # pipeline wrapper
?vb_ard                    # Stage 1 screening
?fit_sem                   # Stage 3 estimation
vignette("biomimic")       # worked example, end to end
help(package = "biomimic") # full index
```

## Data

`neuro_antibodyome` is a 400-feature subset (n = 137: 46 AD cases, 91 controls)
of the ProtoArray data of DeMarshall et al. (2015), GEO accession
[GSE62283](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62283). It
keeps the features selected in the paper plus a reproducible random background
sample. The full 8,323-feature matrix is available from the accession above.

## Citation

```r
citation("biomimic")
```

The manuscript describing the method is in preparation. Until it appears,
please cite the package.

## License

MIT © the package authors. See [LICENSE.md](LICENSE.md).
