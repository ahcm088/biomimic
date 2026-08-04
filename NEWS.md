# biomimic 0.2.0

First public release.

## Stage 1 — VB-ARD screening

* `vb_ard()` ranks all *p* variables by posterior signal-to-noise ratio using a
  variational Bayes factor model with automatic relevance determination priors,
  scaling to thousands of candidates without tuning.
* The Gamma extension places ARD priors on direct-effect coefficients as well,
  flagging variables whose covariate association bypasses the latent factor.
  This automates a step normally done by hand from modification indices.

## Stage 2 — Selection and structure discovery

* `top_k()` and `rank_variables()` retain the top *k* variables.
* `discover_factors()` recovers block structure via SVD with varimax rotation
  for multi-factor models.

## Stage 3 — SEM estimation

* `build_lavaan_model()` generates `lavaan` syntax for four model types: MIMIC,
  multi-group CFA, MIMIC with direct effects, and multi-factor MIMIC.
* `fit_sem()` fits the model with sensible defaults for biological data and
  returns a `biomimic_fit` object whose `print()`, `summary()`, `coef()`,
  `confint()` and `plot()` methods report in biological rather than SEM terms.
* `biomimic()` chains all three stages in one call.

## Post-estimation

* `goodness_of_fit()`, `test_hypotheses()` and `lr_test()` for model assessment
  and nested comparisons.
* `test_invariance()` for configural, metric and scalar measurement invariance.
* `diagnose()`, `residual_analysis()`, `mahalanobis_distance()` for diagnostics.
* `predict_scores()`, `predict_observed()` and `network_analysis()` for
  downstream use.

## Data, plotting and simulation

* `neuro_antibodyome`: 400-feature subset of the neurodegenerative antibodyome
  data (GEO `GSE62283`, AD vs. control), ready for the pipeline.
* Plot helpers for screening, loadings, path diagrams, scores, structural and
  group effects, mediation, residuals, Mahalanobis distances and selection
  recurrence, with `theme_biomimic()` and `palette_biomimic()`.
* `simulate_mimic()`, `run_simulation()`, `sim_metrics()` and
  `sim_screening_metrics()` for method evaluation.
