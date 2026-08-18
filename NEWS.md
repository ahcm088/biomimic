# biomimic 0.3.0

Identification release. The simulation study and all five applications have
been re-run under this version's protocol (ARD on for every K, two-pass
anchoring, bootstrap-calibrated direct-effect thresholds).

## Identification of B versus Gamma

Marginalising the latent scores gives a conditional mean of
`(Lambda B' + Gamma) x`, so the likelihood constrains `B` and `Gamma` only
through their sum: for any `Delta`, the pair `(B + Delta, Gamma - Lambda
Delta')` fits identically. `B` was therefore not identified by the data, and
the split rested entirely on the ARD prior.

* `vb_ard()` gains `gamma_anchor`: indicators whose direct effects are held at
  exactly zero. With at least `K` anchors of linearly independent loadings,
  `B` is identified.
* `select_gamma_anchor()` picks anchors from a first-pass fit — strong loading
  SNR, then weakest direct-effect evidence, then a linear-independence check.
  The second criterion matters: ranking on loading SNR alone can anchor an
  indicator that genuinely has a direct effect, which biases `B`.
* `biomimic()` gains `gamma_anchor`, defaulting to `"auto"` (two-pass
  anchoring). Pass `NULL` for the previous unanchored behaviour.
* `top_k()` forces anchors into the panel and puts them first, so Stage 3 fixes
  its marker loading on an indicator already constrained to be pure.
* `check_identification()` now tests the rank condition (at least `K` pure
  indicators) rather than a fixed `|S| <= p - 2`.

Anchoring fixes identification. It does **not** make the direct-effect SNR
pivotal: under a null with an active latent covariate effect the statistic
still grows with `n` and `p`. A bias/variance decomposition (200 replicates
per cell) locates the cause: it is a **leakage bias**, not variance
miscalibration. The variational SD is well calibrated (slightly conservative;
empirical-to-variational SD ratio 0.85-1.1). What happens at large `p` is that
the fitted latent pathway under-absorbs the covariate effect -- `Lambda` comes
out inflated and `B` deflated, their product roughly preserved -- and the
unabsorbed residue lands in the direct effects of loaded indicators,
proportionally to `B` and to `lambda_j`. With `B = 0` every cell is clean and
flat in `n`; noise indicators are unbiased at every `B`; anchored and
unanchored fits show the same bias.

Consequently, permutation of the covariate cannot calibrate this statistic:
permuting destroys `B` and with it the bias the null needs to contain.
Thresholds for `SNR_Gamma` should come from a parametric bootstrap under
`gamma = 0` with `B` held at its fitted value, per dataset and per covariate
column.

## Multi-factor guidance reversed

* The advice to set `use_ard = FALSE` for K > 1 is withdrawn and reversed:
  ARD should stay ON for every K. A dedicated comparison (150 replicates per
  cell, K = 2) showed ARD-off was the main driver of the multi-factor
  direct-effect pathology (null false-positive rates at the default threshold
  several-fold higher), while lambda-screening F1 is 1.000 under either
  setting. The simulation and application scripts now use ARD on throughout.

## Calibrated direct-effect threshold

* New `calibrate_gamma_threshold()`: parametric-bootstrap calibration of the
  direct-effect SNR threshold. Simulates from the fitted model under the
  global null (`Gamma = 0`, `B`/`Lambda`/`Psi` at their fitted values, real
  `X`), refits, and returns the level-quantile of the per-replicate maximum
  SNR over the panel — a family-wise 5% threshold whose null *contains* the
  leakage bias described above. Refits run on a PSOCK cluster
  (`BIOMIMIC_NCORES` honoured); results are seed-reproducible and identical
  serial or parallel.
* `top_k()` accepts `gamma_threshold = "bootstrap"` and stores the
  calibration in the selection (`$gamma_threshold`, `$gamma_calibration`).
* `biomimic()` now exposes `gamma_threshold` and passes it through to
  `top_k()`. Previously the wrapper hard-wired the fixed 3.0 default and
  rejected any attempt to override it (the argument fell through `...` to
  `vb_ard()` and errored).
* `parallel` added to Imports.

## Fixes

* Removed a dead `if (use_ard) / else` branch in the loading update whose two
  arms computed the same expression, and the parameter it left unused.

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
