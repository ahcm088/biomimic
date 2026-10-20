# biomimic 0.4.0

The package is reorganised as a pipeline from omics features to structural
equation models, with the VB-ARD screening as one front end among several.
The change follows the validation in SLE blood: WGCNA found the interferon
module in four of four cohorts and the VB-ARD factors in two, while the
results that replicated came from the structural models, whichever front end
supplied the panels.

## Modules and selections

* `modules_wgcna()` (signed WGCNA network, soft power by scale-free fit,
  genes ranked by module membership), `modules_from_sets()` (curated gene
  sets) and `modules_ard()` (the factors of a `vb_ard()` fit) return modules
  in one format. `annotate_modules()` names them by known markers.
* `as_selection()` builds a selection from any named list of panels, and
  `select_panels()` from the top genes of chosen modules, so every Stage-3
  function now accepts panels chosen outside the VB-ARD screening. Feature
  names are made safe for `lavaan` and the originals kept as labels.

## Rankings

* `rank_limma()` and `rank_cate()` (CATE, adjusting each feature for hidden
  factors) return the columns of `rank_variables()`.

## Diagnostics

* `global_level_check()` correlates a panel with each sample's whole-array
  mean and re-estimates the group effect adjusted for it; `external_check()`
  correlates panel scores with external measurements; `normalise_arrays()`
  applies per-array median or quantile normalisation; `panel_scores()`
  gives standardised panel means.
* `coordination_diagnostic()` also reports the correlation of the panel
  with the global level, and its print no longer reads a high ratio as a
  coordinated block.

## Structural models

* `fit_bifactor()`: general factor on all indicators, orthogonal specific
  factors on named panels, covariate paths with Wald and likelihood-ratio
  tests, and the share of each panel's common variance on the general factor.
* `fit_mediation()`: latent (panel) or observed mediators, delta-method
  indirect effects and a likelihood-ratio test per mediator.
* `decompose_effects()`: each feature's group effect split into parts
  through each factor and an own part, in parallel.
* `within_between()` and `fit_sem(cluster = )` for repeated measures;
  `fit_riclpm()` fits the random-intercept cross-lagged panel model of
  Hamaker, Kuiper and Grasman (2015).
* `test_invariance()` returns `worst`, the loadings and intercepts whose
  equality constraints the data reject most (univariate score tests). A
  level that does not converge no longer stops the function: its indices
  are NA (`fit_indices$converged`), its verdict "Could not compute", and
  every comparison keeps its row, so verdicts can no longer shift against
  the levels when one likelihood-ratio test fails.
* `lrt_structural()` adds a `label` column with the factor names of the fit.
* `network_analysis()` works with igraph 2.0 (`edge_density()` in place of
  the deprecated `graph.density()`; the partial-correlation matrix is made
  exactly symmetric, as igraph 1.6 and later require).

## Documentation

* New vignette, on simulated data with known structure, that runs the whole
  pipeline.
* `neuro_antibodyome` is now documented as a cautionary example: its AD
  contrast compares disjoint sample collections and disappears after
  per-array normalisation.

# biomimic 0.3.5

* `fit_sem()` now checks the orientation of each factor after fitting. lavaan
  scales a factor by fixing its first (marker) loading to +1, so when the
  marker anti-correlates with most of the panel, every other loading comes out
  negative and the sign of `B` flips: "lower in patients" or "higher in
  patients" then depends on which indicator was listed first. The fit is
  unchanged by the reflection; only the reading of the sign is. Found in the
  SLE validation, where three of ten down-sampled panels returned a positive
  `B` that was the same statement as the other seven's negative one. With the
  new `orient = TRUE` (default), the marker of a misaligned factor is replaced
  by its strongest majority-sign indicator and the model refitted; for
  `type = "direct"` the replacement is drawn only from indicators without a
  direct effect. The result carries an `orientation` table. Fits whose marker
  already sided with the majority --- the usual case --- are unchanged.
* `select_gamma_anchor()` restricts, for `K = 1`, its candidates to indicators
  with the majority loading sign. The anchor goes first in the Stage-2 panel
  and becomes the Stage-3 marker, so this removes the problem at its source
  for panels built by `top_k()`; the check in `fit_sem()` covers panels built
  any other way.
* New `lrt_structural()`: likelihood-ratio test of each structural
  coefficient, by refitting with it fixed at zero, reported next to the Wald
  test. The Wald z of a marker-scaled `B` depends on the marker (Gonzalez and
  Griffin 2001). On the SLE transfer panel, all loadings but two one sign, it
  moved only between p = 0.035 and 0.046 across the 15 possible markers; on a
  methylation panel whose loadings split 9/6 it moved between p = 0.00002 and
  0.48. The likelihood-ratio test gave one value for every marker in both.
* Known, not changed: `fit_sem()`'s default `se = "robust"` makes lavaan
  replace MLR's Huber-White standard errors with `robust.sem` (expected
  information), although the documentation says Huber-White. Estimates are
  unaffected and standard errors move in the third decimal, but the scaled
  likelihood-ratio test depends on the information type. Left as is pending a
  decision, because changing it shifts every standard error already reported.
* `coordination_diagnostic()` reports an empirical p-value per channel, the
  share of random sets at least as co-ordinated as the panel. It is a fair
  test for the gamma channel only: lambda selects on co-ordination and beats
  random sets almost by construction. The documentation and `print()` now
  also say what a high ratio does not mean --- that the panel carries the
  design effect. In a methylation check a cell-composition axis owned the
  lambda channel at ratio 5.5--9.5 while the age signal sat in gamma.

# biomimic 0.3.4

* New `coordination_diagnostic()`, with a `print()` method, reporting for each
  screening channel how strongly the top-`k` selected variables still co-vary
  once the design is accounted for, against the background correlation of
  random variable sets of the same size in the same data. Stage 1 will always
  return a factor, and ranking by loading SNR will always return the variables
  that load on whatever the dominant shared mode happens to be --- biology in
  one dataset, assay structure in another. This separates the two cases before
  a MIMIC is interpreted. Correlations are computed on residuals after
  projecting out `X` by default; `reference` restricts them to a subset of
  samples instead, typically the controls.
* No threshold is imposed and no verdict is returned, deliberately: the ratio
  is a property of the assay as much as of the biology. The useful comparison
  is between channels within one dataset, and the `print()` method gives three
  measured reference points for scale.

# biomimic 0.3.3

* PSOCK workers start with their own default `.libPaths()` and do not inherit
  the master's, so they could not see the temporary library that
  `R CMD build`/`check` installs into and would fail to find biomimic on a
  machine where it is not already installed. The workers are now pointed at
  the master's libraries, so they load the very copy the calling session is
  running. Applies to both `.biomimic_parallel()` and the cluster in
  `run_simulation()`. biomimic's own vignette does not exercise a parallel
  path, so this never broke its build --- unlike the companion package, where
  it did --- but the same latent fault was present.

# biomimic 0.3.2

Parallel-execution robustness; no change to any statistical result.

* Parallel runs now verify that the workers load the same biomimic as the
  calling session, and abort on a version mismatch. PSOCK workers resolve
  the package's internal functions through their own installed copy, so a
  stale install silently ran old code for every replicate; only a missing
  object failed loudly, while a numerical difference would have quietly
  contaminated bootstrap calibration and simulation output. Same-version/
  different-build cases warn once per session, including
  `devtools::load_all()` sessions, whose edits never reach the workers. Set
  `BIOMIMIC_SKIP_VERSION_CHECK` to bypass.
* The worker count is capped at 2 under `R CMD check --as-cran`, which sets
  `_R_CHECK_LIMIT_CORES_` and makes `makeCluster()` error above that.
* `run_simulation()` now reuses `.biomimic_n_cores()` rather than
  duplicating the core-count logic, so both paths follow one policy.

# biomimic 0.3.1

* `run_simulation()` now fits every replicate with `use_ard = TRUE`,
  matching the protocol of the reported simulation study; it previously
  disabled ARD for `K > 1`, the configuration this same study showed to
  be harmful.
* NEWS correction: the 0.3.0 rank-condition check runs in the pre-fit
  identification screen used by `build_lavaan_model()`/`fit_sem()`, not
  in the exported `check_identification()` as previously stated.
* README: the quick-start loadings excerpt now shows consecutive output
  lines (the MGC18216 row was silently skipped).
* Vignette and NEWS wording: `test_invariance()` fits configural, metric,
  scalar and strict models (strict was previously omitted from the
  descriptions), and the direct-effect extension is documented as on by
  default rather than opt-in.

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
* The pre-fit identification screen run by
  `build_lavaan_model()`/`fit_sem()` now tests the rank condition (at
  least `K` pure indicators) rather than a fixed `|S| <= p - 2`.

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
* `test_invariance()` for configural, metric, scalar and strict measurement invariance.
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
