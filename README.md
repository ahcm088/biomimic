# biomimic

<!-- badges: start -->
[![R-CMD-check](https://github.com/ahcm088/biomimic/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ahcm088/biomimic/actions/workflows/R-CMD-check.yaml)
[![License: MIT](https://img.shields.io/badge/license-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R >= 4.1](https://img.shields.io/badge/R-%3E%3D%204.1-blue.svg)](https://cran.r-project.org/)
<!-- badges: end -->

Structural equation models for omics data. `biomimic` takes features through
ranking, module detection and a check of what each module measures, and then
fits the chosen panels as models in [`lavaan`](https://lavaan.ugent.be), with
p-values for every path, likelihood-ratio tests, fit indices and residual
diagnostics.

Omics analyses usually stop at a list: genes ranked by a test, modules of
co-expressed genes. Structural models ask questions that a list cannot
answer. Is a group difference carried by one programme or by two? Does a gene
differ on its own, or only because its programme does? Does a clinical score
follow the programme, and through which path? Does the panel measure the
programme the same way in patients and controls? The obstacle is practical:
SEM software takes the indicators as given, and choosing 15 out of 20,000 is
not something it does. `biomimic` handles that part and leaves the modelling
to `lavaan`.

## Installation

`biomimic` is not on CRAN. Install the development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("ahcm088/biomimic", build_vignettes = TRUE)
```

It requires R >= 4.1.0 and `lavaan` (>= 0.6-17). The front ends use `limma`
and `WGCNA` when installed; `rank_cate()` needs `cate`, which is archived on
CRAN and installs from the
[archive](https://cran.r-project.org/src/contrib/Archive/cate/).

## The pipeline

| Step | Functions | Front ends |
|---|---|---|
| 1. Rank features | `rank_limma()`, `rank_cate()`, `vb_ard()` + `rank_variables()` | limma, CATE, VB-ARD |
| 2. Find modules | `modules_wgcna()`, `modules_from_sets()`, `modules_ard()`, `annotate_modules()` | WGCNA, curated gene sets, VB-ARD factors |
| 3. Check what a module measures | `global_level_check()`, `external_check()`, `normalise_arrays()`, `coordination_diagnostic()` | |
| 4. Model | `select_panels()` or `as_selection()`, then `fit_sem()`, `lrt_structural()`, `test_invariance()`, `decompose_effects()`, `fit_mediation()`, `fit_bifactor()`, `within_between()`, `fit_riclpm()` | lavaan |

VB-ARD is the factor model that `biomimic` introduced: it estimates factors
and each feature's direct effect jointly. It is one front end among several,
not a replacement for the others. In SLE blood it matched CATE as an
adjustment, and found the interferon module in two of four cohorts against
four of four for WGCNA.

## Example

```r
library(biomimic)

# Y: samples x genes; group: 0/1
covs <- data.frame(group = group)
mods <- modules_wgcna(Y)
ann  <- annotate_modules(mods, list(interferon = isg, lymphocyte = lym_markers))

# is the interferon panel just each sample's global level?
global_level_check(Y, head(mods$modules[[ann$module[1]]]$genes, 15), group)

# two factors in one model, with likelihood-ratio tests
sel <- select_panels(mods, ann$module, k = 15, Y = Y, covariates = covs)
fit <- fit_sem(sel, type = "multi_mimic", factor_labels = c("Interferon", "Lymphocyte"))
summary(fit)
lrt_structural(fit)

# do patients and controls share the measurement model, and if not, why not?
panels <- lapply(sel$panels, function(p) unname(sel$labels[p]))
inv <- test_invariance(as_selection(Y, panels[[1]], covs), group_var = "group")
inv$worst$intercepts

# does each gene differ through the programmes or on its own?
decompose_effects(scale(Y), panels, targets = rank_limma(Y, covs)$variable[1:200],
                  group = group)
```

`vignette("biomimic")` runs every step on simulated data with known
structure, including mediation and bifactor models. Each model keeps the
`lavaan` object (`$lavaan_fit`) for anything the wrappers do not cover.

## Data

`neuro_antibodyome` is a 400-feature subset (137 sera, 46 AD cases and 91
controls) of the ProtoArray data of DeMarshall et al. (2015), GEO
[GSE62283](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62283). It
is kept as a cautionary example: cases and controls come from disjoint sample
collections, the panel follows each serum's whole-array level, and the AD
effect disappears after per-array normalisation (`?neuro_antibodyome`).

## Citation

```r
citation("biomimic")
```

A manuscript describing the package is in preparation. Until it appears,
please cite the package.

## License

MIT © the package authors. See [LICENSE.md](LICENSE.md).
