## Submission

This is a new submission. 'biomimic' takes high-dimensional omics data to
structural equation models fitted with 'lavaan': features are ranked
('limma', 'cate' or the package's own variational Bayes screening), grouped
into modules ('WGCNA', curated sets or the screening's factors), checked
against each sample's global level, and modelled (multi-factor MIMIC,
measurement invariance, decomposition, mediation, bifactor and
random-intercept cross-lagged panel models). The version number (0.4.0)
reflects pre-CRAN development; this is the first CRAN submission.

## Test environments

* Local: Windows 10, R 4.4.1

<!-- Maintainer: before the actual submission, please also run and paste the
     results of win-builder (release + devel) and R-hub, e.g.
       devtools::check_win_release(); devtools::check_win_devel()
       rhub::rhub_check()
     CRAN expects at least one clean check on a current R-devel. -->

## R CMD check results

0 errors | 0 warnings | 1 note

```
* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Alexandre Henrique Carvalho Marques <alexandre.marques.088@ufrn.edu.br>'
  New submission
```

The only NOTE is the standard "New submission" message. The URL and
BugReports fields point to the package's public GitHub repository
(https://github.com/ahcm088/biomimic).

## Notes for the reviewer

* Examples that fit a model are wrapped in `\donttest{}` to keep per-example
  run time short; they are exercised by the package's `testthat` suite.
* All packages listed under Suggests (logger, scales, igraph, ggraph,
  corpcor, WGCNA, limma, cate) are used conditionally, guarded by
  `requireNamespace()`; the vignette chunks that use WGCNA and limma run
  only when they are installed.

<!-- Maintainer: 'cate' is archived on CRAN, so the incoming check will note
     a Suggests package not in a mainstream repository. Decide before
     submission: keep it and explain here, or reimplement the part of
     cate::cate() that rank_cate() uses. The local devtools::check() skips
     the incoming step; rerun with remote = TRUE, incoming = TRUE. -->

## Downstream dependencies

There are no reverse dependencies (new package).
