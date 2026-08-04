## Submission

This is a new submission. 'biomimic' implements a three-stage pipeline for
high-dimensional latent-variable modelling: variational Bayes with automatic
relevance determination screens thousands of candidate variables, a top-k
panel is selected, and structural equation models (MIMIC, multi-group CFA,
and MIMIC with direct effects) are fitted with 'lavaan'. The version number
(0.2.0) reflects pre-CRAN development; this is the first CRAN submission.

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
  corpcor) are used conditionally, guarded by `requireNamespace()`.

## Downstream dependencies

There are no reverse dependencies (new package).
