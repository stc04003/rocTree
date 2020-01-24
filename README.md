
**rocTree**
-----------

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active) [![minimal R version](https://img.shields.io/badge/R%3E%3D-3.4.0-6666ff.svg)](https://cran.r-project.org/) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/rocTree)](https://cran.r-project.org/package=rocTree) [![packageversion](https://img.shields.io/badge/Package%20version-1.1.0-orange.svg?style=flat-square)](commits/master) [![Travis-CI Build Status](https://travis-ci.org/stc04003/rocTree.svg?branch=master)](https://travis-ci.org/stc04003/rocTree) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/stc04003/rocTree?branch=master&svg=true)](https://ci.appveyor.com/project/stc04003/rocTree) [![Last-changedate](https://img.shields.io/badge/last%20change-2020--01--24-yellowgreen.svg)](/commits/master)

<!-- README.md is generated from README.Rmd. Please edit that file -->
### ROC-guided survival trees and ensembles

------------------------------------------------------------------------

### Development

The package is under active development.

### Installation

You can install `rocTree` from **GitHub** with:

``` r
## install.packages("devtools")
devtools::install_github("stc04003/rocTree")
```

### Description

The `rocTree` provides implementations to a unified framework for tree-structured analysis with censored survival outcomes. Different from many existing tree building algorithms, the `rocTree` package incorporate time-dependent covariates by constructing a time-invariant partition scheme on the survivor population. The partition-based risk prediction function is constructed using an algorithm guided by the Receiver Operating Characteristic (ROC) curve. Specifically, the generalized time-dependent ROC curves for survival trees show that the target hazard function yields the highest ROC curve. The optimality of the target hazard function motivates us to use a weighted average of the time-dependent area under the curve (AUC) on a set of time points to evaluate the prediction performance of survival trees and to guide splitting and pruning. Moreover, the `rocTree` package also offers a novel ensemble algorithm, where the ensemble is on unbiased martingale estimating equations.

### Online documentations

[Online document](https://www.sychiou.com/rocTree/index.html) includes:

-   Package vignette on [simulating data used in examples](https://www.sychiou.com/rocTree/articles/rocTree-sim.html).
-   Package vignette on [growing time-invariant survival trees](https://www.sychiou.com/rocTree/articles/rocTree-tree.html).
-   Package vignette on [ensemble method](https://www.sychiou.com/rocTree/articles/rocTree-ensemble.html).

Reference
---------

Yifei Sun, Sy Han Chiou, Mei-Cheng Wang. ROC-Guided Survival Trees and Ensembles, (2019+). [doi: 10.1111/biom.13213](https://www.ncbi.nlm.nih.gov/pubmed/31880315).

Disclaimer
----------

The `rocTree` package does not implement the works proposed by Drs. Hossain, Hassan, and Bailey (reference below), though they share similar names.

Hossain, MM; Hassan, MR; Bailey, J, ROC-tree: A novel decision tree induction algorithm based on receiver operating characteristics to classify gene expression data, (2008), 130, 2008, 2 pp. 455--465
