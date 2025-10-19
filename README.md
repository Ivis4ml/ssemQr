## Introduction

`ssemQr` is a package that ultilizes the Proximal Alternating Linearized Maximal to solve the
non-convex non-smooth sparse structrual equation model for identification trans-eQTLs via joint eQTL mapping and
inferring Gene Regulatory Network(GRN).

## Installation

To install `ssemQr`, you need a C++ compiler such as `g++` or `clang++` with C++11 feature,
and for Windows users, the [Rtools](https://cran.r-project.org/bin/windows/Rtools/index.html)
software is needed (unless you can configure the toolchain by yourself).

The installation follows the typical way of R packages on Github:

  ```r
library(devtools)
install_github("Ivis4ml/ssemQr")
```

`ssemQr` package will be uploaded on CRAN later. Then you can install it via CRAN

```r
install.packages("ssemQr")
```

## Vignette
[ssemQr-introduction](https://github.com/Ivis4ml/ssemQr/blob/master/inst/doc/ssemQr.pdf)


## Citation
Xin Zhou, Xiaodong Cai, Joint eQTL mapping and inference of gene regulatory network improves power of detecting both cis- and trans-eQTLs, Bioinformatics, Volume 38, Issue 1, 1 January 2022, Pages 149–156, https://doi.org/10.1093/bioinformatics/btab609


## Acknowledgement
Special thanks must go to Yilun Zhang <3100105044@zju.edu.cn> for his contribution to the numerical optimization for the code and mathematical parts of the research.

Special thanks must go to David López Carbonell <davidlc@unizar.es> for his contribution to the his PR[https://github.com/Ivis4ml/ssemQr/pull/3], which is a comprehensive refactoring of the C++ and R interface for L2 regression functions has significantly improved the ssemQr package by upgrading to double-precision arithmetic for better numerical accuracy, modernizing the Rcpp integration with native exports to replace outdated .Call() wrappers, and enhancing overall code quality and maintainability. These improvements not only make the codebase cleaner and more aligned with current best practices, but also ensure the package is more numerically stable and easier to maintain going forward