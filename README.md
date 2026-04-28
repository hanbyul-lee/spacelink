
<!-- README.md is generated from README.Rmd. Please edit that file -->

# **Spacelink**

<!-- badges: start -->
<!-- badges: end -->

**Spacelink** is a unified statistical framework for detecting and
prioritizing SVGs at both global tissue and cell-type resolution.
Spacelink employs an adaptive multi-kernel model to capture spatial
variance across diverse length scales, and its cell-type specific
version introduces a data-driven gating strategy to correct for spatial
colocalization, designed to improve the specificity for cell types that
are weakly represented in mixed spots relative to more abundant
colocalizing cell types. To summarize spatial variability, we define
**Effective Spatial Variability (ESV)**, a metric which integrates variance
magnitude of each component kernel and its corresponding spatial scale
into a single interpretable score directly suited for genetic analyses.

![](man/figures/spacelink_schematic.png)

## Installation

System requirements:
- R ≥ 4.0.0
- Supported: Windows, macOS, Linux

You can install the development version of spacelink from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hanbyul-lee/spacelink")
```

The following dependency package may need to be installed manually from CRAN:

``` r
install.packages("gaston")
```

The other dependency packages are as follows. If automatic installation fails, you can manually install them from CRAN:

```r
install.packages(c("Rcpp", "RcppML", "pracma", "RANN", "fields", "Matrix"))
```


## Documentation

| Page | Description |
| ------------------------ | ------------------------------ |
| [Installation](https://hanbyul-lee.github.io/spacelink/articles/installation.html) | Setup |
| [Quick Start](https://hanbyul-lee.github.io/spacelink/articles/quick_start.html) | Get started with a minimal working example |
| [Overview](https://hanbyul-lee.github.io/spacelink/articles/overview.html) | Descriptions of main functions |
| [Spacelink Global SVG Vignette](https://hanbyul-lee.github.io/spacelink/articles/spacelink_vignette.html) | Examples and guides for Spacelink analysis at the global tissue level |
| [Spacelink ct-SVG Vignette](https://hanbyul-lee.github.io/spacelink/articles/spacelink_ctSVG_vignette.html) | Examples and guides for Spacelink (ct-SVG) analysis at the cell-type-specific level |
| [Illustration on CosMx Data](https://hanbyul-lee.github.io/spacelink/articles/illustration_on_cosmx_data.html) | Application of Spacelink on a large-scale single-cell resolution dataset |
| [Illustration on Visium Data](https://hanbyul-lee.github.io/spacelink/articles/illustration_on_visium_data.html) | Application of Spacelink on a medium-scale spot resolution dataset |
| [Disease Informativeness Evaluation](https://hanbyul-lee.github.io/spacelink/articles/disease_informativeness_evaluation.html) | Evaluation of ESV disease informativeness using PoPS (Polygenic Priority Score) |
| [Lengthscale Estimation](https://hanbyul-lee.github.io/spacelink/articles/lengthscale_estimation.html) | Comparison of length-scale estimation performance |
| [Runtime & Memory Usage](https://hanbyul-lee.github.io/spacelink/articles/runtime_memory_usage.html) | Benchmarks of computational time and memory requirements across different dataset sizes |
