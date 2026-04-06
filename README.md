
<!-- README.md is generated from README.Rmd. Please edit that file -->

# **Spacelink**

<!-- badges: start -->
<!-- badges: end -->

`spacelink` is a unified statistical framework for detecting and
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

You can install the development version of spacelink from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("hanbyul-lee/spacelink")
```

## Documentation

| Page | Description |
| ----------------- | ------------------------------ |
| [Installation]() | Setup |
| [Overview]() | Descriptions of main functions |
| [Spacelink Workflow]() | Examples and guides for Spacelink analysis at the global tissue level |
| [Spacelink (ct-SVG) Workflow]() | Examples and guides for Spacelink (ct-SVG) analysis at the cell-type-specific level |
| [Runtime & Memory Usage]() | Benchmarks of computational time and memory requirements across different dataset sizes |
| [Illustration on CosMx Data]() | Application of Spacelink on a large-scale single-cell resolution dataset |
| [Illustration on Visium Data]() | Application of Spacelink on a medium-scale spot resolution dataset |
| [Disease Informativeness Evaluation]() | Evaluation of ESV disease informativeness using PoPS (Polygenic Priority Score) |
| [SVG Ranking Metric Comparison]() | Comparison of ESV with other metrics for ranking genes by spatial variability |
