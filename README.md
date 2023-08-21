# Introducing the Orchard Plot for Meta-analysis
[![R-CMD-check](https://github.com/daniel1noble/orchaRd/actions/workflows/check-standard.yaml/badge.svg)](https://github.com/daniel1noble/orchaRd/actions/workflows/check-standard.yaml)
[![test-coverage](https://github.com/daniel1noble/orchaRd/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/daniel1noble/orchaRd/actions/workflows/test-coverage.yaml)
[![codecov](https://codecov.io/gh/daniel1noble/orchaRd/branch/main/graph/badge.svg?token=KqQLvcGfLv)](https://codecov.io/gh/daniel1noble/orchaRd)
[![Ask Us Anything\ !](https://img.shields.io/badge/Ask%20me-anything-1abc9c.svg)](https://github.com/daniel1noble/orchaRd/issues/new)
![Open Source Love](https://badges.frapsoft.com/os/v2/open-source.svg?v=103)
[![DOI](https://zenodo.org/badge/358763476.svg)](https://zenodo.org/badge/latestdoi/358763476)
[![Github Downloads](https://img.shields.io/github/downloads/daniel1noble/orchaRd/total.svg)]()

# Citing orchaRd

To cite `orchaRd` 2.0 in publications one can use the following reference:

Shinichi Nakagawa, Malgorzata Lagisz, Rose E. O'Dea, Patrice Pottier, Joanna Rutkowska, Alistair M. Senior, Yefeng Yang, Daniel W.A. Noble. 2023. orchaRd 2.0: An R package for visualizing meta-analyses with orchard plots. **Methods in Ecology and Evolution**, https://doi.org/10.1111/2041-210X.14152
(preprint = **EcoEvoRxiv**, https://doi.org/10.32942/X2QC7). 

For earlier versions please cite:

Nakagawa, S., Lagisz, M., O'Dea, R. E., Rutkowska, J., Yang, Y., Noble, D. W., & Senior, A. M. (2019). The Orchard Plot: Cultivating Forest Plots for Use in Ecology, Evolution and Beyond. *Research Synthesis Methods* https://doi.org/10.1002/jrsm.1424 12: 4-12 (preprint = *EcoEvoRxiv* https://doi.org/10.32942/osf.io/epqa7)

# Installation
To install `orchaRd` use the following code in R:

```
install.packages("pacman")
pacman::p_load(devtools, tidyverse, metafor, patchwork, R.rsp, emmeans)

devtools::install_github("daniel1noble/orchaRd", force = TRUE)
library(orchaRd)
```

# How to use?
We detail how to use the `orchaRd` package in the [vignette](https://daniel1noble.github.io/orchaRd/). 

# Issues with orchaRd 2.0?
Please note that orchaRd 2.0 is still under active development and testing. If you use it, you should check that the results are what you expect. We do have a number of tests already in place, but there may still be situations where it fails. If you find a bug or a situation that doesn't match your expectations let us know by lodging an [issue](https://github.com/daniel1noble/orchaRd/issues) on GitHub. 
