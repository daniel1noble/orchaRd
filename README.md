## Introducing the Orchard Plot for Meta-analysis

# Citing orchaRd

To cite `orchaRd` in publications one can use the following reference:

Nakagawa, S., Lagisz, M., O'Dea, R. E., Rutkowska, J., Yang, Y., Noble, D. W., & Senior, A. M. (2019). The Orchard Plot: Cultivating Forest Plots for Use in Ecology, Evolution and Beyond. *EcoEvoRxiv* https://doi.org/10.32942/osf.io/epqa7 (submitted to *Research Synthesis Methods*)

# Installation

To install `orchaRd` use the following code in R:

```
install.packages("devtools")
install.packages("tidyverse")
install.packages("metafor")
install.packages("patchwork")
install.packages("R.rsp")

devtools::install_github("daniel1noble/orchard_plot", force = TRUE, build_vignettes = TRUE)

library(orchaRd)
library(patchwork)
library(tidyverse)
library(metafor)
```
