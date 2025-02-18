library(devtools)
library(metafor)
library(tidyverse)


load_all()

# Let's see if by changhin the model data, the VCV should be computed again or not

VCV <- vcalc(lnrr_vi,
             cluster = fish$paper_ID,
             obs     = fish$es_ID,
             data    = fish,
             rho     = 0.5)

res <- rma.mv(lnrr, VCV,
              random = ~ 1 | paper_ID,
              data = fish)

res2 <- rma.mv(lnrr, lnrr_vi,
              random = ~ 1 | paper_ID,
              data = fish)



res

res2

# Check how to output with the data used so it can be used for orchard plot?
