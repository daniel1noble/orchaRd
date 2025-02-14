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

# Both models had different results.

head(res$vi)

head(res2$vi)

identical(res$vi, res2$vi)

res$vi - res2$vi

# And different sampling variances


# Let's check if we can update the data from a model that uses VCV

# 1: Meta analysis removing the first paper

fish_1 <- fish %>% filter(paper_ID != "p047")

VCV1 <- vcalc(lnrr_vi,
             cluster = paper_ID,
             obs     = es_ID,
             data    = fish_1,
             rho     = 0.5)

res1 <- rma.mv(lnrr, VCV1,
              random = ~ 1 | paper_ID,
              data = fish_1)

mod_results(res1, group = "paper_ID")


# 2: Meta analysis but updating the model

# --------------------------------

# test!

vcv_args <- list(vi      = "lnrr_vi",
                 cluster = "paper_ID",
                 obs     = "es_ID",
                 rho     = 0.5)



loo <- leave_one_out(res, group = "paper_ID", vcalc = vcv_args)

loo

loo |>
  filter(left_out == "p047") 
