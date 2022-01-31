# test
rm(list=ls())
install.packages("devtools")
install.packages("tidyverse")
install.packages("metafor")
install.packages("patchwork")
install.packages("R.rsp")

devtools::install_github("daniel1noble/orchaRd", force = TRUE, build_vignettes = TRUE)
#remotes::install_github("rvlenth/emmeans", dependencies = TRUE, build_opts = "")

library(orchaRd)
library(metafor)
library(emmeans)
library(tidyverse)

# Data
data(fish)
warm_dat <- fish

# The Model
model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 | es_ID), mods = ~ experimental_design + trait.type + deg_dif + treat_end_days, method = "REML", test = "t", data = warm_dat,                               control=list(optimizer="optim", optmethod="Nelder-Mead"))

# Example 1: Overall marginal means for each level of experimental design
  # Two step with marginal means
  overall <- marginal_means(model, mod = "experimental_design", group = "group_ID",
                            data = warm_dat)
  orchard_plot(overall, xlab = "lnRR", trunk.size = 2, branch.size = 2, twig.size = 0.5,
               angle = 45)

  # Directly with model
  orchard_plot(model, xlab = "lnRR", data = warm_dat, mod = "experimental_design", group = "group_ID", trunk.size = 2, branch.size = 2, twig.size = 0.5, angle = 45, marginal = TRUE)

# Example 2: Overall marginal mean across all designs
  # Two step with marginal means
  overall1.1 <- marginal_means(model, group = "group_ID", data = warm_dat)
  orchard_plot(overall1.1, xlab = "lnRR", trunk.size = 2, branch.size = 1.2, twig.size = 2)

  # Directly with model
  orchard_plot(model, group = "group_ID", data = warm_dat, xlab = "lnRR", trunk.size = 2, branch.size = 1.2, twig.size = 2, marginal = TRUE)

# Example 3: Marginalised overall mean for each temp
  # Two step with marginal means
  overall2 <- marginal_means(model, group = "group_ID", mod = "1", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif",  weights = "prop", data = warm_dat)
  orchard_plot(overall2, xlab = "lnRR", group = "group_ID", condition.lab = "Temparature")

  # Directly with model
  orchard_plot(model, xlab = "lnRR", group = "group_ID", condition.lab = "Temparature", mod = "1", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif",  weights = "prop", data = warm_dat, marginal = TRUE)


# Example 4: Marginalised means within each trait category as opposed to design
# Two step with marginal means
across_trait <- marginal_means(model,  mod = "trait.type", group = "group_ID", data = warm_dat)
orchard_plot(across_trait, xlab = "lnRR")

# With model directly
orchard_plot(model, xlab = "lnRR", mod = "trait.type", group = "group_ID", data = warm_dat, marginal = TRUE)


# Example 5: Marginalised means for each trait category by different temperature differences
  # Two step with marginal means
  across_trait_by_degree_diff <- marginal_means(model, mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", group = "group_ID", data = warm_dat)
  orchard_plot(across_trait_by_degree_diff,  xlab = "lnRR",  data = warm_dat)

  # With model directly
  orchard_plot(model, mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", group = "group_ID", data = warm_dat, xlab = "lnRR", marginal = TRUE)



  # Example 6: Marginalised means for each trait category by different temperature differences at treatment end days held at 10
  # Two step with marginal means
across_trait_by_degree_diff_at_treat_end_days10 <- marginal_means(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15),  treat_end_days = 10), group = "group_ID", by = "deg_dif")
orchard_plot(across_trait_by_degree_diff_at_treat_end_days10, xlab = "lnRR")

# With model directly
orchard_plot(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15),  treat_end_days = 10), group = "group_ID", by = "deg_dif", xlab = "lnRR", marginal = TRUE)


# Example 7: Marginalised means for each trait category by different temperature differences averaging only treatment end days 10 & 50
# Two step with marginal means
across_trait_by_degree_diff_at_treat_end_days10And50 <- marginal_means(model, data = warm_dat, at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), mod = "trait.type", by = "deg_dif", group = "group_ID")
orchard_plot(across_trait_by_degree_diff_at_treat_end_days10And50, xlab = "lnRR")

# With model directly
orchard_plot(model, data = warm_dat, at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), mod = "trait.type", by = "deg_dif", group = "group_ID", xlab = "lnRR", marginal = TRUE)

# Example 8: Marginalised means for each trait category by different temperature differences averaging only treatment end days 10 & 50 with equal weights
# Two step with marginal means
across_trait_by_treat_end_days10And50 <- marginal_means(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days", weights = "equal", group = "group_ID")
orchard_plot(across_trait_by_treat_end_days10And50, xlab = "lnRR")

# With model directly
orchard_plot(model, data = warm_dat, mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)), by = "treat_end_days", weights = "equal", group = "group_ID", xlab = "lnRR", marginal = TRUE)

# Example 9: Marginalised means for each trait category assuming heteroscedastic error within each
# Model
model_het <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
                             random = list(~1 | group_ID, ~1 + trait.type| es_ID),
                             mods = ~ trait.type + deg_dif, method = "REML",
                             test = "t", rho = 0, struc = "HCS", data = warm_dat,
                             control=list(optimizer="optim", optmethod="Nelder-Mead"))

# Two step process
HetModel <- marginal_means(model_het,
                           group = "group_ID",
                           mod = "trait.type",
                           at = list(deg_dif = c(5, 10, 15)),
                           by = "deg_dif", weights = "prop", data = warm_dat)
orchard_plot(HetModel, xlab = "lnRR")

# With model directly
orchard_plot(model_het, group = "group_ID",
             mod = "trait.type",
             at = list(deg_dif = c(5, 10, 15)),
             by = "deg_dif", weights = "prop", data = warm_dat,
             xlab = "lnRR", marginal = TRUE)
