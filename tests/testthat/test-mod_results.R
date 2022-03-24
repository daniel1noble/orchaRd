context("Checking mod_results function..")

# Simple eklof data
data(eklof)
eklof<-metafor::escalc(measure="ROM", n1i=N_control, sd1i=SD_control,
m1i=mean_control, n2i=N_treatment, sd2i=SD_treatment, m2i=mean_treatment,
data=eklof)
# Add the unit level predictor
eklof$Datapoint<-as.factor(seq(1, dim(eklof)[1], 1))
# fit a MLMR - accouting for some non-independence
eklof_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~ Grazer.type, random=list(~1|ExptID,
~1|Datapoint), data=eklof)
results <- orchaRd::mod_results(eklof_MR, mod = "Grazer.type", group = "ExptID", data=eklof)

testthat::test_that("Checking mod_results output for eklof dataset ..", {

    testthat::expect_equal(
    dim(results)[[2]], 2,
    info = "mod_results output for eklof has correct dimensions")

  testthat::expect_equal(round(results$mod_table[1,2],2), round(-0.8095289, 2), info = "checking mod_results output for eklof calculates correct mean for Amphipod group...")

  testthat::expect_equal(round(results$mod_table[1,4],2), round(-0.20206548, 2), info = "checking mod_results output for eklof calculates correct upper confidence interval for Amphipod group...")

  testthat::expect_equal(round(results$mod_table[1,5],2), round(-3.021706, 2), info = "checking mod_results output for eklof calculates correct lower predcition interval for Amphipod group...")
})


# Fish example demonstrating marginalised means
data(fish)
warm_dat <- fish
model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,random = list(~1 | group_ID, ~1 | es_ID),
mods = ~ experimental_design + trait.type + deg_dif + treat_end_days, method = "REML", test = "t", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))

# Intercept only
 overall <- mod_results(model, group = "group_ID", data = warm_dat)

# Moderators
 #Trait type
across_trait <- mod_results(model, group = "group_ID", mod = "trait.type", data = warm_dat)

# Trait by degrees
across_trait_by_degree_diff <- mod_results(model, group = "group_ID",
mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", data = warm_dat)

# trait by condition
across_trait_by_treat_end_days10And50 <- mod_results(model, group = "group_ID",
mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)),
by = "treat_end_days", data = warm_dat)


# Fish data example with a heteroscedastic error
  model_het <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 + trait.type| es_ID), mods = ~ trait.type + deg_dif, method = "REML", test = "t", rho = 0, struc = "HCS", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))

  HetModel <- mod_results(model_het, group = "group_ID", mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", weights = "prop", data = warm_dat)
  orchard_plot(HetModel, xlab = "lnRR", data = warm_dat)

# tests for fish
testthat::test_that("Checking mod_results output for fish dataset ..", {

  testthat::expect_equal(
    dim(overall$mod_table)[[1]], 1,
    info = "mod_results output for fish 'overall' estimates has correct dimensions")

  testthat::expect_equal(
    dim(across_trait$mod_table)[[1]], 4,
    info = "mod_results output for fish 'across_trait' estimates has correct dimensions")

  testthat::expect_equal(
    dim(across_trait_by_degree_diff$mod_table)[[1]], 12,
    info = "mod_results output for fish 'across_trait_by_degree_diff' estimates has correct dimensions")


  testthat::expect_equal(
    dim(across_trait_by_treat_end_days10And50$mod_table)[[1]], 8,
    info = "mod_results output for fish 'across_trait_by_treat_end_days10And50' estimates has correct dimensions")

  testthat::expect_equal(
    dim(HetModel$mod_table)[[1]], 12,
    info = "mod_results output for fish 'HetModel' estimates has correct dimensions")

  testthat::expect_equal(round(HetModel$mod_table[2,3],2), round(0.736814215, 2), info = "checking mod_results output for HetModel calculates correct mean for Life-history group...")

  testthat::expect_equal(round(overall$mod_table[1,4],2), round(0.04375941, 2), info = "checking mod_results output for overall calculates correct upper confidence interval for Amphipod group...")

  testthat::expect_equal(round(across_trait_by_treat_end_days10And50$mod_table[1,6],2), round(-0.4187192, 2), info = "checking mod_results output for across_trait_by_treat_end_days10And50 calculates correct lower predcition interval for Behaviour group...")
})

