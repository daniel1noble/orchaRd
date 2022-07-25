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
 overall <- orchaRd::mod_results(model, group = "group_ID", data = warm_dat)

# Moderators
 #Trait type
across_trait <- orchaRd::mod_results(model, group = "group_ID", mod = "trait.type", data = warm_dat)

# Trait by degrees
across_trait_by_degree_diff <- orchaRd::mod_results(model, group = "group_ID",
mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", data = warm_dat)

# trait by condition
across_trait_by_treat_end_days10And50 <- orchaRd::mod_results(model, group = "group_ID",
mod = "trait.type", at = list(deg_dif = c(5, 10, 15), treat_end_days = c(10, 50)),
by = "treat_end_days", data = warm_dat)


# Fish data example with a heteroscedastic error
  model_het <- metafor::rma.mv(yi = lnrr, V = lnrr_vi, random = list(~1 | group_ID, ~1 + trait.type| es_ID), mods = ~ trait.type + deg_dif, method = "REML", test = "t", rho = 0, struc = "HCS", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))

  HetModel <- orchaRd::mod_results(model_het, group = "group_ID", mod = "trait.type", at = list(deg_dif = c(5, 10, 15)), by = "deg_dif", weights = "prop", data = warm_dat)
  orchard_plot(HetModel, xlab = "lnRR", data = warm_dat)


  ## English example
  data(english)

  # We need to calculate the effect sizes, in this case d
  english <- metafor::escalc(measure = "SMD", n1i = NStartControl, sd1i = SD_C, m1i = MeanC,
                    n2i = NStartExpt, sd2i = SD_E, m2i = MeanE,
                    var.names=c("SMD","vSMD"),
                    data = english)

  english_MA <- metafor::rma.mv(yi = SMD, V = vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)
  eng_mod_results <- orchaRd::mod_results(english_MA, group = "StudyNo", data = english)


  english_MR0 <- metafor::rma.mv(yi = SMD, V = vSMD, mods = ~ ManipType,
                        random = list(~ 1 | StudyNo, ~ 1 | EffectID), data = english)

  english_MR0_noint <- metafor::rma.mv(yi = SMD, V = vSMD, mods = ~ -1 + ManipType,
                                 random = list(~ 1 | StudyNo, ~ 1 | EffectID),
                                 data = english)

  # Again, we can create a table of results
  res2 <- orchaRd::mod_results(english_MR0, mod = "ManipType", data = english, group = "StudyNo")

  # Check supressing level works
  res3 <- orchaRd::mod_results(english_MR0, mod = "ManipType", data = english, group = "StudyNo",
                               at = list(ManipType = "Quality"), subset = TRUE)

# tests for fish
testthat::test_that("Checking mod_results output for fish dataset ..", {

  testthat::expect_equal(
    dim(overall$mod_table)[[1]], 1,
    info = "mod_results output for fish 'overall' object has correct dimensions")

  testthat::expect_equal(
    dim(across_trait$mod_table)[[1]], 4,
    info = "mod_results output for fish 'across_trait' object has correct dimensions")

  testthat::expect_equal(
    dim(across_trait_by_degree_diff$mod_table)[[1]], 12,
    info = "mod_results output for fish 'across_trait_by_degree_diff' object has correct dimensions")


  testthat::expect_equal(
    dim(across_trait_by_treat_end_days10And50$mod_table)[[1]], 8,
    info = "mod_results output for fish 'across_trait_by_treat_end_days10And50' object has correct dimensions")

  testthat::expect_equal(
    dim(HetModel$mod_table)[[1]], 12,
    info = "mod_results output for fish 'HetModel' object has correct dimensions")

  testthat::expect_equal(
    dim(eng_mod_results$mod_table)[[1]], 1,
    info = "mod_results output for english 'eng_mod_results' object has correct dimensions")

  testthat::expect_equal(
    dim(res2$mod_table)[[1]], 2,
    info = "mod_results output for english 'res2' object has correct dimensions")

  testthat::expect_equal(
    dim(res3$mod_table)[[1]], 1,
    info = "mod_results output for english 'res3' to check subsetting is working and object has correct dimensions")

  testthat::expect_equal(round(res2$mod_table[1,2],2), round(res3$mod_table[1,2], 2), info = "checking mod_results output for english res2 object that has two levels matches the subsetted 'res3' object that only has a single level, Quality...")

  testthat::expect_equal(round(HetModel$mod_table[2,3],2), round(0.736814215, 2), info = "checking mod_results output for HetModel calculates correct mean for Life-history group...")

  testthat::expect_equal(round(overall$mod_table[1,4],2), round(0.04375941, 2), info = "checking mod_results output for overall calculates correct upper confidence interval for Amphipod group...")

  testthat::expect_equal(round(across_trait_by_treat_end_days10And50$mod_table[1,6],2), round(-0.4187192, 2), info = "checking mod_results output for across_trait_by_treat_end_days10And50 calculates correct lower predcition interval for Behaviour group...")


  testthat::expect_warning(orchaRd::mod_results(english_MR0_noint, mod = "ManipType", data = english, group = "StudyNo"))
})


### random-effects model; note that this is obviously not the right model for the data, but it's just testing this class type
mod_uni <- metafor::rma(lnrr, lnrr_vi, data=warm_dat, method="REML")

mod_rs_uni <- orchaRd::mod_results(mod_uni, group = "group_ID", data=warm_dat)

plot_uni <- orchaRd::orchard_plot(mod_uni, group = "group_ID", data=warm_dat, xlab = "Effect Size")

testthat::test_that("Checking mod_results output for rma and rma.uni dataset ..", {

  testthat::expect_equal(
    as.numeric(mod_uni$b), mod_rs_uni$mod_table$estimate,
    info = "mod_results output for rma and rma.uni works fine and matches mod_results...")

})

