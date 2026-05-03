context("Checking r2_ml function..")

## English example
data(english)

# We need to calculate the effect sizes, in this case d
english <- metafor::escalc(measure = "SMD", n1i = NStartControl, sd1i = SD_C, m1i = MeanC,
                           n2i = NStartExpt, sd2i = SD_E, m2i = MeanE,
                           var.names=c("SMD","vSMD"),
                           data = english)

english_MA <- metafor::rma.mv(yi = SMD, V = vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)

# Switch order of obs to check it works still correctly
english_MA_obs <- metafor::rma.mv(yi = SMD, V = vSMD, random = list(~ 1 | EffectID, ~ 1 | StudyNo), data = english)


    english_MA_R2 <- orchaRd::r2_ml(english_MA, data = english)
english_MA_R2_obs <- orchaRd::r2_ml(english_MA_obs, data = english)

english_MA_robust <- metafor::robust(english_MA, cluster = english$StudyNo)

english_MA_R2_robu <- orchaRd::r2_ml(english_MA_robust, data = english)

testthat::test_that("Checking r2_ml function..", {
  testthat::expect_equal(length(english_MA_R2), 2,
                        info = "Checking r2 object has correct dimensions")

  testthat::expect_equal(length(english_MA_R2_robu), 2,
                         info = "Checking r2 robust object has correct dimensions")

  testthat::expect_equal(round(as.vector(english_MA_R2),2),
                         round(c(0.00, 0.45 ), 2),
                         info = "Checking R2 estimates are calculated correctly...")

  testthat::expect_equal(round(as.vector(english_MA_R2),2),
                         round(as.vector(english_MA_R2_obs),2),
                         info = "Checking R2_conditional estimates are calculated correctly and do not depend on order of random effect...")

  testthat::expect_error(r2_ml(english_MA_robust, data = english, boot = 10),
                         info = "Checking R2 estimates correctly throw errow...")
})

# --- Bootstrap path on lim data (matches the i2_ml fixture style) ---
options(warn = -1)
data(lim)
lim$vi <- (1/sqrt(lim$N - 3))^2
lim_MR <- metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1,
                          random=list(~1|Article, ~1|Datapoint), data=lim)

set.seed(42)
lim_R2_boot <- orchaRd::r2_ml(lim_MR, data=lim, boot=2)

testthat::test_that("r2_ml bootstrap path returns expected structure", {
  testthat::expect_equal(dim(lim_R2_boot), c(2, 3),
    info = "r2_ml boot returns 2 rows (marginal, conditional) x 3 quantile columns")
  testthat::expect_equal(colnames(lim_R2_boot), c("Est.", "2.5%", "97.5%"),
    info = "r2_ml boot column names follow the Est./2.5%/97.5% convention")
  testthat::expect_equal(rownames(lim_R2_boot),
    c("R2_marginal", "R2_conditional"),
    info = "r2_ml boot row names match the underlying R2_calc output")
})

testthat::test_that("R2_calc helper computes the same values as r2_ml's no-boot path", {
  testthat::expect_equal(orchaRd::R2_calc(lim_MR), orchaRd::r2_ml(lim_MR, data=lim))
  testthat::expect_error(orchaRd::R2_calc(stats::lm(yi ~ Phylum, data=lim)),
    info = "R2_calc errors when given a non-metafor model object")
})

testthat::test_that("r2_ml rejects unsupported model types", {
  testthat::expect_error(orchaRd::r2_ml(stats::lm(yi ~ Phylum, data=lim), data=lim),
    info = "r2_ml errors on non-metafor model classes")
})

testthat::test_that("r2_ml rejects models with heterogeneous variance (tau2 > 0)", {
  data(fish)
  m_tau2 <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
                            mods = ~ trait.type + deg_dif,
                            random = list(~1 | group_ID, ~1 + trait.type | es_ID),
                            struc = "HCS", data = fish, method = "REML",
                            control = list(optimizer = "optim", optmethod = "Nelder-Mead"))
  testthat::expect_error(orchaRd::r2_ml(m_tau2, data = fish),
    "heterogeneous variance", fixed = TRUE)
})
