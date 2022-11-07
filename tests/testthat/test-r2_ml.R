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
