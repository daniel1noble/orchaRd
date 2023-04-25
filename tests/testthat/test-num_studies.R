context("Checking num_studies function...")

data(english)

# We need to calculate the effect sizes, in this case d
english <- metafor::escalc(measure = "SMD", n1i = NStartControl, sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E, m2i = MeanE,
                           var.names=c("SMD","vSMD"),
                           data = english)

english_MA <- metafor::rma.mv(yi = SMD, V = vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID), data = english)

english_ModRes <- orchaRd::mod_results(english_MA, group = "StudyNo")

overall_english <- orchaRd::num_studies(english_ModRes$data, mod = moderator, group = stdy)

testthat::test_that("Checking num_studies output ...", {

  testthat::expect_equal(
    overall_english$Num_Studies, length(unique(english_ModRes$data$stdy)),
    info = "Check that unique number of studies matches num_studies output...")

})
