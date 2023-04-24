context("Checking get_raw_data function...")

data(english)

# We need to calculate the effect sizes, in this case d
english <- metafor::escalc(measure = "SMD", n1i = NStartControl, sd1i = SD_C, m1i = MeanC, n2i = NStartExpt, sd2i = SD_E, m2i = MeanE,
                  var.names=c("SMD","vSMD"),
                  data = english)

english_MA <- metafor::rma.mv(yi = SMD, V = vSMD, random = list( ~ 1 | StudyNo, ~ 1 | EffectID),data = english)


english_test1 <- orchaRd::get_data_raw(english_MA, mod = "1", group = "StudyNo")

testthat::test_that("Checking get_raw_data output for english dataset ..", {

  testthat::expect_equal(
    dim(english_test1)[[1]], dim(english_test1)[[1]],
    info = "Check that english raw data dimensions match the dimensions of get_raw_data output...")

  testthat::expect_equal(
    names(english_test1), c("yi", "vi","moderator", "stdy", "type"),
    info = "Check that get_raw_data output for english has the correct column names...")


})
