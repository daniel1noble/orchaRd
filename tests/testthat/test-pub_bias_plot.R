context("Checking pub_bias_plot..")

options(warn = -1)

# --- Shared fixtures ---
data(english)
english <- metafor::escalc(measure = "SMD",
                           n1i = NStartControl, sd1i = SD_C, m1i = MeanC,
                           n2i = NStartExpt,    sd2i = SD_E, m2i = MeanE,
                           var.names = c("SMD", "vSMD"),
                           data = english)

# Multi-level meta-analytic model (intercept-only) for the base orchard plot.
english_MA <- metafor::rma.mv(yi = SMD, V = vSMD,
                              random = list(~1 | StudyNo, ~1 | EffectID),
                              test = "t", data = english)

# Step 1: fixed-effect model + robust correction (Yang et al. 2023)
english_FE <- metafor::rma(yi = SMD, vi = vSMD, data = english,
                           test = "t", method = "FE")
english_FE_robust <- metafor::robust(english_FE, cluster = english$StudyNo,
                                     clubSandwich = TRUE)

# Modified Egger model used for the optional v_model branch
english_egger <- metafor::rma.mv(yi = SMD, V = vSMD, mods = ~vSMD,
                                 random = list(~1 | StudyNo, ~1 | EffectID),
                                 test = "t", data = english)

# Base orchard plot to decorate
base_plot <- orchaRd::orchard_plot(english_MA, group = "StudyNo",
                                   xlab = "SMD")

testthat::test_that("pub_bias_plot returns a ggplot when given only fe_model", {
  p <- orchaRd::pub_bias_plot(base_plot, english_FE_robust)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("pub_bias_plot returns a ggplot when given fe_model and v_model", {
  p <- orchaRd::pub_bias_plot(base_plot, english_FE_robust, english_egger)
  testthat::expect_s3_class(p, "ggplot")
})

testthat::test_that("pub_bias_plot errors when fe_model is not intercept-only", {
  testthat::expect_error(
    orchaRd::pub_bias_plot(base_plot, english_egger),
    "intercept only", fixed = TRUE,
    info = "fe_model with multiple coefficients is rejected"
  )
})

testthat::test_that("pub_bias_plot honors col / plotadj / textadj / *.size args", {
  # Non-default styling args still produce a valid ggplot — exercises the same
  # code paths but ensures the geom builders accept the parameters.
  p <- orchaRd::pub_bias_plot(base_plot, english_FE_robust,
                              col = c("forestgreen", "purple"),
                              plotadj = -0.1, textadj = 0.1,
                              branch.size = 1.5, trunk.size = 4)
  testthat::expect_s3_class(p, "ggplot")
})

# --- Direct tests on the helpers ---

testthat::test_that("get_ints_dat returns expected list shape for both 'br' and 'bc' types", {
  br <- get_ints_dat(english_FE_robust, type = "br")
  bc <- get_ints_dat(english_FE_robust, type = "bc")

  testthat::expect_length(br, 2)
  testthat::expect_s3_class(br[[1]], "data.frame")
  testthat::expect_named(br[[1]], c("name", "pred", "ci.lb", "ci.ub"))
  testthat::expect_match(br[[2]], "Bias Robust", fixed = TRUE)
  testthat::expect_match(bc[[2]], "Bias Corrected", fixed = TRUE)
})

testthat::test_that("get_ints_dat rejects unknown 'type' values via match.arg", {
  testthat::expect_error(get_ints_dat(english_FE_robust, type = "nope"))
})

testthat::test_that("geom_pub_stats_yang and geom_pub_stats_naka return ggplot layer lists", {
  br_dat <- get_ints_dat(english_FE_robust, type = "br")
  bc_dat <- get_ints_dat(english_egger,     type = "bc")

  layers_y <- geom_pub_stats_yang(br_dat)
  layers_n <- geom_pub_stats_naka(bc_dat)

  testthat::expect_type(layers_y, "list")
  testthat::expect_type(layers_n, "list")
  testthat::expect_length(layers_y, 3)
  testthat::expect_length(layers_n, 3)
  # Each layer should attach cleanly onto the base plot
  testthat::expect_s3_class(base_plot + layers_y, "ggplot")
  testthat::expect_s3_class(base_plot + layers_n, "ggplot")
})
