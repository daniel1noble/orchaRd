context("Checking cvh1_ml / cvh2_ml / m1_ml / m2_ml functions..")

options(warn = -1)
data(lim)
lim$vi <- (1/sqrt(lim$N - 3))^2

# Two model fixtures so the bootstrap branches that split on
# is.null(mods_formula) are both exercised.
lim_MR <- metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1,
                          random=list(~1|Article, ~1|Datapoint), data=lim)
lim_M0 <- metafor::rma.mv(yi=yi, V=vi,
                          random=list(~1|Article, ~1|Datapoint), data=lim)

# Robust wrapper to confirm the class-check branch accepts robust.rma.
lim_robust <- metafor::robust(lim_MR, cluster=lim$Article)

# Non-metafor object for the wrong-class error branch.
not_a_model <- stats::lm(yi ~ Phylum, data=lim)

# --- Point estimates (no-boot) ---
cvh1_pt <- orchaRd::cvh1_ml(lim_MR)
cvh2_pt <- orchaRd::cvh2_ml(lim_MR)
m1_pt   <- orchaRd::m1_ml(lim_MR)
m2_pt   <- orchaRd::m2_ml(lim_MR)

# --- Bootstrap (with mods) ---
set.seed(42); cvh1_b  <- orchaRd::cvh1_ml(lim_MR, boot=2)
set.seed(42); cvh2_b  <- orchaRd::cvh2_ml(lim_MR, boot=2)
set.seed(42); m1_b    <- orchaRd::m1_ml(lim_MR,   boot=2)
set.seed(42); m2_b    <- orchaRd::m2_ml(lim_MR,   boot=2)

# --- Bootstrap (no mods) -- exercises the is.null(mods_formula) TRUE branch ---
set.seed(42); cvh1_b0 <- orchaRd::cvh1_ml(lim_M0, boot=2)
set.seed(42); cvh2_b0 <- orchaRd::cvh2_ml(lim_M0, boot=2)
set.seed(42); m1_b0   <- orchaRd::m1_ml(lim_M0,   boot=2)
set.seed(42); m2_b0   <- orchaRd::m2_ml(lim_M0,   boot=2)

testthat::test_that("cvh1_ml point estimates and structure", {
  testthat::expect_equal(length(cvh1_pt), 3,
    info = "cvh1_ml returns Total + one entry per random effect")
  testthat::expect_equal(round(as.vector(cvh1_pt), 2),
    round(c(0.9976791, 0.7542212, 0.6530803), 2),
    info = "cvh1_ml point estimates match reference values")
  testthat::expect_equal(names(cvh1_pt),
    c("CVH1_Total", "CVH1_Article", "CVH1_Datapoint"),
    info = "cvh1_ml names follow CVH1_<level> convention")
})

testthat::test_that("cvh1_ml bootstrap structure (with and without mods)", {
  testthat::expect_equal(dim(cvh1_b),  c(3, 3),
    info = "cvh1_ml boot returns 3 rows x 3 quantile columns when mods present")
  testthat::expect_equal(colnames(cvh1_b), c("Est.", "2.5%", "97.5%"),
    info = "cvh1_ml boot has expected column names")
  testthat::expect_equal(dim(cvh1_b0), c(3, 3),
    info = "cvh1_ml boot returns same shape for intercept-only models")
})

testthat::test_that("cvh1_ml accepts robust.rma and rejects non-metafor objects", {
  testthat::expect_silent(orchaRd::cvh1_ml(lim_robust))
  testthat::expect_error(orchaRd::cvh1_ml(not_a_model),
    info = "cvh1_ml errors on non-metafor model classes")
})

testthat::test_that("cvh2_ml point estimates and structure", {
  testthat::expect_equal(length(cvh2_pt), 3,
    info = "cvh2_ml returns Total + one entry per random effect")
  testthat::expect_equal(round(as.vector(cvh2_pt), 2),
    round(c(0.9953635, 0.5688496, 0.4265139), 2),
    info = "cvh2_ml point estimates match reference values")
  testthat::expect_equal(names(cvh2_pt),
    c("CVH2_Total", "CVH2_Article", "CVH2_Datapoint"))
})

testthat::test_that("cvh2_ml bootstrap structure (with and without mods)", {
  testthat::expect_equal(dim(cvh2_b),  c(3, 3))
  testthat::expect_equal(colnames(cvh2_b), c("Est.", "2.5%", "97.5%"))
  testthat::expect_equal(dim(cvh2_b0), c(3, 3))
})

testthat::test_that("cvh2_ml rejects non-metafor objects", {
  testthat::expect_error(orchaRd::cvh2_ml(not_a_model))
})

testthat::test_that("m1_ml point estimates and structure", {
  testthat::expect_equal(length(m1_pt), 3)
  testthat::expect_equal(round(as.vector(m1_pt), 2),
    round(c(0.4994191, 0.3775487, 0.3269195), 2),
    info = "m1_ml point estimates match reference values")
  testthat::expect_equal(names(m1_pt),
    c("M1_Total", "M1_Article", "M1_Datapoint"))
})

testthat::test_that("m1_ml bootstrap structure (with and without mods)", {
  testthat::expect_equal(dim(m1_b),  c(3, 3))
  testthat::expect_equal(colnames(m1_b), c("Est.", "2.5%", "97.5%"))
  testthat::expect_equal(dim(m1_b0), c(3, 3))
})

testthat::test_that("m1_ml rejects non-metafor objects", {
  testthat::expect_error(orchaRd::m1_ml(not_a_model))
})

testthat::test_that("m2_ml point estimates and structure", {
  testthat::expect_equal(length(m2_pt), 3)
  testthat::expect_equal(round(as.vector(m2_pt), 2),
    round(c(0.4988382, 0.2850857, 0.2137525), 2),
    info = "m2_ml point estimates match reference values")
  testthat::expect_equal(names(m2_pt),
    c("M2_Total", "M2_Article", "M2_Datapoint"))
})

testthat::test_that("m2_ml bootstrap structure (with and without mods)", {
  testthat::expect_equal(dim(m2_b),  c(3, 3))
  testthat::expect_equal(colnames(m2_b), c("Est.", "2.5%", "97.5%"))
  testthat::expect_equal(dim(m2_b0), c(3, 3))
})

testthat::test_that("m2_ml rejects non-metafor objects", {
  testthat::expect_error(orchaRd::m2_ml(not_a_model))
})

testthat::test_that("ml_* helpers compute the same values as their wrappers", {
  testthat::expect_equal(orchaRd::ml_cvh1(lim_MR), cvh1_pt)
  testthat::expect_equal(orchaRd::ml_cvh2(lim_MR), cvh2_pt)
  testthat::expect_equal(orchaRd::ml_m1(lim_MR),   m1_pt)
  testthat::expect_equal(orchaRd::ml_m2(lim_MR),   m2_pt)
})

# --- Heterogeneous-variance models (tau2 > 0) are deliberately rejected
#     with an explicit message; build one fixture and check all four wrappers. ---
data(fish)
m_tau2 <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
                          mods = ~ trait.type + deg_dif,
                          random = list(~1 | group_ID, ~1 + trait.type | es_ID),
                          struc = "HCS", data = fish, method = "REML",
                          control = list(optimizer = "optim", optmethod = "Nelder-Mead"))

testthat::test_that("cvh/m wrappers reject models with heterogeneous variance (tau2 > 0)", {
  testthat::expect_error(orchaRd::cvh1_ml(m_tau2), "heterogeneous variance", fixed = TRUE)
  testthat::expect_error(orchaRd::cvh2_ml(m_tau2), "heterogeneous variance", fixed = TRUE)
  testthat::expect_error(orchaRd::m1_ml(m_tau2),   "heterogeneous variance", fixed = TRUE)
  testthat::expect_error(orchaRd::m2_ml(m_tau2),   "heterogeneous variance", fixed = TRUE)
})

# --- Bootstrap path with NAs in yi exercises the model$not.na branch
#     (data <- data[model$not.na, ]). ---
data(lim)
lim_na <- lim
lim_na$vi <- (1/sqrt(lim_na$N - 3))^2
lim_na$yi[1:5] <- NA
lim_MR_na <- metafor::rma.mv(yi = yi, V = vi, mods = ~ Phylum - 1,
                             random = list(~1 | Article, ~1 | Datapoint),
                             data = lim_na)

testthat::test_that("cvh1_ml bootstrap respects model$not.na for missing-data models", {
  set.seed(42)
  out <- orchaRd::cvh1_ml(lim_MR_na, boot = 2)
  testthat::expect_equal(dim(out), c(3, 3))
})
