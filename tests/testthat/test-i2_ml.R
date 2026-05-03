context("Checking i2_ml function..")

options(warn = -1)
data(lim)
lim$vi<-(1/sqrt(lim$N - 3))^2

# Check two types of model specification. One with formula and one without
lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data = lim)

lim_MR2<-metafor::rma.mv(yi~1, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data = lim)

# Check point estimates
 I2_lim_2 <- orchaRd::i2_ml(lim_MR)
I2_lim_22 <- orchaRd::i2_ml(lim_MR2)

# Checking boot
 I2_lim_2_boot <- orchaRd::i2_ml(lim_MR,  boot =10)
I2_lim_22_boot <- orchaRd::i2_ml(lim_MR2,  boot = 10)

# Check it works with robust models
lim_MR_robust <- metafor::robust(lim_MR, cluster=lim$Article)
I2_lim_robust <- orchaRd::i2_ml(lim_MR_robust)

testthat::test_that("Checking i2 has correct dimensions..", {
  testthat::expect_equal(
    length(I2_lim_2), 3,
    info = "i2 has correct dimensions")

  testthat::expect_equal(
    dim(I2_lim_2_boot)[1], 3,
    info = "i2 boot has correct dimensions")

  testthat::expect_equal(
    dim(I2_lim_22_boot)[1], 3,
    info = "i2 boot works with mod formula specification")

  testthat::expect_equal(
    length(I2_lim_robust), 3,
    info = "i2 on robust object has correct dimensions")

testthat::expect_equal(round(as.vector(I2_lim_2),2), round(c(84.32574, 48.19211, 36.13363), 2), info = "checking i2 estimates calculated correctly...")
})

# --- Cover the matrix-method bootstrap branches (with and without mods)
#     plus the missing-data (model$not.na) branch and explicit error paths. ---

set.seed(123)
I2_matrix_pt   <- orchaRd::i2_ml(lim_MR,  method = "matrix")
I2_matrix_pt0  <- orchaRd::i2_ml(lim_MR2, method = "matrix")

set.seed(123)
I2_matrix_b    <- orchaRd::i2_ml(lim_MR,  method = "matrix", boot = 2)
I2_matrix_b0   <- orchaRd::i2_ml(lim_MR2, method = "matrix", boot = 2)

# Missing-data model exercises the data <- data[model$not.na, ] branch
lim_na    <- lim
lim_na$yi[1:5] <- NA
lim_MR_na <- metafor::rma.mv(yi = yi, V = vi, mods = ~ Phylum - 1,
                             random = list(~1 | Article, ~1 | Datapoint),
                             data = lim_na)
set.seed(123)
I2_lim_na_boot <- orchaRd::i2_ml(lim_MR_na, boot = 2)

testthat::test_that("i2_ml matrix method works as point estimate and under bootstrap", {
  testthat::expect_length(I2_matrix_pt, 3)
  testthat::expect_length(I2_matrix_pt0, 3)
  testthat::expect_equal(dim(I2_matrix_b),  c(3, 3))
  testthat::expect_equal(dim(I2_matrix_b0), c(3, 3))
  testthat::expect_equal(colnames(I2_matrix_b), c("Est.", "2.5%", "97.5%"))
})

testthat::test_that("i2_ml respects model$not.na when data has missing values", {
  testthat::expect_equal(dim(I2_lim_na_boot), c(3, 3))
})

testthat::test_that("i2_ml / matrix_i2 / ratio_i2 reject non-metafor models", {
  bogus <- stats::lm(yi ~ Phylum, data = lim)
  testthat::expect_error(orchaRd::i2_ml(bogus))
  testthat::expect_error(orchaRd::matrix_i2(bogus))
  testthat::expect_error(orchaRd::ratio_i2(bogus))
})

testthat::test_that("i2_ml rejects models with heterogeneous variance (tau2 > 0)", {
  data(fish)
  m_tau2 <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,
                            mods = ~ trait.type + deg_dif,
                            random = list(~1 | group_ID, ~1 + trait.type | es_ID),
                            struc = "HCS", data = fish, method = "REML",
                            control = list(optimizer = "optim", optmethod = "Nelder-Mead"))
  testthat::expect_error(orchaRd::i2_ml(m_tau2),
    "heterogeneous variance", fixed = TRUE)
})
