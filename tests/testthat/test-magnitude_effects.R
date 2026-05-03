context("Checking magnitude_effects + .safe_lnM_* SAFE bootstrap functions..")

# Tuned-down min_kept / chunk_init to keep test runtime negligible.
# Defaults are min_kept = 2000 / chunk_init = 4000 which is excessive for tests.

testthat::test_that("magnitude_effects (independent samples) returns expected structure", {
  res <- orchaRd::magnitude_effects(x1bar = 10, x2bar = 7, sd1 = 2, sd2 = 3,
                                    n1 = 30, n2 = 30,
                                    min_kept = 50, chunk_init = 200, seed = 42)
  testthat::expect_s3_class(res, "data.frame")
  testthat::expect_named(res,
    c("lnM_SAFE", "var_lnM_SAFE", "kept", "total", "attempts", "status"))
  testthat::expect_true(is.finite(res$lnM_SAFE))
  testthat::expect_gt(res$var_lnM_SAFE, 0)
  testthat::expect_gte(res$kept, 50)
  testthat::expect_identical(res$status, "ok")
})

testthat::test_that("magnitude_effects (dependent samples) returns expected structure", {
  res <- orchaRd::magnitude_effects(x1bar = 10, x2bar = 7, sd1 = 2, sd2 = 3,
                                    n1 = 30, n2 = 30, paired = TRUE, r = 0.5,
                                    min_kept = 50, chunk_init = 200, seed = 42)
  testthat::expect_s3_class(res, "data.frame")
  testthat::expect_named(res,
    c("lnM_SAFE", "var_lnM_SAFE", "kept", "total", "attempts", "status"))
  testthat::expect_true(is.finite(res$lnM_SAFE))
  testthat::expect_gt(res$var_lnM_SAFE, 0)
})

testthat::test_that("magnitude_effects errors on missing/inconsistent paired inputs", {
  testthat::expect_error(
    orchaRd::magnitude_effects(x1bar = 10, x2bar = 7, sd1 = 2, sd2 = 3,
                               n1 = 30, n2 = 30, paired = TRUE, r = NULL),
    "correlation 'r'", fixed = TRUE,
    info = "paired = TRUE without r should error")
  testthat::expect_error(
    orchaRd::magnitude_effects(x1bar = 10, x2bar = 7, sd1 = 2, sd2 = 3,
                               n1 = 30, n2 = 25, paired = TRUE, r = 0.5),
    "n1 and n2 must be equal", fixed = TRUE,
    info = "paired = TRUE with unequal n should error")
})

testthat::test_that("magnitude_effects results are reproducible under the same seed", {
  res1 <- orchaRd::magnitude_effects(x1bar = 10, x2bar = 7, sd1 = 2, sd2 = 3,
                                     n1 = 30, n2 = 30,
                                     min_kept = 50, chunk_init = 200, seed = 99)
  res2 <- orchaRd::magnitude_effects(x1bar = 10, x2bar = 7, sd1 = 2, sd2 = 3,
                                     n1 = 30, n2 = 30,
                                     min_kept = 50, chunk_init = 200, seed = 99)
  testthat::expect_equal(res1$lnM_SAFE, res2$lnM_SAFE)
  testthat::expect_equal(res1$var_lnM_SAFE, res2$var_lnM_SAFE)
})

testthat::test_that(".safe_lnM_indep handles 'no usable draws' status when MSB cannot exceed MSW", {
  # When the means are equal and the variance is large, MSB > MSW will rarely
  # (if ever) hold; force an early break by setting patience_noaccept = 1 and
  # capping max_draws so we don't loop forever.
  res <- orchaRd::.safe_lnM_indep(x1bar = 5, x2bar = 5, sd1 = 10, sd2 = 10,
                                  n1 = 5, n2 = 5,
                                  min_kept = 1e6, chunk_init = 50,
                                  max_draws = 200, patience_noaccept = 1, seed = 1)
  testthat::expect_true(res$status %in% c("ok", "no_usable_draws"),
    info = "status field is one of the documented values")
})

testthat::test_that(".safe_lnM_indep called directly returns same shape as the wrapper", {
  res <- orchaRd::.safe_lnM_indep(x1bar = 10, x2bar = 7, sd1 = 2, sd2 = 3,
                                  n1 = 30, n2 = 30,
                                  min_kept = 50, chunk_init = 200, seed = 42)
  testthat::expect_named(res,
    c("lnM_SAFE", "var_lnM_SAFE", "kept", "total", "attempts", "status"))
})

testthat::test_that(".safe_lnM_dep called directly returns same shape as the wrapper", {
  res <- orchaRd::.safe_lnM_dep(x1bar = 10, x2bar = 7, sd1 = 2, sd2 = 3,
                                n = 30, r = 0.5,
                                min_kept = 50, chunk_init = 200, seed = 42)
  testthat::expect_named(res,
    c("lnM_SAFE", "var_lnM_SAFE", "kept", "total", "attempts", "status"))
})
