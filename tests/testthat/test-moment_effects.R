
test_that("moment_effects returns correctly named columns", {
  set.seed(1)
  x1 <- rnorm(50)
  x2 <- rnorm(50)

  out_skew <- moment_effects(x1, x2, type = "skew")
  expect_s3_class(out_skew, "data.frame")
  expect_named(out_skew, c("d_skew", "d_skew_v"))
  expect_true(is.finite(out_skew$d_skew_v))
  expect_gt(out_skew$d_skew_v, 0)

  out_kurt <- moment_effects(x1, x2, type = "kurt")
  expect_s3_class(out_kurt, "data.frame")
  expect_named(out_kurt, c("d_kurt", "d_kurt_v"))
  expect_true(is.finite(out_kurt$d_kurt_v))
  expect_gt(out_kurt$d_kurt_v, 0)
})

test_that("type argument must be one of 'skew' or 'kurt'", {
  set.seed(2)
  x1 <- rnorm(10)
  x2 <- rnorm(10)

  expect_error(moment_effects(x1, x2, type = "mean"))
})

test_that("swapping groups flips the sign of the effect", {
  set.seed(3)
  x1 <- rnorm(80)
  x2 <- rnorm(80)

  d1 <- moment_effects(x1, x2, type = "skew")$d_skew
  d2 <- moment_effects(x2, x1, type = "skew")$d_skew
  expect_equal(d1, -d2, tolerance = 1e-12)

  k1 <- moment_effects(x1, x2, type = "kurt")$d_kurt
  k2 <- moment_effects(x2, x1, type = "kurt")$d_kurt
  expect_equal(k1, -k2, tolerance = 1e-12)
})

test_that("jackknife variance increases with an extreme outlier", {
  set.seed(4)
  x <- rnorm(60)
  x_out <- c(x, 25)  # add a large outlier

  v_skew_base <- .sv_jack(x, type = "skew")$var
  v_skew_out  <- .sv_jack(x_out, type = "skew")$var
  expect_gt(v_skew_out, v_skew_base)

  v_kurt_base <- .sv_jack(x, type = "kurt")$var
  v_kurt_out  <- .sv_jack(x_out, type = "kurt")$var
  expect_gt(v_kurt_out, v_kurt_base)
})

test_that("known differences from Pearson distributions are recovered (skew ≈ -1, kurt ≈ 0)", {
  skip_if_not_installed("PearsonDS")
  set.seed(982)

  # m1: skew=0, kurtosis(excess)=3  ; m2: skew=1, kurtosis(excess)=3
  m1 <- c(mean = 0, variance = 1, skewness = 0, kurtosis = 3)
  m2 <- c(mean = 0, variance = 1, skewness = 1, kurtosis = 3)

  x1 <- PearsonDS::rpearson(3000, moments = m1)
  x2 <- PearsonDS::rpearson(3000, moments = m2)

  # Difference in skew should be ~ -1; difference in kurtosis ~ 0
  d_skew <- moment_effects(x1, x2, type = "skew")$d_skew
  d_kurt <- moment_effects(x1, x2, type = "kurt")$d_kurt

  expect_equal(d_skew, -1, tolerance = 0.15)
  expect_equal(d_kurt,  0,  tolerance = 0.20)

  # Variances should be positive and finite
  v_skew <- moment_effects(x1, x2, type = "skew")$d_skew_v
  v_kurt <- moment_effects(x1, x2, type = "kurt")$d_kurt_v
  expect_gt(v_skew, 0)
  expect_gt(v_kurt, 0)
  expect_true(is.finite(v_skew))
  expect_true(is.finite(v_kurt))
})

test_that("known differences from Pearson distributions are recovered (kurt ≈ -2, skew ≈ 0)", {
  skip_if_not_installed("PearsonDS")
  set.seed(981)

  # m1: skew=0, kurtosis(excess)=2 ; m2: skew=0, kurtosis(excess)=4
  # -> d_kurt should be about 2 - 4 = -2 ; d_skew ≈ 0
  m1 <- c(mean = 0, variance = 1, skewness = 0, kurtosis = 2)
  m2 <- c(mean = 0, variance = 1, skewness = 0, kurtosis = 4)

  x1 <- PearsonDS::rpearson(3000, moments = m1)
  x2 <- PearsonDS::rpearson(3000, moments = m2)

  d_skew <- moment_effects(x1, x2, type = "skew")$d_skew
  d_kurt <- moment_effects(x1, x2, type = "kurt")$d_kurt

  expect_equal(d_skew, 0,   tolerance = 0.12)
  expect_equal(d_kurt, -2,  tolerance = 0.20)
})