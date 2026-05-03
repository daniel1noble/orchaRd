
test_that("cor_diff returns expected structure and positive variance (numeric path)", {
  out <- cor_diff(cor1 = 0.5, cor2 = 0.3, n1 = 100, n2 = 120)
  expect_s3_class(out, "data.frame")
  expect_named(out, c("zr_diff", "var_zr_diff"))
  expect_true(is.finite(out$zr_diff))
  expect_true(is.finite(out$var_zr_diff))
  expect_gt(out$var_zr_diff, 0)
})

test_that("cor_diff matches Fisher z arithmetic with numeric inputs", {
  r1 <- 0.5; r2 <- 0.3; n1 <- 100; n2 <- 120
  out <- cor_diff(cor1 = r1, cor2 = r2, n1 = n1, n2 = n2)
  expect_equal(out$zr_diff, atanh(r1) - atanh(r2), tolerance = 1e-12)
  expect_equal(out$var_zr_diff, 1/(n1 - 3) + 1/(n2 - 3), tolerance = 1e-12)
})

test_that("swapping groups flips the sign but not the variance", {
  a <- cor_diff(0.2, 0.6, 80, 150)
  b <- cor_diff(0.6, 0.2, 150, 80)
  expect_equal(a$zr_diff, -b$zr_diff, tolerance = 1e-12)
  expect_equal(a$var_zr_diff, b$var_zr_diff, tolerance = 1e-12)
})

test_that("missing sample sizes with numeric correlations throws a clear error", {
  expect_error(
  cor_diff(cor1 = 0.2, cor2 = 0.3, n1 = 50, n2 = NULL),
  "Provide either \\(cor1, cor2, n1, n2\\) or two data frames x1 and x2",
  perl = TRUE
)
})

test_that("data-frame interface (MASS::mvrnorm): known correlations recovered and variance uses nrow", {
  skip_if_not_installed("MASS")
  set.seed(123)

  # Helper to make a 2-col data frame with target correlation rho
  make_df <- function(n, rho) {
    Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
    xy <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
    data.frame(x = xy[, 1], y = xy[, 2])
  }

  n1 <- 600; n2 <- 800
  rho1 <- 0.65; rho2 <- 0.25
  df1 <- make_df(n1, rho1)
  df2 <- make_df(n2, rho2)

  # Function path using data frames
  out_df <- cor_diff(x1 = df1, x2 = df2, cor1 = 0, cor2 = 0, n1 = 0, n2 = 0)

  # Numeric path computed from the realized sample correlations and true sample sizes
  r1_hat <- cor(df1$x, df1$y)
  r2_hat <- cor(df2$x, df2$y)
  out_num <- cor_diff(cor1 = r1_hat, cor2 = r2_hat, n1 = n1, n2 = n2)

  # z-diff should agree between interfaces
  expect_equal(out_df$zr_diff, out_num$zr_diff, tolerance = 1e-10)

  # Variance MUST be 1/(n1-3) + 1/(n2-3); this catches using length(x) instead of nrow(x)
  expect_equal(out_df$var_zr_diff, 1/(n1 - 3) + 1/(n2 - 3), tolerance = 1e-10)
})

test_that("pairwise.complete.obs works with NAs (MASS::mvrnorm) and variance stays positive", {
  skip_if_not_installed("MASS")
  set.seed(456)

  make_df <- function(n, rho) {
    Sigma <- matrix(c(1, rho, rho, 1), nrow = 2)
    xy <- MASS::mvrnorm(n = n, mu = c(0, 0), Sigma = Sigma)
    data.frame(x = xy[, 1], y = xy[, 2])
  }

  n <- 500
  df1 <- make_df(n, 0.5)
  df2 <- make_df(n, 0.5)

  # Introduce some missingness
  df1$x[sample.int(n, 30)] <- NA
  df1$y[sample.int(n, 20)] <- NA
  df2$x[sample.int(n, 25)] <- NA
  df2$y[sample.int(n, 25)] <- NA

  out <- cor_diff(x1 = df1, x2 = df2, cor1 = 0, cor2 = 0, n1 = 0, n2 = 0)
  expect_true(is.finite(out$zr_diff))
  expect_true(is.finite(out$var_zr_diff))
  expect_gt(out$var_zr_diff, 0)
})

test_that("data-frame inputs must have exactly two columns", {
  df1_bad <- data.frame(a = rnorm(10))
  df2_bad <- data.frame(a = rnorm(10), b = rnorm(10), c = rnorm(10))
  expect_error(
    cor_diff(x1 = df1_bad, x2 = df2_bad, cor1 = 0.2, cor2 = 0.3, n1 = 10, n2 = 10),
    "must each have exactly two columns",
    fixed = TRUE
  )
})

test_that("non-data-frame x1/x2 inputs are rejected", {
  # Passing vectors instead of data frames triggers the !is.data.frame branch
  expect_error(
    cor_diff(x1 = rnorm(10), x2 = rnorm(10),
             cor1 = 0.2, cor2 = 0.3, n1 = 10, n2 = 10),
    "both must be data frames", fixed = TRUE
  )
})

test_that("cor_diff errors when sample sizes are too small for Fisher z variance", {
  expect_error(
    cor_diff(cor1 = 0.5, cor2 = 0.3, n1 = 3, n2 = 50),
    "Sample sizes must be > 3", fixed = TRUE
  )
  expect_error(
    cor_diff(cor1 = 0.5, cor2 = 0.3, n1 = 50, n2 = 2),
    "Sample sizes must be > 3", fixed = TRUE
  )
})

test_that("cor_diff errors when correlations are not finite", {
  expect_error(
    cor_diff(cor1 = NA_real_, cor2 = 0.3, n1 = 50, n2 = 50),
    "must be finite", fixed = TRUE
  )
  expect_error(
    cor_diff(cor1 = 0.5, cor2 = Inf, n1 = 50, n2 = 50),
    "must be finite", fixed = TRUE
  )
})