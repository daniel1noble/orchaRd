test_that("PI matches metafor for sigma2-only model (no G/H structures)", {
  skip_on_cran()
  library(metafor)
  data("dat.konstantopoulos2011")
  df <- dat.konstantopoulos2011

  m <- rma.mv(yi ~ 1, vi, random = list(~1|study, ~1|district, ~1|school), data = df)
  orc <- mod_results(m, mod = "1", group = "study")

  pred <- predict(m)
  # Use mean of metafor PIs as the single-prediction baseline
  expect_equal(orc$mod_table$lowerPR, pred$pi.lb[1], tolerance = 0.01)
  expect_equal(orc$mod_table$upperPR, pred$pi.ub[1], tolerance = 0.01)
})

test_that("PI matches metafor for CS random slope model (single tau2)", {
  skip_on_cran()
  library(metafor)
  data("dat.konstantopoulos2011")
  df <- dat.konstantopoulos2011

  m <- rma.mv(yi ~ year, vi,
              random = list(~year|study, ~1|district, ~1|school),
              rho = 0, data = df, method = "REML")
  expect_length(m$tau2, 1)  # CS structure: single tau2

  orc <- mod_results(m, mod = "year", group = "study")
  # Compare at year = min(year) -- first prediction in orchaRd's grid
  yr1 <- orc$mod_table$name[1]
  pred <- predict(m, newmods = yr1, addx = TRUE)

  expect_equal(orc$mod_table$lowerPR[1], pred$pi.lb[1], tolerance = 0.01)
  expect_equal(orc$mod_table$upperPR[1], pred$pi.ub[1], tolerance = 0.01)
})

test_that("PI matches metafor for model with both tau2 and gamma2 (both scalar)", {
  skip_on_cran()
  library(metafor)
  data("dat.konstantopoulos2011")
  df <- dat.konstantopoulos2011

  m <- rma.mv(yi ~ year, vi,
              random = list(~year|study, ~year|district, ~1|school),
              rho = 0, phi = 0, data = df, method = "REML")
  expect_length(m$tau2, 1)
  expect_length(m$gamma2, 1)

  orc <- mod_results(m, mod = "year", group = "study")
  yr1 <- orc$mod_table$name[1]
  pred <- predict(m, newmods = yr1, addx = TRUE)

  expect_equal(orc$mod_table$lowerPR[1], pred$pi.lb[1], tolerance = 0.01)
  expect_equal(orc$mod_table$upperPR[1], pred$pi.ub[1], tolerance = 0.01)
})

test_that("PI uses weighted average for HCS model (vector tau2), no zigzag", {
  skip_on_cran()
  library(metafor)
  data("dat.konstantopoulos2011")
  df <- dat.konstantopoulos2011

  m <- suppressWarnings(rma.mv(yi ~ year, vi,
              random = list(~year|study, ~1|district, ~1|school),
              struct = "HCS", data = df, method = "REML"))
  expect_true(length(m$tau2) > 1)  # HCS: multiple tau2 values

  orc <- suppressWarnings(mod_results(m, mod = "year", group = "study"))

  # Key test: PI width should vary smoothly (due to SE changes),
  # NOT zigzag due to vector recycling
  pi_upper <- orc$mod_table$upperPR - orc$mod_table$estimate
  diffs <- diff(pi_upper)
  # With vector recycling, adjacent diffs would oscillate wildly.
  # With weighted average, the variation is smooth (only SE changes).
  # Check that the range of PI widths is small (< 20% of mean width)
  expect_lt(diff(range(pi_upper)) / mean(pi_upper), 0.2)

  # Also check that all PI values are finite (no NA, no Inf)
  expect_true(all(is.finite(orc$mod_table$lowerPR)))
  expect_true(all(is.finite(orc$mod_table$upperPR)))
})

test_that("PI is scalar (not vector-recycled) for each prediction row", {
  skip_on_cran()
  library(metafor)
  data("dat.konstantopoulos2011")
  df <- dat.konstantopoulos2011

  m <- suppressWarnings(rma.mv(yi ~ year, vi,
              random = list(~year|study, ~1|district, ~1|school),
              struct = "HCS", data = df, method = "REML"))
  orc <- suppressWarnings(mod_results(m, mod = "year", group = "study"))

  # The weighted average tau_val should be a scalar, so the PI half-width

  # (upper - estimate) should differ from row to row only because of SE,
  # not because of different tau2 values being recycled.
  pi_half <- orc$mod_table$upperPR - orc$mod_table$estimate
  se_vals <- orc$mod_table$upperCL - orc$mod_table$estimate  # CI half-width from SE only
  # PI variation should be proportional to SE variation
  expect_equal(
    cor(pi_half, se_vals),
    1.0,
    tolerance = 0.01
  )
})

test_that("PI uses h.levels.k (not g.levels.k) for gamma2 weights", {
  skip_on_cran()
  library(metafor)

  # Directly test the weighted_var function with h.levels.k vs g.levels.k
  # by constructing a mock model and verifying the PI formula components.
  data("dat.konstantopoulos2011")
  df <- dat.konstantopoulos2011

  m <- tryCatch(
    suppressWarnings(rma.mv(yi ~ year, vi,
              random = list(~year|study, ~year|district, ~1|school),
              struct = c("HCS", "HCS"), data = df, method = "REML",
              control = list(iter.max = 2000, eval.max = 4000))),
    error = function(e) NULL
  )
  skip_if(is.null(m), "Model did not converge on this platform")

  # Only run if model has multiple gamma2 values
  skip_if(length(m$gamma2) <= 1 || is.null(m$h.levels.k),
          "Model does not have multiple gamma2 values")

  expected_gamma <- sum(m$gamma2 * m$h.levels.k) / sum(m$h.levels.k)
  wrong_gamma <- sum(m$gamma2 * m$g.levels.k) / sum(m$g.levels.k)

  # The values must differ for this test to be meaningful
  skip_if(abs(expected_gamma - wrong_gamma) < 1e-10,
          "gamma weights are identical — cannot distinguish")

  orc <- suppressWarnings(mod_results(m, mod = "year", group = "study"))

  # Extract effective gamma_val from PI by back-computing:
  # PI_half = t * sqrt(SE^2 + sigmas + tau_val + gamma_val)
  # So: gamma_val = (PI_half / t)^2 - SE^2 - sigmas - tau_val
  tmp <- orc$mod_table
  # Use an approximate t-stat (large df)
  t_stat <- stats::qt(0.975, 100)
  pi_half <- (tmp$upperPR[1] - tmp$estimate[1])
  # Back-compute SE from CI: CI_half = t * SE
  ci_half <- tmp$upperCL[1] - tmp$estimate[1]
  se2 <- (ci_half / t_stat)^2
  sigmas <- sum(m$sigma2)
  tau_val <- sum(m$tau2 * m$g.levels.k) / sum(m$g.levels.k)

  recovered_gamma <- (pi_half / t_stat)^2 - se2 - sigmas - tau_val

  # Should be closer to expected (h.levels.k) than wrong (g.levels.k)
  dist_correct <- abs(recovered_gamma - expected_gamma)
  dist_wrong <- abs(recovered_gamma - wrong_gamma)
  expect_lt(dist_correct, dist_wrong)
})

test_that("PI handles rma.uni models (gamma2 is NULL)", {
  skip_on_cran()
  library(metafor)
  data("dat.konstantopoulos2011")
  df <- dat.konstantopoulos2011

  m <- rma(yi, vi, data = df)
  expect_null(m$gamma2)

  orc <- mod_results(m, mod = "1", group = "study")
  pred <- predict(m)

  expect_equal(as.numeric(orc$mod_table$lowerPR), as.numeric(pred$pi.lb), tolerance = 0.01)
  expect_equal(as.numeric(orc$mod_table$upperPR), as.numeric(pred$pi.ub), tolerance = 0.01)
})

test_that("PI handles rma.mv with no G/H (only sigma2 terms)", {
  skip_on_cran()
  library(metafor)
  data("dat.konstantopoulos2011")
  df <- dat.konstantopoulos2011

  m <- rma.mv(yi ~ 1, vi, random = list(~1|study, ~1|district), data = df)
  # This model has sigma2 but tau2 = NULL, gamma2 = NULL
  # (tau2 and gamma2 only exist for ~inner|outer terms)

  orc <- mod_results(m, mod = "1", group = "study")

  pred <- predict(m)
  expect_equal(orc$mod_table$lowerPR, pred$pi.lb[1], tolerance = 0.01)
  expect_equal(orc$mod_table$upperPR, pred$pi.ub[1], tolerance = 0.01)
})
