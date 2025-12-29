# tests/testthat/test-lnM-safe.R

testthat::test_that(".MSb computes correct MS between groups (indep and paired)", {
  x1 <- 10
  x2 <- 7
  n1 <- 20
  n2 <- 30

  # independent
  expect_equal(
    .MSb(x1, x2, n1, n2, paired = FALSE),
    ((n1 * n2) / (n1 + n2)) * (x1 - x2)^2
  )

  # paired uses n1 as n (by design of the function)
  expect_equal(
    .MSb(x1, x2, n1, n2, paired = TRUE),
    (n1 / 2) * (x1 - x2)^2
  )

  # symmetry in swapping groups for independent case
  expect_equal(
    .MSb(x1, x2, n1, n2, paired = FALSE),
    .MSb(x2, x1, n2, n1, paired = FALSE)
  )
})

testthat::test_that(".MSw computes correct MS within groups (indep and paired)", {
  sd1 <- 2
  sd2 <- 3
  n1 <- 20
  n2 <- 30

  # independent (pooled variance)
  expect_equal(
    .MSw(sd1, sd2, n1, n2, paired = FALSE),
    ((n1 - 1) * sd1^2 + (n2 - 1) * sd2^2) / (n1 + n2 - 2)
  )

  # paired
  expect_equal(
    .MSw(sd1, sd2, n1, n2, paired = TRUE),
    (sd1^2 + sd2^2) / 2
  )

  # non-negativity
  expect_gte(.MSw(sd1, sd2, n1, n2, paired = FALSE), 0)
  expect_gte(.MSw(sd1, sd2, n1, n2, paired = TRUE), 0)
})

testthat::test_that(".lnM matches manual computation and is symmetric under swapping groups", {
  x1 <- 10
  x2 <- 7
  sd1 <- 2
  sd2 <- 3
  n1 <- 20
  n2 <- 30

  n0 <- (2 * n1 * n2) / (n1 + n2)
  sw2 <- .MSw(sd1, sd2, n1, n2, paired = FALSE)
  sb2 <- (.MSb(x1, x2, n1, n2, paired = FALSE) - sw2) / n0
  lnM_manual <- log(sqrt(sb2)) - log(sqrt(sw2))

  expect_equal(.lnM(x1, x2, sd1, sd2, n1, n2), lnM_manual)

  # swapping groups should not change lnM because it depends on (x1-x2)^2 and pooled terms
  expect_equal(
    .lnM(x1, x2, sd1, sd2, n1, n2),
    .lnM(x2, x1, sd2, sd1, n2, n1)
  )
})

testthat::test_that(".v_lnM returns finite positive variance in typical cases (independent)", {
  x1 <- 10
  x2 <- 7
  sd1 <- 2
  sd2 <- 3
  n1 <- 40
  n2 <- 50

  v <- .v_lnM(x1, x2, sd1, sd2, n1, n2, r = NULL)

  expect_type(v, "double")
  expect_length(v, 1)
  expect_true(is.finite(v))
  expect_gt(v, 0)
})

testthat::test_that(".v_lnM returns finite positive variance in typical cases (dependent)", {
  x1 <- 10
  x2 <- 7
  sd1 <- 2
  sd2 <- 3
  n  <- 60
  r  <- 0.4

  v <- .v_lnM(x1, x2, sd1, sd2, n, n, r = r)

  expect_type(v, "double")
  expect_length(v, 1)
  expect_true(is.finite(v))
  expect_gt(v, 0)
})
