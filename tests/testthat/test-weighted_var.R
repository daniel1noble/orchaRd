context("Checking weighted_var function...")

x <-  c(3, 5, 6, 7)
n <- c(20, 40, 56, 6)

baseFunc <- weighted.mean(x, n)

    orch <- orchaRd::weighted_var(x, n)

    testthat::test_that("Checking weighted_var output for is being calculated correctly ..", {

      testthat::expect_equal(
        baseFunc, orch,
        info = "Check that weighted_var output matches weighted mean...")
    })
