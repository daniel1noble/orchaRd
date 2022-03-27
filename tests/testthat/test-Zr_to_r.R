context("Checking Zr_to_r function...")

# Random correlations, r
set.seed(87)
r <- runif(10, min =-1, max = 1)

# Transform to Zr
Zr <- data.frame(Zr = 0.5*(log( (1+r) / (1-r))))

testthat::test_that("Checking Zr_to_r output ...", {
  testthat::expect_equal(
    as.vector(Zr_to_r(Zr)), r,
    info = "Check that Zr is correctly converted back to r...")
})
