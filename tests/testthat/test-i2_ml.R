context("Checking i2_ml function..")

data(lim)
lim$vi<-(1/sqrt(lim$N - 3))^2
lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article, ~1|Datapoint), data=lim)
I2_lim_2 <- i2_ml(lim_MR, data=lim)

testthat::test_that("Checking i2 has correct dimensions..", {
  testthat::expect_equal(
    length(I2_lim_2), 3,
    info = "i2 has correct dimensions")

testthat::expect_equal(round(as.vector(I2_lim_2),2), round(c(84.32574, 48.19211, 36.13363), 2), info = "checking i2 estaimtes calculated correctly...")
})