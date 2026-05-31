context("Checking catepillars function...")

data(eklof)
eklof<-metafor::escalc(measure="ROM", n1i=N_control, sd1i=SD_control,
                       m1i=mean_control, n2i=N_treatment, sd2i=SD_treatment, m2i=mean_treatment,
                       data=eklof)
# Add the unit level predictor
eklof$Datapoint<-as.factor(seq(1, dim(eklof)[1], 1))
# fit a MLMR - accouting for some non-independence
eklof_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~ Grazer.type, random=list(~1|ExptID,
                                                                       ~1|Datapoint), data=eklof)

plot1 <- caterpillars(eklof_MR, xlab = "Grazing", group = "Grazer.type")


testthat::test_that("Checking caterpillars output ...", {

  testthat::expect_equal(
    plot1$labels$y, "",
    info = "ylab for eklof is empty...")

  testthat::expect_equal(
    plot1$labels$size, NULL,
    info = "Check precision label doesn't exist in caterpillars plot...")

  testthat::expect_error(caterpillars(eklof_MR, xlab = "Grazing", group = "Grazer.type", angle = 45))

  testthat::expect_error(caterpillars(eklof_MR, xlab = "Grazing"))

})


# Regression test for the tanh back-transformation bug: when correlations are
# analysed as Fisher's z and back-transformed with transfm = "tanh", the
# per-effect-size confidence intervals must be computed on the z scale and then
# transformed. Previously the z-scale variance was added to the already
# tanh-transformed estimate, which mixed scales and produced correlations
# outside [-1, 1].
data(lim)
lim$vi <- 1 / (lim$N - 3)
lim_MR <- metafor::rma.mv(yi = yi, V = vi, mods = ~Phylum - 1,
                          random = list(~1 | Article, ~1 | Datapoint),
                          data = lim)
results_lim <- mod_results(lim_MR, mod = "Phylum", group = "Article")

testthat::test_that("caterpillars tanh transform keeps correlations within [-1, 1]", {

  p <- caterpillars(results_lim, mod = "Phylum", group = "Article",
                    xlab = "Correlation coefficient", transfm = "tanh")
  d <- p$data

  # Point estimates and both CI bounds must all lie within the valid range.
  testthat::expect_true(all(d$yi >= -1 & d$yi <= 1),
                        info = "Transformed point estimates exceed [-1, 1]")
  testthat::expect_true(all(d$lower >= -1 & d$lower <= 1),
                        info = "Lower CI bounds exceed [-1, 1]")
  testthat::expect_true(all(d$upper >= -1 & d$upper <= 1),
                        info = "Upper CI bounds exceed [-1, 1]")

  # The transform must work the same way whether a model object or an orchard
  # object of results is supplied.
  p_mod <- caterpillars(lim_MR, mod = "Phylum", group = "Article",
                        xlab = "r", transfm = "tanh")
  testthat::expect_true(all(p_mod$data$lower >= -1 & p_mod$data$upper <= 1),
                        info = "Model-object path produces out-of-range bounds")
})

testthat::test_that("caterpillars tanh CI bounds equal tanh of z-scale interval", {

  # Independently reconstruct the expected bounds: build the CI on the raw
  # (Fisher's z) scale, then apply tanh. This is the statistically correct
  # ordering and what the plot should now use.
  raw <- mod_results(lim_MR, mod = "Phylum", group = "Article")$data
  expected_lower <- tanh(raw$yi - stats::qnorm(0.975) * sqrt(raw$vi))
  expected_upper <- tanh(raw$yi + stats::qnorm(0.975) * sqrt(raw$vi))

  p <- caterpillars(results_lim, mod = "Phylum", group = "Article",
                    xlab = "Correlation coefficient", transfm = "tanh")
  # caterpillars reorders rows within moderator, so compare on sorted values.
  testthat::expect_equal(sort(p$data$lower), sort(expected_lower))
  testthat::expect_equal(sort(p$data$upper), sort(expected_upper))
})

testthat::test_that("caterpillars with transfm = 'none' returns raw-scale CIs", {

  p <- caterpillars(lim_MR, mod = "Phylum", group = "Article", xlab = "Zr")
  d <- p$data

  expected_lower <- d$yi - stats::qnorm(0.975) * sqrt(d$vi)
  expected_upper <- d$yi + stats::qnorm(0.975) * sqrt(d$vi)
  testthat::expect_equal(d$lower, expected_lower)
  testthat::expect_equal(d$upper, expected_upper)
})

