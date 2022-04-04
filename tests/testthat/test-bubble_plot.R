context("Checking bubble_plot function...")

options(warn=-1)
data(lim)
lim[, "year"] <- as.numeric(lim$year)
lim$vi<- 1/(lim$N - 3)
model<-metafor::rma.mv(yi=yi, V=vi, mods= ~Environment*year,
                       random=list(~1|Article,~1|Datapoint), data=na.omit(lim))

options(warn=0)
test <- orchaRd::mod_results(model, mod = "year", group = "Article",
                    data = lim, weights = "prop", by = "Environment")
plot1 <- orchaRd::bubble_plot(test, mod = "year", legend.pos = "top.left")
plot2 <- orchaRd::bubble_plot(test, mod = "year", xlab = "Year", legend.pos = "top.left")


testthat::test_that("Checking bubble_plot output ...", {

  testthat::expect_equal(
    is.null(plot1$facet), FALSE,
    info = "bubble_plot correctly facetting...")

  testthat::expect_equal(
    plot2$labels$x, "Year",
    info = "Check label of x-axis is correctly renamed...")

  testthat::expect_equal(
    plot1$labels$fill, "condition",
    info = "Check that condition variable is present in fill...")

  testthat::expect_error(orchaRd::bubble_plot(test, mod = "year", legend.pos = "top.left", angle = 45))

})
