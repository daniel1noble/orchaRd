context("Checking orchard_plot function...")

data(eklof)
eklof<-metafor::escalc(measure="ROM", n1i=N_control, sd1i=SD_control,
                       m1i=mean_control, n2i=N_treatment, sd2i=SD_treatment, m2i=mean_treatment,
                       data=eklof)
# Add the unit level predictor
eklof$Datapoint<-as.factor(seq(1, dim(eklof)[1], 1))
# fit a MLMR - accouting for some non-independence
eklof_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~ Grazer.type, random=list(~1|ExptID,
                                                                       ~1|Datapoint), data=eklof)

plot1 <- orchard_plot(eklof_MR, xlab = "Grazing", group = "Grazer.type")

plot2 <- orchard_plot(eklof_MR, xlab = "Grazing", group = "Grazer.type", angle = 45)


testthat::test_that("Checking orchard_plot output ...", {

  testthat::expect_equal(
    ggplot2::get_labs(plot1)$y, "Grazing",
    info = "xlab for eklof is correct...")

  testthat::expect_equal(
    ggplot2::get_labs(plot1)$ymin, "lowerCL",
    info = "Check that ymin is 'lowerCL' for eklof...")

  testthat::expect_equal(
    plot1$theme$axis.text.y$angle, 90,
    info = "Check that angle for eklof is defaulting to 90...")

  testthat::expect_equal(
    plot2$theme$axis.text.y$angle, 45,
    info = "Check that angle for eklof is changes to 45 when done...")

  testthat::expect_error(orchard_plot(eklof_MR, xlab = "Grazing", data=eklof))

})


testthat::test_that("point.size controls the size scale range in orchard_plot", {
  p_default <- orchard_plot(eklof_MR, xlab = "Grazing", group = "ExptID")
  p_small   <- orchard_plot(eklof_MR, xlab = "Grazing", group = "ExptID",
                            point.size = c(0.5, 2))

  # Extract the mapped size values from the built plot data
  default_sizes <- range(ggplot2::ggplot_build(p_default)$data[[1]]$size)
  small_sizes   <- range(ggplot2::ggplot_build(p_small)$data[[1]]$size)

  # The small range plot should have smaller max size
  testthat::expect_lt(max(small_sizes), max(default_sizes))
})


testthat::test_that("upper=FALSE keeps lowercase moderator names", {
  data(lim)
  lim$Phylum <- tolower(lim$Phylum)
  lim$vi <- (1/sqrt(lim$N - 3))^2
  lim_MR <- metafor::rma.mv(yi = yi, V = vi, mods = ~ Phylum,
    random = list(~1 | Article, ~1 | Datapoint), data = lim)

  p <- orchard_plot(lim_MR, mod = "Phylum", group = "Article",
    xlab = "Zr", transfm = "tanh", N = "N", upper = FALSE)
  built <- ggplot2::ggplot_build(p)

  # Moderator levels in the data should be lowercase
  mod_levels <- levels(built$plot$data$moderator)
  testthat::expect_true(all(mod_levels == tolower(mod_levels)))
})


testthat::test_that("tree.order works with original-case category names", {
  data(lim)
  lim$Phylum <- tolower(lim$Phylum)
  lim$vi <- (1/sqrt(lim$N - 3))^2
  lim_MR <- metafor::rma.mv(yi = yi, V = vi, mods = ~ Phylum,
    random = list(~1 | Article, ~1 | Datapoint), data = lim)

  new_order <- rev(sort(unique(lim_MR$data$Phylum)))

  # Should not error with lowercase tree.order and default upper=TRUE
  testthat::expect_no_error(
    orchard_plot(lim_MR, mod = "Phylum", group = "Article",
      xlab = "Zr", transfm = "tanh", N = "N", tree.order = new_order)
  )

  # Should not error with lowercase tree.order and upper=FALSE
  testthat::expect_no_error(
    orchard_plot(lim_MR, mod = "Phylum", group = "Article",
      xlab = "Zr", transfm = "tanh", N = "N",
      tree.order = new_order, upper = FALSE)
  )
})

