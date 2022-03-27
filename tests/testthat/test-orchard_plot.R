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

plot1 <- orchard_plot(eklof_MR, xlab = "Grazing", group = "Grazer.type", data=eklof)

plot2 <- orchard_plot(eklof_MR, xlab = "Grazing", group = "Grazer.type", data=eklof, angle = 45)


testthat::test_that("Checking orchard_plot output ...", {

  testthat::expect_equal(
    plot1$labels$y, "Grazing",
    info = "xlab for eklof is correct...")

  testthat::expect_equal(
    plot1$labels$size, "Precision (1/SE)",
    info = "Check precision label for eklof is correct...")

  testthat::expect_equal(
    plot1$labels$ymin, "lowerCL",
    info = "Check that ymin is 'lowerCL' for eklof...")

  testthat::expect_equal(
    plot1$theme$axis.text.y$angle, 90,
    info = "Check that angle for eklof is defaulting to 90...")

  testthat::expect_equal(
    plot2$theme$axis.text.y$angle, 45,
    info = "Check that angle for eklof is changes to 45 when done...")

  testthat::expect_error(orchard_plot(eklof_MR, xlab = "Grazing", group = "Grazer.type"))

  testthat::expect_error(orchard_plot(eklof_MR, xlab = "Grazing", data=eklof))

})

