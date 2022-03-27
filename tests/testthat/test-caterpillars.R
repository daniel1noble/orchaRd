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

plot1 <- caterpillars(eklof_MR, xlab = "Grazing", group = "Grazer.type", data=eklof)


testthat::test_that("Checking caterpillars output ...", {

  testthat::expect_equal(
    plot1$labels$y, "",
    info = "ylab for eklof is empty...")

  testthat::expect_equal(
    plot1$labels$size, NULL,
    info = "Check precision label doesn't exist in caterpillars plot...")

  testthat::expect_error(caterpillars(eklof_MR, xlab = "Grazing", group = "Grazer.type", data=eklof, angle = 45))

  testthat::expect_error(caterpillars(eklof_MR, xlab = "Grazing", group = "Grazer.type"))

  testthat::expect_error(caterpillars(eklof_MR, xlab = "Grazing", data=eklof))

})

