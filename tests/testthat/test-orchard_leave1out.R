test_that("orchard_leave1out creates a ggplot object", {
  data(eklof)

  # Calculate the effect size
  eklof <- metafor::escalc(measure = "ROM",
                           n1i = N_control,
                           sd1i = SD_control,
                           m1i = mean_control, 
                           n2i = N_treatment,
                           sd2i = SD_treatment,
                           m2i = mean_treatment,
                           var.names = c("lnRR", "vlnRR"),
                           data = eklof)
  
  # Add the observation level factor
  eklof$Datapoint <- as.factor(seq(1, dim(eklof)[1], 1))
  
  # Also, we can get the sample size, which we can use for weighting if we would like
  eklof$N <- rowSums(eklof[, c("N_control", "N_treatment")])

  eklof$Author <- paste(eklof$First.author, eklof$Publication.year, sep = ", ")

  # fit a meta-regression with the intercept (contrast)
  res <- metafor::rma.mv(yi = lnRR,
                         V = vlnRR,
                         mods = ~ Grazer.type,
                         random = list(~1 | ExptID,
                                       ~1 | Datapoint),
                         data = eklof)

  p <- orchard_leave1out(res, group = "Author", xlab = "lnRR")

  expect_no_error(orchard_leave1out(res, group = "Author", xlab = "lnRR"))
  expect_s3_class(p, "ggplot")
})
