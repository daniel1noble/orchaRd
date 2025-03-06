test_that("transform_data applies the correct transformations", {
  x <- c(-1, 0, 1)
  
  # No transformation
  expect_equal(transform_data(x, "none"), x)
  
  # Hyperbolic tangent transformation
  expect_equal(transform_data(x, "tanh"), tanh(x))
  
  # Inverse logit transformation
  expect_equal(transform_data(x, "invlogit"), exp(x) / (1 + exp(x)))
  
  # Percentage relative change transformation
  expect_equal(transform_data(x, "percentr"), (exp(x) - 1) * 100)
  
  # Percentage transformation
  expect_equal(transform_data(x, "percent"), exp(x) * 100)
})

test_that("internal transformation functions work independently", {
  x <- c(-2, -1, 0, 1, 2)
  
  expect_equal(transf_tanh(x), tanh(x))
  expect_equal(transf_invlogit(x), exp(x) / (1 + exp(x)))
  expect_equal(transf_percentr(x), (exp(x) - 1) * 100)
  expect_equal(transf_percent(x), exp(x) * 100)
})

test_that("transform_data throws an error for invalid transformation", {
  x <- c(-1, 0, 1)
  
  expect_error(transform_data(x, transf = "foo"))
})


test_that("transform_data works identically to how trasnf was handled by orchard_plot", {
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

  eklof$Datapoint <- as.factor(seq(1, dim(eklof)[1], 1))
  eklof$N <- rowSums(eklof[, c("N_control", "N_treatment")])

  # fit a meta-regression with the intercept (contrast)
  eklof_MR0 <- metafor::rma.mv(yi = lnRR,
                               V = vlnRR,
                               mods = ~Grazer.type,
                               random = list(~1 | ExptID,
                                             ~1 | Datapoint),
                               data = eklof)

  results <- mod_results(eklof_MR0, group = "ExptID")

  data_trim <- results$data
  mod_table <- results$mod_table
	numeric_cols <- sapply(mod_table, is.numeric)

  old_method_table <- mod_table
  new_method_table <- mod_table

  # Inverse logit transformation
  old_method_table[, numeric_cols] <- lapply(mod_table[, numeric_cols], function(x) metafor::transf.ilogit(x))
  new_method_table[, numeric_cols] <- transform_data(mod_table[, numeric_cols], "invlogit")
  expect_equal(old_method_table, new_method_table)
  expect_equal(metafor::transf.ilogit(data_trim$yi),
               transform_data(data_trim$yi, "invlogit"))


  # Hyperbolic tangent transformation
  old_method_table[, numeric_cols] <- lapply(mod_table[, numeric_cols], function(x) tanh(x))
  new_method_table[, numeric_cols] <- transform_data(mod_table[, numeric_cols], "tanh")
  expect_equal(old_method_table, new_method_table)
  expect_equal(tanh(data_trim$yi),
               transform_data(data_trim$yi, "tanh"))

  # Percentage relative change transformation
  old_method_table[, numeric_cols] <- lapply(mod_table[, numeric_cols], function(x) (exp(x) - 1) * 100)
  new_method_table[, numeric_cols] <- transform_data(mod_table[, numeric_cols], "percentr")
  expect_equal(old_method_table, new_method_table)
  expect_equal((exp(data_trim$yi) - 1) * 100,
               transform_data(data_trim$yi, "percentr"))

  # Percentage transformation
  old_method_table[, numeric_cols] <- lapply(mod_table[, numeric_cols], function(x) exp(x) * 100)
  new_method_table[, numeric_cols] <- transform_data(mod_table[, numeric_cols], "percent")
  expect_equal(old_method_table, new_method_table)
  expect_equal(exp(data_trim$yi) * 100,
               transform_data(data_trim$yi, "percent"))
})

