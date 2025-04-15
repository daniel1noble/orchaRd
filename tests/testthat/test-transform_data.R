test_that("transform_data applies the correct transformations", {
  x <- c(-1, 0, 1)
  
  # No transformation
  expect_equal(transform_data(x, transfm ="none"), x)
  
  # Hyperbolic tangent transformation
  expect_equal(transform_data(x, transfm ="tanh"), tanh(x))
  
  # Inverse logit transformation
  expect_equal(transform_data(x, transfm ="invlogit"), exp(x) / (1 + exp(x)))
  
  # Percentage relative change transformation
  expect_equal(transform_data(x, transfm ="percentr"), (exp(x) - 1) * 100)
  
  # Percentage transformation
  expect_equal(transform_data(x, transfm ="percent"), exp(x) * 100)
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
  new_method_table[, numeric_cols] <- transform_data(mod_table[, numeric_cols], transfm ="invlogit")
  expect_equal(old_method_table, new_method_table)
  expect_equal(metafor::transf.ilogit(data_trim$yi),
               transform_data(data_trim$yi, transfm ="invlogit"))


  # Hyperbolic tangent transformation
  old_method_table[, numeric_cols] <- lapply(mod_table[, numeric_cols], function(x) tanh(x))
  new_method_table[, numeric_cols] <- transform_data(mod_table[, numeric_cols], transfm ="tanh")
  expect_equal(old_method_table, new_method_table)
  expect_equal(tanh(data_trim$yi),
               transform_data(data_trim$yi, transfm ="tanh"))

  # Percentage relative change transformation
  old_method_table[, numeric_cols] <- lapply(mod_table[, numeric_cols], function(x) (exp(x) - 1) * 100)
  new_method_table[, numeric_cols] <- transform_data(mod_table[, numeric_cols], transfm ="percentr")
  expect_equal(old_method_table, new_method_table)
  expect_equal((exp(data_trim$yi) - 1) * 100,
               transform_data(data_trim$yi, transfm ="percentr"))

  # Percentage transformation
  old_method_table[, numeric_cols] <- lapply(mod_table[, numeric_cols], function(x) exp(x) * 100)
  new_method_table[, numeric_cols] <- transform_data(mod_table[, numeric_cols], transfm ="percent")
  expect_equal(old_method_table, new_method_table)
  expect_equal(exp(data_trim$yi) * 100,
               transform_data(data_trim$yi, transfm ="percent"))
})

## Checking the inverse transformation of Freeman Tukey
datap <- data.frame(  
  StudyLabel = c("Study A", "Study B", "Study C", "Study D", "Study E",  
                 "Study F", "Study G", "Study H", "Study I", "Study J"), 
  Infected = c(20, 35, 15, 8, 50, 90, 22, 30, 10, 190), #Infected STDs Cases
  TotalSample = c(100, 200, 1000, 100, 300, 1500, 120, 180, 80, 2100), #Total Study Sample
  PublicationYear = c(2013, 2015, 2016, 2016, 2018, 2019, 2021, 2021, 2023, 2023) #Study Publication Year
)  
datap_rma <- escalc(measure = "PFT",   
                   xi = datap$Infected,   
                   ni = datap$TotalSample,   
                   data = datap,   
                   slab = datap$StudyLabel)

# Check transformation works on the raw data, gives same results as metafor.
 expect_equal(transf_inv_ft(datap_rma$yi, datap_rma$TotalSample), metafor::transf.ipft.hm(datap_rma$yi, targs = list(ni = datap_rma$TotalSample)))


## TO DO: Need to test also mod table conversion. Something like this, but note that we would be transforming to median not mean because we are transforming the estimates from the models, so we need to adjust in orchard to get the mean rather than median.:
# Run model
meta_rma <- rma(yi = datap_rma$yi,   
               vi = datap_rma$vi, 
               method = "REML",   
               test = "knha",
               data = datap_rma,  
               verbose = TRUE,
               digits = 3)

meta_rma_new <- predict(meta_rma, transf = transf.ipft.hm, targ = list(ni = datap$TotalSample))  #Back-transformation

mod_table <- mod_results(meta_rma, mod = "1", group = "StudyLabel")$mod_table
mod_table <- transf_inv_ft(as.numeric(mod_table[1,-1]), datap$TotalSample)