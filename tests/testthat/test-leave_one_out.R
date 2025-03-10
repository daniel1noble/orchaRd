# Maybe move this to orchaRd/tests/testthat/setup.R so it can be sourced for all test

library(metafor) 

data(eklof)
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
eklof$paper_ID <- with(eklof, paste(First.author, Publication.year, sep = ", "))

eklof_res <- rma.mv(yi = lnRR,
                    V = vlnRR,
                    mods = ~Grazer.type,
                    random = list(~1 | ExptID,
                                  ~1 | Datapoint),
                    data = eklof)


test_that(".run_leave1out is effectively leaving-out cases and running the expected models", {
  # Creates a model leaving out the studies manually
  # Then compares the parameters with the ones returned by
  # .run_leave1out()

  test_outputs <- .run_leave1out(eklof_res, group = "paper_ID")

  unique_ids <- unique(eklof[["paper_ID"]])

  for (i in seq_along(unique_ids)) { 
    test_model <- test_outputs[[i]]

    tmp_dat <- subset(eklof, paper_ID != unique_ids[i])
    true_model <- metafor::rma.mv(yi     = lnRR,
                                  V      = vlnRR,
                                  mods   = ~ Grazer.type,
                                  random = list(~1 | ExptID,
                                                ~1 | Datapoint),
                                  data   = tmp_dat)

    expect_identical(true_model$k, test_model$k)
    expect_identical(true_model$beta, test_model$beta)
    expect_identical(true_model$ci.ub, test_model$ci.ub)
    expect_identical(true_model$ci.lb, test_model$ci.lb)
  }
})


test_that(".run_leave1out leaves out the correct number of observations", {
  # Create data with known two grouping variables: 'study_ID' and 'species'
  set.seed(123)

  studies_ids <- c(1, 2, 3)
  studies_obs <- c(2, 5, 6)
  id_col <- rep(studies_ids, times = studies_obs)

  species_names <- c("Foo ficticium", "Bar magnificum")
  species_obs <- c(8, 5)
  spp_col <- rep(species_names, times = species_obs)

  mock_data <- data.frame(study_ID = id_col,
                          species = spp_col,
                          yi = rnorm(length(id_col)),
                          vi = abs(rnorm(length(id_col))))

  mock_model <- rma.mv(yi, vi, data = mock_data)

  # Test grouping by 'study_ID'
  # Data has 3 studies with 2, 5 and 6 observations each.
  test_results_study <- .run_leave1out(mock_model, group = "study_ID")
  expect_equal(test_results_study[[1]]$k, (5 + 6))
  expect_equal(test_results_study[[2]]$k, (2 + 6))
  expect_equal(test_results_study[[3]]$k, (2 + 5))

  # Test grouping by 'species'
  # Data has 2 species with 8 and 5 observations
  test_results_spp <- .run_leave1out(mock_model, group = "species")
  expect_equal(test_results_spp[[1]]$k, 5)
  expect_equal(test_results_spp[[2]]$k, 8)
})


test_that(".get_estimates get the correct values", {
  # Run toy model and a leave-one-out using 'species'
    
  set.seed(123)

  studies_ids <- c(1, 2, 3)
  studies_obs <- c(2, 5, 6)
  id_col <- rep(studies_ids, times = studies_obs)

  species_names <- c("Foo ficticium", "Bar magnificum")
  species_obs <- c(8, 5)
  spp_col <- rep(species_names, times = species_obs)

  mock_data <- data.frame(study_ID = id_col,
                          species = spp_col,
                          yi = rnorm(length(id_col)),
                          vi = abs(rnorm(length(id_col))))
  mock_model <- rma.mv(yi, vi, data = mock_data)

  # Results from .get_estimates
  outputs <- .run_leave1out(mock_model, group = "species")
  test_mod_tables <- .get_estimates(outputs, group = "species")

  ## Results from mod_results manually excluding species
  no_foo <- subset(mock_data, species != "Foo ficticium")
  no_foo_mod <- rma.mv(yi,
                       vi,
                       data = no_foo)

  no_bar <- subset(mock_data, species != "Bar magnificum")
  no_bar_mod <- rma.mv(yi,
                       vi,
                       data = no_bar)

  true_mod_tables <- rbind(mod_results(no_foo_mod, group = "species")$mod_table,
                           mod_results(no_bar_mod, group = "species")$mod_table)

  # Don't use the first column because .get_estimates() put species names 
  # while mod_results() only puts 'Intrcpt'.
  expect_equal(test_mod_tables[, -1], true_mod_tables[, -1])
})

test_that(".get_effectsizes gets the correct data", {
    set.seed(123)

    studies_ids <- c(1, 2, 3)
    studies_obs <- c(2, 5, 6)
    id_col <- rep(studies_ids, times = studies_obs)

    species_names <- c("Foo ficticium", "Bar magnificum")
    species_obs <- c(8, 5)
    spp_col <- rep(species_names, times = species_obs)

    mock_data <- data.frame(study_ID = id_col,
                            species = spp_col,
                            yi = rnorm(length(id_col)),
                            vi = abs(rnorm(length(id_col))))

    mock_model <- rma.mv(yi, vi, data = mock_data)
    leave1out_models <- .run_leave1out(mock_model, group = "species")

    test_output <- .get_effectsizes(leave1out_models, group = "species")

    ## Create mod_results manually excluding "Bar magnificum"
    no_bar <- subset(mock_data, species != "Bar magnificum")
    no_bar_mod <- rma.mv(yi, vi, data = no_bar)

    true_effectsizes <- mod_results(no_bar_mod, group = "species")$data

    # Get result from .get_effectsizes() for the model
    # without 'Bar magnificum'.
    # Subset the output using 'moderator', because for 
    # leave-one-out this column has the name of the element
    # left-out.
    test_effectsizes <- subset(test_output, moderator == "Bar magnificum")

    ## Use c() to remove attributes attached
    expect_equal(c(true_effectsizes$yi), test_effectsizes$yi)
    expect_equal(c(true_effectsizes$vi), test_effectsizes$vi)
})


test_that("leave_one_out gives the same result that metafor::leave1out", {
  # Create mock data for the test
  mock_data <- data.frame(study = paste0("Study_", 1:10),
                          yi = c(0.2, -0.1, 0.4, 0.3, 0.5, -0.2, 0.1, 0.6, -0.3, 0.25),
                          vi = c(0.05, 0.07, 0.06, 0.08, 0.04, 0.09, 0.03, 0.05, 0.1, 0.06))

  res <- metafor::rma(yi, vi, data = mock_data)

  # Leave One Out:
  loo_metafor <- metafor::leave1out(res)
  loo_orchard <- leave_one_out(res, group = "study")   
  
  # Metafor output has a lot of column, we only use the estimate, ci.lb and ci.ub
  expect_equal(round(loo_metafor$estimate, 4), round(loo_orchard$mod_table$estimate, 4))
  expect_equal(round(loo_metafor$ci.lb, 4), round(loo_orchard$mod_table$lowerCL, 4))
  expect_equal(round(loo_metafor$ci.ub, 4), round(loo_orchard$mod_table$upperCL, 4))
})


test_that("leave_one_out accepts numeric ID as group and works like metafor::leave1out", {
  # Here 'group' is a numeric variable with ID numbers.
  mock_data2 <- data.frame(study_id = 1:10,
                          yi = c(0.2, -0.1, 0.4, 0.3, 0.5, -0.2, 0.1, 0.6, -0.3, 0.25),
                          vi = c(0.05, 0.07, 0.06, 0.08, 0.04, 0.09, 0.03, 0.05, 0.1, 0.06))

  res <- metafor::rma(yi, vi, data = mock_data2)

  # Leave One Out:
  loo_metafor <- metafor::leave1out(res)
  loo_orchard <- leave_one_out(res, group = "study_id")   
  
  # Metafor output has a lot of column, we only use the estimate, ci.lb and ci.ub
  expect_equal(round(loo_metafor$estimate, 4), round(loo_orchard$mod_table$estimate, 4))
  expect_equal(round(loo_metafor$ci.lb,    4), round(loo_orchard$mod_table$lowerCL,  4))
  expect_equal(round(loo_metafor$ci.ub,    4), round(loo_orchard$mod_table$upperCL,  4))
})
