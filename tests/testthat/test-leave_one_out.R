# Setup code for all tests
library(metafor)

# Create a common setup function to avoid repetition
setup_eklof_data <- function() {
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
  
  return(list(data = eklof, model = eklof_res))
}

# Create a common setup for mock data to avoid repetition
create_mock_data <- function() {
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
  
  return(mock_data)
}

test_that(".run_leave1out correctly leaves out studies and produces expected model parameters", {
  # Compare the results of .run_leave1out() with the model parameters
  # calculated doing a leave-one-out with a for-loop.

  eklof_setup <- setup_eklof_data()
  eklof <- eklof_setup$data
  eklof_res <- eklof_setup$model

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


test_that(".run_leave1out leaves out the correct number of observations for different grouping variables", {
  # Compare the number of observations left out with .run_leave1out() with the
  # known number of observations left out.
  mock_data <- create_mock_data()
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


test_that("leave_one_out() works correctly with robust.rma models", {
  # Compares the model parameters returned by leave_one_out() using robust
  eklof_setup <- setup_eklof_data()
  eklof_data <- eklof_setup$data

  # -----------------------------------------------------------
  # Get the leave_one_out() results with robust_args
  eklof_model_base <- metafor::rma.mv(yi = lnRR,
                                      V  = vlnRR,
                                      random = list(~ 1 | ExptID,
                                                    ~ 1 | Datapoint),
                                      data = eklof_data)

  test_results <- leave_one_out(eklof_model_base,
                                group = "paper_ID",
                                robust_args = list(cluster = "paper_ID"))

  # -----------------------------------------------------------
  # Create 2 subsets, fit the model and use robust
  no_Nelson <- subset(eklof_data, paper_ID != "Nelson, 1981")
  no_Moksnes <- subset(eklof_data, paper_ID != "Moksnes, 2008")


  no_Nelson_base <- metafor::rma.mv(yi = lnRR,
                                    V  = vlnRR,
                                    random = list(~ 1 | ExptID,
                                                  ~ 1 | Datapoint),
                                    data = no_Nelson)
  no_Moksnes_base <- metafor::rma.mv(yi = lnRR,
                                     V  = vlnRR,
                                     random = list(~ 1 | ExptID,
                                                   ~ 1 | Datapoint),
                                     data = no_Moksnes)

  rob_no_Nelson <- robust(no_Nelson_base, cluster = no_Nelson$paper_ID)
  rob_no_Moksnes <- robust(no_Moksnes_base, cluster = no_Moksnes$paper_ID)

  # -----------------------------------------------------------
  # Compare the results

  loo_no_Nelson <- subset(test_results$mod_table, name == "Nelson, 1981") 
  loo_no_Moksnes <- subset(test_results$mod_table, name == "Moksnes, 2008") 

  # First compare with raw model. Confidence intervals must be different
  expect_error(expect_equal(no_Nelson_base$ci.lb, loo_no_Nelson$lowerCL))
  expect_error(expect_equal(no_Nelson_base$ci.ub, loo_no_Nelson$upperCL))

  # Then compare with the robust models
  expect_equal(rob_no_Nelson$beta[1], loo_no_Nelson$estimate)
  expect_equal(rob_no_Nelson$ci.lb, loo_no_Nelson$lowerCL)
  expect_equal(rob_no_Nelson$ci.ub, loo_no_Nelson$upperCL)

  expect_equal(rob_no_Moksnes$beta[1], loo_no_Moksnes$estimate)
  expect_equal(rob_no_Moksnes$ci.lb, loo_no_Moksnes$lowerCL)
  expect_equal(rob_no_Moksnes$ci.ub, loo_no_Moksnes$upperCL)
})


test_that(".get_estimates returns the correct values", {
  mock_data <- create_mock_data()
  mock_model <- rma.mv(yi, vi, data = mock_data)

  # Results from .get_estimates
  outputs <- .run_leave1out(mock_model, group = "species")
  test_mod_tables <- .get_estimates(outputs, group = "species")

  ## Results from mod_results manually excluding species
  no_foo <- subset(mock_data, species != "Foo ficticium")
  no_foo_mod <- rma.mv(yi, vi, data = no_foo)

  no_bar <- subset(mock_data, species != "Bar magnificum")
  no_bar_mod <- rma.mv(yi, vi, data = no_bar)

  true_mod_tables <- rbind(mod_results(no_foo_mod, group = "species")$mod_table,
                           mod_results(no_bar_mod, group = "species")$mod_table)

  # Don't use the first column because .get_estimates() put species names 
  # while mod_results() only puts 'Intrcpt'.
  expect_equal(test_mod_tables[, -1], true_mod_tables[, -1])
})


test_that(".get_effectsizes returns the correct data structure", {
  mock_data <- create_mock_data()
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


test_that("leave_one_out produces identical results to metafor::leave1out for simple models", {
  # Create mock data for the test
  mock_data <- data.frame(study = paste0("Study_", 1:10),
                          yi = c(0.2, -0.1, 0.4, 0.3, 0.5, -0.2, 0.1, 0.6, -0.3, 0.25),
                          vi = c(0.05, 0.07, 0.06, 0.08, 0.04, 0.09, 0.03, 0.05, 0.1, 0.06))

  res <- metafor::rma(yi, vi, data = mock_data)

  # Leave One Out:
  loo_metafor <- metafor::leave1out(res)
  loo_orchard <- leave_one_out(res, group = "study")   
  
  # Metafor output has a lot of columns, we only use the estimate, ci.lb and ci.ub
  expect_equal(round(loo_metafor$estimate, 4), round(loo_orchard$mod_table$estimate, 4))
  expect_equal(round(loo_metafor$ci.lb, 4), round(loo_orchard$mod_table$lowerCL, 4))
  expect_equal(round(loo_metafor$ci.ub, 4), round(loo_orchard$mod_table$upperCL, 4))
})


test_that(".prune_tree correctly removes branches from a phylogenetic tree", {
  # Create a tree with 4 species
  full_spp_names <- c("Festuca_arundinacea",
                      "Lolium_perenne",
                      "Solanum_tuberosum",
                      "Nassella_neesiana")

  # Read phylo tree created with those species
  tree <- ape::read.tree(test_path("dummy_tree.rda"))

  # Remove one species
  removed_species <- "Solanum_tuberosum"
  sliced_spp_names <- full_spp_names[full_spp_names != removed_species]
  pruned_tree <- .prune_tree(tree, sliced_spp_names)

  expect_s3_class(pruned_tree, "phylo")
  expect_false(removed_species %in% pruned_tree$tip.label)
  expect_equal(length(pruned_tree$tip.label), length(sliced_spp_names))
})


test_that(".create_tmp_phylo_matrix creates the correct phylogenetic variance-covariance matrix", {
  full_spp_names <- c("Festuca_arundinacea",
                      "Lolium_perenne",
                      "Solanum_tuberosum",
                      "Nassella_neesiana")
  # Read phylo tree created with those species
  tree <- ape::read.tree(test_path("dummy_tree.rda"))
  
  # Compute phylo matrix for full dataset
  tree <- ape::compute.brlen(tree)
  full_phylo_matrix <- ape::vcv(tree, corr = TRUE)
  
  # --------------------------------------------------------------
  # Remove one species and compute the phylo matrix manually
  short_spp_names <- c("Festuca_arundinacea",
                       "Lolium_perenne",
                       "Nassella_neesiana")
  # Get tree for the reduced dataset 
  short_tree <- .prune_tree(tree, short_spp_names)

  # Compute phylo matrix
  short_tree <- ape::compute.brlen(short_tree)
  short_phylo_matrix <- ape::vcv(short_tree, corr = TRUE)
  
  # --------------------------------------------------------------
  # Test .create_tmp_phylo_matrix
  # It must be: 
  #   - Different than full_phylo_matrix
  #   - Equal to short_phylo_matrix
  dat <- data.frame(spp_names = short_spp_names)
  phylo_args <- list(tree = tree, species_colname = "spp_names")
  test_matrix <- .create_tmp_phylo_matrix(dat, phylo_args)
  
  expect_false(identical(full_phylo_matrix, test_matrix))
  expect_true(identical(short_phylo_matrix, test_matrix))
})


test_that(".run_leave1out correctly handles phylogenetic arguments", {
  # Create mock data using four species
  full_spp_names <- c("Festuca_arundinacea",
                      "Lolium_perenne",
                      "Solanum_tuberosum",
                      "Nassella_neesiana")

  # Read phylo tree created with those species
  tree <- ape::read.tree(test_path("dummy_tree.rda"))

  # Create mock data
  # Create 5 observations for each species
  set.seed(123)
  mock_data <- data.frame(spp_names = rep(full_spp_names, times = 5),
                          yi = rnorm(20),
                          vi = abs(rnorm(20)))

  # Add new column to link species with phylo matrix
  mock_data$phylo <- mock_data$spp_names

  # Compute the full matrix
  tree <- ape::compute.brlen(tree)
  phylo_matrix <- ape::vcv(tree, corr = TRUE)

  # --------------------------------------------------------------
  # Create short data set without "Solanum_tuberosum"
  short_spp_names <- c("Festuca_arundinacea",
                       "Lolium_perenne",
                       "Nassella_neesiana")
  short_tree <- .prune_tree(tree, short_spp_names)
  short_tree <- ape::compute.brlen(short_tree)
  short_phylo_matrix <- ape::vcv(short_tree, corr = TRUE)

  short_data <- subset(mock_data, spp_names != "Solanum_tuberosum")

  true_results <- metafor::rma.mv(yi, vi,
                                  random = list(~ 1 | spp_names,
                                                ~ 1 | phylo),
                                  R = list(phylo = short_phylo_matrix),
                                  data = short_data)

  # --------------------------------------------------------------
  # Create the full model to leave-one-out
  mock_model <- metafor::rma.mv(yi, vi,
                         random = list(~ 1 | spp_names,
                                       ~ 1 | phylo),
                         R = list(phylo = phylo_matrix),
                         data = mock_data)

  # Run leave-one-out
  loo_results <- leave_one_out(mock_model,
                               group = "spp_names",
                               phylo_args = list(tree = tree, species_colname = "phylo"))

  test_results <- subset(loo_results$mod_table, name == "Solanum_tuberosum")

  # --------------------------------------------------------------
  expect_equal(c(true_results$beta), test_results$estimate, tolerance = 1e-4)
  expect_equal(c(true_results$ci.lb), test_results$lowerCL, tolerance = 1e-4)
  expect_equal(c(true_results$ci.ub), test_results$upperCL, tolerance = 1e-4)
})

test_that("leave_one_out() throw errors", {
  mock_data <- create_mock_data()
  mock_model <- rma.mv(yi, vi, data = mock_data)
  
  # Test with non-existent grouping variable
  expect_error(
    leave_one_out(mock_model, group = "nonexistent_group"),
  )
  
  # Test with NULL model
  expect_error(
    leave_one_out(NULL, group = "study_ID")
  )
})
