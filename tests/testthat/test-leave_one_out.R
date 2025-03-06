test_that("leave_one_out gives estimates and CI like a for loop", {
  res <- metafor::rma.mv(lnrr, lnrr_vi,
                         random = list(~ 1 | es_ID,
                                       ~ 1 | paper_ID),
                         data = fish)

  papers_ids <- unique(fish$paper_ID)

  # Initialice data frame to save the results from the models
  loo_loop <- data.frame(left_out = papers_ids,
                         estimate = NA,
                         lowerCL  = NA,
                         upperCL  = NA)

  for (i in seq_along(papers_ids)) {
    # Subset the data an rerun the model
    paper_id <- papers_ids[i]
    fish_sub <- fish[fish$paper_ID != paper_id, ]
    res_sub <- metafor::rma.mv(lnrr, lnrr_vi,
                               random = list(~ 1 | es_ID,
                                             ~ 1 | paper_ID),
                               data = fish_sub)

    loo_loop$estimate[i] <- res_sub$b[1]
    loo_loop$lowerCL[i]  <- res_sub$ci.lb[1]
    loo_loop$upperCL[i]  <- res_sub$ci.ub[1]
  }

  loo_orchard <- leave_one_out(res, group = "paper_ID")

  expect_equal(round(loo_loop$estimate, 4), round(loo_orchard$mod_table$estimate, 4))
  expect_equal(round(loo_loop$lowerCL, 4), round(loo_orchard$mod_table$lowerCL, 4))
  expect_equal(round(loo_loop$upperCL, 4), round(loo_orchard$mod_table$upperCL, 4))
})


test_that("our leave_one_out gives the same result that metafor::leave1out", {
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
