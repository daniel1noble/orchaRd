test_that("leave_one_out gives the same result as a for loop", {
  res <- metafor::rma.mv(lnrr, lnrr_vi,
                         random = list(~ 1 | es_ID,
                                       ~ 1 | paper_ID),
                         data = fish)

  papers_ids <- unique(fish$paper_ID)

  # Initialice data frame to save the results from the models
  loo_loop <- data.frame(left_out = papers_ids,
                        b        = NA,
                        ci_lb    = NA,
                        ci_ub    = NA)

  for (i in seq_along(papers_ids)) {
    # Subset the data an rerun the model
    paper_id <- papers_ids[i]
    fish_sub <- fish[fish$paper_ID != paper_id, ]
    res_sub <- metafor::rma.mv(lnrr, lnrr_vi,
                               random = list(~ 1 | es_ID,
                                             ~ 1 | paper_ID),
                               data = fish_sub)

    loo_loop$b[i]     <- res_sub$b[1]
    loo_loop$ci_lb[i] <- res_sub$ci.lb[1]
    loo_loop$ci_ub[i] <- res_sub$ci.ub[1]
  }

  loo_mine <- leave_one_out(res, group = "paper_ID")

  expect_equal(names(loo_loop), names(loo_mine))
  expect_equal(loo_loop, loo_mine)
})


test_that("our leave_one_out gives the same result that metafor::leave1out", {
  # Create mock data for the test
  mock_data <- data.frame(study = paste0("Study_", 1:10),
                          yi = c(0.2, -0.1, 0.4, 0.3, 0.5, -0.2, 0.1, 0.6, -0.3, 0.25),
                          vi = c(0.05, 0.07, 0.06, 0.08, 0.04, 0.09, 0.03, 0.05, 0.1, 0.06))

  res <- metafor::rma(yi, vi, data = mock_data)

  # Leave One Out:
  loo_metafor <- metafor::leave1out(res)
  loo_mine <- leave_one_out(res, group = "study")   
  
  # Metafor output has a lot of column, we only use the estimate, ci.lb and ci.ub
  expect_equal(
               round(cbind(loo_metafor$estimate, loo_metafor$ci.lb, loo_metafor$ci.ub), 4),
               round(cbind(loo_mine$b, loo_mine$ci_lb, loo_mine$ci_ub), 4)
               )
})
