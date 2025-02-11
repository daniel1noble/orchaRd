test_that("Checking that leave_one_out gives the same result as a for loop", {
  res <- metafor::rma.mv(lnrr, lnrr_vi,
                         random = list(~ 1 | es_ID,
                                       ~ 1 | paper_ID),
                         data = fish)

  papers_ids <- unique(fish$paper_ID)

  # Initialice data frame to save the results from the models
  res_all <- data.frame(left_out = papers_ids,
                        b        = NA,
                        ci_ub    = NA,
                        ci_lb    = NA)

  for (i in seq_along(papers_ids)) {
    # Subset the data an rerun the model
    paper_id <- papers_ids[i]
    fish_sub <- fish[fish$paper_ID != paper_id, ]
    res_sub <- metafor::rma.mv(lnrr, lnrr_vi,
                               random = list(~ 1 | es_ID,
                                             ~ 1 | paper_ID),
                               data = fish_sub)

    res_all$b[i] <- res_sub$b[1]
    res_all$ci_ub[i] <- res_sub$ci.ub[1]
    res_all$ci_lb[i] <- res_sub$ci.lb[1]
  }

  expect_equal(leave_one_out(res, group = "paper_ID"), res_all)
})
