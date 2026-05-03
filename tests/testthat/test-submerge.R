context("Checking sub_merge function..")

data(fish)
warm_dat <- fish
model <- metafor::rma.mv(yi = lnrr, V = lnrr_vi,random = list(~1 | group_ID, ~1 | es_ID),
                         mods = ~ experimental_design + trait.type + deg_dif + treat_end_days, method = "REML", test = "t", data = warm_dat, control=list(optimizer="optim", optmethod="Nelder-Mead"))

# Intercept only
overall <- orchaRd::mod_results(model, group = "group_ID", data = warm_dat)

# Moderators
#Trait type
across_trait <- orchaRd::mod_results(model, group = "group_ID", mod = "trait.type", data = warm_dat)

merge <- orchaRd::submerge(across_trait, overall)

testthat::test_that("Checking submerge function for fish dataset ..", {

  testthat::expect_equal(
    dim(merge$mod_table)[[1]], 5,
    info = "Check that merging of 'overall' (1) and 'across_trait' (4) for fish have dimensions pf 5...")

})

# --- mix = TRUE renames moderator levels by appending the source-table index. ---
mix_merge <- orchaRd::submerge(across_trait, overall, mix = TRUE)

testthat::test_that("submerge(mix = TRUE) re-labels moderator levels by source index", {
  testthat::expect_equal(nrow(mix_merge$mod_table), 5)
  # Every moderator label is suffixed with '1' (from across_trait) or '2'
  # (from overall); the last character of the factor label encodes the source.
  labels <- as.character(mix_merge$data$moderator)
  suffixes <- substr(labels, nchar(labels), nchar(labels))
  testthat::expect_true(all(suffixes %in% c("1", "2")))
  testthat::expect_true(any(suffixes == "1") && any(suffixes == "2"))
  testthat::expect_s3_class(mix_merge, "orchard")
})
