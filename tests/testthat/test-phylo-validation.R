test_that(".validate_phylo_args accepts phylo objects via inherits()", {
  skip_if_not_installed("ape")

  dummy_tree <- ape::rtree(5)
  dummy_tree$tip.label <- paste0("sp", seq_len(5))

  mock_model <- list(Rfix = c(species = TRUE))
  class(mock_model) <- "mock_model"

  phylo_args <- list(
    tree = dummy_tree,
    species_colname = "species"
  )

  expect_no_error(.validate_phylo_args(mock_model, phylo_args))
})
