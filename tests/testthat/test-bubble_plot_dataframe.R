context("Checking bubble_plot with data.frame input (PR #51 feature)...")

# Setup: prepare raw data from lim dataset
data(lim)
lim[, "year"] <- suppressWarnings(as.numeric(lim$year))
lim$vi <- 1 / (lim$N - 3)
lim_clean <- na.omit(lim)

# --- Basic data.frame plotting ---

test_that("bubble_plot works with a raw data.frame", {
  plt <- bubble_plot(
    object = lim_clean,
    yi = "yi",
    vi = "vi",
    stdy = "Article",
    mod = "year",
    group = "Article",
    xlab = "Year",
    ylab = "Effect size"
  )
  expect_s3_class(plt, "ggplot")
})

test_that("bubble_plot with data.frame produces correct axis labels", {
  plt <- bubble_plot(
    object = lim_clean,
    yi = "yi",
    vi = "vi",
    stdy = "Article",
    mod = "year",
    group = "Article",
    xlab = "Publication Year",
    ylab = "Zr"
  )
  expect_equal(plt$labels$x, "Publication Year")
  expect_equal(plt$labels$y, "Zr")
})

test_that("bubble_plot with data.frame does NOT have model fit geom_smooth layers", {
  plt <- bubble_plot(
    object = lim_clean,
    yi = "yi",
    vi = "vi",
    stdy = "Article",
    mod = "year",
    group = "Article"
  )
  # geom_smooth layers come from model fit lines (CI, PI, estimate).
  # A data-only plot should have none.
  layer_types <- vapply(plt$layers, function(l) class(l$geom)[1], character(1))
  expect_false("GeomSmooth" %in% layer_types,
    info = "Data-only bubble plot should not contain model fit smooth lines")
})

test_that("bubble_plot with model DOES have model fit geom_smooth layers", {
  model <- metafor::rma.mv(
    yi = yi, V = vi, mods = ~ year,
    random = list(~1 | Article, ~1 | Datapoint), data = lim_clean
  )
  plt <- bubble_plot(model, mod = "year", group = "Article")
  layer_types <- vapply(plt$layers, function(l) class(l$geom)[1], character(1))
  expect_true("GeomSmooth" %in% layer_types,
    info = "Model-based bubble plot should contain smooth lines")
})

# --- Conditioning variable (by) ---

test_that("bubble_plot with data.frame and by argument creates facets", {
  plt <- bubble_plot(
    object = lim_clean,
    yi = "yi",
    vi = "vi",
    stdy = "Article",
    mod = "year",
    group = "Article",
    by = "Environment"
  )
  expect_s3_class(plt, "ggplot")
  # Should have facet_wrap due to condition
  expect_false(is.null(plt$facet))
  expect_true(inherits(plt$facet, "FacetWrap"))
})

# --- Input validation ---

test_that("bubble_plot errors when yi is missing for data.frame", {
  expect_error(
    bubble_plot(
      object = lim_clean,
      vi = "vi",
      stdy = "Article",
      mod = "year",
      group = "Article"
    ),
    "'yi' must be specified"
  )
})

test_that("bubble_plot errors when vi is missing for data.frame", {
  expect_error(
    bubble_plot(
      object = lim_clean,
      yi = "yi",
      stdy = "Article",
      mod = "year",
      group = "Article"
    ),
    "'vi' must be specified"
  )
})

test_that("bubble_plot errors when stdy is missing for data.frame", {
  expect_error(
    bubble_plot(
      object = lim_clean,
      yi = "yi",
      vi = "vi",
      mod = "year",
      group = "Article"
    ),
    "'stdy' must be specified"
  )
})

test_that("bubble_plot errors when mod is missing for data.frame", {
  expect_error(
    bubble_plot(
      object = lim_clean,
      yi = "yi",
      vi = "vi",
      stdy = "Article",
      group = "Article"
    ),
    "'mod' must be specified"
  )
})

test_that("bubble_plot errors when a column name doesn't exist", {
  expect_error(
    bubble_plot(
      object = lim_clean,
      yi = "nonexistent_col",
      vi = "vi",
      stdy = "Article",
      mod = "year",
      group = "Article"
    ),
    "not found in the data.frame"
  )
})

test_that("bubble_plot errors when by column doesn't exist", {
  expect_error(
    bubble_plot(
      object = lim_clean,
      yi = "yi",
      vi = "vi",
      stdy = "Article",
      mod = "year",
      group = "Article",
      by = "nonexistent_col"
    ),
    "not found in the data.frame"
  )
})

# --- Existing model-based behavior is not broken ---

test_that("bubble_plot still works with model objects after data.frame feature", {
  model <- metafor::rma.mv(
    yi = yi, V = vi, mods = ~ Environment * year,
    random = list(~1 | Article, ~1 | Datapoint), data = lim_clean
  )
  test_res <- orchaRd::mod_results(model, mod = "year", group = "Article",
                                   weights = "prop", by = "Environment")

  # From orchard object
  plt1 <- bubble_plot(test_res, group = "Article", mod = "year",
                      legend.pos = "top.left")
  expect_s3_class(plt1, "ggplot")

  # From model directly
  plt2 <- bubble_plot(model, mod = "year", group = "Article",
                      by = "Environment", legend.pos = "top.left")
  expect_s3_class(plt2, "ggplot")
})

# --- Handles NAs gracefully ---

test_that("bubble_plot with data.frame drops NA rows", {
  lim_with_na <- lim_clean
  lim_with_na$yi[1:3] <- NA
  plt <- bubble_plot(
    object = lim_with_na,
    yi = "yi",
    vi = "vi",
    stdy = "Article",
    mod = "year",
    group = "Article"
  )
  expect_s3_class(plt, "ggplot")
  # Check that the data used has fewer rows
  plot_data <- plt$layers[[1]]$data
  expect_true(nrow(plot_data) < nrow(lim_with_na))
})

# --- Edge cases: type validation ---

test_that("bubble_plot errors when yi column is not numeric", {
  bad <- lim_clean
  bad$yi_char <- as.character(bad$yi)
  expect_error(
    bubble_plot(
      object = bad,
      yi = "yi_char", vi = "vi", stdy = "Article",
      mod = "year", group = "Article"
    ),
    "must be numeric"
  )
})

test_that("bubble_plot errors when vi column is not numeric", {
  bad <- lim_clean
  bad$vi_char <- as.character(bad$vi)
  expect_error(
    bubble_plot(
      object = bad,
      yi = "yi", vi = "vi_char", stdy = "Article",
      mod = "year", group = "Article"
    ),
    "must be numeric"
  )
})

test_that("bubble_plot errors when moderator column is not numeric", {
  expect_error(
    bubble_plot(
      object = lim_clean,
      yi = "yi", vi = "vi", stdy = "Article",
      mod = "Environment", group = "Article"
    ),
    "must be numeric"
  )
})

# --- Edge cases: non-finite values ---

test_that("bubble_plot drops Inf values in yi gracefully", {
  bad <- lim_clean
  bad$yi[1:2] <- Inf
  bad$yi[3] <- -Inf
  plt <- bubble_plot(
    object = bad,
    yi = "yi", vi = "vi", stdy = "Article",
    mod = "year", group = "Article"
  )
  expect_s3_class(plt, "ggplot")
  plot_data <- plt$layers[[1]]$data
  expect_true(nrow(plot_data) == nrow(lim_clean) - 3)
})

test_that("bubble_plot drops Inf values in vi gracefully", {
  bad <- lim_clean
  bad$vi[1] <- Inf
  plt <- bubble_plot(
    object = bad,
    yi = "yi", vi = "vi", stdy = "Article",
    mod = "year", group = "Article"
  )
  expect_s3_class(plt, "ggplot")
  plot_data <- plt$layers[[1]]$data
  expect_true(nrow(plot_data) == nrow(lim_clean) - 1)
})

test_that("bubble_plot errors when all rows are NA/non-finite", {
  bad <- lim_clean
  bad$yi <- NA_real_
  expect_error(
    bubble_plot(
      object = bad,
      yi = "yi", vi = "vi", stdy = "Article",
      mod = "year", group = "Article"
    ),
    "No complete observations"
  )
})

# --- Edge cases: zero/negative vi ---

test_that("bubble_plot warns when vi contains zero or negative values", {
  bad <- lim_clean
  bad$vi[1] <- 0
  bad$vi[2] <- -0.01
  expect_warning(
    bubble_plot(
      object = bad,
      yi = "yi", vi = "vi", stdy = "Article",
      mod = "year", group = "Article"
    ),
    "zero or negative"
  )
})

# --- Edge cases: numeric by column ---

test_that("bubble_plot errors when by column is numeric", {
  expect_error(
    bubble_plot(
      object = lim_clean,
      yi = "yi", vi = "vi", stdy = "Article",
      mod = "year", group = "Article",
      by = "N"
    ),
    "must be categorical"
  )
})

# --- Edge cases: transfm with data.frame ---

test_that("bubble_plot warns when transfm is used with data.frame", {
  expect_warning(
    bubble_plot(
      object = lim_clean,
      yi = "yi", vi = "vi", stdy = "Article",
      mod = "year", group = "Article",
      transfm = "tanh"
    ),
    "transfm.*ignored"
  )
})

# --- Edge cases: single observation ---

test_that("bubble_plot works with a single-row data.frame", {
  one_row <- lim_clean[1, , drop = FALSE]
  plt <- bubble_plot(
    object = one_row,
    yi = "yi", vi = "vi", stdy = "Article",
    mod = "year", group = "Article"
  )
  expect_s3_class(plt, "ggplot")
})

# --- Edge cases: tibble input ---

test_that("bubble_plot works with tibble input", {
  skip_if_not_installed("tibble")
  tbl <- tibble::as_tibble(lim_clean)
  plt <- bubble_plot(
    object = tbl,
    yi = "yi", vi = "vi", stdy = "Article",
    mod = "year", group = "Article"
  )
  expect_s3_class(plt, "ggplot")
})

# --- Edge cases: condition.order with data.frame + by ---

test_that("bubble_plot with data.frame respects condition.order", {
  envs <- sort(unique(lim_clean$Environment))
  rev_order <- rev(envs)
  plt <- bubble_plot(
    object = lim_clean,
    yi = "yi", vi = "vi", stdy = "Article",
    mod = "year", group = "Article",
    by = "Environment",
    condition.order = rev_order
  )
  expect_s3_class(plt, "ggplot")
  # The condition factor levels in the plot data should match the reversed order
  plot_data <- plt$layers[[1]]$data
  expect_equal(levels(plot_data$condition), rev_order)
})

# --- Edge cases: NAs only in the by/condition column ---

test_that("bubble_plot drops rows where only the by column is NA", {
  bad <- lim_clean
  bad$Environment[1:3] <- NA
  plt <- bubble_plot(
    object = bad,
    yi = "yi", vi = "vi", stdy = "Article",
    mod = "year", group = "Article",
    by = "Environment"
  )
  expect_s3_class(plt, "ggplot")
  plot_data <- plt$layers[[1]]$data
  expect_true(nrow(plot_data) == nrow(lim_clean) - 3)
})

# --- Edge cases: NAs in moderator column ---

test_that("bubble_plot drops rows where moderator is NA", {
  bad <- lim_clean
  bad$year[1:5] <- NA
  plt <- bubble_plot(
    object = bad,
    yi = "yi", vi = "vi", stdy = "Article",
    mod = "year", group = "Article"
  )
  expect_s3_class(plt, "ggplot")
  plot_data <- plt$layers[[1]]$data
  expect_true(nrow(plot_data) == nrow(lim_clean) - 5)
})
