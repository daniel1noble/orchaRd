## Tests for whitespace trimming in moderator/factor levels
## Ensures that leading/trailing whitespace in data columns does not cause
## mismatches, errors, or silent failures throughout the orchaRd pipeline.

# ---------- Shared fixture: model with whitespace-polluted moderator ----------

.make_ws_model <- function() {
  data(lim, package = "orchaRd")
  lim$Phylum <- paste0(" ", lim$Phylum, " ")
  lim$vi <- 1 / (lim$N - 3)
  metafor::rma.mv(
    yi = yi, V = vi, mods = ~ Phylum - 1,
    random = list(~1 | Article, ~1 | Datapoint), data = lim
  )
}

model_ws <- .make_ws_model()

# ---------- .trimws_col ----------

test_that(".trimws_col trims whitespace from character vector", {
  x <- c("  foo  ", "bar ", " baz")
  expect_equal(orchaRd:::.trimws_col(x), c("foo", "bar", "baz"))
})

test_that(".trimws_col trims whitespace from factor levels", {
  x <- factor(c("  A  ", "B ", " C"))
  result <- orchaRd:::.trimws_col(x)
  expect_equal(sort(levels(result)), c("A", "B", "C"))
  expect_equal(as.character(result), c("A", "B", "C"))
})

test_that(".trimws_col passes through numeric vectors unchanged", {
  x <- c(1.5, 2.3, 3.7)
  expect_identical(orchaRd:::.trimws_col(x), x)
})

# ---------- mod_results ----------

test_that("mod_results works with whitespace in categorical moderator", {
  skip_on_cran()
  res <- orchaRd::mod_results(model_ws, mod = "Phylum",
    group = "Article", N = "N")
  expect_s3_class(res, "orchard")
  # mod_table names should be trimmed and capitalised
  expect_false(any(grepl("^\\s|\\s$", res$mod_table$name)),
    info = "mod_table$name should have no leading/trailing whitespace")
  # data moderator should be trimmed
  expect_false(any(grepl("^\\s|\\s$", res$data$moderator)),
    info = "data$moderator should have no leading/trailing whitespace")
})

# ---------- get_data_raw ----------

test_that("get_data_raw produces clean moderator values from whitespace data", {
  skip_on_cran()
  raw <- orchaRd::get_data_raw(model_ws, mod = "Phylum", group = "Article",
    N = "N")
  expect_false(any(grepl("^\\s|\\s$", raw$moderator)),
    info = "get_data_raw moderator should have no leading/trailing whitespace")
})

# ---------- orchard_plot ----------

test_that("orchard_plot renders with whitespace in moderator levels", {
  skip_on_cran()
  p <- orchaRd::orchard_plot(model_ws, mod = "Phylum", group = "Article",
    xlab = "Zr", N = "N")
  expect_s3_class(p, "gg")
})

test_that("orchard_plot with mod.order works despite whitespace in data", {
  skip_on_cran()
  # Use clean mod.order values (no whitespace) — should still work
  # because the package trims internally
  clean_order <- c("Chordata", "Arthropoda", "Mollusca", "Echinodermata",
                   "Nematoda", "Platyhelminthes", "Rotifera")
  p <- orchaRd::orchard_plot(model_ws, mod = "Phylum", group = "Article",
    xlab = "Zr", N = "N", mod.order = clean_order)
  expect_s3_class(p, "gg")
})

test_that("orchard_plot with upper=FALSE works despite whitespace in data", {
  skip_on_cran()
  data(lim, package = "orchaRd")
  lim$Phylum <- paste0(tolower(lim$Phylum), " ")
  lim$vi <- 1 / (lim$N - 3)
  model_lower_ws <- metafor::rma.mv(
    yi = yi, V = vi, mods = ~ Phylum - 1,
    random = list(~1 | Article, ~1 | Datapoint), data = lim
  )
  p <- orchaRd::orchard_plot(model_lower_ws, mod = "Phylum",
    group = "Article", xlab = "Zr", N = "N", upper = FALSE)
  expect_s3_class(p, "gg")
})

# ---------- caterpillars ----------

test_that("caterpillars plot renders with whitespace in moderator levels", {
  skip_on_cran()
  res_ws <- orchaRd::mod_results(model_ws, mod = "Phylum",
    group = "Article", N = "N")
  p <- orchaRd::caterpillars(res_ws, group = "Article", xlab = "Zr")
  expect_s3_class(p, "gg")
})

# ---------- bubble_plot (model path) ----------

test_that("bubble_plot works with whitespace in condition/by variable", {
  skip_on_cran()
  data(lim, package = "orchaRd")
  lim$Environment <- paste0("  ", lim$Environment, "  ")
  lim$year <- as.numeric(lim$year)
  lim$vi <- 1 / (lim$N - 3)
  model_bp <- metafor::rma.mv(
    yi = yi, V = vi, mods = ~ Environment * year,
    random = list(~1 | Article, ~1 | Datapoint), data = na.omit(lim)
  )
  res_bp <- orchaRd::mod_results(model_bp, mod = "year",
    group = "Article", weights = "prop", by = "Environment")
  p <- orchaRd::bubble_plot(res_bp, group = "Article", mod = "year",
    xlab = "Year")
  expect_s3_class(p, "gg")
})

# ---------- bubble_plot (data-frame path) ----------

test_that("bubble_plot data-frame path handles whitespace in by column", {
  skip_on_cran()
  data(lim, package = "orchaRd")
  lim_raw <- na.omit(lim)
  lim_raw$vi <- 1 / (lim_raw$N - 3)
  lim_raw$year <- as.numeric(lim_raw$year)
  # Inject whitespace into the by-variable
  lim_raw$Environment <- paste0("  ", lim_raw$Environment, "  ")
  p <- orchaRd::bubble_plot(
    object = lim_raw, yi = "yi", vi = "vi", stdy = "Article",
    mod = "year", group = "Article", by = "Environment",
    xlab = "Year", ylab = "Zr"
  )
  expect_s3_class(p, "gg")
})
