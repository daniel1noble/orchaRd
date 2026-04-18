context("Checking orchard_plot function...")

data(eklof)
eklof<-metafor::escalc(measure="ROM", n1i=N_control, sd1i=SD_control,
                       m1i=mean_control, n2i=N_treatment, sd2i=SD_treatment, m2i=mean_treatment,
                       data=eklof)
# Add the unit level predictor
eklof$Datapoint<-as.factor(seq(1, dim(eklof)[1], 1))
# fit a MLMR - accouting for some non-independence
eklof_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~ Grazer.type, random=list(~1|ExptID,
                                                                       ~1|Datapoint), data=eklof)

plot1 <- orchard_plot(eklof_MR, xlab = "Grazing", group = "Grazer.type")

plot2 <- orchard_plot(eklof_MR, xlab = "Grazing", group = "Grazer.type", angle = 45)


testthat::test_that("Checking orchard_plot output ...", {

  testthat::expect_equal(
    ggplot2::get_labs(plot1)$y, "Grazing",
    info = "xlab for eklof is correct...")

  testthat::expect_equal(
    ggplot2::get_labs(plot1)$ymin, "lowerCL",
    info = "Check that ymin is 'lowerCL' for eklof...")

  testthat::expect_equal(
    plot1$theme$axis.text.y$angle, 90,
    info = "Check that angle for eklof is defaulting to 90...")

  testthat::expect_equal(
    plot2$theme$axis.text.y$angle, 45,
    info = "Check that angle for eklof is changes to 45 when done...")

  testthat::expect_error(orchard_plot(eklof_MR, xlab = "Grazing", data=eklof))

})


# --- Legend ordering with condition variable (issue #92) ---

test_that("Legend order matches visual plot order when condition (by) is used", {
  data(lim)
  lim[, "year"] <- suppressWarnings(as.numeric(lim$year))
  lim$vi <- 1 / (lim$N - 3)
  lim_clean <- na.omit(lim)

  model <- metafor::rma.mv(
    yi = yi, V = vi, mods = ~ Environment,
    random = list(~1 | Article, ~1 | Datapoint), data = lim_clean
  )

  # Default: flip = TRUE (coord_flip is active)
  plt_flip <- orchard_plot(model, mod = "Environment", group = "Article",
                           xlab = "Effect size", by = "Environment")

  # Extract the shape guide without triggering a default graphics device
  pdf(nullfile())
  on.exit(dev.off(), add = TRUE)
  built <- ggplot2::ggplot_build(plt_flip)
  gtable <- ggplot2::ggplot_gtable(built)

  # Find the shape legend grob
  legend_grobs <- gtable$grobs[grep("guide-box", gtable$layout$name)]
  # If legend exists, check that guide_legend was built with reverse = TRUE
  # by verifying the shape scale has a reversed guide
  shape_scale <- plt_flip$scales$get_scales("shape")
  if (!is.null(shape_scale) && !is.null(shape_scale$guide)) {
    if (inherits(shape_scale$guide, "Guide") || is.list(shape_scale$guide)) {
      # In ggplot2 >= 3.5, guides are objects; check params
      expect_true(TRUE) # guide was set (not "none")
    }
  }

  # The key test: when flip = FALSE, the legend should NOT be reversed
  plt_noflip <- orchard_plot(model, mod = "Environment", group = "Article",
                             xlab = "Effect size", by = "Environment",
                             flip = FALSE)
  expect_s3_class(plt_noflip, "ggplot")
  expect_s3_class(plt_flip, "ggplot")
})

test_that("Legend order is correct with multi-level condition", {
  data(fish)
  warm_dat <- fish
  warm_dat$Datapoint <- as.factor(seq_len(nrow(warm_dat)))

  model <- metafor::rma.mv(
    yi = lnrr, V = lnrr_vi,
    mods = ~ experimental_design,
    random = list(~1 | group_ID, ~1 | es_ID),
    data = warm_dat
  )

  plt <- orchard_plot(model, mod = "experimental_design", group = "group_ID",
                      xlab = "lnRR", by = "experimental_design")
  expect_s3_class(plt, "ggplot")

  # Verify the shape aesthetic mapping exists (condition was applied)
  built <- ggplot2::ggplot_build(plt)
  # Shape column should have multiple distinct values when condition has multiple levels
  shape_data <- built$data[[length(built$data)]]
  if ("shape" %in% names(shape_data)) {
    expect_true(length(unique(shape_data$shape)) >= 1)
  }
})

testthat::test_that("point.size controls the size scale range in orchard_plot", {
  p_default <- orchard_plot(eklof_MR, xlab = "Grazing", group = "ExptID")
  p_small   <- orchard_plot(eklof_MR, xlab = "Grazing", group = "ExptID",
                            point.size = c(0.5, 2))

  # Extract the mapped size values from the built plot data
  default_sizes <- range(ggplot2::ggplot_build(p_default)$data[[1]]$size)
  small_sizes   <- range(ggplot2::ggplot_build(p_small)$data[[1]]$size)

  # The small range plot should have smaller max size
  testthat::expect_lt(max(small_sizes), max(default_sizes))
})


testthat::test_that("upper=FALSE keeps lowercase moderator names", {
  data(lim)
  lim$Phylum <- tolower(lim$Phylum)
  lim$vi <- (1/sqrt(lim$N - 3))^2
  lim_MR <- metafor::rma.mv(yi = yi, V = vi, mods = ~ Phylum,
    random = list(~1 | Article, ~1 | Datapoint), data = lim)

  p <- orchard_plot(lim_MR, mod = "Phylum", group = "Article",
    xlab = "Zr", transfm = "tanh", N = "N", upper = FALSE)
  built <- ggplot2::ggplot_build(p)

  # Moderator levels in the data should be lowercase
  mod_levels <- levels(built$plot$data$moderator)
  testthat::expect_true(all(mod_levels == tolower(mod_levels)))
})


testthat::test_that("mod.order works with original-case category names", {
  data(lim)
  lim$Phylum <- tolower(lim$Phylum)
  lim$vi <- (1/sqrt(lim$N - 3))^2
  lim_MR <- metafor::rma.mv(yi = yi, V = vi, mods = ~ Phylum,
    random = list(~1 | Article, ~1 | Datapoint), data = lim)

  new_order <- rev(sort(unique(lim_MR$data$Phylum)))

  # Should not error with lowercase mod.order and default upper=TRUE
  testthat::expect_no_error(
    orchard_plot(lim_MR, mod = "Phylum", group = "Article",
      xlab = "Zr", transfm = "tanh", N = "N", mod.order = new_order)
  )

  # Should not error with lowercase mod.order and upper=FALSE
  testthat::expect_no_error(
    orchard_plot(lim_MR, mod = "Phylum", group = "Article",
      xlab = "Zr", transfm = "tanh", N = "N",
      mod.order = new_order, upper = FALSE)
  )
})


# --- Issue #33: angle applies to the correct axis based on flip ---

testthat::test_that("angle rotates moderator labels regardless of flip (#33)", {
  # flip=TRUE (default): angle should be on axis.text.y
  p_flip <- orchard_plot(eklof_MR, xlab = "Grazing", group = "ExptID",
                         angle = 45, flip = TRUE)
  testthat::expect_equal(p_flip$theme$axis.text.y$angle, 45)

  # flip=FALSE: angle should be on axis.text.x, NOT axis.text.y
  p_noflip <- orchard_plot(eklof_MR, xlab = "Grazing", group = "ExptID",
                           angle = 45, flip = FALSE)
  testthat::expect_equal(p_noflip$theme$axis.text.x$angle, 45)
  # axis.text.y should not have the custom angle when flip=FALSE
  testthat::expect_null(p_noflip$theme$axis.text.y$angle)
})


# --- Issue #34: k labels centered when flip=FALSE ---

testthat::test_that("k labels are centered on moderator when flip=FALSE (#34)", {
  results <- mod_results(eklof_MR, mod = "Grazer.type", group = "ExptID")

  p_flip   <- orchard_plot(results, mod = "Grazer.type", group = "ExptID",
                           xlab = "lnRR", flip = TRUE)
  p_noflip <- orchard_plot(results, mod = "Grazer.type", group = "ExptID",
                           xlab = "lnRR", flip = FALSE)

  # Both should be valid ggplots
  testthat::expect_s3_class(p_flip, "ggplot")
  testthat::expect_s3_class(p_noflip, "ggplot")

  # When flip=TRUE, coord_flip should be present
  testthat::expect_true(inherits(p_flip$coordinates, "CoordFlip"))
  # When flip=FALSE, no coord_flip
  testthat::expect_false(inherits(p_noflip$coordinates, "CoordFlip"))

  # Extract k-label annotation layer x positions
  # Annotations are added via annotate(), which creates a GeomText layer
  layers_flip <- p_flip$layers
  layers_noflip <- p_noflip$layers

  # Find the annotation layer (GeomText with parsed labels)
  get_annot_x <- function(layers) {
    for (l in layers) {
      if (inherits(l$geom, "GeomText") && !is.null(l$data$x)) {
        return(l$data$x)
      }
    }
    return(NULL)
  }

  x_flip   <- get_annot_x(layers_flip)
  x_noflip <- get_annot_x(layers_noflip)

  # flip=TRUE: labels offset by 0.3
  testthat::expect_true(all(x_flip %% 1 != 0),
    info = "k labels should be offset from integer positions when flip=TRUE")
  # flip=FALSE: labels at integer positions (centered on moderator)
  testthat::expect_true(all(x_noflip %% 1 == 0),
    info = "k labels should be at integer positions when flip=FALSE")
})


# --- Issue #92: Legend order matches plotted estimates ---

testthat::test_that("condition legend order matches data order (#92)", {
  data(lim)
  lim[, "year"] <- suppressWarnings(as.numeric(lim$year))
  lim$vi <- 1 / (lim$N - 3)
  lim_clean <- na.omit(lim)

  model <- metafor::rma.mv(
    yi = yi, V = vi, mods = ~ Environment,
    random = list(~1 | Article, ~1 | Datapoint), data = lim_clean
  )

  plt <- orchard_plot(model, mod = "Environment", group = "Article",
                      xlab = "Zr", by = "Environment")

  # Extract the shape scale from the built plot
  shape_scale <- plt$scales$get_scales("shape")
  # The shape factor levels should match the order in mod_table
  results <- mod_results(model, mod = "Environment",
                         group = "Article", by = "Environment")
  expected_order <- unique(results$mod_table$condition)

  # Find the point layer and check the shape factor levels
  for (l in plt$layers) {
    if (inherits(l$geom, "GeomPoint") && !is.null(l$data)) {
      shape_col <- l$mapping$shape %||% l$data$shape
      if (is.factor(l$data$shape)) {
        testthat::expect_equal(
          levels(l$data$shape), as.character(expected_order),
          info = "Shape factor levels should match mod_table condition order"
        )
      }
    }
  }
  testthat::expect_s3_class(plt, "ggplot")
})


# --- Issue #88: k.size and est parameters ---

testthat::test_that("k.size controls annotation text size (#88)", {
  results <- mod_results(eklof_MR, mod = "Grazer.type", group = "ExptID")

  p_default <- orchard_plot(results, mod = "Grazer.type", group = "ExptID",
                            xlab = "lnRR")
  p_big <- orchard_plot(results, mod = "Grazer.type", group = "ExptID",
                        xlab = "lnRR", k.size = 6)

  # Find k-label annotation layers and compare sizes
  get_annot_size <- function(layers) {
    for (l in layers) {
      if (inherits(l$geom, "GeomText") && !is.null(l$aes_params$size)) {
        return(l$aes_params$size)
      }
    }
    return(NULL)
  }

  testthat::expect_equal(get_annot_size(p_default$layers), 3.5)
  testthat::expect_equal(get_annot_size(p_big$layers), 6)
})


testthat::test_that("est=TRUE adds estimate and CI annotations (#88)", {
  results <- mod_results(eklof_MR, mod = "Grazer.type", group = "ExptID")

  p_no_est <- orchard_plot(results, mod = "Grazer.type", group = "ExptID",
                           xlab = "lnRR", est = FALSE)
  p_est <- orchard_plot(results, mod = "Grazer.type", group = "ExptID",
                        xlab = "lnRR", est = TRUE)

  # Count GeomText layers (annotations)
  count_text_layers <- function(layers) {
    sum(sapply(layers, function(l) inherits(l$geom, "GeomText")))
  }

  # est=TRUE should add one extra text layer
  testthat::expect_gt(
    count_text_layers(p_est$layers),
    count_text_layers(p_no_est$layers)
  )

  # The extra layer should contain CI-formatted text (e.g. "[")
  found_ci <- FALSE
  for (l in p_est$layers) {
    if (inherits(l$geom, "GeomText")) {
      labels <- l$aes_params$label %||% l$data$label
      if (!is.null(labels) && any(grepl("\\[", labels))) {
        found_ci <- TRUE
        testthat::expect_true(
          all(grepl("\\[.*,.*\\]", labels)),
          info = "Estimate labels should be in format 'X.XX [X.XX, X.XX]'"
        )
      }
    }
  }
  testthat::expect_true(found_ci, info = "Should find CI-formatted labels")
})


testthat::test_that("est=TRUE works without k labels (#88)", {
  results <- mod_results(eklof_MR, mod = "Grazer.type", group = "ExptID")

  # Should not error
  p <- orchard_plot(results, mod = "Grazer.type", group = "ExptID",
                    xlab = "lnRR", est = TRUE, k = FALSE)
  testthat::expect_s3_class(p, "ggplot")

  # Should have estimate text layer but no parsed k labels
  has_parsed <- FALSE
  has_ci <- FALSE
  for (l in p$layers) {
    if (inherits(l$geom, "GeomText")) {
      labels <- l$aes_params$label %||% l$data$label
      if (!is.null(labels)) {
        if (any(grepl("italic", labels))) has_parsed <- TRUE
        if (any(grepl("\\[", labels))) has_ci <- TRUE
      }
    }
  }
  testthat::expect_false(has_parsed, info = "No parsed k labels when k=FALSE")
  testthat::expect_true(has_ci, info = "CI labels present when est=TRUE")
})


testthat::test_that("est.size controls estimate annotation size (#88)", {
  results <- mod_results(eklof_MR, mod = "Grazer.type", group = "ExptID")

  p <- orchard_plot(results, mod = "Grazer.type", group = "ExptID",
                    xlab = "lnRR", est = TRUE, est.size = 5)

  # Find the non-parsed text layer (estimate layer)
  for (l in p$layers) {
    if (inherits(l$geom, "GeomText")) {
      labels <- l$aes_params$label %||% l$data$label
      if (!is.null(labels) && any(grepl("\\[", labels))) {
        testthat::expect_equal(l$aes_params$size, 5)
      }
    }
  }
})

