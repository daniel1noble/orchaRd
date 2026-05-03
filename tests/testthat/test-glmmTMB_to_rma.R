skip_if_not_installed("glmmTMB")
skip_if_not_installed("lme4")

library(glmmTMB)

set.seed(101)

# glmmTMB in some environments documents equalto()/propto() without exporting
# helper functions. A simple identity wrapper makes the special formula terms
# evaluable while glmmTMB still interprets them internally.
if (!exists("equalto", mode = "function")) {
  equalto <- function(term, group, mat) term
}

if (!exists("propto", mode = "function")) {
  propto <- function(term, group, mat) term
}

make_simple_tmb_data <- function() {
  study <- factor(rep(seq_len(15), each = 3))
  moderator <- factor(rep(c("A", "B", "C"), length.out = length(study)))
  x_mod <- rnorm(length(study))
  study_re <- rnorm(nlevels(study), 0, sqrt(0.15))[study]
  residual <- rnorm(length(study), 0, sqrt(0.2))
  y <- 0.2 + 0.3 * x_mod + c(A = 0, B = 0.25, C = -0.15)[moderator] + study_re + residual
  vi <- runif(length(study), 0.01, 0.08)

  data.frame(
    y = y,
    vi = vi,
    study = study,
    moderator = moderator,
    x_mod = x_mod
  )
}

test_that("glmmTMB_to_rma creates an rma-compatible object for core orchaRd flows", {
  dat <- make_simple_tmb_data()

  fit_tmb <- glmmTMB(
    y ~ moderator + x_mod + (1 | study),
    data = dat,
    REML = TRUE
  )

  converted <- glmmTMB_to_rma(
    fit_tmb,
    yi = "y",
    vi = "vi",
    data = dat,
    measure = "GEN"
  )

  expect_s3_class(converted, "rma.mv")
  expect_equal(class(converted), c("rma.mv", "rma"))
  expect_identical(converted$backend, "glmmTMB")
  expect_equal(nrow(converted$b), 4)
  expect_equal(ncol(converted$vb), 4)
  expect_equal(length(converted$sigma2), 2)
  expect_equal(tail(converted$s.names, 1), "Residual")
  expect_equal(tail(converted$s.nlevels.f, 1), converted$k)

  overall <- mod_results(converted, mod = "1", group = "study")
  categorical <- mod_results(converted, mod = "moderator", group = "study")
  continuous <- mod_results(converted, mod = "x_mod", group = "study")

  expect_equal(nrow(overall$mod_table), 1)
  expect_equal(nrow(categorical$mod_table), 3)
  expect_equal(nrow(continuous$mod_table), 100)

  expect_s3_class(
    orchard_plot(converted, mod = "moderator", group = "study", xlab = "Effect size"),
    "ggplot"
  )
  expect_s3_class(
    bubble_plot(converted, mod = "x_mod", group = "study"),
    "ggplot"
  )
  expect_s3_class(
    caterpillars(converted, mod = "moderator", group = "study", xlab = "Effect size"),
    "ggplot"
  )

  i2_ratio <- i2_ml(converted, method = "ratio")
  r2_vals <- r2_ml(converted, data = dat)

  expect_length(i2_ratio, 3)
  expect_true(all(is.finite(i2_ratio)))
  expect_length(r2_vals, 2)
  expect_true(all(is.finite(r2_vals)))
})

test_that("converted glmmTMB objects fail clearly for unsupported refit workflows", {
  dat <- make_simple_tmb_data()

  fit_tmb <- glmmTMB(
    y ~ moderator + x_mod + (1 | study),
    data = dat,
    REML = TRUE
  )

  converted <- glmmTMB_to_rma(
    fit_tmb,
    yi = "y",
    vi = "vi",
    data = dat
  )

  expect_error(
    i2_ml(converted, method = "ratio", boot = 5),
    "not supported for objects converted with glmmTMB_to_rma",
    fixed = FALSE
  )

  expect_error(
    r2_ml(converted, data = dat, boot = 5),
    "not supported for objects converted with glmmTMB_to_rma",
    fixed = FALSE
  )

  expect_error(
    leave_one_out(converted, group = "study"),
    "leave_one_out\\(\\) is not supported",
    perl = TRUE
  )
})

test_that("glmmTMB_to_rma matches metafor for intercept-only equalto models", {
  set.seed(42)
  k_studies <- 12
  k_per_study <- 3
  study <- rep(seq_len(k_studies), each = k_per_study)
  id <- seq_len(length(study))
  vi <- rbeta(length(id), 2, 20)
  V_mat <- metafor::vcalc(vi, cluster = study, obs = id, rho = 0.5)
  u <- rnorm(k_studies, 0, sqrt(0.2))[study]
  m <- rnorm(length(id), 0, sqrt(0.25))
  e <- as.numeric(MASS::mvrnorm(1, rep(0, length(id)), Sigma = V_mat))
  y <- u + m + e

  dat <- data.frame(
    y = y,
    vi = vi,
    study = factor(study),
    id = id,
    obs = factor(id),
    g = 1L
  )

  fit_rma <- metafor::rma.mv(
    yi = y,
    V = V_mat,
    random = list(~1 | study, ~1 | id),
    data = dat,
    method = "REML"
  )

  fit_tmb <- glmmTMB(
    y ~ 1 + (1 | study) + equalto(0 + obs | g, V_mat),
    data = dat,
    REML = TRUE
  )

  converted <- glmmTMB_to_rma(
    fit_tmb,
    yi = "y",
    vi = "vi",
    data = dat,
    V = V_mat,
    measure = "GEN"
  )

  expect_equal(as.numeric(converted$b), as.numeric(fit_rma$b), tolerance = 1e-3)
  expect_equal(sqrt(diag(converted$vb)), sqrt(diag(fit_rma$vb)), tolerance = 1e-3)
  expect_equal(unname(converted$sigma2[1:2]), unname(fit_rma$sigma2), tolerance = 1e-3)
  expect_equal(unname(i2_ml(converted, method = "ratio")), unname(i2_ml(fit_rma, method = "ratio")), tolerance = 1e-3)
  expect_equal(unname(i2_ml(converted, method = "matrix")), unname(i2_ml(fit_rma, method = "matrix")), tolerance = 1e-3)

  converted_res <- mod_results(converted, mod = "1", group = "study")$mod_table
  metafor_res <- mod_results(fit_rma, mod = "1", group = "study")$mod_table

  expect_equal(converted_res$estimate, metafor_res$estimate, tolerance = 1e-3)
  expect_equal(converted_res$lowerCL, metafor_res$lowerCL, tolerance = 1e-3)
  expect_equal(converted_res$upperCL, metafor_res$upperCL, tolerance = 1e-3)
  expect_equal(converted_res$lowerPR, metafor_res$lowerPR, tolerance = 1e-3)
  expect_equal(converted_res$upperPR, metafor_res$upperPR, tolerance = 1e-3)
})

test_that("glmmTMB_to_rma matches metafor for categorical and continuous equalto models", {
  set.seed(99)
  k_studies <- 18
  k_per_study <- 3
  study <- rep(seq_len(k_studies), each = k_per_study)
  id <- seq_len(length(study))
  moderator <- factor(rep(c("Control", "Treatment"), length.out = length(study)))
  x_mod <- rnorm(length(study))

  vi <- rbeta(length(id), 2, 20)
  V_mat <- metafor::vcalc(vi, cluster = study, obs = id, rho = 0.5)
  u <- rnorm(k_studies, 0, sqrt(0.12))[study]
  m <- rnorm(length(id), 0, sqrt(0.18))
  e <- as.numeric(MASS::mvrnorm(1, rep(0, length(id)), Sigma = V_mat))
  y <- 0.1 + 0.35 * x_mod + c(Control = 0, Treatment = 0.4)[moderator] + u + m + e

  dat <- data.frame(
    y = y,
    vi = vi,
    study = factor(study),
    id = id,
    moderator = moderator,
    x_mod = x_mod,
    obs = factor(id),
    g = 1L
  )

  fit_rma <- metafor::rma.mv(
    yi = y,
    V = V_mat,
    mods = ~ moderator + x_mod,
    random = list(~1 | study, ~1 | id),
    data = dat,
    method = "REML"
  )

  fit_tmb <- glmmTMB(
    y ~ moderator + x_mod + (1 | study) + equalto(0 + obs | g, V_mat),
    data = dat,
    REML = TRUE
  )

  converted <- glmmTMB_to_rma(
    fit_tmb,
    yi = "y",
    vi = "vi",
    data = dat,
    V = V_mat
  )

  expect_equal(as.numeric(converted$b), as.numeric(fit_rma$b), tolerance = 5e-3)
  expect_equal(sqrt(diag(converted$vb)), sqrt(diag(fit_rma$vb)), tolerance = 5e-3)
  expect_equal(unname(converted$sigma2[1:2]), unname(fit_rma$sigma2), tolerance = 5e-3)

  categorical_conv <- mod_results(converted, mod = "moderator", group = "study")$mod_table
  categorical_rma <- mod_results(fit_rma, mod = "moderator", group = "study")$mod_table
  continuous_conv <- mod_results(converted, mod = "x_mod", group = "study")$mod_table
  continuous_rma <- mod_results(fit_rma, mod = "x_mod", group = "study")$mod_table

  expect_equal(categorical_conv$estimate, categorical_rma$estimate, tolerance = 2e-3)
  expect_equal(categorical_conv$lowerPR, categorical_rma$lowerPR, tolerance = 2e-3)
  expect_equal(continuous_conv$estimate, continuous_rma$estimate, tolerance = 3e-3)
  expect_equal(continuous_conv$upperCL, continuous_rma$upperCL, tolerance = 1.5e-2)
})

test_that("glmmTMB_to_rma keeps propto variance components in phylogenetic models", {
  set.seed(222)
  k_sp <- 10
  k_ep <- 2
  tree <- ape::rtree(k_sp, tip.label = paste0("sp", seq_len(k_sp)))
  phylo_cor <- ape::vcv(tree, corr = TRUE)

  study <- rep(seq_len(k_sp), each = k_ep)
  species <- paste0("sp", rep(seq_len(k_sp), each = k_ep))
  id <- seq_len(length(study))

  vi <- rbeta(length(id), 2, 20)
  V_mat <- metafor::vcalc(vi, cluster = study, obs = id, rho = 0.5)
  rownames(V_mat) <- colnames(V_mat) <- as.character(id)

  obs <- factor(id, levels = rownames(V_mat))
  phylo <- factor(species, levels = rownames(phylo_cor))

  study_re <- rnorm(k_sp, 0, sqrt(0.05))[study]
  phylo_re <- as.numeric(MASS::mvrnorm(1, rep(0, k_sp), Sigma = phylo_cor * 0.2))[as.integer(phylo)]
  residual <- rnorm(length(id), 0, sqrt(0.1))
  error <- as.numeric(MASS::mvrnorm(1, rep(0, length(id)), Sigma = V_mat))
  y <- 0.3 + study_re + phylo_re + residual + error

  dat <- data.frame(
    y = y,
    vi = vi,
    study = factor(study),
    species = species,
    id = id,
    phylo = phylo,
    obs = obs,
    g = 1L
  )

  fit_rma <- metafor::rma.mv(
    yi = y,
    V = V_mat,
    random = list(~1 | study, ~1 | species, ~1 | id),
    R = list(species = phylo_cor),
    data = dat,
    method = "REML"
  )

  fit_tmb <- glmmTMB(
    y ~ 1 +
      (1 | study) +
      propto(0 + phylo | g, phylo_cor) +
      equalto(0 + obs | g, V_mat),
    data = dat,
    REML = TRUE
  )

  converted <- glmmTMB_to_rma(
    fit_tmb,
    yi = "y",
    vi = "vi",
    data = dat,
    V = V_mat
  )

  expect_equal(converted$s.names, c("study", "phylo", "Residual"))
  expect_equal(unname(converted$sigma2[1:3]), unname(fit_rma$sigma2), tolerance = 3e-3)

  converted_res <- mod_results(converted, mod = "1", group = "study")$mod_table
  metafor_res <- mod_results(fit_rma, mod = "1", group = "study")$mod_table

  expect_equal(converted_res$estimate, metafor_res$estimate, tolerance = 3e-3)
  expect_equal(converted_res$lowerPR, metafor_res$lowerPR, tolerance = 1e-2)
})

# --- Error-branch and edge-case coverage for glmmTMB_to_rma() ---

test_that("glmmTMB_to_rma rejects non-glmmTMB objects", {
  expect_error(
    glmmTMB_to_rma(stats::lm(mpg ~ wt, data = mtcars), yi = "mpg", vi = "wt"),
    "must be a fitted glmmTMB object", fixed = TRUE
  )
})

test_that("glmmTMB_to_rma rejects invalid 'test' values", {
  dat <- make_simple_tmb_data()
  fit <- glmmTMB(y ~ moderator + (1 | study), data = dat, REML = TRUE)
  expect_error(
    glmmTMB_to_rma(fit, yi = "y", vi = "vi", data = dat, test = "wald"),
    "test", fixed = TRUE
  )
})

test_that("glmmTMB_to_rma errors when yi or vi columns are missing from data", {
  dat <- make_simple_tmb_data()
  fit <- glmmTMB(y ~ moderator + (1 | study), data = dat, REML = TRUE)
  expect_error(
    glmmTMB_to_rma(fit, yi = "not_there", vi = "vi", data = dat),
    "yi", fixed = TRUE
  )
  expect_error(
    glmmTMB_to_rma(fit, yi = "y", vi = "missing_vi", data = dat),
    "vi", fixed = TRUE
  )
})

test_that("glmmTMB_to_rma uses model$frame when data is omitted (numeric yi/vi inputs)", {
  dat <- make_simple_tmb_data()
  fit <- glmmTMB(y ~ moderator + (1 | study), data = dat, REML = TRUE)
  # data = NULL forces the model$frame fallback (line 55). Passing yi/vi as
  # numeric vectors also exercises the as.numeric() branches at lines 58-59.
  out <- glmmTMB_to_rma(fit, yi = dat$y, vi = dat$vi)
  expect_s3_class(out, "rma.mv")
  expect_equal(out$k, nrow(dat))
})

test_that("glmmTMB_to_rma test='t' warns when no study column can be guessed and ddf is unset", {
  dat <- make_simple_tmb_data()
  # Rename 'study' so the auto-detect heuristic finds nothing.
  names(dat)[names(dat) == "study"] <- "panel"
  fit <- glmmTMB(y ~ moderator + (1 | panel), data = dat, REML = TRUE)
  expect_warning(
    glmmTMB_to_rma(fit, yi = "y", vi = "vi", data = dat, test = "t"),
    "Could not find study column", fixed = TRUE
  )
})

test_that("glmmTMB_to_rma test='t' uses an explicit study_col when supplied", {
  dat <- make_simple_tmb_data()
  fit <- glmmTMB(y ~ moderator + (1 | study), data = dat, REML = TRUE)
  out <- glmmTMB_to_rma(fit, yi = "y", vi = "vi", data = dat,
                        test = "t", study_col = "study")
  # ddf is filled per parameter; should match k_studies - p
  expect_false(is.null(out$ddf))
  expect_equal(unique(unlist(out$ddf)), nlevels(dat$study) - out$p)
})
