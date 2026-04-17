#' Convert a glmmTMB model to an rma.mv-compatible object
#'
#' Wraps a fitted \pkg{glmmTMB} model in a light-weight object with class
#' \code{c("rma.mv", "rma")} so it can be passed to \pkg{orchaRd} functions
#' that operate on \pkg{metafor} model objects. The main use-case is
#' multilevel meta-analysis fitted with \pkg{glmmTMB}, including models that
#' use structured covariance terms such as \code{equalto()} and
#' \code{propto()}.
#'
#' The returned object is intended for analytical summaries and plotting
#' functions that only consume the fitted model components. Workflows that
#' require refitting a \pkg{metafor} model, such as \code{leave_one_out()} or
#' bootstrap paths inside some heterogeneity helpers, are not supported for
#' converted objects.
#'
#' @param model A fitted \code{glmmTMB} model object.
#' @param yi Name of the effect-size column in \code{data} (character), or a
#'   numeric vector of effect sizes.
#' @param vi Name of the sampling-variance column in \code{data} (character), or
#'   a numeric vector of sampling variances.
#' @param data Optional data frame. Defaults to \code{model$frame}.
#' @param V Optional known sampling-error variance-covariance matrix. Supply the
#'   same matrix used in \code{equalto()} if you want
#'   \code{orchaRd::i2_ml(method = "matrix")} to work on the converted object.
#' @param measure Character string labelling the effect-size metric. Stored as an
#'   attribute on \code{yi}. Defaults to \code{"GEN"}.
#' @param test Either \code{"z"} or \code{"t"}.
#' @param ddf Optional denominator degrees of freedom used when
#'   \code{test = "t"}.
#' @param study_col Optional study-level grouping variable used for automatic
#'   \code{ddf} calculation when \code{test = "t"} and \code{ddf = NULL}.
#'
#' @return A list of class \code{c("rma.mv", "rma")} containing the slots that
#'   \pkg{orchaRd} and \pkg{metafor} summary helpers expect.
#' @export
glmmTMB_to_rma <- function(model,
                           yi,
                           vi,
                           data = NULL,
                           V = NULL,
                           measure = NULL,
                           test = "z",
                           ddf = NULL,
                           study_col = NULL) {

  if (!inherits(model, "glmmTMB")) {
    stop("'model' must be a fitted glmmTMB object.", call. = FALSE)
  }

  if (!test %in% c("z", "t")) {
    stop("'test' must be \"z\" or \"t\".", call. = FALSE)
  }

  if (is.null(data)) {
    data <- model$frame
  }

  yi_vec <- if (is.character(yi)) data[[yi]] else as.numeric(yi)
  vi_vec <- if (is.character(vi)) data[[vi]] else as.numeric(vi)

  if (is.null(yi_vec)) {
    stop("'yi' column not found in data.", call. = FALSE)
  }
  if (is.null(vi_vec)) {
    stop("'vi' column not found in data.", call. = FALSE)
  }

  attr(yi_vec, "measure") <- if (!is.null(measure)) measure else "GEN"

  beta <- glmmTMB::fixef(model)$cond
  beta_names <- names(beta)
  beta_names[beta_names == "(Intercept)"] <- "intrcpt"
  b_mat <- matrix(beta, ncol = 1L, dimnames = list(beta_names, NULL))
  vb_mat <- as.matrix(stats::vcov(model)$cond)
  rownames(vb_mat) <- colnames(vb_mat) <- beta_names

  ff_full <- model$modelInfo$allForm$formula
  ff_fixed <- lme4::nobars(ff_full)
  formula_mods <- stats::as.formula(
    paste("~", paste(deparse(ff_fixed[[3L]]), collapse = ""))
  )

  not_na <- !is.na(yi_vec) & !is.na(vi_vec)
  X_mat <- stats::model.matrix(formula_mods, data = data[not_na, , drop = FALSE])
  x_colnames <- colnames(X_mat)
  x_colnames[x_colnames == "(Intercept)"] <- "intrcpt"
  colnames(X_mat) <- x_colnames

  re_struc <- model$modelInfo$reStruc$condReStruc
  vc <- glmmTMB::VarCorr(model)$cond
  vc_keys <- names(vc)
  re_fmls <- names(re_struc)

  sigma2_parts <- numeric(0L)
  s_names <- character(0L)
  s_nlevels <- numeric(0L)

  for (i in seq_along(re_struc)) {
    block_code <- re_struc[[i]]$blockCode

    if (block_code == 15L) {
      next
    }

    vc_key <- vc_keys[i]
    mat <- as.matrix(vc[[vc_key]])
    re_parts <- trimws(strsplit(re_fmls[i], "\\|")[[1L]])
    lhs <- re_parts[1L]
    rhs <- if (length(re_parts) > 1L) re_parts[2L] else NULL

    if (block_code == 11L) {
      vars <- mat[1L, 1L]
      lab <- trimws(gsub("^0\\s*\\+\\s*", "", lhs))
      n_levels <- if (!is.null(lab) && lab %in% names(data)) {
        length(unique(data[[lab]][not_na]))
      } else {
        nrow(mat)
      }
    } else {
      vars <- diag(mat)
      if (length(vars) == 0L) {
        vars <- mat[1L, 1L]
      }
      lab <- if (!is.null(rhs)) rhs else lhs
      n_levels <- if (!is.null(rhs) && rhs %in% names(data)) {
        length(unique(data[[rhs]][not_na]))
      } else {
        nrow(mat)
      }
    }

    sigma2_parts <- c(sigma2_parts, vars)
    s_names <- c(
      s_names,
      if (length(vars) == 1L) {
        lab
      } else {
        paste0(lab, ".", colnames(mat))
      }
    )
    s_nlevels <- c(s_nlevels, rep.int(n_levels, length(vars)))
  }

  sigma_resid <- tryCatch(stats::sigma(model)^2, error = function(e) 0)
  sigma2_parts <- c(sigma2_parts, sigma_resid)
  s_names <- c(s_names, "Residual")
  s_nlevels <- c(s_nlevels, sum(not_na))
  names(sigma2_parts) <- s_names

  k <- sum(not_na)
  p <- nrow(b_mat)

  if (test == "t" && is.null(ddf)) {
    if (!is.null(study_col) && study_col %in% names(data)) {
      k_studies <- length(unique(data[[study_col]][not_na]))
    } else {
      study_guess <- intersect(
        c("study", "Study", "STUDY", "cluster", "Cluster"),
        names(data)
      )
      if (length(study_guess) > 0L) {
        k_studies <- length(unique(data[[study_guess[1L]]][not_na]))
      } else {
        warning(
          "Could not find study column for ddf calculation. Set 'ddf' or ",
          "'study_col' explicitly.",
          call. = FALSE
        )
        k_studies <- k
      }
    }
    ddf <- replicate(p, k_studies - p, simplify = FALSE)
  }

  out <- list(
    b = b_mat,
    beta = b_mat,
    vb = vb_mat,
    yi = yi_vec,
    vi = vi_vec,
    data = data,
    not.na = not_na,
    subset = NULL,
    X = X_mat,
    V = V,
    sigma2 = sigma2_parts,
    tau2 = 0,
    gamma2 = 0,
    s.names = s_names,
    s.nlevels.f = s_nlevels,
    g.levels.k = NULL,
    struct = NULL,
    formula.mods = formula_mods,
    test = test,
    ddf = ddf,
    k = k,
    p = p,
    backend = "glmmTMB",
    backend_call = model$call
  )

  class(out) <- c("rma.mv", "rma")
  out
}
