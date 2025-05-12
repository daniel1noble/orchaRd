#' Transform effect sizes in an orchard object
#'
#' Applies a transformation to the effect sizes and their confidence/prediction intervals from a model result object returned by \code{\link{mod_results}}.
#'
#' @param results An object of class \code{orchard}, the output from \code{\link{mod_results}}. 
#' @param transfm A character string specifying the transformation to apply. Options include:
#'   - \code{"none"}: No transformation (default)
#'   - \code{"tanh"}: Hyperbolic tangent transformation
#'   - \code{"invlogit"}: Inverse logit transformation
#'   - \code{"percentr"}: Percentage relative change transformation
#'   - \code{"percent"}: Percentage transformation
#'   - \code{"inv_ft"}: Inverse Freeman-Tukey (double arcsine) transformation for proportions (use with caution)
#' @param n_transfm A numeric vector of sample sizes, used only when \code{transfm = "inv_ft"}. Defaults to \code{NULL}.
#'
#' @return A modified list object of class \code{orchard}, where the effect size estimates and their intervals in both \code{mod_table} and \code{data} are transformed accordingly.
#'
#' @keywords internal
transform_mod_results <- function(results, transfm, n_transfm) {
  mod_table <- results$mod_table
  data <- results$data

  numeric_cols <- c("estimate", "lowerCL", "upperCL", "lowerPR", "upperPR")
  mod_table[, numeric_cols] <- transform_data(mod_table[, numeric_cols],
                                              n = n_transfm,
                                              transfm = transfm)
  data$yi <- transform_data(data$yi, n = n_transfm, transfm = transfm)

  results <- list(mod_table = mod_table, data = data)
  class(results) <- c("orchard", "data.frame")

  return(results)
}


#' Apply a Transformation to Vector
#'
#' This function applies a specified transformation to a numeric vector.
#'
#' @param x A numeric vector.
#' @param n Used with "inv_ft". A numeric value representing the sample size for proportion.
#' @param transfm A character string specifying the transformation to apply. 
#'   Options include:
#'   - "none": No transformation (default)
#'   - "tanh": Hyperbolic tangent transformation
#'   - "invlogit": Inverse logit transformation
#'   - "percentr": Percentage relative change transformation
#'   - "percent": Percentage transformation
#'   - "inv_ft": Inverse Freeman-Tukey (double arcsine) transformation for proportions (use with caution)
#'
#' @return A numeric vector with the applied transformation.
#' @examples
#' \dontrun{
#' transform_data(c(-1, 0, 1), transfm = "tanh")
#' transform_data(c(-1, 0, 1), transfm = "invlogit")
#' transform_data(c(-1, 0, 1), transfm = "none")
#' }
#' @keywords internal
transform_data <- function(x, n = NULL, transfm = c("none", "tanh", "invlogit", "percentr", "percent", "inv_ft")) {
  tryCatch({
    transfm <- match.arg(transfm)
  }, error = function(e) {
    stop(sprintf(
                 "Invalid transformation '%s'. transf must be one of: %s",
                 transfm, paste(c("none", "tanh", "invlogit", "percentr", "percent", "inv_ft"), collapse = ", ")),
         call. = FALSE)
  })

  if (transfm == "none") {
    return(x)
  }

  if(is.null(n) && transfm == "inv_ft") {
    stop("Sample size for each proportion, 'n', must be provided for 'inv_ft' transformation using the n_transfm argument.")
  }

  transf_func <- switch(transfm,
                       "tanh"     = transf_tanh(x),
                       "invlogit" = transf_invlogit(x),
                       "percentr" = transf_percentr(x),
                       "percent"  = transf_percent(x),
                       "inv_ft"   = transf_inv_ft(x, n))

  return(transf_func)
}

#' @title Hyperbolic Tangent Transformation
#'
#' @param x A numeric vector.
#' @return The hyperbolic tangent of x.
#' @keywords internal
transf_tanh <- function(x) {
  tanh(x)
}

#' @title Inverse Logit Transformation
#'
#' @param x A numeric vector.
#' @return The inverse logit of x, computed as exp(x) / (1 + exp(x)).
#' @keywords internal
transf_invlogit <- function(x) {
  exp(x) / (1 + exp(x))
}

#' @title Percentage Relative Change Transformation
#'
#' @param x A numeric vector.
#' @return The percentage relative change transformation, computed as (exp(x) - 1) * 100.
#' @keywords internal
transf_percentr <- function(x) {
  (exp(x) - 1) * 100
}

#' @title Percentage Transformation
#'
#' @param x A numeric vector.
#' @return The percentage transformation, computed as exp(x) * 100.
#' @keywords internal
transf_percent <- function(x) {
  exp(x) * 100
}

#' @title Inverse Freeman-Tukey (Double Arcsine) Transformation for proportions
#'
#' @param x A numeric vector containing effect sizes.
#' @param n A numeric value representing the sample size for proportion.
#' @return Back transformed proportion.
#' @keywords internal
transf_inv_ft <- function(x, n) {
   nhm <- 1/(mean(1/n, na.rm = TRUE))
     z <- suppressWarnings(1/2 * (1 - sign(cos(2 * x)) * sqrt(1 - 
        (sin(2 * x) + (sin(2 * x) - 1/sin(2 * x))/nhm)^2)))
     z <- ifelse(is.nan(z), NA_real_, z)
     z[x > transf_ift(1, nhm)] <- 1
     z[x < transf_ift(0, nhm)] <- 0
    return(z)
}

#' @title Freeman-Tukey (Double Arcsine) Transformation for proportions
#'
#' @param x A numeric vector containing proportions.
#' @param n A numeric value representing the sample size for proportion.
#' @return Back transformed proportion.
#' @keywords internal
transf_ift  <- function(x, n){
     x <- x * n
     z <- 1/2 * (asin(sqrt(x/(n + 1))) + asin(sqrt((x + 1)/(n + 1))))
    return(z)
}
