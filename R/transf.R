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

  if(is.null(n) && transfm == "inv_ft") {
    stop("Sample size for each proportion, 'n', must be provided for 'inv_ft' transformation using the n_transfm argument.")
  }

  if (transfm == "none") {
    return(x)
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
     z <- ifelse(is.nan(z), NA, z)
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