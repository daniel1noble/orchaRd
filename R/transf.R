#' Apply a Transformation to Vector
#'
#' This function applies a specified transformation to a numeric vector.
#'
#' @param x A numeric vector.
#' @param transfm A character string specifying the transformation to apply. 
#'   Options include:
#'   - "none": No transformation (default)
#'   - "tanh": Hyperbolic tangent transformation
#'   - "invlogit": Inverse logit transformation
#'   - "percentr": Percentage relative change transformation
#'   - "percent": Percentage transformation
#'
#' @return A numeric vector with the applied transformation.
#' @examples
#' \dontrun{
#' transform_data(c(-1, 0, 1), "tanh")
#' transform_data(c(-1, 0, 1), "invlogit")
#' }
#' @keywords internal
transform_data <- function(x, transfm = c("none", "tanh", "invlogit", "percentr", "percent")) {
  tryCatch({
    transfm <- match.arg(transfm)
  }, error = function(e) {
    stop(sprintf(
                 "Invalid transformation '%s'. transf must be one of: %s",
                 transfm, paste(c("none", "tanh", "invlogit", "percentr", "percent"), collapse = ", ")),
         call. = FALSE)
  })

  if (transfm == "none") {
    return(x)
  }

  transf_func <- switch(transfm,
                       "tanh"     = transf_tanh,
                       "invlogit" = transf_invlogit,
                       "percentr" = transf_percentr,
                       "percent"  = transf_percent)

  do.call(transf_func, list(x))
}

#' Hyperbolic Tangent Transformation
#'
#' @param x A numeric vector.
#' @return The hyperbolic tangent of x.
#' @keywords internal
transf_tanh <- function(x) {
  tanh(x)
}

#' Inverse Logit Transformation
#'
#' @param x A numeric vector.
#' @return The inverse logit of x, computed as exp(x) / (1 + exp(x)).
#' @keywords internal
transf_invlogit <- function(x) {
  exp(x) / (1 + exp(x))
}

#' Percentage Relative Change Transformation
#'
#' @param x A numeric vector.
#' @return The percentage relative change transformation, computed as (exp(x) - 1) * 100.
#' @keywords internal
transf_percentr <- function(x) {
  (exp(x) - 1) * 100
}

#' Percentage Transformation
#'
#' @param x A numeric vector.
#' @return The percentage transformation, computed as exp(x) * 100.
#' @keywords internal
transf_percent <- function(x) {
  exp(x) * 100
}

