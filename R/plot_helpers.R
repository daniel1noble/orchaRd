#' Colour blind palette
#' @keywords intern
.colour_blind_palette <- c(
  "#88CCEE", "#CC6677", "#DDCC77", "#117733",
  "#332288", "#AA4499", "#44AA99", "#999933",
  "#882255", "#661100", "#6699CC", "#888888",
  "#E69F00", "#56B4E9", "#009E73", "#F0E442",
  "#0072B2", "#D55E00", "#CC79A7", "#999999"
)


#' Get the output from mod_results 
#' 
#' Take the inputs passed to plotting functions and make sure to convert them
#' into the results from \code{mod_results}. 
#'
#' @keywords internal
.get_results <- function(object, mod, group, N, by, at, weights) {
  # Get results to be used. Only accepts metafor or orchard object objects
  if (any(class(object) %in% c("robust.rma", "rma.mv", "rma", "rma.uni"))) {
    results <- orchaRd::mod_results(
      object, 
      mod = mod,  
      group = group, 
      N = N,
      by = by, 
      at = at, 
      weights = weights
    )
  } else if (any(class(object) %in% c("orchard"))) {
    results <- object
  } else {
    stop("object argument must be a 'metafor' model or the output from 'mod_results'", .call = FALSE)
  }

  return(results)
}


#' Set scale for bubble plot.
#' @keywords internal
.get_size_scale <- function(N, dat) {
  if (is.null(N)) {
    scale <- (1 / sqrt(dat[, "vi"]))
    scale_legend <- "Precision (1/SE)"
  }
  else {
    if (!(N %in% names(dat))) {
      stop(sprintf(
        "The column %s doesn't exist in the data. If you use the output from 'mod_results()', make sure that you called 'mod_results(..., N = \"%s\")'",
        N, N), call. = FALSE)
    }
    scale <- dat$N
    scale_legend <- paste0("Sample Size ($\\textit{N}$)")
  } 
  return(list(scale = scale, scale_legend = scale_legend))
}


#' Compute k and g labels for each condition.
#'
#' Takes data from mod_results and compute the K (number of effect sizes)
#' and G (number of groups). 
#' NOTE: It expects the condition from data to be a factor
#' 
#' @keywords internal
.get_kg_labels <- function(dat) {
  cond <- dat$condition

  K <- by(dat, cond, function(x) length(x[, "yi"]))
  G <- by(dat, cond, function(x) length(unique(x[, "stdy"])))

  # Return a table with k, g and condition.
  # Condition must be a factor with the same levels
  # and order than the data used
  kg_labels <- data.frame(
    K = as.vector(K),
    G = as.vector(G),
    condition = factor(levels(cond),
                       levels = levels(cond),
                       ordered = is.ordered(cond))
  )

  return(kg_labels)
}


