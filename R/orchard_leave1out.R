#' Orchard Plot for Leave-One-Out Analysis
#'
#' Performs a leave-one-out analysis on a meta-analytic model and produces an orchard plot.
#'
#' @param model A meta-analytic model from the metafor package.
#' @param mod Character string specifying the model version (default "1").
#' @param group Character string naming the column in model$data to omit iteratively.
#' @param ylab Optional label for the y-axis.
#' @param vcalc_args Optional list of parameters for computing the variance-covariance matrix.
#' @param alpha Numeric value for plot element transparency (default 0.1).
#' @param angle Numeric value for x-axis label angle (default 0).
#' @param g Logical; whether to apply grouping in the plot (default FALSE).
#' @param ... Additional arguments passed to \code{orchard_plot}.
#'
#' @return A ggplot2 object displaying the leave-one-out analysis with reference lines
#'   for the original model's confidence limits.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @examples
#' \dontrun{
#'   res <- metafor::rma.mv(lnrr, lnrr_vi, random = ~ 1 | paper_ID, data = fish)
#'   orchard_leave1out(res, group = "paper_ID", ylab = "Study left out")
#' }
#'
#' @export
orchard_leave1out <- function(model,
                              mod = "1", 
                              group,
                              ylab = NULL,
                              vcalc_args = NULL,
                              angle = 0,
                              g = FALSE,
                              transfm = c("none", "tanh", "invlogit", "percent", "percentr"),
                              ...) {

  transfm <- match.arg(transfm)

  # Extract original model estimates (point estimate and confidence limits)
  orig_table <- mod_results(model, group = group)$mod_table
  orig_results <- orig_table[, c("estimate", "lowerCL", "upperCL"), drop = FALSE]

  # Apply transformation if needed
  if (transfm != "none") {
    orig_results <- transform_data(orig_results, transf = transfm)
  }
  
  # Run leave-one-out analysis. Each iteration omits one element from 'group'
  leave1out_results <- leave_one_out(model = model, group = group, vcalc_args = vcalc_args)

  # Plot each result and include reference lines for the original model's confidence limits.
  p <- orchard_plot(leave1out_results,
                    angle = angle,
                    transfm = transfm,
                    g = g,
                    ...)

  p <- p +
    ggplot2::xlab(ylab) +    # xlab uses ylab? Yes. It is confusing, but because of the coord_flip it is like this.
    ggplot2::geom_hline(yintercept = orig_results$lowerCL, linetype = 2, color = "red") +
    ggplot2::geom_hline(yintercept = orig_results$upperCL, linetype = 2, color = "red")

  return(p)
}
