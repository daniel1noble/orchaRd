#' Orchard Plot for Leave-One-Out Analysis
#'
#' Plots the output of a leave-one-out analysis for a meta-analytic model.
#' The plot displays the effect sizes, overall effect estimate (with confidence and prediction intervals),
#' and optionally, the effect sizes left out ("ghost points") in each iteration. In addition, the original model's
#' 95\% confidence limits are shown as reference lines.
#'
#' @param leave1out Output from \code{leave_one_out()}, containing the leave-one-out analysis results.
#' @param ylab Optional label for the y-axis, which lists the elements left out.
#' @param ci_lines Logical indicating whether to add horizontal reference lines (default is \code{TRUE})
#'   representing the original model's 95\% confidence intervals.
#' @param ci_lines_color Character string specifying the color for the confidence interval lines (default is \code{"red"}).
#' @param ghost_points Logical indicating whether to add "ghost points" for the effect sizes that were left out.
#'   These are plotted as empty points (default is \code{TRUE}).
#' @param angle Numeric value specifying the angle of the x-axis labels in degrees (default is 0).
#' @param g Logical indicating whether to return the plot as a grob object (default is \code{FALSE}).
#' @param transfm Character string specifying the transformation to apply to the data. 
#'   Options are: "none" (default, no transformation), "tanh" (hyperbolic tangent),
#'   "invlogit" (inverse logit), "percent" (percentage), or "percentr" (percentage range).
#' @param ... Additional arguments passed to \code{orchard_plot}.
#'
#' @return A ggplot2 object displaying the leave-one-out analysis with reference lines
#'   for the original model's confidence limits.
#'
#' @details
#' This function is useful to see how sensitive is the overall effect size estimated by a meta-analytic model to exclusion of each element from a group.
#' 
#' The y-axis shows the elements from the group. For each one, the plot shows the result from the model that left-out that element.
#' For each model, the plot shows the effect sizes, overall effect estimate with confidence and prediction intervals, and "ghost points". 
#' This "ghost points" are the effect sizes left out and they are shown as empty points. They can be ommited setting 'ghost_points' to FALSE.
#' 
#' The plot also shows 95% confidence interval for the estimated overall effect. This is shown
#' as red dotted lines. The arguments 'ci_lines' and 'ci_lines_color' allow to remove this lines or select their color.
#'
#' The function internally sets a color palette for the plot. If more than 20 elements are in the \code{group}, a viridis
#' palette is used. Otherwise, a color-blind friendly palette is applied by default.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @examples
#' \dontrun{
#'   # Fit a multivariate meta-analytic model using the metafor package
#'   res <- metafor::rma.mv(lnrr, lnrr_vi, random = ~ 1 | paper_ID, data = fish)
#'   loo_res <- orchaRd::leave_one_out(res, group = "paper_ID")
#'   orchard_leave1out(leave1out = loo_res, xlab = "lnRR")
#' }
#'
#' @export
orchard_leave1out <- function(leave1out,
                              ylab = NULL,
                              ci_lines = TRUE,
                              ci_lines_color = "red",
                              ghost_points = TRUE,
                              angle = 0,
                              g = FALSE,
                              transfm = c("none", "tanh", "invlogit", "percent", "percentr"),
                              ...) {
  transfm <- match.arg(transfm)

  # Check that it is the output of leave_one_out()
  if (!identical(names(leave1out), c("mod_table", "data", "orig_mod_results"))) {
    stop("The 'leave1out' argument must be the output of 'leave_one_out()'")
  }

  # Extract original model estimates and confidence intervals for plotting confidente interval lines
  orig_table <- leave1out$orig_mod_results$mod_table
  orig_results <- orig_table[, c("estimate", "lowerCL", "upperCL"), drop = FALSE]

  # Apply transformation if needed
  if (transfm != "none") {
    orig_results <- lapply(orig_results, function(x) transform_data(x, transf = transfm))
  }
  
  # Set colors for the plot. Check if 20 colours are enough 
  color_palette <- .set_color_palette(leave1out)

  # Create the base plot using orchard_plot() 
  p <- orchard_plot(leave1out,
                    angle = angle,
                    transfm = transfm,
                    cb = FALSE,
                    ...)

  # Add layers to the basic orchard_plot
  if (ci_lines) {
    p <- p + .add_ci_lines(orig_results, ci_lines_color)
  }

  if (ghost_points) {
    p <- p + .add_ghost_points(leave1out$data, transfm)
  }

  p <- p +
    ggplot2::xlab(ylab) +    # xlab uses ylab? Yes. It is confusing, but because of the coord_flip it is like this.
    color_palette$color_layer +
    color_palette$fill_layer

  return(p)
}


#' Set Color Palette for Orchard Plot
#'
#' Internal function to determine the appropriate color palette based on the number of
#' elements in the leave-one-out results.
#'
#' When the number of rows in the leave-one-out results exceeds 20, a viridis color
#' scale is used. Otherwise, a color-blind friendly manual palette is applied.
#'
#' @param loo_results Output from \code{leave_one_out()}, which contains the leave-one-out analysis results.
#'
#' @return A list with two components:
#'   \item{color_layer}{A \pkg{ggplot2} scale object for mapping colors.}
#'   \item{fill_layer}{A \pkg{ggplot2} scale object for mapping fill colors.}
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.set_color_palette <- function(loo_results) {
  color_palette <- list(color_layer = NULL, fill_layer = NULL)

  # Are 20 colours enough? Use viridis if not
  if (nrow(loo_results$mod_table) > 20) {
    color_palette$color_layer <- ggplot2::scale_color_viridis_d()
    color_palette$fill_layer <- ggplot2::scale_fill_viridis_d()
  } else {
    color_blind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733",
                             "#332288", "#AA4499", "#44AA99", "#999933", 
                             "#882255", "#661100", "#6699CC", "#888888",
                             "#E69F00", "#56B4E9", "#009E73", "#F0E442",
                             "#0072B2", "#D55E00", "#CC79A7", "#999999")
    color_palette$color_layer <- ggplot2::scale_color_manual(values = color_blind_palette)
    color_palette$fill_layer <- ggplot2::scale_fill_manual(values = color_blind_palette)
  }

  return(color_palette)
}



#' Add Confidence Interval Lines to Orchard Plot
#'
#' Internal function that creates a \pkg{ggplot2} layer to add horizontal lines
#' corresponding to the original model's 95% confidence limits.
#'
#' These lines help visualize how the leave-one-out effect size estimates compare to the
#' overall model estimate.
#'
#' @param orig_results A data frame containing the original model's estimates. It must
#'   include the columns \code{lowerCL} and \code{upperCL} representing the lower and upper
#'   confidence limits, respectively.
#' @param ci_lines_color Character string specifying the color for the confidence interval lines.
#'
#' @return A \pkg{ggplot2} layer (via \code{geom_hline}) that can be added to a plot.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.add_ci_lines <- function(orig_results, ci_lines_color) {
  ggplot2::geom_hline(yintercept = c(orig_results$lowerCL, orig_results$upperCL),
                      linetype = 2,
                      color = ci_lines_color) 
}


#' Add Ghost Points to Orchard Plot
#'
#' Internal function to add ghost points to an orchard plot.
#'
#' Ghost points represent the effect sizes that were omitted in each leave-one-out iteration.
#' They are plotted as empty points with reduced opacity to indicate their auxiliary role.
#'
#' @param data Data from the leave-one-out analysis, specifically the \code{data} component
#'   containing effect sizes (yi), variances (vi), and grouping variable (stdy).
#' @param transfm Character string indicating the transformation to be applied to the data.
#'   See \code{orchard_leave1out()} for available transformation options.
#'
#' @return A \pkg{ggplot2} layer (via \code{geom_quasirandom}) displaying the ghost points.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.add_ghost_points <- function(data, transfm) {

  # Names are weird because it is the only way to make it work with orchard_plot() without changes.
  # 'stdy' is the name of the grouping variable. So each category from 'stdy' is a different group
  # Plotting the effect sizes using 'stdy' as x-axis will plot the points left-out in each iteration.
  removed_points <- unique(data[, c("yi", "vi", "stdy")])
  removed_points$stdy <- as.factor(removed_points$stdy)
  removed_points$scale <- 1/sqrt(removed_points$vi)

  if (!is.null(transfm)) {
    # NOTE: For some reason 'scale' should not be transformed. It makes the points too small.
    # I should inspect how it is done in orchard_plot()
    numeric_col <- c("yi", "vi")
    removed_points[, numeric_col] <- transform_data(removed_points[, numeric_col], transfm)
  }
  
  ggbeeswarm::geom_quasirandom(data = removed_points,
                               ggplot2::aes(x = stdy, y = yi, size = scale),
                               alpha = 0.3,
                               shape = 21)
}
