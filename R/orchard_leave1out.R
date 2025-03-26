#' Orchard Plot for Leave-One-Out Analysis
#'
#' This function iteratively omits each element defined in a specified grouping variable from a meta-analytic model.
#' For each iteration, a new meta-analytic model is fitted and the resulting effect sizes are visualized using
#' an orchard plot. The plot displays the effect sizes, overall effect estimate (with confidence and prediction intervals),
#' and optionally, the effect sizes left out ("ghost points") in each iteration. In addition, the original model's
#' 95\% confidence limits are shown as reference lines.
#'
#' @param model A meta-analytic model object from the \pkg{metafor} package.
#' @param loo_output Optional data frame output from \code{leave_one_out()}. If provided, the function uses
#'   these pre-computed results instead of re-running the leave-one-out analysis.
#' @param group Character string naming the column in \code{model$data} that indicates the grouping variable.
#'   The analysis iteratively omits each unique element of this variable.
#' @param ylab Optional label for the y-axis, which lists the elements left out.
#' @param vcalc_args Optional list of arguments for the variance-covariance calculation using 
#'   metafor's vcalc function. Must include 'vi', 'cluster', 'obs', and 'rho' elements.
#' @param robust_args Optional list of arguments for robust variance estimation using
#'   metafor's robust function. Must include a 'cluster' element.
#' @param phylo_args Optional list of arguments for phylogenetic matrix calculation using
#'   ape's vcv function. Must include 'tree' and 'species_colname' elements. 'tree' is a phylogenetic
#'   tree object and 'species_colname' is the name of the column in the model data that is linked to the phylo matrix.
#' @param ci_lines Logical indicating whether to add horizontal reference lines (default is \code{TRUE})
#'   representing the original model's 95\% confidence intervals.
#' @param ci_lines_color Character string specifying the color for the confidence interval lines (default is \code{"red"}).
#' @param ghost_points Logical indicating whether to add "ghost points" for the effect sizes that were left out.
#'   These are plotted as empty points (default is \code{TRUE}).
#' @param ... Additional arguments passed to \code{orchard_plot}.
#'
#' @return A ggplot2 object displaying the leave-one-out analysis with reference lines
#'   for the original model's confidence limits.
#'
#' @details
#' This function is useful to see how sensitive are the estimates from 'model' to the absence of each element from \code{group}.
#' 
#' It fits multiple meta-analytic models using \code{model} as a template. In each iteration it leaves out one element of \code{group}.
#' The y-axis shows the elements from \code{group}. Each y-label indicates the element from \code{group} that was left-out in that model.
#' For each model, the plot shows the effect sizes, overall effect estimate with confidence and prediction intervals, and "ghost points". 
#' This "ghost points" are the effect sizes left out and they are shown as empty points. They can be ommited setting 'ghost_points' to FALSE
#' 
#' The plot also shows 95% confidence interval for the estimated overall effect by \code{model}. This is shown
#' as red dotted lines. The arguments 'ci_lines' and 'ci_lines_color' allow to remove this lines or select their color.
#'
#' The function can be slow because it has to fit various meta-analytic models. 
#' This can be a problem if you need to repeatedly call \code{orchard_leave1out()}, 
#' for example when trying to make aesthetic changes to the plot. 
#' To avoid this problem, you can pass the result from a \code{orchaRd::leave_one_out()} so
#' the function doesn't have to re-run it every time. Note that you still has to pass \code{model} and \code{group}.
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
#'
#'   # Generate an orchard plot using leave-one-out analysis
#'   orchard_leave1out(model = res, group = "paper_ID", xlab = "lnRR")
#'
#'   # Alternatively, precompute the leave-one-out results for faster re-plotting
#'   loo_res <- orchaRd::leave_one_out(res, group = "paper_ID")
#'   orchard_leave1out(model = res, loo_output = loo_res, group = "paper_ID", xlab = "lnRR")
#' }
#'
#' @export
orchard_leave1out <- function(model,
                              loo_output = NULL,
                              group,
                              ylab = NULL,
                              vcalc_args = NULL,
                              robust_args = NULL,
                              phylo_args = NULL,
                              ci_lines = TRUE,
                              ci_lines_color = "red",
                              ghost_points = TRUE,
                              angle = 0,
                              g = FALSE,
                              transfm = c("none", "tanh", "invlogit", "percent", "percentr"),
                              ...) {

  transfm <- match.arg(transfm)

  # Extract original model estimates (estimates and confidence intervals)
  orig_table <- mod_results(model, group = group)$mod_table
  orig_results <- orig_table[, c("estimate", "lowerCL", "upperCL"), drop = FALSE]

  # Apply transformation if needed
  if (transfm != "none") {
    orig_results <- transform_data(orig_results, transf = transfm)
  }
  
  # If not leave-one-out outputs is provided, run leave-one-out analysis. 
  # Each iteration omits one element from 'group'
  if (is.null(loo_output)) { 
    loo_output <- leave_one_out(model = model,
                                group = group,
                                vcalc_args = vcalc_args,
                                robust_args = robust_args,
                                phylo_args = phylo_args)
  }

  # Set colors for the plot. Check if 20 colours are enough 
  color_palette <- .set_color_palette(loo_output)

  # Create the base plot using orchard_plot() and then
  # add lines for the original model's 95% confidence interval 
  p <- orchard_plot(loo_output,
                    angle = angle,
                    transfm = transfm,
                    cb = FALSE,
                    ...)

  # Add layers to the basic orchard_plot
  if (ci_lines) {
    p <- p + .add_ci_lines(orig_results, ci_lines_color)
  }

  if (ghost_points) {
    p <- p + .add_ghost_points(model, group, transfm)
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
#' @param loo_results A list or data frame containing the leave-one-out analysis results,
#'   including a \code{mod_table} component with one row per leave-one-out model.
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
#' @param model A meta-analytic model object from the \pkg{metafor} package.
#' @param group Character string specifying the column name in \code{model$data} that
#'   identifies the grouping variable. Each unique value in this variable is omitted in turn.
#' @param transfm Character string with the transformation to be applied
#'
#' @return A \pkg{ggplot2} layer (via \code{geom_quasirandom}) displaying the ghost points.
#'
#' @details
#' The function extracts the removed points from the original model's data, converts the grouping
#' variable to a factor, and computes a scaling factor (the inverse square root of the variance)
#' for the point sizes.
#'
#' @author Facundo Decunta - fdecunta@agro.uba.ar
#'
#' @keywords internal
.add_ghost_points <- function(model, group, transfm) {
  # Add the removed points from each study in red
  removed_points <- mod_results(model, group = group)$data 

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
