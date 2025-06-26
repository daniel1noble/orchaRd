#' @title orchard_plot
#' @description Using a \pkg{metafor} model object of class \code{rma} or \code{rma.mv}, or a results table of class \code{orchard}, it creates an orchard plot from mean effect size estimates for all levels of a given categorical moderator, and their corresponding confidence and prediction intervals.
#' @param object model object of class \code{rma.mv}, \code{rma}, or \code{orchard} table of model results.
#' @param mod the name of a moderator. Defaults to \code{"1"} for an intercept-only model. Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or any grouping variable one wishes to present sample sizes for. Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param by Character vector indicating the name that predictions should be conditioned on for the levels of the moderator.
#' @param at List of levels one wishes to predict at for the corresponding varaibles in 'by'. Used when one wants marginalised means. This argument can also be used to suppress levels of the moderator when argument \code{subset = TRUE}. Provide a list as follows: \code{list(mod = c("level1", "level2"))}.
#' @param weights Used when one wants marginalised means. How to marginalize categorical variables. The default is \code{weights = "prop"}, which weights moderator level means based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. In the case of sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param xlab The effect size measure label.
#' @param N The name of the column in the data specifying the sample size so that each effect size estimate is scaled to the sample size, N. Defaults to \code{NULL}, so that precision is used for scaling each raw effect size estimate instead of sample size.
#' @param alpha The level of transparency for effect sizes represented in the orchard plot.
#' @param angle The angle of y labels. The default is 90 degrees.
#' @param cb If \code{TRUE}, it uses 20 colour blind friendly colors.
#' @param k If \code{TRUE}, it displays k (number of effect sizes) on the plot.
#' @param g If \code{TRUE}, it displays g (number of grouping levels for each level of the moderator) on the plot.
#' @param transfm If set to \code{"tanh"}, a tanh transformation will be applied to effect sizes, converting Zr to a correlation or pulling in extreme values for other effect sizes (lnRR, lnCVR, SMD).  \code{"invlogit"} can be used to convert lnRR to the inverse logit scale. \code{"percentr"} can convert to the percentage change scale when using response ratios and \code{"percent"} can convert to the percentage change scale of an log transformed effect size. Defaults to \code{"none"}.
#' @param n_transfm The vector of sample sizes for each effect size estimate. This is used when \code{transfm = "inv_ft"}. Defaults to NULL.
#' @param condition.lab Label for the condition being marginalized over.
#' @param tree.order Order in which to plot the groups of the moderator when it is a categorical one. Should be a vector of equal length to number of groups in the categorical moderator, in the desired order (bottom to top, or left to right for flipped orchard plot)
#' @param trunk.size Size of the mean, or central point.
#' @param branch.size Size of the confidence intervals.
#' @param twig.size Size of the prediction intervals.
#' @param legend.pos Where to place the legend. To remove the legend, use \code{legend.pos = "none"}.
#' @param k.pos Where to put k (number of effect sizes) on the plot. Users can specify the exact position or they can use specify \code{"right"}, \code{"left"},  or \code{"none"}. Note that numeric values (0, 0.5, 1) can also be specified and this would give greater precision.
#' @param refline.pos Where to put the reference line. defaults to 0.
#' @param colour Colour of effect size shapes. By default, effect sizes are colored according to the \code{mod} argument. If \code{TRUE}, they are colored according to the grouping variable
#' @param fill If \code{TRUE}, effect sizes will be filled with colours. If \code{FALSE}, they will not be filled with colours.
#' @param weights Used when one wants marginalised means. How to marginalize categorical variables. The default is \code{weights = "prop"}, which weights moderator level means based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. In the case of sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param upper Logical, defaults to \code{TRUE}, indicating that the first letter of the character string for the moderator variable should be capitalized.
#' @param flip Logical, defaults to \code{TRUE}, indicating whether the plot should be flipped.
#' @return Orchard plot
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' data(eklof)
#' eklof<-metafor::escalc(measure="ROM", n1i=N_control, sd1i=SD_control,
#' m1i=mean_control, n2i=N_treatment, sd2i=SD_treatment, m2i=mean_treatment,
#' data=eklof)
#' # Add the unit level predictor
#' eklof$Datapoint<-as.factor(seq(1, dim(eklof)[1], 1))
#' # fit a MLMR - accounting for some non-independence
#' eklof_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~ Grazer.type-1,
#' random=list(~1|ExptID, ~1|Datapoint), data=eklof)
#' results <- mod_results(eklof_MR, mod = "Grazer.type", group = "ExptID")
#' orchard_plot(results, mod = "Grazer.type",
#' group = "ExptID", xlab = "log(Response ratio) (lnRR)")
#' # or
#' orchard_plot(eklof_MR, mod = "Grazer.type", group = "ExptID",
#' xlab = "log(Response ratio) (lnRR)")
#'
#' # Example 2
#' data(lim)
#' lim$vi<- 1/(lim$N - 3)
#' lim_MR<-metafor::rma.mv(yi=yi, V=vi, mods=~Phylum-1, random=list(~1|Article,
#' ~1|Datapoint), data=lim)
#' orchard_plot(lim_MR, mod = "Phylum", group = "Article",
#' xlab = "Correlation coefficient", transfm = "tanh", N = "N")
#' }
#' @export

orchard_plot <- function(
  object,
  mod = "1",
  group,
  xlab,
  N = NULL,
  alpha = 0.5,
  angle = 90,
  cb = TRUE,
  k = TRUE,
  g = TRUE,
  tree.order = NULL,
  trunk.size = 0.5,
  branch.size = 1.2,
  twig.size = 0.5,
  transfm = c("none", "tanh", "invlogit", "percent", "percentr", "inv_ft"),
  n_transfm = NULL,
  condition.lab = NULL,
  legend.pos = c("bottom.right", "bottom.left", "top.right", "top.left",
                  "top.out", "bottom.out", "none"), 
  k.pos = c("right", "left", "none"),
  refline.pos = 0,
  colour = FALSE,
  fill = TRUE,
  weights = "prop",
  by = NULL,
  at = NULL,
  upper = TRUE,
  flip = TRUE
) {
  ## evaluate choices, if not specified it takes the first choice
  transfm <- match.arg(NULL, choices = transfm)
  legend.pos <- match.arg(NULL, choices = legend.pos)
  k.pos <- match.arg(NULL, choices = k.pos)

  results <- .get_results(object, mod, group, N, by, at, weights)

  if (transfm != "none") {
    results <- transform_mod_results(results, transfm, n_transfm)
  }

  if (!is.null(tree.order)) {
    results <- .order_tree(results, tree.order)
  }

  mod_table <- results$mod_table
  data_trim <- results$data

  # Make sure levels from data and mod table are equally ordered
  # TODO: This should be done by mod_results: return factors with ordered levels
  data_trim$moderator <- factor(data_trim$moderator, levels = mod_table$name, labels = mod_table$name)

  size_scale <- .get_size_scale(N, data_trim)
  data_trim$scale <- size_scale$scale
  scale_legend <- size_scale$scale_legend
 
  if ("condition" %in% names(mod_table)) {
    condition.lab <- if (!is.null(condition.lab)) condition.lab else "Condition"
  }

  plt <- .base_orchard_plot(data_trim, colour, fill, alpha) +
    .orcd_theme(angle) +
    .orcd_reference_line(alpha, refline.pos) +
    .orcd_conf_intervals(mod_table, branch.size) +
    .orcd_pred_intervals(mod_table, trunk.size, twig.size) +
    .orcd_point_estimates(mod_table, colour, trunk.size) +
    .orcd_legends(legend.pos, scale_legend, condition.lab) +
    .orcd_axis_labels(xlab) +
    .orcd_k_and_g(k, k.pos, g, data_trim, mod_table) +
    .orcd_colour_blind_palette(cb) +
    if (flip) ggplot2::coord_flip()

  return(plt)
}


#' Create base for orchard plot
#'
#' @keywords internal
.base_orchard_plot <- function(data_trim, colour, fill, alpha) {

  # Defaults:  colour == FALSE ; fill == TRUE
  effectsize_color <- if (colour) as.factor(data_trim$stdy) else data_trim$moderator
  effectsize_fill <- if (fill) effectsize_color else NULL

  p <- ggplot2::ggplot() +
    ggbeeswarm::geom_quasirandom(
      data = data_trim,
      ggplot2::aes(
        y = yi,
        x = moderator,
        size = scale,
        colour = effectsize_color,
        fill = effectsize_fill
      ),
      alpha = alpha,
      shape = 21
    )
  return(p)
}

#' Theme for orchard plot
#'
#' @keywords internal
.orcd_theme <- function(angle) {
  list(
    ggplot2::theme_bw(),
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 9),
      legend.direction = "horizontal",
      legend.background = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 10,
                                          colour = "black",
                                          hjust = 0.5,
                                          angle = angle)
      )
  )
}


#' Add reference line for orchard plot
#'
#' @keywords internal
.orcd_reference_line <- function(alpha, refline.pos) {
  ggplot2::geom_hline(
    yintercept = refline.pos,
    linetype = 2,
    colour = "black",
    alpha = alpha
  )
}


#' Add confidence intervals
#'
#' @keywords internal
.orcd_conf_intervals <- function(mod_table, branch.size) {
  ggplot2::geom_linerange(
    data = mod_table,
    ggplot2::aes(x = name, ymin = lowerCL, ymax = upperCL),
    position = ggplot2::position_dodge2(width = .set_width(mod_table)),
    size = branch.size
  )
}


#' Add prediction intervals
#'
#' @keywords internal
.orcd_pred_intervals <- function(mod_table, trunk.size, twig.size) {
  if (twig.size == "none" || twig.size == 0)
    return(ggplot2::geom_blank())

  ggplot2::geom_linerange(
    data = mod_table,
    ggplot2::aes(y = estimate, x = name, min = lowerPR, max = upperPR),
    position = ggplot2::position_dodge2(width = .set_width(mod_table)),
    size = trunk.size,
    linewidth = twig.size
  )
}


#' Add point estimates for orchard plot
#'
#' Note: Point estimate previously used 'geom_pointrange',
#' plotting the estimate and prediction interval at the same time.
#' However, that was less flexible.
#' Now PI are plotted with 'geom_linerange' and point estimate with 'geom_point'.
#'
#' @keywords internal
.orcd_point_estimates <- function(mod_table, colour, trunk.size) {

  estimates_shape <- factor(1)    # Dummy var. Same shape for all
  shapes_values   <- 21           # Filled circles by default
  estimates_fill  <- if (colour) factor(1) else mod_table$name  # Default: colour == FALSE

  if ("condition" %in% names(mod_table)) {
    condition_no    <- nlevels(factor(mod_table$condition))
    shapes_values   <- 20 + (1:condition_no)           
    estimates_shape <- as.factor(mod_table$condition)  
  } 

  # To make points look the same as before, multiply size and stroke by those factors (found them by try and error)
  points_layer <- ggplot2::geom_point(
    data = mod_table,
    ggplot2::aes(y = estimate,
                 x = name,
                 fill = estimates_fill,
                 shape = estimates_shape),
      position = ggplot2::position_dodge2(width = .set_width(mod_table)), 
      size = trunk.size * 4,            # 4 times the trunk.size to make it look like before
      stroke = 1 + (trunk.size * 0.07) # Increase the stroke by 7% the size of the point. 
  )

  shapes_layer <- ggplot2::scale_shape_manual(values = shapes_values)

  return(list(points_layer, shapes_layer))
}


#' Add legends to orchard plot
#'
#' @keywords internal
.orcd_legends <- function(legend.pos, scale_legend, condition.lab) {
  # Define legend position and add the legend
  legend_position <- switch(legend.pos,
    "bottom.right" = ggplot2::theme(legend.position = "inside", legend.justification = c(1, 0)),
    "bottom.left"  = ggplot2::theme(legend.position = "inside", legend.justification = c(0, 0)),
    "top.right"    = ggplot2::theme(legend.position = "inside", legend.justification = c(1, 1)),
    "top.left"     = ggplot2::theme(legend.position = "inside", legend.justification = c(0, 1)),
    "top.out"      = ggplot2::theme(legend.position = "top"),
    "bottom.out"   = ggplot2::theme(legend.position = "bottom"),
    "none"         = ggplot2::theme(legend.position = "none")
  )

  size_legend <- ggplot2::guide_legend(title = latex2exp::TeX(scale_legend))
  shape_legend <- if (is.null(condition.lab)) {
    "none" 
  } else {
    ggplot2::guide_legend(title = latex2exp::TeX(condition.lab))
  }

  legend_layer <- ggplot2::guides(
    size = size_legend,
    shape = shape_legend,
    colour = "none",
    fill = "none"
    )

  return(list(legend_position, legend_layer))
}


#' Add axis labels to orchard plot
#'
#' @keywords internal
.orcd_axis_labels <- function(xlab) {
  ggplot2::labs(y = xlab, x = "")
}


#' Add k and g to base orchard plot
#'
#' @param k.pos Position for the annotations. Can be "right", "left", or a numeric value specifying the y-coordinate.
#' @param g Logical. If \code{TRUE}, adds the g values in parentheses after k.
#' @param data_trim A data frame containing the yi values used to calculate annotation positions.
#' @param mod_table A data frame containing the k (and optionally g) values for each group.
#'
#' @return ggplot object with added text annotations.
#' @keywords internal
.orcd_k_and_g <- function(k, k.pos, g, data_trim, mod_table) {
  # Early return if no need of k-g labels
  if (k == FALSE || k.pos == "none") return()

  mod_table$K <- as.vector(by(data_trim, data_trim[, "moderator"], function(x) length(x[, "yi"])))
  mod_table$g <- as.vector(num_studies(data_trim, moderator, stdy)[, 2])
  group_no <- nlevels(factor(mod_table$name))

  # Define x position
  x_pos <- seq(1, group_no, 1) + 0.3

  # Define y position based on k.pos parameter
  if (k.pos == "right") {
    y_pos <- max(data_trim$yi) * 1.1
    hjust_val <- "right"
  } else if (k.pos == "left") {
    y_pos <- min(data_trim$yi) * 1.1
    hjust_val <- "left"
  } else {
    y_pos <- k.pos # Manually defined by user
    hjust_val <- 0.5 # 0.5 because this is centered in ggplot
  }

  # Determine label based on whether g should be included
  if (g) {
    label_text <- paste("italic(k)==", mod_table$K[1:group_no], "~", "(", mod_table$g[1:group_no], ")")
  } else {
    label_text <- paste("italic(k)==", mod_table$K[1:group_no])
  }

  ggplot2::annotate(
    "text",
    x = x_pos,
    y = y_pos,
    label = label_text,
    parse = TRUE,
    hjust = hjust_val,
    size = 3.5
  )
}


#' Set a colour-blind-friendly palette
#'
#' @keywords internal
.orcd_colour_blind_palette <- function(cb) {
  if (cb) {
    return(
      list(
        ggplot2::scale_fill_manual(values = .colour_blind_palette),
        ggplot2::scale_colour_manual(values = .colour_blind_palette)
      )
    )
  }
}


#' @title Zr_to_r
#' @description Converts Zr back to r (Pearson's correlation coefficient)
#' @param df data frame of results of class 'orchard'
#' @return A data frame containing all the model results including mean effect size estimate, confidence and prediction intervals with estimates converted back to r
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @export
Zr_to_r <- function(df){
	return(sapply(df, tanh))
}


#' Reorder the tree
#'
#' Reorder the position of the moderators in the plot by reordering
#' the factos in the model table and the data.
#'
#' @keywords internal
.order_tree <- function(results, tree.order) {
  mod_table <- results$mod_table
  data_trim <- results$data

  # If length of tree order does not match number of categories
  # in categorical moderator, then stop with error.
  if (length(tree.order)!= nlevels(factor(data_trim[,'moderator']))) {
    stop("Length of 'tree.order' does not equal number of categories in moderator")
  }

  data_trim$moderator <- factor(data_trim$moderator, levels = tree.order, labels = tree.order)
  mod_table <- mod_table %>% dplyr::arrange(factor(name, levels = tree.order))

  return(list(mod_table = mod_table, data = data_trim))
}


#' Set width for ggplot2::position_dodge2()
#'
#' When the model use has some condition, the plot needs separation for
#' estimates from the same moderator.
#' This function adds space when a condition is present.
#' 
#' @keywords internal
.set_width <- function(mod_table) {
  if ("condition" %in% names(mod_table)) {
    dodge_width <- 0.35  # Add spacing for conditions
  } else {
    dodge_width <- 0.75  # Default in ggplot
  }
  return(dodge_width)
}
