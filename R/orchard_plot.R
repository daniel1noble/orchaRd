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
  condition.lab = "Condition",
  legend.pos = c("bottom.right", "bottom.left", "top.right", "top.left",
                  "top.out", "bottom.out", "none"), 
  k.pos = c("right", "left", "none"),
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

  if (any(class(object) %in% c("robust.rma", "rma.mv", "rma", "rma.uni"))) {
    # TODO: Test this. I think that the old code has a useless if-else 
    #       to decide how to ise mod
    results <- orchaRd::mod_results(object, mod, group, N, by = by, at = at,
                                    weights = weights, upper = upper)
  } 
  
  if (any(class(object) %in% c("orchard"))) {
    results <- object
  }

  # -----------------------------------------
  # Prepare data

  # Get model table and data.
  mod_table <- results$mod_table
  data_trim <- results$data

  # Transform data if needed
  if (transfm != "none") {
    numeric_cols <- c("estimate", "lowerCL", "upperCL", "lowerPR", "upperPR")
    mod_table[, numeric_cols] <- transform_data(mod_table[, numeric_cols],
                                                n = n_transfm,
                                                transfm = transfm)
    data_trim$yi <- transform_data(data_trim$yi, n = n_transfm, transfm = transfm)
  }

  # making sure factor names match
  data_trim$moderator <- factor(data_trim$moderator, levels = mod_table$name, labels = mod_table$name)

  # Reorder data if needed
  # If tree.order isn't equal to NULL, and length of tree order does not 
  # match number of categories in categorical moderator, then stop with error.
  if (!is.null(tree.order) & length(tree.order)!= nlevels(data_trim[,'moderator'])) {
    stop("Length of 'tree.order' does not equal number of categories in moderator")
  }

  # If tree.order isn't equal to NULL but passes above check, then reorder mod table
  # according to custom order if there is one.
  if (!is.null(tree.order)) {
    data_trim$moderator<-factor(data_trim$moderator, levels = tree.order, labels = tree.order)
    mod_table <- mod_table %>% dplyr::arrange(factor(name, levels = tree.order))
  }


  # Set scale used for points sizes and it's legend.
  # If N is not null, it is the name of the column that has the sample size
  if (is.null(N) == FALSE) {
    data_trim$scale <- data_trim$N
    legend <- paste0("Sample Size ($\\textit{N}$)") # we want to use italic
  } else {
    # This is the default
    data_trim$scale <- (1 / sqrt(data_trim[, "vi"]))
    legend <- "Precision (1/SE)"
  }
  size_legend <- ggplot2::guide_legend(title = latex2exp::TeX(legend))

  # ----------------------------------------
  # Set color and fill
  # setting fruit colour
  if (colour == TRUE) {
    color <- as.factor(data_trim$stdy)
    color2 <- NULL
  } else {
    color <- data_trim$mod
    color2 <- mod_table$name
  }

  # whether we fill fruit or not
  if (fill == TRUE) {
    fill <- color
  } else {
    fill <- NULL
  }

  # ----------------------------------------
  # Build orchard plot by layers:
  #   1. Effect sizes as beeswarm and bubbles
  #   2. Horizontal reference line (at 0 by default)
  #   3. Confidence interval
  #   4. Prediction interval
  #   5. Point estimate
  #   6. Basic theme, labels and legend
  #   7. Colors and flip

  # 1,2 and 7 are common for every plot
  # 3 to 6 are different if conditions are used (i.e., `at` and `by` options)

  # ----------------------------------------
  # Parse options for plots with or without conditions:

  # If the plot has conditions it must add:
  # - Space between lines from different conditions
  # - Different shapes for each element in 'condition'
  # - Shapes for all conditions:                        
  # - Labels for conditions:                           

  if (names(mod_table)[2] == "condition") {
    condition_no <- length(unique(mod_table[, "condition"])) # How many conditions?

    lines_position <- ggplot2::position_dodge2(width = 0.3) # Add spacing for multiple lines
    points_shapes <- as.factor(mod_table$condition) # Different shapes for each condition
    shapes_values <- 20 + (1:condition_no) # One shape for each condition
    shape_legend <- ggplot2::guide_legend(title = latex2exp::TeX(condition.lab))
  } else {
    lines_position <- ggplot2::position_dodge2(width = 0.75) # Set to default width in ggplot2
    points_shapes <- factor(1) # Dummy variable, same shape for all
    shapes_values <- 21 # Only one shape
    shape_legend <- "none"
  }


  plt <- .base_orchard_plot(data_trim, color, fill, alpha) %>%
    .add_reference_line(alpha) %>% 
    .add_conf_intervals(mod_table, lines_position, branch.size) %>%
    .add_pred_intervals(mod_table, lines_position, trunk.size, twig.size) %>%
    .add_point_estimates(mod_table, color2, points_shapes, lines_position, trunk.size, shapes_values) %>%
    .add_orchard_theme() %>%
    .add_legends(legend.pos, size_legend, shape_legend) %>%
    .add_axis_labels(xlab, angle)

  if (k == TRUE && k.pos != "none") {
    plt <- .add_k_and_g(plt, k.pos, g, data_trim, mod_table)
  }
  if (cb == TRUE) {
    plt <- .add_colour_blind_palette(plt)
  }
  if (flip) {
    plt <- plt + ggplot2::coord_flip()
  }

  return(plt)
}


#' Create base for orchard plot
#'
#' @keywords internal
.base_orchard_plot <- function(data_trim, color, fill, alpha) {
  p <- ggplot2::ggplot() +
    ggbeeswarm::geom_quasirandom(
      data = data_trim,
      ggplot2::aes(
        y = yi,
        x = moderator,
        size = scale,
        colour = color,
        fill = fill
      ),
      alpha = alpha,
      shape = 21
    )
  return(p)
}


#' Add reference line for orchard plot
#'
#' @keywords internal
.add_reference_line <- function(plt, alpha) {
  plt +
    ggplot2::geom_hline(
      yintercept = 0,
      linetype = 2,
      colour = "black",
      alpha = alpha
    )
}


#' Add confidence intervals
#'
#' @keywords internal
.add_conf_intervals <- function(plt, mod_table, lines_position, branch.size) {
  plt +
    ggplot2::geom_linerange(
      data = mod_table,
      ggplot2::aes(
        x = name,
        ymin = lowerCL,
        ymax = upperCL
      ),
      position = lines_position,
      size = branch.size
    )
}


#' Add prediction intervals
#'
#' @keywords internal
.add_pred_intervals <- function(plt, mod_table, lines_position, trunk.size, twig.size) {
  plt +
    ggplot2::geom_linerange(
      data = mod_table,
      ggplot2::aes(y = estimate, x = name, min = lowerPR, max = upperPR),
      position = lines_position,
      size = trunk.size,
      linewidth = twig.size
    )
}


#' Add point estimates for orchard plot
#'
#' @keywords internal
.add_point_estimates <- function(plt, mod_table, color2, points_shapes, lines_position, trunk.size, shapes_values) {
  # NOTE: Point estimate previously used in 'geom_pointrange',
  # plotting the point estimate and prediction interval at the same time.
  # However, that was less flexible.
  # Now PI are plotted with 'geom_linerange' and point estimate with 'geom_point'.
  #
  # To make points look the same as before, 'geom_point' use
  # some factors multiplying size and stroke.
  plt +
    ggplot2::geom_point(
      data = mod_table,
      ggplot2::aes(y = estimate, x = name, fill = color2, shape = points_shapes),
      position = lines_position,
      size = trunk.size * 4,
      stroke = 1 + trunk.size * 0.07   # Increase the stroke of the point 7% the size of the point. I think it's nice that way
    ) +
    ggplot2::scale_shape_manual(values = shapes_values)
}


#' Theme for orchard plot
#'
#' @keywords internal
.add_orchard_theme <- function(plt) {
  plt + ggplot2::theme_bw() +
    ggplot2::theme(
      legend.title = ggplot2::element_text(size = 9),
      legend.direction = "horizontal",
      legend.background = ggplot2::element_blank()
    )
}


#' Add legends to the plot
#'
#' @keywords internal
.add_legends <- function(plt, legend.pos, size_legend, shape_legend) {
  # Define legend position and add the legend
  plt <- switch(legend.pos,
    "bottom.right" = plt + ggplot2::theme(legend.position = "inside", legend.justification = c(1, 0)),
    "bottom.left"  = plt + ggplot2::theme(legend.position = "inside", legend.justification = c(0, 0)),
    "top.right"    = plt + ggplot2::theme(legend.position = "inside", legend.justification = c(1, 1)),
    "top.left"     = plt + ggplot2::theme(legend.position = "inside", legend.justification = c(0, 1)),
    "top.out"      = plt + ggplot2::theme(legend.position = "top"),
    "bottom.out"   = plt + ggplot2::theme(legend.position = "bottom"),
    "none"         = plt + ggplot2::theme(legend.position = "none"),
    plt
  )

  plt <- plt +
    ggplot2::guides(
      size   = size_legend,
      shape  = shape_legend,
      colour = "none",
      fill   = "none"
    )

  return(plt)
}


#' Add axis labels to orchard plot
#'
#' @keywords internal
.add_axis_labels <- function(plt, xlab, angle) {
  plt + ggplot2::labs(y = xlab, x = "") +
    ggplot2::theme(axis.text.y = ggplot2::element_text(
      size = 10,
      colour = "black",
      hjust = 0.5,
      angle = angle
    ))
}


#' Add k and g to base orchard plot
#'
#' @param plot A ggplot object representing the base orchard plot.
#' @param k.pos Position for the annotations. Can be "right", "left", or a numeric value specifying the y-coordinate.
#' @param g Logical. If \code{TRUE}, adds the g values in parentheses after k.
#' @param data_trim A data frame containing the yi values used to calculate annotation positions.
#' @param mod_table A data frame containing the k (and optionally g) values for each group.
#'
#' @return ggplot object with added text annotations.
#' @keywords internal
.add_k_and_g <- function(plt, k.pos, g, data_trim, mod_table) {
  # Add in total effect sizes for each level
  mod_table$K <- as.vector(by(data_trim, data_trim[, "moderator"], function(x) length(x[, "yi"])))
  # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
  mod_table$g <- as.vector(num_studies(data_trim, moderator, stdy)[, 2])
  # Number of groups in a moderator & data points
  group_no <- length(unique(mod_table[, "name"]))

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

  plt <- plt + ggplot2::annotate(
    "text",
    x = x_pos,
    y = y_pos,
    label = label_text,
    parse = TRUE,
    hjust = hjust_val,
    size = 3.5
  )

  return(plt)
}


#' Set colour blind friendly palette
#'
#' @keywords internal
.add_colour_blind_palette <- function(plt) {
  colour_blind_palette <- c(
    "#88CCEE", "#CC6677", "#DDCC77", "#117733",
    "#332288", "#AA4499", "#44AA99", "#999933",
    "#882255", "#661100", "#6699CC", "#888888",
    "#E69F00", "#56B4E9", "#009E73", "#F0E442",
    "#0072B2", "#D55E00", "#CC79A7", "#999999"
  )

  plt +
    ggplot2::scale_fill_manual(values = colour_blind_palette) +
    ggplot2::scale_colour_manual(values = colour_blind_palette)
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
