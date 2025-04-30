#' @title bubble_plot
#' @description Using a \pkg{metafor} model object of class \code{rma} or \code{rma.mv}, or a results table of class \code{orchard}, the \code{bubble_plot} function creates a bubble plot from slope estimates. In cases when a model includes interaction terms, this function creates panels of bubble plots.
#' @param object Model object of class \code{rma}, \code{rma.mv}, or \code{orchard} table of model results
#' @param mod The name of a continuous moderator, to be plotted on the x-axis of the bubble plot.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or any grouping variable one wishes to present sample sizes for. Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param by Character vector indicating the name that predictions should be conditioned on for the levels of the moderator.
#' @param at List of levels one wishes to predict at for the corresponding variables in \code{by}. Used when one wants marginalised means. This argument can also be used to suppress levels of the moderator when argument \code{subset = TRUE}. Provide a list as follows: \code{list(mod = c("level1", "level2"))}.
#' @param weights How to marginalize categorical variables; used when one wants marginalised means. The default is \code{weights = "prop"}, which weights means for moderator levels based on their proportional representation in the data. For example, if \code{"sex"} is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal when, for example, males and females are typically roughly equally prevalent in a population. In cases such as these, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param transfm If set to \code{"tanh"}, a tanh transformation will be applied to effect sizes, converting Zr to a correlation or pulling in extreme values for other effect sizes (lnRR, lnCVR, SMD).  \code{"invlogit"} can be used to convert lnRR to the inverse logit scale. \code{"percentr"} can convert to the percentage change scale when using response ratios and \code{"percent"} can convert to the percentage change scale of an log transformed effect size. Defaults to \code{"none"}.
#' @param n_transfm The vector of sample sizes for each effect size estimate. This is used when \code{transfm = "inv_ft"}. Defaults to NULL.
#' @param xlab Moderator label.
#' @param ylab Effect size measure label.
#' @param k.pos The position of effect size number, k.
#' @param N The vector of sample size which an effect size is based on. Defaults to precision (the inverse of sampling standard error).
#' @param alpha The level of transparency for pieces of fruit (effect size).
#' @param cb If \code{TRUE}, it uses a colourblind-friendly palette of 20 colours (do not make this \code{TRUE}, when colour = \code{TRUE}).
#' @param k If \code{TRUE}, it displays k (number of effect sizes) on the plot.
#' @param g If \code{TRUE}, it displays g (number of grouping levels for each level of the moderator) on the plot.
#' @param est.lwd Size of the point estimate.
#' @param ci.lwd Size of the confidence interval.
#' @param pi.lwd Size of the prediction interval.
#' @param est.col Colour of the point estimate.
#' @param ci.col Colour of the confidence interval.
#' @param pi.col Colour of the prediction interval.
#' @param condition.nrow Number of rows to plot condition variable.
#' @param condition.order Order of the levels of the condition variable in the plot. Defaults to NULL.
#' @param legend.pos Where to place the legend, or not to include a legend ("none").
#'
#' @return Bubble plot
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' data(lim)
#' lim[, "year"] <- as.numeric(lim$year)
#' lim$vi <- 1 / (lim$N - 3)
#' model <- metafor::rma.mv(
#'   yi = yi, V = vi, mods = ~ Environment * year,
#'   random = list(~ 1 | Article, ~ 1 | Datapoint), data = na.omit(lim)
#' )
#' test <- orchaRd::mod_results(model, mod = "year", group = "Article", data = lim, weights = "prop", by = "Environment")
#' orchaRd::bubble_plot(test, mod = "year", group = "Article", legend.pos = "top.left")
#' # Or just using model directly
#' orchaRd::bubble_plot(model, mod = "year", legend.pos = "top.left", group = "Article", weights = "prop", by = "Environment")
#' }
#' @export

# TODO - make poly works for bubble???
# TODO - write to https://github.com/rvlenth/emmeans/issues (missing combinations or interaction not allowed)

bubble_plot <- function(
  object,
  mod,
  group = NULL,
  xlab = "Moderator",
  ylab = "Effect size",
  N = "none",
  alpha = 0.5,
  cb = TRUE,
  k = TRUE,
  g = FALSE,
  transfm = c("none", "tanh", "invlogit", "percent", "percentr"),
  n_transfm = NULL,
  est.lwd = 1,
  ci.lwd = 0.5,
  pi.lwd = 0.5,
  est.col = "black",
  ci.col = "black",
  pi.col = "black",
  legend.pos = c("top.left", "top.right", "bottom.right", "bottom.left",
                 "top.out", "bottom.out", "none"),
  k.pos = c("top.right", "top.left",
            "bottom.right", "bottom.left",
            "none"),
  condition.nrow = 2,
  condition.order = NULL,
  weights = "prop",
  by = NULL,
  at = NULL
) {
  transfm <- match.arg(NULL, choices = transfm)
  legend.pos <- match.arg(NULL, choices = legend.pos)
  k.pos <- match.arg(NULL, choices = k.pos)

  if (missing(object)) stop("'object' must be provided. See ?bubble_plot", .call = FALSE)
  if (missing(group)) {
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?bubble_plot")
  }

  # Get results to be used. Only accepts metafor or orchard object objects
  if (any(class(object) %in% c("robust.rma", "rma.mv", "rma", "rma.uni"))) {
    results <- orchaRd::mod_results(
      object, 
      mod = mod,  
      group = group, 
      by = by, 
      at = at, 
      weights = weights
    )
  } else if (any(class(object) %in% c("orchard"))) {
    results <- object
  } else {
    stop("object argument must be a 'metafor' model or the output from 'mod_results'", .call = FALSE)
  }

  # Transform results if needed
  if (transfm != "none") {
    results <- .apply_transformation(results, n_transfm, transfm)
    mod_table <- results$mod_table
    data_trim <- results$data_trim
  }

  # Set condition in mod_table and data_trim.
  # If no condition, it uses a dummy variable. Make it easier to handle
  # reordering, facets and all that, because always reference to 'condition'.
  results <- .set_condition(results, condition.order)
  mod_table <- results$mod_table
  data_trim <- results$data

  # Define scale: size of points
  scale <- .set_scale(N, data_trim)
  data_trim$scale <- scale$scale
  size_legend <- scale$size_legend

  # Get k/g labels for the plot based on the processed data_trim.
  kg_labels <- .get_kg_labels(data_trim)

  # Note: the bbp (bubble plot) prefix is to avoid clashes with other functions
  plt <- .base_bubble_plot(data_trim, alpha) +
    .bbp_theme() +
    .bbp_pred_interval(mod_table, pi.lwd, pi.col) +
    .bbp_conf_interval(mod_table, ci.lwd, ci.col) +
    .bbp_estimate_line(mod_table, est.lwd, est.col) +
    .bbp_axis_labels(xlab, ylab) +
    .bbp_legends(legend.pos, size_legend) +
    .bbp_kg_labels(k, g, k.pos, kg_labels) + 
    .bbp_facets(data_trim, condition.nrow, condition.order) +
    .bbp_colors(data_trim, cb)

  return(plt)
}


#' Apply transformation to results
#' @keywords internal
.apply_transformation <- function(results, n_transfm, transfm) {
  mod_table <- results$mod_table
  data_trim <- results$data

  numeric_cols <- c("estimate", "lowerCL", "upperCL", "lowerPR", "upperPR")
  mod_table[, numeric_cols] <- transform_data(mod_table[, numeric_cols],
                                              n = n_transfm,
                                              transfm = transfm)
  data_trim$yi <- transform_data(data_trim$yi, n = n_transfm, transfm = transfm)

  return(list(mod_table = mod_table, data_trim = data_trim))
}


#' Set the condition variable for plotting.
#' @keywords internal
.set_condition <- function(results, condition.order) {
  # There are 2 cases:
  #  - No condition. 
  #  - Categorical condition
  # If the condition is cuantitative it is not possible to plot bubbles (points)

  mod_table <- results$mod_table
  data <- results$data

  # This is not necessary but helps with readability
  if (is.null(mod_table$condition)) {
    condition_type <- "none"
  } else if (is.character(mod_table$condition) || is.factor(mod_table$condition)) {
    condition_type <- "categorical"
  } else {
    stop("The condition must be categorical", call. = FALSE)
  }

  if (condition_type == "none") {
    data$condition <- factor(1)
    mod_table$condition <- factor(1)
  } else if (condition_type == "categorical") {
    if (is.null(condition.order)) {
      condition.order <- levels(factor(data$condition))
    } 
    mod_table$condition <- factor(mod_table$condition,
                                  levels = condition.order,
                                  labels = condition.order,
                                  ordered = TRUE)
    data$condition <- factor(data$condition,
                             levels = condition.order,
                             labels = condition.order,
                             ordered = TRUE)
  } 

  results <- list(
    mod_table = mod_table,
    data = data
  )

  return(results)
}


#' Set scale for bubble plot.
#' @keywords internal
.set_scale <- function(N, data_trim) {
  if (any(N != "none")) {
    scale <- data_trim$N
    size_legend <- paste0("Sample Size ($\\textit{N}$)")
  } else {
    scale <- (1 / sqrt(data_trim[, "vi"]))
    size_legend <- "Precision (1/SE)"
  }
  return(list(scale = scale, size_legend = size_legend))
}


#' Compute k and g labels for each condition.
#' @keywords internal
.get_kg_labels <- function(data_trim) {
  cond <- data_trim$condition

  K <- by(data_trim, cond, function(x) length(x[, "yi"]))
  G <- by(data_trim, cond, function(x) length(unique(x[, "stdy"])))

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


#' Create the base bubble plot layer.
#' @keywords internal
.base_bubble_plot <- function(data_trim, alpha) {
  p <- ggplot2::ggplot() + 
    ggplot2::geom_point(
      data = data_trim,
      ggplot2::aes(x = moderator,
                   y = yi,
                   size = scale,
	                 fill = condition),
		shape = 21,
		alpha = alpha)

  return(p)
}


#' Bubble plot ggplot2 theme settings.
#' @keywords internal
.bbp_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      legend.direction = "horizontal", 
      legend.background = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 10, colour = "black", hjust = 0.5, angle = 90),
      strip.text = ggplot2::element_text(size = 11)  # facet-label font size
  )
   
}


#' Add prediction interval lines to bubble plot.
#' @keywords internal
.bbp_pred_interval <- function(mod_table, pi.lwd, pi.col) {
  # Note: Multiple ggplot2::geoms must be passed as a list
  list(
    ggplot2::geom_smooth(data = mod_table,
      ggplot2::aes(x = moderator,
                   y = lowerPR),
        method = "loess",
        formula = y ~ x,
        se = FALSE,
        lty = "dotted",
        lwd = pi.lwd,
        colour = pi.col),
    ggplot2::geom_smooth(data = mod_table,
      ggplot2::aes(x = moderator,
                   y = upperPR),
        method = "loess",
        formula = y ~ x,
        se = FALSE,
        lty = "dotted",
        lwd = pi.lwd,
        colour = pi.col)
  )
}


#' Add confidence interval lines to bubble plot.
#' @keywords internal
.bbp_conf_interval <- function(mod_table, ci.lwd, ci.col) {
  # Note: Multiple ggplot2::geoms must be passed as a list
  list(
    ggplot2::geom_smooth(data = mod_table,
      ggplot2::aes(x = moderator,
                   y = lowerCL),
  		  method = "loess",
  		  formula = y ~ x,
  		  se = FALSE,
  		  lty = "dashed",
  		  lwd = ci.lwd,
  		  colour = ci.col),
    ggplot2::geom_smooth(data = mod_table,
      ggplot2::aes(x = moderator,
                   y = upperCL),
  		  method = "loess",
  		  formula = y ~ x,
  		  se = FALSE,
  		  lty = "dashed",
  		  lwd = ci.lwd,
  		  colour = ci.col)
    )
}


#' Add estimate line to bubble plot.
#' @keywords internal
.bbp_estimate_line <- function(mod_table, est.lwd, est.col) {
  ggplot2::geom_smooth(data = mod_table,
    ggplot2::aes(x = moderator, 
    y = estimate),
                 method = "loess",
                 formula = y ~ x,
                 se = FALSE,
                 lwd = est.lwd,
                 colour = est.col)
}


#' Add axis labels to bubble plot.
#' @keywords internal
.bbp_axis_labels <- function(xlab, ylab) {
  ggplot2::labs(x = xlab, y = ylab, parse = TRUE)
}


#' Add and position legends on bubble plot.
#' @keywords internal
.bbp_legends <- function(legend.pos, size_legend) {
  # Returns two ggplot layers:
  # One modifies the legends positions
  # The other are the legends
  legend_position_layer <- switch(legend.pos,
    "bottom.right" = ggplot2::theme(legend.position = "inside", legend.justification = c(1, 0)),
    "bottom.left"  = ggplot2::theme(legend.position = "inside", legend.justification = c(0, 0)),
    "top.right"    = ggplot2::theme(legend.position = "inside", legend.justification = c(1, 1)),
    "top.left"     = ggplot2::theme(legend.position = "inside", legend.justification = c(0, 1)),
    "top.out"      = ggplot2::theme(legend.position = "top"),
    "bottom.out"   = ggplot2::theme(legend.position = "bottom"),
    "none"         = ggplot2::theme(legend.position = "none")
  )

  legends_layer <- ggplot2::guides(
    size = ggplot2::guide_legend(title = latex2exp::TeX(size_legend)),
    colour = "none",
    fill = "none"
    )

  return(list(legend_position_layer, legends_layer))
}


#' Add k and g annotation labels to bubble plot.
#' @keywords internal
.bbp_kg_labels <- function(k, g, k.pos, kg_labels) {
  if (k == FALSE || k.pos == "none") return()

  # Build label string
  if (g) {
    label_str <- paste0("italic(k)==", kg_labels$K, " ~ (", kg_labels$G, ")")
  } else {
    label_str <- paste0("italic(k)==", kg_labels$K)
  }
  kg_labels$label_str <- label_str

  # If left, anchor at -inf. At right, at +inf. Same for bottom and top
  x_val <- if (k.pos %in% c("top.left", "bottom.left")) -Inf else Inf
  y_val <- if (k.pos %in% c("top.left", "top.right")) Inf else -Inf

  # Set justification
  hjust_val <- if (k.pos %in% c("top.left", "bottom.left")) -0.5 else 1.5
  vjust_val <- if (k.pos %in% c("top.left", "top.right")) 2.5 else -1.5

  ggplot2::geom_text(
    data = kg_labels,
    ggplot2::aes(x = x_val, y = y_val, label = label_str),
    parse = TRUE,
    hjust = hjust_val,
    vjust = vjust_val
  )
}


#' Add facetting to bubble plot by condition.
#' @keywords internal
.bbp_facets <- function(data_trim, condition.nrow, condition.order) {
  condition <- data_trim$condition

  if (nlevels(condition) == 1) return()

  if (is.character(condition) || is.factor(condition)) {
    ggplot2::facet_wrap(ggplot2::vars(condition), nrow = condition.nrow) 
  }
}


#' Set fill colors for bubble plot points.
#' @keywords internal
.bbp_colors <- function(data_trim, cb) {
  # # colour blind friendly colours with grey
  cbpl <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733",
            "#332288", "#AA4499", "#44AA99", "#999933",
            "#882255", "#661100", "#6699CC", "#888888",
            "#E69F00", "#56B4E9", "#009E73", "#F0E442",
            "#0072B2", "#D55E00", "#CC79A7", "#999999")

  # Boolean variable to make it readable
  has_conditions <- if (nlevels(data_trim$condition) > 1) TRUE else FALSE

  if (!has_conditions) {
    ggplot2::scale_fill_manual(values = "grey90")
  } else if (has_conditions && cb == TRUE) {
    ggplot2::scale_fill_manual(values = cbpl) 
  } else {
    ggplot2::scale_fill_hue()  # Default in ggplot2 for discrete variables
  }

}
