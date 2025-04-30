#' @title bubble_plot
#' @description Using a \pkg{metafor} model object of class \code{rma} or \code{rma.mv}, or a results table of class \code{orchard}, the \code{bubble_plot} function creates a bubble plot from slope estimates. In cases when a model includes interaction terms, this function creates panels of bubble plots.
#' @param object Model object of class \code{rma}, \code{rma.mv}, or \code{orchard} table of model results
#' @param mod The name of a continuous moderator, to be plotted on the x-axis of the bubble plot.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species, or any grouping variable one wishes to present sample sizes for. Not needed if an \code{orchard_plot} is provided with a \code{mod_results} object of class \code{orchard}.
#' @param by Character vector indicating the name that predictions should be conditioned on for the levels of the moderator.
#' @param at List of levels one wishes to predict at for the corresponding variables in \code{by}. Used when one wants marginalised means. This argument can also be used to suppress levels of the moderator when argument \code{subset = TRUE}. Provide a list as follows: \code{list(mod = c("level1", "level2"))}.
#' @param weights How to marginalize categorical variables; used when one wants marginalised means. The default is \code{weights = "prop"}, which weights means for moderator levels based on their proportional representation in the data. For example, if \code{"sex"} is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal when, for example, males and females are typically roughly equally prevalent in a population. In cases such as these, you can give the moderator levels equal weight using \code{weights = "equal"}.
#' @param transfm If set to \code{"tanh"}, a tanh transformation will be applied to effect sizes, converting Zr to a correlation or pulling in extreme values for other effect sizes (lnRR, lnCVR, SMD).  \code{"invlogit"} can be used to convert lnRR to the inverse logit scale. \code{"percentr"} can convert to the percentage change scale when using response ratios and \code{"percent"} can convert to the percentage change scale of an log transformed effect size. Defaults to \code{"none"}.
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
#' @param legend.pos Where to place the legend, or not to include a legend ("none").
#' @param condition.levels Order of the levels of the condition variable in the order to plot. Defaults to NULL.
#'
#' @return Orchard plot
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
#' orchaRd::bubble_plot(test, mod = "year", data = lim, group = "Article", legend.pos = "top.left")
#' # Or just using model directly
#' orchaRd::bubble_plot(model, mod = "year", legend.pos = "top.left", data = lim, group = "Article", weights = "prop", by = "Environment")
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
  # condition.lab = "Condition",
  weights = "prop",
  by = NULL,
  at = NULL,
  condition.levels = NULL
) {
  transfm <- match.arg(NULL, choices = transfm)
  legend.pos <- match.arg(NULL, choices = legend.pos)
  k.pos <- match.arg(NULL, choices = k.pos)
  # facet <- match.arg(NULL, choices = facet)

  if (missing(group)) {
    stop("Please specify the 'group' argument by providing the name of the grouping variable. See ?bubble_plot")
  }

  if (is.numeric(by)) {
    k <- FALSE
    g <- FALSE
  }


  if (any(class(object) %in% c("robust.rma", "rma.mv", "rma", "rma.uni"))) {
    if (mod != "1") {
      results <- orchaRd::mod_results(object, mod, group,
        by = by, at = at, weights = weights
      )
    } else {
      results <- orchaRd::mod_results(object,
        mod = "1", group,
        by = by, at = at, weights = weights
      )
    }
  }

  if (any(class(object) %in% c("orchard"))) {
    results <- object
  }

  # ---------------------------
  # Prepare data for plotting 

  # Get model table and raw data
  mod_table <- results$mod_table
  data_trim <- results$data

  # Define scale: size of points
  if (any(N != "none")) {
    data_trim$scale <- data_trim$N
    size_legend <- paste0("Sample Size ($\\textit{N}$)") # we want to use italic
  } else {
    data_trim$scale <- (1 / sqrt(data_trim[, "vi"]))
    size_legend <- "Precision (1/SE)"
  }

  # Transform data if needed
  if (transfm != "none") {
    numeric_cols <- c("estimate", "lowerCL", "upperCL", "lowerPR", "upperPR")
    mod_table[, numeric_cols] <- transform_data(mod_table[, numeric_cols],
                                                n = n_transfm,
                                                transfm = transfm)
    data_trim$yi <- transform_data(data_trim$yi, n = n_transfm, transfm = transfm)
  }


  # NOTE: All this sections is just for creating labels:
  # k, g and condition
  #
  # Don't think it should be here. It can be computed inside
  # a function that takes care of those labels

  if (is.null(data_trim$condition) == TRUE) {
    # the number of effect sizes
    effect_num <- nrow(data_trim)
    # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
    group_num <- length(unique(data_trim$stdy))
    dat_text <- data.frame(K = effect_num, G = group_num)
  } else {
    if (!is.null(condition.levels)) {
      data_trim$condition <- factor(data_trim$condition, levels = condition.levels, labels = condition.levels, ordered = TRUE)
    } else {
      # making sure factor names match
      data_trim$condition <- factor(data_trim$condition, levels = mod_table$condition, labels = mod_table$condition)
    }

    effect_num <- as.vector(by(data_trim, data_trim[, "condition"], function(x) base::length(x[, "yi"])))
    # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
    group_num <- as.vector(by(data_trim, data_trim[, "condition"], function(x) base::length(base::unique(x[, "stdy"]))))
    dat_text <- data.frame(K = effect_num,
			   G = group_num,
			   condition = as.vector(base::levels(data_trim$condition)),
                           stringsAsFactors = FALSE)
    dat_text$condition <- factor(dat_text$condition,
                                 levels = levels(data_trim$condition),
                                 ordered = is.ordered(data_trim$condition))
  }


  if (is.null(data_trim$condition) == TRUE) {
    fill <- factor(1)
    fill_palette <- ggplot2::scale_fill_manual(values = "grey90")
    # Remove guides for fill!
  } else if (is.character(data_trim$condition) == TRUE || is.factor(data_trim$condition) == TRUE) {
    fill <- data_trim$condition
    fill_palette <- ggplot2::scale_fill_viridis_d()
    # TODO: Don't use viridis, set colour blind
  }


  # Note: the bbp (bubble plot) prefix is to avoid clashes with other functions
  plt <- .base_bubble_plot(data_trim, fill, fill_palette, alpha) +
	  .bbp_theme() +
	  .bbp_pred_interval(mod_table, pi.lwd, pi.col) +
	  .bbp_conf_interval(mod_table, ci.lwd, ci.col) +
	  .bbp_estimate_line(mod_table, est.lwd, est.col) +
	  .bbp_axis_labels(xlab, ylab) +
	  .bbp_legends(legend.pos, size_legend) +
	  .bbp_k_and_g(k, g, k.pos, dat_text) + 
	  .bbp_facets(data_trim$condition, condition.nrow, condition.levels) 



  return(plt)

  # putting colors in
  # # colour blind friendly colours with grey
  cbpl <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")


  if(cb == TRUE){
    plot <- plot +
      ggplot2::scale_fill_manual(values=cbpl) +
      ggplot2::scale_colour_manual(values=cbpl)
  }

  return(plot)
}


.base_bubble_plot <- function(data_trim, fill, fill_palette, alpha) {
  p <- ggplot2::ggplot() + 
    ggplot2::geom_point(
      data = data_trim,
      ggplot2::aes(
        x = moderator,
	y = yi,
	size = scale,
	fill = fill
      ),
      shape = 21,
      alpha = alpha
    ) +
    fill_palette

  return(p)
}


.bbp_theme <- function() {
  ggplot2::theme_bw() +
    ggplot2::theme(
      legend.direction = "horizontal", 
      legend.background = ggplot2::element_blank(),
      axis.text.y = ggplot2::element_text(size = 10, colour = "black", hjust = 0.5, angle = 90)
    )
}


.bbp_pred_interval <- function(mod_table, pi.lwd, pi.col) {
  # Multiple geoms should be passed as a list
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


.bbp_conf_interval <- function(mod_table, ci.lwd, ci.col) {
  # Multiple geoms should be passed as a list:
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


.bbp_axis_labels <- function(xlab, ylab) {
  ggplot2::labs(x = xlab,
	       	y = ylab,
	       	parse = TRUE) 
}


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

  # TODO: The legends in the old version show the circles filled, not blank!

  legends_layer <- ggplot2::guides(
    size = ggplot2::guide_legend(title = latex2exp::TeX(size_legend)),
    colour = "none",
    fill = "none"
    )

  return(list(legend_position_layer, legends_layer))
}


.bbp_facets <- function(condition, condition.nrow, condition.levels) {
  if (!is.null(condition.levels)) {
    condition <- factor(condition, levels = condition.levels)
  }

  if (is.character(condition) == TRUE || is.factor(condition) == TRUE) {
    ggplot2::facet_wrap(ggplot2::vars(condition), nrow = condition.nrow) 
  }
}


.bbp_k_and_g <- function(k, g, k.pos, dat_text) {
  if (k == FALSE || k.pos == "none") {
    return()
  } 

  # Build label string
  if (g) {
    label_str <- paste0("italic(k)==", dat_text$K, " ~ (", dat_text$G, ")")
  } else {
    label_str <- paste0("italic(k)==", dat_text$K)
  }

  # If left, anchor at -inf. At right, at +inf. Same for bottom and top
  x_val <- if (k.pos %in% c("top.left", "bottom.left")) -Inf else Inf
  y_val <- if (k.pos %in% c("top.left", "top.right")) Inf else -Inf

  # Set justification
  hjust_val <- if (k.pos %in% c("top.left", "bottom.left")) -0.5 else 2
  vjust_val <- if (k.pos %in% c("top.left", "top.right")) 2.5 else -1.5
  if (g) {
    # Adjust a little bit the position when adding g because of extra space needed
    hjust_val <- hjust_val - 0.5
    vjust_val <- vjust_val - 0.5
  }

  ggplot2::geom_text(
    data = dat_text,
    mapping = ggplot2::aes(x = x_val, y = y_val),
    inherit.aes = FALSE  # Need this for silencing ggplot warnings 
    label = label_str,
    parse = TRUE,
    hjust = hjust_val,
    vjust = vjust_val
  )
}
