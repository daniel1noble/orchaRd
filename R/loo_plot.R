#' Leave-One-Out plot
#'
#' Create a leave-one-out sensitivity plot for a meta-analytic model.
#'
#' @param model A metafor model from class \code{rma}, \code{rma.mv}, or \code{rma.uni}.
#' @param loo_dataframe A data frame containing leave-one-out results. It must include the columns:
#'   \code{left_out} (paper identifier), \code{b} (effect size), \code{ci_lb}, and \code{ci_ub}.
#' @param order A character string specifying the order for the \code{left_out} factor.
#'   Options are \code{"none"}, \code{"alphabetic"}, \code{"ascending"}, or \code{"descending"}.
#'   \code{"ascending"} and \code{"descending"} order the factor by the value of \code{b}.
#' @param labels Optional. A data frame mapping the raw codes in \code{left_out} to friendlier names.
#'   It must have two columns:
#'   \itemize{
#'     \item \code{left_out}: the codes as they appear in \code{loo_dataframe}.
#'     \item \code{label}: the friendly label you want to display.
#'   }
#'
#' @return A ggplot2 object.
#' @export
#'
#' @examples
#' \dontrun{
#'   # Load necessary libraries
#'   library(metafor)
#'   library(dplyr)
#'
#'   # Calculate effect sizes using the metafor package
#'   dat <- metafor::escalc(measure = "RR", ai = tpos, bi = tneg, ci = cpos, di = cneg, data = dat.bcg)
#'
#'   # Fit the meta-analytic model
#'   res <- metafor::rma(yi, vi, data = dat)
#'
#'   # Create a mapping of trial codes to friendly labels.
#'   # This step could also be done in a spreadsheet and read in as a CSV.
#'   labels_df <- dat.bcg %>%
#'     dplyr::distinct(trial, .keep_all = TRUE) %>%
#'     dplyr::mutate(left_out = trial,
#'                   label = paste(author, year, sep = ", ")) %>%
#'     dplyr::select(left_out, label)
#'
#'   # Run the leave-one-out analysis (assuming a leave_one_out function exists)
#'   res_loo <- leave_one_out(res, group = "trial")
#'
#'   # Create the leave-one-out plot using custom labels and descending order.
#'   loo_plot(res, res_loo, order = "descending", labels = labels_df)
#' }

loo_plot <- function(model, loo_dataframe, 
                     order = c("none", "alphabetic", "ascending", "descending"),
                     labels = NULL) {
  
  order <- match.arg(order)
  
  # If a 'labels' data frame is provided, merge it with the main data.
  if (!is.null(labels)) {
    if (!is.data.frame(labels) || !all(c("left_out", "label") %in% names(labels))) {
      stop("'labels' must be a data frame with columns 'left_out' and 'label'. See ?loo_plot")
    }
    # Merge the labels data frame and replace left_out values with the friendly labels.
    loo_dataframe <- merge(loo_dataframe, labels, by = "left_out", all.x = TRUE)
    loo_dataframe$left_out <- as.character(loo_dataframe$label)
  } else {
    # Default cleaning if no mapping is provided: replace underscores with spaces and convert to title case.
    loo_dataframe$left_out <- as.character(loo_dataframe$left_out)
    loo_dataframe$left_out <- gsub("_", " ", loo_dataframe$left_out)
    loo_dataframe$left_out <- tools::toTitleCase(loo_dataframe$left_out)
  }
  
  # Ensure left_out is a factor for ordering.
  loo_dataframe$left_out <- factor(loo_dataframe$left_out)
  
  # Apply ordering if requested.
  if (order == "alphabetic") {
    loo_dataframe$left_out <- factor(loo_dataframe$left_out,
                                     levels = sort(levels(loo_dataframe$left_out)))
  } else if (order == "ascending") {
    new_levels <- unique(loo_dataframe$left_out[order(loo_dataframe$b)])
    loo_dataframe$left_out <- factor(loo_dataframe$left_out, levels = new_levels)
  } else if (order == "descending") {
    new_levels <- unique(loo_dataframe$left_out[order(-loo_dataframe$b)])
    loo_dataframe$left_out <- factor(loo_dataframe$left_out, levels = new_levels)
  }
  
  # Create the plot.
  p <- ggplot2::ggplot(loo_dataframe, ggplot2::aes(x = left_out,
                                                    y = b,
                                                    ymin = ci_lb,
                                                    ymax = ci_ub)) +
    # Add overall effect and its confidence intervals.
    ggplot2::geom_hline(yintercept = model$b[1], linetype = "solid") +
    ggplot2::geom_hline(yintercept = model$ci.lb[1], linetype = "dotted") +
    ggplot2::geom_hline(yintercept = model$ci.ub[1], linetype = "dotted") +
    # Plot leave-one-out effect sizes with their confidence intervals.
    ggplot2::geom_pointrange(alpha = 0.7, color = "darkorange") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1),
                   panel.grid = ggplot2::element_blank()) +
    ggplot2::coord_flip()
  
  return(p)
}

