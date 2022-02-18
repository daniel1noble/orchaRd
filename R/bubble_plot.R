#' @title bubble_plot
#' @description Using a metafor model object of class rma or rma.mv or a results table of class orchard, it creates a an bubble plot from slope estimates or panels of bubble plot in cases when a model includes interaction terms.
#' @param object model object of class 'rma.mv', 'rma' or 'orchard' table of model results
#' @param mod the name of a contionous moderator.
#' @param group The grouping variable that one wishes to plot beside total effect sizes, k. This could be study, species or whatever other grouping variable one wishes to present sample sizes. Not needed of a orchard_plot is provided with a mod_results object of class 'orchard'.
#' @param data The data frame used to fit the rma.mv model object. Not needed of a orchard_plot is provided with a mod_results object of class 'orchard'.
#' @param by Used when one wants marginalised means. The 'condition' variable that one wishes to have the mean for the moderator vary.
#' @param at Used when one wants marginalised means. The 'condition' that one wishes to calculate the means at, but is not presented in output
#' @param weights Used when one wants marginalised means. How to marginalize categorical variables. The default is weights = "prop", which wights means for moderator levels based on their proportional representation in the data. For example, if "sex" is a moderator, and males have a larger sample size than females, then this will produce a weighted average, where males are weighted more towards the mean than females. This may not always be ideal. IN the case if sex, for example, males and females are roughly equally prevalent in a population. As such, you can give the moderator levels equal weight using weights = "equal".
#' @param xlab Moderator label.
#' @param xlab Effect size measure label.
#' @param N  The vector of sample size which an effect size is based on. If default, we use precision (the inverse of sampling standard error)
#' @param alpha The level of transparency for pieces of fruit (effect size)
#' @param cb If TRUE, it uses 20 colour blind friendly colors
#' @param k If TRUE, it displays k (number of effect sizes) on the plot
#' @param g If TRUE, it displays g (number of grouping levels for each level of the moderator) on the plot
#' @param condition.lab Label for the condition being marginalized over.
#' @return Orchard plot
#' @author Shinichi Nakagawa - s.nakagawa@unsw.edu.au
#' @author Daniel Noble - daniel.noble@anu.edu.au
#' @examples
#' \dontrun{
#' data(eklof)
#'
#' }
#' @export

bubble_plot <- function(object, mod, group, data,
                         xlab = "Moderator", ylab = "Effect size", N = "none",
                         alpha = 0.5, cb = FALSE, k = TRUE, g = TRUE,
                         condition.lab = "Condition",
                         weights = "prop", by = NULL, at = NULL)
{

  if(any(class(object) %in% c("rma.mv", "rma"))){

    if(mod != "1"){
      results <-  orchaRd::mod_results(object, mod, group, data,
                                       by = by, at = at, weights = weights)
    } else {
      results <-  orchaRd::mod_results(object, mod = "1", group, data,
                                       by = by, at = at, weights = weights)
    }
  }

  if(any(class(object) %in% c("orchard"))) {
    results <- object
  }

  mod_table <- results$mod_table

  data <- results$data
  #data_trim$moderator <- factor(data_trim$moderator, levels = mod_table$name, labels = mod_table$name)

  data$scale <- (1/sqrt(data[,"vi"]))
  legend <- "Precision (1/SE)"
  #legend <- "Precision (1/SE)"

  if(any(N != "none")){
    data$scale <- N
    legend <- paste0("Sample Size (", "N",")") # we want to use italic
  }


  # Add in total effect sizes for each level
  mod_table$K <- length(data[,"moderator"])

  # Add in total levels of a grouping variable (e.g., study ID) within each moderator level.
  mod_table$g <-length(unique(data[,"stdy"]))

  #data_no <- nrow(data)

  # colour blind friendly colours with grey
  cbpl <- c("#88CCEE", "#CC6677", "#DDCC77", "#117733", "#332288", "#AA4499", "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#999999")


  plot <-ggplot() +
    geom_point(data = data, aes(x = moderator, y = yi, size = scale, colour = condition, fill = condition), shape = 21, alpha = alpha) +
    #geom_ribbon(aes(ymin = ymin, ymax = ymax), fill = "#0072B2")  + # not quite sure why this does not work
    geom_smooth(data = mod_table, aes(x = moderator, y = lowerPR), method =  "loess", se = FALSE, lty =  "dotted", lwd = 0.25, colour = "#0072B2") +
    geom_smooth(data = mod_table, aes(x = moderator, y = upperPR), method =  "loess", se = FALSE, lty = "dotted", lwd = 0.25, colour = "#0072B2") +
    geom_smooth(data = mod_table, aes(x = moderator, y = lowerCL), method =  "loess", se = FALSE,lty = "dotted", lwd = 0.25, colour ="#D55E00") +
    geom_smooth(data = mod_table, aes(x = moderator, y = upperCL), method =  "loess", se = FALSE, lty ="dotted", lwd = 0.25, colour ="#D55E00") +
    geom_smooth(data = mod_table, aes(x = moderator, y = estimate), method =  "loess", se = FALSE, lty ="dashed", lwd = 0.5, colour ="black") +
    facet_grid(rows = vars(condition)) +
    labs(x = xlab, y = ylab, size = legend, parse = TRUE) +
    guides(fill = "none", colour = "none") +
    # themses
    theme_bw() +
    theme(legend.position= c(1, 1), legend.justification = c(1, 1)) +
    theme(legend.direction="horizontal") +
    #theme(legend.background = element_rect(fill = "white", colour = "black")) +
    theme(legend.background = element_blank()) +
    theme(axis.text.y = element_text(size = 10, colour ="black", hjust = 0.5, angle = 90))



  # putting colors in
  if(cb == TRUE){
    plot <- plot +
      ggplot2::scale_fill_manual(values=cbpl) +
      ggplot2::scale_colour_manual(values=cbpl)
  }



  return(plot)
}
